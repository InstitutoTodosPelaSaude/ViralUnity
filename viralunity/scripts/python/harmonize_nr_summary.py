#!/usr/bin/env python3
"""Harmonize a contig-track taxa summary with the per-contig NR verdict table.

Refactored from the REVISA nr_validation prototype for in-pipeline use:

* Keys are ``sample|contig`` (single run), not ``RUN|sample|contig``.
* Works for both contig tracks — pass the track's own per-sample
  ``(contig, taxid)`` krona map (``*.diamond.supported.krona_input.tsv`` for
  diamond_contigs, ``*.output.krona.txt`` for kraken2_contigs).
* It is a *taxonomic filter*: species rows the NR consensus confidently calls
  non-viral (``nr_is_virus == False``) are REMOVED (``nr_pass = False``) so the
  downstream bleed/negative-control statistics see only NR-validated taxa. Rows
  that are ambiguous / have no NR data, and all non-species rows, are kept
  (``nr_pass = NA``). The informational columns ``nr_is_virus``,
  ``nr_species_correct``, ``nr_correct_species`` are appended to every row.

Species-level RefSeq-vs-NR disagreements are surfaced (not filtered) in a
separate ``*_nr_flags.tsv`` (reason ``misid_novel`` / ``misid_known``).
"""

import argparse
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import pandas as pd

try:
    from viralunity.scripts.python.taxonomy import get_lineage, load_taxdump
except ImportError:
    from taxonomy import get_lineage, load_taxdump

NA = "NA"
VIRAL_PHYLUM_SUFFIX = "viricota"


def majority_bool(n_true: int, n_total: int):
    """True/False/'ambiguous' by strict majority; None when n_total == 0."""
    if n_total == 0:
        return None
    if n_true * 2 > n_total:
        return True
    if n_true * 2 < n_total:
        return False
    return "ambiguous"


def _bool_str(v) -> str:
    if v is None:
        return NA
    if v is True:
        return "True"
    if v is False:
        return "False"
    return str(v)


def load_nr_verdicts(nr_table_path: str) -> Dict[Tuple[str, str], dict]:
    """Index the NR top-species verdict table by (sample, contig).

    Each qseqid is ``sample|contig`` (combine_contigs re-heading). One row per
    contig (filter_top_species_hits emits one representative per query).
    """
    verdicts: Dict[Tuple[str, str], dict] = {}
    if not os.path.exists(nr_table_path) or os.path.getsize(nr_table_path) == 0:
        return verdicts
    df = pd.read_csv(nr_table_path, sep="\t", dtype=str).fillna(NA)
    for _, row in df.iterrows():
        qseqid = str(row.get("qseqid", ""))
        bits = qseqid.split("|", 1)
        if len(bits) != 2:
            continue
        sample, contig = bits
        phylum = str(row.get("phylum", NA))
        verdicts[(sample, contig)] = {
            "species": str(row.get("species", NA)),
            "phylum": phylum,
            "is_viral": phylum != NA and phylum.lower().endswith(VIRAL_PHYLUM_SUFFIX),
            "pident": str(row.get("pident", NA)),
            "evalue": str(row.get("evalue", NA)),
            "bitscore": str(row.get("bitscore", NA)),
        }
    return verdicts


def build_species_contigs(
    krona_files: List[Tuple[str, str]],
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
    name_map: Dict[str, str],
):
    """Map each contig's track taxid up to its species; index contigs by
    (sample, species_taxid). ``krona_files`` is a list of (sample, path) where
    each path is a ``contig<TAB>taxid`` TSV."""
    species_cache: Dict[str, Tuple] = {}

    def climb_species(taxid: str):
        if taxid in species_cache:
            return species_cache[taxid]
        sp_taxid, sp_name = None, None
        if taxid and taxid != "0" and taxid in parent_map:
            for tid in get_lineage(taxid, parent_map):
                if rank_map.get(tid) == "species":
                    sp_taxid, sp_name = tid, name_map.get(tid, NA)
                    break
        species_cache[taxid] = (sp_taxid, sp_name)
        return sp_taxid, sp_name

    contig_info: Dict[Tuple[str, str], dict] = {}
    species_contigs: Dict[Tuple[str, str], List[str]] = defaultdict(list)
    for sample, path in krona_files:
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 2:
                    continue
                contig, hit_taxid = cols[0], cols[1]
                sp_taxid, sp_name = climb_species(hit_taxid)
                contig_info[(sample, contig)] = {
                    "hit_taxid": hit_taxid,
                    "sp_name": sp_name,
                }
                if sp_taxid is not None:
                    species_contigs[(sample, sp_taxid)].append(contig)
    return contig_info, species_contigs


def aggregate_species_row(sample, species_taxid, species_name, contigs, verdicts):
    """Compute (nr_is_virus, nr_species_correct, nr_correct_species) for a row."""
    n_with_nr = n_viral = 0
    viral_species = []
    for contig in contigs:
        rec = verdicts.get((sample, contig))
        if rec is None:
            continue
        n_with_nr += 1
        if rec["is_viral"]:
            n_viral += 1
            viral_species.append(rec["species"])

    is_virus = majority_bool(n_viral, n_with_nr)
    species_correct = None
    correct_species = NA
    if is_virus is True:
        matches = [s for s in viral_species if s != NA and s.lower() == species_name.lower()]
        species_correct = majority_bool(len(matches), len(viral_species))
        if species_correct is False:
            mism = [s for s in viral_species if s != NA and s.lower() != species_name.lower()]
            if mism:
                correct_species = Counter(mism).most_common(1)[0][0]
    return {
        "nr_is_virus": _bool_str(is_virus),
        "nr_species_correct": _bool_str(species_correct) if is_virus is True else NA,
        "nr_correct_species": correct_species,
    }


def _nr_pass(nr_is_virus: str):
    """True keep / False drop / NA keep, from the nr_is_virus verdict string."""
    if nr_is_virus == "True":
        return "True"
    if nr_is_virus == "False":
        return "False"
    return NA


def harmonize(
    summary_path: str,
    nr_table_path: str,
    krona_files: List[Tuple[str, str]],
    output_path: str,
    dropped_path: str,
    flags_path: str,
    taxdump_dir: str,
) -> Tuple[int, int]:
    """Append NR columns + nr_pass to the summary, drop confidently-non-viral
    species rows, and write the dropped sidecar and the nr_flags file. Returns
    (n_kept, n_dropped)."""
    parent_map, rank_map, name_map = load_taxdump(
        os.path.join(taxdump_dir, "nodes.dmp"), os.path.join(taxdump_dir, "names.dmp")
    )
    verdicts = load_nr_verdicts(nr_table_path)
    contig_info, species_contigs = build_species_contigs(
        krona_files, parent_map, rank_map, name_map
    )

    if not os.path.exists(summary_path) or os.path.getsize(summary_path) == 0:
        open(output_path, "w").close()
        open(dropped_path, "w").close()
        _write_flags(flags_path, [])
        return 0, 0

    df = pd.read_csv(summary_path, sep="\t", dtype=str).fillna(NA)
    sample_summary_species = defaultdict(set)
    for _, r in df.iterrows():
        if r["rank"] == "species" and r["name"] != NA:
            sample_summary_species[r["sample"]].add(r["name"].lower())

    nr_is_virus_col, sc_col, cs_col, pass_col = [], [], [], []
    for _, row in df.iterrows():
        if row["rank"] != "species":
            nr_is_virus_col.append(NA)
            sc_col.append(NA)
            cs_col.append(NA)
            pass_col.append(NA)
            continue
        contigs = species_contigs.get((row["sample"], str(row["taxid"])), [])
        agg = aggregate_species_row(row["sample"], row["taxid"], row["name"], contigs, verdicts)
        nr_is_virus_col.append(agg["nr_is_virus"])
        sc_col.append(agg["nr_species_correct"])
        cs_col.append(agg["nr_correct_species"])
        pass_col.append(_nr_pass(agg["nr_is_virus"]))

    df["nr_is_virus"] = nr_is_virus_col
    df["nr_species_correct"] = sc_col
    df["nr_correct_species"] = cs_col
    df["nr_pass"] = pass_col

    keep_mask = df["nr_pass"] != "False"
    kept = df[keep_mask]
    dropped = df[~keep_mask]
    kept.to_csv(output_path, sep="\t", index=False)
    dropped.to_csv(dropped_path, sep="\t", index=False)

    _write_flags(
        flags_path,
        _build_flags(verdicts, contig_info, sample_summary_species),
    )
    return len(kept), len(dropped)


FLAG_COLS = [
    "sample",
    "contig",
    "reason",
    "refseq_species",
    "nr_species",
    "nr_species_in_sample_refseq",
    "nr_phylum",
    "pident",
    "evalue",
    "bitscore",
]


def _build_flags(verdicts, contig_info, sample_summary_species) -> List[list]:
    rows = []
    for (sample, contig), rec in sorted(verdicts.items()):
        if not rec["is_viral"] or rec["species"] == NA:
            continue
        info = contig_info.get((sample, contig))
        refseq_sp = info["sp_name"] if info and info["sp_name"] else NA
        if refseq_sp != NA and rec["species"].lower() == refseq_sp.lower():
            continue  # NR agrees with the refseq species -> not a flag
        in_sample = rec["species"].lower() in sample_summary_species.get(sample, set())
        reason = "misid_known" if in_sample else "misid_novel"
        rows.append(
            [
                sample,
                contig,
                reason,
                refseq_sp,
                rec["species"],
                "True" if in_sample else "False",
                rec["phylum"],
                rec["pident"],
                rec["evalue"],
                rec["bitscore"],
            ]
        )
    return rows


def _write_flags(flags_path: str, rows: List[list]) -> None:
    with open(flags_path, "w") as out:
        out.write("\t".join(FLAG_COLS) + "\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Harmonize a contig summary with NR verdicts.")
    p.add_argument("--summary", required=True)
    p.add_argument("--nr", required=True, help="NR top_species_hit_lca TSV")
    p.add_argument("--output", required=True)
    p.add_argument("--dropped", required=True)
    p.add_argument("--flags", required=True)
    p.add_argument("--taxdump-dir", required=True)
    p.add_argument("--krona", action="append", required=True, metavar="SAMPLE=PATH")
    return p.parse_args()


def _parse_krona_pairs(items):
    pairs = []
    for item in items:
        sample, sep, path = item.partition("=")
        if not sep:
            raise ValueError(f"--krona must be SAMPLE=PATH, got: {item}")
        pairs.append((sample, path))
    return pairs


def main() -> None:
    if "snakemake" in globals():
        samples = list(snakemake.params.samples)  # noqa: F821
        krona = list(snakemake.input.krona)  # noqa: F821
        krona_files = list(zip(samples, krona))
        kept, dropped = harmonize(
            summary_path=str(snakemake.input.summary),  # noqa: F821
            nr_table_path=str(snakemake.input.nr),  # noqa: F821
            krona_files=krona_files,
            output_path=str(snakemake.output.summary),  # noqa: F821
            dropped_path=str(snakemake.output.dropped),  # noqa: F821
            flags_path=str(snakemake.output.flags),  # noqa: F821
            taxdump_dir=str(snakemake.params.taxdump),  # noqa: F821
        )
    else:
        args = parse_args()
        kept, dropped = harmonize(
            summary_path=args.summary,
            nr_table_path=args.nr,
            krona_files=_parse_krona_pairs(args.krona),
            output_path=args.output,
            dropped_path=args.dropped,
            flags_path=args.flags,
            taxdump_dir=args.taxdump_dir,
        )
    print(f"[harmonize_nr_summary] kept={kept} dropped={dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
