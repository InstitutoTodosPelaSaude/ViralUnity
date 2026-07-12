#!/usr/bin/env python3
"""Add per-taxon largest-contig statistics to a diamond_contigs taxa summary.

For each taxon row, report the length of the largest de novo viral contig
assigned to that taxon, that contig's median per-position sequencing depth, and
what fraction of the reference genome its length spans:

  * ``largest_contig_bp``               — length (bp) of the largest assigned contig
  * ``largest_contig_ref_coverage_pct`` — largest_contig_bp / genome_length_bp * 100
  * ``largest_contig_median_depth``     — median depth over that single largest contig

These are a cheap proxy for "how much of the genome is covered, and how well",
reusing the viral read-remap BAM the diamond_contigs track already produces
(``samtools depth -a`` -> per-position depth; contig length = number of
positions). A contig is assigned to a taxon by climbing its leaf taxid to the
family/genus/species ancestor, mirroring how ``summarize_krona_taxa`` and
``harmonize_nr_summary`` aggregate hits up the lineage.

``largest_contig_ref_coverage_pct`` is a preliminary genome-completeness estimate
(e.g. a 7 kb largest contig for a 10 kb virus ≈ 70%), meant to help decide
whether a reference assembly is worth attempting. It reuses the
``genome_length_bp`` column already added by the RPKM step (the median RefSeq
genome length per ``(rank, taxid)`` node), so it is populated on the same tracks
and under the same ``--viral-genomes`` gate. It is reported raw and *uncapped*:
a value above 100% means the largest contig exceeds the median reference length
(over-assembly, a chimeric contig, or a strain longer than the RefSeq median).
It is ``NA`` where no contig is assigned or no reference length is available, and
is approximate at family/genus ranks (median length across a diverse node).

Caveat: contig length is a proxy for, not a direct measure of, the fraction of
the reference genome covered — a long contig can still be a partial or chimeric
assembly. Only available when ``--viral-genomes`` (RPKM) is enabled, which is
also where the diamond_contigs viral remap runs.
"""

import argparse
import os
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import pandas as pd

try:
    from viralunity.scripts.python.taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump
except ImportError:
    from taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump

NA = "NA"


def _ref_coverage_pct(contig_len: int, genome_length_raw) -> object:
    """Return largest_contig_bp / genome_length_bp * 100 (raw, uncapped, 2 dp).

    ``NA`` when the reference length is missing, non-numeric, NaN, or <= 0.
    Values > 100 are legitimate (the largest contig exceeds the median reference
    length) and are reported as-is.
    """
    try:
        genome_length = float(genome_length_raw)
    except (TypeError, ValueError):
        return NA
    if genome_length != genome_length or genome_length <= 0:  # NaN or non-positive
        return NA
    return round(contig_len / genome_length * 100.0, 2)


def load_contig_depths(depth_path: str) -> Dict[str, Tuple[int, float]]:
    """Parse a ``samtools depth -a`` file into per-contig (length, median_depth).

    The file has ``contig<TAB>pos<TAB>depth`` with one line per position (``-a``
    emits zero-depth positions too), so the number of lines for a contig is its
    length and the median is taken over every position.
    """
    per_contig: Dict[str, List[int]] = defaultdict(list)
    if not depth_path or not os.path.exists(depth_path) or os.path.getsize(depth_path) == 0:
        return {}
    with open(depth_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            try:
                depth = int(parts[2])
            except ValueError:
                continue
            per_contig[contig].append(depth)
    return {c: (len(vals), float(statistics.median(vals))) for c, vals in per_contig.items()}


def build_taxon_best_contig(
    krona_files: List[Tuple[str, str]],
    depth_by_sample: Dict[str, Dict[str, Tuple[int, float]]],
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
) -> Dict[Tuple[str, str, str], Tuple[int, float]]:
    """Map (sample, rank, taxid) -> (largest_contig_bp, its median depth).

    Each contig's leaf taxid is climbed to its family/genus/species ancestors;
    the contig contributes to every ancestor row at a rank of interest. The
    winner per (sample, rank, taxid) is the contig with the greatest length.
    """
    best: Dict[Tuple[str, str, str], Tuple[int, float]] = {}
    for sample, path in krona_files:
        depths = depth_by_sample.get(sample, {})
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 2:
                    continue
                contig, taxid = cols[0], cols[1]
                stat = depths.get(contig)
                if stat is None:
                    continue  # no depth data for this contig (unmapped / empty bam)
                length, med_depth = stat
                for tid in get_lineage(taxid, parent_map):
                    rank = rank_map.get(tid)
                    if rank not in RANKS_OF_INTEREST:
                        continue
                    key = (sample, rank, tid)
                    prev = best.get(key)
                    if prev is None or length > prev[0]:
                        best[key] = (length, med_depth)
    return best


def add_contig_stats(
    summary_df: pd.DataFrame,
    best: Dict[Tuple[str, str, str], Tuple[int, float]],
) -> pd.DataFrame:
    """Return *summary_df* with largest_contig_bp / largest_contig_median_depth."""
    required = {"sample", "rank", "taxid"}
    missing = required - set(summary_df.columns)
    if missing:
        raise ValueError(
            f"Summary is missing required column(s) {sorted(missing)}; "
            f"columns present: {list(summary_df.columns)}"
        )
    out = summary_df.copy()
    has_glen = "genome_length_bp" in out.columns
    lengths, covs, depths = [], [], []
    for _, row in out.iterrows():
        key = (str(row["sample"]), str(row["rank"]), str(row["taxid"]))
        stat = best.get(key)
        if stat is None:
            lengths.append(NA)
            covs.append(NA)
            depths.append(NA)
        else:
            lengths.append(stat[0])
            covs.append(_ref_coverage_pct(stat[0], row["genome_length_bp"]) if has_glen else NA)
            depths.append(stat[1])
    out["largest_contig_bp"] = lengths
    out["largest_contig_ref_coverage_pct"] = covs
    out["largest_contig_median_depth"] = depths
    return out


def run(
    summary_path: str,
    output_path: str,
    depth_files: List[Tuple[str, str]],
    krona_files: List[Tuple[str, str]],
    taxdump_dir: str,
) -> None:
    parent_map, rank_map, _ = load_taxdump(
        os.path.join(taxdump_dir, "nodes.dmp"), os.path.join(taxdump_dir, "names.dmp")
    )
    depth_by_sample = {sample: load_contig_depths(path) for sample, path in depth_files}
    best = build_taxon_best_contig(krona_files, depth_by_sample, parent_map, rank_map)

    if not os.path.exists(summary_path) or os.path.getsize(summary_path) == 0:
        open(output_path, "w").close()
        return

    summary_df = pd.read_csv(summary_path, sep="\t", dtype=str)
    out = add_contig_stats(summary_df, best)
    out.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    sys.stderr.write(f"[add_contig_stats_to_summary] wrote {len(out)} rows to {output_path}\n")


def _parse_pairs(items: Optional[List[str]]) -> List[Tuple[str, str]]:
    pairs = []
    for item in items or []:
        sample, sep, path = item.partition("=")
        if not sep:
            raise ValueError(f"expected SAMPLE=PATH, got: {item}")
        pairs.append((sample, path))
    return pairs


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Add largest-contig size + median depth per taxon.")
    p.add_argument("--summary", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--taxdump-dir", required=True)
    p.add_argument("--depth", action="append", required=True, metavar="SAMPLE=PATH")
    p.add_argument("--krona", action="append", required=True, metavar="SAMPLE=PATH")
    return p.parse_args()


def main() -> None:
    if "snakemake" in globals():
        samples = list(snakemake.params.samples)  # noqa: F821
        depth_files = list(zip(samples, list(snakemake.input.depth)))  # noqa: F821
        krona_files = list(zip(samples, list(snakemake.input.krona)))  # noqa: F821
        run(
            summary_path=str(snakemake.input.summary),  # noqa: F821
            output_path=str(snakemake.output.summary),  # noqa: F821
            depth_files=depth_files,
            krona_files=krona_files,
            taxdump_dir=str(snakemake.params.taxdump),  # noqa: F821
        )
    else:
        args = parse_args()
        run(
            summary_path=args.summary,
            output_path=args.output,
            depth_files=_parse_pairs(args.depth),
            krona_files=_parse_pairs(args.krona),
            taxdump_dir=args.taxdump_dir,
        )


if __name__ == "__main__":
    main()
