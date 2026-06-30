#!/usr/bin/env python3
"""Build a per-taxon representative genome-length table from a RefSeq viral FASTA index.

For each sequence in a samtools-faidx .fai file the taxid is looked up via a
genome2taxid mapping (accession<TAB>taxid, no header), then the length is walked
up the lineage and accumulated at every family/genus/species ancestor.  The final
output records the **median** length across all RefSeq sequences under each node,
which is used as the representative genome length for RPKM computation.

Output TSV columns:
  rank, taxid, name, genome_length_bp (median), n_genomes
"""

import argparse
import os
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from viralunity.scripts.python.taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump


def _parse_fai(fai_path: str) -> Dict[str, int]:
    """Return {accession_key: length} from a samtools faidx .fai file.

    The accession key is the first whitespace-delimited token of the sequence
    name, matching the format used in genome2taxid.tsv.
    """
    lengths: Dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            # name column may be "NC_001461.1 Dengue virus 1, ..."
            accession = parts[0].split()[0]
            try:
                lengths[accession] = int(parts[1])
            except ValueError:
                sys.stderr.write(f"WARNING: skipping malformed .fai line: {line.rstrip()}\n")
    return lengths


def _parse_genome2taxid(path: str) -> Dict[str, str]:
    """Return {accession: taxid} from a two-column no-header TSV."""
    mapping: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                sys.stderr.write(f"WARNING: skipping malformed genome2taxid line: {line.rstrip()}\n")
                continue
            accession, taxid = parts[0].strip(), parts[1].strip()
            mapping[accession] = taxid
    return mapping


def build_genome_length_table(
    fai_path: str,
    genome2taxid_path: str,
    nodes_dmp: str,
    names_dmp: Optional[str] = None,
) -> List[Tuple[str, str, str, int, int]]:
    """Compute median genome lengths per (rank, taxid) node.

    Parameters
    ----------
    fai_path:
        Path to the viral genomes ``.fai`` (samtools faidx output).
    genome2taxid_path:
        Path to the ``genome2taxid.tsv`` (accession<TAB>taxid, no header).
    nodes_dmp:
        Path to NCBI taxdump ``nodes.dmp``.
    names_dmp:
        Optional path to NCBI taxdump ``names.dmp``.  Used only for the name
        column in the output; the length roll-up itself uses only nodes.

    Returns
    -------
    List of ``(rank, taxid, name, median_length_bp, n_genomes)`` tuples, one
    per (rank, taxid) node that has at least one associated genome.
    """
    parent_map, rank_map, name_map = load_taxdump(nodes_dmp, names_dmp)
    fai_lengths = _parse_fai(fai_path)
    g2t = _parse_genome2taxid(genome2taxid_path)

    # Accumulate lengths at each lineage node whose rank is in RANKS_OF_INTEREST.
    node_lengths: Dict[Tuple[str, str], List[int]] = defaultdict(list)

    n_no_taxid = 0
    n_no_lineage = 0
    for accession, length in fai_lengths.items():
        taxid = g2t.get(accession)
        if taxid is None:
            n_no_taxid += 1
            continue

        lineage = get_lineage(taxid, parent_map)
        found_rank = False
        for anc in lineage:
            rk = rank_map.get(anc)
            if rk in RANKS_OF_INTEREST:
                node_lengths[(rk, anc)].append(length)
                found_rank = True

        if not found_rank:
            n_no_lineage += 1

    if n_no_taxid:
        sys.stderr.write(
            f"[build_genome_length_table] {n_no_taxid} accessions had no taxid in genome2taxid.\n"
        )
    if n_no_lineage:
        sys.stderr.write(
            f"[build_genome_length_table] {n_no_lineage} accessions had no family/genus/species ancestor.\n"
        )

    rows: List[Tuple[str, str, str, int, int]] = []
    for (rk, taxid), lengths in node_lengths.items():
        median_len = int(statistics.median(lengths))
        name = name_map.get(taxid, "NA")
        rows.append((rk, taxid, name, median_len, len(lengths)))

    return rows


def write_table(rows: List[Tuple[str, str, str, int, int]], out_path: str) -> None:
    """Write the table to a TSV file."""
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("rank\ttaxid\tname\tgenome_length_bp\tn_genomes\n")
        for rank, taxid, name, length, n in rows:
            fh.write(f"{rank}\t{taxid}\t{name}\t{length}\t{n}\n")
    sys.stderr.write(
        f"[build_genome_length_table] wrote {len(rows)} rows to {out_path}\n"
    )


def run_cli() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Build per-taxon median genome-length table from a viral RefSeq FASTA "
            "index for RPKM normalisation."
        )
    )
    ap.add_argument("--fai", required=True, help="Path to samtools faidx .fai of the viral FASTA.")
    ap.add_argument(
        "--genome2taxid",
        required=True,
        help="Two-column no-header TSV: accession<TAB>taxid.",
    )
    ap.add_argument(
        "--taxdump-dir",
        required=True,
        help="Directory containing nodes.dmp (and optionally names.dmp).",
    )
    ap.add_argument("--out", required=True, help="Output TSV path.")
    args = ap.parse_args()

    nodes_dmp = os.path.join(args.taxdump_dir, "nodes.dmp")
    names_dmp = os.path.join(args.taxdump_dir, "names.dmp")
    names_dmp = names_dmp if os.path.exists(names_dmp) else None

    rows = build_genome_length_table(
        fai_path=args.fai,
        genome2taxid_path=args.genome2taxid,
        nodes_dmp=nodes_dmp,
        names_dmp=names_dmp,
    )
    write_table(rows, args.out)


def run_snakemake() -> None:
    rows = build_genome_length_table(
        fai_path=str(snakemake.input.fai),
        genome2taxid_path=str(snakemake.input.genome2taxid),
        nodes_dmp=os.path.join(str(snakemake.params.taxdump), "nodes.dmp"),
        names_dmp=os.path.join(str(snakemake.params.taxdump), "names.dmp"),
    )
    write_table(rows, str(snakemake.output[0]))


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
