#!/usr/bin/env python3
"""Reduce a per-query DIAMOND-blastx-vs-nr table (21 columns: 12 BLAST + 9
lineage ranks, no header, grouped by query, up to 10 hits/query best-first) to
ONE representative hit per query plus an LCA consensus classification.

Ported from the validated REVISA nr_validation prototype. Per query:

1. Representative hit: the first hit with a non-"NA" species; else the best hit.
2. LCA consensus over ALL of the query's hits: climb species -> genus -> ... ->
   domain; the first rank where the most-common non-NA value reaches
   ``threshold`` (default 0.5 of the query's hits, NA counting against) wins.
3. Viral flag: the representative's phylum ends in ``viricota`` (ICTV convention).

Outputs a 23-column table (21 + consensus_rank + consensus_taxon) with a header,
and a viruses-only subset. This is the per-(sample,contig) NR verdict consumed by
harmonize_nr_summary.py.
"""

import argparse
import sys
from collections import Counter
from typing import List, Tuple

NA = "NA"
VIRAL_PHYLUM_SUFFIX = "viricota"

BLAST_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]
TAX_COLS = ["domain", "kingdom", "realm", "phylum", "class", "order", "family", "genus", "species"]
ALL_COLS = BLAST_COLS + TAX_COLS
N_COLS = len(ALL_COLS)

SPECIES_IDX = ALL_COLS.index("species")
PHYLUM_IDX = ALL_COLS.index("phylum")

RANKS_SPECIFIC_TO_GENERAL = [
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "realm",
    "kingdom",
    "domain",
]
RANK_IDX = {rank: ALL_COLS.index(rank) for rank in RANKS_SPECIFIC_TO_GENERAL}
OUT_HEADER = ALL_COLS + ["consensus_rank", "consensus_taxon"]


def choose_representative(hits: List[List[str]]) -> List[str]:
    """First hit with a non-NA species, else the best (first) hit."""
    for hit in hits:
        if hit[SPECIES_IDX] != NA:
            return hit
    return hits[0]


def is_viral(row: List[str]) -> bool:
    """True if the row's phylum is a viral phylum (name ends in 'viricota')."""
    phylum = row[PHYLUM_IDX]
    return phylum != NA and phylum.lower().endswith(VIRAL_PHYLUM_SUFFIX)


def lca_consensus(hits: List[List[str]], threshold: float) -> Tuple[str, str]:
    """Climb ranks from species upward; return (rank, taxon) of first consensus."""
    n = len(hits)
    for rank in RANKS_SPECIFIC_TO_GENERAL:
        idx = RANK_IDX[rank]
        counts = Counter(hit[idx] for hit in hits if hit[idx] != NA)
        if not counts:
            continue
        taxon, count = counts.most_common(1)[0]
        if count / n >= threshold:
            return rank, taxon
    return "unclassified", NA


def process_file(
    in_path: str,
    out_path: str,
    viral_path: str,
    threshold: float = 0.5,
) -> Tuple[int, int]:
    """Collapse the grouped input to one row/query; write the full and
    viruses-only tables (both with a header). Returns (queries, viral_queries)."""
    queries = 0
    viral = 0

    def _flush(group, out, vout):
        nonlocal queries, viral
        if not group:
            return
        queries += 1
        rep = choose_representative(group)
        rank, taxon = lca_consensus(group, threshold)
        record = "\t".join(rep + [rank, taxon]) + "\n"
        out.write(record)
        if is_viral(rep):
            viral += 1
            vout.write(record)

    header = "\t".join(OUT_HEADER) + "\n"
    skipped = 0
    with open(in_path) as inp, open(out_path, "w") as out, open(viral_path, "w") as vout:
        out.write(header)
        vout.write(header)
        current_q = None
        group: List[List[str]] = []
        for line in inp:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != N_COLS:
                skipped += 1
                continue
            if cols[0] != current_q:
                _flush(group, out, vout)
                current_q = cols[0]
                group = []
            group.append(cols)
        _flush(group, out, vout)
    if skipped:
        print(
            f"[filter_top_species_hits] WARNING: skipped {skipped} row(s) not exactly "
            f"{N_COLS} columns (upstream annotation malformed?); {queries} queries kept.",
            file=sys.stderr,
        )
    return queries, viral


def _derive_viral_path(output: str) -> str:
    if output.endswith(".tsv"):
        return output[: -len(".tsv")] + ".viruses_only.tsv"
    return output + ".viruses_only.tsv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pick representative hit + LCA consensus per query.")
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--viral-output", default=None)
    p.add_argument("--consensus-threshold", type=float, default=0.5)
    return p.parse_args()


def main() -> None:
    if "snakemake" in globals():
        in_path = str(snakemake.input[0])  # noqa: F821
        out_path = str(snakemake.output.table)  # noqa: F821
        viral_path = str(snakemake.output.viral)  # noqa: F821
        threshold = float(getattr(snakemake.params, "consensus_threshold", 0.5))  # noqa: F821
    else:
        args = parse_args()
        in_path = args.input
        out_path = args.output
        viral_path = args.viral_output or _derive_viral_path(args.output)
        threshold = args.consensus_threshold

    queries, viral = process_file(in_path, out_path, viral_path, threshold)
    print(
        f"[filter_top_species_hits] queries={queries} viral={viral} -> {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
