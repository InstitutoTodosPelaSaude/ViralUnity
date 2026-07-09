#!/usr/bin/env python3

import os
import sys
from collections import defaultdict

try:
    from viralunity.scripts.python.taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump
except ImportError:
    # Running inside Snakemake's `script:` directive (per-rule conda env without
    # the viralunity package installed); the script's own dir is on sys.path.
    from taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump


def load_diamond_reads(diamond_tax_file):
    reads = {}
    with open(diamond_tax_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            contig_id = parts[0]
            try:
                mapped_reads = int(parts[-1])
            except ValueError:
                continue  # skip header or malformed lines
            reads[contig_id] = mapped_reads
    return reads


def summarize_krona(krona_file, parent_map, rank_map, diamond_reads=None):
    contig_counts = defaultdict(int)
    read_counts = defaultdict(int)
    totals_per_rank = defaultdict(int)

    with open(krona_file) as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                sys.stderr.write(f"WARNING: skipping malformed krona line: {line.rstrip()}\n")
                continue
            contig_id, taxid = parts[0], parts[1]
            lineage = get_lineage(taxid, parent_map)
            seen_ranks = set()

            for tid in lineage:
                r = rank_map.get(tid)
                if r in RANKS_OF_INTEREST and r not in seen_ranks:
                    key = (r, tid)
                    contig_counts[key] += 1
                    totals_per_rank[r] += 1

                    if diamond_reads is not None:
                        read_counts[key] += diamond_reads.get(contig_id, 0)

                    seen_ranks.add(r)

    return contig_counts, read_counts, totals_per_rank


def run(
    krona,
    diamond_tax,
    taxdump_nodes,
    taxdump_names,
    sample,
    classifier,
    unit,
    output,
):
    parent_map, rank_map, name_map = load_taxdump(taxdump_nodes, taxdump_names)

    diamond_reads = None
    if diamond_tax:
        diamond_reads = (
            load_diamond_reads(diamond_tax)
            if os.path.exists(diamond_tax) and os.path.getsize(diamond_tax) > 0
            else {}
        )

    contig_counts, read_counts, totals_per_rank = summarize_krona(
        krona, parent_map, rank_map, diamond_reads
    )

    header = [
        "sample",
        "tool",
        "mode",
        "rank",
        "taxid",
        "name",
        "count",
        "percent",
        "source",
    ]

    if diamond_reads is not None:
        header.append("mapped_reads")

    with open(output, "w") as out:
        out.write("\t".join(header) + "\n")

        for rank in RANKS_OF_INTEREST:
            for (r, taxid), count in sorted(contig_counts.items()):
                if r != rank:
                    continue

                total = totals_per_rank[r]
                percent = (count / total * 100) if total else 0.0
                taxname = name_map.get(taxid, "NA")

                fields = [
                    sample,
                    classifier,
                    unit,
                    r,
                    taxid,
                    taxname,
                    str(count),
                    f"{percent:.4f}",
                    krona,
                ]

                if diamond_reads is not None:
                    fields.append(str(read_counts.get((r, taxid), 0)))

                out.write("\t".join(fields) + "\n")


# ---- Snakemake entrypoint ----
if "snakemake" in globals():
    run(
        krona=snakemake.input.krona,
        diamond_tax=snakemake.input.get("annotated", None),
        taxdump_nodes=os.path.join(snakemake.params.taxdump, "nodes.dmp"),
        taxdump_names=os.path.join(snakemake.params.taxdump, "names.dmp"),
        sample=snakemake.params.sample,
        classifier=snakemake.params.tool,
        unit=snakemake.params.mode,
        output=snakemake.output[0],
    )
elif __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--krona", required=True)
    parser.add_argument("--diamond-tax")
    parser.add_argument("--taxdump-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--classifier", required=True)
    parser.add_argument("--unit", default="contigs")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    run(
        krona=args.krona,
        diamond_tax=args.diamond_tax,
        taxdump_nodes=os.path.join(args.taxdump_dir, "nodes.dmp"),
        taxdump_names=os.path.join(args.taxdump_dir, "names.dmp"),
        sample=args.sample,
        classifier=args.classifier,
        unit=args.unit,
        output=args.output,
    )
