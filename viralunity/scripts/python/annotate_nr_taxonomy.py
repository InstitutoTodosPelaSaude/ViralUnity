#!/usr/bin/env python3
"""Append NCBI lineage ranks to a DIAMOND-blastx-vs-nr outfmt-6 table.

Input rows are the 12 standard BLAST columns plus a trailing ``staxids`` field
(as emitted by ``diamond blastx --outfmt 6 ... staxids``). For each row we
resolve the subject taxid to its lineage and emit the 12 BLAST columns followed
by 9 rank-name columns:

    domain kingdom realm phylum class order family genus species

(NCBI ``superkingdom`` is reported as ``domain``; viruses populate ``realm``
rather than ``domain``.) Missing ranks are the literal string ``NA``. The 21-
column, header-less, query-grouped output is exactly what
``filter_top_species_hits.py`` consumes. This replaces the prototype's separate
accession->taxid step by using DIAMOND's own ``staxids``.
"""

import argparse
import os
import sys
from typing import Dict, List

try:
    from viralunity.scripts.python.taxonomy import get_lineage, load_taxdump
except ImportError:
    from taxonomy import get_lineage, load_taxdump

NA = "NA"
RANK_ORDER = [
    "domain",
    "kingdom",
    "realm",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]
# NCBI rank name -> our output rank (superkingdom is reported as domain).
_RANK_ALIAS = {"superkingdom": "domain"}


def resolve_lineage_ranks(
    taxid: str,
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
    name_map: Dict[str, str],
) -> List[str]:
    """Return the 9 rank names (RANK_ORDER) for ``taxid`` ('NA' where absent).

    ``taxid`` may be a ';'-separated list (DIAMOND staxids); the first is used.
    """
    taxid = str(taxid).split(";")[0].strip()
    result = {rank: NA for rank in RANK_ORDER}
    if not taxid or taxid in ("0", "nan") or taxid not in parent_map:
        return [NA] * len(RANK_ORDER)
    for node in get_lineage(taxid, parent_map):
        raw_rank = rank_map.get(node, "no rank")
        out_rank = _RANK_ALIAS.get(raw_rank, raw_rank)
        if out_rank in result and result[out_rank] == NA:
            result[out_rank] = name_map.get(node, NA) if name_map else NA
    return [result[rank] for rank in RANK_ORDER]


def annotate(
    in_path: str,
    out_path: str,
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
    name_map: Dict[str, str],
) -> int:
    """Write the 21-column annotated table. Returns rows written. Empty/missing
    input yields an empty output."""
    if not os.path.exists(in_path) or os.path.getsize(in_path) == 0:
        open(out_path, "w").close()
        return 0
    written = 0
    with open(in_path) as inp, open(out_path, "w") as out:
        for line in inp:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 13:
                continue
            blast = cols[:12]
            staxids = cols[12]
            ranks = resolve_lineage_ranks(staxids, parent_map, rank_map, name_map)
            out.write("\t".join(blast + ranks) + "\n")
            written += 1
    return written


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Append lineage ranks to diamond-vs-nr TSV.")
    p.add_argument("--input", required=True, help="diamond outfmt6 + staxids TSV")
    p.add_argument("--output", required=True, help="21-column annotated TSV")
    p.add_argument("--nodes", required=True, help="taxdump nodes.dmp")
    p.add_argument("--names", required=True, help="taxdump names.dmp")
    return p.parse_args()


def main() -> None:
    if "snakemake" in globals():
        taxdump = str(snakemake.params.taxdump)  # noqa: F821
        in_path = str(snakemake.input[0])  # noqa: F821
        out_path = str(snakemake.output[0])  # noqa: F821
        nodes = os.path.join(taxdump, "nodes.dmp")
        names = os.path.join(taxdump, "names.dmp")
    else:
        args = parse_args()
        in_path, out_path, nodes, names = args.input, args.output, args.nodes, args.names

    parent_map, rank_map, name_map = load_taxdump(nodes, names)
    n = annotate(in_path, out_path, parent_map, rank_map, name_map)
    print(f"[annotate_nr_taxonomy] annotated {n} rows -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
