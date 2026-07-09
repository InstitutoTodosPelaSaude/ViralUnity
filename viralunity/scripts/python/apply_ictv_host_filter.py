#!/usr/bin/env python3
"""Taxonomic filter: keep only viruses whose lineage is a known vertebrate-
infecting clade (ICTV-derived allowlist); drop everything else (bacteriophages,
archaeal viruses, plant/fungal/algal/protist/invertebrate-only viruses).

This is a *taxonomic* filter and runs BEFORE the bleed / negative-control
filters, so it removes rows (rather than flagging them): the statistical filters
downstream then compute their thresholds on taxonomically-valid taxa only. A row
is kept iff its taxid or any ancestor is in the allowlist. Removed rows are
written to a ``*.dropped.tsv`` sidecar with a ``drop_reason`` column for audit.

The allowlist is a plain text file of NCBI taxids (one per line, ``#`` comments
and blank lines ignored), built from the ICTV Virus Metadata Resource by
``build_ictv_vertebrate_taxids.py``. It is user-editable and the filter is
off by default.
"""

import argparse
import os
import sys
from typing import Dict, Set, Tuple

import pandas as pd

try:
    from viralunity.scripts.python.taxonomy import get_lineage, load_taxdump
except ImportError:
    from taxonomy import get_lineage, load_taxdump

DROP_REASON = "not_vertebrate_virus"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Keep only vertebrate-infecting virus taxa (ICTV allowlist)."
    )
    parser.add_argument("--summary", required=True, help="Input taxa summary TSV")
    parser.add_argument("--output", required=True, help="Filtered summary TSV")
    parser.add_argument("--dropped", required=True, help="Dropped-rows audit TSV")
    parser.add_argument("--allowlist", required=True, help="Vertebrate-virus taxids file")
    parser.add_argument("--nodes", required=True, help="taxdump nodes.dmp")
    parser.add_argument("--names", default=None, help="taxdump names.dmp (optional)")
    return parser.parse_args()


def load_allowlist(path: str) -> Set[str]:
    """Read a taxid allowlist file (one taxid per line; ``#`` comments/blank
    lines ignored; anything after whitespace on a line is discarded)."""
    allow: Set[str] = set()
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            allow.add(line.split()[0])
    return allow


def lineage_allowed(taxid: str, allowset: Set[str], parent_map: Dict[str, str]) -> bool:
    """True iff ``taxid`` or any ancestor is in ``allowset``. Taxid 0 (unclassified)
    and taxids absent from the taxdump are never allowed."""
    taxid = str(taxid)
    if taxid in ("0", "", "nan") or taxid not in parent_map:
        return False
    return any(node in allowset for node in get_lineage(taxid, parent_map))


def filter_summary(
    in_path: str,
    out_path: str,
    dropped_path: str,
    allowset: Set[str],
    parent_map: Dict[str, str],
) -> Tuple[int, int]:
    """Split the summary into kept (vertebrate-virus lineage) and dropped rows.

    Returns ``(n_kept, n_dropped)``. An empty/missing input yields empty outputs.
    """
    if not os.path.exists(in_path) or os.path.getsize(in_path) == 0:
        open(out_path, "w").close()
        open(dropped_path, "w").close()
        return 0, 0

    df = pd.read_csv(in_path, sep="\t", dtype=str)
    if df.empty:
        df.to_csv(out_path, sep="\t", index=False)
        dropped = df.copy()
        dropped["drop_reason"] = pd.Series(dtype=str)
        dropped.to_csv(dropped_path, sep="\t", index=False)
        return 0, 0

    if "taxid" not in df.columns:
        raise ValueError(
            f"Input summary {in_path!r} has no 'taxid' column; cannot apply the "
            f"ICTV host filter. Columns present: {list(df.columns)}"
        )

    keep_mask = df["taxid"].apply(lambda t: lineage_allowed(t, allowset, parent_map))
    kept = df[keep_mask]
    dropped = df[~keep_mask].copy()
    dropped["drop_reason"] = DROP_REASON

    kept.to_csv(out_path, sep="\t", index=False)
    dropped.to_csv(dropped_path, sep="\t", index=False)
    return len(kept), len(dropped)


def run(
    summary: str,
    output: str,
    dropped: str,
    allowlist_file: str,
    nodes_dmp: str,
    names_dmp: str = None,
) -> Tuple[int, int]:
    allowset = load_allowlist(allowlist_file)
    parent_map, _rank, _name = load_taxdump(nodes_dmp, names_dmp)
    kept, n_dropped = filter_summary(summary, output, dropped, allowset, parent_map)
    print(
        f"[apply_ictv_host_filter] allowlist={len(allowset)} taxids "
        f"kept={kept} dropped={n_dropped} -> {output}",
        file=sys.stderr,
    )
    return kept, n_dropped


def main() -> None:
    if "snakemake" in globals():
        taxdump = str(snakemake.params.taxdump)  # noqa: F821
        run(
            summary=str(snakemake.input.summary),  # noqa: F821
            output=str(snakemake.output.summary),  # noqa: F821
            dropped=str(snakemake.output.dropped),  # noqa: F821
            allowlist_file=str(snakemake.params.allowlist),  # noqa: F821
            nodes_dmp=os.path.join(taxdump, "nodes.dmp"),
            names_dmp=os.path.join(taxdump, "names.dmp"),
        )
    else:
        args = parse_args()
        run(
            summary=args.summary,
            output=args.output,
            dropped=args.dropped,
            allowlist_file=args.allowlist,
            nodes_dmp=args.nodes,
            names_dmp=args.names,
        )


if __name__ == "__main__":
    main()
