#!/usr/bin/env python3
"""
Filter a per-sample krona_input.tsv (query<TAB>leaf_taxid) using the per-sample
'Pass' set extracted from the fully-filtered taxa-summary table (the last file in
the cumulative chain, e.g. `*_taxa_summary_RPM.bleed[.neg][.ictv].tsv`).

Matching is *lineage-aware*: a krona_input row is kept iff any ancestor of its
leaf taxid whose rank is family/genus/species appears in the Pass set built for
the row's (sample, tool, mode). This is the inverse of the lineage walk that
`summarize_krona_taxa.py` uses to build the per-rank summaries, so by
construction every krona_input row that contributed to a passing summary row
is kept and nothing else is.
"""

import argparse
import os
import sys
from typing import Dict, Optional, Set, Tuple

import pandas as pd

try:
    from viralunity.scripts.python.taxonomy import (  # noqa: F401 (re-exported)
        RANKS_OF_INTEREST,
        get_lineage,
    )
    from viralunity.scripts.python.taxonomy import load_taxdump as _load_taxdump_full
except ImportError:
    # Running inside Snakemake's `script:` directive (per-rule conda env without
    # the viralunity package installed); the script's own dir is on sys.path.
    from taxonomy import (  # noqa: F401 (re-exported)
        RANKS_OF_INTEREST,
        get_lineage,
    )
    from taxonomy import load_taxdump as _load_taxdump_full


def load_taxdump(
    nodes_dmp: str, names_dmp: Optional[str] = None
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Read nodes.dmp and return (parent_map, rank_map).

    The names_dmp argument is accepted for API parity; the lineage walk itself
    does not need taxon names. Delegates to
    :func:`viralunity.scripts.python.taxonomy.load_taxdump`.
    """
    parent, rank, _name = _load_taxdump_full(nodes_dmp, names_dmp)
    return parent, rank


def build_pass_taxids(
    summary_path: str,
    sample: str,
    tool: str,
    mode: str,
) -> Set[str]:
    """
    Extract the per-sample pass-set of taxids from the fully-filtered summary file.

    - Always requires bleed_pass == True.
    - If a `neg_pass` column is present, only drop explicit False; NA (taxon
      absent from negative controls) is treated as a pass so the filter stays
      conservative.
    """
    df = pd.read_csv(summary_path, sep="\t", dtype={"taxid": str})

    if "taxid" not in df.columns:
        raise ValueError(f"Summary {summary_path} is missing required 'taxid' column.")
    if "bleed_pass" not in df.columns:
        raise ValueError(
            f"Summary {summary_path} is missing required 'bleed_pass' column "
            "(run apply_max_rpm_bleed_filter.py first)."
        )

    mask = (
        (df.get("sample") == sample)
        & (df.get("tool") == tool)
        & (df.get("mode") == mode)
        & df["bleed_pass"].astype(bool)
    )

    if "neg_pass" in df.columns:
        # True → keep, False → drop, NA → keep (conservative).
        # neg_pass is set by add_negative_control_enrichment.py (fold-enrichment /
        # log10-ratio / z-score gate) — it replaced the old Poisson-based filter.
        # NA means zero negative controls were configured, so no filtering is applied.
        mask &= df["neg_pass"].where(df["neg_pass"].notna(), True).astype(bool)

    return set(df.loc[mask, "taxid"].astype(str))


def filter_krona_input(
    krona_input_path: str,
    output_path: str,
    pass_taxids: Set[str],
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
) -> Tuple[int, int]:
    """
    Stream `krona_input_path` and write rows whose leaf taxid has any
    family/genus/species ancestor in `pass_taxids`.

    Returns (kept, total) row counts.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    keep_cache: Dict[str, bool] = {}
    total = 0
    kept = 0

    with open(krona_input_path) as fin, open(output_path, "w") as fout:
        if not pass_taxids:
            return 0, sum(1 for line in fin if line.strip())

        for line in fin:
            if not line.strip():
                continue
            total += 1

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            tid = parts[1]

            cached = keep_cache.get(tid)
            if cached is None:
                cached = _lineage_passes(tid, pass_taxids, parent_map, rank_map)
                keep_cache[tid] = cached

            if cached:
                fout.write(line)
                kept += 1

    return kept, total


def _lineage_passes(
    taxid: str,
    pass_taxids: Set[str],
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
) -> bool:
    if taxid == "0" or taxid not in parent_map:
        return False
    for anc in get_lineage(taxid, parent_map):
        if rank_map.get(anc) in RANKS_OF_INTEREST and anc in pass_taxids:
            return True
    return False


def run(
    summary: str,
    krona_input: str,
    output: str,
    sample: str,
    tool: str,
    mode: str,
    nodes_dmp: str,
    names_dmp: Optional[str] = None,
) -> None:
    pass_taxids = build_pass_taxids(summary, sample=sample, tool=tool, mode=mode)
    parent_map, rank_map = load_taxdump(nodes_dmp, names_dmp)
    kept, total = filter_krona_input(
        krona_input_path=krona_input,
        output_path=output,
        pass_taxids=pass_taxids,
        parent_map=parent_map,
        rank_map=rank_map,
    )
    sys.stderr.write(
        f"[filter_krona_by_pass_taxids] sample={sample} tool={tool} mode={mode} "
        f"pass_taxids={len(pass_taxids)} kept={kept}/{total}\n"
    )


def run_cli() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Filter a krona_input.tsv by the per-sample Pass set extracted from "
            "the fully-filtered taxa-summary table, using lineage-aware matching."
        )
    )
    ap.add_argument("--summary", required=True, help="Path to the fully-filtered summary TSV.")
    ap.add_argument("--krona-input", required=True, help="Per-sample krona_input.tsv.")
    ap.add_argument("--out", required=True, help="Output filtered krona_input.tsv.")
    ap.add_argument(
        "--sample", required=True, help="Sample id (must match summary's 'sample' column)."
    )
    ap.add_argument("--tool", required=True, choices=["kraken2", "diamond"])
    ap.add_argument("--mode", required=True, choices=["reads", "contigs"])
    ap.add_argument(
        "--taxdump-dir", required=True, help="Directory with nodes.dmp [and names.dmp]."
    )
    args = ap.parse_args()

    run(
        summary=args.summary,
        krona_input=args.krona_input,
        output=args.out,
        sample=args.sample,
        tool=args.tool,
        mode=args.mode,
        nodes_dmp=os.path.join(args.taxdump_dir, "nodes.dmp"),
        names_dmp=os.path.join(args.taxdump_dir, "names.dmp"),
    )


def run_snakemake() -> None:
    run(
        summary=str(snakemake.input.summary),
        krona_input=str(snakemake.input.krona_input),
        output=str(snakemake.output[0]),
        sample=str(snakemake.params.sample),
        tool=str(snakemake.params.tool),
        mode=str(snakemake.params.mode),
        nodes_dmp=os.path.join(str(snakemake.params.taxdump), "nodes.dmp"),
        names_dmp=os.path.join(str(snakemake.params.taxdump), "names.dmp"),
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
