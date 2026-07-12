#!/usr/bin/env python3
"""Split a combined taxa summary into per-rank tables (the user-facing output).

The filter chain computes on one combined per-track table (kept internally under
``<track>/summaries/full/``). This terminal step splits that table into the browsable
deliverable: one file per taxonomic rank under ``<track>/summaries/<rank>/``, propagating
higher-rank *names* down as columns so each table is self-contained:

  * species table gains ``family`` and ``genus`` columns
  * genus table gains a ``family`` column
  * family table is unchanged

It also guarantees a ``final_species`` column exists on every track (tracks
without NR validation — reads tracks, NR-off contig tracks — get
``final_species = name``), so the confirmed-taxonomy column is universal.
"""

import argparse
import os
import sys
from typing import Dict, List, Tuple

import pandas as pd

try:
    from viralunity.scripts.python.taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump
except ImportError:
    from taxonomy import RANKS_OF_INTEREST, get_lineage, load_taxdump

NA = "NA"
# Higher-rank name columns propagated into each rank's table.
HIGHER_RANKS: Dict[str, Tuple[str, ...]] = {
    "family": (),
    "genus": ("family",),
    "species": ("family", "genus"),
}


def _insert_after(cols: List[str], anchor: str, new_cols: List[str]) -> List[str]:
    """Return column order with *new_cols* inserted right after *anchor*."""
    cols = [c for c in cols if c not in new_cols]
    if anchor in cols:
        at = cols.index(anchor) + 1
    else:
        at = len(cols)
    return cols[:at] + new_cols + cols[at:]


def ensure_final_species(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``final_species`` = coalesce(nr_correct_species, name) if not present."""
    if "final_species" in df.columns or "name" not in df.columns:
        return df
    df = df.copy()
    has_corr = "nr_correct_species" in df.columns

    def _coalesce(row):
        if has_corr:
            cs = row["nr_correct_species"]
            if cs is not None and str(cs) not in ("", NA):
                return cs
        return row["name"]

    df["final_species"] = df.apply(_coalesce, axis=1)
    anchor = "nr_correct_species" if has_corr else "name"
    return df[_insert_after(list(df.columns), anchor, ["final_species"])]


def add_higher_rank_names(
    df: pd.DataFrame,
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
    name_map: Dict[str, str],
) -> pd.DataFrame:
    """Add ``family``/``genus`` ancestor-name columns (positioned after ``name``).

    Each row's taxid is climbed once; the family and genus ancestor names are
    recorded. Rows are filtered per-rank later, so both columns are computed for
    every row and the irrelevant ones dropped at write time.
    """
    df = df.copy()
    cache: Dict[str, Tuple[str, str]] = {}

    def _ancestors(taxid: str) -> Tuple[str, str]:
        if taxid in cache:
            return cache[taxid]
        fam = gen = NA
        if taxid and taxid != "0":
            for tid in get_lineage(taxid, parent_map):
                r = rank_map.get(tid)
                if r == "family" and fam == NA:
                    fam = name_map.get(tid, NA)
                elif r == "genus" and gen == NA:
                    gen = name_map.get(tid, NA)
        cache[taxid] = (fam, gen)
        return fam, gen

    fams, gens = [], []
    for taxid in df["taxid"].astype(str):
        fam, gen = _ancestors(taxid)
        fams.append(fam)
        gens.append(gen)
    df["family"] = fams
    df["genus"] = gens
    return df[_insert_after(list(df.columns), "name", ["family", "genus"])]


def split_by_rank(
    df: pd.DataFrame,
    out_paths: Dict[str, str],
    parent_map: Dict[str, str],
    rank_map: Dict[str, str],
    name_map: Dict[str, str],
) -> None:
    """Write one file per rank in *out_paths*, with higher-rank columns propagated."""
    if "rank" not in df.columns or "taxid" not in df.columns:
        raise ValueError(
            f"Summary must contain 'rank' and 'taxid'; columns present: {list(df.columns)}"
        )
    df = ensure_final_species(df)
    df = add_higher_rank_names(df, parent_map, rank_map, name_map)

    for rank, path in out_paths.items():
        os.makedirs(os.path.dirname(path), exist_ok=True)
        sub = df[df["rank"] == rank].copy()
        # Keep only the higher-rank name columns relevant to this rank.
        drop = [c for c in ("family", "genus") if c not in HIGHER_RANKS.get(rank, ())]
        sub = sub.drop(columns=[c for c in drop if c in sub.columns])
        sub.to_csv(path, sep="\t", index=False, na_rep="NA")


def run(summary_path: str, out_paths: Dict[str, str], taxdump_dir: str) -> None:
    parent_map, rank_map, name_map = load_taxdump(
        os.path.join(taxdump_dir, "nodes.dmp"), os.path.join(taxdump_dir, "names.dmp")
    )
    if not os.path.exists(summary_path) or os.path.getsize(summary_path) == 0:
        for path in out_paths.values():
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").close()
        return
    df = pd.read_csv(summary_path, sep="\t", dtype=str)
    split_by_rank(df, out_paths, parent_map, rank_map, name_map)
    sys.stderr.write(f"[split_summary_by_rank] wrote {len(out_paths)} per-rank tables\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Split a combined taxa summary into per-rank tables.")
    p.add_argument("--summary", required=True)
    p.add_argument("--taxdump-dir", required=True)
    for rank in RANKS_OF_INTEREST:
        p.add_argument(f"--{rank}", required=True, help=f"output path for the {rank} table")
    return p.parse_args()


def main() -> None:
    if "snakemake" in globals():
        out_paths = {
            rank: str(getattr(snakemake.output, rank)) for rank in RANKS_OF_INTEREST
        }  # noqa: F821
        run(
            summary_path=str(snakemake.input.summary),  # noqa: F821
            out_paths=out_paths,
            taxdump_dir=str(snakemake.params.taxdump),  # noqa: F821
        )
    else:
        args = parse_args()
        out_paths = {rank: getattr(args, rank) for rank in RANKS_OF_INTEREST}
        run(summary_path=args.summary, out_paths=out_paths, taxdump_dir=args.taxdump_dir)


if __name__ == "__main__":
    main()
