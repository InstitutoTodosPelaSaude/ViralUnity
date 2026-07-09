#!/usr/bin/env python3
"""Taxonomic filter (kraken2 tracks): drop taxa supported by too few DISTINCT
minimizers, KrakenUniq-style.

Kraken2 run with ``--report-minimizer-data`` reports, per taxon, the total
minimizer count and the number of *distinct* minimizers. A high total with few
distinct minimizers is the classic signature of a spurious hit (a handful of
low-complexity/conserved k-mers matched many times), whereas a genuine taxon
accumulates many distinct minimizers across its genome. We drop a taxon when its
distinct-minimizer count is below ``min_distinct`` and/or its duplication ratio
(total / distinct) exceeds ``max_dup``.

Like the other taxonomic filters this runs BEFORE bleed/negative-control and
removes rows (writing a ``*.dropped.tsv`` audit sidecar). Rows for which no
report data is found are kept (the filter cannot assess them). Applies to
kraken2 tracks only; off by default.
"""

import argparse
import os
import sys
from typing import Dict, Optional, Tuple

import pandas as pd

DROP_REASON = "low_distinct_minimizers"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Drop kraken2 taxa with too few distinct minimizers."
    )
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dropped", required=True)
    parser.add_argument("--report", action="append", required=True, metavar="SAMPLE=PATH")
    parser.add_argument("--min-distinct", type=int, default=None)
    parser.add_argument("--max-duplication", type=float, default=None)
    return parser.parse_args()


def parse_report_minimizers(path: str) -> Dict[str, Tuple[int, int]]:
    """Parse a kraken2 ``--report-minimizer-data`` report.

    Returns ``{taxid: (n_minimizers, n_distinct_minimizers)}``. The
    minimizer-augmented report has 8 tab-separated columns:
    pct, clade_reads, taxon_reads, n_minimizers, n_distinct, rank, taxid, name.
    Lines without the two minimizer columns (plain 6-column reports) are skipped.
    """
    counts: Dict[str, Tuple[int, int]] = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return counts
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            try:
                n_min = int(parts[3])
                n_distinct = int(parts[4])
            except ValueError:
                continue
            taxid = parts[6].strip()
            counts[taxid] = (n_min, n_distinct)
    return counts


def passes_minimizers(
    n_distinct: int,
    n_total: int,
    min_distinct: Optional[int],
    max_dup: Optional[float],
) -> bool:
    """True iff the taxon clears the (optional) distinct-count and duplication
    criteria. A None/non-positive threshold disables that criterion."""
    if min_distinct and n_distinct < min_distinct:
        return False
    if max_dup and max_dup > 0:
        dup = (n_total / n_distinct) if n_distinct > 0 else float("inf")
        if dup > max_dup:
            return False
    return True


def filter_summary(
    in_path: str,
    out_path: str,
    dropped_path: str,
    counts_by_sample: Dict[str, Dict[str, Tuple[int, int]]],
    min_distinct: Optional[int],
    max_dup: Optional[float],
) -> Tuple[int, int]:
    """Drop rows whose (sample, taxid) minimizer support fails the thresholds.

    Rows with no matching report entry are kept. Returns ``(n_kept, n_dropped)``.
    """
    if not os.path.exists(in_path) or os.path.getsize(in_path) == 0:
        open(out_path, "w").close()
        open(dropped_path, "w").close()
        return 0, 0

    df = pd.read_csv(in_path, sep="\t", dtype=str)
    if df.empty:
        df.to_csv(out_path, sep="\t", index=False)
        df.assign(drop_reason=pd.Series(dtype=str)).to_csv(dropped_path, sep="\t", index=False)
        return 0, 0

    def _keep(row) -> bool:
        sample_counts = counts_by_sample.get(str(row["sample"]))
        if not sample_counts:
            return True
        entry = sample_counts.get(str(row["taxid"]))
        if entry is None:
            return True
        n_total, n_distinct = entry
        return passes_minimizers(n_distinct, n_total, min_distinct, max_dup)

    keep_mask = df.apply(_keep, axis=1)
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
    reports,
    samples,
    min_distinct: Optional[int],
    max_dup: Optional[float],
) -> Tuple[int, int]:
    counts_by_sample = {
        sample: parse_report_minimizers(report) for sample, report in zip(samples, reports)
    }
    kept, n_dropped = filter_summary(
        summary, output, dropped, counts_by_sample, min_distinct, max_dup
    )
    print(
        f"[apply_minimizer_filter] min_distinct={min_distinct} max_dup={max_dup} "
        f"kept={kept} dropped={n_dropped} -> {output}",
        file=sys.stderr,
    )
    return kept, n_dropped


def _parse_pairs(items):
    samples, paths = [], []
    for item in items:
        sample, sep, path = item.partition("=")
        if not sep:
            raise ValueError(f"--report must be SAMPLE=PATH, got: {item}")
        samples.append(sample)
        paths.append(path)
    return samples, paths


def main() -> None:
    if "snakemake" in globals():
        samples = list(snakemake.params.samples)  # noqa: F821
        reports = list(snakemake.input.reports)  # noqa: F821
        md = getattr(snakemake.params, "min_distinct", None)  # noqa: F821
        mx = getattr(snakemake.params, "max_duplication", None)  # noqa: F821
        run(
            summary=str(snakemake.input.summary),  # noqa: F821
            output=str(snakemake.output.summary),  # noqa: F821
            dropped=str(snakemake.output.dropped),  # noqa: F821
            reports=reports,
            samples=samples,
            min_distinct=int(md) if md else None,
            max_dup=float(mx) if mx else None,
        )
    else:
        args = parse_args()
        samples, reports = _parse_pairs(args.report)
        run(
            summary=args.summary,
            output=args.output,
            dropped=args.dropped,
            reports=reports,
            samples=samples,
            min_distinct=args.min_distinct,
            max_dup=args.max_duplication,
        )


if __name__ == "__main__":
    main()
