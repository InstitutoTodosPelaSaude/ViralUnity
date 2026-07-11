#!/usr/bin/env python3
import argparse
from typing import List, Optional

import pandas as pd


def infer_group_cols(df: pd.DataFrame, extra_group_cols: Optional[List[str]] = None) -> List[str]:
    """
    Decide what defines a 'taxon' for max-RPM purposes.

    Always include rank+taxid if present.
    Optionally include tool/mode if they exist in the table (keeps categories separate).
    Optionally include user-provided extra columns (e.g. database, classifier, etc.).
    """
    group_cols = []

    # Keep tool/mode separation if present (matches your earlier design: apply per category)
    for c in ["tool", "mode"]:
        if c in df.columns:
            group_cols.append(c)

    # Rank and taxid identify the taxon at a given level
    if "rank" in df.columns:
        group_cols.append("rank")
    if "taxid" in df.columns:
        group_cols.append("taxid")
    else:
        raise ValueError("Input must contain 'taxid' column to identify taxa.")

    if extra_group_cols:
        for c in extra_group_cols:
            if c not in df.columns:
                raise ValueError(f"Requested group column '{c}' not found in input columns.")
            group_cols.append(c)

    return group_cols


def apply_bleed_filter(
    df: pd.DataFrame,
    rpm_col: str = "rpm",
    rpkm_col: str = "rpkm",
    fraction: float = 0.005,
    rpm_floor: float = 1.0,
    rpkm_floor: float = 0.1,
    metric: str = "auto",
    group_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Add max-metric bleed-through filter columns.

    The comparison metric is chosen per taxon group: ``rpkm`` when the ``rpkm``
    column is present and the group has a non-NA value, else ``rpm`` (this is
    the ``metric="auto"`` default; ``"rpm"``/``"rpkm"`` force one). Because the
    bleed test is a WITHIN-taxon ratio (``value >= fraction * max`` across
    samples of the *same* taxid), and genome length is constant within a taxon,
    switching rpm→rpkm rescales every value in a group by the same constant and
    leaves the pass/fail decision unchanged. What differs is the *floor* gate:
    ``rpkm`` values sit on a different scale, so a metric-appropriate floor
    (``rpkm_floor`` vs ``rpm_floor``) decides whether the filter is applied.

    Columns added:
    - bleed_metric: "rpkm" or "rpm" — the metric used for that taxon group
    - bleed_max: maximum metric value for the group across samples
    - bleed_threshold: fraction * bleed_max (only when bleed_max >= floor)
    - bleed_applied: whether the filter was applied for that taxon
    - bleed_pass: True if value >= bleed_threshold OR the filter was not applied
    """
    if "sample" not in df.columns:
        raise ValueError("Input must contain 'sample' column.")
    if rpm_col not in df.columns:
        raise ValueError(
            f"Input missing RPM column '{rpm_col}'. Did you run add_RPM_to_summary first?"
        )
    if metric not in ("auto", "rpm", "rpkm"):
        raise ValueError(f"metric must be 'auto', 'rpm', or 'rpkm'; got {metric!r}.")

    out = df.copy()
    if group_cols is None:
        group_cols = infer_group_cols(out)

    has_rpkm = rpkm_col in out.columns
    if metric == "rpkm":
        if not has_rpkm:
            raise ValueError(f"metric='rpkm' but no '{rpkm_col}' column present.")
        out["bleed_metric"] = "rpkm"
    elif metric == "rpm":
        out["bleed_metric"] = "rpm"
    elif has_rpkm:
        # auto: rpkm when the group has any non-NA rpkm, else rpm.
        out["bleed_metric"] = out.groupby(group_cols, dropna=False)[rpkm_col].transform(
            lambda s: "rpkm" if s.notna().any() else "rpm"
        )
    else:
        out["bleed_metric"] = "rpm"

    def _metric_value(row):
        if row["bleed_metric"] == "rpkm" and has_rpkm:
            val = row[rpkm_col]
            if pd.notna(val):
                return float(val)
        return float(row[rpm_col])

    out["_bleed_val"] = out.apply(_metric_value, axis=1)
    out["_bleed_floor"] = out["bleed_metric"].map(
        {"rpkm": float(rpkm_floor), "rpm": float(rpm_floor)}
    )

    # Max metric value per taxon group across samples.
    out["bleed_max"] = out.groupby(group_cols, dropna=False)["_bleed_val"].transform("max")

    out["bleed_applied"] = out["bleed_max"] >= out["_bleed_floor"]

    out["bleed_threshold"] = 0.0
    mask = out["bleed_applied"]
    out.loc[mask, "bleed_threshold"] = out.loc[mask, "bleed_max"].astype(float) * float(fraction)

    # Pass if above threshold, OR if not applied (bleed_max below floor).
    out["bleed_pass"] = True
    out.loc[mask, "bleed_pass"] = out.loc[mask, "_bleed_val"].astype(float) >= out.loc[
        mask, "bleed_threshold"
    ].astype(float)

    out = out.drop(columns=["_bleed_val", "_bleed_floor"], errors="ignore")
    return out


def run_cli():
    ap = argparse.ArgumentParser(
        description="Apply max-RPM bleed-through filter to an RPM-augmented taxa summary TSV."
    )
    ap.add_argument(
        "--in",
        dest="inp",
        required=True,
        help="Input TSV (must include sample,taxid,rpm).",
    )
    ap.add_argument("--out", required=True, help="Output TSV with bleed filter columns added.")
    ap.add_argument("--rpm-col", default="rpm", help="Name of RPM column (default: rpm).")
    ap.add_argument("--rpkm-col", default="rpkm", help="Name of RPKM column (default: rpkm).")
    ap.add_argument(
        "--metric",
        default="auto",
        choices=["auto", "rpm", "rpkm"],
        help="Metric to bleed-filter on (default: auto = rpkm when available per group, else rpm).",
    )
    ap.add_argument(
        "--fraction",
        type=float,
        default=0.005,
        help="Threshold fraction of the group's max metric (default: 0.005).",
    )
    ap.add_argument(
        "--rpm-floor",
        type=float,
        default=1.0,
        help="If a rpm-metric group's max < this floor, do NOT apply the filter (default: 1.0).",
    )
    ap.add_argument(
        "--rpkm-floor",
        type=float,
        default=0.1,
        help="If a rpkm-metric group's max < this floor, do NOT apply the filter (default: 0.1).",
    )
    ap.add_argument(
        "--group-cols",
        default=None,
        help="Optional comma-separated columns to define groups. If omitted, uses tool/mode (if present) + rank + taxid.",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t")

    group_cols = None
    if args.group_cols:
        group_cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]

    out = apply_bleed_filter(
        df,
        rpm_col=args.rpm_col,
        rpkm_col=args.rpkm_col,
        fraction=args.fraction,
        rpm_floor=args.rpm_floor,
        rpkm_floor=args.rpkm_floor,
        metric=args.metric,
        group_cols=group_cols,
    )
    out.to_csv(args.out, sep="\t", index=False)


def run_snakemake():
    inp = str(snakemake.input[0])
    outp = str(snakemake.output[0])

    rpm_col = getattr(snakemake.params, "rpm_col", "rpm")
    rpkm_col = getattr(snakemake.params, "rpkm_col", "rpkm")
    metric = getattr(snakemake.params, "metric", "auto")
    fraction = float(getattr(snakemake.params, "fraction", 0.005))
    rpm_floor = float(getattr(snakemake.params, "rpm_floor", 1.0))
    rpkm_floor = float(getattr(snakemake.params, "rpkm_floor", 0.1))

    group_cols = getattr(snakemake.params, "group_cols", None)
    if group_cols is not None:
        # allow passing list or comma-separated string
        if isinstance(group_cols, str):
            group_cols = [c.strip() for c in group_cols.split(",") if c.strip()]
        elif not isinstance(group_cols, list):
            raise ValueError(
                "snakemake.params.group_cols must be a list or comma-separated string."
            )

    df = pd.read_csv(inp, sep="\t")
    out = apply_bleed_filter(
        df,
        rpm_col=rpm_col,
        rpkm_col=rpkm_col,
        fraction=fraction,
        rpm_floor=rpm_floor,
        rpkm_floor=rpkm_floor,
        metric=metric,
        group_cols=group_cols,
    )
    out.to_csv(outp, sep="\t", index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
