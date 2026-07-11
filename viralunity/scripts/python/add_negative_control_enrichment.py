#!/usr/bin/env python3
"""Add negative-control enrichment metrics to a ViralUnity taxa summary.

Replaces the previous Poisson-based ``apply_negative_background_filter.py``
with interpretable statistics that are valid for overdispersed, batch-specific
viral contamination:

  * fold_enrichment  = (sample_metric + pc) / (control_mean + pc)
  * log10_ratio       = log10((sample_metric + pc) / (control_mean + pc))
  * z_score          = (sample_metric - control_mean) / control_sd
                       [only when n_controls >= 2 AND control_sd > 0]

The *decision metric* is ``rpkm`` when that column is present and non-NA for
the taxon, else ``rpm``.  This is recorded in the ``neg_metric`` column.

Pass/fail logic (``neg_pass``):
  * n_controls == 0  → neg_pass = NA  (no filter; bleed-only mode)
  * n_controls >= 2  → neg_pass = (z_score >= z_score_threshold);
                       if control_sd == 0, falls back to log10_ratio gate
  * n_controls == 1  → neg_pass = (log10_ratio >= log10_ratio_threshold)
  * taxa absent from controls → control_mean = 0, computed against pseudocount

Control statistics use *zero-fill*: a taxon not detected in a given control
contributes a metric of 0 for that control, so control_mean / control_sd /
control_median are computed over all ``n_negative_controls`` values (not only
the controls where the taxon appears). This keeps the z-score denominator
consistent with ``n_controls`` and prevents a contaminant seen in one of many
controls from overstating the control baseline.

The downstream lineage-aware Krona filter (filter_krona_by_pass_taxids.py)
treats NA ``neg_pass`` as keep (conservative), so the existing contract is
preserved.

Output column additions:
  is_negative_control, n_negative_controls, neg_metric, control_mean,
  control_sd, control_median, control_max, fold_enrichment, log10_ratio,
  z_score, log10_ratio_threshold_used, z_score_threshold_used,
  enrichment_pseudocount, neg_decision, neg_pass
"""

import argparse
import math
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Sequence

import pandas as pd

try:
    from viralunity.scripts.python.apply_max_rpm_bleed_filter import infer_group_cols
except ImportError:
    # Running inside Snakemake's `script:` directive (per-rule conda env without
    # the viralunity package installed); the script's own dir is on sys.path.
    from apply_max_rpm_bleed_filter import infer_group_cols

# ─────────────────────────────────────────────────────────────────────────────
# Core enrichment functions (unit-testable, no Snakemake dependencies)
# ─────────────────────────────────────────────────────────────────────────────


def calculate_fold_enrichment(
    sample_metric: float,
    control_mean: float,
    pseudocount: float = 1.0,
) -> float:
    """Return (sample + pc) / (control_mean + pc)."""
    return (sample_metric + pseudocount) / (control_mean + pseudocount)


def calculate_log10_ratio(
    sample_metric: float,
    control_mean: float,
    pseudocount: float = 1.0,
) -> float:
    """Return log10((sample + pc) / (control_mean + pc))."""
    numerator = sample_metric + pseudocount
    denominator = control_mean + pseudocount
    return math.log10(numerator / denominator)


def calculate_z_score(
    sample_metric: float,
    control_metrics: Sequence[float],
    min_controls: int = 2,
) -> Optional[float]:
    """Return the z-score of *sample_metric* relative to *control_metrics*.

    Returns None when:
      * fewer than *min_controls* values are supplied, OR
      * the control standard deviation is 0.
    """
    if len(control_metrics) < min_controls:
        return None
    mean = statistics.mean(control_metrics)
    sd = statistics.stdev(control_metrics)  # sample stdev (n-1)
    if sd == 0:
        return None
    return (sample_metric - mean) / sd


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def _build_ctrl_stats(
    control_rows: pd.DataFrame,
    group_cols: List[str],
    n_controls: int,
) -> Dict[tuple, dict]:
    """Return a dict mapping group-key tuples to per-taxon control statistics.

    A taxon absent from a given negative control contributes a metric of 0 for
    that control (zero-fill), so control_mean/SD/median are computed over all
    ``n_controls`` values rather than only the controls where the taxon was
    detected. This keeps the z-score denominator (the number of control values)
    consistent with the ``n_controls`` the pass/fail decision branch keys off,
    and stops a contaminant seen in a single control from setting control_mean
    to that lone value.

    Uses an explicit Python loop to avoid pandas MultiIndex surprises when
    ``groupby().apply()`` returns a Series from the aggregation function.
    """
    metric_lists: Dict[tuple, List[float]] = defaultdict(list)
    for _, row in control_rows.iterrows():
        key = tuple(row[c] for c in group_cols)
        val = row["_metric"]
        if not (isinstance(val, float) and math.isnan(val)):
            metric_lists[key].append(float(val))

    stats: Dict[tuple, dict] = {}
    for key, vals in metric_lists.items():
        # Zero-fill the controls where this taxon was not detected.
        padded = vals + [0.0] * max(0, n_controls - len(vals))
        stats[key] = {
            "control_mean": statistics.mean(padded),
            "control_sd": statistics.stdev(padded) if len(padded) >= 2 else None,
            "control_median": statistics.median(padded),
            "control_max": max(padded),
            "_ctrl_vals": padded,
        }
    return stats


def _add_na_enrichment_cols(
    out: pd.DataFrame,
    pseudocount: float,
    z_thresh: float,
    l2r_thresh: float,
) -> None:
    """In-place: add all output columns as NA for the zero-control case."""
    out["is_negative_control"] = False
    out["n_negative_controls"] = 0
    out["control_mean"] = pd.NA
    out["control_sd"] = pd.NA
    out["control_median"] = pd.NA
    out["control_max"] = pd.NA
    out["fold_enrichment"] = pd.NA
    out["log10_ratio"] = pd.NA
    out["z_score"] = pd.NA
    out["enrichment_pseudocount"] = pseudocount
    out["z_score_threshold_used"] = z_thresh
    out["log10_ratio_threshold_used"] = l2r_thresh
    out["neg_decision"] = "none"
    out["neg_pass"] = pd.NA


# ─────────────────────────────────────────────────────────────────────────────
# Main enrichment logic
# ─────────────────────────────────────────────────────────────────────────────


def apply_negative_control_enrichment(
    df: pd.DataFrame,
    negatives: List[str],
    pseudocount: float = 1.0,
    z_score_threshold: float = 3.0,
    log10_ratio_threshold: float = 1.0,
    group_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Add enrichment, log10-ratio, z-score, and neg_pass columns to *df*.

    Parameters
    ----------
    df:
        Taxa summary TSV loaded as a DataFrame.  Must contain the columns
        produced by the add_RPM/add_RPKM steps: ``sample``, ``taxid``,
        ``rpm``, and optionally ``rpkm``.
    negatives:
        List of sample IDs to treat as negative controls.
    pseudocount:
        Additive pseudocount for fold-enrichment and log10-ratio.
    z_score_threshold:
        Minimum z-score required for ``neg_pass = True`` when
        ``n_negative_controls >= 2``.
    log10_ratio_threshold:
        Minimum log10-ratio required for ``neg_pass = True`` when
        ``n_negative_controls == 1`` or z-score is undefined.
    group_cols:
        Columns that define a unique taxon group. Defaults to the standard
        ``infer_group_cols`` result: (tool, mode if present) + rank + taxid.
    """
    if "sample" not in df.columns:
        raise ValueError("Input must contain a 'sample' column.")
    if "rpm" not in df.columns:
        raise ValueError("Input missing 'rpm' column. Run add_RPM_to_summary.py first.")

    out = df.copy()

    if group_cols is None:
        group_cols = infer_group_cols(out)

    # ── Determine decision metric per taxon group ─────────────────────────────
    # Use 'rpkm' if the column exists and has at least one non-NA value for the
    # group; otherwise fall back to 'rpm'. Use transform so no merge is needed.
    has_rpkm_col = "rpkm" in out.columns
    if has_rpkm_col:
        out["neg_metric"] = out.groupby(group_cols, dropna=False)["rpkm"].transform(
            lambda s: "rpkm" if s.notna().any() else "rpm"
        )
    else:
        out["neg_metric"] = "rpm"

    # ── Resolve the scalar metric for every row ───────────────────────────────
    def _metric_value(row):
        if row["neg_metric"] == "rpkm" and has_rpkm_col:
            val = row["rpkm"]
            if pd.notna(val):
                return float(val)
        return float(row["rpm"])

    out["_metric"] = out.apply(_metric_value, axis=1)

    # ── No negative controls ──────────────────────────────────────────────────
    if not negatives:
        _add_na_enrichment_cols(out, pseudocount, z_score_threshold, log10_ratio_threshold)
        out = out.drop(columns=["_metric"], errors="ignore")
        return out

    # ── Validate negative-control labels ─────────────────────────────────────
    neg_set = set(negatives)
    out["is_negative_control"] = out["sample"].isin(neg_set)

    neg_samples_present = sorted(out.loc[out["is_negative_control"], "sample"].unique())
    if not neg_samples_present:
        raise ValueError(
            "None of the provided negative controls appear in the input table's " "'sample' column."
        )

    n_controls = len(neg_samples_present)
    out["n_negative_controls"] = n_controls

    # ── Per-taxon control statistics ──────────────────────────────────────────
    control_rows = out[out["is_negative_control"]].copy()
    ctrl_stats = _build_ctrl_stats(control_rows, group_cols, n_controls)

    # Attach control statistics to each row via explicit lookup. A taxon absent
    # from *all* controls has n_controls zeros: mean/median 0 and (for >= 2
    # controls) SD 0, which routes it to the log10-ratio fallback below.
    out["control_mean"] = 0.0
    out["control_sd"] = 0.0 if n_controls >= 2 else pd.NA
    out["control_median"] = 0.0
    out["control_max"] = 0.0
    out["_ctrl_vals"] = pd.NA

    for idx, row in out.iterrows():
        key = tuple(row[c] for c in group_cols)
        if key in ctrl_stats:
            s = ctrl_stats[key]
            out.at[idx, "control_mean"] = s["control_mean"]
            out.at[idx, "control_sd"] = s["control_sd"]
            out.at[idx, "control_median"] = s["control_median"]
            out.at[idx, "control_max"] = s["control_max"]
            out.at[idx, "_ctrl_vals"] = s["_ctrl_vals"]

    # ── Enrichment metrics ────────────────────────────────────────────────────
    out["enrichment_pseudocount"] = pseudocount
    out["z_score_threshold_used"] = z_score_threshold
    out["log10_ratio_threshold_used"] = log10_ratio_threshold

    out["fold_enrichment"] = out.apply(
        lambda r: calculate_fold_enrichment(r["_metric"], r["control_mean"], pseudocount),
        axis=1,
    )
    out["log10_ratio"] = out.apply(
        lambda r: calculate_log10_ratio(r["_metric"], r["control_mean"], pseudocount),
        axis=1,
    )

    # z-score: only meaningful for non-control rows
    out["z_score"] = pd.NA
    non_ctrl_mask = ~out["is_negative_control"]
    for idx in out[non_ctrl_mask].index:
        ctrl_vals = out.at[idx, "_ctrl_vals"]
        # Normalise: a taxon absent from all controls → n_controls zeros.
        if not isinstance(ctrl_vals, list):
            ctrl_vals = [0.0] * n_controls
        z = calculate_z_score(out.at[idx, "_metric"], ctrl_vals, min_controls=2)
        if z is not None:
            out.at[idx, "z_score"] = z

    # ── neg_pass decision ─────────────────────────────────────────────────────
    out["neg_decision"] = pd.NA
    out["neg_pass"] = pd.NA

    if n_controls >= 2:
        for idx in out[non_ctrl_mask].index:
            z = out.at[idx, "z_score"]
            if pd.notna(z):
                out.at[idx, "neg_decision"] = "z_score"
                out.at[idx, "neg_pass"] = bool(float(z) >= z_score_threshold)
            else:
                # z undefined (control SD == 0 or taxon absent from controls)
                l2r = out.at[idx, "log10_ratio"]
                out.at[idx, "neg_decision"] = "log10_ratio_fallback"
                out.at[idx, "neg_pass"] = bool(float(l2r) >= log10_ratio_threshold)
    else:
        # n_controls == 1: log10-ratio gate
        for idx in out[non_ctrl_mask].index:
            l2r = out.at[idx, "log10_ratio"]
            out.at[idx, "neg_decision"] = "log10_ratio"
            out.at[idx, "neg_pass"] = bool(float(l2r) >= log10_ratio_threshold)

    # Tidy up internal columns
    out = out.drop(columns=["_metric", "_ctrl_vals"], errors="ignore")

    return out


# ─────────────────────────────────────────────────────────────────────────────
# CLI / Snakemake entry points
# ─────────────────────────────────────────────────────────────────────────────


def run_cli() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Add negative-control enrichment metrics (fold-enrichment, log10-ratio, z-score) "
            "and a neg_pass gate to a ViralUnity taxa summary TSV."
        )
    )
    ap.add_argument(
        "--in", dest="inp", required=True, help="Input TSV (must contain sample, taxid, rpm)."
    )
    ap.add_argument("--out", required=True, help="Output TSV.")
    ap.add_argument(
        "--negatives",
        default="",
        help="Comma-separated negative-control sample IDs (empty = no filter).",
    )
    ap.add_argument(
        "--pseudocount",
        type=float,
        default=1.0,
        help="Pseudocount for fold-enrichment/log10-ratio (default: 1.0).",
    )
    ap.add_argument(
        "--z-score-threshold",
        type=float,
        default=3.0,
        help="Min z-score for neg_pass when n_controls >= 2 (default: 3.0).",
    )
    ap.add_argument(
        "--log10-ratio-threshold",
        type=float,
        default=1.0,
        help="Min log10-ratio for neg_pass when n_controls == 1 (default: 1.0).",
    )
    ap.add_argument(
        "--group-cols",
        default=None,
        help="Comma-separated grouping columns (default: infer from data).",
    )
    args = ap.parse_args()

    negatives = [x.strip() for x in args.negatives.split(",") if x.strip()]
    group_cols = None
    if args.group_cols:
        group_cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]

    df = pd.read_csv(args.inp, sep="\t")
    out = apply_negative_control_enrichment(
        df,
        negatives=negatives,
        pseudocount=args.pseudocount,
        z_score_threshold=args.z_score_threshold,
        log10_ratio_threshold=args.log10_ratio_threshold,
        group_cols=group_cols,
    )
    out.to_csv(args.out, sep="\t", index=False, na_rep="NA")
    sys.stderr.write(f"[add_negative_control_enrichment] wrote {len(out)} rows to {args.out}\n")


def run_snakemake() -> None:
    inp = str(snakemake.input[0])
    outp = str(snakemake.output[0])

    negatives = getattr(snakemake.params, "negatives", None) or []
    if isinstance(negatives, str):
        negatives = [x.strip() for x in negatives.split(",") if x.strip()]
    if not isinstance(negatives, list):
        raise ValueError("snakemake.params.negatives must be a list or comma-separated string.")

    pseudocount = float(getattr(snakemake.params, "pseudocount", 1.0))
    z_score_threshold = float(getattr(snakemake.params, "z_score_threshold", 3.0))
    log10_ratio_threshold = float(getattr(snakemake.params, "log10_ratio_threshold", 1.0))

    group_cols = getattr(snakemake.params, "group_cols", None)
    if isinstance(group_cols, str):
        group_cols = [c.strip() for c in group_cols.split(",") if c.strip()]

    df = pd.read_csv(inp, sep="\t")
    out = apply_negative_control_enrichment(
        df,
        negatives=negatives,
        pseudocount=pseudocount,
        z_score_threshold=z_score_threshold,
        log10_ratio_threshold=log10_ratio_threshold,
        group_cols=group_cols,
    )
    out.to_csv(outp, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
