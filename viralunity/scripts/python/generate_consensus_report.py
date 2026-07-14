#!/usr/bin/env python3
"""Generate a self-contained interactive HTML report for a ViralUnity consensus run.

Single source of truth for both delivery paths:
  - Snakemake ``script:`` directive (``run_snakemake``), and
  - the ``viralunity report`` CLI (``run_cli`` / imported ``write_report``).

The report visualizes ``<output>/assembly/assembly_stats_summary.csv`` plus the
per-base coverage tables under ``<output>/assembly/[<segment>/]coverage_stats/``.
Segmented mode is detected by the presence of a ``segment`` column in the stats
CSV. Every coverage path is reconstructed from ``output_dir`` + the CSV rows (via
``resolve_basewise_path``) so the CLI path, which has no ``snakemake`` object,
resolves inputs exactly like the workflow does.
"""

import argparse
import json
import logging
import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs
from plotly.subplots import make_subplots

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")
TEMPLATE_NAME = "report_template.html.j2"

# The percentage_above_*x columns are dropped from the displayed table; only the
# breadth metric (horizontal_coverage) is kept, still formatted as a percentage.
PERCENTAGE_ABOVE_COLS = [
    "percentage_above_10x",
    "percentage_above_100x",
    "percentage_above_1000x",
]
PERCENTAGE_COLS = PERCENTAGE_ABOVE_COLS + ["horizontal_coverage"]
INT_COLS = [
    "number_of_reads",
    "number_of_trim_paired_reads",
    "number_of_mapped_reads",
]

# Human-readable table headers, replacing the raw snake_case column names.
COLUMN_LABELS = {
    "sample_name": "Sample",
    "segment": "Segment",
    "number_of_reads": "Total reads",
    "number_of_trim_paired_reads": "QC-passed reads",
    "number_of_mapped_reads": "Mapped reads",
    "average_depth": "Mean depth",
    "horizontal_coverage": "Genome coverage",
}

FONT_FAMILY = 'system-ui, -apple-system, "Segoe UI", sans-serif'
GRID_COLOR = "rgba(137,135,129,0.25)"  # semi-transparent gray: reads on light + dark

# Pre-validated, colourblind-safe categorical palette (dataviz skill reference
# instance: worst adjacent CVD dE 24.2 light / 10.3 dark). Fixed order, never
# cycled; a run with >8 samples falls back to emphasis (single hue + hover).
PALETTE_LIGHT = [
    "#2a78d6",
    "#1baf7a",
    "#eda100",
    "#008300",
    "#4a3aa7",
    "#e34948",
    "#e87ba4",
    "#eb6834",
]
PALETTE_DARK = [
    "#3987e5",
    "#199e70",
    "#c98500",
    "#008300",
    "#9085e9",
    "#e66767",
    "#d55181",
    "#d95926",
]
SERIES_1_LIGHT = "#2a78d6"
EMPHASIS_MUTED = "#898781"  # recessive ink for >8-series line plots
MAX_CATEGORICAL = 8

FIG_WIDTH = 900
FIG_HEIGHT = 420


# --------------------------------------------------------------------------- #
# Data loading / shape
# --------------------------------------------------------------------------- #
def read_stats_summary(csv_path: str) -> pd.DataFrame:
    """Read the assembly stats summary CSV (has a header row)."""
    return pd.read_csv(csv_path)


def is_segmented(df: pd.DataFrame) -> bool:
    """Segmented runs carry a ``segment`` column (one row per sample x segment)."""
    return "segment" in df.columns


def resolve_basewise_path(output_dir: str, sample: str, segment: Optional[str] = None) -> str:
    """Reconstruct the per-base coverage path for a sample (+ segment).

    Unsegmented: ``<out>/assembly/coverage_stats/<sample>.table_cov_basewise.txt``
    Segmented:   ``<out>/assembly/<segment>/coverage_stats/<sample>.table_cov_basewise.txt``
    """
    base = os.path.join(output_dir, "assembly")
    if segment:
        base = os.path.join(base, segment)
    return os.path.join(base, "coverage_stats", f"{sample}.table_cov_basewise.txt")


def load_basewise_table(path: str) -> pd.DataFrame:
    """Load a per-base coverage table (tab-separated, no header, 3 columns).

    Missing files log a warning and return an empty frame so a partial output
    directory still renders for the samples that do have data.
    """
    cols = ["reference_id", "position", "depth"]
    if not os.path.exists(path):
        logger.warning("Coverage file not found, skipping: %s", path)
        return pd.DataFrame(columns=cols)
    return pd.read_csv(path, sep=r"\s+", header=None, names=cols, engine="python")


# --------------------------------------------------------------------------- #
# Run metadata
# --------------------------------------------------------------------------- #
# NOTE: this module runs inside Snakemake's ``--use-conda`` env (envs/report.yaml),
# which does not install the viralunity package, so the config-key/data-type
# strings below are inlined rather than imported from viralunity.constants (they
# mirror ``ConfigKeys.DATA``/``ConfigKeys.SCHEME`` and ``DataType``).
_ILLUMINA = "illumina"
_NANOPORE = "nanopore"


def build_report_metadata(config: Optional[dict], output_dir: str) -> Dict[str, object]:
    """Derive report run-metadata from the run config, or infer from the output dir.

    The config YAML that drives the whole pipeline already carries everything the
    report needs:

    * ``data`` (``illumina``/``nanopore``) -> platform, library layout, and whether
      a QC step ran. Illumina consensus runs are paired-end and run a fastp QC
      step; Nanopore runs are single-end and run no QC.
    * ``scheme`` (a primer BED path, or ``"NA"``) -> primer scheme.

    When no config is available (``viralunity report`` on an older output dir with
    no config to hand), fall back to inferring the platform from the presence of a
    ``qc/`` directory (fastp writes one on Illumina; Nanopore has none). This is a
    last-resort heuristic; a declared ``data`` field always wins.

    Returns ``{platform, library_layout, primer_scheme, qc_performed}``.
    """
    platform: Optional[str] = None
    primer_scheme: Optional[str] = None
    if config:
        data = config.get("data")
        if data in (_ILLUMINA, _NANOPORE):
            platform = data
        scheme = config.get("scheme")
        if scheme and str(scheme).strip().upper() != "NA":
            primer_scheme = str(scheme)
    if platform is None:
        platform = _ILLUMINA if os.path.isdir(os.path.join(output_dir, "qc")) else _NANOPORE
    paired = platform == _ILLUMINA
    return {
        "platform": platform,
        "library_layout": "paired" if paired else "single",
        "primer_scheme": primer_scheme,
        "qc_performed": paired,
    }


# --------------------------------------------------------------------------- #
# Transforms
# --------------------------------------------------------------------------- #
def downsample_min_pool(
    positions: np.ndarray, depths: np.ndarray, max_points: int = 2000
) -> Tuple[np.ndarray, np.ndarray]:
    """Downsample a coverage track to <= ``max_points`` using min-pooling.

    The position range is split into ``max_points`` contiguous bins; per bin the
    (position, depth) pair with the MINIMUM depth is kept. Keeping the minimum
    (not the mean) preserves coverage dips instead of averaging them away.
    """
    positions = np.asarray(positions)
    depths = np.asarray(depths)
    n = len(positions)
    if n <= max_points:
        return positions, depths

    out_pos: List = []
    out_depth: List = []
    for idx in np.array_split(np.arange(n), max_points):
        if len(idx) == 0:
            continue
        j = idx[int(np.argmin(depths[idx]))]
        out_pos.append(positions[j])
        out_depth.append(depths[j])
    return np.array(out_pos), np.array(out_depth)


def dedupe_and_sum_reads(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per sample for the global reads view.

    ``number_of_reads`` / ``number_of_trim_paired_reads`` are identical across a
    sample's segment rows -> take the first. ``number_of_mapped_reads`` is
    per-segment -> sum across segments. Unsegmented input is a no-op per sample.
    """
    agg = {
        "number_of_reads": "first",
        "number_of_trim_paired_reads": "first",
        "number_of_mapped_reads": "sum",
    }
    return df.groupby("sample_name", as_index=False, sort=False).agg(agg)


# --------------------------------------------------------------------------- #
# Display formatting
# --------------------------------------------------------------------------- #
def _fmt_pct(value) -> str:
    try:
        return f"{float(value) * 100:.1f}%"
    except (TypeError, ValueError):
        return str(value)


def _fmt_depth(value) -> str:
    try:
        return f"{float(value):.1f}"
    except (TypeError, ValueError):
        return str(value)


def _fmt_int(value) -> str:
    try:
        return f"{int(round(float(value)))}"
    except (TypeError, ValueError):
        return str(value)


def build_stats_table_html(df: pd.DataFrame) -> str:
    """Render the stats summary as a sortable HTML table.

    The ``percentage_above_*x`` columns are dropped and headers are humanized.
    ``horizontal_coverage`` is shown x100 with a ``%`` suffix; ``average_depth``
    is rounded; read counts render as integers. Sorting is handled by inline JS
    in the template (``data-sort`` carries the raw numeric value).
    """
    columns = [c for c in df.columns if c not in PERCENTAGE_ABOVE_COLS]
    header_cells = "".join(
        f'<th data-col="{i}" onclick="sortTable(this)">{COLUMN_LABELS.get(col, col)}'
        f'<span class="sort-arrow"></span></th>'
        for i, col in enumerate(columns)
    )

    body_rows = []
    for _, row in df.iterrows():
        cells = []
        for col in columns:
            raw = row[col]
            if col in PERCENTAGE_COLS:
                display, sort_val = _fmt_pct(raw), raw
            elif col == "average_depth":
                display, sort_val = _fmt_depth(raw), raw
            elif col in INT_COLS:
                display, sort_val = _fmt_int(raw), raw
            else:
                display, sort_val = str(raw), raw
            cells.append(f'<td data-sort="{sort_val}">{display}</td>')
        body_rows.append("<tr>" + "".join(cells) + "</tr>")

    return (
        '<table class="stats-table" id="stats-table">'
        f"<thead><tr>{header_cells}</tr></thead>"
        f'<tbody>{"".join(body_rows)}</tbody></table>'
    )


# --------------------------------------------------------------------------- #
# Figure builders
# --------------------------------------------------------------------------- #
def _base_layout(**kwargs) -> dict:
    # Transparent surfaces + muted grid so the chart blends into the card in
    # both light and dark; the page's theme toggle relayouts font/grid colour.
    layout = dict(
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        margin=dict(l=60, r=30, t=50, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(
            family='system-ui, -apple-system, "Segoe UI", sans-serif',
            size=13,
            color="#52514e",
        ),
        colorway=PALETTE_LIGHT,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        xaxis=dict(gridcolor="rgba(137,135,129,0.25)", zeroline=False),
        yaxis=dict(gridcolor="rgba(137,135,129,0.25)", zeroline=False),
    )
    layout.update(kwargs)
    return layout


def _log_depth_range(max_value: float) -> list:
    """Log10 axis range for a depth chart, generalizable across run scales.

    Floor is 1 (10^0) so a depth of 0/1 sits on the axis bottom and the 20x/100x
    guides are always inside the range; the top is 1.5x the data max but never
    below 150, so the 100x guide never hugs the ceiling on shallow runs.
    """
    top = max(150.0, float(max_value) * 1.5)
    return [0.0, math.log10(top)]


def _clamp_for_log(depths) -> np.ndarray:
    """Clamp depths to >= 1 so a log axis can render zeros (shown at the floor)."""
    return np.maximum(np.asarray(depths, dtype=float), 1.0)


def _add_depth_guides(fig: go.Figure) -> None:
    """20x / 100x horizontal guide lines with literal '20x'/'100x' annotations."""
    fig.add_hline(
        y=20,
        line_width=1,
        line_dash="dot",
        line_color=EMPHASIS_MUTED,
        annotation_text="20x",
        annotation_position="right",
    )
    fig.add_hline(
        y=100,
        line_width=1,
        line_dash="dot",
        line_color=EMPHASIS_MUTED,
        annotation_text="100x",
        annotation_position="right",
    )


def build_reads_histogram(df: pd.DataFrame) -> go.Figure:
    """Reads per sample as two stacked panels sharing the sample x-axis.

    Top panel: total vs QC-passed reads. These are a subset relationship
    (QC-passed <= total), so they are overlaid — the QC-passed bar sits in front
    of the total bar and the remainder is the reads dropped in QC. Bottom panel:
    mapped reads (a single measure, its own hue).
    """
    samples = df["sample_name"].tolist()
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.14,
        subplot_titles=("Total and QC-passed reads", "Mapped reads"),
    )
    fig.add_bar(
        x=samples,
        y=df["number_of_reads"],
        name="Total reads",
        marker_color=PALETTE_LIGHT[0],
        opacity=0.4,
        row=1,
        col=1,
    )
    fig.add_bar(
        x=samples,
        y=df["number_of_trim_paired_reads"],
        name="QC-passed reads",
        marker_color=PALETTE_LIGHT[0],
        row=1,
        col=1,
    )
    fig.add_bar(
        x=samples,
        y=df["number_of_mapped_reads"],
        name="Mapped reads",
        marker_color=PALETTE_LIGHT[1],
        row=2,
        col=1,
    )
    fig.update_layout(
        width=FIG_WIDTH,
        height=560,
        margin=dict(l=60, r=30, t=60, b=60),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family=FONT_FAMILY, size=13, color="#52514e"),
        colorway=PALETTE_LIGHT,
        barmode="overlay",
        legend=dict(orientation="h", yanchor="bottom", y=1.06, xanchor="left", x=0),
    )
    fig.update_xaxes(gridcolor=GRID_COLOR, zeroline=False)
    fig.update_yaxes(gridcolor=GRID_COLOR, zeroline=False, rangemode="tozero", title_text="Reads")
    fig.update_xaxes(title_text="Sample", row=2, col=1)
    return fig


def build_coverage_bar_chart(df: pd.DataFrame, segmented: bool) -> go.Figure:
    """Log-y bar of average depth (one bar per sample x segment row).

    A single measure -> one hue for every bar. 20x/100x guide lines annotated.
    """
    if segmented:
        labels = [f"{s} | {seg}" for s, seg in zip(df["sample_name"], df["segment"])]
    else:
        labels = df["sample_name"].tolist()
    max_depth = float(df["average_depth"].max()) if len(df) else 100.0
    fig = go.Figure()
    fig.add_bar(x=labels, y=df["average_depth"], marker_color=SERIES_1_LIGHT, name="Average depth")
    fig.update_layout(
        _base_layout(showlegend=False),
        yaxis_type="log",
        yaxis_range=_log_depth_range(max_depth),
        yaxis_title="Mean depth (log)",
        xaxis_title="Sample" + (" x segment" if segmented else ""),
    )
    _add_depth_guides(fig)
    return fig


def build_aggregated_coverage_line_plot(
    per_sample_series: Dict[str, Tuple[np.ndarray, np.ndarray]], title: str
) -> go.Figure:
    """Aggregated coverage: x = genome position, one line per sample.

    <= 8 samples get the fixed categorical palette; more than 8 fall back to a
    single recessive hue (emphasis + hover) so identity is never colour-alone
    past the CVD-safe ceiling (the stats table carries per-sample identity).
    """
    emphasis = len(per_sample_series) > MAX_CATEGORICAL
    max_depth = max(
        (float(np.max(d)) for _, d in per_sample_series.values() if len(d)), default=100.0
    )
    fig = go.Figure()
    for sample, (positions, depths) in per_sample_series.items():
        kwargs = dict(mode="lines", name=sample, line=dict(width=2))
        if emphasis:
            kwargs["line"]["color"] = EMPHASIS_MUTED
            kwargs["opacity"] = 0.6
        fig.add_scatter(x=positions, y=_clamp_for_log(depths), **kwargs)
    fig.update_layout(
        _base_layout(showlegend=not emphasis, hovermode="x unified"),
        title=title,
        yaxis_type="log",
        yaxis_range=_log_depth_range(max_depth),
        yaxis_title="Depth (log)",
        xaxis_title="Genome position",
    )
    _add_depth_guides(fig)
    return fig


def build_sample_detail_plot(sample: str, label: str, positions, depths) -> go.Figure:
    """Single-sample (single-segment) coverage line, used by the lazy path."""
    max_depth = float(np.max(depths)) if len(depths) else 100.0
    fig = go.Figure()
    fig.add_scatter(
        x=list(positions),
        y=list(_clamp_for_log(depths)),
        mode="lines",
        line=dict(width=2, color=SERIES_1_LIGHT),
        name=label or sample,
    )
    title = sample if not label else f"{sample} — {label}"
    fig.update_layout(
        _base_layout(showlegend=False, hovermode="x unified"),
        title=title,
        yaxis_type="log",
        yaxis_range=_log_depth_range(max_depth),
        yaxis_title="Depth (log)",
        xaxis_title="Genome position",
    )
    _add_depth_guides(fig)
    return fig


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def _ordered_unique(values) -> List[str]:
    seen, out = set(), []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _load_coverage_cache(
    output_dir: str, df: pd.DataFrame, segmented: bool
) -> Dict[Tuple[str, Optional[str]], Tuple[np.ndarray, np.ndarray]]:
    """Load + downsample every (sample, segment) coverage track exactly once."""
    cache: Dict[Tuple[str, Optional[str]], Tuple[np.ndarray, np.ndarray]] = {}
    for _, row in df.iterrows():
        sample = row["sample_name"]
        segment = row["segment"] if segmented else None
        table = load_basewise_table(resolve_basewise_path(output_dir, sample, segment))
        if table.empty:
            cache[(sample, segment)] = (np.array([]), np.array([]))
            continue
        pos, depth = downsample_min_pool(table["position"].to_numpy(), table["depth"].to_numpy())
        cache[(sample, segment)] = (pos, depth)
    return cache


def _fig_to_html(fig: go.Figure) -> str:
    return pio.to_html(fig, include_plotlyjs=False, full_html=False)


def render_report(output_dir: str, metadata: Optional[dict] = None) -> str:
    """Build the full self-contained HTML report string for a consensus run.

    ``metadata`` is the run-metadata dict from :func:`build_report_metadata`; when
    omitted it is inferred from the output dir (the CLI-on-an-old-dir path).
    """
    if metadata is None:
        metadata = build_report_metadata(None, output_dir)
    stats_path = os.path.join(output_dir, "assembly", "assembly_stats_summary.csv")
    df = read_stats_summary(stats_path)
    segmented = is_segmented(df)
    samples = _ordered_unique(df["sample_name"].tolist())
    segments = _ordered_unique(df["segment"].tolist()) if segmented else [None]

    cache = _load_coverage_cache(output_dir, df, segmented)

    # Global figures.
    stats_table_html = build_stats_table_html(df)
    reads_fig_html = _fig_to_html(build_reads_histogram(dedupe_and_sum_reads(df)))
    depth_fig_html = _fig_to_html(build_coverage_bar_chart(df, segmented))

    # Aggregated coverage: one panel per segment (segments have different lengths).
    aggregated_panels = []
    for segment in segments:
        series = {
            s: cache[(s, segment)]
            for s in samples
            if len(cache.get((s, segment), (np.array([]),))[0]) > 0
        }
        title = f"Aggregated coverage — {segment}" if segment else "Aggregated coverage"
        fig_html = _fig_to_html(build_aggregated_coverage_line_plot(series, title))
        aggregated_panels.append({"segment": segment or "", "html": fig_html})

    # Per-sample coverage data, embedded once as JSON for the lazy detail section.
    coverage_json: Dict[str, list] = {}
    for sample in samples:
        entries = []
        for segment in segments:
            pos, depth = cache.get((sample, segment), (np.array([]), np.array([])))
            if len(pos) == 0:
                continue
            entries.append(
                {
                    "label": segment or "",
                    "x": [int(p) for p in pos],
                    "y": [float(d) for d in depth],
                }
            )
        coverage_json[sample] = entries

    stats_by_sample = _stats_rows_by_sample(df, segmented)

    # On Illumina (paired-end), total/QC-passed counts are read pairs while mapped
    # reads are counted individually, so the two reads panels are on different
    # units. Driven by the declared library layout, never inferred from ratios.
    paired = metadata["library_layout"] == "paired"
    reads_note = (
        "Top panel: read-pair counts. Bottom panel: individual read counts." if paired else ""
    )

    env = Environment(
        loader=FileSystemLoader(TEMPLATE_DIR),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template(TEMPLATE_NAME)
    return template.render(
        plotly_js=get_plotlyjs(),
        metadata=metadata,
        segmented=segmented,
        samples=samples,
        segments=[s or "" for s in segments],
        stats_table_html=stats_table_html,
        reads_fig_html=reads_fig_html,
        reads_note=reads_note,
        depth_fig_html=depth_fig_html,
        aggregated_panels=aggregated_panels,
        coverage_json=json.dumps(coverage_json),
        stats_by_sample=json.dumps(stats_by_sample),
        palette_light=PALETTE_LIGHT,
        palette_dark=PALETTE_DARK,
    )


def _stats_rows_by_sample(df: pd.DataFrame, segmented: bool) -> Dict[str, list]:
    """Formatted per-sample stat rows for the per-sample detail panel."""
    out: Dict[str, list] = {}
    for _, row in df.iterrows():
        rendered = {}
        for col in df.columns:
            if col in PERCENTAGE_ABOVE_COLS:
                continue
            label = COLUMN_LABELS.get(col, col)
            raw = row[col]
            if col in PERCENTAGE_COLS:
                rendered[label] = _fmt_pct(raw)
            elif col == "average_depth":
                rendered[label] = _fmt_depth(raw)
            elif col in INT_COLS:
                rendered[label] = _fmt_int(raw)
            else:
                rendered[label] = str(raw)
        out.setdefault(row["sample_name"], []).append(rendered)
    return out


def write_report(output_dir: str, dest: str, metadata: Optional[dict] = None) -> None:
    """Render the report and write it to ``dest`` (the shared entry point)."""
    html = render_report(output_dir, metadata)
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    with open(dest, "w") as fh:
        fh.write(html)
    logger.info("Wrote consensus report: %s", dest)


# --------------------------------------------------------------------------- #
# Entry points
# --------------------------------------------------------------------------- #
def run_cli():
    ap = argparse.ArgumentParser(
        description="Generate an interactive HTML report for a consensus output directory."
    )
    ap.add_argument("--input", dest="input_dir", required=True, help="Consensus output directory.")
    ap.add_argument(
        "--output",
        dest="output_path",
        default=None,
        help="Destination HTML path (default: <input>/report.html).",
    )
    ap.add_argument(
        "--config-file",
        dest="config_file",
        default=None,
        help="Run config YAML (the one used to launch the pipeline); supplies "
        "platform/library-layout/primer-scheme metadata. Inferred from the "
        "output dir when omitted.",
    )
    args = ap.parse_args()
    dest = args.output_path or os.path.join(args.input_dir, "report.html")
    metadata = build_report_metadata(load_run_config(args.config_file), args.input_dir)
    write_report(args.input_dir, dest, metadata)


def load_run_config(config_file: Optional[str]) -> Optional[dict]:
    """Load the run config YAML if a path is given (CLI path only).

    Imported lazily so the module still loads under the Snakemake ``report.yaml``
    conda env, which does not ship PyYAML; ``run_snakemake`` uses the pre-parsed
    ``snakemake.config`` instead.
    """
    if not config_file:
        return None
    import yaml

    with open(config_file) as fh:
        return yaml.safe_load(fh)


def run_snakemake():
    output_dir = str(snakemake.params[0])  # noqa: F821
    dest = str(snakemake.output[0])  # noqa: F821
    metadata = build_report_metadata(dict(snakemake.config), output_dir)  # noqa: F821
    write_report(output_dir, dest, metadata)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
