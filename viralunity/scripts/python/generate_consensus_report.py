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
import html
import json
import logging
import math
import os
import re
import urllib.request
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

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
    "number_of_trim_paired_reads": "QC-passed",
    "number_of_mapped_reads": "Mapped %",
    "average_depth": "Mean depth",
    # "Horizontal coverage" (not "Genome coverage") to disambiguate breadth from
    # the segmented heatmap and from depth (spec §3).
    "horizontal_coverage": "Horizontal coverage",
}

# Columns rendered right-aligned with tabular figures (everything numeric); the
# sample/segment identity columns stay left-aligned.
NUMERIC_COLS = set(INT_COLS) | {"average_depth", "horizontal_coverage"} | set(PERCENTAGE_COLS)

FONT_FAMILY = 'system-ui, -apple-system, "Segoe UI", sans-serif'
GRID_COLOR = "rgba(137,135,129,0.25)"  # semi-transparent gray: reads on light + dark

# Pre-validated, colourblind-safe categorical palette (dataviz skill reference
# instance: worst adjacent CVD dE ~27.6 light / ~13.6 dark; see
# palette_validation_test.py). Fixed order, never cycled; a run with >8 samples
# falls back to emphasis (single hue + hover).
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

# Annotation-track hues (gene + primer), deliberately low-chroma and outside the
# categorical palette so a track never reads as a data series or collides with
# the depth-series blue (PALETTE_LIGHT[0]). Pools A/B share the primer hue (the
# lane label distinguishes them). Border hexes; the JS builds a translucent fill
# from each. Validated in palette_validation_test.py.
TRACK_GENE_COLOR = "#5b6b7f"
TRACK_GENE_COLOR_DARK = "#98a3b4"
TRACK_PRIMER_COLOR = "#4f7d6a"
TRACK_PRIMER_COLOR_DARK = "#8bb4a2"

FIG_WIDTH = 900
FIG_HEIGHT = 420


# --------------------------------------------------------------------------- #
# Report parameters (thresholds + colours) — the tweakables in the spec's §9
# --------------------------------------------------------------------------- #
DEFAULT_PASS_THRESHOLD = 0.90
DEFAULT_WARN_THRESHOLD = 0.70
DEFAULT_CHART_COLOR = "#2a78d6"
DEFAULT_COLORBAR_THICKNESS = 14

# Coverage status tiers. Always paired with a dot/label/bar (never colour alone),
# and held constant across light/dark themes (the accent + status hues do not
# flip; only the surrounding surfaces do).
STATUS_PASS_COLOR = "#1b9e5a"
STATUS_WARN_COLOR = "#e0a100"
STATUS_FAIL_COLOR = "#e34948"
TIER_COLORS = {"pass": STATUS_PASS_COLOR, "warn": STATUS_WARN_COLOR, "fail": STATUS_FAIL_COLOR}


@dataclass(frozen=True)
class ReportParams:
    """Tweakable report-generation parameters (thresholds + chart colours).

    ``pass_threshold`` / ``warn_threshold`` are coverage fractions in ``[0, 1]``.
    Every threshold-derived label the report shows (KPI tiles, the "low coverage
    only" filter, status legends) is computed from these, so the literal numbers
    90/70 never appear hardcoded in user-facing text (spec §9).
    """

    pass_threshold: float = DEFAULT_PASS_THRESHOLD
    warn_threshold: float = DEFAULT_WARN_THRESHOLD
    chart_color: str = DEFAULT_CHART_COLOR
    colorbar_thickness: int = DEFAULT_COLORBAR_THICKNESS

    def __post_init__(self):
        # Fail loudly on bad params rather than silently emitting a nonsensical
        # report (e.g. --pass-threshold 90 read as a fraction, or warn > pass
        # inverting the tiers, or a colour Plotly can't parse).
        for name in ("pass_threshold", "warn_threshold"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and 0.0 <= float(v) <= 1.0):
                raise ValueError(f"{name} must be a fraction in [0, 1], got {v!r}")
        if self.warn_threshold > self.pass_threshold:
            raise ValueError(
                f"warn_threshold ({self.warn_threshold}) must be <= pass_threshold "
                f"({self.pass_threshold})"
            )
        if not re.fullmatch(r"#[0-9A-Fa-f]{6}", str(self.chart_color)):
            raise ValueError(f"chart_color must be a #RRGGBB hex colour, got {self.chart_color!r}")
        if not (isinstance(self.colorbar_thickness, int) and self.colorbar_thickness > 0):
            raise ValueError(
                f"colorbar_thickness must be a positive integer, got {self.colorbar_thickness!r}"
            )

    def tier(self, coverage: float) -> str:
        """Status tier for a horizontal-coverage fraction: ``pass``/``warn``/``fail``."""
        try:
            cov = float(coverage)
        except (TypeError, ValueError):
            return "fail"
        if cov >= self.pass_threshold:
            return "pass"
        if cov >= self.warn_threshold:
            return "warn"
        return "fail"

    def tier_color(self, coverage: float) -> str:
        return TIER_COLORS[self.tier(coverage)]

    @property
    def pass_pct_label(self) -> str:
        """e.g. ``"90%"`` — the pass threshold as a compact percent string."""
        return f"{self.pass_threshold * 100:g}%"

    @property
    def warn_pct_label(self) -> str:
        """e.g. ``"70%"`` — the warn threshold as a compact percent string."""
        return f"{self.warn_threshold * 100:g}%"

    def as_client_dict(self) -> dict:
        """Params the client JS needs (thresholds + colours + status hues)."""
        return {
            "passThreshold": self.pass_threshold,
            "warnThreshold": self.warn_threshold,
            "passPctLabel": self.pass_pct_label,
            "warnPctLabel": self.warn_pct_label,
            "chartColor": self.chart_color,
            "colorbarThickness": self.colorbar_thickness,
            "tierColors": TIER_COLORS,
        }


# --------------------------------------------------------------------------- #
# Data loading / shape
# --------------------------------------------------------------------------- #
def read_stats_summary(csv_path: str) -> pd.DataFrame:
    """Read the assembly stats summary CSV (has a header row).

    ``keep_default_na=False`` so a segment literally named ``NA`` (influenza's
    neuraminidase segment) is read as the string ``"NA"``, not parsed to ``NaN``.
    The identity columns are then forced to ``str``: a segment named ``1``..``8``
    (numbered influenza segments, or numeric sample ids) would otherwise infer to
    ``int64`` and break every path/key built from it (``os.path.join`` rejects a
    non-str). Numeric metric columns keep their inferred numeric dtype.
    """
    df = pd.read_csv(csv_path, keep_default_na=False)
    for col in ("sample_name", "segment"):
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df


def is_segmented(df: pd.DataFrame) -> bool:
    """Segmented runs carry a ``segment`` column (one row per sample x segment)."""
    return "segment" in df.columns


def _safe_float(value, default: float = 0.0) -> float:
    """Parse a metric cell to float, degrading to ``default`` on garbage.

    The stats table already coerces its numeric cells defensively; the KPI and
    ordering paths use this so one non-numeric coverage/depth cell (a crafted or
    corrupt CSV) degrades that sample rather than crashing the whole report.
    """
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    return f if math.isfinite(f) else default


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
# Annotation tracks (gene GFF3 + primer BED)
# --------------------------------------------------------------------------- #
# Feature-type and label precedence are reporting choices (documented in
# REPORT.md); genes are drawn by default, falling back to CDS when a contig has
# no gene features.
GFF_FEATURE_TYPES = ("gene",)
GFF_FALLBACK_FEATURE_TYPES = ("CDS",)
GFF_LABEL_KEYS = ("Name", "gene", "gene_name", "locus_tag", "product", "ID")


def _sanitize_contig(s: str) -> str:
    """Mirror the nanopore workflow's contig sanitisation (``/\\|,~`` + space -> ``_``)."""
    return re.sub(r"[/\\|,~ ]", "_", s)


def _accession(s: str) -> str:
    """First whitespace/pipe-delimited token of a contig id (``NC_007373.1|...`` -> ``NC_007373.1``)."""
    return re.split(r"[\s|]", s.strip())[0] if s else s


def _contig_matches(feature_contig: str, coverage_contig: str) -> bool:
    """Whether a BED/GFF3 seqid refers to the coverage contig.

    Tolerant of the nanopore header sanitisation and of accession-only vs
    full-description headers, so an annotation file authored against the raw
    reference still matches a sanitised or pipe-delimited coverage contig.
    """
    if feature_contig == coverage_contig:
        return True
    if _sanitize_contig(feature_contig) == coverage_contig:
        return True
    acc = _accession(feature_contig)
    return bool(acc) and acc == _accession(coverage_contig)


def _parse_gff_attributes(field: str) -> Dict[str, str]:
    """Parse a GFF3 column-9 attribute string (``key=value;...``) into a dict."""
    attrs: Dict[str, str] = {}
    for part in field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, value = part.partition("=")
        attrs[key.strip()] = value.strip()
    return attrs


def _gff_rows_for_contig(path: str, contig: str, feature_types) -> List[dict]:
    """Return raw feature dicts of the requested types on ``contig`` (1-based)."""
    wanted = {t.lower() for t in feature_types}
    features: List[dict] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if not _contig_matches(cols[0], contig) or cols[2].lower() not in wanted:
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            attrs = _parse_gff_attributes(cols[8])
            name = next((attrs[k] for k in GFF_LABEL_KEYS if attrs.get(k)), cols[2])
            strand = cols[6] if cols[6] in ("+", "-") else "."
            features.append({"start": start, "end": end, "name": name, "strand": strand})
    features.sort(key=lambda f: (f["start"], f["end"]))
    return features


def parse_gff3(
    path: str,
    contig: str,
    feature_types=GFF_FEATURE_TYPES,
    fallback_feature_types=GFF_FALLBACK_FEATURE_TYPES,
) -> List[dict]:
    """Parse gene features for ``contig`` from a GFF3 file.

    Returns a list of ``{start, end, name, strand}`` with 1-based inclusive
    coordinates (GFF3 is already 1-based), sorted by start. Features of
    ``feature_types`` are used; if the contig has none, ``fallback_feature_types``
    is tried. A missing file yields ``[]`` with a warning.
    """
    if not os.path.exists(path):
        logger.warning("Gene annotation file not found, skipping: %s", path)
        return []
    features = _gff_rows_for_contig(path, contig, feature_types)
    if not features and fallback_feature_types:
        features = _gff_rows_for_contig(path, contig, fallback_feature_types)
    return features


# Pool lane labels; distinct pools (or amplicon parities) map to A, B, C, ...
_POOL_LANE_LABELS = [f"Pool {chr(ord('A') + i)}" for i in range(12)]


def _primer_side_and_amplicon(name: str):
    """Split an ARTIC-style primer name into (amplicon_key, side).

    ``scheme_7_LEFT`` / ``scheme_7_RIGHT_alt`` -> (``scheme_7``, ``LEFT``/``RIGHT``).
    Returns ``(None, None)`` when no LEFT/RIGHT token is present.
    """
    upper = name.upper()
    for side in ("LEFT", "RIGHT"):
        token = "_" + side
        idx = upper.find(token)
        if idx != -1:
            return name[:idx], side
    return None, None


def parse_primer_bed(path: str, contig: str) -> List[dict]:
    """Parse a primer-scheme BED into ordered lanes for ``contig``.

    Returns ``[{label, features:[{start, end, name}]}]``. ARTIC-style schemes
    are paired LEFT/RIGHT into amplicon spans and split into pool lanes
    (``Pool A``/``Pool B``/...), using the BED pool column when present, else
    amplicon-number parity. If pairing fails (no LEFT/RIGHT tokens or an
    unmatched primer), every primer is drawn individually in a single
    ``Primers`` lane. BED coordinates (0-based half-open) become 1-based
    inclusive. A missing file yields ``[]``.
    """
    if not os.path.exists(path):
        logger.warning("Primer scheme file not found, skipping: %s", path)
        return []

    rows: List[dict] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4 or not _contig_matches(cols[0], contig):
                continue
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError:
                continue
            rows.append(
                {
                    "start": start + 1,  # BED 0-based -> 1-based inclusive
                    "end": end,
                    "name": cols[3],
                    "pool": cols[4].strip() if len(cols) > 4 and cols[4].strip() else None,
                }
            )

    if not rows:
        return []

    amplicons: Dict[str, dict] = {}
    pairable = True
    for row in rows:
        key, side = _primer_side_and_amplicon(row["name"])
        if key is None:
            pairable = False
            break
        amp = amplicons.setdefault(
            key, {"start": row["start"], "end": row["end"], "pool": row["pool"], "sides": set()}
        )
        amp["start"] = min(amp["start"], row["start"])
        amp["end"] = max(amp["end"], row["end"])
        amp["sides"].add(side)
        if amp["pool"] is None:
            amp["pool"] = row["pool"]

    if not pairable or any({"LEFT", "RIGHT"} - amp["sides"] for amp in amplicons.values()):
        # Fallback: one lane, each primer drawn individually.
        return [
            {
                "label": "Primers",
                "features": [
                    {"start": r["start"], "end": r["end"], "name": r["name"]} for r in rows
                ],
            }
        ]

    # Assign each amplicon to a pool lane: distinct pool-column values in
    # first-seen order, else amplicon-number parity (odd -> A, even -> B). If an
    # amplicon's LEFT/RIGHT primers disagree on pool (a malformed scheme), the
    # first non-empty pool seen for that amplicon wins.
    pool_order: List[str] = []
    lanes: Dict[str, List[dict]] = {}
    for key, amp in amplicons.items():
        if amp["pool"] is not None:
            pool_id = amp["pool"]
        else:
            match = "".join(ch for ch in key if ch.isdigit())
            pool_id = "odd" if (match and int(match) % 2 == 1) else "even"
        if pool_id not in pool_order:
            pool_order.append(pool_id)
        lanes.setdefault(pool_id, []).append(
            {"start": amp["start"], "end": amp["end"], "name": key}
        )

    result = []
    for i, pool_id in enumerate(pool_order):
        label = _POOL_LANE_LABELS[i] if i < len(_POOL_LANE_LABELS) else f"Pool {i + 1}"
        features = sorted(lanes[pool_id], key=lambda f: f["start"])
        result.append({"label": label, "features": features})
    return result


# --------------------------------------------------------------------------- #
# Annotation from the run config + an NCBI fallback for a missing gene GFF3
# --------------------------------------------------------------------------- #
_EFETCH_FT = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={acc}&rettype=ft&retmode=text"
)


def _config_scheme_path(config: Optional[dict]) -> Optional[str]:
    """Primer-scheme BED path from the run config (``scheme``), or None/"NA"."""
    if not config:
        return None
    v = config.get("scheme")
    return str(v) if v and str(v).strip().upper() != "NA" else None


def _config_gene_path(config: Optional[dict], segment: Optional[str]) -> Optional[str]:
    """Gene-annotation GFF3 path from the config (``gene_annotation``).

    Accepts a plain path (unsegmented) or a ``{segment: path}`` dict (segmented,
    mirroring ``--segmented-gene-annotation``). Returns None when absent or "NA".
    """
    if not config:
        return None
    v = config.get("gene_annotation")
    if not v:
        return None
    if isinstance(v, dict):
        p = v.get(segment)
        return str(p) if p and str(p).strip().upper() != "NA" else None
    s = str(v).strip()
    return s if s and s.upper() != "NA" else None


def _parse_ncbi_feature_table(text: str, feature_types) -> List[dict]:
    """Parse an NCBI ``rettype=ft`` feature table into gene-like features.

    Returns ``{start, end, name, strand}`` (1-based inclusive, span across all
    intervals of a multi-exon feature), keeping features of ``feature_types``.
    The label prefers the ``gene`` qualifier, then ``locus_tag``, then ``product``.
    """
    wanted = {t.lower() for t in feature_types}

    def _is_coord(s: str) -> bool:
        try:
            int(s.lstrip("<>"))
            return True
        except ValueError:
            return False

    feats: List[dict] = []
    cur: Optional[dict] = None
    for line in text.splitlines():
        if not line or line.startswith(">Feature"):
            continue
        cols = line.split("\t")
        if (
            len(cols) >= 3
            and cols[0]
            and cols[1]
            and cols[2]
            and _is_coord(cols[0])
            and _is_coord(cols[1])
        ):
            a, b = int(cols[0].lstrip("<>")), int(cols[1].lstrip("<>"))
            cur = {
                "type": cols[2].strip().lower(),
                "lo": min(a, b),
                "hi": max(a, b),
                "strand": "+" if a <= b else "-",
                "q": {},
            }
            feats.append(cur)
        elif (
            cur
            and len(cols) >= 2
            and cols[0]
            and cols[1]
            and _is_coord(cols[0])
            and _is_coord(cols[1])
            and (len(cols) < 3 or not cols[2])
        ):
            a, b = int(cols[0].lstrip("<>")), int(cols[1].lstrip("<>"))
            cur["lo"], cur["hi"] = min(cur["lo"], a, b), max(cur["hi"], a, b)
        elif cur and len(cols) >= 5 and cols[0] == "":
            key, val = cols[3].strip().lower(), cols[4].strip()
            if key and val and key not in cur["q"]:
                cur["q"][key] = val

    out = []
    for f in feats:
        if f["type"] not in wanted:
            continue
        name = f["q"].get("gene") or f["q"].get("locus_tag") or f["q"].get("product") or f["type"]
        out.append({"start": f["lo"], "end": f["hi"], "name": name, "strand": f["strand"]})
    out.sort(key=lambda x: (x["start"], x["end"]))
    return out


# A nuccore accession looks like 1-2 letters, optional "_", digits, ".version"
# (MN908947.3, NC_007373.1, KP164568.1). Contig ids that don't match (de-novo
# contigs, custom names) are skipped so we never burn a network timeout on an id
# NCBI could never resolve.
_ACCESSION_RE = re.compile(r"^[A-Z]{1,2}_?\d{4,}\.\d+$")


def fetch_ncbi_gene_features(accession: str, timeout: int = 12) -> List[dict]:
    """Fetch gene features for a nuccore accession from NCBI (efetch feature table).

    Network- and parse-failures degrade to ``[]`` with a warning — the report is
    never blocked by an unreachable NCBI or an unannotated accession. Genes are
    preferred; a record with no gene features falls back to CDS. An id that does
    not look like a nuccore accession is skipped without a request.
    """
    if not accession or not _ACCESSION_RE.match(accession):
        if accession:
            logger.warning("Skipping NCBI annotation fetch for non-accession contig %r", accession)
        return []
    url = _EFETCH_FT.format(acc=quote(accession))
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:  # noqa: S310 (https, fixed host)
            text = resp.read().decode("utf-8", "replace")
    except Exception as exc:  # network / HTTP / decode
        logger.warning("NCBI annotation fetch failed for %s: %s", accession, exc)
        return []
    genes = _parse_ncbi_feature_table(text, GFF_FEATURE_TYPES)
    if not genes:
        genes = _parse_ncbi_feature_table(text, GFF_FALLBACK_FEATURE_TYPES)
    if not genes:
        logger.warning("NCBI returned no gene/CDS features for %s", accession)
    return genes


def _write_gff3_cache(path: str, contig: str, features: List[dict]) -> None:
    """Cache fetched gene features as a GFF3 next to the run (best-effort).

    Persisting under ``<output>/annotation/`` means a re-render (or the automatic
    report) reuses the fetch instead of hitting NCBI again. A read-only output
    dir just skips the cache.
    """
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n")
            for f in features:
                # Strip GFF3-structural chars from the label so a product name with
                # ';', a tab or a newline can't corrupt the cached record on re-read.
                name = re.sub(r"[;\t\r\n]+", " ", str(f["name"])).strip()
                fh.write(
                    f"{contig}\tNCBI\tgene\t{f['start']}\t{f['end']}\t.\t"
                    f"{f.get('strand', '.')}\t.\tName={name}\n"
                )
    except OSError as exc:
        logger.warning("Could not cache fetched annotation to %s: %s", path, exc)


def resolve_annotation_path(output_dir: str, kind: str, segment: Optional[str] = None) -> str:
    """Reconstruct the staged annotation path for a track ``kind``.

    The pipeline stages annotation files under ``<output>/annotation/`` (see the
    ``stage_*`` rules). Primer BED is a single file; gene GFF3 is per-segment in
    segmented mode.
    """
    base = os.path.join(output_dir, "annotation")
    if kind == "primer":
        return os.path.join(base, "primer_scheme.bed")
    if segment:
        return os.path.join(base, f"{segment}.gene_annotation.gff3")
    return os.path.join(base, "gene_annotation.gff3")


def _resolve_gene_features(
    output_dir: str,
    segment: Optional[str],
    contig: str,
    config: Optional[dict],
    fetch_missing: bool,
) -> List[dict]:
    """Gene features for a segment, in source-preference order.

    1. the staged ``<output>/annotation/[<seg>.]gene_annotation.gff3`` (pipeline);
    2. the config's ``gene_annotation`` path (a run whose annotation was never
       staged — e.g. re-reporting an older output dir);
    3. an NCBI fetch by the contig's accession (only when ``fetch_missing``),
       cached back to the staged path so later renders skip the network.
    """
    staged = resolve_annotation_path(output_dir, "gene", segment)
    if os.path.exists(staged):
        feats = parse_gff3(staged, contig)
        if feats:
            return feats
    cfg_path = _config_gene_path(config, segment)
    if cfg_path and os.path.exists(cfg_path):
        feats = parse_gff3(cfg_path, contig)
        if feats:
            return feats
    if fetch_missing:
        feats = fetch_ncbi_gene_features(_accession(contig))
        if feats:
            _write_gff3_cache(staged, contig, feats)
            return feats
    return []


def _resolve_primer_path(output_dir: str, config: Optional[dict]) -> Optional[str]:
    """Primer BED path: the staged copy, else the config ``scheme`` path."""
    staged = resolve_annotation_path(output_dir, "primer")
    if os.path.exists(staged):
        return staged
    cfg_path = _config_scheme_path(config)
    if cfg_path and os.path.exists(cfg_path):
        return cfg_path
    return None


def build_annotation_model(
    output_dir: str,
    segments: List[Optional[str]],
    contig_by_segment: Dict[Optional[str], str],
    config: Optional[dict] = None,
    fetch_missing: bool = False,
) -> dict:
    """Assemble the per-segment annotation-track model.

    For each segment with a known coverage contig, resolve gene features (staged
    GFF3 → config ``gene_annotation`` → optional NCBI fetch) and a primer BED
    (staged → config ``scheme``), and build ordered lanes matched to the contig.
    Returns ``{"by_segment": {seg_label: {"lanes": [...]}}, "has_genes": bool,
    "has_primers": bool}``. Segments with no drawable feature are omitted, so a
    run with no annotation anywhere yields empty lanes and false flags.
    """
    by_segment: Dict[str, dict] = {}
    has_genes = False
    has_primers = False
    primer_path = _resolve_primer_path(output_dir, config)

    for segment in segments:
        contig = contig_by_segment.get(segment)
        if not contig:
            continue
        lanes: List[dict] = []

        gene_features = _resolve_gene_features(output_dir, segment, contig, config, fetch_missing)
        if gene_features:
            lanes.append({"kind": "gene", "label": "Annotation", "features": gene_features})
            has_genes = True

        if primer_path:
            for lane in parse_primer_bed(primer_path, contig):
                if lane["features"]:
                    lanes.append({"kind": "primer", **lane})
                    has_primers = True

        if lanes:
            by_segment[segment or ""] = {"lanes": lanes}

    return {"by_segment": by_segment, "has_genes": has_genes, "has_primers": has_primers}


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
        logger.warning(
            "No library layout declared (no config passed); inferred platform=%s "
            "from the %spresence of a qc/ directory. The mapping-rate denominator "
            "is doubled only for paired-end (Illumina) runs, so a misidentified "
            "platform scales every mapping rate by 2x. Pass --config-file to avoid "
            "guessing.",
            platform,
            "" if platform == _ILLUMINA else "absence of a ",
        )
    paired = platform == _ILLUMINA
    return {
        "platform": platform,
        "library_layout": "paired" if paired else "single",
        "primer_scheme": primer_scheme,
        "qc_performed": paired,
    }


def _render_config_value(value) -> str:
    """Escaped HTML for one config value (scalars, lists, and nested dicts)."""
    if value is None or value == "":
        return '<span class="cfg-null">—</span>'
    if isinstance(value, dict):
        if not value:
            return '<span class="cfg-null">—</span>'
        rows = "".join(
            f'<div class="cfg-row"><span class="cfg-key">{html.escape(str(k))}</span>'
            f'<span class="cfg-val">{_render_config_value(v)}</span></div>'
            for k, v in value.items()
        )
        return f'<div class="cfg-nested">{rows}</div>'
    if isinstance(value, (list, tuple)):
        if not value:
            return '<span class="cfg-null">—</span>'
        items = "".join(f"<li>{_render_config_value(v)}</li>" for v in value)
        return f'<ul class="cfg-list">{items}</ul>'
    return html.escape(str(value))


def build_config_panel_html(config: Optional[dict]) -> str:
    """Render the full run config as escaped, grouped key–value blocks.

    Everything in the config is shown (a memory of exactly what was run); the
    per-rule ``*_cpus``/``*_ram`` keys are split into a separate "Resources"
    group. Returns ``""`` when no config is available (the panel is then omitted).
    """
    if not config:
        return ""
    resources = {k: v for k, v in config.items() if k.endswith(("_cpus", "_ram"))}
    params = {k: v for k, v in config.items() if k not in resources}

    def section(title: str, data: dict) -> str:
        if not data:
            return ""
        rows = "".join(
            f'<div class="cfg-row"><span class="cfg-key">{html.escape(str(k))}</span>'
            f'<span class="cfg-val">{_render_config_value(v)}</span></div>'
            for k, v in data.items()
        )
        return (
            f'<h3 class="cfg-section">{html.escape(title)}</h3><div class="cfg-block">{rows}</div>'
        )

    return section("Parameters", params) + section("Resources", resources)


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
def _json_for_script(obj) -> str:
    """JSON safe to embed inside an HTML ``<script>`` element.

    ``json.dumps`` does not escape ``<``/``>``/``&``, so a data value containing
    ``</script>`` (e.g. a crafted sample id) would otherwise terminate the script
    element during HTML parsing. Escaping them as ``\\uXXXX`` is valid JSON and
    parses back to the original characters, but cannot break out of the element.
    """
    return json.dumps(obj).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def _fmt_pct(value) -> str:
    try:
        return f"{float(value) * 100:.1f}%"
    except (TypeError, ValueError):
        return str(value)


def _fmt_depth(value) -> str:
    # Group the integer part (a comma every 3 digits) but keep the single
    # fractional digit ungrouped. Python's ``,`` format spec is locale-
    # independent (always a comma), so the separator never renders as a decimal
    # point under a de-DE-style ambient locale.
    try:
        return f"{float(value):,.1f}"
    except (TypeError, ValueError):
        return str(value)


def _fmt_int(value) -> str:
    # Thousands-grouped display for read counts. ``,`` is locale-independent;
    # the raw value is kept separately as the numeric sort key (``data-sort``).
    try:
        return f"{int(round(float(value))):,}"
    except (TypeError, ValueError):
        return str(value)


def _mapping_rate(mapped, qc_passed, paired: bool) -> Optional[float]:
    """Mapped-read rate (%) = mapped / QC-passed, unit-reconciled by layout.

    Mapped reads are counted individually; on a paired-end (Illumina) run the
    QC-passed count is read *pairs*, so the denominator is doubled to bring both
    onto the individual-read scale. Layout comes from the declared metadata, never
    from the mapped/total ratio. Returns ``None`` when the inputs are unusable.
    """
    try:
        mapped_f = float(mapped)
        qc = float(qc_passed)
    except (TypeError, ValueError):
        return None
    denom = qc * 2 if paired else qc
    if denom <= 0:
        return None
    return mapped_f / denom * 100.0


def _fmt_rate(pct: Optional[float]) -> str:
    """Format a mapping rate, keeping small on-target fractions legible."""
    if pct is None or not math.isfinite(pct):
        return "—"  # em dash
    if pct >= 1:
        return f"{pct:.1f}%"
    if pct >= 0.01:
        return f"{pct:.2f}%"
    if pct > 0:
        return f"{pct:.3f}%"
    return "0%"


def _coverage_cell_html(raw, params: ReportParams) -> str:
    """Coverage cell: tinted value + status dot are handled elsewhere; this cell
    carries the value and an inline proportion bar (width = coverage %)."""
    try:
        frac = float(raw)
    except (TypeError, ValueError):
        return f'<td class="num">{html.escape(str(raw))}</td>'
    tier = params.tier(frac)
    pct = max(0.0, min(100.0, frac * 100.0))
    display = html.escape(_fmt_pct(raw))
    return (
        f'<td class="num cov-cell" data-sort="{html.escape(str(raw))}">'
        f'<span class="cov-val cov-{tier}">{display}</span>'
        f'<span class="cov-bar cov-bar-{tier}" style="width:{pct:.1f}%"></span></td>'
    )


_ROLLUP_HEADERS = [
    ("Sample", "col-freeze"),
    ("Total reads", "th-num"),
    ("QC-passed", "th-num"),
    ("Mapped %", "th-num"),
    ("Mean depth", "th-num"),
    ("Horizontal coverage", "th-num"),
]


def _segmented_stats_table_html(df: pd.DataFrame, paired: bool, params: ReportParams) -> str:
    """Roll-up stats table for segmented runs (spec §3): one summary row per sample.

    Instead of ~96×8 flat rows, each sample gets one summary row — coverage = the
    **mean across its segments**, depth = the **median across its segments**, and
    Mapped % = the **whole-sample rate** (Σ mapped ÷ QC-passed, denominator doubled
    when paired; per the user's definition). A caret expands the row to its
    per-segment subrows (indented, each with its own status dot + per-segment
    Mapped %). Samples sort worst-mean-coverage-first. The segment selector chips
    (rendered in the template) drive the All ⇄ single-segment view client-side.
    """
    header_cells = []
    for i, (label, extra) in enumerate(_ROLLUP_HEADERS):
        is_cov = i == 5
        aria = "ascending" if is_cov else "none"
        asc = ' data-asc="true"' if is_cov else ""
        arrow = "▲" if is_cov else ""
        header_cells.append(
            f'<th class="{extra}" data-col="{i}" role="button" tabindex="0" '
            f'aria-sort="{aria}"{asc} onclick="sortTable(this)">'
            f'{html.escape(label)}<span class="sort-arrow">{arrow}</span></th>'
        )

    # Per-sample roll-up, preserving first-seen sample order before the worst-first
    # sort so ties are stable.
    samples = _ordered_unique(df["sample_name"].tolist())
    rollups = []
    for s in samples:
        rows = df[df["sample_name"] == s]
        covs = pd.to_numeric(rows["horizontal_coverage"], errors="coerce").dropna().tolist()
        depths = pd.to_numeric(rows["average_depth"], errors="coerce").dropna().tolist()
        qc = rows["number_of_trim_paired_reads"].iloc[0]
        total = rows["number_of_reads"].iloc[0]
        sum_mapped = pd.to_numeric(rows["number_of_mapped_reads"], errors="coerce").sum()
        rollups.append(
            {
                "sample": s,
                "mean_cov": float(np.mean(covs)) if covs else 0.0,
                "median_depth": float(np.median(depths)) if depths else 0.0,
                "total": total,
                "qc": qc,
                "rate": _mapping_rate(sum_mapped, qc, paired),
                "n_seg": len(rows),
                "rows": rows,
            }
        )
    rollups.sort(key=lambda r: r["mean_cov"])  # worst mean coverage first

    body_rows = []
    for sid, r in enumerate(rollups):
        s_esc = html.escape(str(r["sample"]))
        tier = params.tier(r["mean_cov"])
        rate_tip = html.escape(
            f'sum of segment mapped reads / {"2× " if paired else ""}QC-passed reads'
        )
        summary = (
            f'<tr class="seg-summary" data-role="summary" data-sid="{sid}" '
            f'data-sample="{s_esc}" data-segment="" data-coverage="{r["mean_cov"]:.6f}" '
            f'tabindex="0" role="button" aria-expanded="false">'
            f'<td class="col-freeze cell-sample">'
            f'<span class="seg-caret" aria-hidden="true"></span>'
            f'<span class="cov-dot cov-dot-{tier}" aria-hidden="true"></span>{s_esc}'
            f'<span class="seg-tag">{r["n_seg"]} segments</span></td>'
            f'<td class="num" data-sort="{html.escape(str(r["total"]))}">{html.escape(_fmt_int(r["total"]))}</td>'
            f'<td class="num" data-sort="{html.escape(str(r["qc"]))}">{html.escape(_fmt_int(r["qc"]))}</td>'
            f'<td class="num" data-sort="{html.escape(str(r["rate"] if r["rate"] is not None else -1))}" '
            f'title="{rate_tip}">{html.escape(_fmt_rate(r["rate"]))}</td>'
            f'<td class="num" data-sort="{r["median_depth"]}">{html.escape(_fmt_depth(r["median_depth"]))}</td>'
            + _coverage_cell_html(r["mean_cov"], params)
            + "</tr>"
        )
        body_rows.append(summary)
        # Per-segment subrows (hidden until the summary row is expanded).
        for _, seg_row in r["rows"].iterrows():
            seg = html.escape(str(seg_row.get("segment", "")))
            try:
                seg_cov = float(seg_row.get("horizontal_coverage"))
            except (TypeError, ValueError):
                seg_cov = 0.0
            seg_tier = params.tier(seg_cov)
            seg_rate = _mapping_rate(
                seg_row.get("number_of_mapped_reads"),
                seg_row.get("number_of_trim_paired_reads"),
                paired,
            )
            seg_tip = html.escape(
                f'{_fmt_int(seg_row.get("number_of_mapped_reads"))} mapped reads / '
                f'{"2× " if paired else ""}QC-passed reads'
            )
            body_rows.append(
                f'<tr class="seg-subrow row-hidden" data-role="subrow" data-parent="{sid}" '
                f'data-sample="{s_esc}" data-segment="{seg}" data-coverage="{seg_cov:.6f}">'
                f'<td class="col-freeze cell-sample cell-subrow">'
                f'<span class="cov-dot cov-dot-{seg_tier}" aria-hidden="true"></span>{seg}</td>'
                f'<td class="num cell-muted">—</td>'
                f'<td class="num cell-muted">—</td>'
                f'<td class="num" data-sort="{html.escape(str(seg_rate if seg_rate is not None else -1))}" '
                f'title="{seg_tip}">{html.escape(_fmt_rate(seg_rate))}</td>'
                f'<td class="num" data-sort="{html.escape(str(seg_row.get("average_depth")))}">'
                f'{html.escape(_fmt_depth(seg_row.get("average_depth")))}</td>'
                + _coverage_cell_html(seg_cov, params)
                + "</tr>"
            )

    return (
        '<table class="stats-table stats-rollup" id="stats-table">'
        f"<thead><tr>{''.join(header_cells)}</tr></thead>"
        f'<tbody>{"".join(body_rows)}</tbody></table>'
    )


def build_stats_table_html(
    df: pd.DataFrame, paired: bool = True, params: Optional[ReportParams] = None
) -> str:
    """Render the stats summary as a filterable, sortable data table (spec §3).

    The ``percentage_above_*x`` columns are dropped and headers are humanized;
    ``horizontal_coverage`` (labelled "Horizontal coverage") is shown x100 with a
    ``%`` suffix, coloured by status tier, and carries an inline proportion bar; a
    status **dot** precedes the sample name. ``average_depth`` is rounded; read
    counts are grouped; ``number_of_mapped_reads`` becomes a **mapping rate**
    (``mapped / QC-passed``, denominator doubled when ``paired``), raw count in the
    tooltip. Rows sort **worst-coverage-first** by default (server-side) so
    problems surface without interaction; every numeric column stays client-
    sortable. Each row carries ``data-sample``/``data-segment``/``data-coverage``
    for the search box and the "low coverage only" filter. Numeric cells are
    right-aligned with tabular figures.

    Segmented runs (a ``segment`` column) roll up to one summary row per sample
    with expandable per-segment subrows — see :func:`_segmented_stats_table_html`.
    """
    if params is None:
        params = ReportParams()
    if is_segmented(df):
        return _segmented_stats_table_html(df, paired, params)
    columns = [c for c in df.columns if c not in PERCENTAGE_ABOVE_COLS]
    # Default sort = coverage ascending (worst first). Stable so ties keep input
    # order. NaN/garbage coverage sinks to the top (treated as worst).
    if "horizontal_coverage" in df.columns:
        order = pd.to_numeric(df["horizontal_coverage"], errors="coerce")
        df = df.assign(_covsort=order).sort_values(
            "_covsort", ascending=True, kind="stable", na_position="first"
        )
    cov_idx = columns.index("horizontal_coverage") if "horizontal_coverage" in columns else -1

    header_cells = []
    for i, col in enumerate(columns):
        # The coverage column starts pre-sorted ascending; its header reflects that
        # so a click toggles to descending (worst<->best).
        is_cov = i == cov_idx
        aria = "ascending" if is_cov else "none"
        asc = ' data-asc="true"' if is_cov else ""
        arrow = "▲" if is_cov else ""
        cls = " ".join(
            c
            for c in ("th-num" if col in NUMERIC_COLS else "", "col-freeze" if i == 0 else "")
            if c
        )
        header_cells.append(
            f'<th class="{cls}" data-col="{i}" role="button" tabindex="0" '
            f'aria-sort="{aria}"{asc} onclick="sortTable(this)">'
            f"{html.escape(str(COLUMN_LABELS.get(col, col)))}"
            f'<span class="sort-arrow">{arrow}</span></th>'
        )

    body_rows = []
    for _, row in df.iterrows():
        sample = html.escape(str(row.get("sample_name", "")))
        segment = html.escape(str(row.get("segment", ""))) if "segment" in columns else ""
        try:
            cov_frac = float(row.get("horizontal_coverage"))
        except (TypeError, ValueError):
            cov_frac = 0.0
        tier = params.tier(cov_frac)
        cells = []
        for i, col in enumerate(columns):
            raw = row[col]
            if col == "sample_name":
                # Status dot + name, frozen as the first column.
                cells.append(
                    f'<td class="col-freeze cell-sample"><span class="cov-dot cov-dot-{tier}" '
                    f'aria-hidden="true"></span>{html.escape(str(raw))}</td>'
                )
                continue
            if col == "number_of_mapped_reads":
                rate = _mapping_rate(raw, row.get("number_of_trim_paired_reads"), paired)
                sort_val = html.escape(str(rate if rate is not None else -1))
                tip = html.escape(
                    f'{_fmt_int(raw)} mapped reads / {"2× " if paired else ""}QC-passed reads'
                )
                cells.append(
                    f'<td class="num" data-sort="{sort_val}" title="{tip}">'
                    f"{html.escape(_fmt_rate(rate))}</td>"
                )
                continue
            if col == "horizontal_coverage":
                cells.append(_coverage_cell_html(raw, params))
                continue
            # The numeric formatters fall back to str(value) on a non-numeric
            # value, so a corrupt/crafted stats CSV could otherwise inject markup
            # through a "numeric" column. Escape every emitted display value and
            # data-sort attribute; for clean numbers this is a no-op.
            if col in PERCENTAGE_COLS:
                display, sort_val = html.escape(_fmt_pct(raw)), html.escape(str(raw))
            elif col == "average_depth":
                display, sort_val = html.escape(_fmt_depth(raw)), html.escape(str(raw))
            elif col in INT_COLS:
                display, sort_val = html.escape(_fmt_int(raw)), html.escape(str(raw))
            else:
                display = sort_val = html.escape(str(raw))
            cls = "num" if col in NUMERIC_COLS else ""
            cells.append(f'<td class="{cls}" data-sort="{sort_val}">{display}</td>')
        body_rows.append(
            f'<tr data-sample="{sample}" data-segment="{segment}" '
            f'data-coverage="{cov_frac:.6f}">' + "".join(cells) + "</tr>"
        )

    return (
        '<table class="stats-table" id="stats-table">'
        f"<thead><tr>{''.join(header_cells)}</tr></thead>"
        f'<tbody>{"".join(body_rows)}</tbody></table>'
    )


# --------------------------------------------------------------------------- #
# Figure builders
# --------------------------------------------------------------------------- #
def _log_depth_range(max_value: float) -> list:
    """Log10 axis range for a depth chart, generalizable across run scales.

    Floor is 1 (10^0) so a depth of 0/1 sits on the axis bottom and the 20x/100x
    guides are always inside the range; the top is 1.5x the data max but never
    below 150, so the 100x guide never hugs the ceiling on shallow runs.
    """
    top = max(150.0, float(max_value) * 1.5)
    return [0.0, math.log10(top)]


def throughput_series(df: pd.DataFrame, paired: bool) -> pd.DataFrame:
    """Decompose per-sample reads into the 3 stacked throughput series (spec §4).

    Reconciled to individual-read units: on paired-end (Illumina) runs the total
    and QC-passed counts are read *pairs*, so both are doubled; mapped reads are
    already counted individually (and pre-summed across segments by
    :func:`dedupe_and_sum_reads`). Nanopore is single-end — no doubling, and
    ``removed`` is 0 because there is no QC step (total == QC-passed).

    Returns the input frame with ``_total``/``_mapped``/``_qc_unmapped``/
    ``_removed`` columns, sorted by total reads ascending.
    """
    factor = 2 if paired else 1
    total = pd.to_numeric(df["number_of_reads"], errors="coerce").fillna(0.0) * factor
    qc = pd.to_numeric(df["number_of_trim_paired_reads"], errors="coerce").fillna(0.0) * factor
    mapped = pd.to_numeric(df["number_of_mapped_reads"], errors="coerce").fillna(0.0).clip(lower=0)
    # mapped is a subset of QC-passed, which is a subset of total; clamp the
    # deltas so rounding/quirks never draw a negative segment.
    qc_unmapped = (qc - mapped).clip(lower=0)
    removed = (total - qc).clip(lower=0)
    out = df.assign(_total=total, _mapped=mapped, _qc_unmapped=qc_unmapped, _removed=removed)
    return out.sort_values("_total", ascending=True, kind="stable")


def build_reads_histogram(
    df: pd.DataFrame, paired: bool = True, segmented: bool = False
) -> go.Figure:
    """Sequencing throughput per sample as horizontal 3-series stacked bars (spec §4).

    One row per sample (``df`` is the :func:`dedupe_and_sum_reads` output), plot
    height growing with sample count so 96 samples scroll instead of crushing into
    hairlines. The three stacked series — **Mapped reads** (accent), **QC-passed,
    unmapped** (green), **Removed by QC** (neutral grey) — sum to the sample's
    total reads (unit-reconciled; see :func:`throughput_series`). Bars sort by
    total reads ascending. The Absolute ⇄ Percent toggle (client-side ``barnorm``)
    normalises each sample to 100% for a fair QC-loss comparison across depths.
    """
    work = throughput_series(df, paired)
    samples = work["sample_name"].tolist()
    n = len(samples)
    height = max(200, 46 + n * (18 if segmented else 14))
    fig = go.Figure()
    fig.add_bar(
        y=samples,
        x=work["_mapped"],
        name="Mapped reads",
        orientation="h",
        marker_color=PALETTE_LIGHT[0],
        hovertemplate="%{y}<br>Mapped: %{x:,.0f}<extra></extra>",
    )
    fig.add_bar(
        y=samples,
        x=work["_qc_unmapped"],
        name="QC-passed, unmapped",
        orientation="h",
        marker_color=PALETTE_LIGHT[1],
        hovertemplate="%{y}<br>QC-passed, unmapped: %{x:,.0f}<extra></extra>",
    )
    fig.add_bar(
        y=samples,
        x=work["_removed"],
        name="Removed by QC",
        orientation="h",
        marker_color=EMPHASIS_MUTED,
        hovertemplate="%{y}<br>Removed by QC: %{x:,.0f}<extra></extra>",
    )
    fig.update_layout(
        autosize=True,
        height=height,
        barmode="stack",
        bargap=0.25,
        margin=dict(l=60, r=30, t=54, b=44),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family=FONT_FAMILY, size=13, color="#52514e"),
        # Fixed-width legend entries spaced out so the three labels don't crowd.
        legend=dict(
            orientation="h", yanchor="bottom", y=1.01, xanchor="left", x=0, traceorder="normal"
        ),
        xaxis=dict(title_text="Reads", gridcolor=GRID_COLOR, zeroline=False, rangemode="tozero"),
        yaxis=dict(
            type="category",
            gridcolor="rgba(0,0,0,0)",
            zeroline=False,
            automargin=True,
            title_text="",
        ),
    )
    return fig


# --------------------------------------------------------------------------- #
# Coverage heatmap (spec §5) — the centrepiece that replaces N overlaid lines
# --------------------------------------------------------------------------- #
HEATMAP_BINS = 150


def _bin_depths(positions, depths, genome_len: int, n_bins: int = HEATMAP_BINS) -> List:
    """Mean depth per genome bin (spec §5A), aligned to a fixed ``[1, genome_len]``.

    Bins with no sampled position become ``None`` (rendered as a gap, honest about
    missing data). Fed the already-downsampled cache arrays, so the means inherit
    the min-pool bias that keeps dropouts visible — the right emphasis for a
    dropout-at-a-glance overview.
    """
    if genome_len <= 0 or len(positions) == 0:
        return [None] * n_bins
    positions = np.asarray(positions)
    depths = np.asarray(depths, dtype=float)
    edges = np.linspace(1, genome_len + 1, n_bins + 1)
    idx = np.clip(np.digitize(positions, edges) - 1, 0, n_bins - 1)
    sums = np.zeros(n_bins)
    counts = np.zeros(n_bins)
    np.add.at(sums, idx, depths)
    np.add.at(counts, idx, 1)
    return [float(sums[i] / counts[i]) if counts[i] > 0 else None for i in range(n_bins)]


def _bin_centers(genome_len: int, n_bins: int = HEATMAP_BINS) -> List[float]:
    edges = np.linspace(1, genome_len + 1, n_bins + 1)
    return [float((edges[i] + edges[i + 1]) / 2) for i in range(n_bins)]


def _log10_or_none(v) -> Optional[float]:
    """log10 of a depth for the log heatmap; ``None`` for 0/None (honest gap)."""
    if v is None or v <= 0:
        return None
    return float(math.log10(v))


def sample_overall_coverage(
    df: pd.DataFrame, segmented: bool, lengths: Dict[Tuple[str, Optional[str]], int]
) -> Dict[str, float]:
    """Per-sample whole-genome horizontal coverage (for worst-first ordering).

    Unsegmented: the row's ``horizontal_coverage``. Segmented: length-weighted
    across segments (same weighting as the KPI tiles), so the ordering matches the
    headline figures. Falls back to an equal-weight mean when no track lengths are
    available (a pruned output dir).
    """
    samples = _ordered_unique(df["sample_name"].tolist())
    if not segmented:
        by = df.drop_duplicates("sample_name").set_index("sample_name")
        return {s: _safe_float(by.loc[s, "horizontal_coverage"]) for s in samples}
    segments = _ordered_unique(df["segment"].tolist())
    seg_len = {seg: max((lengths.get((s, seg), 0) for s in samples), default=0) for seg in segments}
    if sum(seg_len.values()) == 0:
        seg_len = {seg: 1 for seg in segments}
    genome = sum(seg_len.values()) or 1
    hc = {
        (r["sample_name"], r["segment"]): _safe_float(r["horizontal_coverage"])
        for _, r in df.iterrows()
    }
    return {
        s: sum(hc.get((s, seg), 0.0) * seg_len[seg] for seg in segments) / genome for s in samples
    }


def worst_first_order(overall: Dict[str, float], samples: List[str]) -> List[str]:
    """Sample order, worst overall coverage first (stable on ties)."""
    return sorted(samples, key=lambda s: (overall.get(s, 0.0), samples.index(s)))


def build_heatmap_model(
    df: pd.DataFrame,
    cache: Dict[Tuple[str, Optional[str]], Tuple[np.ndarray, np.ndarray]],
    lengths: Dict[Tuple[str, Optional[str]], int],
    segmented: bool,
    order: List[str],
) -> dict:
    """Assemble the coverage-heatmap model (spec §5), embedded as JSON for the client.

    ``position`` holds one binned depth grid per segment key (``""`` unsegmented):
    ``x`` bin centres, ``zNatural`` and ``zLog`` matrices (samples × bins, rows in
    ``order``). ``grid`` (segmented only) is the per-segment horizontal-coverage %
    matrix (samples × segments) for the "All segments" mode. The client renders
    natural/log depth (mode A) or the coverage grid (mode B) from these arrays —
    plain lists, so they dodge Plotly's typed-array serialisation.
    """
    segments = _ordered_unique(df["segment"].tolist()) if segmented else [None]

    position: Dict[str, dict] = {}
    for seg in segments:
        genome_len = max((lengths.get((s, seg), 0) for s in order), default=0)
        x = _bin_centers(genome_len) if genome_len > 0 else []
        z_nat, z_log = [], []
        for s in order:
            pos, depth = cache.get((s, seg), (np.array([]), np.array([])))
            binned = _bin_depths(pos, depth, genome_len)
            z_nat.append(binned)
            z_log.append([_log10_or_none(v) for v in binned])
        position[seg or ""] = {"x": x, "zNatural": z_nat, "zLog": z_log, "genomeLen": genome_len}

    grid = None
    if segmented:
        seg_labels = [s for s in segments]
        hc = {
            (r["sample_name"], r["segment"]): float(r["horizontal_coverage"])
            for _, r in df.iterrows()
        }
        z = [[hc.get((s, seg), None) for seg in seg_labels] for s in order]
        # Coverage as a percentage (0-100) so the fixed colour scale reads in %.
        z = [[(v * 100.0 if v is not None else None) for v in row] for row in z]
        grid = {"segments": [s or "" for s in seg_labels], "z": z}

    return {"segmented": segmented, "samples": order, "position": position, "grid": grid}


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


def _kpi_block(
    samples: List[str],
    breadth: Dict[str, float],
    depth: Dict[str, float],
    params: ReportParams,
) -> dict:
    """One KPI card set (spec §2): sample count, #>=pass, #<warn, median coverage, mean depth.

    ``pass_count`` / ``below_warn`` are computed from the tweakable thresholds so
    the tile labels can restate them without hardcoding 90/70.
    """
    depths = [depth[s] for s in samples if s in depth]
    breadths = [breadth[s] for s in samples if s in breadth]
    return {
        "samples": len(samples),
        "pass_count": sum(1 for s in samples if breadth.get(s, 0.0) >= params.pass_threshold),
        "below_warn": sum(1 for s in samples if breadth.get(s, 0.0) < params.warn_threshold),
        "median_coverage": float(np.median(breadths)) if breadths else 0.0,
        "mean_depth": float(np.mean(depths)) if depths else 0.0,
    }


def _kpi_display(block: dict) -> dict:
    """Formatted strings for one KPI block (server-rendered tile fallback)."""
    samples = block["samples"]
    pass_pct = (block["pass_count"] / samples * 100) if samples else 0.0
    return {
        "samples": _fmt_int(block["samples"]),
        "pass_count": _fmt_int(block["pass_count"]),
        "pass_sub": f"{pass_pct:.0f}% of run",
        "below_warn": _fmt_int(block["below_warn"]),
        "below_warn_nonzero": block["below_warn"] > 0,
        "median_coverage": _fmt_pct(block["median_coverage"]),
        "mean_depth": _fmt_depth(block["mean_depth"]) + "×",  # × suffix
    }


def _kpi_samples_subtitle(
    segments: List[Optional[str]], segmented: bool, metadata: Optional[dict]
) -> str:
    """Subtitle under the 'Samples analyzed' tile (spec §2).

    Segmented: "N segments each". Amplicon: the primer-scheme basename if the run
    declared one, else the generic "amplicon".
    """
    if segmented:
        n = len([s for s in segments if s is not None]) or len(segments)
        return f"{n} segments each"
    scheme = (metadata or {}).get("primer_scheme")
    if scheme:
        return os.path.basename(str(scheme))
    return "amplicon"


def build_kpi_summary(
    df: pd.DataFrame,
    lengths: Dict[Tuple[str, Optional[str]], int],
    segmented: bool,
    params: ReportParams,
) -> dict:
    """Global (and, for segmented runs, per-segment) KPI figures for the top tiles.

    Reports the sample count, the number of samples at/above the pass threshold,
    the number below the warn threshold, the median horizontal coverage, and the
    run-level mean depth. For segmented runs, a sample's breadth and depth are
    **length-weighted across segments** (weight = each segment's coverage-track
    length), so both read as whole-genome figures; a sample that drops a segment
    is correctly penalised rather than having the segment ignored.

    Returns ``{"global": {...}, "per_segment": {seg: {...}}}``.
    """
    samples = _ordered_unique(df["sample_name"].tolist())

    if not segmented:
        by = df.drop_duplicates("sample_name").set_index("sample_name")
        breadth = {s: _safe_float(by.loc[s, "horizontal_coverage"]) for s in samples}
        depth = {s: _safe_float(by.loc[s, "average_depth"]) for s in samples}
        return {"global": _kpi_block(samples, breadth, depth, params), "per_segment": {}}

    segments = _ordered_unique(df["segment"].tolist())
    # Canonical per-segment length = longest track seen for that segment.
    seg_len = {seg: max((lengths.get((s, seg), 0) for s in samples), default=0) for seg in segments}
    if sum(seg_len.values()) == 0:
        # No coverage tracks available for any segment (e.g. the CLI was pointed
        # at an output dir whose per-segment coverage files were pruned) while the
        # stats CSV still carries valid per-segment breadth/depth. Fall back to an
        # equal-weight mean across segments so the whole-genome KPIs reflect the
        # CSV rather than collapsing to 0%/0x.
        seg_len = {seg: 1 for seg in segments}
    genome_len = sum(seg_len.values()) or 1
    hc = {
        (r["sample_name"], r["segment"]): _safe_float(r["horizontal_coverage"])
        for _, r in df.iterrows()
    }
    dp = {
        (r["sample_name"], r["segment"]): _safe_float(r["average_depth"]) for _, r in df.iterrows()
    }

    g_breadth, g_depth = {}, {}
    for s in samples:
        g_breadth[s] = sum(hc.get((s, seg), 0.0) * seg_len[seg] for seg in segments) / genome_len
        g_depth[s] = sum(dp.get((s, seg), 0.0) * seg_len[seg] for seg in segments) / genome_len

    per_segment = {}
    for seg in segments:
        s_with = [s for s in samples if (s, seg) in hc]
        per_segment[seg] = _kpi_block(
            s_with,
            {s: hc[(s, seg)] for s in s_with},
            {s: dp[(s, seg)] for s in s_with},
            params,
        )
    return {
        "global": _kpi_block(samples, g_breadth, g_depth, params),
        "per_segment": per_segment,
    }


def _load_coverage_cache(output_dir: str, df: pd.DataFrame, segmented: bool) -> Tuple[
    Dict[Tuple[str, Optional[str]], Tuple[np.ndarray, np.ndarray]],
    Dict[Tuple[str, Optional[str]], int],
    Dict[Optional[str], str],
]:
    """Load + downsample every (sample, segment) coverage track exactly once.

    Returns ``(cache, lengths, contig_by_segment)``: ``cache`` maps each key to
    the downsampled ``(positions, depths)`` for plotting; ``lengths`` maps it to
    the true pre-downsample track length (row count), used for length-weighting
    the whole-genome KPI figures on segmented runs; ``contig_by_segment`` maps
    each segment to its coverage contig name (column 1 of the coverage table),
    the join key for aligning annotation features.
    """
    cache: Dict[Tuple[str, Optional[str]], Tuple[np.ndarray, np.ndarray]] = {}
    lengths: Dict[Tuple[str, Optional[str]], int] = {}
    contig_by_segment: Dict[Optional[str], str] = {}
    for _, row in df.iterrows():
        sample = row["sample_name"]
        segment = row["segment"] if segmented else None
        table = load_basewise_table(resolve_basewise_path(output_dir, sample, segment))
        lengths[(sample, segment)] = len(table)
        if table.empty:
            cache[(sample, segment)] = (np.array([]), np.array([]))
            continue
        if segment not in contig_by_segment:
            contig_by_segment[segment] = str(table["reference_id"].iloc[0])
        pos, depth = downsample_min_pool(table["position"].to_numpy(), table["depth"].to_numpy())
        cache[(sample, segment)] = (pos, depth)
    return cache, lengths, contig_by_segment


def _fig_to_html(fig: go.Figure) -> str:
    # default_width="100%" makes the wrapper div fill the card (no fixed 900px
    # that would force sideways scrolling on narrow viewports); responsive:true
    # keeps the figure sized to that div as it changes.
    return pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        default_width="100%",
        config={"responsive": True},
    )


def render_report(
    output_dir: str,
    metadata: Optional[dict] = None,
    config: Optional[dict] = None,
    params: Optional[ReportParams] = None,
    fetch_annotation: bool = False,
) -> str:
    """Build the full self-contained HTML report string for a consensus run.

    ``metadata`` is the run-metadata dict from :func:`build_report_metadata`; when
    omitted it is inferred from the output dir (the CLI-on-an-old-dir path).
    ``config`` is the full run config YAML (as a dict); when given it powers the
    run-parameters panel and supplies annotation source paths, otherwise both are
    omitted. ``params`` carries the tweakable thresholds/colours (defaults
    preserve historical behaviour). ``fetch_annotation`` allows an NCBI lookup for
    a missing gene GFF3 (off in the workflow, on for ``viralunity report``).
    """
    if metadata is None:
        metadata = build_report_metadata(config, output_dir)
    if params is None:
        params = ReportParams()
    # On Illumina (paired-end), total/QC-passed counts are read pairs while mapped
    # reads are counted individually; the mapping-rate denominator is doubled to
    # reconcile the units. Driven by the declared library layout, never by ratios.
    paired = metadata["library_layout"] == "paired"
    stats_path = os.path.join(output_dir, "assembly", "assembly_stats_summary.csv")
    df = read_stats_summary(stats_path)
    segmented = is_segmented(df)
    samples = _ordered_unique(df["sample_name"].tolist())
    segments = _ordered_unique(df["segment"].tolist()) if segmented else [None]

    cache, coverage_lengths, contig_by_segment = _load_coverage_cache(output_dir, df, segmented)
    kpi_summary = build_kpi_summary(df, coverage_lengths, segmented, params)
    kpi_samples_sub = _kpi_samples_subtitle(segments, segmented, metadata)
    annotation_model = build_annotation_model(
        output_dir, segments, contig_by_segment, config, fetch_annotation
    )

    # Worst-coverage-first sample order, shared by the heatmap y-axis and the
    # by-sample list so both surface problem samples first.
    overall_coverage = sample_overall_coverage(df, segmented, coverage_lengths)
    order = worst_first_order(overall_coverage, samples)
    # Server-rendered by-sample list rows (worst-first), each with a status dot +
    # coverage badge so a reviewer can jump straight to a bad sample.
    sample_list = [
        {
            "sample": s,
            "tier": params.tier(overall_coverage.get(s, 0.0)),
            "pct": _fmt_pct(overall_coverage.get(s, 0.0)),
        }
        for s in order
    ]

    # Global figures.
    stats_table_html = build_stats_table_html(df, paired, params)
    reads_fig_html = _fig_to_html(
        build_reads_histogram(dedupe_and_sum_reads(df), paired, segmented)
    )
    heatmap_model = build_heatmap_model(df, cache, coverage_lengths, segmented, order)

    # Per-sample coverage data, embedded once as JSON for the lazy detail section.
    coverage_json: Dict[str, list] = {}
    for sample in samples:
        entries = []
        for segment in segments:
            pos, depth = cache.get((sample, segment), (np.array([]), np.array([])))
            if len(pos) == 0:
                continue
            # Compute both axis ranges once, here, and hand them to the client so
            # the JS scale toggle consumes them as data instead of re-deriving them
            # (kills the Python<->JS _log_depth_range duplication and gives the
            # linear axis an explicit data-fit range rather than autorange).
            md = float(np.max(depth)) if len(depth) else 0.0
            entries.append(
                {
                    "label": segment or "",
                    "x": [int(p) for p in pos],
                    "y": [float(d) for d in depth],
                    "logRange": _log_depth_range(md),
                    "linRange": [0.0, md * 1.08 if md > 0 else 1.0],
                }
            )
        coverage_json[sample] = entries

    # The by-sample panel is only meaningful when at least one coverage track was
    # found (a fully pruned dir has none).
    has_coverage = any(entries for entries in coverage_json.values())

    env = Environment(
        loader=FileSystemLoader(TEMPLATE_DIR),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template(TEMPLATE_NAME)
    return template.render(
        plotly_js=get_plotlyjs(),
        metadata=metadata,
        kpi_global=_kpi_display(kpi_summary["global"]),
        kpi_json=_json_for_script(kpi_summary),
        kpi_samples_sub=kpi_samples_sub,
        config_panel_html=build_config_panel_html(config),
        segmented=segmented,
        samples=samples,
        segments=[s or "" for s in segments],
        has_coverage=has_coverage,
        stats_table_html=stats_table_html,
        reads_fig_html=reads_fig_html,
        heatmap_json=_json_for_script(heatmap_model),
        sample_list=sample_list,
        coverage_json=_json_for_script(coverage_json),
        annotation_json=_json_for_script(annotation_model["by_segment"]),
        palette_light=PALETTE_LIGHT,
        palette_dark=PALETTE_DARK,
        track_genes=f"{TRACK_GENE_COLOR},{TRACK_GENE_COLOR_DARK}",
        track_primers=f"{TRACK_PRIMER_COLOR},{TRACK_PRIMER_COLOR_DARK}",
        params=params,
        params_json=_json_for_script(params.as_client_dict()),
    )


def write_report(
    output_dir: str,
    dest: str,
    metadata: Optional[dict] = None,
    config: Optional[dict] = None,
    params: Optional[ReportParams] = None,
    fetch_annotation: bool = False,
) -> None:
    """Render the report and write it to ``dest`` (the shared entry point)."""
    html_str = render_report(output_dir, metadata, config, params, fetch_annotation)
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    with open(dest, "w") as fh:
        fh.write(html_str)
    logger.info("Wrote consensus report: %s", dest)


# --------------------------------------------------------------------------- #
# Entry points
# --------------------------------------------------------------------------- #
def report_params_from_config(config: Optional[dict], **overrides) -> ReportParams:
    """Build :class:`ReportParams` from optional config keys + explicit overrides.

    The pipeline config may carry ``report_pass_threshold`` /
    ``report_warn_threshold`` / ``report_chart_color`` /
    ``report_colorbar_thickness`` (all optional); explicit ``overrides`` (from CLI
    flags) win over config, which wins over the dataclass defaults. ``None``
    overrides are ignored so an unset CLI flag falls through to config/default.
    """
    config = config or {}
    fields = {
        "pass_threshold": config.get("report_pass_threshold"),
        "warn_threshold": config.get("report_warn_threshold"),
        "chart_color": config.get("report_chart_color"),
        "colorbar_thickness": config.get("report_colorbar_thickness"),
    }
    for key, value in overrides.items():
        if value is not None:
            fields[key] = value
    return ReportParams(**{k: v for k, v in fields.items() if v is not None})


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
    ap.add_argument(
        "--pass-threshold",
        dest="pass_threshold",
        type=float,
        default=None,
        help="Coverage fraction for the green/pass tier (default 0.90). Drives the "
        "'>=90%% coverage' KPI and the status dots/bars.",
    )
    ap.add_argument(
        "--warn-threshold",
        dest="warn_threshold",
        type=float,
        default=None,
        help="Coverage fraction for the amber/warn tier (default 0.70). Drives the "
        "'Below 70%%' KPI and the 'low coverage only' table filter.",
    )
    ap.add_argument(
        "--chart-color",
        dest="chart_color",
        default=None,
        help="Accent hex colour for the heatmap scale, by-sample line, and mapped-reads bar.",
    )
    ap.add_argument(
        "--colorbar-thickness",
        dest="colorbar_thickness",
        type=int,
        default=None,
        help="Heatmap colour-bar thickness in px (default 14).",
    )
    ap.add_argument(
        "--fetch-annotation",
        dest="fetch_annotation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="When a run has no gene annotation staged or in its config, fetch it "
        "from NCBI by the reference accession (default on). --no-fetch-annotation "
        "keeps the report fully offline.",
    )
    args = ap.parse_args()
    dest = args.output_path or os.path.join(args.input_dir, "report.html")
    config = load_run_config(args.config_file)
    metadata = build_report_metadata(config, args.input_dir)
    params = report_params_from_config(
        config,
        pass_threshold=args.pass_threshold,
        warn_threshold=args.warn_threshold,
        chart_color=args.chart_color,
        colorbar_thickness=args.colorbar_thickness,
    )
    write_report(args.input_dir, dest, metadata, config, params, args.fetch_annotation)


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
    config = dict(snakemake.config)  # noqa: F821
    metadata = build_report_metadata(config, output_dir)
    params = report_params_from_config(config)
    write_report(output_dir, dest, metadata, config, params)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
