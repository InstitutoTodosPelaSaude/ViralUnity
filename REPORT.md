# Consensus report conventions

The `viralunity report` / `generate_html_report` output is a single self-contained
interactive HTML file built by `viralunity/scripts/python/generate_consensus_report.py`
and `templates/report_template.html.j2`. These conventions exist so the data-visualization
mistakes fixed in the v1.3.4 UX pass don't quietly return. Follow them when touching the
report; the test suite (below) enforces most of them.

## Depth axes are linear by default
Coverage/depth charts render **linear** with the y-range fit to the data. A log axis makes
Plotly's default ticks (1, 10, 100) read as linear and silently misjudges every value by
orders of magnitude. A per-chart **Linear / Log10** toggle is available for genuinely
wide-dynamic-range cases.

- **Never flip an axis with `Plotly.relayout({'yaxis.type': …})`.** On the ~2000-point
  coverage lines that hits a Plotly slow path (~6 s freeze). Rebuild the figure with
  `Plotly.react` / `Plotly.newPlot` and a freshly built layout instead. The snapshot test
  guards this (asserts `Plotly.react`, forbids a quoted `yaxis.type` relayout key).
- Threshold guides (20×, 100×) sit at their **true value** (20, 100), never `10^value`, and
  the range is fit to the data — never a raw autorange that a stray annotation can poison.
- The Log10 axis range is computed **once in Python** (`_log_depth_range`) and passed to the
  client as `logRange` per coverage entry. Do not re-derive it in JS.

## Zero coverage is honest, never clamped
Never substitute a fake depth to satisfy a log axis. Plot 0 on the baseline (linear) or break
the line at true zeros (`null` + `connectgaps:false`) on log. Collapsing "no coverage" into
"depth 1" erases exactly the dropouts a coverage plot exists to show.

## Read counts are reconciled to one unit via declared metadata
Total/QC-passed are read **pairs** on paired-end (Illumina) runs; mapped reads are counted
individually. The mapping rate is `mapped / QC-passed`, with the denominator **doubled when
the run is paired-end**. Library layout comes from the declared run-metadata field
(`build_report_metadata`, from the config YAML's `data`/`scheme`), **never** inferred from the
mapped/total ratio. The table shows the rate as the primary value; the raw count stays in the
cell tooltip.

## One number-formatting helper
All counts route through the locale-pinned formatters (`_fmt_int` / `_fmt_depth`, Python's
`,` spec — always a comma, never a locale decimal point). The visible string is grouped; the
raw value stays in `data-sort` so column sorting is unaffected.

## Fixed, validated categorical palette
`PALETTE_LIGHT` / `PALETTE_DARK` are a fixed 8-hue order (never cycled); >8 samples fall back
to a single recessive hue plus hover so identity is never colour-alone. Dark mode uses a
**dark-specific derivation**, not the light hexes reused. `test/scripts/palette_validation_test.py`
checks colour-blind separation and dark-mode contrast in CI before any reserve slot ships.

## Self-contained, responsive, accessible
No external JS/CSS/font requests — inline everything (Plotly is vendored once). Figures are
responsive (`autosize` + `default_width="100%"`, CSS caps the width), so narrow viewports
reflow instead of scrolling sideways. Sortable headers are keyboard-operable with a live
`aria-sort`; encodings are never colour-only (legend + labels). Both light and dark themes
must keep working on anything you touch.

## Progressive disclosure
Run-level cards (assembly stats, sequencing throughput, aggregated coverage) lead. Per-sample
and per-segment detail live in collapsed, lazily rendered `<details>` accordions ("By sample",
"By segment") so a many-sample plate stays fast to open.

## Annotation tracks (genes + primers)
Optional gene (GFF3) and primer-scheme (BED) tracks are drawn as rectangles on a second y-axis
beneath the depth line, sharing the genome x-axis, on all three coverage views. They are Plotly
traces (not layout shapes) so they survive the Linear/Log10 rebuild and stay fixed while the depth
axis toggles. The files are staged into `<output>/annotation/` by the workflow (column 1 is
sanitized on Nanopore to match the sanitized coverage contig names), and the report matches features
to each segment's contig by exact name. Genes / Primers chips toggle each track; both default on and
appear only when that data is present. Track hues are low-chroma, sourced from Python, and kept
distinct from the depth-series blue (validated in `palette_validation_test.py`).

Two reporting choices are baked into the GFF3 parser (`parse_gff3`) and open to revision:
- **Feature type:** `gene` features are drawn, falling back to `CDS` on a contig with no genes.
- **Label:** the first present of `Name → gene → gene_name → locus_tag → product → ID`.

## Tests
- `test/scripts/generate_consensus_report_test.py` — pure helpers (metadata, formatting,
  mapping rate, figure builders).
- `test/report_generation_test.py` — end-to-end render against `test/fixtures/report/`.
- `test/report_snapshot_test.py` — HTML-content snapshot + the relayout-freeze JS guard.
- `test/scripts/palette_validation_test.py` — palette CVD/contrast validation.

Fixtures live under `test/fixtures/report/` (Illumina paired, Nanopore single-end with a known
zero-coverage run, segmented Influenza-like). Keep them tiny and committed.
