# Empirical before → after analysis

Acceptance evidence for the report UX / data-viz pass, produced by **generating real
reports** and extracting checkable signals from them — not by eyeballing. It renders the
committed fixtures under `test/fixtures/report/` with both the pre-UX generator (git
`6277c03`) and the current generator, then diffs concrete signals.

Reproduce (also writes the rendered before/after HTML for browser comparison):

```bash
python report-ux-improvements/empirical_analysis.py --out /tmp/report_ux
# open /tmp/report_ux/before_unsegmented.html vs after_unsegmented.html, etc.
```

## Signal table

Rendered from the two shared fixtures (before + after) plus the new Nanopore fixture
(after-only). `True`/`False` are presence of each signal in the generated HTML.

| signal | before unseg | after unseg | before seg | after seg | after nanopore |
|---|:--:|:--:|:--:|:--:|:--:|
| counts grouped with thousands separators | ✗ | ✓ | ✗ | ✓ | ✓ |
| mapped shown as a rate (`Mapped %`) | ✗ | ✓ | ✗ | ✓ | ✓ |
| baked "Depth (log)" axis title | ✓ | ✗ | ✓ | ✗ | ✗ |
| standalone average-depth chart | ✓ | ✗ | ✓ | ✗ | ✗ |
| card titled "Sequencing throughput" | ✗ | ✓ | ✗ | ✓ | ✓ |
| honest-zero gap (`connectgaps:false`) | ✗ | ✓ | ✗ | ✓ | ✓ |
| scale toggle rebuilds via `Plotly.react` | ✗ | ✓ | ✗ | ✓ | ✓ |
| `yaxis.type` relayout flip (the freeze) | ✗ | ✗ | ✗ | ✗ | ✗ |
| "By sample" accordion | ✗ | ✓ | ✗ | ✓ | ✓ |
| "By segment" accordion | ✗ | ✓ | ✗ | ✓ | ✓ |
| fixed 900px chart width | ✓ | ✗ | ✓ | ✗ | ✗ |
| min depth in embedded coverage | 4.0 | 4.0 | 4.0 | 4.0 | **0.0** |

## What the numbers show, finding by finding

- **Thousands separators (F1).** Pre-UX reports print raw integers; after, every count is
  grouped (`10000` → `10,000`, `462000` → `462,000`). Sorting is unaffected — the raw value
  stays in `data-sort`.
- **Mapping rate (F3b) + unit reconciliation (F3/F6).** The `Mapped %` column is computed
  from the declared library layout. Extracted cells from the generated tables:
  - Illumina paired, unsegmented: `sample-A` = **47.4 %** (`9000 / (2 × 9500)`), `sample-B` =
    **33.3 %** (`1200 / (2 × 1800)`). The ×2 denominator is what makes the funnel honest.
  - Nanopore single-end: `barcode05` = **55.5 %** (`256410 / 462000`, no doubling) — matching
    the value the review measured on the real Nanopore run.
- **Linear-default depth (F2) + no log title (F2d).** The baked "Depth (log)" axis title is
  gone; depth is linear by default with a client-side Log10 toggle.
- **Honest zeros (F4).** The Nanopore fixture carries a deliberate depth-0 run at positions
  1–5; it survives into the embedded coverage data (**min depth 0.0**, not clamped to 1), and
  the log path breaks the line at zeros (`connectgaps:false`) instead of inventing a depth.
- **No relayout freeze (F2f).** The scale toggle rebuilds via `Plotly.react`; no report emits
  a `yaxis.type` relayout — the guard against the ~6 s tab freeze.
- **Average-depth chart removed (F3d)** and **reads chart renamed to "Sequencing
  throughput" (F3c)**, both confirmed present/absent as expected.
- **Responsive (F5).** The fixed 900px chart wrapper is gone (`default_width="100%"` +
  `responsive:true`), so charts reflow into the card.
- **Progressive disclosure (F5c/F5d).** Both "By sample" and "By segment" accordions are built
  client-side from the embedded per-sample/per-segment data.

## Notes

- The two shared fixtures have no zero-coverage positions (min depth 4.0), so the honest-zero
  win is demonstrated by the purpose-built Nanopore fixture (min depth 0.0 preserved).
- `report size` is ~4.9 MB in both before and after — dominated by the vendored Plotly bundle,
  unchanged by these edits; the report stays a single self-contained file.
