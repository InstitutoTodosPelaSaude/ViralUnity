# Task: implement the reviewed report UX / data-viz improvements in the ViralUnity generator

You are working in the ViralUnity repository. Your job is to move a reviewed set of
UX and data-visualization improvements into the **HTML consensus report generator**
so that every report ViralUnity produces has them — not as a client-side patch, but
baked into the source that builds the report.

Everything you need is in the **`report-ux-improvements/`** folder at the repo root:
- `README.md` — orientation, the report's inferred architecture, and the **gotchas**.
- `RECOMMENDATIONS.md` — the full findings (issue → why → source-side fix).
- `patch_reports.py` — the **reference implementation** (the desired end state as
  small, commented vanilla-JS/CSS functions).
- `reference-reports/*.html` — the **patched "after" reports** to diff against.

Read all four before writing any code.

## Guiding principle: fix at the right layer

The reference implementation applies fixes **client-side at page load** because the
reviewer only had rendered HTML. **You have the source, so fix at the source
layer** — that is almost always cleaner and permanent:

- **Number formatting, the mapping-rate column, "honest zeros", linear-default
  axes, removing the average-depth chart, renaming sections, the throughput
  chart's % panel** → do these in the **Python / Plotly / template** code that
  builds the figure or table. Do *not* ship a JS post-processor that rewrites the
  DOM after load.
- **Genuinely interactive behavior** — the Linear/Log10 toggles, the collapsed
  `<details>` accordions with lazy rendering, the segment/sample selectors — stays
  as emitted client-side JS. Port those from `patch_reports.py` as real template
  JS, keeping the report a **single self-contained HTML file** (no external assets).

Use `patch_reports.py` as the *spec* for the exact desired markup/behavior;
translate it idiomatically into the generator rather than copying it verbatim.

## Workflow

1. **Discover, don't guess.** Find where reports are generated: the Python that
   computes stats and builds the Plotly figures, the HTML/Jinja template, and any
   client-side JS. Map the three layers (README describes what to expect). Identify
   the `_log_depth_range` duplication. Summarize what you found before proceeding.
2. **Plan and confirm scope.** Propose an ordered plan (see priority below) and
   confirm it with the maintainer before large changes. Prefer many small,
   reviewable commits over one big one.
3. **Add the enabler first:** a run/library-metadata field
   (`platform`, `library_layout`, `primer_scheme`, `qc_performed`) threaded into
   the report context. Several later items depend on it.
4. **Implement in priority order** (below), one change per commit.
5. **Verify every change end-to-end.** After each change, **regenerate the reports
   from a small fixture dataset** and compare against `reference-reports/` — open
   in a browser, check the specific behavior, and confirm nothing else regressed.
   A change isn't done until you've seen the regenerated report do the right thing.
6. **Lock it in with tests** (see below) so these classes of issue can't silently
   return.

## Priority order

1. **Metadata field** (enabler).
2. **P1 — the reader hits these first:** thousands separators (locale-pinned);
   depth charts linear-by-default with a Log10 toggle and data-fit range; honest
   zero-coverage; reconcile read units to one scale via the library-layout field.
3. **P2:** mapping-rate % (table + throughput chart's lower panel); rename
   "Reads per sample" → "Sequencing throughput"; distinct hues + legend; responsive
   charts; keyboard/AT-accessible sortable headers; progressive-disclosure
   accordions ("By sample" + "By segment", collapsed, lazy).
4. **P3:** remove the average-depth chart; dark-mode palette derivation +
   validation; spacing/dividers; drop redundant embedded plot titles.

(Full detail and rationale for each is in `RECOMMENDATIONS.md`.)

## Respect the gotchas (from README — these cost real debugging time)

- **Pin the number-format locale**; keep display vs sort key separate.
- **Never flip a Plotly axis type with `relayout`** — rebuild with `react`/`newPlot`
  (the `relayout` path froze the tab ~6 s). This matters wherever you emit a
  scale toggle.
- **Threshold labels belong at the true value (20, 100), not `10^value`**; set an
  explicit data-fit range, don't trust autorange with off-scale annotations.
- **Infer nothing about library layout from read-count ratios** — use the field.
- **Compute `_log_depth_range` once** (Python), pass to the client.
- **Never clamp depth `0 → 1`** to satisfy a log axis.

## Tests & conventions to add

- A **single number-formatting helper** in the generator, used everywhere.
- A **visual-regression test** from a tiny committed fixture (one Illumina, one
  Nanopore, one segmented sample with a known zero-coverage run and a known
  paired ×2 case): snapshot the rendered HTML + a headless screenshot of each
  chart, **plus a toggle-latency assertion** (guards the `relayout` freeze).
- **Palette validation** in CI (light + dark) for color-blind separation + contrast.
- A short **`REPORT.md`** documenting the conventions (linear-default depth; zeros
  as zero/gap; read counts in one unit via the metadata field; fixed validated
  palette order) so contributors follow them.

## Quality bar

- Match the repo's existing code style, structure, and naming.
- Keep each report a **single self-contained HTML file** — no external JS/CSS/font
  requests (inline everything).
- Don't regress existing outputs for platforms not covered by the four sample
  reports; think about single-sample runs, many-sample plates, and both segmented
  and unsegmented genomes.
- Preserve accessibility (keyboard operation, `aria-sort`, non-color-only encoding)
  and both light/dark themes on anything you touch.

When in doubt about intended behavior, open the corresponding `reference-reports/`
file — it is the ground truth for the desired end state.
