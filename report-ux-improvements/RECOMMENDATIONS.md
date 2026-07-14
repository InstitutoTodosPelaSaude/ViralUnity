# ViralUnity consensus report — visual QA review and development recommendations

**Scope of review:** the four rendered report previews (`real_illumina.html`,
`real_nanopore.html`, `real_influenza.html`, `real_guaroa.html`), covering
Illumina paired-end, Nanopore single-end, and segmented genomes (Influenza 8
segments, Guaroa 3). All four share one byte-identical template (CSS block +
~5.6 KB custom JS tail); only the embedded data and the pre-baked Plotly JSON
differ.

**What this document is.** We only had rendered HTML, not the pipeline source
that generates it. The four report files in this folder have been patched
client-side as a *working reference implementation* of the fixes; this document
translates each fix into the equivalent **source-side change** for the report
generator, and then proposes process changes so these classes of issue stop
recurring. Wherever a fix was applied to the rendered files, it is marked
**[applied to previews]**.

A companion interactive page (built during review) shows current-vs-proposed
treatments on real data and is the design reference for the deeper redesigns
called out below.

---

## How the report is built (inferred)

From the rendered output, the generator appears to have three layers:

1. **A Python step** that computes statistics and builds the aggregate figures
   (reads-per-sample, average-depth, aggregated-coverage) with Plotly, then
   serializes them to JSON embedded in the page.
2. **An HTML template** (the shared shell: header, cards, `<style>`, the stats
   table, and a `data-palette` attribute on `<body>`).
3. **Client-side JS** that renders the per-sample coverage plots lazily
   (`coverageLayout` / `showSample`) and handles theming and table sorting.

Two observations drive several recommendations below:

- A JS comment reads *"must match Python's `_log_depth_range`."* The log-axis
  range logic is **duplicated** in Python and JS. Duplicated logic that must be
  kept in sync by hand is a standing bug source.
- The reads-per-sample chart already carries the caption *"Top panel: read-pair
  counts. Bottom panel: individual read counts."* — so the maintainers already
  **know** the two panels are in different units. The problem is that they are
  presented on two separate axes rather than reconciled to one, which is what
  makes them hard to compare (see Finding 3).

---

## Findings and fixes

Priority: **P1** = misleads or blocks reading the data · **P2** = clarity /
accessibility · **P3** = polish / robustness.

### Finding 1 — Numbers have no thousands separators · P1 · [applied to previews]

Every count in both the main statistics table and the per-sample table is
printed as a raw integer (`348018`, `5256777`, `21440`). At 6–7 digits with no
grouping these are genuinely hard to read and compare at a glance.

**Why it matters.** These tables are the first thing a reviewer scans to decide
whether a sample passed. Misreading `4381806` vs `438186` is a real QC error,
not a cosmetic one.

**Fix applied to previews.** A single formatter groups digits with a comma and
is applied to the main table on load and inside `renderSampleStats`; the raw
value is preserved in the `data-sort` attribute so column sorting is unaffected.

**Source-side fix.** Format at generation time in **one** place. Two gotchas we
hit and that the source must handle:

- **Pin the locale.** A bare `toLocaleString()` / `{:n}` renders `46341` as
  `46.341` in some locales — which reads as a *decimal*. Use an explicit
  grouping locale (`toLocaleString('en-US')`, or Python
  `f"{n:,}"`), never the ambient one.
- **Keep display and sort separate.** The visible string is grouped; the sort
  key stays numeric. (The template already uses `data-sort` for exactly this —
  keep feeding it the raw number.)
- Do not group the fractional part of mean-depth values, and leave
  non-numbers (sample IDs, `S1`, `99.6%`) untouched.

### Finding 2 — "Depth" charts default to a log axis with no visual cue · P1 · [applied to previews]

Three charts put depth on a log axis: per-sample coverage (client-side) and,
labelled "(log scale)", the average-depth bar chart and the aggregated-coverage
profile. Plotly's default log ticks print as plain numbers (1, 10, 100, 1000),
so unless the reader notices the uneven tick spacing the axis reads as linear
and every value is silently misjudged by orders of magnitude.

The **average-depth chart is the worst case: it is a bar chart on a log axis.**
Bar length is only an honest magnitude cue against a linear, zero-anchored
baseline. On a log axis a 40 % real difference (3,245× vs 2,319×) reads as a
~6 % height difference — the primary visual channel actively misleads.

**Why it matters.** Depth and its uniformity are core QC signals. A scale that
distorts them undermines the report's main job.

**Decision taken (with the maintainer during review):** default both depth
charts to **linear**, with an explicit **Log10 toggle** for the genuinely
wide-dynamic-range cases (e.g. a low-on-target segmented/metagenomic library
where the mapped bar would otherwise vanish). Linear-by-default also removes the
need for the zero-clamping hack in Finding 4.

**Fix applied to previews.** The client-side coverage chart is rebuilt
linear-first with a per-card Linear/Log10 toggle. The aggregated-coverage profile
is flipped to linear with an explicit **data-driven range** (`[0, max × 1.08]`,
not raw autorange — see Finding 2c) plus the same toggle, and the axis
titles/heading no longer hard-code "(log)"/"(log scale)". (The average-depth chart
that this originally also fixed was later **removed entirely** — Finding 3d.)

**Source-side fix.** Generate these charts **linear by default**, with the
y-range fit to the data. If you keep a log option, add real log tick formatting
(`10⁰, 10¹, 10²` or `1×/10×/100×`) so the scale is unmistakable, and label the
axis "Depth (log₁₀)" only when log is actually active. For a one-value-per-sample
aggregate, prefer a **dot/lollipop plot** over log-scaled bars if the range is
wide — dots don't rely on baseline-anchored length.

### Finding 2b — Average-depth chart isn't segment-aware ~~· P2 · [applied to previews]~~ · SUPERSEDED by Finding 3d (chart removed)

*Originally: for segmented genomes the baked average-depth bar chart plotted every
sample×segment combination as one flat row of bars, so you couldn't compare
samples within a segment. This was first fixed by rebuilding it segment-aware —
but the chart was subsequently **removed entirely** (Finding 3d) as low-value,
since the mean-depth numbers live in the stats tables. Retained here only to
record the reasoning.*

**Source-side fix (if the chart is ever reinstated).** Emit it faceted or
filterable **by segment** (x = sample, one panel/selection per segment), from the same
per-segment table the report already computes.

### Finding 2c — A latent annotation bug blows up the axis range when linear · P1 · [applied to previews]

The baked depth charts place their "20×"/"100×" threshold *labels* at
`y = 10^20` and `y = 10^100` (the code appears to raise 10 to the threshold
value instead of using the value). On the original fixed **log** range these were
off-screen and unnoticed — but any **autorange** (as when switching to linear)
pulls `1e100` into the range calculation and the entire profile collapses to a
flat line near zero. This is the "y-axis range far too wide / depth looks near
zero" symptom.

**Fix applied to previews.** Threshold annotations are re-anchored to their true
depth (20 and 100), and the range is set **explicitly from the data in both
linear and log modes** (never autorange, which the bogus annotation poisons
before the re-anchor lands).

**Source-side fix.** Place the threshold annotations at `y = 20` and `y = 100`
(not `10^20` / `10^100`). This is a one-line bug in the figure builder; fixing it
at the source removes the need for the runtime re-anchor.

### Finding 2d — Log-axis y-caption and tick labels · P2 · [applied to previews]

Two readability problems on the log depth axes, both echoing the original
"overlaying natural and log scales" report: (a) the y-axis title repeated the
scale ("Depth (log₁₀)") even though a Linear/Log10 **toggle** sits right beside
the chart, and (b) Plotly's default log axis labels minor ticks as bare `2` / `5`
(meaning 2,000 / 5,000), easily misread as linear values.

**Fix applied to previews.** The y-axis caption is just **"Depth"** ("Mean depth"
for the bar chart) in both modes — the toggle states the scale — and log ticks
are pinned to **powers of ten** (`dtick = 1`: 1, 10, 100, 1k, 10k) so there are
no ambiguous minor labels.

**Source-side fix.** Same: don't encode the scale in the axis title when a scale
control is present, and set `dtick=1` (or an explicit `10⁰,10¹,…` tick format) on
log axes.

### Finding 2e — Aggregated-coverage: redundant embedded title + per-sample legend that won't scale · P2 · [applied to previews]

The aggregated-coverage figure carried an embedded plot title ("Aggregated
coverage — S") that duplicates the card heading and the Segment selector, and a
**legend/direct-labelled line per sample**. The per-sample legend is fine for two
samples but does not scale — a run of 20–50 samples would produce an unreadable
tangle of near-identical lines and a legend taller than the plot.

**Fix applied to previews.** The embedded title and the legend are removed from
the aggregated chart (top margin reclaimed). The per-sample "Coverage by sample"
plots keep their titles, because for segmented genomes those stack one plot per
segment and the title carries the segment identity (not redundant there).

**Source-side fix.** Don't emit a plot title that repeats the section heading.
More importantly, **rethink the many-sample aggregated view**: overlaying N raw
per-sample lines is a form that breaks down past a handful of samples. For large
N, show a **summary band** — median depth with an interquartile (or min–max)
ribbon across samples — so the "typical coverage and its spread along the genome"
reads at a glance regardless of sample count, with individual lines available on
demand (hover / a drill-down), not all at once.

### Finding 2f — Toggling a chart's scale via `Plotly.relayout` froze the tab (~6 s) · P1 · [applied to previews]

**This one is worth carrying into the real implementation.** Flipping an axis
between log and linear with `Plotly.relayout(gd, {'yaxis.type': …})` on the
~2,000-point aggregated-coverage line figures hit a Plotly slow path that took
**~6 seconds per toggle** and made the browser report the page as unresponsive.
It is not a data-volume problem: the identical 2,000 points render in **~13 ms**
via `Plotly.react` / `Plotly.newPlot`. The slow path is specific to an
incremental `relayout` that changes axis type (and, to a lesser degree, any
`relayout` while the axis is linear); a full re-render side-steps it entirely.

**Fix applied to previews.** The aggregated-coverage scale toggle now rebuilds the
figure with `Plotly.react(gd, gd.data, layout)` instead of mutating it with
`relayout`. Toggling is now 7–70 ms (single panel to 8-panel segmented) with no
freeze. The per-sample coverage plots were already safe because they re-plot with
`Plotly.newPlot` on every change; the average-depth chart uses `Plotly.react`.

**Source-side fix.** For any interactive control that changes an axis type or
otherwise re-scales a Plotly figure (a Dash callback, an ipywidgets handler, a
hand-rolled toggle), prefer **`react`/`newPlot` with a freshly built layout over
`relayout` that flips `type`**, especially on line traces with more than a few
hundred points. If a build ever adds such a control, include a **toggle-latency
check** in the visual-regression suite (Finding "process #3") so a multi-second
freeze can't ship unnoticed.

### Finding 3 — Reads-per-sample: units split across two axes; identity by opacity only · P1/P2 · [partly applied]

Two separate issues in one chart:

- **(P2, applied)** "Total reads" and "QC-passed reads" were the *same* blue,
  distinguished only by `opacity: 0.4` vs `1.0`, with no legend — identity was
  carried entirely by a small subplot annotation. This fails for color-blind
  readers and anyone screenshotting the panel. *Fix applied:* distinct hues
  (blue / green / orange), full opacity, and a legend.

- **(P1, source-side)** Total/QC-passed are plotted on one axis (read **pairs**,
  ~348k) and Mapped on a second axis (individual **reads**, ~643k). The report
  captions this but never reconciles it, so the three quantities that form one
  funnel — sequenced → QC-passed → mapped — cannot be compared on one scale.

**Verified on the real data.** Doubling QC-passed lands just above Mapped for
both Illumina samples (330,568 × 2 = 661,136 vs 643,219 mapped → ~2.7 % unmapped;
248,936 × 2 = 497,872 vs 466,724 → ~6.3 %). This confirms Total/QC-passed are
fragment (pair) counts and Mapped is an individual-read count. Nanopore is
single-end and its Total already equals QC-passed exactly — no doubling applies.

**Source-side fix.** Put all three quantities in **one unit on one axis**,
driven by library layout (see Finding 6): for paired-end, multiply
Total/QC-passed by 2 (or divide Mapped by 2 — pick one and label it) so the
funnel is directly readable. The design reference for this is a **three-dot
lollipop per sample** (one row per sample, dots for Total / QC-passed / Mapped,
the connecting line = reads lost at each step). It stays readable at 50+ samples
(one row each) and, for segmented viruses, a **segment switch** re-points only
the Mapped dot (Total/QC-passed are sample-level and hold still). This redesign
needs the per-sample/per-segment data the generator already has and belongs in
source, not a runtime patch on baked JSON — see the companion page for the
worked interaction.

**Important caveat for the implementer.** The ×2 correction must be driven by a
**declared library-layout field**, never inferred from the mapped/total ratio.
It reconciles cleanly here only because these are high-on-target amplicon runs; a
metagenomic library legitimately shows a low mapped ratio for unrelated reasons,
which is indistinguishable by ratio alone from "forgot to double."

### Finding 3b — Mapped reads shown as a raw count, not a rate · P2 · [applied to previews]

A raw mapped-read count (`643,219`, `33,112`, `426`) is hard to interpret without
mentally dividing by the input. The useful QC signal is the **mapping rate** —
what fraction of QC-passed reads aligned.

**Fix applied to previews.** The "Mapped reads" column (main table and per-sample
table) now shows a **percentage**, computed on the fly as
`mapped / QC-passed × 100`, made unit-consistent by doubling the QC-passed
denominator when the run is paired-end. Paired vs single-end is detected per row
by `Total > QC-passed` (ViralUnity's Illumina path filters, so Total > QC; its
Nanopore path runs no filter, so Total == QC). The raw count is preserved in the
cell's hover tooltip, and the column sorts by rate. Real results: Illumina
`sample-4117` → **97.3 %**; Nanopore `barcode05` (single-end) → **55.5 %**;
Influenza S1 (paired, off-target-heavy) → **0.36 %**. The **lower panel of the
"Sequencing throughput" chart** (Finding 3c) was also switched from an absolute
mapped-read count to this same mapping-rate percentage (0–100 axis), so the chart
and table agree.

**Source-side fix.** Compute the mapping rate at generation time from the
**declared library layout** (Finding 6), not the `Total > QC` heuristic — the
heuristic is correct for the current pipeline but would misclassify an unfiltered
paired run or a filtered single-end run. Show the rate as the primary value and
keep the absolute count available (tooltip or a secondary column).

### Finding 3c — "Reads per sample" naming and the mapped panel · P2 · [applied to previews]

Two small clarity items on the reads chart: the card was titled "Reads per
sample" (describes the x-axis, not what the chart is for), and its lower panel
plotted **absolute mapped reads** on a second read-count axis — inconsistent with
the table's mapping-rate column and hard to compare across samples.

**Fix applied to previews.** The card is retitled **"Sequencing throughput"**, and
the lower panel now shows **% mapped** on a fixed 0–100 axis (per Finding 3b), so
chart and table tell the same story. Series get distinct hues + a legend (from
Finding 3's P2 half). The vertical gap between the two stacked subplots was also
widened (subplot domains 0.14 → 0.22 apart) so the panels read as clearly
separate. (Note: the chart's mapped % is *sample-level* — total mapped across
segments ÷ QC-passed — whereas the table shows per-segment %; same formula,
different granularity, matching each view's axis.)

**Source-side fix.** Rename the section, emit the mapped panel as a rate, and set
a comfortable `vertical_spacing` on the two-row subplot.

### Finding 3d — "Average depth" chart removed · P3 · [applied to previews]

The stand-alone "Average depth" bar chart restated a single number already present
in every stats table (the `Mean depth` column) and, per the maintainer, added
little. It has been **removed** from the previews. The mean-depth values remain in
the assembly / per-sample / per-segment tables.

**Source-side fix.** Drop the average-depth figure from the generator; keep the
`Mean depth` column.

### Finding 4 — Zero-coverage positions are clamped to depth 1 · P1 · [applied to previews]

The client-side coverage code did `y.map(v => Math.max(1, v))` so the log axis
could render zeros. This collapses "no coverage at all" and "covered at depth 1"
into the same point — erasing exactly the dropouts a coverage plot exists to
show. (Real example: `sample-4117` has true zero runs at positions 1–46 and
29,837–29,890.)

**Fix applied to previews.** Linear mode (now default) plots raw depth, so zeros
sit honestly on the baseline. Log mode breaks the line at true zeros
(`null` + `connectgaps:false`) so no-coverage reads as a gap, never as a false
depth of 1.

**Source-side fix.** Never substitute a fake value to satisfy a log axis. Plot
zeros as zero (linear) or as an explicit gap / separate "no coverage" marker
(log). This also removes the duplicated `_log_depth_range` clamp logic.

### Finding 5 — Fixed 900 px chart width forces horizontal scrolling · P2 · [applied to previews]

The baked charts hard-code `width: 900px` (and the client-side coverage plots
used a fixed `FIG_WIDTH`), so on any viewport narrower than 900 px the card
scrolls sideways and content (including the legend's third item) is clipped.

**Fix applied to previews.** Charts are made responsive (`autosize`, wrapper
width neutralized to `max-width` + `100%`, `responsive: true`).

**Source-side fix.** Emit Plotly figures with `autosize` / `responsive` and CSS
width, not a pixel width. For the segmented "Coverage by sample" and
"Aggregated coverage" sections, lay panels out in a responsive grid
(`repeat(auto-fit, minmax(260px, 1fr))`) instead of a single fixed-width column
— an 8-segment sample currently produces ~3,360 px of vertical scroll before any
cross-segment comparison is possible.

### Finding 5b — "Coverage by sample" spacing: table→plot and segment→segment · P2 · [applied to previews]

Two cramped seams in the "Coverage by sample" card: the per-sample stats table
butted directly against the first coverage plot's title, and (for segmented
viruses) each per-segment plot's x-axis label ran straight into the next plot's
title with **zero gap** — hard to tell where one segment ended and the next began.

**Fix applied to previews.** A 28 px gap now separates the stats table from the
first plot (relying on the table's own last-row hairline, no extra line), and a
32 px gap + theme-aware hairline divider (`1px var(--border)`) separates
consecutive segment plots (the first has none), so each segment reads as its own
panel.

**Source-side fix.** Whatever renders the per-segment stack should emit the gap +
divider (or, better, the responsive grid from Finding 5). Trivial once the panels
are real layout children rather than bare appended `<div>`s.

### Finding 5c — Progressive disclosure: per-sample detail behind an accordion · P2 · [applied to previews]

A report reads best when run-level information leads and per-sample detail is
opt-in: most readers want "did the run work?" answered at a glance, and only some
drill into a specific sample's coverage. Previously all per-sample coverage plots
rendered inline and always, adding length (and render cost) most readers didn't
ask for.

**Fix applied to previews.** The per-sample card is now a native
`<details>`/`<summary>` accordion titled **"By sample"**, **collapsed by
default**, under a light-gray **"Details"** eyebrow — so the run-level cards
(assembly stats, sequencing throughput, aggregated coverage) are what you see
first. `<details>` is keyboard-operable and screen-reader-friendly with no
JavaScript. The per-sample plots are **lazy-rendered on first expand** — if a
reader never opens the section, that Plotly work never runs.

**Source-side fix.** Emit the per-sample section inside `<details>` (collapsed),
and defer building/serializing its figures until expand (a small `toggle`-event
handler, as done here). This also trims the initial payload and time-to-first-paint
for runs with many samples — for a 96-sample plate, deferring per-sample coverage
is a large, free win.

### Finding 5d — A second "By segment" view · P3 · [applied to previews]

The per-sample view answers "how did this sample do?"; the complementary
grouping — all samples' coverage for one segment (or, for an unsegmented genome,
simply all samples at once) — had no view.

**Fix applied to previews.** A second accordion, **"By segment"**, sits beside
"By sample" under the "Details" eyebrow, with the same mechanics (collapsed, lazy,
Linear/Log10 toggle, hairline dividers). Its behaviour adapts to the genome:

- **Segmented** (Influenza, Guaroa): a **Segment selector** picks a segment and
  shows its coverage across **every sample**, plus a focused stats row per sample.
  Surfaces cross-sample outliers directly — e.g. Influenza segment S3 mapped
  0.50 % in one sample vs 9.4 % in the other, side by side.
- **Unsegmented** (Illumina, Nanopore): one whole-genome "segment", so **no
  selector** — it just stacks every sample's coverage for quick all-samples
  comparison. The stats table is omitted there (it would duplicate Assembly
  statistics).

**Source-side fix.** Emit both groupings from the same per-sample/per-segment data
the generator already has (`COVERAGE[sample][segment]`, `STATS[sample][segment]`).
Only render the per-segment view when more than one segment exists.

### Finding 6 — No library-layout / platform / kit metadata in the report · P2 (enabler)

The report carries no field recording sequencing platform, library layout
(paired vs single end), or primer scheme. That single missing field is what
forces Findings 3 to be a guess and blocks correct per-platform behavior
(the QC step isn't even run for Nanopore, yet nothing in the report says so).

**Source-side fix.** Thread a small run-metadata block from the pipeline inputs
into the report context: `platform`, `library_layout`, `primer_scheme`,
`qc_performed`. Drive the ×2 read-unit reconciliation, the "no QC step run"
note, and any platform-specific captions from it. This is the highest-leverage
data change — several findings collapse into "read the metadata field."

### Finding 7 — Sortable table headers are not keyboard/screen-reader accessible · P2 · [applied to previews]

Sortable `<th>` used a bare `onclick` with no `role`, `tabindex`, `aria-sort`,
or key handling, so sorting was mouse-only and silent to assistive tech.

**Fix applied to previews.** Headers get `role="button"`, `tabindex="0"`,
Enter/Space activation, and a live `aria-sort` that updates on each sort.

**Source-side fix.** Emit those attributes in the template, and keep `aria-sort`
in sync inside the existing `sortTable`.

### Finding 8 — The categorical palette isn't validated, and fails in dark mode · P3

`<body data-palette="…">` declares 8 categorical hues, but only 2 (blue, green)
are drawn today, and the same hex values are reused for light and dark themes.
Running the six-check color procedure:

- **Light:** passes the lightness, chroma, and color-blind-separation checks;
  three slots fall below the 3:1 surface-contrast floor (mitigated only because
  every plotted value also appears in the stats table).
- **Dark:** **fails** — 4 of 8 slots fall outside the dark lightness band and one
  drops below contrast. The two hues actually in use pass in both modes, so
  there's no urgent visual break today, but the reserve slots are not safe to
  start using in dark mode as-is.

**Source-side fix.** Derive a **dark-mode-specific** set of steps (don't reuse
the light hexes), and run a palette validator in CI (see below) so any new
series color is checked before it ships. A quick lightness-only correction is
*not* enough — we tried it and it broke color-blind separation between two
adjacent slots; the fix needs the full check loop.

---

## Process recommendations (adjusting the development course)

These target the *causes*, so the findings above don't recur.

1. **One number-formatting helper, used everywhere.** Finding 1 exists because
   values are emitted ad hoc. Add a single `format_count()` / `format_depth()`
   in the generator (locale pinned, separator explicit) and route every rendered
   number through it. New columns then inherit correct formatting for free.

2. **Kill the Python↔JS logic duplication.** The `_log_depth_range` value (and
   any other quantity computed in both languages) should be computed **once** in
   Python and passed to the client as data, not re-implemented in JS "to match."
   Duplicated numeric logic silently drifts.

3. **A visual-regression / snapshot test on the template.** Generate reports
   from a tiny fixed fixture dataset (one Illumina, one Nanopore, one segmented
   sample with a known zero-coverage run and a known ×2 case) and snapshot the
   rendered HTML + a headless screenshot of each chart in CI. This is what would
   have caught the log-scale bar chart, the clamped zeros, and the missing
   thousands separators before release. Keep the fixture in the repo. **Include a
   toggle-latency assertion** for any interactive control (Finding 2f) — a chart
   interaction that takes seconds is a shipping defect a static snapshot misses.

4. **Prefer `react`/`newPlot` over type-flipping `relayout`.** Finding 2f showed a
   `relayout` that changes `yaxis.type` freezing the tab for ~6 s where a full
   re-render is ~13 ms. Make "rebuild the figure, don't mutate its axis type" the
   house rule for any scale/axis control.

5. **Palette validation in CI.** Feed `data-palette` (light and dark variants)
   through a color-blindness / contrast validator on every build; fail the build
   on a hard violation. Prevents Finding 8 from reaching users when the reserve
   slots start being used.

6. **An accessibility checklist in PR review** for anything touching the
   template: keyboard operation of every control, `aria-sort` on sortable
   headers, legend/label not color-only, and a color-blind check of chart hues.

7. **Document the conventions once, in the repo.** A short `REPORT.md`: depth
   axes are linear by default with an optional log toggle; zero coverage is shown
   as zero/gap and never clamped; read counts are reconciled to a single unit via
   the library-layout field; the categorical palette order is fixed and
   validated. New contributors then follow the convention instead of
   rediscovering these traps.

---

## Priority order for the pipeline backlog

1. **Finding 2c (annotation `10^value` bug)** — a one-line source fix that
   removes a latent axis-range blow-up; do it first, it's nearly free.
2. **Finding 6 (metadata field)** — unblocks 3 / 3b and several captions;
   smallest data change, largest downstream simplification.
3. **Finding 1 (number formatting)**, **Finding 2 (linear-default depth)**, and
   **Finding 3b (mapping rate)** — what the reader hits first and most often.
4. **Finding 4 (honest zeros)** and **Finding 3 (unit reconciliation + lollipop)**
   — correctness of the coverage and funnel stories.
5. **Findings 5, 7 (responsive + a11y)** — broaden who can read the report and
   on what device.
6. **Finding 8 (dark palette) + process items** — hardening so none of the
   above regress.

---

## Files touched in this review

- `real_illumina.html`, `real_nanopore.html`, `real_influenza.html`,
  `real_guaroa.html` — patched in place with the client-side reference fixes for
  every `[applied to previews]` finding above (number formatting; linear-default
  depth + Log10 toggle; honest zeros; aggregated-coverage range/title/legend fix;
  reads chart → "Sequencing throughput" with a % mapped panel; average-depth chart
  removed; responsive charts; a11y; "Details" section with **By sample** / **By
  segment** accordions, collapsed + lazy). Originals preserved as `*.html.bak`.
- **The reference implementation is a single, self-contained Python patcher** that
  makes exactly two string substitutions (a CSS block and the ~5.6 KB custom JS
  tail) — the four templates are byte-identical, so one patch covers all. That
  script is the clearest artifact to port from: every fix is a small, commented
  vanilla-JS function. (Ask the reviewer for `patch_reports.py`.)
- The full reads-chart redesign (three-dot lollipop + segment switch, Finding 3)
  is demonstrated on the companion interactive page, not on the baked charts.
- The four `.bak` files can be deleted once the source-side changes land; they
  are only a rollback safety net for the preview patch.
