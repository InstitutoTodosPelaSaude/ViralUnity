# Report UX / data-viz improvements — handoff kit

This folder is a self-contained handoff for improving the **HTML consensus reports**
that ViralUnity generates. A UX / data-visualization review was done against four
rendered reports (Illumina paired-end, Nanopore single-end, and two segmented
genomes — Influenza with 8 segments, Guaroa with 3). Every fix was first
implemented **client-side on the rendered HTML** as a working reference; the job
now is to move those fixes into the **report generator** in this repo so they are
produced on every run.

## What's in this folder

| File | What it is |
|---|---|
| **`PROMPT.md`** | The task prompt. Start here — it tells you how to do the work well. |
| **`RECOMMENDATIONS.md`** | The detailed findings: each issue → why it matters for QC → the source-side fix, with real numbers. |
| **`patch_reports.py`** | The **reference implementation**. A self-contained Python script that makes two string substitutions (a CSS block + the ~5.6 KB custom JS tail) to the byte-identical report templates. **Every fix is a small, commented vanilla-JS function in here** — the clearest thing to translate into the generator. (Reference only; its absolute paths point at the reviewer's machine — do not run it here.) |
| **`reference-reports/`** | The four **patched "after" reports**. Open them in a browser to see the target; diff a freshly generated report against them. |

## How the report appears to be built (inferred from output)

Three layers — locate each in this repo first:

1. **Python stats + Plotly figures** → embedded JSON (reads-per-sample,
   aggregated-coverage, etc.).
2. **HTML template** → shell, `<style>`, the assembly-statistics table, a
   `data-palette` on `<body>`.
3. **Client-side JS** → lazy per-sample coverage plots (`coverageLayout` /
   `showSample`), theming, table sort. A comment references **`_log_depth_range`**,
   which is **duplicated in Python and JS** (a bug source — see gotcha 5).

The four report templates are **byte-identical** apart from embedded data.

## The single highest-leverage change

The report has **no field for platform / library layout (paired vs single-end) /
primer scheme / whether QC ran**. That gap forces the mapping-rate math to be a
heuristic and blocks platform-aware behavior. Add that metadata to the report
context first — several other items collapse into "read the field."

## Gotchas that will bite (reproduce before you doubt them)

1. **Locale-pinned number formatting.** A bare `toLocaleString()` / `{:n}` renders
   `46341` as `46.341` in de-DE etc. — reads as a *decimal*. Use `en-US` /
   `f"{n:,}"`. Keep the display string and the sort key separate.
2. **Never flip a Plotly axis type with `relayout`.**
   `Plotly.relayout(gd, {'yaxis.type': …})` on the ~2000-point coverage line
   figures froze the tab **~6 s per toggle**; `Plotly.react`/`newPlot` with a fresh
   layout is **~15 ms**. Rebuild the figure; don't mutate its axis type. Add a
   toggle-latency check to CI.
3. **Threshold-annotation bug.** The baked depth charts place "20×"/"100×" labels
   at `y = 10^value` (`1e20`, `1e100`). Harmless on a fixed log range, but any
   **autorange** (e.g. switching to linear) pulls `1e100` into the range and the
   profile collapses to a flat line. Place them at the true value; set an explicit
   data-fit range.
4. **Don't infer library layout from the numbers.** "`Total > QC` ⇒ paired" is
   correct for the current pipeline but misclassifies an unfiltered-paired or
   filtered-single-end run. Drive it from the declared library-layout field.
5. **Kill the Python↔JS `_log_depth_range` duplication** — compute once in Python,
   pass to the client.
6. **Honest zeros** — never clamp depth `0 → 1` for a log axis; plot 0 on linear,
   break the line (`null`) on log.

See `RECOMMENDATIONS.md` for the full reasoning behind each, and `PROMPT.md` for
how to execute.
