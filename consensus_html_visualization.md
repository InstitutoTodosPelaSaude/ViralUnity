# Implementation spec — interactive HTML consensus reports (v1.3.4)

**Audience:** an agent picking this up fresh, with no memory of how this spec was produced.
**Status:** design complete and reviewed; **no code written yet**. This document is the
authoritative source — it supersedes any older/partial notes on this topic.

**Repo:** ViralUnity. 

Read `CLAUDE.md` at the repo root first — it explains that the Python layer is a thin
CLI/orchestration shim and the substantive work lives in the Snakemake `.smk` files under
`viralunity/scripts/`. Follow its editing rules exactly (don't rename `validate_args` /
`generate_config_file` / `run_snakemake_workflow`, keep per-rule conda envs, edit `rule all` +
`organize_files` together when adding outputs, the versioning policy, etc.).

## Skills / workflow expectations (if your harness supports skills)

- Follow TDD: write failing tests for pure functions first, then implement.
- Run a verification pass before declaring done — actually generate reports from the two example
  output dirs and inspect them, don't just rely on unit tests.

---

## 1. Goal

Generate a single self-contained interactive HTML report visualizing the results of a ViralUnity
**consensus** run. Scope is the four `consensus_{illumina,nanopore}[_segmented].smk` entry points
only (not `meta`). The report must work for both unsegmented viruses (e.g. SARS-CoV-2) and
segmented viruses (e.g. Influenza A, 8 segments).

## 2. Branch & version

- Create and work on a new branch `feature/interactive-visualizations-v1.3.4`, branched off
  `feature/filters-v1.3.2` (the branch this design work happened on — check it exists and is what
  you expect before branching; if `feature/filters-v1.3.2` has since merged to `main`, branch from
  `main` instead and note the deviation).
- The version for this feature (all interactive visualizations) is **v1.3.4**.
- Per `CLAUDE.md`'s versioning policy this is arguably MINOR-sized work, but the version number
  for this feature was explicitly fixed at v1.3.4 by direct instruction — don't renumber it.

## 3. Two example output dirs to validate against

Gitignored working dirs at the repo root — use for manual validation only, **not** as test
fixtures (do not commit anything from these into `test/`):
- Unsegmented (SARS-CoV-2): `outputs/test_consensus_illumina/`
- Segmented (Influenza A, segments S1–S8): `outputs/test_consensus_influenza/`

Both were inspected directly and confirmed to contain real data (see §4).

## 4. Data sources — exact, verified formats

### 4.1 Global stats table

`<output>/assembly/assembly_stats_summary.csv`.

Unsegmented header (verified against `outputs/test_consensus_illumina/assembly/assembly_stats_summary.csv`):
```
sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage
sample-4117,348018,330568,643219,3245.2445239608064,0.9959535832525165,0.9959535832525165,0.9475972310470522,0.9959535832525165
```

Segmented header inserts `segment` as column 2 (one row per sample × segment; verified against
`outputs/test_consensus_influenza/assembly/assembly_stats_summary.csv`):
```
sample_name,segment,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage
sample-sample-1177,S1,5256777,4661046,33112,1893.439982913285,0.9794959419051688,0.9700982486117044,0.6569841947885519,0.9756514310123878
```

**Detect segmented mode by the presence of the `segment` column.**

**The `percentage_above_*x` and `horizontal_coverage` columns are fractions in [0,1], not already
×100.** The report must multiply by 100 and add a `%` suffix when displaying them; round
`average_depth` for display (e.g. to 1 decimal or nearest integer — your call, be consistent).

This file is produced by `viralunity/scripts/python/calculate_assembly_stats.py` (per-row) and
concatenated by a per-entry-point `rule unify_assembly_statistics_reports` (see §7).

### 4.2 Per-base coverage

`<output>/assembly/coverage_stats/{sample}.table_cov_basewise.txt` (unsegmented) or
`<output>/assembly/{segment}/coverage_stats/{sample}.table_cov_basewise.txt` (segmented).

Verified format: **whitespace-separated (tab), no header, 3 columns**: `reference_id  position  depth`.

Unsegmented example (`outputs/test_consensus_illumina/assembly/coverage_stats/sample-4117.table_cov_basewise.txt`, 29,903 rows):
```
MN908947.3	1	0
MN908947.3	2	0
```

Segmented example (`outputs/test_consensus_influenza/assembly/S1/coverage_stats/sample-sample-1177.table_cov_basewise.txt`, 2,341 rows):
```
NC_007373.1|Influenza_A|H3N2|segment_1	1	0
NC_007373.1|Influenza_A|H3N2|segment_1	2	0
```

Note: segment reference IDs can contain `|` characters — treat `reference_id` as an opaque string,
don't split on any assumed delimiter. Files can be ~30k rows (SARS-CoV-2 genome length) —
**downsample for plotting** (see §6.4).

Produced by `rule calculate_coverage_basewise` in `viralunity/scripts/rules/stats.smk` (shared by
all 4 entry points via `include:`), using `bedtools genomecov -d -ibam {input.bam}`. Output path
template: `config['output'] + "assembly/" + SEGMENT_WILDCARD + "coverage_stats/{sample}.table_cov_basewise.txt"`
where `SEGMENT_WILDCARD` is `""` (unsegmented) or `"{segment}/"` (segmented) — set per entry point.

### 4.3 Sample naming

Sample names appear **exactly** as in the pipeline's files — e.g. `sample-4117`, and in the
influenza fixture literally `sample-sample-1177` (a double `sample-` prefix, because the fixture's
own input sample id already started with `sample-`). **Do not strip, reformat, or "clean up" any
prefix.** Treat `sample_name` (from the CSV) and the `{sample}` in coverage filenames as opaque,
identical strings — they are the same wildcard value.

## 5. Report layout (top = global, bottom = per-sample)

Visual language borrows from https://dsystem.sinapse.org.br/ (the "Sinapse UI Theme Layer"): OKLCH
color tokens, semantic `background/foreground/card/border/primary` roles, light/dark mode, rounded
cards, clean typography, generous spacing, card-based sections, a light/dark toggle. (This is a
styling reference to imitate the *feel* of, not a resource to fetch from at runtime — the report
must be offline/self-contained, so hand-author CSS matching that aesthetic rather than linking to
the site.)

### 5.1 Global section

- **Stats table**: render `assembly_stats_summary.csv` as a clean, sortable table (client-side
  sort via a small inline `<script>`, not a second library). Format percentages as `%` (×100),
  round depths.
- **Reads histogram**: grouped bar chart per sample of total reads (`number_of_reads`),
  QC/trimmed reads (`number_of_trim_paired_reads`), and mapped reads (`number_of_mapped_reads`).
  In segmented mode, `number_of_reads`/`number_of_trim_paired_reads` are identical across a
  sample's segment rows — **dedupe to one value per sample** (take any row, e.g. `.first()`);
  `number_of_mapped_reads` is per-segment, so **sum it across segments** for this global view (or
  stack by segment — either satisfies the requirement, sum is simpler and was the design choice
  made; note if you deviate).
- **Coverage (depth) bar chart**: `average_depth` per sample (per sample × segment in segmented
  mode — no dedupe here, one bar per row of the summary CSV). **Log y-axis**; draw light
  horizontal guide lines at **y=20** and **y=100** (the "20x"/"100x" marks — the end-to-end test
  will grep the rendered HTML for the literal substrings `"20x"` and `"100x"`, so label the guide
  lines with that exact text, e.g. as an annotation).
- **Aggregated coverage line plot**: x = genome position, one line per sample, hover shows depth.
  Same 20x/100x horizontal guide lines. **In segmented mode, produce one aggregated plot per
  segment**, selectable via a dropdown/tabs (segments have different lengths, so they cannot share
  one x-axis).

### 5.2 Per-sample section (below global)

Must be viewable **one sample at a time** — do NOT stack all samples' plots on the page (it must
stay short/lightweight). Use a sample `<select>` dropdown. For the selected sample show: its stats
row, and its per-base coverage line plot (one plot per segment in segmented mode), with the same
20x/100x guide lines.

### 5.3 UI mechanism for dropdowns/tabs (decided — do not redesign)

**Pre-render every figure statically at HTML-build time (in Python, via Plotly); toggle visibility
with a small inline `<script>` driven by `<select onchange>`.** Do **not** use Plotly's
`updatemenus`/animation frames for this — those are built for switching traces/frames *within one
figure*, not for swapping between independently-built figures with different x-axis ranges
(segments have different lengths), and using them would fight the "only the selected item's plot
is active" requirement.

Concretely:
- Sample panels: `<div class="sample-panel" data-sample="{name}">` each containing that sample's
  coverage-line div(s) (one per segment if segmented). Only one panel has `display:block`; the
  rest `display:none` via CSS, flipped by the dropdown's `onchange` handler
  (`document.querySelectorAll('.sample-panel').forEach(p => p.style.display = 'none');
  document.querySelector('[data-sample="'+val+'"]').style.display = 'block';` — reuse the same
  generic function for the segment dropdown, parameterized by a CSS class name/data attribute).
- Segment panels for the *global aggregated* plot: same pattern, `.segment-panel[data-segment]`.
- **Every Plotly figure must set fixed pixel `width`/`height` in `update_layout` (not
  `autosize`).** This is what makes Plotly render correctly inside a `display:none` container at
  page-parse time — Plotly sizes from the layout config, not the container's computed box, so no
  `Plotly.Plots.resize()` call-on-reveal is needed. Skipping this causes the classic "chart renders
  as 0×0 inside a hidden div" bug.

## 6. Self-contained + charting requirements

- The HTML must be fully self-contained and open offline: **no CDN links**; all JS/CSS embedded
  inline.
- Use Plotly via `plotly.py`. Embed the Plotly library inline **exactly once**: build the document
  with a Jinja2 template, inline `plotly.min.js` once (via `plotly.offline.get_plotlyjs()`), and
  render each individual figure with `plotly.io.to_html(fig, include_plotlyjs=False,
  full_html=False)`. Do not let each figure carry its own copy of `plotly.js`.
- **Downsampling**: reduce per-base coverage to at most ~2000 points per line before plotting,
  using **min-pooling** within bins — i.e. split the position range into ~2000 equal-width bins
  and, per bin, keep the row (position AND depth together) with the **minimum depth**, not the
  mean. This preserves coverage dips instead of averaging them away. Function signature suggestion:
  `downsample_min_pool(positions: np.ndarray, depths: np.ndarray, max_points: int = 2000) ->
  tuple[np.ndarray, np.ndarray]`.

## 7. Snakemake structure — verified facts (read before editing any `.smk` file)

All four consensus entry points (`viralunity/scripts/consensus_illumina.smk`,
`consensus_illumina_segmented.smk`, `consensus_nanopore.smk`, `consensus_nanopore_segmented.smk`):

- Each defines `SEGMENT_WILDCARD` near the top: `""` for unsegmented, `"{segment}/"` for
  segmented; segmented variants additionally define `SEGMENTS = config["reference"]` (a dict:
  segment name → fasta path) and build `rule all` targets via `expand(..., segment=SEGMENTS.keys())`.
- Each `include:`s `rules/stats.smk`, which defines `rule calculate_coverage_basewise` (the
  per-base coverage rule, output template in §4.2) — shared, not duplicated.
- Each entry point **duplicates its own** `rule unify_assembly_statistics_reports` (writes
  `config['output'] + "assembly/assembly_stats_summary.csv"`, with the header hard-coded in a
  shell heredoc — segmented variants' header includes `segment`) and its own `rule organize_files`
  (the final target, output `config['output'] + "benchmark.tsv"`). These two rules are the closest
  existing precedent for the new report rule — **duplicate the new rule per entry point the same
  way**, rather than trying to factor it into one shared `rules/*.smk` file, because the `expand()`
  shape (with/without `segment=`) differs between segmented and unsegmented variants.
- `<output>` (the pipeline's top-level output dir, i.e. what the task spec calls `<output>`) is
  `config['output']`; `benchmark.tsv` sits directly at `<output>/benchmark.tsv`, alongside
  `<output>/samples/`, `<output>/assembly/`. The new report should be written to
  `<output>/report.html` — same directory as `benchmark.tsv`.
- Conditional `rule all` targets already exist, e.g. (both Illumina variants only — Nanopore has
  no isnv step):
  ```python
  config['output'] + "isnvs/isnvs_summary.tsv" if config.get("run_isnv", False) else [],
  ```
  This is the pattern to copy for the new `generate_html_report` toggle (§9).
- A representative `script:`-directive rule for calling convention reference:
  `rule calculate_assembly_statistics` in `viralunity/scripts/rules/consensus_illumina_common.smk`
  calls `viralunity/scripts/python/calculate_assembly_stats.py` via
  `main(snakemake.input, snakemake.output[0], snakemake.params[0], getattr(snakemake.wildcards,
  "segment", None))` — i.e. params/wildcards are readily accessible from the injected `snakemake`
  global.
- Conda envs with pandas already exist: `viralunity/scripts/envs/utils.yaml` (`pandas=2.0` +
  `seqtk`, `biopython`, `samtools`, `bcftools`) — but it doesn't have `plotly`/`jinja2`, and is used
  by many lightweight shell rules that shouldn't need a plotting stack, so **add a new env** rather
  than extending it (see §8).

### 7.1 New rule to add (duplicated in all 4 entry-point `.smk` files)

```python
rule generate_html_report:
    conda:
        "envs/report.yaml"
    input:
        stats_summary = rules.unify_assembly_statistics_reports.output.unified_stats_summary,
        basewise = expand(rules.calculate_coverage_basewise.output.table_cov, sample=config["samples"])
        # segmented variants: expand(..., sample=config["samples"], segment=SEGMENTS.keys())
    output:
        report = config['output'] + "report.html"
    params:
        output_dir = config['output']
    script:
        "python/generate_consensus_report.py"
```

`input.basewise` exists purely so Snakemake schedules/tracks staleness against the upstream
coverage rule — the script itself must reconstruct every basewise file path independently from
`params.output_dir` plus the rows of the stats CSV (via a `resolve_basewise_path` helper, §8), not
from `snakemake.input` directly. **This is the load-bearing design reason the same core module
serves both the Snakemake path and the CLI path unmodified** — the CLI path has no `snakemake`
object at all and must resolve paths the same way from a plain directory argument.

Also add to `rule all` in each entry point:
```python
config['output'] + "report.html" if config.get("generate_html_report", True) else [],
```

### 7.2 Resource lists — do NOT touch

`viralunity/constants.py`'s `ResourceDefaults.CONSENSUS_ILLUMINA_RULES` /
`CONSENSUS_NANOPORE_RULES` lists only contain genuinely heavy rules (`perform_qc`, `map_reads`,
`trim_primer_sequences`, `detect_isnv`, `infer_consensus_sequence`). Lightweight utility rules
(`organize_files`, `unify_assembly_statistics_reports`, `rename_sequences`) are deliberately absent
— they get no dedicated `{rule}_cpus`/`{rule}_ram` config keys. `generate_html_report` is
lightweight; follow the same precedent and do **not** add it to these lists.

## 8. Core module and file layout

### 8.1 `viralunity/scripts/python/generate_consensus_report.py` (new)

Single source of truth imported by **both** delivery paths. Follow the established dual-mode
pattern used by ~19 existing files in this same directory (e.g.
`viralunity/scripts/python/apply_max_rpm_bleed_filter.py`):
```python
def run_cli():
    ...  # argparse entry, used when invoked directly / from the CLI wrapper
def run_snakemake():
    ...  # reads snakemake.input / snakemake.output / snakemake.params
if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
```
Both `viralunity/scripts/` and `viralunity/scripts/python/` already have `__init__.py` files (this
was verified — `find viralunity/scripts -name __init__.py` returns both), so a normal Python import
`from viralunity.scripts.python.generate_consensus_report import write_report` works fine from the
new CLI module. **No CLI module in this repo does that today** (verified: no `scripts.python`
imports exist anywhere under `viralunity/*.py` as of this writing) — this is a deliberate, new but
low-risk pattern; call it out in your commit message / PR description as intentional.

Suggested pure functions (design them for independent unit-testability — each should take
in-memory data, not file paths, wherever the logic itself doesn't need I/O):

- `read_stats_summary(csv_path: str) -> pd.DataFrame`
- `is_segmented(df: pd.DataFrame) -> bool` — `"segment" in df.columns`
- `resolve_basewise_path(output_dir: str, sample: str, segment: str | None = None) -> str` —
  reconstructs the exact path from §4.2/§7.
- `load_basewise_table(path: str) -> pd.DataFrame` — columns `["reference_id", "position",
  "depth"]`, no header; if the file is missing, log a warning and return an empty frame rather than
  raising, so a partial output directory still renders for the samples that do have data.
- `downsample_min_pool(positions: np.ndarray, depths: np.ndarray, max_points: int = 2000) ->
  tuple[np.ndarray, np.ndarray]` — see §6.
- `dedupe_and_sum_reads(df: pd.DataFrame) -> pd.DataFrame` — per `sample_name`: `.first()` for
  `number_of_reads`/`number_of_trim_paired_reads`, `.sum()` for `number_of_mapped_reads`.
- `build_stats_table_html(df: pd.DataFrame) -> str`
- `build_reads_histogram(df: pd.DataFrame) -> plotly.graph_objects.Figure`
- `build_coverage_bar_chart(df: pd.DataFrame, segmented: bool) -> Figure` (log y-axis, `add_hline`
  at y=20 and y=100 with `"20x"`/`"100x"` annotations)
- `build_aggregated_coverage_line_plot(per_sample_series: dict, title: str) -> Figure` (one call
  per segment in segmented mode)
- `build_sample_detail_plot(sample: str, per_segment_series: dict) -> Figure`
- `render_report(output_dir: str) -> str` — orchestrates everything above, renders the Jinja2
  template, returns the full HTML string
- `write_report(output_dir: str, dest: str) -> None` — the function both delivery paths ultimately
  call

**Every figure must set fixed pixel `width`/`height`** (§5.3) — this is not optional, it's required
for the hide/show mechanism to work.

### 8.2 `viralunity/scripts/python/templates/report_template.html.j2` (new file + new `templates/` dir)

Single Jinja2 template. Contains:
- One `<script>{{ plotly_js }}</script>` block, populated once from `plotly.offline.get_plotlyjs()`.
- One generic inline `<script>` for client-side table sorting (stats table) and the
  sample/segment `<select>` toggle logic (§5.3) — reuse one small parameterized JS function for
  both dropdowns, don't duplicate it.
- Placeholders for each figure's `plotly.io.to_html(fig, include_plotlyjs=False,
  full_html=False)` output, and for the CSS implementing the card-based, light/dark-toggle look
  described in §5 (hand-authored, inline `<style>`, no external stylesheet).

### 8.3 `viralunity/scripts/envs/report.yaml` (new)

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - pandas=2.0
  - plotly
  - jinja2
```
(match the pinning style of `viralunity/scripts/envs/utils.yaml` for the exact version syntax used
elsewhere in this repo — check that file before finalizing pins.)

### 8.4 `viralunity/viralunity_report.py` (new) — standalone CLI subcommand

Model this on `viralunity/viralunity_create_samplesheet.py` — a **single file**, no
`_orchestrator.run_pipeline` (this isn't a pipeline launcher, it's a one-shot operation reading an
existing directory and writing one HTML file):
```python
@click.command("report")
@click.option("--input", "input_dir", required=True, help="Existing consensus output directory.")
@click.option("--output", "output_path", default=None, help="Destination HTML path (default: <input>/report.html).")
def report(input_dir, output_path):
    dest = output_path or os.path.join(input_dir, "report.html")
    try:
        write_report(input_dir, dest)
    except Exception as e:
        raise click.ClickException(str(e))
```
Import `write_report` from `viralunity.scripts.python.generate_consensus_report`.

Register in `viralunity/viralunity_cli.py`: add `from viralunity.viralunity_report import report`
near the other subcommand imports, and `cli.add_command(report)` alongside the existing
`cli.add_command(...)` lines (that file is a `click.group()` — read it, it's short, ~60 lines).

## 9. Config toggle — `generate_html_report` (default `True`)

Decided: this ships as an opt-out config toggle, not an always-on step (matches the existing
`run_isnv`-style pattern and lets users disable it on very large runs). Wire through all 6 touch
points per `CLAUDE.md`'s rule ("If you add a new pipeline option, you touch four places... plus the
`.smk` files"):

1. **`viralunity/constants.py`**: add `GENERATE_HTML_REPORT = "generate_html_report"` to the
   `ConfigKeys` class.
2. **`viralunity/config_generator.py`**: add `generate_html_report: bool = True` parameter to both
   `add_illumina_settings(...)` and `add_nanopore_settings(...)`, with
   `self._set(ConfigKeys.GENERATE_HTML_REPORT, generate_html_report, self.SECTION_PARAMETERS)` (or
   whatever the local section variable is named — check the existing `run_isnv` line right next to
   it for the exact idiom).
3. **`viralunity/viralunity_consensus_cli.py`**: add a two-state boolean flag —
   `click.option("--generate-html-report/--no-generate-html-report", default=True, help="Generate an interactive HTML report at the end of the run.")`
   — to whichever shared option group both the illumina and nanopore subcommands pull from (the
   file has a `_COMMON_OPTIONS`/`_add_common_options` decorator pattern shared across both
   subcommands — read the file, it's ~424 lines, and look at where `--create-config-only` is
   defined for the exact shared-options mechanism to hook into).
4. **`viralunity/viralunity_consensus.py`**: forward `generate_html_report=args.get("generate_html_report", True)`
   into the `ConfigGenerator` call, right next to the existing `run_isnv=args.get("run_isnv", False)`
   line.
5. **`viralunity/validators.py`**: no new validator function needed — this is a plain boolean with
   no cross-flag interaction to check (unlike e.g. `run_denovo_assembly`, which interacts with
   other metagenomics flags).
6. **The 4 `.smk` entry points**: `config.get("generate_html_report", True)` conditional in
   `rule all` (§7.1) — this is the actual gate; the config-generator/CLI touch points above just
   control whether that key ends up `False` in the generated YAML.

## 10. Packaging fix — required, not in the original ask

`pyproject.toml`'s `[tool.setuptools.package-data]` section currently reads:
```toml
[tool.setuptools.package-data]
"viralunity.scripts" = ["*.smk", "**/*.smk", "envs/*.yaml"]
```
This does **not** glob `.j2` template files. Add `"python/templates/*.j2"` to that list — otherwise
the new Jinja2 template silently goes missing from a non-editable `pip install` (breaks the
published/PyPI package even though local editable installs and tests won't notice, since they read
straight from the source tree). Don't skip this.

Also, per the original task: add `plotly` and `jinja2` directly to `[project.dependencies]` in
`pyproject.toml` (not a new optional extra) — the CLI path needs them unconditionally regardless of
whether Snakemake/conda is involved. Check the current dependency list first
(`pyyaml>=6.0, click>=8.0, biopython>=1.81, snakemake>=7.32,<8, pandas>=1.5, pulp<2.8` as of this
writing) and append to it.

## 11. Tests (TDD — write these first, then implement until green)

### 11.1 `test/scripts/generate_consensus_report_test.py` (new)

Follow the established convention in `test/scripts/*_test.py`: stdlib `unittest.TestCase`,
package-absolute imports (`from viralunity.scripts.python.generate_consensus_report import
is_segmented, dedupe_and_sum_reads, ...`), **inline** fixtures (small `pd.DataFrame`s built via
helper functions, `io.StringIO` for CSV text, small `numpy` arrays) — do not create checked-in
files for these, matching how e.g. `apply_max_rpm_bleed_filter_test.py` and
`add_contig_stats_to_summary_test.py` are structured.

Cover at minimum:
- `is_segmented`: `True`/`False` on inline DataFrames with/without a `segment` column.
- `dedupe_and_sum_reads`: a 2-segment synthetic DataFrame for one sample → assert
  `number_of_reads`/`number_of_trim_paired_reads` are deduped (not doubled) and
  `number_of_mapped_reads` is the sum of the two segment rows.
- `downsample_min_pool`: build a large synthetic `positions`/`depths` array (e.g. 20,000 points)
  with a known narrow, deep dip planted somewhere in the middle → assert (a) the output has
  ≤`max_points` points, and (b) the planted dip's minimum depth value is still present in the
  downsampled output (i.e. it wasn't averaged away).
- `resolve_basewise_path`: assert the unsegmented path (no `segment` arg) and segmented path (with
  `segment` arg) match the exact templates in §4.2/§7.
- `read_stats_summary` / `build_stats_table_html`: small inline CSV via `io.StringIO`, assert
  percentage columns render as `%`-formatted values (×100) and depth is rounded.

### 11.2 `test/report_generation_test.py` (new, end-to-end)

Uses small, **checked-in** fixture files under `test/fixtures/report/` — a deliberate deviation
from the inline-only convention in `test/scripts/`, justified because this test needs a realistic
*nested directory tree* (`assembly/coverage_stats/...` and `assembly/{segment}/coverage_stats/...`)
matching `resolve_basewise_path`'s exact layout; that's more legible checked in than assembled
file-by-file inside the test body via `tempfile.TemporaryDirectory()`. Create:
- A tiny `assembly_stats_summary.csv` (2–3 samples, unsegmented).
- 2–3 small `table_cov_basewise.txt` files (a few hundred rows each is plenty — they don't need to
  be realistic genome lengths, just enough rows to exercise the downsampler at a small scale, or
  simply not trigger it at all if under 2000 rows — both cases are fine to cover).
- A second small fixture set exercising the **segmented** case (a `segment` column in the CSV,
  and coverage files nested under a `{segment}/coverage_stats/` dir).

Assertions on the generated HTML:
- Output is a single file.
- Contains **no** `http://` or `https://` substring anywhere (grep the whole file).
- Contains each fixture sample name (verbatim, including any odd prefixes you chose in the
  fixtures — pick simple names like `sample-A`/`sample-B` to keep this assertion trivial).
- Contains the literal substrings `"20x"` and `"100x"`.
- Contains `"Plotly.newPlot"` (or equivalent Plotly-generated marker) at least once per expected
  figure — count how many figures your design produces for the fixture and assert that count, or
  at minimum assert it appears more than zero times.

### 11.3 Dryrun suite

**No `test/dryrun_configs/*.yaml` edits are required.** Once `rule all` in each `.smk` entry point
references `config.get("generate_html_report", True)` (default `True`), `report.html` becomes a
real DAG target automatically — existing dryrun YAMLs don't set this key, so they get the default
(on). `snakemake -n` (used by `test/viralunity_dryrun_test.py`) only resolves the DAG; it never
executes rule bodies, so it doesn't need the `report.yaml` conda env or the actual Python script to
work, just for the `script:` path referenced in the rule to exist on disk. Run
`pytest test/viralunity_dryrun_test.py -v` after the `.smk` edits to confirm all 4 consensus
workflows (plus their segmented variants) still dry-run clean.

If you want explicit coverage of the `generate_html_report: false` path too, add one new
`consensus_illumina__no_report.yaml`-style variant to `test/dryrun_configs/` — optional, not
required for the suite to pass.

### 11.4 Full local loop

```bash
conda activate viralunity   # or your env name; package installed via pip install -e .
python -m pytest test/ -q                          # or python -m unittest discover ./test -p *test.py, per CLAUDE.md
black --check viralunity/ test/ && ruff check viralunity/ test/
pytest test/viralunity_dryrun_test.py -v
```
If the dryrun suite fails with `AttributeError: module 'pulp' has no attribute 'list_solvers'`, run
`pip install 'pulp<3'` (Snakemake 7.32 is incompatible with newer pulp — known issue, documented in
`CLAUDE.md`).

## 12. Manual validation (required — do not skip)

```bash
viralunity report --input outputs/test_consensus_illumina --output /tmp/r1.html
viralunity report --input outputs/test_consensus_influenza --output /tmp/r2.html
```
Open both in a browser with no network (or just grep for zero `http(s)://` substrings) and
confirm, for each:
- Sortable global stats table with correctly formatted percentages.
- Reads histogram (grouped bars, one group per sample).
- Log-scale depth bar chart with visible 20x/100x guide lines.
- Aggregated coverage line plot — for the influenza fixture, confirm the segment dropdown/tabs
  actually switch between 8 distinct plots with different x-axis extents.
- Per-sample selector — confirm switching samples actually swaps the visible plot(s) (open dev
  tools, watch the DOM/`display` styles toggle) and that only one sample's plot(s) render at a
  time (page doesn't visibly contain all samples' charts rendered simultaneously).

Then run one real (small) consensus workflow end-to-end (illumina, unsegmented, `--create-config-only`
off) and confirm `<output>/report.html` is produced automatically by the Snakemake rule, and that
its content matches what the CLI produces for the same output directory.

## 13. Docs

- New `docs/report.md`: MyST/Sphinx style, **no YAML frontmatter** (matches `docs/output.md`'s
  style — it starts directly with `# Output Layout`; this new page should start directly with
  `# Interactive HTML Report`). Describe what each chart shows, and how to regenerate via
  `viralunity report --input <output_dir> --output report.html` without rerunning the pipeline.
- `docs/index.md`: add `report` to the `{toctree}` list (currently: `installation, tutorial/index,
  usage, commands, architecture, embedding, output, notes, citation` — insert `report` in a
  sensible spot, e.g. right after `output`).
- `docs/output.md`: add `report.html` to the consensus pipeline's ASCII tree diagram and its "Key
  files" table (both already exist near the top of that file for the consensus section).

## 14. Versioning

- `viralunity/__init__.py`: bump `__version__` from `"1.3.2"` to `"1.3.3"`.
- `Dockerfile`: bump the `LABEL version="1.3.2"` line to `"1.3.3"` (line ~2, per `RELEASING.md`'s
  instructions — check that file for the full release checklist, though **do not** actually tag or
  push a release; that's the maintainer's call).
- `CHANGELOG.md`: add a new `## [1.3.3] - <today's date>` entry (follow the existing entries'
  `### Added`/`### Changed`/`### Fixed` subsection style) describing the interactive HTML report
  feature and the `generate_html_report` toggle.
- Do **not** push a release tag — leave tagging/PyPI publishing to the maintainer, per explicit
  instruction.

## 15. Suggested implementation order (TDD)

1. Create branch `feature/interactive-visualizations-v1.3.4` off `feature/filters-v1.3.2` (check
   it still exists / hasn't merged first — if it has, branch off `main` instead and note the
   deviation in your final summary).
2. Write failing unit tests (§11.1) against not-yet-created `generate_consensus_report.py`.
3. Implement the pure functions until green, in dependency order: parsing → dedupe/sum →
   downsample → figure builders.
4. Add checked-in fixtures (§11.2) + the failing end-to-end test.
5. Implement `render_report`/`write_report` + the Jinja2 template + toggle JS + CSS until the e2e
   test passes.
6. Add the `run_cli()`/`run_snakemake()` dual-mode wrapper.
7. Add `viralunity_report.py` and register it in `viralunity_cli.py`; manually smoke-test against
   both example output dirs (§12, first two commands only at this point).
8. Wire the 6 config touch points (§9) + the 4 `.smk` rule/`rule all` edits (§7.1) +
   `envs/report.yaml` (§8.3) + the `pyproject.toml` dependency and package-data fixes (§10).
9. Run `make test-dryrun`; fix any DAG issues.
10. Write/link the docs page (§13), bump versions and CHANGELOG (§14).
11. Full loop: `make lint`, `make test`, `make test-dryrun` (§11.4) — all green.
12. Final full manual validation per §12, including the real end-to-end pipeline run.

## 16. Definition of done

- `viralunity report --input outputs/test_consensus_illumina --output /tmp/r1.html` and the
  influenza equivalent both produce single self-contained HTML files that open offline and match
  the layout in §5.
- Running a consensus workflow (any of the 4 entry points) auto-produces `<output>/report.html`
  when `generate_html_report` is left at its default (`True`), and does **not** produce it when
  `--no-generate-html-report` is passed.
- `make lint`, `make test`, and `make test-dryrun` all pass.
- The new docs page is linked from `docs/index.md` and builds without Sphinx errors (if you have a
  way to build the docs locally — check for a `docs/Makefile` or similar; if none exists, at least
  visually confirm the MyST syntax is well-formed).
- Version bumped to 1.3.3 in both `viralunity/__init__.py` and the `Dockerfile`, with a matching
  `CHANGELOG.md` entry, and **no** git tag pushed.
