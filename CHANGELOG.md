# Changelog

All notable changes to this project will be documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aspires to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The release process is documented in [RELEASING.md](RELEASING.md).

## [1.4.0] - 2026-07-16

Adds an interactive, self-contained HTML report for consensus runs, designed to
stay legible from 1 to ~96 samples (≈96×8 for segmented viruses).

### Added

- Interactive consensus HTML report: a single offline-openable `report.html`
  visualizing a consensus run —
  - **Summary KPI tiles**: samples analyzed (with an "N segments each"/amplicon
    subtitle), samples ≥ the pass threshold (with "% of run"), samples below the
    warn threshold (tile turns red when non-zero), median horizontal coverage, and
    mean depth. Thresholds are configurable and every label is derived from them.
  - **Assembly-statistics data table**: a sample/segment search box, a "low
    coverage only" filter, a live row count, worst-coverage-first default sort, and
    visual pass/fail encoding (status dot + tinted coverage cell + inline bar).
    Segmented runs roll up to one summary row per sample (mean coverage, median
    depth, whole-sample mapped %) that expands to per-segment rows, plus a
    segment-focus chip row.
  - **Sequencing throughput**: horizontal stacked bars (mapped / QC-passed-unmapped
    / removed-by-QC) summing to total reads in reconciled units, with an
    Absolute/Percent toggle; height grows with sample count.
  - **Coverage heatmap** (replacing per-sample line overlays): samples worst-first
    on the y-axis; a position-depth mode (binned genome position × depth, with a
    Natural/Log10 toggle and the annotation band) and, for segmented runs, an "All
    segments" grid of per-segment horizontal coverage %.
  - **By-sample panel**: a searchable, worst-first sample list (status dot +
    coverage badge) beside the selected sample's depth plot with 20x/100x guides;
    segmented samples default to an "All (concatenated)" view plus a per-segment
    selector.
  Depth charts default to honest zeros (no coverage reads as a gap, never a false
  depth of 1). Charts use a colourblind-safe palette and support light/dark mode
  and print/PDF. Built with Plotly, inlined once so the file has no external
  dependencies.
- Report-generation parameters, as `viralunity report` flags (or config keys):
  `--pass-threshold` (0.90), `--warn-threshold` (0.70), `--chart-color`, and
  `--colorbar-thickness`.
- Run-parameters slide-over drawer showing the full run config, HTML-escaped and
  behind a focus-trapped modal.
- Optional annotation tracks beneath the by-sample and position-mode coverage
  views: a **gene track** from a new `--gene-annotation FILE` /
  `--segmented-gene-annotation SEGMENT=PATH` GFF3 input, and a **primer-scheme
  track** (amplicons paired into pool rows) from the existing `--primer-scheme`
  BED. Both are drawn as rectangles aligned to the genome x-axis, with Genes /
  Primers on-off toggles. The files are staged into `<output>/annotation/` so the
  report and `viralunity report` keep working after the run. Absent when no
  annotation is provided.
- `viralunity report --input <consensus-output-dir> [--output report.html]`: a
  standalone subcommand to (re)generate the report from an existing output
  directory without rerunning the pipeline.
- `--generate-html-report / --no-generate-html-report` (default on) on both
  `viralunity consensus illumina` and `viralunity consensus nanopore`: the
  report is produced automatically at the end of every consensus run and can be
  opted out of on very large runs.

### Changed

- New runtime dependencies `plotly` and `jinja2`; the packaged `.j2` template is
  now shipped as package data.

### Fixed

- A segment literally named `NA` (influenza's neuraminidase segment) is no longer
  parsed as `NaN` when reading the stats summary, which previously crashed
  coverage-path resolution on real NA-segment runs.

## [1.3.3] - 2026-07-14

Documentation-only release: the user-facing docs are brought in line with the
1.3.2 CLI and behaviour, and the changelog is tidied. No code or workflow changes.

### Documentation

- Removed the stale Portuguese-docs language switcher from `docs/index.md` (the
  `docs-pt/` tree was retired in a prior PR).
- Corrected doc/code drift in `docs/`: the `--evalue` default (`1e-10`), the
  `--method`/`--source` reference-assembly flags (no default; required with
  `--run-reference-assembly`), the subcommand count (six — adds `setup`), and
  replaced references to a non-existent `get-databases ictv-vertebrate-taxids`
  command with the `build_ictv_vertebrate_taxids.py` script.
- Documented 1.3.2 behaviour that had not reached the docs: the
  `absent_from_controls` negative-control decision, the per-rank
  `summaries/{family,genus,species}/` + `summaries/full/` output layout, and the
  `.ctgstats` step in the taxa-summary filter chain.
- Added the PyPI install path (`pip install viralunity`) to the installation docs
  and documented the auto-generated per-rule `--<rule>-cpus` / `--<rule>-ram`
  options.

### Changed

- `CHANGELOG.md` housekeeping: folded the stale `[Unreleased]` section into
  `[1.2.0]` (where the `viralunity setup` subcommand, `--conda-prefix`, and the
  conda-26.x fix actually shipped), dropped entries already recorded under
  `[1.3.0]`, and corrected an erroneous `get-databases ictv-vertebrate-taxids`
  command reference — restoring reverse-chronological order.

## [1.3.2] - 2026-07-11

Empirically-motivated refinements to the metagenomics contamination filters,
grounded in the REVISA lab-contamination episode (mayaro + HIV libraries seen
across negative controls).

### Changed (breaking)

- **Per-rank output layout.** Each track's taxa summary now lives under a
  `taxonomic_assignments/<track>/summaries/` tree: one table per taxonomic rank
  under `summaries/{family,genus,species}/` (the user-facing deliverable, most
  users want `species/`), with higher-rank names propagated down (species tables
  gain `family`+`genus` columns; genus tables gain `family`). The combined
  cumulative filter chain lives internally under `summaries/full/`, and the
  per-sample intermediate tables under `summaries/per_sample_summaries/`.
  Downstream scripts that read the old flat `<track>/<track>_taxa_summary_*.tsv`
  path must switch to `<track>/summaries/species/` (or `<track>/summaries/full/`
  for the combined table).
- **Negative-control log-ratio is now log10, not log2.** The `log2_ratio` column,
  `neg_decision` values, `--log2-ratio-threshold` flag and config key are renamed
  to their `log10` equivalents. The default threshold stays `1.0` but now marks a
  **10-fold** enrichment over controls (log10 = 1) rather than the old 2-fold
  (log2 = 1). This gate drives `neg_pass` when exactly one control is present or
  the control variance is zero; with ≥2 controls the primary z-score gate
  (unchanged) decides, so the tightening applies to that fallback path.
- **Removed the `source` column** from taxa summaries (it held an internal
  krona-input path of no user interest).

### Added

- **`final_species` column** — the confirmed species call, coalescing
  `nr_correct_species` (NR's correction where it disagreed) with the original
  `name`; present on every track.
- **Aggregate (pooled) negative-control filter** — treats all controls as one
  pooled library (raw reads pooled, weighting each control by its library size)
  and reports `pooled_control_metric`, `agg_fold_enrichment`, `agg_log10_ratio`,
  and `agg_fold_enrichment_10x/100x_pass`. Complementary to the z-score; catches
  widespread, high-variance contamination that inflates per-control SD.
- **Convenience pass-flag columns** — `fold_enrichment_10x_pass`,
  `fold_enrichment_100x_pass`, `neg_pass_5`, `neg_pass_10` (`>=` inclusive; NA
  where the underlying statistic is NA).
- **Largest-contig statistics** (both contig tracks, with `--viral-genomes`):
  `largest_contig_bp`, `largest_contig_ref_coverage_pct`, and
  `largest_contig_median_depth`, a cheap genome-fraction/coverage proxy. Each
  contig track remaps the host-filtered reads to its own viral contigs and runs
  `samtools depth -a` for per-contig depth (kraken2_contigs extracts its viral
  contigs by lineage; diamond_contigs reuses its existing viral remap).
  `largest_contig_ref_coverage_pct` = `largest_contig_bp / genome_length_bp * 100`
  is a preliminary genome-completeness estimate (raw/uncapped, so >100% flags a
  contig longer than the median reference; NA without a reference length).

### Changed

- **Bleed filter is metric-aware.** It now uses RPKM when `--viral-genomes` is
  supplied (per taxon), else RPM, recorded in a new `bleed_metric` column; the
  group-max column `max_rpm` is renamed `bleed_max`. Because bleed is a
  within-taxon ratio the pass/fail is unchanged by the metric; the application
  floor is now metric-specific (`bleed_rpm_floor` 1.0, new `bleed_rpkm_floor`
  0.1 — both YAML-only).
- **Taxa absent from all negative controls now auto-pass the enrichment gates**
  (science-behaviour change — flagged for PI review). Previously a taxon never
  seen in any control still got a finite `fold_enrichment = (sample + pc)/(0 + pc)`
  that tracked only the sample's own abundance, so low/moderate-abundance
  control-absent taxa could *fail* `fold_enrichment_100x_pass` /
  `agg_fold_enrichment_100x_pass` (pervasive on contig tracks, where RPKM is
  tiny). Now, when a taxon's `control_max == 0`, the meaningless ratios
  (`fold_enrichment`, `log10_ratio`, `agg_fold_enrichment`, `agg_log10_ratio`)
  are blanked to NA and every enrichment gate — the four fold-enrichment
  10x/100x pass flags **and** the primary `neg_pass` — is set True, with
  `neg_decision = "absent_from_controls"`. Absence from controls is itself the
  evidence the taxon is not control-borne contamination. The z-score gates
  (`neg_pass_5`/`neg_pass_10`) stay NA (a z-score is undefined without control
  signal). This changes which taxa survive the lineage-aware Krona "filtered"
  output and reference-assembly selection.
- **ICTV host-filter docstrings corrected.** The filter runs *last* in the chain
  (after bleed/negative-control), not before as three docstrings claimed. The
  behaviour is unchanged and numerically equivalent to running it earlier, since
  bleed/neg are per-taxon ratios; only the comments were wrong.

## [1.3.1] - 2026-07-10

### Added

- **`viralunity get-databases nr`** — download and configure the NCBI **nr** protein
  database used by `meta --run-nr-validation`, closing the one remaining bring-your-own
  gap in `get-databases`. By default it fetches NCBI's preformatted, md5-verified BLAST+
  nr volumes with the official `update_blastdb.pl` tool and runs `diamond prepdb` so
  DIAMOND (and its `staxids`) are ready. Two bring-your-own modes cover HPC/air-gapped
  setups: `--from-blastdb <prefix>` registers an existing nr, and `--from-fasta <faa>`
  builds a native `nr.dmnd` via `diamond makedb` (with optional
  `--taxonmap`/`--taxonnodes`/`--taxonnames` for taxonomy). Also supports `--source`
  (ncbi/aws/gcp), `--threads`, and `--skip-prepdb`. nr is intentionally left out of
  `get-databases all` because of its size (100+ GB).

## [1.3.0] - 2026-07-09

### Added

- **ICTV vertebrate-virus taxonomic filter** (`viralunity meta --run-ictv-host-filter`)
  — optional false-positive removal that drops hits outside a vertebrate-infecting
  allowlist built from the ICTV VMR sheet (`build_ictv_vertebrate_taxids.py`), matched
  lineage-aware. Requires `--ictv-vertebrate-taxids-file`.
- **NR validation of viral contigs** (`--run-nr-validation`) — combines the viral
  contigs, searches them against a DIAMOND NR database, assigns an LCA verdict, and
  harmonizes the result back into the contig tracks as an `nr_pass` column (dropping
  contigs that fail). New options: `--nr-diamond-database`, `--nr-evalue`,
  `--nr-max-target-seqs`, `--nr-sensitivity`, `--nr-consensus-threshold`.
- **Aggregated contig DIAMOND search** — optional combine → search → split path that
  reduces redundant DIAMOND invocations across samples.
- A dedicated taxonomic-filter stage that runs **before** the bleed / negative-control
  filters, with each track's taxa summary flowing through one cumulative, consistently
  named chain: `…_taxa_summary[_RPM|_RPKM][.nr].bleed[.neg][.ictv].tsv`.

### Changed

- **Negative-control enrichment now zero-fills absent controls.** Control mean/SD are
  computed over *all* negative controls (a taxon undetected in a control counts as 0),
  and the z-score gate uses the same denominator as the pass/fail decision. This changes
  which taxa pass the contamination filter relative to 1.2.0 (which used only the controls
  where a taxon was detected). No division by zero — the z-score is `NA` when the control
  SD is 0.
- `make test-empirical` runs inside the `viralunity` conda env, so it no longer fails
  under a Python 3.12 base environment.
- Pinned the `genome_selection` conda env (`python=3.11`, `pandas=2.0`, `blast=2.16`).
- Assembly-stats read counts are computed as `lines // 4` instead of counting `+`
  separator lines (which miscounted `+<id>` separators and Phred-10 quality lines).

### Fixed

- `convert_diamond_output_to_krona_input.py` crashed on CLI use (`argparse` was never
  imported).
- Reference selection silently produced an empty `reference_targets.tsv` when a summary
  contained a blank taxid (float coercion); taxids are now read as strings.
- Hardened taxonomy post-processing (`annotate_diamond_taxonomy`, `summarize_krona_taxa`)
  and removed error-swallowing patterns in the dehosting/alignment/consensus shell rules.
- `create-samplesheet` writes CSV via `csv.writer` and sanitizes sample names.

### Engineering

- Input hardening for service use: `run_name` sanitization and numeric-parameter range
  validation (threads, thresholds, e-value); `wildcard_constraints` on every workflow.
- Observability: `onerror`/`onsuccess` handlers and per-rule logs on the Clair3 and
  MultiQC rules.
- Run-provenance manifest now records the config SHA-256 and a post-run completion status.
- CI: `mypy` cleared to zero errors and promoted to a **gating** check; new unit tests
  take every metagenomics metric-path script off 0% coverage.

### Removed

- Nothing user-facing. (A KrakenUniq-style minimizer filter was prototyped during
  development and removed before release.)

## [1.2.0] - 2026-07-07

### Added

- **Published to PyPI** — ViralUnity is now installable with `pip install viralunity`.
  Releases are built and uploaded automatically via GitHub Actions trusted publishing
  (OIDC) when a `vX.Y.Z` tag is pushed; see [RELEASING.md](RELEASING.md). Note that
  conda/mamba is still required at runtime, since Snakemake builds the per-rule tool
  environments via `--use-conda`.
- **RPKM normalisation** — when `--viral-genomes` (RefSeq FASTA) and `--viral-taxids`
  (genome2taxid TSV) are provided, a per-taxon genome-length table is computed (median
  genome length at family, genus, and species level across all accessions under each node).
  RPKM is then emitted alongside RPM in a new `*_taxa_summary_RPKM.tsv` intermediate file
  for every classifier × mode track (kraken2/diamond × reads/contigs). RPKM is `NA` for
  taxa with no matching genome length.
- **Negative-control enrichment filter** — replaces the previous Poisson-based filter
  (`apply_negative_background_filter.py`) with interpretable statistics:
  - `fold_enrichment` and `log2_ratio` (always computed, configurable `--enrichment-pseudocount`).
  - `z_score` (when n_controls ≥ 2 and control SD > 0).
  - `neg_metric`: `rpkm` when available for that taxon, else `rpm`.
  - Tiered `neg_pass` gate: z-score (≥2 controls), log2-ratio (1 control or SD = 0 fallback),
    or NA (0 controls; treated as *keep* by the lineage-aware Krona filter).
  - All metrics and thresholds recorded in `*_RPM.bleed.neg.tsv` for full reproducibility.
- Three new `viralunity meta` CLI options: `--enrichment-pseudocount` (default `1.0`),
  `--z-score-threshold` (default `3.0`), `--log2-ratio-threshold` (default `1.0`).
- `viralunity/scripts/python/taxonomy.py` — shared NCBI taxonomy utilities
  (`load_taxdump`, `get_lineage`, `RANKS_OF_INTEREST`) extracted from
  `summarize_krona_taxa.py` and `filter_krona_by_pass_taxids.py` to eliminate code
  duplication; 16 unit tests added.
- `viralunity/scripts/python/build_genome_length_table.py` — builds the per-taxon median
  genome-length table from `.fai` + genome2taxid + taxdump; 18 unit tests.
- `viralunity/scripts/python/add_rpkm_to_summary.py` — merges genome lengths and computes
  RPKM; 17 unit tests.
- `viralunity/scripts/python/add_negative_control_enrichment.py` — new enrichment filter;
  49 unit tests covering all decision tiers, RPKM/RPM metric selection, SD = 0 fallback,
  absent-from-controls zero-background assumption, and the NA-as-keep contract.
- **`viralunity setup` subcommand** — pre-builds every per-rule conda environment into a
  shared cache at install time, so `viralunity consensus` / `viralunity meta` runs never
  have to materialize envs on the hot path. Options: `--conda-prefix PATH`, `--pipelines
  [consensus-illumina|consensus-nanopore|meta-illumina|meta-nanopore|all]` (repeatable,
  default `all`), `--threads INT`, `--dry-run`. Segmented workflow variants intentionally
  share envs with their non-segmented counterparts and are not separate selections.
- New `--conda-prefix` option on both `viralunity consensus` (illumina + nanopore) and
  `viralunity meta` (illumina + nanopore) so pipeline runs reuse the cache produced by
  `viralunity setup`. Picks up `$VIRALUNITY_CONDA_PREFIX` if set; otherwise defaults to
  `~/.cache/viralunity/conda-envs/`.
- `ConfigGenerator.write_skeleton(pipeline, data_type, config_path, placeholder_dir)`
  classmethod + `SKELETON_PLACEHOLDERS` map: emits a minimal YAML config sufficient for
  Snakemake to parse a workflow with `conda_create_envs_only=True`, paired with the set of
  empty input files the DAG build needs to touch.
- Docs: `docs/installation.md` "Troubleshooting" subsection (conda 26.x / bioconda shards
  404 symptom + escape hatch) and "First-time environment setup" subsection;
  `docs/tutorial/setup.md` step 1b "Build per-rule environments"; `docs/commands.md`
  `viralunity setup` section plus the `--conda-prefix` flag on the consensus and meta
  shared-options tables.
- Tests: `--conda-prefix` wiring on both CLIs, orchestrator forwarding, and skeleton +
  setup CLI coverage. +13 tests total.

### Changed

- `--negative-p-threshold` CLI option **removed** (previously set the Poisson p-value
  threshold for negative-control filtering). Replace workflow calls with the three new
  enrichment options above. Existing YAML configs that contain `negative_p_threshold` will
  have the key silently ignored by Snakemake (it is no longer read by any rule).
- The RPM denominator (and therefore RPKM) always reflects the read count *after* the
  dehosting step: dehosted reads (when dehosting is on), post-QC reads (Illumina, dehosting
  off), or raw reads (Nanopore, dehosting off). This is the correct normalisation base and
  was already the implicit behaviour — it is now documented and tested explicitly.
- `environment.yml`: pinned `conda` and `conda-libmamba-solver` to `>=24,<26` with an
  inline comment naming the lift condition (bioconda publishes shards OR conda's shards
  path tolerates 404s gracefully) — the immediate workaround for the
  `CreateCondaEnvironmentException` (see "Fixed" below).
- Default conda env cache path for `viralunity consensus` / `viralunity meta` is now
  `~/.cache/viralunity/conda-envs/` (override with `--conda-prefix` or
  `$VIRALUNITY_CONDA_PREFIX`). Snakemake previously materialized envs into
  `<workdir>/.snakemake/conda/` per working directory, so existing users' per-workdir
  caches become orphaned on first post-upgrade run and rebuild once into the shared cache.
  Subsequent runs are faster and survive across workdirs.
- `viralunity/_orchestrator.py:run_workflow` now forwards `args.get("conda_prefix")` to the
  `snakemake(...)` call. Absent key resolves to `None`, preserving the previous per-workdir
  behaviour for direct callers and tests that do not set the key.

### Removed

- `viralunity/scripts/python/apply_negative_background_filter.py` — the Poisson-based
  negative-control filter script. Replaced by `add_negative_control_enrichment.py`.

### Fixed

- **Reference assembly now selects from the post-filter taxa tables.** The
  `select_references_meta` checkpoint previously read the *raw* counts table
  (`*_taxa_summary.tsv`), so contaminants suppressed by the cross-sample filters could
  still trigger reference-guided consensus runs. It now reads the
  negative-control-filtered table (`*_RPM.bleed.neg.tsv`) when `--negative-controls` is
  set, else the bleed-filtered table (`*_RPM.bleed.tsv`), dropping taxa that explicitly
  fail `bleed_pass`/`neg_pass` (NA is kept; older outputs fall back to raw counts).
  New `--summary-suffix` arg on `select_reference_genomes.py`; 13 unit tests.
- **Shared-module imports work inside per-rule conda envs.** The taxonomy refactor made
  four `script:` scripts import `from viralunity.scripts.python...`, which is not importable
  in the per-rule conda envs and crashed the metagenomics rules at runtime. Added a
  sibling-module import fallback (`summarize_krona_taxa`, `filter_krona_by_pass_taxids`,
  `add_negative_control_enrichment`, `build_genome_length_table`).
- **Negative-control IDs are matched after `sample-` prefixing.** `--negative-controls`
  takes raw sample IDs, but the summary `sample` column is prefixed; the enrichment step
  matched the raw IDs and aborted with "None of the provided negative controls appear...".
  IDs are now prefixed to match, with fast-fail validation that each names a real sample.
- `viralunity setup` pre-builds the conda envs gated by every optional pipeline flag —
  the skeleton config previously left `run_isnv`, `run_denovo_assembly`, `run_diamond_*`,
  `run_reference_assembly`, and the meta-nanopore polish flags off, so Snakemake's DAG walk
  pruned the matching rules and their envs were never materialized. `setup --pipelines
  consensus-illumina` now builds 4 envs instead of 3 (adds `consensus.yaml`); the meta
  variants pick up `assembly.yaml`, `genome_selection.yaml`, and (nanopore) `medaka.yaml`.
  `SKELETON_PLACEHOLDERS` grew the reference-assembly placeholders needed for the expanded
  DAG.
- **Critical:** first pipeline run on a fresh install with conda 26.x +
  conda-libmamba-solver 26.x aborted on `Creating conda environment .../qc.yaml` with a
  `CreateCondaEnvironmentException`. Root cause is upstream: conda-libmamba-solver's
  "repodata shards" optimization queries `repodata_shards.msgpack.zst` from bioconda, which
  does not publish one; the resulting HTTP 404 trips conda's `RepodataIsEmpty` constructor,
  which calls `response.json()` on the non-JSON 404 body and crashes inside conda's own
  error path. Fixed by the `environment.yml` toolchain pin (see "Changed") and structurally
  neutralized by `viralunity setup` (see "Added"): containerless pipeline runs no longer
  touch the host solver after the cache is pre-warmed.

## [1.1.0] - 2026-05-21

### Added

- New `viralunity get-databases` subcommand group with downloaders for
  kraken2, krona, taxdump, virus-genome, host-genome, deacon-index, and
  diamond databases, plus an `all` aggregator and a `clean-protein-fasta`
  utility.
- New `viralunity build-deacon-index` subcommand to build a Deacon
  minimizer index from a host FASTA.
- DIAMOND-based protein-level classification (reads + contigs) for both
  Illumina and Nanopore metagenomics workflows.
- Reference-assembly checkpoint workflow that builds per-family
  consensuses from filtered Kraken2/DIAMOND hits.
- Filtered Krona reports aligned with bleed and negative-control taxa
  summaries.
- Per-rule conda environments (`viralunity/scripts/envs/*.yaml`) and
  `resources:` declarations on memory-heavy rules; cluster execution is
  now feasible.
- Real Snakemake dry-run integration tests
  (`test/viralunity_dryrun_test.py` + `test/dryrun_configs/`).
- Path-resolution unit tests (`test/viralunity_path_resolution_test.py`).
- Sphinx documentation under `docs/` (EN) and `docs-pt/` (PT) with full
  CLI flag tables; ReadTheDocs configuration.
- `pyproject.toml` migration (PEP 621): `[project]` metadata, dev extras
  (`black`, `pytest`, `ruff`), `[tool.black]`, `[tool.ruff]`,
  `[tool.mypy]`, `[tool.pytest.ini_options]` configurations.
- CI hardening: `push: main` trigger, dedicated `lint` job (black --check,
  ruff check), dedicated `dryrun` job; CI now runs on PRs *and* on
  pushes to `main`.
- Three previously-hardcoded tool flags are now config-driven keys
  (`minimap2_consensus_align_flags`, `diamond_max_target_seqs`,
  `kraken2_extra_flags`) with the historical defaults preserved.
- Community files: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` (Contributor
  Covenant 2.1), `.github/ISSUE_TEMPLATE/`, `.github/PULL_REQUEST_TEMPLATE.md`.
- `RELEASING.md` documenting the version-bump and tag workflow.
- `--clair3-model` flag on `viralunity meta nanopore`, mirroring the
  consensus CLI option. Only consulted by the reference-assembly
  path's `infer_consensus_sequence` rule; falls back to
  `r1041_e82_400bps_sup_v500` when omitted.
- English `docs/output.md` documenting consensus + meta output
  directory layouts (parity with PT `docs-pt/saida.md`); fixes the
  previously-dead `output` entry in the EN toctree.

### Changed

- Migrated from `argparse` to Click for all CLI surfaces; shared options
  are stacked via `_add_common_options` and `_add_resource_options`
  decorators.
- Granular Snakemake rule organisation: 21 files under
  `viralunity/scripts/rules/` (vs 5 on `main`).
- Aligned dependency pinning between `environment.yml` and
  `pyproject.toml` on `>=X.Y` semantics, with `<8` upper bound on
  Snakemake (its 8.x release dropped the Python API this package calls).
- LICENSE copyright extended to `2021-2026`.
- README restored to a usable shape: install + quick-start blocks for
  each subcommand, with the ReadTheDocs link kept as the canonical
  source.
- Dockerfile `LABEL description` now matches `viralunity._description`
  ("A pipeline for viral metagenomics analysis.") instead of the
  generic "Docker image for viralunity" stub.
- `viralunity meta` no longer passes `forceall=True, lock=False,
  workdir=os.getcwd()` to `snakemake(...)`. The pipeline now respects
  Snakemake's caching and matches the consensus pipeline's invocation;
  the workdir defaults to cwd, which is what `resolve_path_args`
  already prepares paths against.
- `viralunity meta nanopore`'s `--run-polish-racon`,
  `--run-polish-medaka`, and `--medaka-model` flags now actually
  populate the generated YAML config. Previously they were accepted
  by Click but `generate_config_file` never invoked
  `add_nanopore_settings`, so the keys silently vanished before
  reaching the workflow.
- Docs review for v1.1.0:
  - `docs/usage.md` `get-databases` tree expanded with
    `clean-protein-fasta`, `virus-genome`, and `all`.
  - `docs/installation.md` gained a "Development install" section
    with `pip install -e ".[dev]"` and a `CONTRIBUTING.md` link.
  - `docs/commands.md` and `docs-pt/comandos.md` document the three
    YAML-only configuration keys (`minimap2_consensus_align_flags`,
    `diamond_max_target_seqs`, `kraken2_extra_flags`).
  - `docs-pt/comandos.md` adds the missing Illumina level-0
    `create-samplesheet` example (parity with EN).
  - `docs/notes.md` + `docs-pt/notas.md`: the dynamic-reference-
    assembly strategy comparison and required-databases tables now
    correctly state that both `taxid` and `similarity` strategies
    use `--viral-taxids` and `--taxdump`.
  - Sphinx copyright in `docs/conf.py` and `docs-pt/conf.py` bumped
    from "2025" to "2021-2026" to match `LICENSE`.

### Fixed

- **Critical**: `viralunity/validators.py` raised `AdaptersNotFoundError`
  without importing it; any Illumina run providing a non-existent
  `--adapters` path triggered a `NameError` instead of the intended
  `AdaptersNotFoundError`. Existing tests mocked the validator so this
  did not surface in CI.
- Dockerfile `LABEL version` bumped from `1.0.3` to `1.1.0` to match
  `viralunity/__init__.py:__version__`.
- Replaced `print("No samples were provided.")` in
  `viralunity/viralunity_meta.py` with `logger.warning(...)`.
- Added `from e` exception chaining to all custom-exception re-raises
  inside `except ... as e` blocks in `viralunity/validators.py` and
  `viralunity/config_generator.py`.
- Replaced deprecated `pd.read_csv(..., delim_whitespace=True)` with
  `sep=r"\s+"` in
  `viralunity/scripts/python/calculate_assembly_stats.py`.
- Removed `subprocess.Popen(..., shell=True)` from
  `calculate_assembly_stats.py`; now uses `subprocess.run([...])` with
  list args so the path is not interpreted by a shell.
- Added `set -euo pipefail` to every multi-line shell block across the
  rule files so piped commands fail fast instead of letting an upstream
  error get swallowed by a downstream `samtools sort`.
- `make test` now depends on `make install` so a fresh clone runs
  cleanly.
- Removed leftover `print(os.environ.get("PATH"))` debug calls from
  `test/viralunity_consensus_test.py` and `test/viralunity_meta_test.py`.
- `metagenomics_nanopore.smk`'s `organize_files` benchmark-aggregation
  rule wrote literal `\t` strings into the first column of
  `benchmark.tsv`. The shell block had been declared as a Python
  *raw* triple-quoted string (`r"""..."""`), so the `\\t` escapes in
  the embedded awk no longer collapsed to tabs. Reverted to a
  regular triple-quoted string, matching the Illumina counterpart.
- `select_reference_genomes.py` (taxid strategy): a missing or empty
  `--viral-taxids` previously produced a silent header-only
  `reference_targets.tsv`. Now hard-fails early with an actionable
  message that points users at
  `viralunity get-databases virus-genome`.
- `select_reference_genomes.py` (taxid strategy): also index each
  genome at its species-rank ancestor when building the
  taxid → accessions map. Handles cross-database taxid drift such
  as the NCBI promotion of "Betacoronavirus pandemicum" to a species
  taxid (3418604) above the historical SARS-CoV-2 strain taxid
  (2697049) that NCBI Datasets / RefSeq still uses for
  `NC_045512.2`.
- `extract_reference_fasta` rule: now produces the `.fai` sibling
  alongside the per-sample reference FASTA via `samtools faidx`, so
  the downstream `clair3` and `bcftools consensus` calls can find
  the index they require.
- `infer_consensus_sequence` rule: `clair3_model` fallback bumped
  from `r1041_e82_400bps_sup_v420` (not shipped with the clair3
  conda env we install) to `r1041_e82_400bps_sup_v500` (matches the
  consensus CLI default and the env's installed model).

### Refactored

- Renamed the custom `FileNotFoundError` exception to
  `ViralUnityFileNotFoundError` so it no longer shadows Python's
  builtin.
- Extracted a shared `viralunity._subprocess.run_command` helper out of
  the duplicated `_run` helpers in
  `viralunity_get_databases_cli.py` and
  `viralunity_build_deacon_index_cli.py`.
- Consolidated the two near-identical `data_report.jsonl` parsers in
  `viralunity_get_databases_cli.py` into a single parametric function
  driven by a `key_extractor` callable.
- Extracted shared consensus rules (`generate_multiqc_report`,
  `calculate_assembly_statistics`, `align_consensus_to_reference_genome`)
  into `rules/consensus_illumina_common.smk` and
  `rules/consensus_nanopore_common.smk`; the four top-level snakefiles
  now keep only the segment-specific rules.
- Lifted shared pipeline orchestration (`ConfigGenerator` preamble, the
  `snakemake(...)` invocation, the `main` try/except skeleton) into
  `viralunity/_orchestrator.py`.
- Expanded one-line docstrings on all consensus/meta Click handlers and
  added type hints to all 18 Click handler signatures.
- Polished helper scripts (`calculate_assembly_stats.py`,
  `rename_sequences.py`): module + function docstrings, type hints,
  renamed the `input` parameter that shadowed the Python builtin.
- Removed two no-op `args.update({k: v for k, v in kwargs.items() if k
  not in args})` lines in `viralunity_meta_cli.py` that obscured the
  intent of `_build_meta_args`.

## [1.1.0] — branch-vs-main overview

See git history for the substantial changes on this branch versus
`main`: per-rule conda environments, resources declarations, set -euo
pipefail sweep, dryrun integration tests, sphinx docs, argparse→Click
CLI migration, new `get-databases` / `build-deacon-index` subcommands,
metagenomics expansion (DIAMOND, reference assembly, dehosting,
kraken2 reads/contigs).
