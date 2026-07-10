# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

ViralUnity is a Python package whose only job at runtime is to validate inputs, write a YAML config, and launch one of six Snakemake workflows for viral HTS analysis. All of the actual bioinformatics — alignment, classification, assembly, consensus calling — lives in the Snakemake files under `viralunity/scripts/`. Treat the Python layer as a CLI + orchestration shim; treat the `.smk` rule files as the substantive code.

## Common commands

Setup is conda-based (the per-rule tool stacks live in `viralunity/scripts/envs/*.yaml` and are pulled in by Snakemake's `--use-conda`):

```bash
conda env create -n viralunity -f environment.yml
conda activate viralunity
pip install -e ".[dev]"
```

Development loop:

```bash
make test            # python -m unittest discover ./test -p *test.py
make test-dryrun     # pytest test/viralunity_dryrun_test.py -v   (snakemake -n on every workflow)
make lint            # black --check viralunity/ test/  +  ruff check viralunity/ test/
make format          # black + ruff --fix
```

Run a single test file or test:

```bash
python -m unittest test.viralunity_consensus_test -v
python -m unittest test.viralunity_consensus_test.Test_Validate.test_validate_ok
pytest test/viralunity_dryrun_test.py -v -k metagenomics_illumina
```

Pipeline invocations (the user runs these against real data; for editing, prefer `--create-config-only` plus the dryrun tests to avoid 30+ minute runs):

```bash
viralunity create-samplesheet --input <runs-dir> --output samples.csv
viralunity consensus  illumina --sample-sheet ... --reference ... --config-file ... --output ...
viralunity consensus  nanopore --sample-sheet ... --reference ... --config-file ... --output ...
viralunity meta       illumina --sample-sheet ... --kraken2-database ... --config-file ... --output ...
viralunity meta       nanopore --sample-sheet ... --kraken2-database ... --config-file ... --output ...
viralunity get-databases all --path databases/
viralunity build-deacon-index --input host.fasta
```

If `pytest` on the dryrun suite fails with `AttributeError: module 'pulp' has no attribute 'list_solvers'`, pin `pip install 'pulp<3'` — Snakemake 7.32 is incompatible with newer pulp.

## Architecture

### CLI → orchestrator → Snakemake

The single console entry point `viralunity` (declared in `pyproject.toml`'s `[project.scripts]`) lives in `viralunity/viralunity_cli.py` and is a Click group. Each subcommand has its own `*_cli.py` module that defines Click options, then hands a plain `args` dict to a `main()` in a non-CLI module (`viralunity_consensus.py`, `viralunity_meta.py`, etc.).

Both `consensus` and `meta` `main()` functions are thin wrappers around `_orchestrator.run_pipeline(...)` in `viralunity/_orchestrator.py`. The orchestrator runs four callbacks in order: `resolve_paths` → `validate` → `generate_config` → `run_workflow_fn`, with shared error handling and a `--create-config-only` short-circuit. The pipeline-specific modules still own `validate_args`, `generate_config_file`, and `run_snakemake_workflow` so existing test patches at those module-level names keep working — do not move those names back into `_orchestrator` without updating the tests.

`run_snakemake_workflow` picks the workflow file by formatting a path: `viralunity/scripts/{consensus,metagenomics}_{illumina,nanopore}[_segmented].smk`. Segmentation is selected by detecting that `args["reference"]` is a dict (built from repeatable `--segmented-reference SEGMENT=PATH` options) rather than a string. There is no segmented variant of the meta workflows.

### Config file is the contract

`viralunity/config_generator.py` (`ConfigGenerator`) writes the YAML that Snakemake reads. The keys it emits are the same strings hard-coded in the `.smk` files (e.g., `config["samples"]`, `config["output"]`, `config["run_denovo_assembly"]`). If you add a new pipeline option, you touch four places: the Click option in the relevant `*_cli.py`, the `validators.py` checks, a `ConfigGenerator.add_*` setter, and the rule(s) that read it. Constants for these key names live in `viralunity/constants.py` (`ConfigKeys`) but the `.smk` files reference the raw strings directly, so renaming a key means grepping `viralunity/scripts/` too.

`ConfigGenerator.add_resource_settings(args, rule_list)` emits one `{rule}_cpus` / `{rule}_ram` pair per rule. The rule lists are declared as class attributes on `ResourceDefaults` in `constants.py` (`CONSENSUS_ILLUMINA_RULES`, `META_SHARED_RULES`, etc.). When you add a computationally heavy rule, add it to the right list so its resources land in the generated config.

Three tool-level flags are config-only (no CLI), defaulted to the historical values for backwards compatibility:
- `minimap2_consensus_align_flags` (consensus, both data types)
- `diamond_max_target_seqs` (meta)
- `kraken2_extra_flags` (meta)

Users edit them by re-running with `--create-config-only`, editing the YAML, and rerunning Snakemake directly. Don't promote them to CLI flags without discussion.

### Snakemake workflows

Each top-level `.smk` in `viralunity/scripts/` is small — it sets up wildcard helpers, computes `_all_inputs()` based on which optional steps are enabled (`run_denovo`, `run_k2_reads`, `run_diamond_reads`, `has_negative_controls`, `run_reference_assembly`, …), and `include:`s rule modules from `viralunity/scripts/rules/`. The rule modules are split by phase (`qc_illumina.smk`, `alignment_*.smk`, `metagenomics_dehost_*.smk`, `metagenomics_kraken2_{reads,contigs}_*.smk`, `metagenomics_diamond_{reads,contigs}_*.smk`, `metagenomics_assembly_*.smk`, `metagenomics_reference_assembly.smk`, `stats.smk`, `consensus_*_common.smk`).

A few cross-cutting conventions to know before editing rules:

- **Per-rule conda envs.** Every rule has `conda: "envs/<name>.yaml"` (relative to the workflow file). Adding a new tool means either reusing an env or adding a YAML there. Snakemake's `--use-conda` is enabled in `_orchestrator.run_workflow`.
- **Outputs are computed conditionally.** The top-level `rule all` calls a `_all_inputs()` function that appends targets based on config booleans. Each optional branch (Kraken2 reads, Diamond contigs, reference assembly, negative-control filter, …) shows up both there and inside `organize_files`'s `expand(... if <flag> else [])` inputs. Whenever you add an optional step, edit both places.
- **`organize_files` is the symlink terminus.** It is the last rule before `benchmark.tsv` and creates the per-sample `samples/<sample>/...` symlinks that users actually browse. New per-sample outputs need a `ln -sf` block there to be discoverable.
- **Krona "raw vs filtered" pairs.** Every classifier × source (kraken2_reads, kraken2_contigs, diamond_reads, diamond_contigs) emits two Krona HTMLs: `*.krona.html` (raw) and `*.filtered.krona.html` (post bleed / negative-control filter). The filter is lineage-aware and lives in `viralunity/scripts/python/filter_krona_by_pass_taxids.py`.
- **Reference sanitization (nanopore).** The nanopore consensus pipeline sanitizes reference FASTA headers (replacing `/ \ | , ~` and spaces with `_`) before use because clair3 makes per-contig directories from the seq IDs. Don't bypass this.

### Reference assembly (meta)

When `--run-reference-assembly` is set, `metagenomics_reference_assembly.smk` adds a Snakemake checkpoint that reads the bleed-filtered taxa summary TSVs, picks reference accessions, and triggers per-`{family}_{accession}` (`ref_key`) consensus assembly using the same alignment/consensus rules as `viralunity consensus`. Two strategies live behind `--reference-selection-strategy`:

- `taxid` (default): exact taxid lookup in `--viral-taxids`, with a species-level fallback resolved through `--taxdump`. Fast; assumes the classifier DB and `--viral-genomes` come from the same RefSeq release.
- `similarity`: blastn from de novo contigs against `--viral-genomes`, requires `--run-denovo-assembly` and a prebuilt BLAST index alongside the FASTA.

Both validate the resulting taxids against `--taxdump` to confirm they trace to a target family.

### Sample sheets

CSV with no header. Illumina has 3 columns (`sample_id,R1,R2`), Nanopore has 2 (`sample_id,fastq`). `create-samplesheet` builds them by scanning a run directory; the parser (`viralunity/validators.py:validate_sample_sheet`) keys off column count, not the data-type flag of the file itself, so a malformed Nanopore CSV that happens to have 3 columns will be silently parsed as Illumina.

Sample names are prefixed with `sample-` inside the generated YAML by `ConfigGenerator.add_samples` — the `.smk` files refer to `sample-<id>` everywhere, but users only ever see `<id>` in their inputs/outputs. Tests and config inspection should expect the prefixed form.

## Tests

The unittest suite under `test/` covers Python-layer behavior (CLI parsing, validators, config generation, path resolution). The Snakemake dryrun suite (`test/viralunity_dryrun_test.py`) runs `snakemake -n` against every workflow + a YAML in `test/dryrun_configs/`; the placeholder fixture `create_dryrun_placeholders.sh` writes empty input files so paths resolve. CI runs lint, the unit suite, and the dryruns in three separate jobs.

When changing rule wiring, run `make test-dryrun` — it catches missing inputs, broken `expand` patterns, and circular dependencies that the Python suite cannot see.

## Editing notes

- Don't rename `validate_args`, `generate_config_file`, or `run_snakemake_workflow` in `viralunity_consensus.py`/`viralunity_meta.py` — multiple tests patch those exact module-level names.
- The `viralunity/scripts/python/*.py` scripts run inside Snakemake's `script:` directive, which injects a `snakemake` global. Ruff would otherwise flag F821; `pyproject.toml` already silences this for that path. Don't add `from snakemake import snakemake` to those files.
- `viralunity` is a published console script (`pyproject.toml: [project.scripts]`). After a fresh checkout, `pip install -e .` is required before `viralunity --help` works — `make install-dev` does this.
- Version is single-sourced from `viralunity/__init__.py:__version__`. Bumping it for release also requires editing the `Dockerfile` LABEL — see `RELEASING.md`.
- `my_results/` and `my_test_data/` are the user's local working directories, not test fixtures; they are gitignored. Use `test/dryrun_configs/` for additions to the regression suite.

## Versioning policy

Follow Semantic Versioning (`MAJOR.MINOR.PATCH`) with these conventions:

- **PATCH** (`x.y.Z`) — small code adjustments: bug fixes, minor tweaks, docs/CI/packaging changes that don't add user-facing features.
- **MINOR** (`x.Y.0`) — larger blocks of changes: a group of significant features or a substantial enhancement. Reset PATCH to 0.
- **MAJOR** (`X.0.0`) — reserved for exceptional circumstances and bumped **manually by the maintainer only**. Do not bump MAJOR on your own; if a change seems to warrant it, propose it and let the user decide.

When cutting any release, follow the full `RELEASING.md` flow (bump `__version__` + `Dockerfile` LABEL, update `CHANGELOG.md`, tag `vX.Y.Z`) — pushing the tag triggers the PyPI publish workflow.
