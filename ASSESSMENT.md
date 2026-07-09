# ViralUnity — Independent Engineering Assessment

*Senior/staff-level code review of the ViralUnity Snakemake + Python viral-metagenomics
pipeline. This is an assessment, not a refactor: no pipeline code was modified. Every claim is
backed by an executed command or a `file:line` citation, with fact separated from judgment.*

## Evidence base (what was actually run)

- **Tooling present:** snakemake 7.32.0, ruff 0.15.20, mypy 2.1.0, pytest 9.0.3, pytest-cov 7.1.0,
  conda 25.3.0, mamba 2.1.1, black 26.5.1.
  **Absent:** `snakefmt`, `singularity`/`apptainer`, and the bio binaries `diamond`/`blastn`/`samtools`
  — so no real rule execution was possible; Snakemake was exercised via `--lint` and dry-run only.
- `ruff check viralunity/ test/` → **All checks passed**. `black --check viralunity/ test/` → **pass, 78 files**.
  (`ruff format --check` reports 10 diffs, but the project's formatter of record is **black**, per
  `[tool.black]` in `pyproject.toml` and CI — so that is not a finding.)
- `mypy viralunity/` → **12 errors in 8 files** (CI runs this `continue-on-error`; documented backlog).
- `pytest -m "not empirical" --cov=viralunity --cov-branch` → **373 passed, 4 deselected; 55% total coverage**.
- `snakemake --lint` + `-n` dry-runs on consensus_illumina/nanopore + metagenomics illumina/nanopore +
  variants (`__nr`, `__ictv`, `__fullchain`) → **all exit 0; DAGs resolve; no ambiguous-rule warnings**.
- Direct source reads verifying every load-bearing claim (run_name handling, dead code, per-rule log
  coverage, wildcard ambiguity, taxid dtype coercion, the argparse CLI crash, the assembly-stats read
  counter, and the negative-control enrichment math).

> **Note for the PI:** finding **H1** (negative-control enrichment statistics) is a *science-behavior*
> question, not an infrastructure bug. It is documented here for review and should not be changed
> without sign-off.

---

## Executive summary

**Maturity: Functional (advanced beta) — not yet production-ready for an unattended public-health
webservice.**

ViralUnity is a genuinely well-engineered research pipeline that is far past prototype: a
pip-installable package, CHANGELOG/SemVer, a ReadTheDocs Sphinx site, a multi-job CI matrix
(ruff + black + mypy-advisory + pytest-coverage + dry-run + pip-audit/bandit + docker + conda-smoke),
per-rule pinned conda environments, a run-provenance manifest with per-input SHA256, a typed/coded
error hierarchy, safe subprocess handling, and 373 passing unit tests. Its lint and format gates are
green and its core scientific metrics are, on inspection, correct.

It falls short of *production-ready* on four axes:

1. **Input trust** — the generated config is validated only by hand-written presence/consistency
   checks (no JSON schema, no numeric range checks), and `run_name` reaches every output path
   unsanitized.
2. **Observability** — 60% of executing rules (112 of 186) declare no per-rule `log:`.
3. **Scientific-code assurance** — several analysis scripts on the metagenomics metric path have zero
   unit tests and contain real defects (a CLI crash, a fragile read counter, a silent taxid-coercion
   path), plus one methodology inconsistency in the negative-control filter.
4. **Reproducibility** — conda envs pin only minor versions, one env floats `pandas` entirely, and the
   `conda-lock.yml` the Makefile references is not committed.

None of these is an *active* correctness break in the shipped default wiring — they are latent risks
that matter precisely because this codebase is being hardened toward a service.

---

## Findings (by severity)

| # | Sev | Dimension | Location | Finding | Recommended fix |
|---|-----|-----------|----------|---------|-----------------|
| H1 | High | Correctness (science) | `viralunity/scripts/python/add_negative_control_enrichment.py:109-134,288,296-306` | `control_mean`/`SD` are computed only over controls where the taxon is **present**, not zero-filled over all controls; the z-score requires ≥2 *present* values while the z-vs-log2 decision branch keys off *total* `n_controls`. A contaminant seen in 1 of 5 controls gets `control_mean` = that single value and silently falls back to the log2 gate. This biases the contamination filter toward dropping real low-prevalence taxa. **Science call → PI sign-off.** | Decide and document the "absent-in-control = 0 vs no-data" contract; make the z-score gate and the decision branch use the same denominator. |
| H2 | High | Flexibility / Security | `viralunity/validators.py:66-105` (dead), `viralunity_meta_cli.py:40`, `viralunity_consensus_cli.py:71`, `config_generator.py:83`, `provenance.py:97` | `run_name` is written verbatim into `config["output"]` and thus into every rule's output path and the DAG, with **no validation**. `sanitize_identifier`/`ensure_within_base` exist for exactly this purpose but have **zero callers**. A `run_name` containing `../`, spaces, or glob/shell metacharacters is a path-traversal / wildcard hazard. | Call `sanitize_identifier` on `run_name` (and reuse it in the sample-sheet path); add `wildcard_constraints`. |
| H3 | High | Error handling / Rule design | 112 of 186 executing rules; e.g. `rules/consensus_nanopore.smk:7` (Clair3 `infer_consensus_sequence`), all `script:` post-processing rules, `rules/metagenomics_multiqc_illumina.smk:1` | **60% of executing rules declare no `log:`.** Nanopore consensus/variant calling (Clair3) has `benchmark:` but no `log:`; the entire pandas post-processing chain routes stderr only to Snakemake's master log, so a failed per-sample transform is hard to trace. The Illumina common module logs these rules but the Nanopore twin (`rules/consensus_nanopore_common.smk`) does not — inconsistent. | Add `log:` to executing rules, prioritising Clair3 and the `summarize_*`/`add_*`/`apply_*` chain; mirror Illumina/Nanopore. |
| M1 | Med | Flexibility | package-wide; `_orchestrator.py:56-63`, `validators.py:289-443` | **No JSON-schema config validation** (no `schemas/`, no `snakemake.utils.validate`, no `jsonschema`). The hand-written checks are broad on presence/mutual-exclusion but have **no type/range/enum checks**: negative `--threads`, `minimum_coverage=-5`, `af_threshold=5.0`, `threads_total=0` all pass validation and reach the workflow. There is no config-key contract, so key drift is caught only deep in Snakemake. | Add a JSON schema validated via `snakemake.utils.validate`; add numeric bounds. |
| M2 | Med | Correctness | `viralunity/scripts/python/convert_diamond_output_to_krona_input.py:3-4,131` | The CLI branch calls `argparse.ArgumentParser()` but the module imports only `gzip`/`os` → `NameError` on any CLI use. The Snakemake `script:` path is unaffected. Module has 0% test coverage. | `import argparse`; add a smoke test. |
| M3 | Med | Correctness | `viralunity/scripts/python/select_reference_genomes.py:239` | Summaries are read without `dtype={"taxid": str}` (unlike `add_rpkm_to_summary.py` and `filter_krona_by_pass_taxids.py`, which force it). A single NaN taxid in any concatenated summary types the whole column `float64`; taxids then stringify as `"3418604.0"` and never match the `genome2taxid` keys → **silent empty `reference_targets.tsv`** with no error. | Read with `dtype={"taxid": str}`; assert non-empty selection where expected. |
| M4 | Med | Rule design | `rules/consensus_illumina_common.smk:36-56,79-86`; `rules/metagenomics_multiqc_illumina.smk`; `organize_files` (`metagenomics_illumina.smk:167-320`) | Undeclared shell I/O breaks incremental/parallel correctness: the MultiQC rule scans a directory whose fastp inputs aren't declared; `align_consensus_to_reference_genome` writes 3 undeclared outputs and reads `*renamed.fasta` via a glob rather than its declared input; `organize_files` builds a large symlink tree and `find`s other rules' outputs while declaring only `benchmark.tsv`. | Declare real inputs/outputs (or `directory()`); extract `organize_files` logic to a tested Python `script:`. |
| M5 | Med | Correctness | `viralunity/scripts/python/annotate_diamond_taxonomy.py:66` | `mapped_reads = parts[-1]` blindly takes the last column. This is correct only because the upstream rule appends the count; fed a raw outfmt6 it silently writes `bitscore` as reads, and `summarize_krona_taxa.py:23-25` then drops all counts silently. 0% coverage. | Reference the column by name/known index; validate the header. |
| M6 | Med | Correctness | `viralunity/scripts/python/calculate_assembly_stats.py:34` | Reads are counted by `line.rstrip("\n") == "+"`. This breaks on spec-legal `+<read_id>` separators (count → 0) and miscounts single-base reads whose quality char is `+` (Phred 10). Inconsistent with `add_RPM_to_summary.py`'s `lines // 4`. Affects the reported `number_of_reads`/`number_of_trim_reads`. | Count `lines // 4` (or use the fastp JSON / seqkit). |
| M7 | Med | Reproducibility | `viralunity/scripts/envs/*.yaml`; `environment.yml`; Makefile `lock` target | Envs pin only **minor** versions (`samtools=1.23`, `kraken2=2.1`, `diamond=2.1`); `genome_selection.yaml` floats `pandas` and `python>=3.9` entirely. **`conda-lock.yml` is not committed** despite the Makefile instructing it. Results are not bit-reproducible over time. | Commit `conda-lock.yml`; pin exact patch versions or consume the lockfile at runtime; pin `pandas` in `genome_selection.yaml`. |
| M8 | Med | Testing | coverage run; `test/scripts/` | Zero unit tests for `add_RPM_to_summary`, `annotate_diamond_taxonomy`, `convert_diamond_output_to_krona_input`, `filter_diamond_by_idxstats`, `filter_taxids`, `summarize_krona_taxa` (all 0% line coverage); `select_reference_genomes` at 15%. Findings M2/M3/M5 all sit in untested modules. No `.tests/`-style end-to-end run in CI (dry-run only); no `--generate-unit-tests`. | Add unit tests for the untested metric-path scripts; keep the empirical suite but wire a tiny fixture E2E into CI. |
| M9 | Med | Error handling | `rules/metagenomics_dehost_illumina.smk:126`; `rules/consensus_illumina.smk:96,103`; `rules/alignment_illumina.smk:28-30` | Error-swallowing shell: `gzip -dc "$1" 2>/dev/null \|\| true` means a corrupt host-filtered FASTQ silently yields an empty merge; `tabix … \|\| touch {index}` masks a real tabix failure with an empty/invalid index; the alignment pipeline's `2> {log}` binds only to the last command, so `minimap2`/`samtools view` stderr is lost (pipefail still aborts the run). | Remove blanket `\|\| true`; capture per-stage stderr; fail loudly on decompression errors. |
| M10 | Med | Reproducibility | `_orchestrator.py:117-119`; `provenance.py:74` | The manifest is written **before** the run (records intent, not result), stores the config **path, not a hash/copy** (later edits silently diverge), and omits tool/DB/git versions (acknowledged in the module docstring but not collected anywhere). | Hash/copy the config into the manifest; write/patch a completion record post-run; capture tool + DB versions. |
| M11 | Med | Code quality | `mypy viralunity/` | 12 mypy errors in 8 files (implicit-Optional defaults, `Sequence[object]`→`list`, missing annotations, a `Callable` return-type mismatch). CI runs mypy `continue-on-error` (documented ~14-item backlog), so it is non-gating. | Clear the backlog and flip mypy to gating. |
| L1 | Low | Wildcards | tree-wide; `rules/metagenomics_reference_assembly.smk:50` | **No `wildcard_constraints` anywhere.** Snakemake's default wildcard regex is greedy `.+` and crosses `/`; `assembly/{ref_key}/references/{sample}.fasta` has two unconstrained wildcards. **Latent, not active:** dry-runs resolved cleanly with no ambiguity warnings, and sample names are regex-sanitized at parse time (`^[A-Za-z0-9][A-Za-z0-9._-]*$`, no `/`). Defense-in-depth. | Add `wildcard_constraints: sample=r"[^/]+", ref_key=r"[^/]+", segment=r"[^/]+"`. |
| L2 | Low | Hygiene | repo root; `.gitignore` | Large untracked, **un-ignored** working dirs risk accidental commit: `revisa5_run/` (**12G**), `results/` (558M), `my_test_chikv/` (290M), `my_test_data.tar.gz` (51M). `.gitignore` covers `my_results/`/`scratch/` but not these. Loose working `.md`/`.csv` files sit in the root. | Extend `.gitignore` (`results/`, `revisa5_run/`, `*.tar.gz`, `my_test_*`); move scratch docs out of the tree. |
| L3 | Low | Maintainability | rule wiring; Makefile | Minor defects: `perform_qc` declares `conda:` twice (`qc_illumina.smk:9,37`); `validators.py:377`'s `.gz` sub-clause is logically redundant (reduces to `not os.path.isfile(taxids)`); `collect_reference_assemblies` is the one executing rule with no `conda:` (`metagenomics_reference_assembly.smk`); the Makefile `run-consensus` target is a copy-paste of `run-meta` (passes kraken2/krona args to consensus; `exmaple` typo); `rpm_floor=1.0`/`min_mapped=1` are hardcoded in params rather than config. | Housekeeping. |
| L4 | Low | Provenance / robustness | `viralunity_create_samplesheet.py:107-118,263-269`; `summarize_krona_taxa.py:40`; `filter_taxids.py:56-57` | `create-samplesheet` doesn't sanitize names and writes CSV via f-strings (a comma in a path corrupts the row); `summarize_krona_taxa` does an unguarded 2-field unpack; `filter_taxids` treats `taxid_column` as a position into a sliced column Index (works by coincidence for the shipped configs). | Use `csv.writer` + `sanitize_identifier`; guard the unpack; clarify the column-index semantics. |

---

## Per-dimension notes

**1. Correctness & scientific validity.** No coordinate-system (BED/GFF/VCF/SAM) interval math exists
in the Python layer — the highest-risk category is empty (BLAST `qstart/qend/sstart/send` in
`filter_top_species_hits.py` are passed through verbatim, never used in arithmetic). Verified-correct
by direct read: RPM (`mapped/total*1e6`, zero-denominator guarded, with a unit-consistent denominator =
the merged host-filtered FASTQ that also feeds classification), RPKM (`rpm*1000/genome_length_bp`, taxid
forced to `str` on both sides of the merge), the bleed filter, `bedtools genomecov -d` coverage stats,
the lineage walk (root-terminating and cycle-safe), the ICTV token match (splits on `,`/`;` to avoid the
`invertebrate ⊃ vertebrate` substring trap), and `filter_krona`'s boolean handling. Defects concentrate
in the untested metric-path scripts (H1, M2, M3, M5, M6).

**2. Snakemake rule design.** `--lint` findings are style-only (`+` path composition; "mixed rules and
functions") plus **false-positive** "absolute path" warnings that misread a `sed 's/[\/|,~ ]/_/g'`
expression (`consensus_nanopore.smk:13`) and a `"/" + track` mid-path concatenation
(`metagenomics_illumina.smk:49`). No `run:` directives anywhere (good — logic is in `shell:`/`script:`).
Every rule but one references a per-rule `conda:` env. The single checkpoint, `select_references_meta`,
is structurally correct: its aggregator calls `checkpoints.select_references_meta.get()`, guards
empty/missing with `EmptyDataError`, and de-dupes with a `set()`; the taxonomic tracks correctly use
static `expand(...)`, so the data-dependent DAG is confined (appropriately) to reference assembly. Main
weaknesses: log coverage (H3), undeclared I/O (M4), and embedded awk/sed that belongs in tested scripts.

**3. Wildcards.** No `wildcard_constraints` and no `ruleorder`. Mitigated in practice by sample-name
sanitization plus clean DAG resolution; recommend adding constraints as defense-in-depth (L1).

**4. Flexibility & configurability.** Samples are driven by a sample sheet and thresholds/paths by the
config — good. But there is no schema validation and no numeric bounds (M1), and some relative DB-path
defaults resolve against the current working directory.

**5. Python code quality.** ruff and black are clean. The code is idiomatic, uses streaming parsers, and
has no `shell=True`. Weak spots: 12 mypy errors (M11); a few oversized functions
(`validate_metagenomics_requirements` ~150 lines; `config_generator.add_metagenomics_settings` ~40
parameters; `viralunity_get_databases_cli.py` at 862 LOC); and the enrichment script uses
`iterrows`/`apply` (three row-wise passes) — fine at taxa-table scale, but a vectorization opportunity.

**6. Error handling & robustness.** A strong orchestration skeleton: typed `ViralUnityError` with stable
codes, clean user-facing messages versus `logger.exception` for the genuinely unexpected, and
best-effort provenance that never aborts an analysis. Gaps: no `onerror`/`onsuccess`; missing per-rule
logs (H3); error-swallowing shell one-liners (M9).

**7. Edge cases.** Handled: empty FASTQ (`n_lines==0 → 0`), zero-denominator RPM, empty consensus (a mock
VCF is emitted), empty idxstats/summary (touch / empty table), and malformed `.fai`/genome2taxid lines
(warn + skip). Not handled or fragile: NaN taxid (M3), non-standard FASTQ `+` lines (M6), input-format
drift (M5), and paths with spaces in the unquoted Illumina shell rules and `organize_files` `ln -sf`.

**8. Testing & coverage.** 373 tests with real fixture-based assertions (e.g. taxdump-driven harmonize
tests) — not superficial; ~712 asserts across ~5,900 test LOC. Overall 55% branch coverage;
orchestration/CLI/validators are well covered, but 6 metric scripts sit at 0% and `select_reference_genomes`
at 15% (M8). Good dry-run matrix (13 configs) plus an opt-in empirical suite (SARS-CoV-2 with a real
download); no in-CI end-to-end run.

**9. Reproducibility.** Per-rule conda envs with `--use-conda`, a docker image (non-root, pinned base),
a provenance manifest with per-input SHA256, and PyPI trusted-publishing. Gaps: minor-only version pins,
a floating `pandas`, and an **uncommitted lockfile** (M7); the config is recorded by path not hash and
the manifest is written pre-run (M10); no `--use-singularity` path.

**10. Data & performance.** `benchmark:` is on all heavy rules (excellent); `threads:` and
`resources: mem_mb` are on the compute-heavy rules. `temp()` is used but inconsistently — host-filtered
FASTQs, MEGAHIT contigs, Kraken2/DIAMOND tables, and viral BAMs are not `temp()`, and
`extract_viral_contigs.ids` is `temp()` in Illumina but persistent in Nanopore. **No `protected()`** on
final deliverables, so a stray `--forceall` overwrites published results.

**11. Documentation.** Strong: a Sphinx site (install / usage / commands / architecture / output /
tutorials), ReadTheDocs, a CHANGELOG (Keep-a-Changelog + SemVer), RELEASING.md, CONTRIBUTING, and
issue/PR templates. No Snakemake `report()` HTML report. The Makefile `run-consensus` example is wrong
(L3).

**12. Maintainability & hygiene.** No TODO/FIXME, no secrets, and no absolute/user-specific paths in the
package. The main issue is working-tree clutter not covered by `.gitignore` (L2) plus minor dead code
(L3/L4).

---

## Quick wins (high impact, low effort)

- `import argparse` in `convert_diamond_output_to_krona_input.py` (M2).
- Add `dtype={"taxid": str}` to `select_reference_genomes.py:239` (M3).
- Call `sanitize_identifier(run_name)` — the helper already exists (H2).
- Add `wildcard_constraints` (L1); commit `conda-lock.yml` and pin `pandas` in `genome_selection.yaml` (M7).
- Fix read counting to `lines // 4` in `calculate_assembly_stats.py` (M6).
- Extend `.gitignore`; remove the double `conda:`; fix the Makefile `run-consensus` target (L2/L3).

## Structural (larger investment)

- A JSON schema with numeric-bounds validation via `snakemake.utils.validate` (M1).
- Systematic `log:` coverage plus `onerror`/`onsuccess` handlers (H3).
- Unit tests for the 6 untested metric scripts, and a tiny in-CI end-to-end fixture (M8).
- Resolve H1 (negative-control statistics) with the PI, then encode the decided contract and add tests.
- Declare rule I/O honestly / extract `organize_files` to a tested script (M4).
- Clear the mypy backlog and gate on it (M11).

## Genuinely well done (preserve in any refactor)

- The safe subprocess layer (list-form args, timeouts, no `shell=True`), the zip-slip guard, and the
  guarded `rmtree`.
- The typed/coded `ViralUnityError` hierarchy and the clean orchestrator error routing.
- Strong sample-sheet parsing (sanitized, deduplicated, paired/single validated, exact column counts).
- The provenance manifest with per-input SHA256; correct checkpoint usage; broad `benchmark:` coverage.
- A mature CI (lint + format + mypy + coverage + dry-run + security + docker + conda-smoke) and real
  fixture-based tests.
- Verified-correct core metrics (RPM/RPKM/bleed/coverage/lineage/ICTV) and a deliberately
  coordinate-math-free design.
