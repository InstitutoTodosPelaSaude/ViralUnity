# Handoff — v1.3.2 metagenomics filter refinements

**Branch:** `feature/filters-v1.3.2` (branched off `origin/main` @ `840ca31`, which is the
merged 1.3.1 release). **Version bumped to `1.3.2`.** Not pushed / not PR'd / not tagged.

This document is a continuation point: if the session breaks, everything needed to resume,
verify, and finish the release is here. Companion docs on the branch: `IMPLEMENTATION_PLAN.md`
(decision record) and the planning file `/home/gevop/.claude/plans/sharded-hugging-acorn.md`.
Source spec: `metagenomcis_filters_issue.md` (repo root, untracked, do not delete).

## Why this work exists

The REVISA project (`/home/gevop/projects/REVISA`) ran many metagenomic runs of clinical
samples through ViralUnity. Its negative controls revealed the lab was contaminated with
**mayaro** and **HIV** libraries from other projects, giving an empirical testbed for the
metagenomics contamination filters. The PI proposed eight refinements; this branch
implements them as v1.3.2. Empirical anchor found during planning: HIV
(*Lentivirus humimdef1*) shows negative-control RPKM `[276, 1988, 1991, 17144]` — a 62×
spread — in `REVISA/selected_sample_outputs/REVISA1_diamond_contigs_taxa_summary_RPKM.nr.bleed.neg.tsv`.

## Status

- **All 8 refinements implemented + 1 later extension (kraken2 contig stats).**
- **466 unit tests + 14 Snakemake dryruns pass; `black`/`ruff` clean.**
- Working tree clean. **Nothing pushed/tagged** — tagging triggers the PyPI publish workflow,
  which is the maintainer's call (see `RELEASING.md`).

## Hard operational constraints (still in force)

1. Work only on `feature/filters-v1.3.2`; never commit to `main`.
2. **Do NOT launch/re-run/trigger any REVISA analysis** until the PI confirms **REVISA run 5**
   has finished. Coding, reading existing outputs, and unit/dryrun tests are fine.
3. After run 5 finishes **and the PI confirms**: run the test suite, then the **toy
   sars-cov-2** end-to-end run; report real results.
4. Only after the toy run passes **and the PI explicitly says go**: run ONLY the latest
   summarization steps on the real REVISA data (heavy upstream compute already done — treat
   it as precious/read-mostly; never overwrite existing results without confirming).

## Commits on this branch (oldest → newest)

| SHA | Commit | What |
|-----|--------|------|
| `2aa9af4` | docs: implementation plan | `IMPLEMENTATION_PLAN.md` (decisions + commit map) |
| `fcb129e` | remove `source` column | dropped from `summarize_krona_taxa.py` |
| `72fcd72` | log2 → log10 (**breaking**) | rename everywhere; default `1.0` now = **10× fold** (was 2×), i.e. stricter neg gate |
| `6e5e025` | pass-flag columns | `fold_enrichment_10x/100x_pass`, `neg_pass_5/10` (`>=`, NA→NA) |
| `f0bd78e` | `final_species` | coalesce(`nr_correct_species`, `name`) in NR harmonize |
| `4514514` | bleed metric-aware | RPKM when available else RPM; `bleed_metric`; `max_rpm`→`bleed_max`; metric-specific floors |
| `41b075f` | bleed unit test | new `apply_max_rpm_bleed_filter_test.py` (was missing) |
| `9f69b31` | aggregate neg-control filter | pooled raw-read control; `pooled_control_metric`, `agg_*` cols |
| `44118c3` | contig stats (diamond_contigs) | `largest_contig_bp` + `largest_contig_median_depth`; `.ctgstats` step |
| `5e3c319` | per-rank split (**breaking**) | user-facing `<track>/{family,genus,species}/`; chain moved to `<track>/chain/` |
| `143be79` | release: 1.3.2 | `__version__` + Dockerfile LABEL + CHANGELOG |
| `b99df7a` | extend contig stats to kraken2_contigs | PI reversed the diamond-only scope; both contig tracks now |

## The eight refinements (+ extension) and the decisions behind them

Numbers below are the proposal numbers from `metagenomcis_filters_issue.md`. Each scientific
fork was signed off by the PI (memory: `scientific-changes-need-user-signoff`).

1. **Bleed filter → RPKM-aware.** *Key finding to remember:* the bleed test is a
   **within-taxon** ratio across samples (`value >= fraction * max`, grouped by
   `tool,mode,rank,taxid`). Genome length is constant within a taxon, so RPM→RPKM is a
   uniform rescale and **`bleed_pass` is unchanged**. The only real effect is the
   *application floor*. So: metric auto-selected per group (`bleed_metric` column), floors
   are metric-specific — `bleed_rpm_floor` (1.0) and **`bleed_rpkm_floor` (0.1, OPEN — see
   below)**, both YAML-only knobs. Column `max_rpm` renamed `bleed_max`.
2. **Aggregate (pooled) negative-control filter.** Complementary to the z-score (does NOT
   change `neg_pass`). Pools **raw reads** across controls: `pooled = Σ(metric·total_reads)/
   Σ(total_reads)`, so a deeply-sequenced control counts more (falls back to equal weights if
   `total_reads` absent). New cols: `pooled_control_metric`, `agg_fold_enrichment`,
   `agg_log10_ratio`, `agg_fold_enrichment_10x/100x_pass`. Catches the HIV-style high-variance
   contamination that inflates per-control SD and drags the z-score down.
3. **Pass-flag columns.** `fold_enrichment_10x_pass`, `fold_enrichment_100x_pass`,
   `neg_pass_5`, `neg_pass_10`. `>=` inclusive; **NA stays NA** (matches the neg_pass
   NA-as-keep contract). Nullable-boolean.
4. **log2 → log10.** Interpretability. Column `log2_ratio`→`log10_ratio`, `neg_decision`
   values `log10_ratio`/`log10_ratio_fallback`, `log10_ratio_threshold_used`, CLI
   `--log10-ratio-threshold`, `ConfigKeys.LOG10_RATIO_THRESHOLD`. **Default `1.0` now means a
   10-fold enrichment** (log10=1) vs the old 2-fold (log2=1) — a deliberate, PI-approved
   stricter default.
5. **`final_species`.** `coalesce(nr_correct_species, name)` — NR's corrected species when NR
   confidently disagreed, else the original name. Positioned right of `nr_correct_species` in
   NR harmonize; the per-rank split (proposal 7) guarantees it on non-NR tracks too
   (`= name`).
6. **Contig stats.** `largest_contig_bp` + `largest_contig_median_depth` per taxon, via a
   `.ctgstats` chain step. **Now on BOTH contig tracks** (PI reversed the initial
   diamond-only scope). Each track remaps host-filtered reads to *its own* viral contigs and
   runs `samtools depth -a`; median is over the single largest contig assigned to the taxon
   (lineage-climb assignment). Gated on `--viral-genomes` (`compute_rpkm`). *Caveat documented
   in the tutorial:* contig length is a proxy for, not a measure of, covered genome fraction.
7. **Per-rank output reorg (breaking).** User-facing deliverable is now one table per rank
   under `<track>/{family,genus,species}/`, with higher-rank names propagated down (species
   gains `family`+`genus`; genus gains `family`). The combined cumulative chain is computed
   internally under `<track>/chain/`. `select_reference_genomes.py` reads `<track>/chain/`
   (with a fallback to the old flat layout); the Krona filter reads the combined chain file
   unchanged.
8. **Remove `source` column** (was an internal krona-input path; no consumer).

## New Python scripts (all in `viralunity/scripts/python/`, each with a `*_test.py`)

- `split_summary_by_rank.py` — terminal per-rank split; propagates family/genus name columns;
  guarantees `final_species`.
- `add_contig_stats_to_summary.py` — per-taxon largest-contig length + median depth from a
  `samtools depth -a` file + contig→taxid map; classifier-agnostic (shared by both tracks).
- `extract_viral_contigs.py` — selects contigs a classifier called viral (lineage under
  taxid `10239`) so kraken2_contigs can build a light viral-only remap; also has a FASTA
  writer. Used only by the kraken2_contigs track.

## Modified Python / config surface

- `apply_max_rpm_bleed_filter.py` — metric selection, `bleed_metric`, `bleed_max`, rpm/rpkm floors.
- `add_negative_control_enrichment.py` — log10 rename, pass flags (`_threshold_flag`,
  `_add_pass_flag_columns`), aggregate pooled control (`_add_aggregate_control`).
- `harmonize_nr_summary.py` — `_add_final_species`.
- `summarize_krona_taxa.py` — dropped `source`.
- `select_reference_genomes.py` — `resolve_summary_file` now globs `<track>/chain/` first,
  falls back to the flat layout.
- CLI/config: `viralunity_meta_cli.py` (`--log10-ratio-threshold`), `constants.py`
  (`LOG10_RATIO_THRESHOLD`), `config_generator.py`, `viralunity_meta.py`, `validators.py`.
- `viralunity/__init__.py` (1.3.2), `Dockerfile` LABEL (1.3.2).

## Snakemake workflow surface

- Top-level `metagenomics_{illumina,nanopore}.smk`:
  - `_summary_stem` now points into `<track>/chain/`.
  - `_chain_steps` adds `ctgstats` for any `track.endswith("contigs")` when `compute_rpkm`
    (order: `[nr] → [ctgstats] → bleed → [neg] → [ictv]`).
  - Added `per_rank_summary`/`per_rank_summaries`/`_chain_tail`; `_all_inputs()` requests the
    per-rank files (not the combined `final_summary`).
  - Four `split_<track>_by_rank` rules after the includes (guarded by the track flags).
- `rules/metagenomics_diamond_contigs_{illumina,nanopore}.smk`: `depth_of_viral_contigs`
  (samtools depth on the existing viral BAM) + `add_contig_stats_diamond_contigs`.
- `rules/metagenomics_kraken2_contigs_{illumina,nanopore}.smk`: `extract_viral_contigs_kraken2`
  + `remap_and_depth_viral_contigs_kraken2` (writes to `mapping/viral_kraken2/`) +
  `add_contig_stats_kraken2_contigs`.
- All 8 per-track bleed rules gained a `rpkm_floor = config.get("bleed_rpkm_floor", 0.1)` param.
- Chain-base paths in all 8 track rule files relocated under `<track>/chain/`.

## New / changed output schema (for the PI's records)

Fully-loaded diamond_contigs deliverable (species table):
`taxonomic_assignments/diamond_contigs/species/diamond_contigs_species_taxa_summary_RPKM.nr.ctgstats.bleed.neg.ictv.tsv`

New columns vs 1.3.1: `family`, `genus` (propagated, rank-dependent), `final_species`,
`bleed_metric`, `bleed_max` (was `max_rpm`), `fold_enrichment_10x_pass`,
`fold_enrichment_100x_pass`, `neg_pass_5`, `neg_pass_10`, `pooled_control_metric`,
`agg_fold_enrichment`, `agg_log10_ratio`, `agg_fold_enrichment_10x_pass`,
`agg_fold_enrichment_100x_pass`, `largest_contig_bp`, `largest_contig_median_depth`.
Renamed: `log2_ratio`→`log10_ratio`, `log2_ratio_threshold_used`→`log10_ratio_threshold_used`.
Removed: `source`. Layout: combined chain moved to `<track>/chain/`; per-rank tree is new.

## OPEN item needing the PI (non-blocking)

**`bleed_rpkm_floor` default = `0.1`** is a guess (≈ the old RPM floor of 1.0 evaluated at a
10 kb genome, since RPKM ≈ RPM·1000/len). It only changes *which taxa the bleed filter is
applied to*, never the ratio test. Confirm/retune once the toy + REVISA outputs are visible.
It's a YAML-only knob — edit the generated config and rerun Snakemake.

## How to verify (run these when unblocked)

```bash
conda activate viralunity            # env must have the package installed (pip install -e .)
python -m pytest test/ -q            # 466 unit tests + 14 dryruns
# lint (Makefile's `make lint` reinstalls and hits a py3.12 pin; run the tools directly):
black --check viralunity/ test/ && ruff check viralunity/ test/
```

Dryrun coverage note: `test/dryrun_configs/metagenomics_illumina__fullchain.yaml` exercises
both contig tracks + RPKM + NR + neg + ICTV (so `.ctgstats` + per-rank split on both);
`metagenomics_nanopore__ctgstats.yaml` covers the nanopore contig tracks.

### Phase 5 — toy sars-cov-2 (after REVISA run 5 done + PI confirms)
Run the toy dataset end-to-end (see `CLAUDE.md` invocation examples / `docs/tutorial/`),
confirm the per-rank tree and new columns are produced, and report real results.

### Phase 6 — REVISA summarization (after Phase 5 passes + explicit PI go)
Run ONLY the latest summarization steps on the real REVISA data (heavy compute already done).
Confirm scope first; do not re-trigger upstream compute; do not overwrite existing results
without confirming.

## Finishing the release (maintainer)

Per `RELEASING.md`: `__version__` and Dockerfile LABEL are already at 1.3.2 and `CHANGELOG.md`
has the `[1.3.2]` entry. Remaining: open a PR / merge to `main`, then tag `v1.3.2` and push the
tag (pushing the tag triggers the PyPI publish workflow). Left intentionally to the maintainer.

## Gotchas for whoever continues

- `make lint` runs `make install-dev` first, which fails on Python 3.12 (`requires <3.12`).
  Run `black`/`ruff` directly instead; the installed env already works for tests.
- Tests are stdlib `unittest` classes collected by `pytest`; run via
  `python -m pytest test/...`, not `python -m unittest` from the repo root (package-path issues).
- With both contig tracks on, each does its own viral remap → two (light, viral-only) remaps
  per sample. Intentional: keeps per-track depth semantics correct.
- Do not rename `validate_args` / `generate_config_file` / `run_snakemake_workflow`
  (tests patch those names) — see `CLAUDE.md`.
