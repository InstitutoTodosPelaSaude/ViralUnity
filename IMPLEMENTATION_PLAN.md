# IMPLEMENTATION_PLAN — v1.3.2 metagenomics filter refinements

Branch: `feature/filters-v1.3.2` (off `main` @ 1.3.1). Tracks eight PI-approved
refinements to the metagenomics contamination filters, motivated by the REVISA
contamination episode (mayaro + HIV libraries seen across negative controls).

## Decisions needed / locked (PI sign-off obtained in planning)

All eight scientific/output decisions were signed off by the PI before coding. One
**tunable remains open for the PI to confirm during verification** (non-blocking).

| # | Proposal | Decision |
|---|----------|----------|
| 1 | Bleed filter metric | **RPKM when available (per group) else RPM**, recalibrated floor. *NB the pass/fail is mathematically identical to RPM (within-taxon ratio); only the floor gate changes.* |
| 2 | Aggregate negative-control filter | **Pool RAW reads then normalize** (`Σreads/Σtotal_reads·1e6`); complementary to the z-score, does not replace it. |
| 3 | Pass-flag columns | `fold_enrichment_10x_pass`, `fold_enrichment_100x_pass`, `neg_pass_5`, `neg_pass_10`; `>=` inclusive; **NA stays NA**. |
| 4 | log2 → log10 | Rename everywhere; default log10-ratio threshold **1.0 (=10× fold)** — accepted as ~5× stricter than the old log2 default of 1.0 (=2×). |
| 5 | `final_species` | `coalesce(nr_correct_species, name)` for all rows, immediately right of `nr_correct_species`. |
| 6 | Contig filter | `largest_contig_bp` + `largest_contig_median_depth` (true `samtools depth`, median over the single largest contig for the taxon); gated on `--viral-genomes`. **Updated 2026-07-11 (PI reversed the earlier diamond-only scope): now on BOTH contig tracks** — kraken2_contigs extracts its viral contigs by lineage (taxid 10239) and remaps reads for its own depth BAM; diamond_contigs reuses its existing viral remap. |
| 7 | Output reorg | **Clean break** to a per-method/per-rank tree; combined user-facing file removed; family/genus names propagated down; `select_reference_genomes.py` + krona filter rewired. |
| 8 | Remove `source` column | Confirmed no consumer; drop it. |

### OPEN tunable for PI (confirm at toy/REVISA verification)

**RPKM bleed floor default.** RPKM ≈ RPM·(1000/len), so a 10 kb genome's RPKM is ~0.1×
its RPM. New config key `bleed_rpkm_floor` is introduced with a **proposed default of
`0.1`** (≈ the current RPM floor of 1.0 evaluated at 10 kb). The RPM path keeps
`bleed_rpm_floor = 1.0`. Please confirm or retune `bleed_rpkm_floor` once you see the toy
and REVISA outputs — it only affects which taxa the bleed filter is *applied* to, never
the ratio test itself.

## Commit map

Each commit leaves `make test` + `make test-dryrun` + `make lint` green.

- **C1 — remove `source` (P8).** `summarize_krona_taxa.py` + test + `docs/output.md`.
- **C2 — log2 → log10 (P4).** `add_negative_control_enrichment.py` (fn, column,
  `neg_decision` strings, threshold-used column), CLI `--log10-ratio-threshold`,
  `ConfigKeys`, `config_generator`, `viralunity_meta`, `validators`, all track `.smk`
  `config.get`, tests, docs. Default 1.0 = 10×.
- **C3 — pass-flag columns (P3).** `add_negative_control_enrichment.py` + test.
- **C4 — `final_species` (P5).** `harmonize_nr_summary.py` (+ reads-track summary so the
  column exists on every track) + test.
- **C5 — bleed metric + floor (P1).** `apply_max_rpm_bleed_filter.py` (per-group metric
  select, `bleed_metric` col, rpm/rpkm floors) + 6-touch wiring for `bleed_rpkm_floor` +
  **new `apply_max_rpm_bleed_filter_test.py`** + track `.smk` params.
- **C6 — aggregate pooled neg filter (P2).** `add_negative_control_enrichment.py`
  (`pooled_control_metric`, `agg_fold_enrichment`, `agg_log10_ratio`, agg pass flags) +
  test with REVISA HIV fixture `[276,1988,1991,17144]`.
- **C7 — contig size + median depth (P6).** New `depth_of_viral_contigs` rule
  (`samtools depth`) in both contig tracks × both platforms (gated `if compute_rpkm`) +
  new `add_contig_stats_to_summary.py` + test + dryrun config.
- **C8 — output reorg clean break (P7).** New `split_summary_by_rank.py` + terminal
  rules; `_all_inputs`/`organize_files` retargeted; `select_reference_genomes.py` reads
  species file; `filter_krona_by_pass_taxids` reads the combined intermediate; docs +
  dryrun configs + test.
- **C9 — release prep.** `__version__` + `Dockerfile` LABEL → 1.3.2; `CHANGELOG.md`
  `[1.3.2]`; finalize docs.

## Output-change summary (for the PI's records)

| Change | Alters existing outputs? |
|---|---|
| Remove `source` | Yes — column dropped |
| log2→log10 + stricter default | Yes — column renamed; `neg_pass` decisions shift (2×→10×) |
| Pass-flag cols | Additive |
| `final_species` | Additive |
| Bleed metric/floor | pass/fail unchanged; `bleed_metric` added; floor gate may change which taxa are filtered |
| Aggregate neg filter | Additive; `neg_pass` unchanged |
| Contig size/depth | Additive (only with `--viral-genomes`) |
| Rank reorg | **Breaking** — layout replaced by per-rank tree |

## Verification & operational gates (hard constraints)

1. Work only on `feature/filters-v1.3.2`; never commit to `main`.
2. **Do NOT launch/re-run/trigger any REVISA analysis** — REVISA run 5 is still running.
   Coding + reading existing outputs + unit tests only, until the PI confirms run 5 done.
3. After run 5 finishes **and PI confirms**: `make test`, `make test-dryrun`, `make lint`,
   then the **toy sars-cov-2** end-to-end run; report real results.
4. Only after the toy run passes **and the PI explicitly says go**: run ONLY the latest
   summarization steps on real REVISA data. Confirm scope first; never overwrite existing
   results without confirmation.
