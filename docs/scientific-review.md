# Scientific-behavior review (for PI sign-off)

This document lists candidate changes to ViralUnity's **scientific behavior** — places
where the pipeline's numerical output (consensus calls, abundance metrics, contamination
filtering, reference selection) may be wrong, inconsistent across platforms, or
non-reproducible. **Nothing here has been changed in code.** Each item is a proposal for
the maintainer (PI) to accept, reject, or adjust in a later session, because changing any
of these alters reported results and therefore needs domain sign-off for reproducibility
and defensibility in a public-health / genomic-epidemiology setting.

All file:line references were verified against `main` (tip `972516b`) on 2026-07-06.

For each item: **Current behavior**, **Proposed change**, **Rationale**, **Public-health
impact**, **Suggested test** to add when the change is made.

Severity legend: 🔴 affects reported genomes/metrics · 🟠 cross-platform comparability /
reproducibility · 🟡 QC/reporting only.

---

## 1. 🟠 Nanopore vs Illumina low-coverage masking is off-by-one

**Current behavior.** Nanopore masks a position to `N` when depth `<= minimum_depth`
(`viralunity/scripts/rules/consensus_nanopore.smk:61-62`, awk `$3 <= int(minimum_depth)`,
default `minimum_depth=10`). Illumina uses `samtools consensus -d {minimum_depth}`
(`viralunity/scripts/rules/consensus_illumina.smk:66`), where `-d` is an *inclusive*
minimum — a position at depth exactly `minimum_depth` is **called**. So at 10× the
nanopore consensus gets an `N` and the Illumina consensus gets a base, from the same
nominal setting.

**Proposed change.** Change the nanopore comparison to `$3 < int(minimum_depth)` so both
platforms mean "call sites with depth ≥ `minimum_depth`".

**Rationale.** "Sites ≥ minimum_depth are reported" should be one definition, not two.

**Public-health impact.** N-count / genome completeness and any depth-boundary
lineage-defining mutation can differ purely by sequencing platform, undermining
cross-platform comparability in surveillance datasets.

**Suggested test.** Synthetic BAM with a position at exactly `minimum_depth`; assert it is
called (not `N`) on both nanopore and Illumina paths.

---

## 2. 🔴 Paired-end DIAMOND-reads RPM/RPKM denominator is ~2× too large

**Current behavior.** The read-count denominator is the **merged R1+R2** host-filtered
FASTQ (`viralunity/scripts/metagenomics_illumina.smk:12-14`, `get_sample_to_fastq()` →
`*.merged.fastq.gz`; merge concatenates R1 then R2 with no renaming,
`metagenomics_dehost_illumina.smk:109-130`), counted as `lines // 4`
(`add_RPM_to_summary.py:12-33`) — i.e. **per-read ≈ 2× fragments**. For **DIAMOND reads**,
the numerator `count` collapses R1/R2 mates to one id (see Item 3), i.e. **per-fragment**.
Numerator per-fragment ÷ denominator per-read → DIAMOND-reads RPM/RPKM is systematically
~½ the true per-fragment value on paired data.

> Note: **kraken2-reads is NOT affected** — kraken2 runs R1 and R2 as two files without
> `--paired` (`metagenomics_kraken2_reads_illumina.smk:32-34`), emitting one line per read,
> so its `count` is per-read and matches the per-read denominator.

**Proposed change.** Make numerator and denominator use the same unit. Options: (a) count
fragments for the denominator on paired data (e.g. count R1 only, or reads//2), or
(b) count reads (not fragments) for the DIAMOND-reads numerator. Whichever is chosen,
apply it consistently to kraken2-reads too so both classifiers are comparable.

**Rationale.** RPM/RPKM is the primary abundance metric and feeds the bleed filter,
negative-control gate, and reference selection; a hidden per-classifier factor makes
kraken2 and diamond RPMs non-comparable and shifts values across fixed thresholds.

**Public-health impact.** A true positive can drop below a fixed RPM threshold for DIAMOND
but not kraken2 (or vice-versa) on the same sample.

**Suggested test.** Tiny paired FASTQ (known fragment count) + stub diamond/kraken2 hits;
assert both classifiers report the same RPM for the same fragment set.

---

## 3. 🔴 DIAMOND read-id counting collapses mate pairs

**Current behavior.** `viralunity/scripts/python/convert_diamond_output_to_krona_input.py`
stores hits in `read2taxid[read_id]` (dict, line 73) keyed on `qseqid` (line 59), and
builds the read set `ids` (a `set()`, line 79) from `header.lstrip("@").split()[0]`
(lines 84-87). Illumina R1/R2 mates share that token (the `1:N:…`/`2:N:…` flag is after
the space), so both mates collapse to one id; output is one line per unique id.

**Proposed change.** Decide the intended counting unit (fragment vs read) and make it
explicit; if per-read is intended, disambiguate mates (e.g. append `/1`,`/2`). Coordinate
with Item 2.

**Rationale.** A fragment where one mate hits and a fragment where both mates hit currently
count identically, biasing relative abundance; and this is the mechanism behind Item 2.

**Public-health impact.** Skews relative abundance and any mate-coverage-based confidence.

**Suggested test.** Diamond output with both mates hitting one fragment and only one mate
hitting another; assert the resulting counts match the chosen unit definition.

---

## 4. 🟡 Two disagreeing FASTQ read counters

**Current behavior.** `calculate_assembly_stats.py:24-36` counts records by matching the
bare separator line `line.rstrip("\n") == "+"`. Elsewhere reads are counted as `lines//4`
(`add_RPM_to_summary.py:33`) or by header modulo
(`convert_diamond_output_to_krona_input.py`). The `"+"` method miscounts when the separator
carries the id (`+READID` → undercount) or when a 1-base read's quality line is literally
`"+"` (Phred 10 → overcount).

**Proposed change.** Replace the `"+"`-matching counter with `lines // 4` (with the
not-divisible-by-4 warning already used in `add_RPM_to_summary.py`) and share one helper.

**Rationale.** One read-count definition across the codebase.

**Public-health impact.** QC/reporting only — corrupts the per-sample read-tracking numbers
reviewers use to judge sample quality/dropout, but does not change consensus or abundance.

**Suggested test.** FASTQ fixtures with `+READID` separators and a 1-base `+`-quality read;
assert the counter returns the true record count.

---

## 5. 🟠 Cross-platform consensus allele-frequency thresholds differ by default

**Current behavior.** Illumina consensus AF default `0.5` (`consensus_illumina.smk:60,66`,
`samtools consensus -c 0.5`); nanopore consensus AF default `0.7`
(`consensus_nanopore.smk:25,58`, `bcftools filter FORMAT/AF >= 0.7`). Illumina iSNV band is
`0.05 ≤ AF < 0.5` with the **upper 0.5 hardcoded** (`consensus_illumina.smk:21,46`);
nanopore has no iSNV rule. Nanopore variants also gate on defaults `variant_depth=5`,
`variant_quality=20`, `minimum_map_quality=20` (`consensus_nanopore.smk:26-30`).

**Proposed change.** Decide whether the platforms should share a default AF (and document
why they differ if intentional — nanopore's higher error rate is a legitimate reason).
Either align defaults or document the rationale prominently.

**Rationale / impact.** A variant at ~55% frequency is fixed into an Illumina consensus but
left as reference on nanopore, so the same specimen can yield different consensus genomes /
lineage assignments by platform alone. (This may be *intended* given error profiles — hence
sign-off, not an automatic fix.)

**Suggested test.** Document-only, or a fixture asserting the documented per-platform
threshold is applied.

---

## 6. 🟠 Negative-control contamination filter goes near-inert on an empty control

**Current behavior.** The Poisson filter (`apply_negative_background_filter.py`) was
replaced by an enrichment/z-score gate (`add_negative_control_enrichment.py`, wired e.g.
`metagenomics_kraken2_reads_illumina.smk:157-170`). It is numerically safe — all
denominators use a pseudocount (`calculate_fold_enrichment:61`, `calculate_log2_ratio:70`),
and `calculate_z_score` returns `None` when `n<2` or `sd==0` (lines 86-92), so there is **no
division-by-zero or inf**. But when a negative control has ~0 reads, `control_mean≈0` and the
gate falls back to `log2((sample+1)/(0+1)) >= 1`, i.e. a taxon passes when
`sample_metric >= 1` — nearly everything passes. With **zero** controls, `neg_pass = NA`
(kept). So a failed/empty negative control provides essentially no contamination protection
rather than raising an error.

**Proposed change.** Emit a loud warning (and optionally fail/quarantine the run) when a
declared negative control has read count below a minimum, instead of silently degrading to a
permissive gate.

**Rationale.** An empty negative is a common clean-lab outcome *and* a common failure mode;
silently trusting it is the dangerous case.

**Public-health impact.** Contamination may pass unfiltered into reported taxa when the
control that was supposed to catch it is itself empty.

**Suggested test.** Run the enrichment gate with a zero-read control; assert a warning is
raised (and, if adopted, that the run is flagged rather than passing everything).

---

## 7. 🟠 Non-deterministic reference selection (and unbounded blastn)

**Current behavior.** `select_reference_genomes.py:367-377` runs
`blastn … -max_target_seqs 10` and keeps the first qualifying hit per contig relying on
output order (lines 389-433). `-max_target_seqs` is the well-known non-deterministic /
early-termination heuristic — the top N are not guaranteed to be the global best hits — so
the chosen reference can change across runs / DB versions. The blastn call
(`subprocess.run`, line 381) has **no `timeout`**. Separately,
`convert_diamond_output_to_krona_input.py:79,101` iterates a `set()` to write krona-input
rows, so output row order varies run-to-run (counts unaffected; file checksums not stable).

**Proposed change.** For reference selection, sort blastn hits explicitly by
`bitscore`/`pident`/`qcov` and pick deterministically (and consider raising
`-max_target_seqs` or using `-subject_besthit`); add a `timeout=` to the subprocess. For
krona-input, emit rows in a sorted, stable order.

**Rationale / impact.** Auto-selected reference determines the mapping reference and thus
the consensus — a reproducibility hazard for reported genomes; a hung blastn can stall a
run with no bound.

**Suggested test.** Same contigs + DB run twice → identical selected reference; a stable
byte-for-byte krona-input file across runs (fixed `PYTHONHASHSEED`).

---

## Appendix — hardcoded (non-configurable) scientific thresholds

Reviewers should be aware these defaults are baked in (some overridable via config but
undocumented). Consider promoting the load-bearing ones to documented config keys.

| Threshold | Default | Location | Configurable? |
|---|---|---|---|
| RPM bleed floor | `1.0` | all `metagenomics_*` rules (e.g. `metagenomics_kraken2_reads_illumina.smk:150`) | No (hardcoded); `bleed_fraction` 0.005 is |
| iSNV upper band | `0.5` | `consensus_illumina.smk:46` | No (lower `af_isnv_threshold` 0.05 is) |
| Nanopore variant gates | qual 20 / mapq 20 / depth 5 / chunk 50000 | `consensus_nanopore.smk:26-30` | Config-overridable, undocumented |
| Reference reads/contigs count | 100 / 1 | `select_reference_genomes.py:30-31` | CLI |
| BLAST qcov / pident | 80 / 80 | `select_reference_genomes.py:36-37` | CLI |
| DIAMOND evalue / max_target_seqs / sensitivity | 0.001 / 1 / sensitive | `metagenomics_diamond_reads_illumina.smk:16-18` | config-only |
| Kraken2 minimum_hit_group | 4 | `metagenomics_kraken2_reads_illumina.smk:13` | config |
| Enrichment gate | z≥3.0 / log2≥1.0 / pseudocount 1.0 | `add_negative_control_enrichment.py:169-171` | config |
| Coverage completeness breakpoints | 10×/100×/1000× | `calculate_assembly_stats.py:60-62` | No (reporting only) |

---

## How to work through this document

Each item is independent. When you accept one, the corresponding fix belongs in a small
dedicated branch (`fix/science-<item>`), with the "Suggested test" turned into a real test
that pins the new behavior. Items 2 and 3 should be tackled together (same root cause).
Item 5 may be a documentation-only outcome if the per-platform difference is intentional.
