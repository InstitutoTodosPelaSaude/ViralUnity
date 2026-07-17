# Design note: automating the metagenomics → consensus linkage

**Status:** internal design / roadmap — not part of the published documentation
(excluded from the Sphinx build). No code changes are described here; this records
options for future versions.

**Origin:** the REVISA consensus stress test (July 2026) ran
`viralunity consensus illumina` across 305 metagenomics-derived targets / 65
reference groups and compared reference-based consensus against the de novo
assembly. The two pipeline-breaking bugs it found are already fixed (see
`CHANGELOG.md`, 1.4.0 `### Fixed`). This note captures the *larger* findings that
are features rather than fixes — most of them about tightening the automatic
"assemble a consensus for every virus metagenomics found" path.

---

## Where the linkage stands today

ViralUnity can already assemble a consensus after metagenomics via
`viralunity meta --run-reference-assembly`
(`viralunity/scripts/rules/metagenomics_reference_assembly.smk`):

- A Snakemake `checkpoint select_references_meta` reads the fully bleed/neg/ICTV
  filtered taxa-summary TSVs and calls
  `viralunity/scripts/python/select_reference_genomes.py`, which picks a reference
  accession per `(sample, taxon)` and writes `reference_targets.tsv`
  (`sample, ref_key, reference_genome`).
- `rule extract_reference_fasta` pulls each accession out of `--viral-genomes` and
  the workflow reuses the ordinary consensus rules (`alignment_illumina.smk`,
  `consensus_illumina.smk`) keyed by a `{ref_key}/` wildcard, producing one
  `assembly/{ref_key}/consensus/final_consensus/{sample}.consensus.fasta` per
  `(sample, ref_key)` row.
- Two selection strategies exist: `taxid` (default; exact/species-fallback taxid
  lookup) and `similarity` (blastn of de novo contigs against `--viral-genomes`).

Known limits of the current path, relevant below:

- **Verification is family-level only.** Both strategies accept an accession if its
  taxid traces up to a requested *family*; species-rank taxids are used only as a
  lookup bridge (`select_reference_genomes.py`, taxid strategy). There is no check
  that the chosen accession's species matches the classifier's species call.
- **No segmented variant.** The meta path never includes `consensus_*_segmented.smk`;
  it treats every accession as a single-FASTA reference. Segmented viruses
  (influenza, rotavirus, bunyaviruses) cannot be assembled as multi-segment genomes
  after metagenomics.
- **Grouping is per-`(sample, ref_key)`.** Each sample re-extracts its own copy of a
  shared accession and is assembled independently; there is no cross-sample
  collation, and reads map to fastp-trimmed reads (not the dehosted reads produced
  for classification).

---

## Proposed enhancements

### 1. Species-level reference verification

Cross-check the chosen accession's *species* against the classifier's species call,
not just the family, and exclude mismatches. The stress test did exactly this
(BLAST hit's family/species vs the metagenomics call) and excluded targets whose
best reference disagreed at family level, rather than silently assembling against a
wrong-family genome. This tightens the `similarity` and `taxid` strategies in
`select_reference_genomes.py` (reuse `get_taxid_at_rank(..., "species")`, already
present). For `*sp.*` / genus-only calls, fall back to family-level as today.

### 2. Segmented reference selection and assembly after metagenomics

There is currently no segmented meta path. A ready-made approach, validated in the
stress test: once the best hit is known, take **all accessions in the
accession→taxid map that share the hit's strain-level taxid** as the segment set
(this guarantees the segments come from one genome — e.g. the 8 PR8 accessions
under one influenza-A taxid), label them `seg1..segN`, and drive the existing
`consensus_illumina_segmented.smk` / `consensus_nanopore_segmented.smk` machinery.
The meta path already repurposes the segment wildcard slot for `ref_key`, so wiring
a genuine segmented variant is a natural extension rather than a rewrite.

- Validate the segment count against a per-family expectation (influenza 6–8,
  rotavirus 10–12, bunyaviruses 3); if a taxid resolves to an anomalous number of
  accessions, fall back to the single best accession and record why.
- **Caveat to document for users:** the strain a sample is grouped under is decided
  by whichever *single* segment happened to be its largest/most-confident contig, so
  two samples of the same real strain can land in different strain-groups if their
  largest contigs came from different segments. Per-segment reference selection (or
  a consensus-of-segment-strains step) would be more principled for a definitive
  analysis.

### 3. Reference-level sample grouping

Group samples by the *selected reference* rather than assembling each
`(sample, accession)` in isolation. In the stress test this produced
subtype/serotype-consistent runs for free — dengue split into per-serotype
references, influenza A collapsed into H3N2 vs H1N1pdm09 groups — because BLAST of
the largest contig already discriminates strains. Benefits: one shared reference
extraction per group, and the option of a cross-sample alignment / combined report
per reference (the current path collates nothing across samples).

### 4. Consensus-vs-de-novo QC panel

Add an automatic QC comparison (as a report panel or a summary TSV) that flags when
the nearest RefSeq reference is too divergent for reference-based consensus to be
trusted. Stress-test evidence over 159 targets with ≥70 % contig coverage:

- Where both the consensus and the de novo contig confidently call a base, they
  agree at ~100 % median identity — good cross-validation.
- But de novo recovers more genome: the contig called bases the consensus left as
  `N` in 124/159 targets, and 62/159 consensuses recovered **<50 %** of the
  reference — almost all divergent taxa (enteroviruses, pegivirus) whose nearest
  RefSeq genome is only 58–80 % identical, so read mapping (and thus the consensus)
  collapses while de novo assembly still reconstructs the genome.

Surfacing "reference divergence too high — consensus under-recovers vs de novo" per
target would let users catch this automatically.

### 5. Hybrid / contig-as-reference (needs PI sign-off)

For divergent taxa, a reference-based consensus against the nearest RefSeq genome
systematically under-recovers relative to the de novo assembly metagenomics already
produced. Options: scaffold the consensus onto the de novo contig, or use the contig
itself as the reference when its identity to RefSeq is low. **This is a
scientific-method decision, not a mechanical fix — it changes what "the consensus"
means — and should be decided by the PI before any implementation.**

### 6. Minor / documentation-only observations

- **Dehosted reads not carried into reference assembly.** The reference-assembly
  `map_reads` consumes fastp-trimmed reads (`perform_qc` output), not the
  host-depleted reads from `metagenomics_dehost_illumina.smk`. Carrying dehosting
  through to mapping is a possible refinement.
- **"Skip fastp / accept pre-trimmed reads" is *not* needed for the automated path.**
  The stress test observed redundant re-trimming only because it launched a separate
  `viralunity consensus` per reference group; within a single `viralunity meta`
  invocation, `perform_qc` (`qc_illumina.smk`) is keyed by `{sample}` only and its
  trimmed reads are already reused across every `ref_key` target. Noted for
  completeness; low priority.
- **RefSeq DB resolves subtypes/serotypes better than assumed.** A prior assumption
  held the local `virus_genomes` DB to ~one strain per species, so influenza
  subtypes could not be resolved. That is wrong for well-sampled respiratory/arbo
  viruses: BLAST of the largest contig split influenza A into H3N2 and H1N1pdm09,
  influenza B to B/Lee/1940, and dengue into per-serotype references automatically.
  This is what makes reference-level grouping (item 3) worthwhile.

---

## Where these would land

Items 1–3 belong in `select_reference_genomes.py` and
`metagenomics_reference_assembly.smk` (extending the existing `similarity`/`taxid`
strategies and adding a segmented include). Item 4 fits the consensus report
generator (`generate_consensus_report.py`) or a meta-side summary. Item 5 is a
larger, PI-gated change to the assembly strategy itself.
