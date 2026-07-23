# Commands Reference

## `viralunity setup`

Pre-builds the per-rule conda environments declared in the Snakemake workflows into a shared cache directory. Run once after installing ViralUnity (or after upgrading); subsequent `viralunity consensus` / `viralunity meta` runs reuse the cached envs and skip env creation.

```bash
viralunity setup --pipelines all
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--conda-prefix` | `$VIRALUNITY_CONDA_PREFIX` or `~/.cache/viralunity/conda-envs` | Directory where per-rule envs are cached. Reused by every later pipeline run that points at the same prefix. |
| `--pipelines` | `all` | Repeatable. Pick from `consensus-illumina`, `consensus-nanopore`, `meta-illumina`, `meta-nanopore`, or `all`. Segmented variants share envs with their non-segmented counterparts. |
| `--threads` | `4` | Cores given to Snakemake while materializing envs. |
| `--dry-run` | off | Print the envs that would be created and exit without invoking conda. |

### Examples

**Build everything once after install:**

```bash
viralunity setup --pipelines all
```

**Build only what you need:**

```bash
viralunity setup --pipelines consensus-illumina --pipelines meta-illumina
```

**Inspect first:**

```bash
viralunity setup --pipelines consensus-illumina --dry-run
```

**Use a shared cache on a cluster:**

```bash
export VIRALUNITY_CONDA_PREFIX=/shared/viralunity-envs
viralunity setup --pipelines all
```

Every subsequent `viralunity consensus` / `viralunity meta` invocation will pick up `$VIRALUNITY_CONDA_PREFIX` automatically; no per-run flag needed.

---

## `viralunity create-samplesheet`

Before running any pipeline you need a CSV sample sheet. The `create-samplesheet` command generates it automatically from a sequencing run directory.

```bash
viralunity create-samplesheet --input /path/to/run/ --output samples.csv
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | *(required)* | Run directory: contains FASTQ files directly (level 0) or one subdirectory per sample (level 1). |
| `--output` | *(required)* | Output CSV file path. |
| `--level` | `1` | `0` = files in `--input`, `1` = files in subdirectories (one per sample). |
| `--separator` | `-` | Character used to split the file/directory name and extract the sample ID (`-`, `_`, or `.`). |
| `--pattern` | `R1` | Pattern identifying the first read file at level 0 (`R1` for Illumina, `barcode` for Nanopore). |

### Output format

The generated CSV has no header row:

```
# Illumina (3 columns)
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz

# Nanopore (2 columns)
barcode01,/path/to/barcode01.fastq
```

### Examples

**Illumina** — FASTQs organised in one subdirectory per sample (level 1, default):

```bash
viralunity create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv
```

**Illumina** — all FASTQs directly in the run directory (level 0):

```bash
viralunity create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

**Nanopore** — one FASTQ per barcode subdirectory:

```bash
viralunity create-samplesheet \
    --input /path/to/nanopore_run/ \
    --output /path/to/samples_nano.csv \
    --separator _ \
    --pattern barcode
```

---

## `viralunity consensus`

The consensus pipeline takes raw reads to processed consensus genome sequences with a single command. Select the data type as a subcommand (`illumina` or `nanopore`).

### Options — shared (both data types)

| Option | Default | Description |
|--------|---------|-------------|
| `--sample-sheet` | *(required)* | CSV file with sample IDs and file paths. |
| `--config-file` | *(required)* | Path for the YAML config file to be created. |
| `--output` | *(required)* | Base output directory. |
| `--reference` | — | Reference genome FASTA (mutually exclusive with `--segmented-reference`). A multi-record FASTA is treated as a segmented virus and split into one reference per record, with segment names taken from the headers (first `\|`/whitespace token, sanitized). |
| `--single-reference` | off | Treat a multi-record `--reference` as a single reference (align all records together in one pass) instead of auto-splitting it into segments. |
| `--segmented-reference` | — | Per-segment reference: `SEGMENT=PATH` (repeatable). |
| `--primer-scheme` | — | Primer scheme BED file (amplicon sequencing only). |
| `--minimum-coverage` | `20` | Minimum depth for consensus base inclusion. |
| `--minimum-read-length` | `50` | Minimum read length threshold. |
| `--af-threshold` | `0.51` | Min allele frequency to call variant into consensus. |
| `--run-name` | `undefined` | Name for the sequencing run. |
| `--threads` | `1` | Threads per individual task. |
| `--threads-total` | `1` | Total threads for the workflow. |
| `--create-config-only` | off | Only generate the config file; do not run the workflow. |
| `--conda-prefix` | `~/.cache/viralunity/conda-envs` | Cache directory for per-rule conda envs. Picked up from `$VIRALUNITY_CONDA_PREFIX` if set. Pre-warm with `viralunity setup`. |

### Options — Illumina only

| Option | Default | Description |
|--------|---------|-------------|
| `--adapters` | — | Adapter sequences FASTA (fastp QC). |
| `--trim-head` | `0` | Bases to trim from 5′ end. |
| `--trim-tail` | `0` | Bases to trim from 3′ end. |
| `--cut-front-mean-quality` | `10` | fastp cut_front quality threshold. |
| `--cut-tail-mean-quality` | `10` | fastp cut_tail quality threshold. |
| `--cut-right-window-size` | `4` | fastp cut_right window size. |
| `--cut-right-mean-quality` | `15` | fastp cut_right quality threshold. |
| `--af-isnv-threshold` | `0.0` | Min allele frequency for iSNV analysis. |
| `--run-isnv` | off | Run intra-host SNV analysis (LoFreq). |

### Options — Nanopore only

| Option | Default | Description |
|--------|---------|-------------|
| `--chunk-size` | `10000` | Chunk size for clair3 processing. |
| `--clair3-model` | `r1041_e82_400bps_sup_v500` | Clair3 model for variant calling. |
| `--variant-quality` | `20` | Minimum variant quality (clair3). |
| `--variant-depth` | `10` | Minimum alt allele depth (clair3). |
| `--minimum-map-quality` | `30` | Minimum mapping quality (clair3). |

### Examples

**Illumina — single reference:**

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/references/viral_genome.fasta \
    --primer-scheme /path/to/primers.bed \
    --adapters /path/to/adapters.fa \
    --threads 2 \
    --threads-total 4
```

**Illumina — segmented genome, single multi-FASTA** (simplest; one file with all
segments, segment names taken from the FASTA headers):

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_segmented.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/references/all_segments.fasta \
    --threads 2 \
    --threads-total 4
```

**Illumina — segmented genome, one file per segment** (explicit segment names via
repeatable `--segmented-reference`):

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_segmented.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --segmented-reference L=/path/to/L_segment.fasta \
    --segmented-reference S=/path/to/S_segment.fasta \
    --threads 2 \
    --threads-total 4
```

**Illumina — single reference with multiple contigs (fragmented genome)**: when a
genome is only available as several contigs/scaffolds but should be treated as **one**
reference (all contigs aligned together in a single pass, e.g. to report one
whole-assembly horizontal coverage), pass the multi-record FASTA with
`--single-reference` so it is *not* auto-split into segments:

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/fragmented_assembly.fasta \
    --single-reference \
    --threads 2 \
    --threads-total 4
```

Behavior notes for a multi-contig single reference:

- **Coverage/depth statistics** (`assembly_stats_summary.csv`) are aggregated across
  all contigs: `horizontal_coverage` is the fraction of *all* reference positions at
  or above the threshold, and `average_depth` is the whole-assembly mean.
- **Consensus** keeps one record per contig (contigs are never fused).
- In the **HTML report**, the coverage track lays the contigs end-to-end along one
  genome axis (in FASTA order).
- The cross-sample alignment is written **per contig** under
  `consensus/final_consensus/per_contig_alignments/<contig>.fasta` (one MSA each).
  A concatenation of all of them is also kept as `samples_alignment.fasta` (one
  padded block per contig, not a single rectangular alignment) for backward
  compatibility.
- **Known limitation:** a `--gene-annotation` track on a multi-contig single
  reference only renders features for the *first* contig. For per-contig
  annotation, run in segmented mode instead (omit `--single-reference`; a single
  combined `--gene-annotation` is then split per segment automatically).

**Nanopore:**

```bash
viralunity consensus nanopore \
    --sample-sheet /path/to/samplesheet_nano.csv \
    --config-file /path/to/example_nano.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/reference.fasta \
    --threads 4 \
    --threads-total 4
```

**Config only** (generates config without running the workflow):

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --output /path/to/example_output \
    --reference /path/to/reference.fasta \
    --create-config-only
```

---

## `viralunity meta`

The metagenomics pipeline takes raw reads to taxonomic classifications and visualizations. You can use **Kraken2 only**, **Diamond only**, or **both**. Each tool can be run on **reads** and, when assembly is enabled, on **contigs**.

### Pipeline overview (Illumina)

1. **Quality control** — [fastp](https://github.com/OpenGene/fastp): trimming, quality filtering, adapter handling, JSON/HTML reports.
2. **Optional dehosting** — minimap2 (if `--host-reference` is set) or [Deacon](https://github.com/bede/deacon) (if `--deacon-index` is set).
3. **Read classification** — [Kraken2](https://github.com/DerrickWood/kraken2) and/or [DIAMOND](https://github.com/bbuchfink/diamond) blastx on reads (optional each).
4. **Optional de novo assembly** — [MEGAHIT](https://github.com/voutcn/megahit) on host-filtered pairs.
5. **Contig classification** — Kraken2 and/or Diamond on contigs (when assembly is run).
6. **Summaries and filters** — Per-sample taxa tables, RPM normalisation (and RPKM when `--viral-genomes` is supplied), optional max-RPM bleed filter, and negative-control enrichment filter (fold-enrichment, log10-ratio, z-score; replaces the old Poisson filter). [Krona](https://github.com/marbl/Krona) plots (both a raw `*.krona.html` and a post-filter `*.filtered.krona.html` per sample/classifier/mode) and [MultiQC](https://multiqc.info) reports.
7. **Optional dynamic reference assembly** — Automatic reference sequence selection from taxonomic hits and subsequent consensus assembly.

### Pipeline overview (Nanopore)

1. **Optional dehosting** — minimap2 or Deacon.
2. **Read classification** — Kraken2 and/or Diamond on reads.
3. **Optional de novo assembly** — MEGAHIT on host-filtered reads.
4. **Optional polishing** — [Racon](https://github.com/lbcb-sci/racon) and/or [Medaka](https://github.com/nanoporetech/medaka) on MEGAHIT contigs.
5. **Contig classification** — Kraken2 and/or Diamond on contigs.
6. **Summaries and filters** — Same as Illumina.
7. **Optional dynamic reference assembly** — Same as Illumina.

### Options — shared (both data types)

| Option | Default | Description |
|--------|---------|-------------|
| `--sample-sheet` | *(required)* | CSV with sample IDs and file paths. |
| `--config-file` | *(required)* | Path to YAML config file to be generated. |
| `--output` | *(required)* | Base output directory. |
| `--run-name` | `undefined` | Run name; output goes under `{output}/{run_name}/`. |
| `--kraken2-database` | `NA` | Path to Kraken2 database. |
| `--krona-database` | `NA` | Path to Krona taxonomy (required for any classification). |
| `--taxdump` | `NA` | NCBI taxdump dir (`nodes.dmp`, `names.dmp`). |
| `--remove-human-reads` | off | Remove human reads from summaries/Krona. |
| `--remove-unclassified-reads` | off | Remove unclassified reads from summaries/Krona. |
| `--host-reference` | `NA` | Host genome FASTA for minimap2 dehosting. |
| `--deacon-index` | `NA` | Deacon minimizer index for host depletion. |
| `--run-denovo-assembly` | off | Run MEGAHIT de novo assembly. |
| `--run-kraken2-reads/--no-kraken2-reads` | on | Kraken2 on reads. |
| `--run-kraken2-contigs/--no-kraken2-contigs` | on | Kraken2 on contigs. |
| `--run-diamond-reads/--no-diamond-reads` | off | Diamond on reads. |
| `--run-diamond-contigs/--no-diamond-contigs` | off | Diamond on contigs. |
| `--diamond-database` | `NA` | Protein FASTA for Diamond. |
| `--taxids` | `NA` | NCBI taxid mapping for Diamond taxonomy. |
| `--diamond-sensitivity` | `sensitive` | `sensitive` / `mid-sensitive` / `more-sensitive` / `ultra-sensitive`. |
| `--evalue` | `1e-10` | Diamond E-value threshold. |
| `--bleed-fraction` | `0.005` | Max-RPM bleed filter fraction. |
| `--negative-controls` | (empty) | Comma-separated sample IDs used as negative controls. |
| `--enrichment-pseudocount` | `1.0` | Pseudocount for fold-enrichment and log10-ratio vs negative controls. |
| `--z-score-threshold` | `3.0` | Z-score threshold for `neg_pass` (used when ≥2 negative controls are present). |
| `--log10-ratio-threshold` | `1.0` | Log10-ratio threshold for `neg_pass`; the default `1.0` marks a 10-fold enrichment over controls. Used when exactly 1 negative control is present, or when z-score is undefined due to zero control variance. |
| `--minimum-hit-group` | `4` | Kraken2 minimum-hit-group parameter. |
| `--run-ictv-host-filter`/`--no-ictv-host-filter` | off | Keep only vertebrate-infecting viruses (drop phages, plant/fungal/algal/invertebrate-only viruses). See *Taxonomic false-positive filters* below. |
| `--ictv-vertebrate-taxids-file` | `NA` | Vertebrate-virus taxid allowlist (built by the `build_ictv_vertebrate_taxids.py` script — see *Taxonomic false-positive filters* below). Required with `--run-ictv-host-filter`. |
| `--run-nr-validation`/`--no-nr-validation` | off | Re-search de novo viral contigs against NCBI nr and keep only NR-confirmed viral contigs. Contig tracks only; requires `--run-denovo-assembly` and `--run-diamond-contigs`. |
| `--nr-diamond-database` | `NA` | nr database for NR validation: a BLAST+ nr db (searched via `diamond prepdb`) or a native `.dmnd`. |
| `--nr-evalue` | `1e-10` | E-value threshold for the NR DIAMOND search. |
| `--nr-max-target-seqs` | `10` | Top hits kept per contig for the NR LCA consensus. |
| `--nr-sensitivity` | `fast` | DIAMOND sensitivity for the NR search. |
| `--nr-consensus-threshold` | `0.5` | Fraction of a contig's hits that must agree at a rank for the LCA consensus. |
| `--run-reference-assembly`/`--no-run-reference-assembly` | off | Enable reference assembly from filtered taxonomic hits. |
| `--method` | *(none)* | Method used for reference assembly (`kraken2`, `diamond`, `both`). Required when `--run-reference-assembly` is set. |
| `--source` | *(none)* | Source of taxonomy data for reference assembly (`reads`, `contigs`, `both`). Required when `--run-reference-assembly` is set. |
| `--reads-count` | `100` | Minimum reads assigned to a viral family to trigger reference assembly. |
| `--contigs-count` | `1` | Minimum contigs assigned to a viral family to trigger reference assembly. |
| `--families` | `Coronaviridae,...` | Comma-separated list of targeted viral families for reference assembly. |
| `--reference-selection-strategy` | `taxid` | Strategy to associate a reference genome (`taxid` or `similarity`). |
| `--blast-qcov` | `80` | Minimum query coverage for similarity strategy. |
| `--blast-pident` | `80` | Minimum percent identity for similarity strategy. |
| `--viral-genomes` | `NA` | FASTA of viral genomes for reference assembly database. |
| `--viral-taxids` | `NA` | TSV of genome to TaxID mapping for reference assembly. |
| `--threads` | `1` | Threads per task. |
| `--threads-total` | `1` | Total threads for the workflow. |
| `--create-config-only` | off | Only generate the config; do not run the workflow. |
| `--conda-prefix` | `~/.cache/viralunity/conda-envs` | Cache directory for per-rule conda envs. Picked up from `$VIRALUNITY_CONDA_PREFIX` if set. Pre-warm with `viralunity setup`. |

### Options — Illumina only (fastp QC)

| Option | Default | Description |
|--------|---------|-------------|
| `--adapters` | `NA` | Adapter FASTA; `NA` = auto-detect. |
| `--minimum-read-length` | `50` | Minimum read length after trimming. |
| `--trim-head` | `0` | Bases to trim from 5′. |
| `--trim-tail` | `0` | Bases to trim from 3′. |
| `--cut-front-mean-quality` | `20` | fastp cut_front quality threshold. |
| `--cut-tail-mean-quality` | `20` | fastp cut_tail quality threshold. |
| `--cut-right-window-size` | `4` | fastp cut_right window size. |
| `--cut-right-mean-quality` | `20` | fastp cut_right quality threshold. |

### Options — Nanopore only (polishing)

| Option | Default | Description |
|--------|---------|-------------|
| `--run-polish-racon/--no-polish-racon` | off | Racon polishing on MEGAHIT assembly. |
| `--run-polish-medaka/--no-polish-medaka` | off | Medaka polishing on assembly. |
| `--medaka-model` | — | Medaka model (e.g. `r941_min_high_g360`). |

### Examples

**Illumina — Kraken2 on reads (default):**

```bash
viralunity meta illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_meta.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --threads 2 \
    --threads-total 4
```

**Illumina — with adapters and dehosting:**

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --adapters /path/to/adapters.fa \
    --host-reference /path/to/host_genome.fa \
    --threads 4 --threads-total 8
```

**Illumina — full pipeline (assembly + Kraken2 + Diamond, reads and contigs):**

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-diamond-reads --run-diamond-contigs \
    --diamond-database /path/to/proteins.faa \
    --taxids /path/to/protein2taxid.tsv \
    --run-denovo-assembly \
    --host-reference /path/to/host.fa \
    --threads 4 --threads-total 8
```

**Illumina — with reference assembly (Kraken2 hits drive consensus assembly):**

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-reference-assembly \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids databases/virus_genomes/genome2taxid.tsv \
    --threads 4 --threads-total 8
```

**Nanopore — full pipeline (dehosting, assembly, Medaka polishing, Kraken2 + Diamond):**

```bash
viralunity meta nanopore \
    --sample-sheet /path/to/samplesheet_nano.csv \
    --config-file /path/to/config_nano.yml \
    --run-name example_run \
    --output /path/to/output_nano \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --diamond-database /path/to/proteins.faa \
    --taxids /path/to/protein2taxid.tsv \
    --host-reference /path/to/host.fa \
    --run-denovo-assembly \
    --run-kraken2-reads --run-kraken2-contigs \
    --run-diamond-reads --run-diamond-contigs \
    --run-polish-medaka --medaka-model r941_min_high_g360 \
    --threads 4 --threads-total 4
```

### Krona outputs (raw vs. filtered)

For every classifier/mode that runs, the pipeline emits two Krona HTMLs per sample:

- `samples/<sample>/<classifier>_<mode>.krona.html` — built directly from the per-sample taxonomic classification before any cross-sample filtering. This is the raw view of what the classifier reported.
- `samples/<sample>/<classifier>_<mode>.filtered.krona.html` — built from the same krona input, but pruned to taxa that survive the filters recorded in the fully-filtered taxa-summary table (the longest-named file in the cumulative chain, e.g. `_taxa_summary_RPKM.nr.ctgstats.bleed.neg.ictv.tsv`) for that `(sample, tool, mode)`.

Filtering is *lineage-aware*: a contig/read is kept when any ancestor of its leaf taxid at the `family`, `genus`, or `species` rank passes `bleed_pass` (and, when negative controls are configured, also `neg_pass`). This is the inverse of how `summarize_krona_taxa.py` aggregates rows up the lineage, so the filtered Krona shows exactly the contigs/reads that contributed to a passing rank-row. Strain and sub-species hits whose species/genus/family passes are preserved; rows with `taxid==0` are dropped.

**`neg_pass` semantics** — set by the negative-control enrichment filter:
- `neg_pass = NA` — no negative controls configured (zero-control mode); treated as *keep* by the Krona filter.
- `neg_pass = True` with `neg_decision = "absent_from_controls"` — controls *are* configured but the taxon is absent from every one of them (`control_max == 0`). The enrichment ratios (`fold_enrichment`, `log10_ratio`, `agg_fold_enrichment`, `agg_log10_ratio`) would only track the sample's own abundance against the pseudocount, so they are blanked to NA and every enrichment gate auto-passes — absence from the controls is itself the evidence the taxon is not control-borne. The z-score gates (`neg_pass_5`/`neg_pass_10`) stay NA (a z-score is undefined without control signal).
- `neg_pass = True/False` with `neg_decision = "log10_ratio"` — one negative control; gate is `log10_ratio ≥ --log10-ratio-threshold`.
- `neg_pass = True/False` with `neg_decision = "z_score"` — two or more controls; gate is `z_score ≥ --z-score-threshold`.
- `neg_pass = True/False` with `neg_decision = "log10_ratio_fallback"` — two or more controls but all identical (SD = 0); falls back to the log10-ratio gate.

Additional diagnostic columns added by the negative-control step (the `.neg` suffix): `neg_metric` (`rpkm` or `rpm`), `control_mean`, `control_sd`, `control_median`, `control_max`, `fold_enrichment`, `log10_ratio`, `z_score`, `enrichment_pseudocount`, `z_score_threshold_used`, `log10_ratio_threshold_used`, `n_negative_controls`.

Convenience pass/fail flags are also emitted (`>=` inclusive; NA where the underlying statistic is NA): `fold_enrichment_10x_pass` (`fold_enrichment >= 10`), `fold_enrichment_100x_pass` (`fold_enrichment >= 100`), `neg_pass_5` (`z_score >= 5`) and `neg_pass_10` (`z_score >= 10`).

A complementary **aggregate (pooled) control** view is also emitted (diagnostic only — it does not change `neg_pass`): `pooled_control_metric` treats all controls as one pooled library (raw reads pooled, i.e. each control weighted by its library size; absent controls contribute 0), plus `agg_fold_enrichment`, `agg_log10_ratio`, and `agg_fold_enrichment_10x_pass`/`agg_fold_enrichment_100x_pass`. This catches widespread, high-variance contamination that inflates per-control SD and drags the z-score down.

`viralunity/scripts/python/filter_krona_by_pass_taxids.py` also exposes a CLI for standalone use:

```bash
python viralunity/scripts/python/filter_krona_by_pass_taxids.py \
    --summary path/to/diamond_contigs_taxa_summary_RPKM.nr.ctgstats.bleed.neg.ictv.tsv \
    --krona-input path/to/{sample}.diamond.supported.krona_input.tsv \
    --out path/to/{sample}.diamond.supported.filtered.krona_input.tsv \
    --sample {sample} --tool diamond --mode contigs \
    --taxdump path/to/taxdump
```

### Taxonomic false-positive filters (optional)

Open classifiers (kraken2 + viral index, DIAMOND + viral RefSeq) can report taxa that are
artifacts of imperfect taxonomic identification — bacteriophages, giant algal viruses, endogenous
retroviruses, and bacterial/eukaryal contigs mislabelled as viral. Three optional, tunable filters
remove these. They are **off by default**, and each writes a `*.dropped.tsv` audit sidecar next to
its output listing the removed rows and why.

**Order.** The taxa summary flows through **one cumulative chain**; each enabled step appends
exactly one suffix, in this order:

```
_RPM/_RPKM → .nr → .ctgstats → .bleed → .neg → .ictv
```

so a fully-loaded contig track ends at `..._taxa_summary_RPKM.nr.ctgstats.bleed.neg.ictv.tsv` and
the *fully-filtered* table is always the file with the longest name. NR validation and the
`.ctgstats` largest-contig statistics are contig-tracks only (and `.ctgstats` additionally requires
`--viral-genomes`); the ICTV host filter applies to all tracks. The order is **result-neutral** with
respect to the surviving taxa: `.ctgstats`, `.bleed` and `.neg` only add per-taxon columns (they
never remove rows and don't depend on which other taxa are present), while the row-removing steps
(`.nr`, `.ictv`) use criteria independent of those columns. With every filter off the chain is just
`..._RPM.bleed.tsv` (or `_RPKM.bleed.tsv` with `--viral-genomes`), byte-identical in content to
before. Row-removing steps write a `*.dropped.tsv` audit sidecar.

| Filter | Flag | Scope | What it drops |
|--------|------|-------|---------------|
| ICTV host | `--run-ictv-host-filter` + `--ictv-vertebrate-taxids-file` | all four tracks | Taxa not in a vertebrate-infecting virus lineage (phages, plant/fungal/algal/invertebrate-only, giant viruses). Lineage-aware against an ICTV-derived taxid allowlist (built by `build_ictv_vertebrate_taxids.py`, see below). |
| NR validation | `--run-nr-validation` + `--nr-diamond-database` | contig tracks only | Contigs the NR (full `nr`) LCA consensus confidently calls non-viral. Runs one aggregated `diamond blastx` over all samples' de novo viral contigs; requires `--run-denovo-assembly` and `--run-diamond-contigs`. |

NR validation also writes `<track>_nr_flags.tsv` surfacing (not filtering) species-level
RefSeq-vs-NR disagreements — `misid_novel` (an NR virus species absent from that sample's RefSeq
calls, highest surveillance interest) and `misid_known`.

**Building the ICTV allowlist.** The vertebrate-virus taxid allowlist is derived from the ICTV
Virus Metadata Resource (VMR; keep taxa whose "Host source" contains `vertebrates`) and is a
plain, user-editable list of NCBI taxids. Download the current VMR spreadsheet from
<https://ictv.global/vmr>, then:

```bash
python viralunity/scripts/python/build_ictv_vertebrate_taxids.py \
    --vmr VMR_MSLxx.xlsx --names <taxdump>/names.dmp \
    --output databases/ictv/vertebrate_virus_taxids.txt
```

The script maps ICTV family/genus names to taxids via the taxdump `names.dmp`; regenerate it
against a new VMR/taxdump release, or edit the taxid list by hand to tune coverage.

**Aggregated contig search (config-only).** `combine_contig_search: true` in the YAML runs the
per-sample contig DIAMOND search once over all samples' contigs combined (sample-prefixed headers),
splitting the output per sample — often faster (one DB load). Off by default; toggle via
`--create-config-only`, edit the YAML, and rerun. Benchmark against the per-sample path using the
`logs/diamond_contigs/*.benchmark.txt` wall-times before adopting it.

**Not yet in scope (future).** Physically slimmed diamond/kraken2 databases, an empirical
REVISA-trained FP classifier, and a targeted mapping/genotyping track for a curated
virus-of-interest panel are documented directions, not implemented.

---

## `viralunity get-databases`

The `get-databases` command downloads and sets up the reference databases required by the meta pipeline. Each subcommand takes an optional `--path` argument (default: `databases` in the current directory), inside which a named subdirectory is created.

```bash
viralunity get-databases --help
viralunity get-databases kraken2 --help
```

### `kraken2` — Kraken2 pre-built index

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/kraken2/`. |
| `--url` | k2_viral_20240112 | URL of the pre-built index archive. Check [available indexes](https://benlangmead.github.io/aws-indexes/k2). |

```bash
viralunity get-databases kraken2 --path /data/dbs
# database at /data/dbs/kraken2/
# use: --kraken2-database /data/dbs/kraken2
```

### `krona` — Krona taxonomy

Requires the viralunity conda environment to be active (`CONDA_PREFIX` must be set).

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/krona/taxonomy/`. |

```bash
conda activate viralunity
viralunity get-databases krona --path /data/dbs
# taxonomy at /data/dbs/krona/taxonomy/
# use: --krona-database /data/dbs/krona/taxonomy
```

### `taxdump` — NCBI taxdump

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/taxdump/`. |
| `--url` | NCBI taxdump | URL of the taxdump archive. |

```bash
viralunity get-databases taxdump --path /data/dbs
# nodes.dmp, names.dmp, etc. at /data/dbs/taxdump/
# use: --taxdump /data/dbs/taxdump
```

### `diamond` — Diamond protein database

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/diamond/`. |
| `--taxon` | Viruses | NCBI taxon name to download (e.g. `Viruses`, `coronaviridae`). |
| `--refseq/--no-refseq` | `on` | Limit to RefSeq genomes only. |
| `--threads` | `1` | Threads for `diamond makedb`. |
| `--skip-makedb` | off | Download and reformat files only; skip `diamond makedb`. |

Requires the NCBI Datasets CLI (`conda install -c conda-forge ncbi-datasets-cli`).

```bash
viralunity get-databases diamond --path /data/dbs --threads 4
# /data/dbs/diamond/viral.dmnd
# /data/dbs/diamond/protein2taxid.tsv
# use: --diamond-database /data/dbs/diamond/viral.dmnd
#      --taxids /data/dbs/diamond/protein2taxid.tsv
```

NCBI Datasets occasionally bundles a small number of nucleotide CDS records
(with a genome accession such as `NC_xxxxxx.1` as the header) inside the
`protein.faa` it returns for `--include protein`. The `diamond` subcommand
automatically detects and drops these records before running
`diamond makedb` and prints a warning summarizing how many were skipped.

### `clean-protein-fasta` — Strip nucleotide records from a protein FASTA

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | *(required)* | Input protein FASTA file (e.g. `databases/diamond/viral.protein.faa`). |
| `--output` | `<input>.cleaned.faa` | Output FASTA path. If equal to `--input`, the file is replaced atomically. |
| `--no-backup` | off | When replacing the input in place, do not keep a `.with_dna.bak` backup. |
| `--min-dna-length` | `20` | Minimum sequence length to classify a record as DNA. |

Use this command to recover an already-downloaded `viral.protein.faa` that
fails `diamond makedb` with:

```
Error: The sequences are expected to be proteins but only contain DNA letters.
```

```bash
viralunity get-databases clean-protein-fasta \
    --input databases/diamond/viral.protein.faa \
    --output databases/diamond/viral.protein.faa

diamond makedb \
    --in databases/diamond/viral.protein.faa \
    --db databases/diamond/viral.dmnd \
    --threads 4
```

### `nr` — NCBI nr database for NR validation

Downloads and configures the full NCBI **nr** protein database used by
`meta --run-nr-validation` (`--nr-diamond-database`). By default it fetches
NCBI's preformatted, md5-verified BLAST+ nr volumes with the official
`update_blastdb.pl` tool and runs `diamond prepdb` so DIAMOND can search them
directly; BLAST v5 databases carry taxonomy, so the `staxids` the NR-validation
step needs are populated with no extra download.

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates/uses `{path}/diamond/`. |
| `--db-name` | `nr` | BLAST database name to fetch with `update_blastdb.pl`. |
| `--source` | `ncbi` | Download mirror (`ncbi`, `aws`, `gcp`; aws/gcp are often faster). |
| `--threads` | `1` | Threads for `update_blastdb.pl` and `diamond prepdb`. |
| `--skip-prepdb` | off | Download only; let the pipeline run `diamond prepdb` at runtime. |
| `--from-blastdb` | *(none)* | Bring-your-own: register an existing BLAST+ nr prefix; skip download. |
| `--from-fasta` | *(none)* | Bring-your-own: build `nr.dmnd` from a protein FASTA via `diamond makedb`. |
| `--taxonmap` / `--taxonnodes` / `--taxonnames` | *(none)* | With `--from-fasta`: `prot.accession2taxid` + taxdump `nodes.dmp`/`names.dmp` so `staxids` populate. |

> **nr is very large** (100+ GB decompressed) and can take hours to download.
> It is intentionally **not** part of `get-databases all`.

```bash
# Default: download preformatted BLAST+ nr and prepare it for DIAMOND
viralunity get-databases nr --path databases/ --source gcp --threads 8
# use: --nr-diamond-database databases/diamond/nr

# Bring-your-own: register an nr you already downloaded
viralunity get-databases nr --from-blastdb /data/blast/nr

# Bring-your-own: build a native .dmnd from a FASTA (with taxonomy for staxids)
viralunity get-databases nr --from-fasta nr.faa \
    --taxonmap prot.accession2taxid.FULL.gz \
    --taxonnodes taxdump/nodes.dmp --taxonnames taxdump/names.dmp
# use: --nr-diamond-database databases/diamond/nr.dmnd
```

### `virus-genome` — Viral genomes database

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/virus_genomes/`. |
| `--taxon` | Viruses | NCBI taxon name to download (e.g. `Viruses`, `coronaviridae`). |
| `--refseq`/`--no-refseq` | `on` | Limit to RefSeq genomes only. |
| `--skip-makeblastdb` | off | Download and reformat files only; skip `makeblastdb` index creation. |

Requires the NCBI Datasets CLI and BLAST+ (`makeblastdb`).

After downloading and reformatting genome sequences, the command automatically runs `makeblastdb` to build a nucleotide BLAST index alongside the FASTA. This index is required when using the `similarity` reference selection strategy (`--reference-selection-strategy similarity`).

```bash
viralunity get-databases virus-genome --taxon Viruses
# FASTA at databases/virus_genomes/viral.genomes.fasta
# BLAST DB index files at databases/virus_genomes/viral.genomes.fasta.{nhr,nin,nsq,...}
# taxids at databases/virus_genomes/genome2taxid.tsv
# use: --viral-genomes databases/virus_genomes/viral.genomes.fasta
#      --viral-taxids databases/virus_genomes/genome2taxid.tsv
```

### `host-genome` — Host genome download

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/host_genomes/`. |
| `--accession` | *(required)* | NCBI genome accession ID (e.g. `GCA_000001405.29`). |

```bash
viralunity get-databases host-genome --accession GCA_000001405.29
# fasta at databases/host_genomes/GCA_000001405.29.fasta
```

### `deacon-index` — Pre-built Deacon index

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/deacon_indexes/`. |
| `--index-name` | `panhuman-1` | Pre-built index name (`panhuman-1` or `panmouse-1`). |

```bash
viralunity get-databases deacon-index --index-name panhuman-1
# index at databases/deacon_indexes/panhuman-1.idx
# use: --deacon-index databases/deacon_indexes/panhuman-1.idx
```

### Download all databases at once

```bash
viralunity get-databases all --threads 4
```

---

## `viralunity build-deacon-index`

Build a Deacon minimizer index from a host genome FASTA.

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/deacon_indexes/`. |
| `--input` | *(required)* | Input FASTA file to index. |
| `--threads` | `8` | Number of threads to use. |

```bash
viralunity build-deacon-index \
    --input databases/host_genomes/GCA_000001405.29.fasta \
    --threads 4
# output at databases/deacon_indexes/GCA_000001405.29.idx
```

---

## Configuration file overrides

A few tool-level parameters are tunable only through the YAML config file produced by `--config-file` (they are not exposed as CLI flags because they rarely need to change). The defaults preserve the historical behaviour, so most users can ignore this section.

Open the generated YAML config file after running with `--create-config-only` and edit the corresponding key under the `# parameters` section:

| Config key | Default | Effect |
|------------|---------|--------|
| `minimap2_consensus_align_flags` | `-a --sam-hit-only --secondary=no --score-N=0` | Flags passed to `minimap2` when re-aligning the per-sample consensus back to the reference for the final multiple-sequence alignment. Used by the `consensus` pipelines. |
| `diamond_max_target_seqs` | `1` | Value of DIAMOND `--max-target-seqs`. Raise it if you need multiple protein hits per query (will slow the run and increase output size). Used by `meta`. |
| `kraken2_extra_flags` | `--report-minimizer-data` | Extra flags appended to every Kraken2 invocation alongside `--threads` and `--minimum-hit-group`. Set to `""` to drop the minimizer-data column from the report. Used by `meta`. |

Example: to bump DIAMOND to keep up to five protein hits per query and drop Kraken2's minimizer column, edit the YAML config to:

```yaml
diamond_max_target_seqs: 5
kraken2_extra_flags: ""
```

Then re-run the same `viralunity meta ...` invocation (without `--create-config-only`) and pass the edited config to Snakemake.

---

## Per-rule CPU and RAM overrides

Both `viralunity consensus` and `viralunity meta` auto-generate a `--<rule>-cpus` and a `--<rule>-ram` option for each computationally significant Snakemake rule, so you can size individual steps without touching the global `--threads` / `--threads-total`. Every such option defaults to **2** CPUs / **4** GB and is written into the generated YAML as `{rule}_cpus` / `{rule}_ram`.

The exact set depends on the pipeline and data type — run the relevant subcommand's `--help` to list them all. Representative rules (flags replace `_` with `-` and append `-cpus` / `-ram`):

- **consensus** — `map_reads`, `trim_primer_sequences`, `perform_qc` + `detect_isnv` (Illumina), `infer_consensus_sequence` (Nanopore). e.g. `--map-reads-cpus`, `--map-reads-ram`.
- **meta** — `remove_host_reads`, `run_megahit`, `run_kraken2_reads`, `run_kraken2_contigs`, `run_diamond_reads`, `run_diamond_contigs`, `run_diamond_nr`, `perform_qc` (Illumina), `run_racon` / `run_medaka` (Nanopore), … e.g. `--run-megahit-cpus`, `--run-kraken2-reads-ram`.

```bash
viralunity meta illumina ... \
    --run-megahit-cpus 8 --run-megahit-ram 32 \
    --run-kraken2-reads-cpus 8
```
