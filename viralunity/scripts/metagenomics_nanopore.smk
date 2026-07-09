def get_exclude_taxids():
    """TaxIDs to exclude from classification outputs (e.g. human 9606, unclassified 0)."""
    exclude = []
    if config.get("remove_human_reads", False):
        exclude.append("9606")
    if config.get("remove_unclassified_reads", False):
        exclude.append("0")
    return exclude

EXCLUDE_TAXIDS = get_exclude_taxids()

def get_sample_to_fastq():
    """Map each sample to its host-filtered FASTQ for read counting (RPM)."""
    return {s: config["output"] + "host_filtered/" + s + ".filtered.fastq.gz" for s in config["samples"]}

def get_map_input_fastqs(wildcards):
    """Return list of read files for this sample (single-end long reads)."""
    paths = config["samples"][wildcards.sample]
    if isinstance(paths, str):
        paths = paths.strip().split()
    return paths

host_filtering_enabled = (config.get("host_reference", "NA") not in ("NA", "", None)) or (config.get("deacon_index", "NA") not in ("NA", "", None))
dehost_with_deacon = config.get("deacon_index", "NA") not in ("NA", "", None)
run_denovo = config.get("run_denovo_assembly", False)
run_polish_racon = config.get("run_polish_racon", False)
run_polish_medaka = config.get("run_polish_medaka", False)
run_k2_reads = config.get("run_kraken2_reads", True)
run_k2_contigs = config.get("run_kraken2_contigs", True)
run_diamond_reads = config.get("run_diamond_reads", False)
run_diamond_contigs = config.get("run_diamond_contigs", False)
has_negative_controls = bool(config.get("negative_controls", []))
compute_rpkm = bool(config.get("compute_rpkm", False))
combine_contig_search = bool(config.get("combine_contig_search", False))
run_ictv_host_filter = bool(config.get("run_ictv_host_filter", False))
run_nr_validation = bool(config.get("run_nr_validation", False))

diamond_db_input_path = config.get("diamond_database", "NA")
if diamond_db_input_path != "NA":
    diamond_db_is_ready = diamond_db_input_path.endswith(".dmnd")
    diamond_db_file = diamond_db_input_path if diamond_db_is_ready else diamond_db_input_path + ".dmnd"
else:
    diamond_db_is_ready = False
    diamond_db_file = "NA"

def _classification_contigs_name():
    if run_polish_medaka:
        return "polished.fasta"
    if run_polish_racon:
        return "racon.fasta"
    return "final.contigs.fa"

def get_final_contigs(wildcards):
    """Path to contigs used for classification (polished, racon, or raw MEGAHIT)."""
    base = config["output"] + "denovo_assembly/megahit/{sample}/"
    return base.format(sample=wildcards.sample) + _classification_contigs_name()

def all_classification_contigs():
    """Per-sample classification-contig paths, in config['samples'] order, for the
    aggregated combined search."""
    base = config["output"] + "denovo_assembly/megahit/{sample}/"
    name = _classification_contigs_name()
    return [base.format(sample=s) + name for s in config["samples"]]

def get_medaka_assembly_input(wildcards):
    """Assembly input for Medaka (racon output or MEGAHIT)."""
    base = config["output"] + "denovo_assembly/megahit/{sample}/"
    s = wildcards.sample
    if run_polish_racon:
        return base.format(sample=s) + "racon.fasta"
    return base.format(sample=s) + "final.contigs.fa"

def _summary_stem(track):
    return config["output"] + "metagenomics/taxonomic_assignments/" + track + "/" + track + "_taxa_summary"

def _summary_base(track):
    return _summary_stem(track) + ("_RPKM" if compute_rpkm else "_RPM")

def _chain_steps(track):
    """Ordered post-metric processing steps enabled for this track.

    The taxa summary flows through one cumulative chain:
    base -> _RPM/_RPKM -> [.nr] -> .bleed -> [.neg] -> [.ictv]. Each enabled step
    appends exactly one filename suffix. NR validation applies to contig tracks
    only; bleed is always on; neg runs when negative controls are present; the
    ICTV host filter applies to all tracks and runs last.
    """
    steps = []
    if run_nr_validation and track.endswith("contigs"):
        steps.append("nr")
    steps.append("bleed")
    if has_negative_controls:
        steps.append("neg")
    if run_ictv_host_filter:
        steps.append("ictv")
    return steps

def _chain_path(track, upto):
    """Active-metric base + '.'-joined step suffixes up to index ``upto`` (-1 = base only)."""
    steps = _chain_steps(track)
    return _summary_base(track) + "".join("." + s for s in steps[: upto + 1]) + ".tsv"

def chain_input(track, step):
    """Summary a chain step consumes (the prior step's output, or the metric base)."""
    return _chain_path(track, _chain_steps(track).index(step) - 1)

def chain_output(track, step):
    """Summary a chain step produces (base + suffixes up to and including it)."""
    return _chain_path(track, _chain_steps(track).index(step))

def dropped_sidecar(path):
    """Audit path for rows a row-removing step dropped: '<x>.tsv' -> '<x>.dropped.tsv'."""
    return path[:-4] + ".dropped.tsv" if path.endswith(".tsv") else path + ".dropped.tsv"

def final_summary(track):
    """The last file in the chain for a track (the fully-filtered summary)."""
    return _chain_path(track, len(_chain_steps(track)) - 1)

def _all_inputs():
    targets = [
        config["output"] + "benchmark.tsv"
    ]
    if run_k2_reads:
        targets.append(final_summary("kraken2_reads"))
        targets.extend(expand(
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/reports/{sample}.output.filtered.krona.html",
            sample=config["samples"],
        ))
    if run_denovo and run_k2_contigs:
        targets.append(final_summary("kraken2_contigs"))
        targets.extend(expand(
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/reports/{sample}.output.filtered.krona.html",
            sample=config["samples"],
        ))
    if run_diamond_reads:
        targets.append(final_summary("diamond_reads"))
        targets.extend(expand(
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/reports/{sample}.diamond.filtered.krona.html",
            sample=config["samples"],
        ))
    if run_denovo and run_diamond_contigs:
        targets.append(final_summary("diamond_contigs"))
        targets.extend(expand(
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/reports/{sample}.diamond.supported.filtered.krona.html",
            sample=config["samples"],
        ))

    if run_denovo and run_diamond_contigs and run_nr_validation:
        targets.append(config["output"] + "metagenomics/nr_validation/nr_query.nr.top_species_hit_lca.viruses_only.tsv")

    if config.get("run_reference_assembly", False):
        targets.append(config["output"] + "reference_assembly_done.txt")
    return targets

rule all:
    input:
        _all_inputs(),

if (run_diamond_reads or run_diamond_contigs) and not diamond_db_is_ready and diamond_db_input_path != "NA":
    rule create_diamond_db_shared:
        input:
            diamond_db_input_path
        output:
            diamond_db_file
        log:
            config["output"] + "logs/diamond/diamond_makedb.log"
        benchmark:
            config["output"] + "logs/diamond/diamond_makedb.benchmark.log"
        conda:
            "envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output}) $(dirname {log})
            diamond makedb --in {input} --db {output} 2> {log}
            """
include: "rules/metagenomics_genome_lengths.smk"
include: "rules/metagenomics_dehost_nanopore.smk"
include: "rules/metagenomics_kraken2_reads_nanopore.smk"
include: "rules/metagenomics_diamond_reads_nanopore.smk"
include: "rules/metagenomics_assembly_nanopore.smk"
include: "rules/metagenomics_kraken2_contigs_nanopore.smk"
include: "rules/metagenomics_diamond_contigs_nanopore.smk"
include: "rules/metagenomics_nr_validation_nanopore.smk"
if config.get("run_reference_assembly", False):
    include: "rules/metagenomics_reference_assembly.smk"

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        kraken2_reads_reports = expand(rules.run_kraken2_reads.output.report, sample=config["samples"]) if run_k2_reads else [],
        kraken2_reads_krona = expand(rules.create_krona_report_from_kraken2_reads.output, sample=config["samples"]) if run_k2_reads else [],
        kraken2_reads_filtered_krona = expand(rules.create_filtered_krona_report_from_kraken2_reads.output, sample=config["samples"]) if run_k2_reads else [],
        diamond_reads_tsv = expand(rules.run_diamond_reads.output.tsv, sample=config["samples"]) if run_diamond_reads else [],
        diamond_reads_krona = expand(rules.create_krona_report_diamond_reads.output, sample=config["samples"]) if run_diamond_reads else [],
        diamond_reads_filtered_krona = expand(rules.create_filtered_krona_report_diamond_reads.output, sample=config["samples"]) if run_diamond_reads else [],
        host_filtered = expand(rules.remove_host_reads.output.filtered, sample=config["samples"]),
        megahit_contigs = expand(rules.run_megahit.output.contigs, sample=config["samples"]) if run_denovo else [],
        kraken2_contigs_reports = expand(rules.run_kraken2_contigs.output.report, sample=config["samples"]) if run_denovo and run_k2_contigs else [],
        kraken2_contigs_krona = expand(rules.create_krona_report_from_kraken2_contigs.output, sample=config["samples"]) if run_denovo and run_k2_contigs else [],
        kraken2_contigs_filtered_krona = expand(rules.create_filtered_krona_report_from_kraken2_contigs.output, sample=config["samples"]) if run_denovo and run_k2_contigs else [],
        diamond_contigs_tsv = expand(rules.run_diamond_contigs.output.tsv, sample=config["samples"]) if run_denovo and run_diamond_contigs else [],
        diamond_contigs_krona = expand(rules.create_krona_report_diamond_contigs.output, sample=config["samples"]) if run_denovo and run_diamond_contigs else [],
        diamond_contigs_filtered_krona = expand(rules.create_filtered_krona_report_diamond_contigs.output, sample=config["samples"]) if run_denovo and run_diamond_contigs else [],
        summary_k2_reads = final_summary("kraken2_reads") if run_k2_reads else [],
        summary_diamond_reads = final_summary("diamond_reads") if run_diamond_reads else [],
        summary_k2_contigs = final_summary("kraken2_contigs") if run_denovo and run_k2_contigs else [],
        summary_diamond_contigs = final_summary("diamond_contigs") if run_denovo and run_diamond_contigs else [],
    output:
        benchmark = config["output"] + "benchmark.tsv"
    params:
        outdir = config["output"],
        samples = list(config["samples"].keys()),
        run_k2_reads = run_k2_reads,
        run_diamond_reads = run_diamond_reads,
        run_denovo = run_denovo,
        run_k2_contigs = run_k2_contigs,
        run_diamond_contigs = run_diamond_contigs
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            mkdir -p {params.outdir}samples/$sample;
        done

        # Kraken2 Reads
        if [ "{params.run_k2_reads}" = "True" ]; then
            for _file in {input.kraken2_reads_reports}; do
                sample=$(basename $_file .report.txt);
                ln -sf $_file {params.outdir}samples/$sample/kraken2_reads.report.txt;
            done
            for _file in {input.kraken2_reads_krona}; do
                sample=$(basename $_file .output.krona.html);
                ln -sf $_file {params.outdir}samples/$sample/kraken2_reads.krona.html;
            done
            for _file in {input.kraken2_reads_filtered_krona}; do
                sample=$(basename $_file .output.filtered.krona.html);
                ln -sf $_file {params.outdir}samples/$sample/kraken2_reads.filtered.krona.html;
            done
        fi

        # DIAMOND Reads
        if [ "{params.run_diamond_reads}" = "True" ]; then
            for _file in {input.diamond_reads_tsv}; do
                sample=$(basename $_file .diamond.tsv);
                ln -sf $_file {params.outdir}samples/$sample/diamond_reads.tsv;
            done
            for _file in {input.diamond_reads_krona}; do
                sample=$(basename $_file .diamond.krona.html);
                ln -sf $_file {params.outdir}samples/$sample/diamond_reads.krona.html;
            done
            for _file in {input.diamond_reads_filtered_krona}; do
                sample=$(basename $_file .diamond.filtered.krona.html);
                ln -sf $_file {params.outdir}samples/$sample/diamond_reads.filtered.krona.html;
            done
        fi

        # De novo assembly
        if [ "{params.run_denovo}" = "True" ]; then
            for _file in {input.megahit_contigs}; do
                sample=$(basename $(dirname $_file));
                ln -sf $_file {params.outdir}samples/$sample/contigs.fa;
            done
        fi

        # Taxonomy (Contigs)
        if [ "{params.run_denovo}" = "True" ]; then
            if [ "{params.run_k2_contigs}" = "True" ]; then
                for _file in {input.kraken2_contigs_reports}; do
                    sample=$(basename $_file .report.txt);
                    ln -sf $_file {params.outdir}samples/$sample/kraken2_contigs.report.txt;
                done
                for _file in {input.kraken2_contigs_krona}; do
                    sample=$(basename $_file .output.krona.html);
                    ln -sf $_file {params.outdir}samples/$sample/kraken2_contigs.krona.html;
                done
                for _file in {input.kraken2_contigs_filtered_krona}; do
                    sample=$(basename $_file .output.filtered.krona.html);
                    ln -sf $_file {params.outdir}samples/$sample/kraken2_contigs.filtered.krona.html;
                done
            fi
            if [ "{params.run_diamond_contigs}" = "True" ]; then
                for _file in {input.diamond_contigs_tsv}; do
                    sample=$(basename $_file .diamond.tsv);
                    ln -sf $_file {params.outdir}samples/$sample/diamond_contigs.tsv;
                done
                for _file in {input.diamond_contigs_krona}; do
                    sample=$(basename $_file .diamond.supported.krona.html);
                    ln -sf $_file {params.outdir}samples/$sample/diamond_contigs.krona.html;
                done
                for _file in {input.diamond_contigs_filtered_krona}; do
                    sample=$(basename $_file .diamond.supported.filtered.krona.html);
                    ln -sf $_file {params.outdir}samples/$sample/diamond_contigs.filtered.krona.html;
                done
            fi
        fi

        # Host filtered reads (Nanopore)
        for _file in {input.host_filtered}; do
            sample=$(basename $_file .filtered.fastq.gz);
            ln -sf $_file {params.outdir}samples/$sample/host_filtered_reads.fastq.gz;
        done

        # Benchmark aggregation
        echo -e "sample\\ttask\\tseconds\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output}
        find {params.outdir} -name "*.benchmark.txt" | while read -r file; do
            task=$(basename $(dirname $file))
            sample=$(basename "$file" .benchmark.txt)
            
            matched=false
            for s in {params.samples}; do
                if [[ "$sample" == "$s" ]]; then
                    matched=true
                    break
                fi
            done
            
            if [[ "$matched" == "false" ]]; then
                sample="All"
            else
                sample=$(echo "$sample" | sed 's/sample-//')
            fi

            tail -n +2 "$file" | awk -v sample="$sample" -v task="$task" '{{print sample"\\t"task"\\t"$0}}' >> {output}
        done
        """
