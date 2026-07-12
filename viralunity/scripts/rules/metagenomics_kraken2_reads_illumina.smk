rule run_kraken2_reads:
    input:
        filtered_R1 = rules.remove_host_reads.output.filtered_R1,
        filtered_R2 = rules.remove_host_reads.output.filtered_R2,
    output:
        report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.report.txt",
        outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.txt",
    threads: config.get("run_kraken2_reads_cpus", 2)
    resources:
        mem_mb = config.get("run_kraken2_reads_ram", 4) * 1024
    params:
        database = config["kraken2_database"],
        minimum_hit_group = config.get("minimum_hit_group", 4),
        extra_flags = config.get("kraken2_extra_flags", "--report-minimizer-data")
    log:
        config["output"] + "logs/kraken2_reads/{sample}.log"
    benchmark:
        config["output"] + "logs/kraken2_reads/{sample}.benchmark.txt"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.report}) $(dirname {log})
        unc_size_R1=$(gzip -l "{input.filtered_R1}" | awk 'NR==2 {{print $1}}')
        unc_size_R2=$(gzip -l "{input.filtered_R2}" | awk 'NR==2 {{print $1}}')
        if [ "$unc_size_R1" = "0" ] || [ "$unc_size_R2" = "0" ]; then
            echo "WARNING: {input.filtered_R1} or {input.filtered_R2} empty. Creating dummy Kraken2 READS outputs." > {log}
            : > {output.report}
            : > {output.outfile}
        else
            kraken2 --db {params.database} --threads {threads} {params.extra_flags} \
                --minimum-hit-group {params.minimum_hit_group} --report {output.report} \
                --output {output.outfile} {input.filtered_R1} {input.filtered_R2} 2> {log}
        fi
        """

rule create_krona_input_from_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.txt"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt"
    params:
        keep_columns = [1, 2],
        taxid_column = 1,
        exclude_taxids = EXCLUDE_TAXIDS
    conda:
        "../envs/utils.yaml"
    script:
        "../python/filter_taxids.py"

rule create_krona_report_from_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/reports/{sample}.output.krona.html"
    params:
        krona_database = config["krona_database"]
    log:
        config["output"] + "logs/krona_kraken2_reads/{sample}.log"
    benchmark:
        config["output"] + "logs/krona_kraken2_reads/{sample}.benchmark.txt"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        set -euo pipefail
        if [ -s {input} ]; then
            ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
        else
            echo "Empty krona input (kraken2 reads)." >> {log}
            touch {output}
        fi
        """

rule summarize_taxa_kraken2_reads:
    input:
        krona = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt",
        plot = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/reports/{sample}.output.krona.html"
    output:
        temp(config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summary/{sample}.taxa.tsv")
    params:
        taxdump = config["taxdump"],
        tool = "kraken2",
        mode = "reads",
        sample = "{sample}"
    conda:
        "../envs/utils.yaml"
    script:
        "../python/summarize_krona_taxa.py"

rule summarize_taxa_kraken2_reads_all:
    input:
        expand(
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summary/{sample}.taxa.tsv",
            sample=config["samples"]
        )
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summaries/full/kraken2_reads_taxa_summary.tsv"
    conda:
        "../envs/utils.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})
        header_file=""
        for f in {input}; do [ -s "$f" ] && header_file="$f" && break; done
        if [ -n "$header_file" ]; then head -n 1 "$header_file" > {output}; else echo -e "sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource" > {output}; fi
        for f in {input}; do [ -s "$f" ] && tail -n +2 "$f" >> {output}; done
        """

rule add_RPM_to_kraken2_reads_summary:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summaries/full/kraken2_reads_taxa_summary.tsv",
        merged_fastqs = expand(
            config["output"] + "host_filtered/{sample}.merged.fastq.gz",
            sample=config["samples"]
        )
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summaries/full/kraken2_reads_taxa_summary_RPM.tsv"
    params:
        sample_to_fastq = get_sample_to_fastq(),
        reads_col = "count"
    conda:
        "../envs/utils.yaml"
    script:
        "../python/add_RPM_to_summary.py"

if compute_rpkm:
    rule add_rpkm_to_kraken2_reads_summary:
        input:
            summary = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summaries/full/kraken2_reads_taxa_summary_RPM.tsv",
            genome_lengths = config["output"] + "metagenomics/genome_lengths.tsv",
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summaries/full/kraken2_reads_taxa_summary_RPKM.tsv",
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_rpkm_to_summary.py"

if run_ictv_host_filter:
    rule apply_ictv_filter_kraken2_reads:
        input:
            summary = chain_input("kraken2_reads", "ictv")
        output:
            summary = chain_output("kraken2_reads", "ictv"),
            dropped = dropped_sidecar(chain_output("kraken2_reads", "ictv"))
        params:
            allowlist = config.get("ictv_vertebrate_taxids_file", "NA"),
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_ictv_host_filter.py"


rule apply_bleed_filter_kraken2_reads:
    input:
        chain_input("kraken2_reads", "bleed")
    output:
        chain_output("kraken2_reads", "bleed")
    params:
        fraction = config.get("bleed_fraction", 0.005),
        rpm_floor = config.get("bleed_rpm_floor", 1.0),
        rpkm_floor = config.get("bleed_rpkm_floor", 0.1),
        rpm_col = "rpm",
    conda:
        "../envs/utils.yaml"
    script:
        "../python/apply_max_rpm_bleed_filter.py"

if has_negative_controls:
    rule add_negative_control_enrichment_kraken2_reads:
        input:
            chain_input("kraken2_reads", "neg")
        output:
            chain_output("kraken2_reads", "neg")
        params:
            negatives = config.get("negative_controls", []),
            pseudocount = config.get("enrichment_pseudocount", 1.0),
            z_score_threshold = config.get("z_score_threshold", 3.0),
            log10_ratio_threshold = config.get("log10_ratio_threshold", 1.0)
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_negative_control_enrichment.py"

rule make_filtered_krona_input_kraken2_reads:
    input:
        summary = final_summary("kraken2_reads"),
        krona_input = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.filtered.krona.txt"
    params:
        sample = "{sample}",
        tool = "kraken2",
        mode = "reads",
        taxdump = config["taxdump"]
    conda:
        "../envs/utils.yaml"
    script:
        "../python/filter_krona_by_pass_taxids.py"

rule create_filtered_krona_report_from_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.filtered.krona.txt"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/reports/{sample}.output.filtered.krona.html"
    params:
        krona_database = config["krona_database"]
    log:
        config["output"] + "logs/krona_kraken2_reads/{sample}.filtered.log"
    benchmark:
        config["output"] + "logs/krona_kraken2_reads/{sample}.filtered.benchmark.txt"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        r"""
        set -euo pipefail
        if [ -s {input} ]; then
            ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
        else
            echo "Empty filtered krona input (kraken2 reads)." >> {log}
            touch {output}
        fi
        """
