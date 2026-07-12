if run_denovo and run_k2_contigs:
    rule run_kraken2_contigs:
        input:
            fasta = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa"
        output:
            report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.report.txt",
            outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.txt",
        threads: config.get("run_kraken2_contigs_cpus", 2)
        resources:
            mem_mb = config.get("run_kraken2_contigs_ram", 4) * 1024
        params:
            database = config["kraken2_database"],
            minimum_hit_group = config.get("minimum_hit_group", 4),
            extra_flags = config.get("kraken2_extra_flags", "--report-minimizer-data")
        log:
            config["output"] + "logs/kraken2_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/kraken2_contigs/{sample}.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            if [ ! -s {input.fasta} ]; then
                echo "WARNING: {input.fasta} empty. Creating dummy Kraken2 contigs outputs." >> {log}
                touch {output.report} {output.outfile}
            else
                kraken2 --db {params.database} --threads {threads} {params.extra_flags} \
                    --minimum-hit-group {params.minimum_hit_group} --report {output.report} \
                    --output {output.outfile} {input.fasta} 2> {log}
            fi
            """

    rule create_krona_input_from_kraken2_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.txt"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt"
        params:
            keep_columns = [1, 2],
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_taxids.py"

    rule create_krona_report_from_kraken2_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/reports/{sample}.output.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_kraken2_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/krona_kraken2_contigs/{sample}.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty krona input (kraken2 contigs)." >> {log}
                touch {output}
            fi
            """

    rule summarize_taxa_kraken2_contigs:
        input:
            krona = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt",
            plot = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/reports/{sample}.output.krona.html"
        output:
            temp(config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summary/{sample}.taxa.tsv")
        params:
            taxdump = config["taxdump"],
            tool = "kraken2",
            mode = "contigs",
            sample = "{sample}"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/summarize_krona_taxa.py"

    rule summarize_taxa_kraken2_contigs_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summaries/full/kraken2_contigs_taxa_summary.tsv"
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

    rule add_RPM_to_kraken2_contigs_summary:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summaries/full/kraken2_contigs_taxa_summary.tsv",
            merged_fastqs = expand(
                config["output"] + "host_filtered/{sample}.merged.fastq.gz",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summaries/full/kraken2_contigs_taxa_summary_RPM.tsv"
        params:
            sample_to_fastq = get_sample_to_fastq(),
            reads_col = "count"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_RPM_to_summary.py"

    if compute_rpkm:
        rule add_rpkm_to_kraken2_contigs_summary:
            input:
                summary = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summaries/full/kraken2_contigs_taxa_summary_RPM.tsv",
                genome_lengths = config["output"] + "metagenomics/genome_lengths.tsv",
            output:
                config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summaries/full/kraken2_contigs_taxa_summary_RPKM.tsv",
            conda:
                "../envs/utils.yaml"
            script:
                "../python/add_rpkm_to_summary.py"

        rule extract_viral_contigs_kraken2:
            input:
                contigs = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
                krona = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt"
            output:
                fasta = config["output"] + "denovo_assembly/viral_contigs_kraken2/{sample}.viral_contigs.fa",
                ids = temp(config["output"] + "denovo_assembly/viral_contigs_kraken2/{sample}.viral.ids.txt")
            params:
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/extract_viral_contigs.py"

        rule remap_and_depth_viral_contigs_kraken2:
            input:
                fasta = config["output"] + "denovo_assembly/viral_contigs_kraken2/{sample}.viral_contigs.fa",
                R1 = rules.remove_host_reads.output.filtered_R1,
                R2 = rules.remove_host_reads.output.filtered_R2
            output:
                bam = config["output"] + "mapping/viral_kraken2/{sample}.viral.bam",
                bai = config["output"] + "mapping/viral_kraken2/{sample}.viral.bam.bai",
                depth = config["output"] + "mapping/viral_kraken2/{sample}.viral.depth.txt"
            threads: config.get("remap_reads_to_viral_contigs_cpus", 2)
            resources:
                mem_mb = config.get("remap_reads_to_viral_contigs_ram", 4) * 1024
            log:
                config["output"] + "logs/remap_viral_kraken2/{sample}.log"
            conda:
                "../envs/alignment.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.bam}) $(dirname {log})
                if [ ! -s {input.fasta} ]; then
                    echo "No viral kraken2 contigs for {wildcards.sample}; skipping remap." >> {log}
                    touch {output.bam} {output.bai}
                    : > {output.depth}
                else
                    minimap2 -t {threads} -ax sr {input.fasta} {input.R1} {input.R2} | \
                        samtools sort -@ {threads} -o {output.bam} -
                    samtools index -@ {threads} {output.bam}
                    samtools depth -a {output.bam} > {output.depth}
                fi
                """

        rule add_contig_stats_kraken2_contigs:
            input:
                summary = chain_input("kraken2_contigs", "ctgstats"),
                depth = expand(
                    config["output"] + "mapping/viral_kraken2/{sample}.viral.depth.txt",
                    sample=list(config["samples"])
                ),
                krona = expand(
                    config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt",
                    sample=list(config["samples"])
                )
            output:
                summary = chain_output("kraken2_contigs", "ctgstats")
            params:
                samples = list(config["samples"]),
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/add_contig_stats_to_summary.py"

    if run_ictv_host_filter:
        rule apply_ictv_filter_kraken2_contigs:
            input:
                summary = chain_input("kraken2_contigs", "ictv")
            output:
                summary = chain_output("kraken2_contigs", "ictv"),
                dropped = dropped_sidecar(chain_output("kraken2_contigs", "ictv"))
            params:
                allowlist = config.get("ictv_vertebrate_taxids_file", "NA"),
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/apply_ictv_host_filter.py"

    if run_nr_validation:
        rule harmonize_nr_kraken2_contigs:
            input:
                summary = chain_input("kraken2_contigs", "nr"),
                nr = config["output"] + "metagenomics/nr_validation/nr_query.nr.top_species_hit_lca.tsv",
                krona = expand(config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt", sample=list(config["samples"]))
            output:
                summary = chain_output("kraken2_contigs", "nr"),
                dropped = dropped_sidecar(chain_output("kraken2_contigs", "nr")),
                flags = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summaries/full/kraken2_contigs_nr_flags.tsv"
            params:
                samples = list(config["samples"]),
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/harmonize_nr_summary.py"

    rule apply_bleed_filter_kraken2_contigs:
        input:
            chain_input("kraken2_contigs", "bleed")
        output:
            chain_output("kraken2_contigs", "bleed")
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
        rule add_negative_control_enrichment_kraken2_contigs:
            input:
                chain_input("kraken2_contigs", "neg")
            output:
                chain_output("kraken2_contigs", "neg")
            params:
                negatives = config.get("negative_controls", []),
                pseudocount = config.get("enrichment_pseudocount", 1.0),
                z_score_threshold = config.get("z_score_threshold", 3.0),
                log10_ratio_threshold = config.get("log10_ratio_threshold", 1.0)
            conda:
                "../envs/utils.yaml"
            script:
                "../python/add_negative_control_enrichment.py"

    rule make_filtered_krona_input_kraken2_contigs:
        input:
            summary = final_summary("kraken2_contigs"),
            krona_input = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.filtered.krona.txt"
        params:
            sample = "{sample}",
            tool = "kraken2",
            mode = "contigs",
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_krona_by_pass_taxids.py"

    rule create_filtered_krona_report_from_kraken2_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.filtered.krona.txt"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/reports/{sample}.output.filtered.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_kraken2_contigs/{sample}.filtered.log"
        benchmark:
            config["output"] + "logs/krona_kraken2_contigs/{sample}.filtered.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            r"""
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty filtered krona input (kraken2 contigs)." >> {log}
                touch {output}
            fi
            """
