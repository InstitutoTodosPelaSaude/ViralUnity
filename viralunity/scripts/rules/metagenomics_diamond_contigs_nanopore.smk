if run_denovo and run_diamond_contigs:
    # Contig DIAMOND search: either one search per sample (default) or a single
    # aggregated search over all samples' contigs (sample-prefixed headers) that
    # is split back per sample. Toggled by config["combine_contig_search"].
    if combine_contig_search:
        rule combine_diamond_contigs_query:
            input:
                fastas = all_classification_contigs()
            output:
                config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/_combined.contigs.fasta"
            params:
                samples = list(config["samples"])
            conda:
                "../envs/utils.yaml"
            script:
                "../python/combine_contigs.py"

        rule run_diamond_contigs_combined:
            input:
                fasta = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/_combined.contigs.fasta",
                db = diamond_db_file
            output:
                tsv = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/_combined.diamond.tsv"
            threads: config.get("run_diamond_contigs_cpus", 2)
            resources:
                mem_mb = config.get("run_diamond_contigs_ram", 4) * 1024
            log:
                config["output"] + "logs/diamond_contigs/_combined.log"
            benchmark:
                config["output"] + "logs/diamond_contigs/_combined.benchmark.txt"
            params:
                sensitivity = config.get("diamond_sensitivity", "sensitive"),
                evalue = config.get("evalue", 0.001),
                max_target_seqs = config.get("diamond_max_target_seqs", 1)
            conda:
                "../envs/taxonomy.yaml"
            shell:
                r"""
                set -euo pipefail
                if [ ! -s {input.fasta} ]; then
                    echo "WARNING: {input.fasta} empty. Creating dummy combined DIAMOND output." >> {log}
                    touch {output.tsv}
                else
                    diamond blastx --db {input.db} --query {input.fasta} \
                        --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen qseq qseq_translated full_qseq sseq full_sseq \
                        --max-target-seqs {params.max_target_seqs} --evalue {params.evalue} \
                        --{params.sensitivity} --threads {threads} 2> {log}
                fi
                """

        rule split_diamond_contigs:
            input:
                config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/_combined.diamond.tsv"
            output:
                expand(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv", sample=list(config["samples"]))
            params:
                samples = list(config["samples"])
            conda:
                "../envs/utils.yaml"
            script:
                "../python/split_search_output.py"
    else:
        rule run_diamond_contigs:
            input:
                fasta = get_final_contigs,
                db = diamond_db_file
            output:
                tsv = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv"
            threads: config.get("run_diamond_contigs_cpus", 2)
            resources:
                mem_mb = config.get("run_diamond_contigs_ram", 4) * 1024
            log:
                config["output"] + "logs/diamond_contigs/{sample}.log"
            benchmark:
                config["output"] + "logs/diamond_contigs/{sample}.benchmark.txt"
            params:
                sensitivity = config.get("diamond_sensitivity", "sensitive"),
                evalue = config.get("evalue", 0.001),
                max_target_seqs = config.get("diamond_max_target_seqs", 1)
            conda:
                "../envs/taxonomy.yaml"
            shell:
                r"""
                set -euo pipefail
                if [ ! -s {input.fasta} ]; then
                    echo "WARNING: {input.fasta} empty. Creating dummy DIAMOND contigs output." >> {log}
                    touch {output.tsv}
                else
                    diamond blastx --db {input.db} --query {input.fasta} \
                        --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen qseq qseq_translated full_qseq sseq full_sseq \
                        --max-target-seqs {params.max_target_seqs} --evalue {params.evalue} \
                        --{params.sensitivity} --threads {threads} 2> {log}
                fi
                """

    rule extract_viral_contigs:
        input:
            contigs = get_final_contigs,
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv"
        output:
            # NOT temp(): remap_reads_to_viral_contigs reuses this id list
            # downstream (see rules.extract_viral_contigs.output.ids), so unlike
            # the Illumina track it must persist for the rest of the DAG.
            ids = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral.ids.txt",
            fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa"
        threads: 1
        log:
            config["output"] + "logs/viral_contigs/{sample}.log"
        conda:
            "../envs/utils.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.fasta}) $(dirname {log})
            cut -f1 {input.diamond} | sort -u > {output.ids}
            if [ -s {output.ids} ]; then
                seqtk subseq {input.contigs} {output.ids} > {output.fasta}
            else
                echo "No viral contigs for {wildcards.sample}" >> {log}
                touch {output.fasta}
            fi
            """

    rule index_viral_contigs:
        input:
            fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa"
        output:
            index = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.mmi"
        threads: config.get("index_viral_contigs_cpus", 2)
        resources:
            mem_mb = config.get("index_viral_contigs_ram", 4) * 1024
        conda:
            "../envs/alignment.yaml"
        shell:
            """
            set -euo pipefail
            if [ -s {input.fasta} ]; then
                minimap2 -d {output.index} {input.fasta}
            else
                touch {output.index}
            fi
            """

    if run_polish_medaka:

        rule bam_sort_index_idxstats_from_medaka:
            input:
                bam = rules.run_medaka.output.bam,
                viral_ids = rules.extract_viral_contigs.output.ids
            output:
                bam = config["output"] + "mapping/viral/{sample}.viral.bam",
                bai = config["output"] + "mapping/viral/{sample}.viral.bam.bai",
                idxstats = config["output"] + "mapping/viral/{sample}.viral.idxstats.txt",
                idxstats_filtered = config["output"] + "mapping/viral/{sample}.viral.idxstats.filtered.txt"
            threads: config.get("bam_sort_index_idxstats_from_medaka_cpus", 2)
            resources:
                mem_mb = config.get("bam_sort_index_idxstats_from_medaka_ram", 4) * 1024
            log:
                config["output"] + "logs/bam_idxstats/{sample}.log"
            benchmark:
                config["output"] + "logs/bam_idxstats/{sample}.benchmark.txt"
            conda:
                "../envs/medaka.yaml"
            shell:
                """
                set -euo pipefail
                mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
                if [ ! -s "{input.bam}" ]; then
                    echo "No Medaka BAM for {wildcards.sample}." > "{log}"
                    touch "{output.bam}" "{output.bai}" "{output.idxstats}" "{output.idxstats_filtered}"
                elif [ ! -s "{input.viral_ids}" ]; then
                    echo "No viral contigs for {wildcards.sample}." > "{log}"
                    touch "{output.bam}" "{output.bai}" "{output.idxstats}" "{output.idxstats_filtered}"
                else
                    refs=()
                    while IFS= read -r r; do [[ -n "$r" ]] && refs+=("$r"); done < "{input.viral_ids}"
                    samtools view -@ {threads} -b "{input.bam}" "$${{refs[@]}}" | samtools sort -@ {threads} -o "{output.bam}" -
                    samtools index -@ {threads} "{output.bam}"
                    samtools idxstats "{output.bam}" > "{output.idxstats}"
                    awk '$3 > 0 && $1 != "*"' "{output.idxstats}" > "{output.idxstats_filtered}"
                fi
                """

    else:

        rule remap_reads_to_viral_contigs:
            input:
                index = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.mmi",
                fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa",
                reads = rules.remove_host_reads.output.filtered
            output:
                bam = config["output"] + "mapping/viral/{sample}.viral.bam",
                bai = config["output"] + "mapping/viral/{sample}.viral.bam.bai",
                idxstats = config["output"] + "mapping/viral/{sample}.viral.idxstats.txt"
            threads: config.get("remap_reads_to_viral_contigs_cpus", 2)
            resources:
                mem_mb = config.get("remap_reads_to_viral_contigs_ram", 4) * 1024
            log:
                config["output"] + "logs/remap_viral/{sample}.log"
            benchmark:
                config["output"] + "logs/remap_viral/{sample}.benchmark.txt"
            conda:
                "../envs/alignment.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.bam}) $(dirname {log})
                if [ ! -s {input.fasta} ]; then
                    echo "No viral contigs for {wildcards.sample}; skipping remap." >> {log}
                    touch {output.bam} {output.bai} {output.idxstats}
                else
                    minimap2 -t {threads} -ax map-ont {input.index} {input.reads} | \
                        samtools sort -@ {threads} -o {output.bam} -
                    samtools index -@ {threads} {output.bam}
                    samtools idxstats {output.bam} > {output.idxstats}
                fi
                """

    rule depth_of_viral_contigs:
        input:
            bam = config["output"] + "mapping/viral/{sample}.viral.bam",
            bai = config["output"] + "mapping/viral/{sample}.viral.bam.bai"
        output:
            depth = config["output"] + "mapping/viral/{sample}.viral.depth.txt"
        log:
            config["output"] + "logs/remap_viral/{sample}.depth.log"
        conda:
            "../envs/alignment.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.depth}) $(dirname {log})
            if [ ! -s {input.bam} ]; then
                : > {output.depth}
            else
                samtools depth -a {input.bam} > {output.depth} 2> {log}
            fi
            """

    rule diamond_filter_by_idxstats:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv",
            idxstats = config["output"] + "mapping/viral/{sample}.viral.idxstats.txt"
        output:
            filtered = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv"
        params:
            min_mapped = config.get("diamond_min_mapped", 1)
        log:
            config["output"] + "logs/diamond_filter/{sample}.log"
        benchmark:
            config["output"] + "logs/diamond_filter/{sample}.benchmark.txt"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_diamond_by_idxstats.py"

    rule annotate_diamond_taxonomy:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv",
            assembly = config["taxids"]
        output:
            annotated = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tax.tsv"
        log:
            config["output"] + "logs/annotate_tax/{sample}.log"
        benchmark:
            config["output"] + "logs/annotate_tax/{sample}.benchmark.txt"
        params:
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/annotate_diamond_taxonomy.py"

    rule create_krona_input_from_diamond_contigs:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv",
            fasta = get_final_contigs,
            assembly = config["taxids"]
        output:
            krona_input = temp(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.temp.tsv")
        params:
            data_format = "fasta"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/convert_diamond_output_to_krona_input.py"

    rule filter_krona_input_from_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.temp.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv"
        params:
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_taxids.py"

    rule create_krona_report_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/reports/{sample}.diamond.supported.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_diamond_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/krona_diamond_contigs/{sample}.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty krona input (diamond contigs)." >> {log}
                touch {output}
            fi
            """

    rule summarize_taxa_diamond_contigs:
        input:
            krona = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv",
            plot = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/reports/{sample}.diamond.supported.krona.html",
            annotated = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tax.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/summary/{sample}.taxa.tsv"
        params:
            taxdump = config["taxdump"],
            tool = "diamond",
            mode = "contigs",
            sample = "{sample}"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/summarize_krona_taxa.py"

    rule summarize_taxa_diamond_contigs_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/chain/diamond_contigs_taxa_summary.tsv"
        conda:
            "../envs/utils.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output})
            header_file=""
            for f in {input}; do [ -s "$f" ] && header_file="$f" && break; done
            if [ -n "$header_file" ]; then head -n 1 "$header_file" > {output}; else echo -e "sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource\tmapped_reads" > {output}; fi
            for f in {input}; do [ -s "$f" ] && tail -n +2 "$f" >> {output}; done
            """

    rule add_RPM_to_diamond_contigs_summary:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/chain/diamond_contigs_taxa_summary.tsv",
            merged_fastqs = expand(
                config["output"] + "host_filtered/{sample}.filtered.fastq.gz",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/chain/diamond_contigs_taxa_summary_RPM.tsv"
        params:
            sample_to_fastq = get_sample_to_fastq(),
            reads_col = "mapped_reads"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_RPM_to_summary.py"

    if compute_rpkm:
        rule add_rpkm_to_diamond_contigs_summary:
            input:
                summary = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/chain/diamond_contigs_taxa_summary_RPM.tsv",
                genome_lengths = config["output"] + "metagenomics/genome_lengths.tsv",
            output:
                config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/chain/diamond_contigs_taxa_summary_RPKM.tsv",
            conda:
                "../envs/utils.yaml"
            script:
                "../python/add_rpkm_to_summary.py"

        rule add_contig_stats_diamond_contigs:
            input:
                summary = chain_input("diamond_contigs", "ctgstats"),
                depth = expand(
                    config["output"] + "mapping/viral/{sample}.viral.depth.txt",
                    sample=list(config["samples"])
                ),
                krona = expand(
                    config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv",
                    sample=list(config["samples"])
                )
            output:
                summary = chain_output("diamond_contigs", "ctgstats")
            params:
                samples = list(config["samples"]),
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/add_contig_stats_to_summary.py"

    if run_ictv_host_filter:
        rule apply_ictv_filter_diamond_contigs:
            input:
                summary = chain_input("diamond_contigs", "ictv")
            output:
                summary = chain_output("diamond_contigs", "ictv"),
                dropped = dropped_sidecar(chain_output("diamond_contigs", "ictv"))
            params:
                allowlist = config.get("ictv_vertebrate_taxids_file", "NA"),
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/apply_ictv_host_filter.py"

    if run_nr_validation:
        rule harmonize_nr_diamond_contigs:
            input:
                summary = chain_input("diamond_contigs", "nr"),
                nr = config["output"] + "metagenomics/nr_validation/nr_query.nr.top_species_hit_lca.tsv",
                krona = expand(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv", sample=list(config["samples"]))
            output:
                summary = chain_output("diamond_contigs", "nr"),
                dropped = dropped_sidecar(chain_output("diamond_contigs", "nr")),
                flags = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_nr_flags.tsv"
            params:
                samples = list(config["samples"]),
                taxdump = config["taxdump"]
            conda:
                "../envs/utils.yaml"
            script:
                "../python/harmonize_nr_summary.py"

    rule apply_bleed_filter_diamond_contigs:
        input:
            chain_input("diamond_contigs", "bleed")
        output:
            chain_output("diamond_contigs", "bleed")
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
        rule add_negative_control_enrichment_diamond_contigs:
            input:
                chain_input("diamond_contigs", "neg")
            output:
                chain_output("diamond_contigs", "neg")
            params:
                negatives = config.get("negative_controls", []),
                pseudocount = config.get("enrichment_pseudocount", 1.0),
                z_score_threshold = config.get("z_score_threshold", 3.0),
                log10_ratio_threshold = config.get("log10_ratio_threshold", 1.0)
            conda:
                "../envs/utils.yaml"
            script:
                "../python/add_negative_control_enrichment.py"

    rule make_filtered_krona_input_diamond_contigs:
        input:
            summary = final_summary("diamond_contigs"),
            krona_input = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.filtered.krona_input.tsv"
        params:
            sample = "{sample}",
            tool = "diamond",
            mode = "contigs",
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_krona_by_pass_taxids.py"

    rule create_filtered_krona_report_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.filtered.krona_input.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/reports/{sample}.diamond.supported.filtered.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_diamond_contigs/{sample}.filtered.log"
        benchmark:
            config["output"] + "logs/krona_diamond_contigs/{sample}.filtered.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            r"""
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty filtered krona input (diamond contigs)." >> {log}
                touch {output}
            fi
            """
