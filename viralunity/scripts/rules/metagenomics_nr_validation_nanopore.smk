# NR validation (Nanopore): re-search the de novo viral contigs against NCBI nr
# ONCE (all samples combined, sample-prefixed headers), resolve taxonomy from
# DIAMOND staxids, and reduce to one representative hit + LCA consensus per
# contig. Produces the per-(sample,contig) NR verdict table consumed by the
# harmonize step. Contig tracks only; requires denovo + diamond_contigs.
if run_denovo and run_diamond_contigs and run_nr_validation:

    rule build_nr_query:
        input:
            fastas = expand(config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa", sample=list(config["samples"]))
        output:
            config["output"] + "metagenomics/nr_validation/nr_query.fasta"
        params:
            samples = list(config["samples"])
        conda:
            "../envs/utils.yaml"
        script:
            "../python/combine_contigs.py"

    rule run_diamond_nr:
        input:
            query = config["output"] + "metagenomics/nr_validation/nr_query.fasta"
        output:
            tsv = config["output"] + "metagenomics/nr_validation/nr_query.nr.diamond.tsv"
        params:
            db = config.get("nr_diamond_database", "NA"),
            evalue = config.get("nr_evalue", 1e-10),
            max_target_seqs = config.get("nr_max_target_seqs", 10),
            sensitivity = config.get("nr_sensitivity", "fast")
        threads: config.get("run_diamond_nr_cpus", 4)
        resources:
            mem_mb = config.get("run_diamond_nr_ram", 16) * 1024
        log:
            config["output"] + "logs/nr_validation/diamond_nr.log"
        benchmark:
            config["output"] + "logs/nr_validation/diamond_nr.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            if [ ! -s {input.query} ]; then
                echo "Empty NR query; creating empty output." >> {log}
                touch {output.tsv}
            else
                if [[ "{params.db}" != *.dmnd ]]; then
                    diamond prepdb --db {params.db} >> {log} 2>&1 || true
                fi
                diamond blastx --db {params.db} --query {input.query} \
                    --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
                    --max-target-seqs {params.max_target_seqs} --evalue {params.evalue} \
                    --{params.sensitivity} --threads {threads} 2>> {log}
            fi
            """

    rule annotate_nr_taxonomy:
        input:
            config["output"] + "metagenomics/nr_validation/nr_query.nr.diamond.tsv"
        output:
            config["output"] + "metagenomics/nr_validation/nr_query.nr.taxonomy.tsv"
        params:
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/annotate_nr_taxonomy.py"

    rule filter_nr_top_species:
        input:
            config["output"] + "metagenomics/nr_validation/nr_query.nr.taxonomy.tsv"
        output:
            table = config["output"] + "metagenomics/nr_validation/nr_query.nr.top_species_hit_lca.tsv",
            viral = config["output"] + "metagenomics/nr_validation/nr_query.nr.top_species_hit_lca.viruses_only.tsv"
        params:
            consensus_threshold = config.get("nr_consensus_threshold", 0.5)
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_top_species_hits.py"
