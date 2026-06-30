# Optional per-taxon genome-length table, required for RPKM computation.
#
# Rules in this file are only defined when compute_rpkm is True (i.e. the
# user provided --viral-genomes, which implies viral_genomes != "NA"). When
# compute_rpkm is False, the bleed-filter rules in each track file fall back
# to the *_RPM.tsv input directly and no rules here are activated.

if compute_rpkm:

    rule index_viral_genomes_for_rpkm:
        """Index the viral reference FASTA so per-accession lengths are available."""
        input:
            config["viral_genomes"]
        output:
            config["viral_genomes"] + ".fai"
        conda:
            "../envs/utils.yaml"
        shell:
            "samtools faidx {input}"

    rule build_genome_length_table:
        """
        Build a median genome-length table per (rank, taxid) from the RefSeq
        viral genome FASTA index and genome2taxid mapping.  The output feeds
        each per-track RPKM computation step.
        """
        input:
            fai          = config["viral_genomes"] + ".fai",
            genome2taxid = config["viral_taxids"],
        output:
            config["output"] + "metagenomics/genome_lengths.tsv"
        params:
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/build_genome_length_table.py"
