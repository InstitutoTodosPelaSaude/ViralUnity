"""Rules shared between the segmented and non-segmented Nanopore consensus
workflows.

The including snakefile must define:
* ``SEGMENT_WILDCARD``       -- ``"{segment}/"`` (segmented) or ``""`` (single)
* ``REFERENCE``              -- path string / callable taking ``wildcards`` and
                                returning a path (segmented) -- forwarded to
                                ``params.reference``
"""


rule calculate_assembly_statistics:
    conda:
        "../envs/utils.yaml"
    input:
        get_map_input_fastqs,
        get_map_input_fastqs,
        get_map_input_fastqs,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "coverage_stats/{sample}.stats_summary.csv")
    params:
        minimum_depth = config["minimum_depth"]
    script:
        "../python/calculate_assembly_stats.py"


rule align_consensus_to_reference_genome:
    conda:
        "../envs/alignment.yaml"
    input:
        # Literal path (instead of ``rules.unify_assembly_statistics_reports.``)
        # to avoid an include-order cycle: ``unify_assembly_statistics_reports``
        # lives in the top-level snakefile and itself references
        # ``rules.calculate_assembly_statistics`` defined here.
        stats = config['output'] + "assembly/assembly_stats_summary.csv",
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"],
            allow_missing=True
        )
    output:
        aln_consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/samples_alignment.fasta",
        combined = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/consensus.fasta",
        sam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/aln.consensus.sam",
        masked = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/aln.consensus.indelsMasked.fasta"
    params:
        reference = REFERENCE,
        minimap2_flags = config.get("minimap2_consensus_align_flags", "-a --sam-hit-only --secondary=no --score-N=0")
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/align_consensus_to_reference_genome.log"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        cat {params.reference} {input.consensus_files} > {output.combined}
        minimap2 {params.minimap2_flags} {params.reference} {output.combined} -o {output.sam}
        gofasta sam toMultiAlign --pad -s {output.sam} -o {output.aln_consensus}
        sed '/^>/ ! s/-/N/g' {output.aln_consensus} > {output.masked}
        """
