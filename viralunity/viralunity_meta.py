#!/usr/bin/env python

"""
This scripts dynamically creates config files and runs
the viralunity metagenomics snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import logging
import os
from pathlib import Path
from typing import Any, Dict

from viralunity import _orchestrator
from viralunity.constants import DataType, ResourceDefaults
from viralunity.exceptions import ValidationError
from viralunity.validators import (
    META_PATH_ARG_KEYS,
    get_samples_from_args,
    resolve_path_args,
    sanitize_identifier,
    validate_illumina_requirements,
    validate_metagenomics_requirements,
    validate_nanopore_requirements,
)

# Set up logging
logger = logging.getLogger(__name__)


def validate_args(args: Dict[str, Any]) -> Dict[str, list]:
    """Validate all pipeline arguments.

    Args:
        args: Dictionary of pipeline arguments

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If validation fails
    """
    logger.info("Validating pipeline arguments")

    # Reject unsafe run names before they become output-path / DAG components
    # (run_name is spliced into config["output"] and every rule's output path).
    if args.get("run_name"):
        sanitize_identifier(args["run_name"], "run_name")

    # Get and validate samples
    samples = get_samples_from_args(args)

    logger.info(f"Found {len(samples)} samples")

    # Validate metagenomics-specific requirements
    validate_metagenomics_requirements(args)

    # Validate data-type specific requirements
    data_type = args.get("data_type")
    if data_type == DataType.ILLUMINA:
        validate_illumina_requirements(args)
    elif data_type == DataType.NANOPORE:
        validate_nanopore_requirements(args)

    logger.info("All arguments validated successfully")
    return samples


def generate_config_file(samples: Dict[str, list], args: Dict[str, Any]) -> None:
    """Generate configuration file for the pipeline.

    Args:
        samples: Dictionary mapping sample names to file paths
        args: Dictionary of pipeline arguments

    Raises:
        ValidationError: If config generation fails
    """
    data_type = args["data_type"]
    generator = _orchestrator.start_config(args, samples)

    # Add metagenomics-specific settings
    negative = args.get("negative_controls")
    if isinstance(negative, str):
        negative = [x.strip() for x in negative.split(",") if x.strip()]
    elif negative is None:
        negative = []

    # Negative-control IDs are user-facing raw sample IDs, but samples are stored
    # in the config (and referenced throughout the .smk rules and the summary
    # 'sample' column) with a 'sample-' prefix applied by ConfigGenerator.add_samples.
    # Prefix the negatives to match, validating that each names a real sample so a
    # typo fails fast here instead of silently disabling the negative-control filter.
    known_samples = set(samples or {})
    prefixed_negatives = []
    for nc in negative:
        raw = nc[len("sample-") :] if nc.startswith("sample-") else nc
        if known_samples and raw not in known_samples:
            raise ValidationError(
                f"Negative control '{nc}' does not match any sample in the sample sheet. "
                f"Known samples: {sorted(known_samples)}"
            )
        prefixed_negatives.append(f"sample-{raw}")
    negative = prefixed_negatives

    # Normalize Krona database path: if it points to a dir containing 'taxonomy', append it.
    krona_db = args.get("krona_database", "NA")
    if krona_db and krona_db != "NA":
        kp = Path(krona_db)
        if (kp / "taxonomy").is_dir() and kp.name != "taxonomy":
            krona_db = str(kp / "taxonomy")

    generator.add_metagenomics_settings(
        kraken2_database=args["kraken2_database"],
        krona_database=krona_db,
        remove_human_reads=args.get("remove_human_reads", False),
        remove_unclassified_reads=args.get("remove_unclassified_reads", False),
        host_reference=args.get("host_reference", "NA"),
        deacon_index=args.get("deacon_index", "NA"),
        taxdump=args.get("taxdump", "NA"),
        run_denovo_assembly=args.get("run_denovo_assembly", False),
        run_kraken2_reads=args.get("run_kraken2_reads", True),
        run_kraken2_contigs=args.get("run_kraken2_contigs", True),
        run_diamond_reads=args.get("run_diamond_reads", False),
        run_diamond_contigs=args.get("run_diamond_contigs", False),
        taxids=args.get("taxids", "NA"),
        diamond_database=args.get("diamond_database", "NA"),
        diamond_sensitivity=args.get("diamond_sensitivity", "sensitive"),
        evalue=args.get("evalue", 0.001),
        bleed_fraction=args.get("bleed_fraction", 0.005),
        negative_controls=negative,
        minimum_hit_group=args.get("minimum_hit_group", 4),
        diamond_max_target_seqs=args.get("diamond_max_target_seqs", 1),
        kraken2_extra_flags=args.get("kraken2_extra_flags", "--report-minimizer-data"),
        combine_contig_search=args.get("combine_contig_search", False),
        run_ictv_host_filter=args.get("run_ictv_host_filter", False),
        ictv_vertebrate_taxids_file=args.get("ictv_vertebrate_taxids_file", "NA"),
        run_nr_validation=args.get("run_nr_validation", False),
        nr_diamond_database=args.get("nr_diamond_database", "NA"),
        nr_evalue=args.get("nr_evalue", 1e-10),
        nr_max_target_seqs=args.get("nr_max_target_seqs", 10),
        nr_sensitivity=args.get("nr_sensitivity", "fast"),
        nr_consensus_threshold=args.get("nr_consensus_threshold", 0.5),
        compute_rpkm=args.get("viral_genomes", "NA") not in ("NA", "", None),
        enrichment_pseudocount=args.get("enrichment_pseudocount", 1.0),
        z_score_threshold=args.get("z_score_threshold", 3.0),
        log2_ratio_threshold=args.get("log2_ratio_threshold", 1.0),
    )

    generator.add_reference_assembly_settings(
        run_reference_assembly=args.get("run_reference_assembly", False),
        method=args.get("method", "kraken2"),
        source=args.get("source", "reads"),
        reads_count=args.get("reads_count", 100),
        contigs_count=args.get("contigs_count", 1),
        families=args.get(
            "families",
            "Coronaviridae,Orthomyxoviridae,Flaviviridae,Herpesviridae,Papillomaviridae,Paramyxoviridae,Adenoviridae",
        ),
        reference_selection_strategy=args.get("reference_selection_strategy", "taxid"),
        blast_qcov=args.get("blast_qcov", 80),
        blast_pident=args.get("blast_pident", 80),
        viral_genomes=args.get("viral_genomes", "databases/virus_genomes/viral.genomes.fasta"),
        viral_taxids=args.get("viral_taxids", "databases/virus_genomes/genome2taxid.tsv"),
    )

    if data_type == DataType.ILLUMINA:
        generator.add_illumina_settings(
            adapters=args.get("adapters") or "NA",
            minimum_read_length=args.get("minimum_read_length", 50),
            trim_head=args.get("trim_head", 0),
            trim_tail=args.get("trim_tail", 0),
        )
    else:
        generator.add_nanopore_settings(
            run_polish_racon=args.get("run_polish_racon", False),
            run_polish_medaka=args.get("run_polish_medaka", False),
            medaka_model=args.get("medaka_model"),
            clair3_model=args.get("clair3_model"),
        )

    # Add resource settings
    shared_rules = ResourceDefaults.META_SHARED_RULES
    if data_type == DataType.ILLUMINA:
        generator.add_resource_settings(args, shared_rules + ResourceDefaults.META_ILLUMINA_RULES)
    else:
        generator.add_resource_settings(args, shared_rules + ResourceDefaults.META_NANOPORE_RULES)

    # Save config file
    generator.save()

    logger.info(f"Configuration file generated: {args['config_file']}")


def run_snakemake_workflow(args: Dict[str, Any]) -> bool:
    """Run the Snakemake workflow for the metagenomics pipeline."""
    thisdir = os.path.abspath(os.path.dirname(__file__))
    workflow_path = os.path.join(thisdir, "scripts", f"metagenomics_{args['data_type']}.smk")
    return _orchestrator.run_workflow(workflow_path, args)


def main(args: Dict[str, Any]) -> int:
    """Main entry point for the metagenomics pipeline.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    return _orchestrator.run_pipeline(
        args,
        resolve_paths=lambda a: resolve_path_args(a, META_PATH_ARG_KEYS),
        validate=validate_args,
        generate_config=generate_config_file,
        run_workflow_fn=run_snakemake_workflow,
        skip_when_no_samples=True,
    )
