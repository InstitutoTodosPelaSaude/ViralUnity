#!/usr/bin/env python

"""Standalone ``viralunity report`` subcommand.

A one-shot operation (not a pipeline launcher): it reads an existing consensus
output directory and writes a single self-contained interactive HTML report.
The heavy lifting lives in the shared core module
``viralunity.scripts.python.generate_consensus_report``, which is the same code
the Snakemake ``generate_html_report`` rule invokes — so the CLI and the
pipeline produce identical reports for the same output directory.
"""

import logging
import os

import click

from viralunity.scripts.python.generate_consensus_report import (
    build_report_metadata,
    load_run_config,
    write_report,
)

logger = logging.getLogger(__name__)


@click.command("report")
@click.option(
    "--input",
    "input_dir",
    required=True,
    help="Existing consensus output directory (contains assembly/assembly_stats_summary.csv).",
)
@click.option(
    "--output",
    "output_path",
    default=None,
    help="Destination HTML path (default: <input>/report.html).",
)
@click.option(
    "--config-file",
    "config_file",
    default=None,
    help="Run config YAML used to launch the pipeline; supplies platform, "
    "library-layout, and primer-scheme metadata. Inferred from the output "
    "directory when omitted.",
)
def report(input_dir, output_path, config_file):
    """Generate an interactive HTML report from a consensus output directory."""
    if not os.path.isdir(input_dir):
        raise click.ClickException(f"Input directory not found: {input_dir}")
    if config_file and not os.path.isfile(config_file):
        raise click.ClickException(f"Config file not found: {config_file}")
    dest = output_path or os.path.join(input_dir, "report.html")
    try:
        metadata = build_report_metadata(load_run_config(config_file), input_dir)
        write_report(input_dir, dest, metadata)
        logger.info(f"Consensus report generated: {dest}")
    except Exception as e:
        raise click.ClickException(str(e))
