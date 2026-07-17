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
@click.option(
    "--pass-threshold",
    "pass_threshold",
    type=float,
    default=None,
    help="Coverage fraction for the green/pass tier (default 0.90); drives the "
    "'>=90% coverage' KPI and the status dots/bars.",
)
@click.option(
    "--warn-threshold",
    "warn_threshold",
    type=float,
    default=None,
    help="Coverage fraction for the amber/warn tier (default 0.70); drives the "
    "'Below 70%' KPI and the 'low coverage only' filter.",
)
@click.option(
    "--chart-color",
    "chart_color",
    default=None,
    help="Accent #RRGGBB hex colour for the coverage-heatmap colour scale.",
)
@click.option(
    "--colorbar-thickness",
    "colorbar_thickness",
    type=int,
    default=None,
    help="Heatmap colour-bar thickness in px (default 14).",
)
@click.option(
    "--fetch-annotation/--no-fetch-annotation",
    "fetch_annotation",
    default=True,
    help="When a run has no gene annotation staged or in its config, fetch it from "
    "NCBI by the reference accession (default on); --no-fetch-annotation keeps the "
    "report fully offline.",
)
def report(
    input_dir,
    output_path,
    config_file,
    pass_threshold,
    warn_threshold,
    chart_color,
    colorbar_thickness,
    fetch_annotation,
):
    """Generate an interactive HTML report from a consensus output directory."""
    # Imported lazily so the report generator's dependencies (plotly, jinja2)
    # are only required when a report is actually built. A module-scope import
    # here would make them hard requirements for every ``viralunity`` subcommand,
    # even though the pipeline builds reports in a separate ``envs/report.yaml``
    # conda env. See test/cli_import_constraint_test.py.
    try:
        from viralunity.scripts.python.generate_consensus_report import (
            build_report_metadata,
            load_run_config,
            report_params_from_config,
            write_report,
        )
    except ImportError as e:
        raise click.ClickException(
            f"The report generator's dependencies are missing ({e}). "
            "Install them with 'pip install plotly jinja2' (or reinstall "
            "viralunity, which declares them)."
        )
    if not os.path.isdir(input_dir):
        raise click.ClickException(f"Input directory not found: {input_dir}")
    if config_file and not os.path.isfile(config_file):
        raise click.ClickException(f"Config file not found: {config_file}")
    dest = output_path or os.path.join(input_dir, "report.html")
    try:
        config = load_run_config(config_file)
        metadata = build_report_metadata(config, input_dir)
        params = report_params_from_config(
            config,
            pass_threshold=pass_threshold,
            warn_threshold=warn_threshold,
            chart_color=chart_color,
            colorbar_thickness=colorbar_thickness,
        )
        write_report(input_dir, dest, metadata, config, params, fetch_annotation)
        logger.info(f"Consensus report generated: {dest}")
    except Exception as e:
        raise click.ClickException(str(e))
