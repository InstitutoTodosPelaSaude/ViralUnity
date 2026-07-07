"""Top-level click CLI for ViralUnity."""

import click

from viralunity import __program__, __version__
from viralunity.logging_config import configure_logging
from viralunity.viralunity_build_deacon_index_cli import build_deacon_index
from viralunity.viralunity_consensus_cli import consensus
from viralunity.viralunity_create_samplesheet import create_samplesheet
from viralunity.viralunity_get_databases_cli import get_databases
from viralunity.viralunity_meta_cli import meta
from viralunity.viralunity_setup_cli import setup


@click.group()
@click.version_option(version=__version__, prog_name=__program__)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    help="Logging verbosity.",
)
@click.option(
    "--json-logs",
    is_flag=True,
    default=False,
    help="Emit structured JSON logs (one object per line) instead of text.",
)
def cli(log_level: str, json_logs: bool) -> None:
    """ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data.

    \b
    Subcommands:
    * consensus           reference-guided consensus assembly (illumina/nanopore)
    * meta                metagenomic classification and de novo assembly
    * setup               pre-build per-rule conda envs into a shared cache
    * create-samplesheet  generate a sample sheet from a sequencing directory
    * get-databases       download/build reference databases
    * build-deacon-index  build a Deacon index for host depletion

    Run ``viralunity <subcommand> --help`` for subcommand-specific options.
    """
    configure_logging(level=log_level, json_logs=json_logs)


cli.add_command(consensus)
cli.add_command(meta)
cli.add_command(setup)
cli.add_command(create_samplesheet)
cli.add_command(get_databases)
cli.add_command(build_deacon_index)


def main() -> None:
    cli()


if __name__ == "__main__":
    main()
