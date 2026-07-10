# ViralUnity

ViralUnity is a tool for analysing viral high-throughput sequencing data. It is a Python package that orchestrates Snakemake workflows for data quality control, taxonomic assignment, and reference genome assembly. ViralUnity runs on *nix systems and can process entire sequencing runs in minimal time on a regular computer.

> **Full documentation:** <https://viralunity.readthedocs.io/en/latest/>

## Installation

Install the ViralUnity CLI from PyPI:

```bash
pip install viralunity
```

> **conda/mamba is still required at runtime.** ViralUnity orchestrates Snakemake, which
> builds the per-rule tool environments (aligners, classifiers, assemblers) via
> `--use-conda` on first run. Make sure conda or mamba is installed and on your `PATH`;
> you can pre-build those environments up front with `viralunity setup --pipelines all`.

To install from source for development instead:

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity
conda env create -n viralunity -f environment.yml
conda activate viralunity
pip install -e .
```

Per-rule conda environments under `viralunity/scripts/envs/` are managed automatically by Snakemake; the top-level `environment.yml` only installs ViralUnity itself and its core runtime dependencies.

## Quick start

Six top-level subcommands are exposed via the `viralunity` CLI:

```bash
viralunity create-samplesheet --input <runs-dir> --output samples.csv
viralunity get-databases all --path databases/
viralunity setup --pipelines all                # pre-build per-rule conda envs
viralunity consensus illumina --sample-sheet samples.csv --reference ref.fasta --output run/
viralunity meta      illumina --sample-sheet samples.csv --kraken2-database <db> --output run/
viralunity build-deacon-index --input host.fasta --output host.dcn
```

`get-databases all` grabs the four common databases; large optional ones are separate
subcommands (`virus-genome`, `deacon-index`, and `nr` for `meta --run-nr-validation`).

Global options: `--log-level {DEBUG,INFO,WARNING,ERROR}` and `--json-logs`
(e.g. `viralunity --log-level DEBUG meta ...`).

Each subcommand has its own `--help`; the same information is exhaustively documented in the `docs/` Sphinx site (rendered on ReadTheDocs at the link above).

## Tests

```bash
make test
```

This installs the package in editable mode (if not already installed) and runs the `unittest` suite under `test/`. Snakemake dry-run tests live in `test/viralunity_dryrun_test.py` and use `pytest`.

## Citation

A scientific publication describing this pipeline is being prepared. Meanwhile, please cite this repository. Primary references for upstream tools (fastp, MultiQC, Minimap2, Samtools, BCFtools, BEDtools, gofasta, MEGAHIT, Racon, BLAST, Kraken2, Krona, DIAMOND, Clair3, Medaka, Deacon) are listed in the ReadTheDocs site.

## License

MIT — see `LICENSE`.
