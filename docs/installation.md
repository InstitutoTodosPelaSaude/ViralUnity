# Installation

ViralUnity is a Python package that launches Snakemake workflows. The tool dependencies are documented in the conda environment file `environment.yml` and the workflow-specific dependencies are in `viralunity/scripts/envs/`.

## Clone and create the environment

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity/
conda env create -f environment.yml
conda activate viralunity
```

Or with **micromamba** (recommended on macOS with Apple Silicon):

```bash
micromamba env create -f environment.yml --platform osx-64
micromamba activate viralunity
```

```{warning}
On macOS with Apple Silicon (M1 or later), the `viralunity/scripts/envs/clair3.yaml` environment may not install correctly due to compatibility constraints in the clair3 dependencies.
```

## Troubleshooting

### `CreateCondaEnvironmentException` with a `repodata_shards.msgpack.zst` 404

If the first pipeline run aborts on `Creating conda environment .../qc.yaml ...` with a chained `HTTPError: 404` for `repodata_shards.msgpack.zst` followed by a `JSONDecodeError`, you are hitting a known incompatibility between conda 26.x's libmamba shards path and bioconda (which does not publish shards). `environment.yml` already pins `conda<26` and `conda-libmamba-solver<26` to avoid this, so the fix is to **rebuild the `viralunity` env from the pinned file**:

```bash
conda deactivate
conda env remove -n viralunity
conda env create -n viralunity -f environment.yml
conda activate viralunity
pip install -e .
```

If you cannot rebuild the env (for example, on a shared cluster install), the manual escape hatch is to fall back to the classic solver:

```bash
conda config --set solver classic
```

Full background and the remediation plan are documented in [`CONDA_ENV_CREATION_FIX.md`](https://github.com/InstitutoTodosPelaSaude/ViralUnity/blob/main/CONDA_ENV_CREATION_FIX.md).

## Verify the installation

```bash
viralunity --version
viralunity --help
```

## First-time environment setup

After installing ViralUnity, build the per-rule conda environments once with:

```bash
viralunity setup --pipelines all
```

This downloads and resolves every sub-pipeline dependency into `~/.cache/viralunity/conda-envs/` (override with `--conda-prefix PATH`, or set `$VIRALUNITY_CONDA_PREFIX`). Future `viralunity consensus` / `viralunity meta` runs reuse the cached envs and do not need to re-create them per working directory, which both speeds up first runs and isolates env-creation failures from real pipeline runs.

Run `viralunity setup --pipelines consensus-illumina --dry-run` first to inspect what would be built.

## Development install

To work on ViralUnity itself, install the optional `dev` extras (linters and tests):

```bash
pip install -e ".[dev]"
```

See [CONTRIBUTING.md](https://github.com/InstitutoTodosPelaSaude/ViralUnity/blob/main/CONTRIBUTING.md) for the full development workflow.
