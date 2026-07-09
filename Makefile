.PHONY: test test-dryrun test-empirical lint format run-meta run-consensus install install-dev build build-docker run-docker lock typecheck

# Conda env used for the end-to-end suite. viralunity pins python <3.12, so the
# empirical target runs inside this env rather than whatever interpreter invoked
# make (a python 3.12 base env would fail `pip install -e`). Override as needed:
#   make test-empirical VIRALUNITY_ENV=my-env
VIRALUNITY_ENV ?= viralunity

test: install-dev
	python3 -m unittest discover ./test -p *test.py

test-dryrun: install-dev
	pytest test/viralunity_dryrun_test.py -v

# Opt-in end-to-end suite: downloads real data, runs full pipelines. Runs inside
# the $(VIRALUNITY_ENV) conda env (see the variable above) so it works regardless
# of the interpreter that invoked make. Prerequisite: `viralunity setup
# --pipelines all` once to build the per-rule envs. Override the data cache with
# VIRALUNITY_TEST_CACHE=/some/path.
test-empirical:
	@command -v conda >/dev/null 2>&1 || { \
		echo "ERROR: conda not found on PATH. Install conda/miniforge (see README)."; exit 1; }
	@conda run -n $(VIRALUNITY_ENV) python --version >/dev/null 2>&1 || { \
		echo "ERROR: conda env '$(VIRALUNITY_ENV)' not found. Create it with:"; \
		echo "    conda env create -n $(VIRALUNITY_ENV) -f environment.yml"; \
		echo "or override the name: make test-empirical VIRALUNITY_ENV=<name>"; exit 1; }
	conda run --no-capture-output -n $(VIRALUNITY_ENV) python -m pip install -e ".[dev]"
	conda run --no-capture-output -n $(VIRALUNITY_ENV) pytest test/viralunity_empirical_test.py -v -m empirical

lint: install-dev
	black --check viralunity/ test/
	ruff check viralunity/ test/

# Type check (gating in CI; the backlog has been cleared to 0 errors).
typecheck: install-dev
	mypy viralunity/

format: install-dev
	black viralunity/ test/
	ruff check --fix viralunity/ test/

# Regenerate the pinned conda lockfile from environment.yml. Requires network
# access and conda-lock (pip install conda-lock). Commit conda-lock.yml so the
# runtime environment is reproducible across machines and over time.
lock:
	conda-lock lock --file environment.yml --lockfile conda-lock.yml

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

run-meta:
	viralunity meta \
		--data-type illumina \
		--sample-sheet input/viralunity_samplesheet.csv \
		--config-file output/config_meta.yml \
		--run-name test-meta \
		--kraken2-database input/database/kraken2 \
		--krona-database input/database/krona/taxonomy/ \
		--adapters input/references/SARS-CoV-2_RefSeq.fasta \
		--threads 6 \
		--threads-total 6 \
		--output output/test-meta-exmaple

run-consensus:
	viralunity consensus \
		--data-type illumina \
		--sample-sheet input/viralunity_samplesheet.csv \
		--config-file output/config_consensus.yml \
		--run-name test-consensus \
		--reference input/references/SARS-CoV-2_RefSeq.fasta \
		--primer-scheme input/references/scheme.bed \
		--threads 1 \
		--threads-total 1 \
		--output output/test-consensus-example

# Build the sdist + wheel for PyPI (see RELEASING.md). Requires `build` + `twine`.
build:
	python -m build
	twine check dist/*

build-docker:
	docker build -t viralunity/viralunity:latest .

run-docker:
	docker run --rm -i -t viralunity/viralunity:latest