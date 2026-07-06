.PHONY: test test-dryrun test-empirical lint format run-meta run-consensus install install-dev build-docker run-docker lock typecheck

test: install-dev
	python3 -m unittest discover ./test -p *test.py

test-dryrun: install-dev
	pytest test/viralunity_dryrun_test.py -v

# Opt-in end-to-end suite: downloads real data, runs full pipelines.
# Prerequisite: viralunity on PATH + `viralunity setup --pipelines all` once.
# Override cache location with VIRALUNITY_TEST_CACHE=/some/path.
test-empirical: install-dev
	pytest test/viralunity_empirical_test.py -v -m empirical

lint: install-dev
	black --check viralunity/ test/
	ruff check viralunity/ test/

# Advisory type check (not yet gating; ~14 known findings to clear over time).
typecheck: install-dev
	mypy viralunity/ || true

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
		--config-file output/config_meta.yml \
		--run-name test-meta \
		--kraken2-database input/database/kraken2 \
		--krona-database input/database/krona/taxonomy/ \
		--adapters input/references/SARS-CoV-2_RefSeq.fasta \
		--threads 1 \
		--threads-total 1 \
		--output output/test-meta-exmaple

build-docker:
	docker build -t viralunity/viralunity:latest .

run-docker:
	docker run --rm -i -t viralunity/viralunity:latest

conda-build:
	conda build viralunity/meta.yaml