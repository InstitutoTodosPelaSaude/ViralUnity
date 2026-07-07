# Embedding ViralUnity in a service

ViralUnity is designed to be driven by a web service as its bioinformatics core.
This page documents the integration seam: how to call it in-process, the input
contract you must enforce, and what it emits for observability and provenance.
It stops short of prescribing an HTTP API or job queue — those are the caller's
responsibility.

## Calling the pipeline

Each pipeline's `main(args)` is the entry point; both delegate to
`viralunity._orchestrator.run_pipeline`, which runs
`resolve_paths → validate → generate_config → write_run_manifest → run_workflow`
and returns an exit code (`0` success, `1` failure).

```python
from viralunity.logging_config import configure_logging
from viralunity.viralunity_meta import main as meta_main

configure_logging(level="INFO", json_logs=True, run_id=job_id)
exit_code = meta_main({
    "data_type": "illumina",
    "sample_sheet": "/jobs/<job>/samples.csv",
    "config_file": "/jobs/<job>/config.yaml",
    "output": "/jobs/<job>/out",
    "run_name": job_id,          # must be a safe identifier — see below
    "kraken2_database": "...",
    "threads": 8, "threads_total": 8,
})
```

## Per-job isolation (required)

Snakemake writes a `.snakemake/` lock/state directory into the working directory
and per-run outputs under `<output>/<run_name>/`. **Run each job in its own
working directory** (a fresh temp dir per job) so concurrent jobs cannot collide
on locks, logs, or conda-env creation. The in-process `snakemake()` call is not
safe to run concurrently from the same cwd.

## Input contract (enforce on untrusted uploads)

- **Sample ids** and **`run_name`** become shell tokens and filesystem path
  components. The sample-sheet parser rejects unsafe ids; for service-supplied
  values use `validators.sanitize_identifier(value, field=...)`.
- Constrain user-supplied output/config paths to a per-job base with
  `validators.ensure_within_base(path, base)` to block path traversal.
- Sample sheets are validated up-front: duplicate ids, ragged rows, wrong column
  counts, and missing files all raise `SampleSheetError` / `ViralUnityFileNotFoundError`.

## Observability

Call `configure_logging(level, json_logs=True, run_id=<job id>)` once per process
(or per job in a worker). Every log record is stamped with the run id; JSON mode
emits one object per line for ingestion.

## Structured errors

Expected failures raise subclasses of `viralunity.exceptions.ViralUnityError`,
each with a machine-readable `code` and a `to_dict()` payload:

```python
from viralunity.exceptions import ViralUnityError
try:
    ...
except ViralUnityError as e:
    return {"status": "error", **e.to_dict()}   # {"error","code","message"}
```

`run_pipeline` already logs these as `[<code>] <message>` and returns `1`;
catch them earlier if you need the structured payload.

## Provenance

Each run writes `<output>/<run_name>/run_manifest.json` with the ViralUnity
version, a UTC timestamp, the resolved config path, and a `sha256`+size for every
input FASTQ — enough to reproduce a result by record. Persist it alongside job
outputs. Per-rule tool/DB versions are not captured here (they live in the
per-rule conda envs); pin those via `make lock` for full reproducibility.
