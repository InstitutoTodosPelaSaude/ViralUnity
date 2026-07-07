# Architecture

ViralUnity is a thin Python orchestration layer over Snakemake workflows. The
Python code validates inputs, writes a YAML config, and launches one of six
Snakemake workflows; all the bioinformatics lives in the `.smk` rule files.

## Data flow

```
CLI (click)                       viralunity_cli.py  →  <subcommand>_cli.py
   │  parses options, configures logging (--log-level / --json-logs)
   ▼
main(args)                        viralunity_consensus.py / viralunity_meta.py
   │  thin wrapper around…
   ▼
_orchestrator.run_pipeline(...)   viralunity/_orchestrator.py
   │  resolve_paths → validate → generate_config → [run manifest] → run_workflow
   ▼
ConfigGenerator → YAML            viralunity/config_generator.py
   │  emits the exact keys the .smk files read (config["samples"], …)
   ▼
snakemake(workflow.smk, config)   viralunity/scripts/<pipeline>_<datatype>.smk
   │  rule all → _all_inputs() → include: rules/*.smk
   ▼
per-rule conda envs               viralunity/scripts/envs/*.yaml  (--use-conda)
```

## Module map

| Module | Responsibility |
|---|---|
| `viralunity_cli.py` | Top-level click group; wires `--log-level` / `--json-logs`. |
| `<sub>_cli.py` | Per-subcommand click options → plain `args` dict. |
| `viralunity_consensus.py`, `viralunity_meta.py` | Own `validate_args`, `generate_config_file`, `run_snakemake_workflow`; call the orchestrator. |
| `_orchestrator.py` | Shared `run_pipeline` skeleton (resolve → validate → config → manifest → run) with structured error handling. |
| `validators.py` | File/dir existence, sample-sheet parsing, cross-dependency checks, input sanitization. |
| `config_generator.py` | Writes the YAML config (the contract with the `.smk` files). |
| `constants.py` | `ConfigKeys`, `DataType`, `ResourceDefaults` (per-pipeline rule lists). |
| `exceptions.py` | Typed error hierarchy with machine-readable `code`s. |
| `logging_config.py` | Central logging (run id, text/JSON). |
| `provenance.py` | `run_manifest.json` (version, config, input checksums). |
| `_subprocess.py` | `run_command` for external tools (get-databases / build-deacon-index), with timeouts. |
| `scripts/*.smk`, `scripts/rules/*.smk` | The actual pipelines. |
| `scripts/python/*.py` | Helpers run via Snakemake's `script:` directive. |

## The config is the contract

`ConfigGenerator` emits keys that are hard-coded as strings in the `.smk` files
(e.g. `config["samples"]`, `config["output"]`, `config["run_denovo_assembly"]`).
Adding a pipeline option means touching four places: the click option in the
relevant `*_cli.py`, the `validators.py` check, a `ConfigGenerator.add_*` setter,
and the rule(s) that read it. `config["samples"]` maps `sample-<id>` to a **list**
of FASTQ paths (readers tolerate the legacy space-joined string form too).

## Workflow selection

`run_snakemake_workflow` formats a path:
`scripts/{consensus,metagenomics}_{illumina,nanopore}[_segmented].smk`.
Segmentation is chosen when `args["reference"]` is a dict (from repeatable
`--segmented-reference SEGMENT=PATH`). There is no segmented metagenomics variant.
