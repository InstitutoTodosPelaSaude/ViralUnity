"""Run provenance manifest.

For public-health reporting, results must be reproducible-by-record: given an
output, you must be able to recover *which* pipeline version, config, and exact
input files produced it. This module writes a ``run_manifest.json`` into the run
output directory capturing the ViralUnity version, a timestamp, the resolved
config path, and a checksum/size for every input FASTQ.

Tool and database versions are best captured at rule-execution time (inside the
per-rule conda envs) and are intentionally out of scope here; this manifest
covers the orchestration-level provenance that the Python layer can record
reliably.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
from typing import Any, Dict, List, Optional

from viralunity import __version__

MANIFEST_FILENAME = "run_manifest.json"


def _sha256(path: str, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def _describe_input(path: str) -> Dict[str, Any]:
    """Return a provenance record for a single input file."""
    record: Dict[str, Any] = {"path": os.path.abspath(path)}
    if os.path.isfile(path):
        stat = os.stat(path)
        record["size_bytes"] = stat.st_size
        record["sha256"] = _sha256(path)
    else:
        record["missing"] = True
    return record


def build_run_manifest(
    args: Dict[str, Any],
    samples: Dict[str, List[str]],
    *,
    timestamp: Optional[str] = None,
) -> Dict[str, Any]:
    """Build (but do not write) the run-manifest dict.

    Args:
        args: The pipeline argument dict (needs ``output``, ``run_name``,
            ``config_file``, ``data_type``).
        samples: Mapping of sample id -> list of input FASTQ paths.
        timestamp: ISO-8601 timestamp; generated (UTC) if omitted.
    """
    if timestamp is None:
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    sample_inputs = {
        sample: [_describe_input(p) for p in paths] for sample, paths in (samples or {}).items()
    }

    return {
        "viralunity_version": __version__,
        "created_utc": timestamp,
        "run_name": args.get("run_name"),
        "data_type": args.get("data_type"),
        "config_file": (os.path.abspath(args["config_file"]) if args.get("config_file") else None),
        "output": (
            os.path.abspath(os.path.join(args["output"], args.get("run_name", "")))
            if args.get("output")
            else None
        ),
        "sample_count": len(samples or {}),
        "samples": sample_inputs,
    }


def write_run_manifest(
    args: Dict[str, Any],
    samples: Dict[str, List[str]],
    *,
    timestamp: Optional[str] = None,
) -> str:
    """Write the run manifest into ``<output>/<run_name>/run_manifest.json``.

    Returns:
        The path to the written manifest.
    """
    manifest = build_run_manifest(args, samples, timestamp=timestamp)
    run_dir = os.path.join(args["output"], args.get("run_name", ""))
    os.makedirs(run_dir, exist_ok=True)
    manifest_path = os.path.join(run_dir, MANIFEST_FILENAME)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return manifest_path
