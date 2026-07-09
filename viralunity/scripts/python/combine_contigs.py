#!/usr/bin/env python3
"""Combine per-sample contig FASTAs into a single FASTA for one aggregated
similarity search.

Each record header is rewritten as ``>{sample}|{contig_id}`` so the origin
sample can be recovered when the search output is split back per sample (see
``split_search_output.py``). ``contig_id`` is the first whitespace-delimited
token of the original header, matching how megahit / downstream rules refer to
contigs. Missing or empty inputs are skipped so the rule stays robust when a
sample produced no contigs.
"""

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine per-sample contig FASTAs with sample-prefixed headers."
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        metavar="SAMPLE=PATH",
        help="Repeatable. A sample label and its FASTA path, e.g. sample-A=a.fa",
    )
    parser.add_argument("--output", required=True, help="Combined FASTA output path")
    return parser.parse_args()


def combine_contigs(inputs: Iterable[Tuple[str, str]], output_path) -> int:
    """Write a combined FASTA with ``>{sample}|{contig_id}`` headers.

    ``inputs`` is an iterable of ``(sample_label, fasta_path)``. Missing or
    empty files are skipped. Returns the number of records written.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with output_path.open("w") as out:
        for sample, fasta_path in inputs:
            path = Path(fasta_path)
            if not path.exists() or path.stat().st_size == 0:
                continue
            with path.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        contig_id = line[1:].strip().split()[0]
                        out.write(f">{sample}|{contig_id}\n")
                        written += 1
                    else:
                        out.write(line)
    return written


def _parse_input_pairs(items: List[str]) -> List[Tuple[str, str]]:
    pairs = []
    for item in items:
        sample, sep, path = item.partition("=")
        if not sep:
            raise ValueError(f"--input must be SAMPLE=PATH, got: {item}")
        pairs.append((sample, path))
    return pairs


def main() -> None:
    if "snakemake" in globals():
        samples = list(snakemake.params.samples)  # noqa: F821
        fastas = list(snakemake.input.fastas)  # noqa: F821
        inputs = list(zip(samples, fastas))
        output_path = snakemake.output[0]  # noqa: F821
    else:
        args = parse_args()
        inputs = _parse_input_pairs(args.input)
        output_path = args.output

    n = combine_contigs(inputs, output_path)
    print(f"[combine_contigs] wrote {n} records -> {output_path}")


if __name__ == "__main__":
    main()
