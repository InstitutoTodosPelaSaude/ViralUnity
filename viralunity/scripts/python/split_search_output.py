#!/usr/bin/env python3
"""Split an aggregated outfmt-6 similarity-search TSV back into per-sample files.

The combined query used ``{sample}|{contig_id}`` headers (see
``combine_contigs.py``); this reverses that: each row is routed to its sample's
output file with the ``{sample}|`` prefix stripped from the qseqid (column 1) so
downstream rules see the original contig IDs and behave exactly as they did for
per-sample searches. Every sample gets an output file (empty if it had no hits),
so Snakemake targets always resolve. Rows whose prefix is not a known sample are
skipped.
"""

import argparse
from pathlib import Path
from typing import Dict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split an aggregated outfmt-6 TSV into per-sample TSVs."
    )
    parser.add_argument("--combined", required=True, help="Aggregated outfmt-6 TSV")
    parser.add_argument(
        "--output",
        action="append",
        required=True,
        metavar="SAMPLE=PATH",
        help="Repeatable. A sample label and its per-sample output path.",
    )
    parser.add_argument(
        "--delimiter", default="|", help="Sample/contig header delimiter (default: |)"
    )
    return parser.parse_args()


def split_search_output(
    combined_tsv, sample_to_output: Dict[str, str], delimiter: str = "|"
) -> Dict[str, int]:
    """Route rows of ``combined_tsv`` to per-sample outputs, stripping the prefix.

    Returns ``{sample: n_rows}``. Every key in ``sample_to_output`` gets a file,
    empty if no rows matched. A missing/empty combined input produces all-empty
    outputs.
    """
    counts: Dict[str, int] = {sample: 0 for sample in sample_to_output}
    handles = {}
    try:
        for sample, out_path in sample_to_output.items():
            path = Path(out_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            handles[sample] = path.open("w")

        combined_path = Path(combined_tsv)
        if combined_path.exists() and combined_path.stat().st_size > 0:
            with combined_path.open() as fh:
                for raw in fh:
                    line = raw.rstrip("\n")
                    if not line.strip():
                        continue
                    qseqid, sep, rest = line.partition("\t")
                    sample, dsep, contig = qseqid.partition(delimiter)
                    if not dsep or sample not in handles:
                        continue
                    restored = contig + ("\t" + rest if sep else "")
                    handles[sample].write(restored + "\n")
                    counts[sample] += 1
    finally:
        for fh in handles.values():
            fh.close()
    return counts


def _parse_output_pairs(items):
    mapping = {}
    for item in items:
        sample, sep, path = item.partition("=")
        if not sep:
            raise ValueError(f"--output must be SAMPLE=PATH, got: {item}")
        mapping[sample] = path
    return mapping


def main() -> None:
    if "snakemake" in globals():
        samples = list(snakemake.params.samples)  # noqa: F821
        outputs = list(snakemake.output)  # noqa: F821
        sample_to_output = dict(zip(samples, outputs))
        combined = snakemake.input[0]  # noqa: F821
        delimiter = getattr(snakemake.params, "delimiter", "|")  # noqa: F821
    else:
        args = parse_args()
        sample_to_output = _parse_output_pairs(args.output)
        combined = args.combined
        delimiter = args.delimiter

    counts = split_search_output(combined, sample_to_output, delimiter)
    print(f"[split_search_output] rows per sample: {counts}")


if __name__ == "__main__":
    main()
