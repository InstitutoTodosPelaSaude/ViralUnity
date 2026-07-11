#!/usr/bin/env python3
"""Extract the viral subset of de novo contigs from a classifier's contig→taxid map.

Used by the kraken2_contigs track to build the FASTA of contigs classified as
viral, so reads can be remapped to just those contigs for depth measurement
(mirroring how the diamond_contigs track remaps to its diamond-viral contigs).

A contig is kept when its assigned taxid's lineage includes the target root
(default ``10239`` = Viruses). This keeps the remap light — only viral contigs
are mapped, not the whole (often bacterial-dominated) assembly.
"""

import argparse
import os
import sys
from typing import Dict, Set

try:
    from viralunity.scripts.python.taxonomy import get_lineage, load_taxdump
except ImportError:
    from taxonomy import get_lineage, load_taxdump

VIRUSES_TAXID = "10239"


def viral_contig_ids(
    krona_path: str, parent_map: Dict[str, str], root: str = VIRUSES_TAXID
) -> Set[str]:
    """Return the set of contig ids whose taxid lineage includes *root*."""
    viral: Set[str] = set()
    if not krona_path or not os.path.exists(krona_path) or os.path.getsize(krona_path) == 0:
        return viral
    with open(krona_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            contig, taxid = cols[0], cols[1]
            if not taxid or taxid == "0" or taxid not in parent_map:
                continue
            if root in get_lineage(taxid, parent_map):
                viral.add(contig)
    return viral


def _iter_fasta(path: str):
    """Yield (header_id, full_record_text) for each record; id is the first token."""
    header = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header = line[1:].split()[0]
                buf = [line]
            elif header is not None:
                buf.append(line)
    if header is not None:
        yield header, "".join(buf)


def write_viral_contigs(
    contigs_fasta: str,
    krona_path: str,
    taxdump_dir: str,
    out_fasta: str,
    out_ids: str,
    root: str = VIRUSES_TAXID,
) -> int:
    """Write the viral-contig FASTA and its id list. Returns the number kept."""
    parent_map, _, _ = load_taxdump(os.path.join(taxdump_dir, "nodes.dmp"))
    viral = viral_contig_ids(krona_path, parent_map, root)

    os.makedirs(os.path.dirname(out_fasta) or ".", exist_ok=True)
    kept = 0
    with open(out_ids, "w") as ids_out, open(out_fasta, "w") as fa_out:
        if contigs_fasta and os.path.exists(contigs_fasta) and os.path.getsize(contigs_fasta) > 0:
            for header, record in _iter_fasta(contigs_fasta):
                if header in viral:
                    ids_out.write(header + "\n")
                    fa_out.write(record)
                    kept += 1
    return kept


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract viral contigs by classifier taxid.")
    p.add_argument("--contigs", required=True)
    p.add_argument("--krona", required=True, help="contig<TAB>taxid TSV")
    p.add_argument("--taxdump-dir", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-ids", required=True)
    p.add_argument("--root", default=VIRUSES_TAXID)
    return p.parse_args()


def main() -> None:
    if "snakemake" in globals():
        kept = write_viral_contigs(
            contigs_fasta=str(snakemake.input.contigs),  # noqa: F821
            krona_path=str(snakemake.input.krona),  # noqa: F821
            taxdump_dir=str(snakemake.params.taxdump),  # noqa: F821
            out_fasta=str(snakemake.output.fasta),  # noqa: F821
            out_ids=str(snakemake.output.ids),  # noqa: F821
        )
    else:
        args = parse_args()
        kept = write_viral_contigs(
            contigs_fasta=args.contigs,
            krona_path=args.krona,
            taxdump_dir=args.taxdump_dir,
            out_fasta=args.out_fasta,
            out_ids=args.out_ids,
            root=args.root,
        )
    sys.stderr.write(f"[extract_viral_contigs] kept {kept} viral contigs\n")


if __name__ == "__main__":
    main()
