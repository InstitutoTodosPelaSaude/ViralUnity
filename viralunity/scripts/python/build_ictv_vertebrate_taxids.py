#!/usr/bin/env python3
"""Build the vertebrate-infecting-virus taxid allowlist from the ICTV Virus
Metadata Resource (VMR).

The VMR spreadsheet (https://ictv.global/vmr/current) has a controlled "Host
source" vocabulary — Archaea, Bacteria (bacteriophages), Fungi, Invertebrates,
Plants, Protists, Vertebrates, and combinations. We keep every virus family and
genus whose host source includes ``vertebrates`` (this includes arbovirus
families listed as ``invertebrates, vertebrates``), map those ICTV names to NCBI
taxids via a taxdump ``names.dmp``, and write the taxids one-per-line. That file
is the allowlist consumed by ``apply_ictv_host_filter.py``.

The result is a reproducible, user-editable artifact: re-run this to refresh it
against a new VMR / taxdump release. Pure parsing/resolution functions are unit
tested; ``main`` wires download + Excel parsing + name resolution.
"""

import argparse
import sys
from typing import Dict, Iterable, List, Set, Tuple

VMR_CURRENT_URL = "https://ictv.global/vmr/current"


def is_vertebrate_host(host_source) -> bool:
    """True iff the ICTV Host source string lists ``vertebrates`` as a token.

    Token-aware to avoid the substring trap where ``invertebrates`` contains
    ``vertebrates``.
    """
    if not host_source:
        return False
    text = str(host_source).lower()
    tokens = [t.strip() for chunk in text.split(",") for t in chunk.split(";")]
    return "vertebrates" in tokens


def extract_vertebrate_taxa(rows: Iterable[dict]) -> Set[str]:
    """Collect Family and Genus names from VMR rows whose host is vertebrate."""
    names: Set[str] = set()
    for row in rows:
        if not is_vertebrate_host(row.get("Host source")):
            continue
        for key in ("Family", "Genus"):
            value = (row.get(key) or "").strip()
            if value:
                names.add(value)
    return names


def build_name_index(names_dmp: str) -> Dict[str, str]:
    """Map lowercased scientific names -> taxid from a taxdump ``names.dmp``."""
    index: Dict[str, str] = {}
    with open(names_dmp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4 or parts[3] != "scientific name":
                continue
            taxid, name = parts[0], parts[1]
            index[name.lower()] = taxid
    return index


def resolve_names_to_taxids(
    names: Iterable[str], name_index: Dict[str, str]
) -> Tuple[Set[str], Set[str]]:
    """Return (taxids, unresolved_names) resolving ICTV names via ``name_index``."""
    taxids: Set[str] = set()
    unresolved: Set[str] = set()
    for name in names:
        taxid = name_index.get(name.lower())
        if taxid:
            taxids.add(taxid)
        else:
            unresolved.add(name)
    return taxids, unresolved


def _read_vmr_rows(vmr_path: str) -> List[dict]:
    """Read the VMR Excel/CSV into a list of column->value dicts, tolerating
    minor column-name variants (e.g. 'Host Source')."""
    import pandas as pd

    if vmr_path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(vmr_path)
    else:
        df = pd.read_csv(vmr_path, sep=None, engine="python")
    rename = {}
    for col in df.columns:
        low = str(col).strip().lower()
        if low == "host source":
            rename[col] = "Host source"
        elif low == "family":
            rename[col] = "Family"
        elif low == "genus":
            rename[col] = "Genus"
    df = df.rename(columns=rename)
    return df.to_dict(orient="records")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the vertebrate-virus taxid allowlist from the ICTV VMR."
    )
    parser.add_argument("--vmr", required=True, help="Path to a downloaded VMR .xlsx")
    parser.add_argument("--names", required=True, help="taxdump names.dmp")
    parser.add_argument("--output", required=True, help="Output taxid allowlist file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = _read_vmr_rows(args.vmr)
    taxa = extract_vertebrate_taxa(rows)
    index = build_name_index(args.names)
    taxids, unresolved = resolve_names_to_taxids(taxa, index)

    with open(args.output, "w") as out:
        out.write(
            "# Vertebrate-infecting virus taxids (ICTV VMR 'Host source' contains 'vertebrates').\n"
        )
        out.write(f"# Source: {VMR_CURRENT_URL}\n")
        out.write("# Regenerate with build_ictv_vertebrate_taxids.py; edit freely.\n")
        for taxid in sorted(taxids, key=int):
            out.write(f"{taxid}\n")

    print(
        f"[build_ictv_vertebrate_taxids] vertebrate taxa={len(taxa)} "
        f"resolved={len(taxids)} unresolved={len(unresolved)} -> {args.output}",
        file=sys.stderr,
    )
    if unresolved:
        print(
            "[build_ictv_vertebrate_taxids] unresolved names (not in taxdump): "
            + ", ".join(sorted(unresolved)[:20])
            + (" ..." if len(unresolved) > 20 else ""),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
