#!/usr/bin/env python3
"""Shared NCBI taxonomy utilities used across ViralUnity metagenomics scripts.

Centralises the three helpers that were previously copied between
summarize_krona_taxa.py and filter_krona_by_pass_taxids.py.
"""

from typing import Dict, List, Optional, Tuple

# The three ranks at which the metagenomics pipeline summarises hits.
RANKS_OF_INTEREST: Tuple[str, ...] = ("family", "genus", "species")


def load_taxdump(
    nodes_dmp: str,
    names_dmp: Optional[str] = None,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Parse nodes.dmp and optionally names.dmp.

    Parameters
    ----------
    nodes_dmp:
        Path to NCBI taxdump ``nodes.dmp``.
    names_dmp:
        Path to NCBI taxdump ``names.dmp``. Optional; when omitted the
        returned name map is an empty dict.

    Returns
    -------
    parent_map : dict[taxid, parent_taxid]
    rank_map   : dict[taxid, rank_string]
    name_map   : dict[taxid, scientific_name]  (empty when names_dmp is None)
    """
    parent: Dict[str, str] = {}
    rank: Dict[str, str] = {}
    name: Dict[str, str] = {}

    with open(nodes_dmp) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]
            parent[parts[0]] = parts[1]
            rank[parts[0]] = parts[2]

    if names_dmp is not None:
        with open(names_dmp) as fh:
            for line in fh:
                if not line.strip():
                    continue
                parts = [p.strip() for p in line.split("|")]
                if parts[3] == "scientific name":
                    name[parts[0]] = parts[1]

    return parent, rank, name


def get_lineage(taxid: str, parent_map: Dict[str, str]) -> List[str]:
    """Return the full lineage list from *taxid* up to (and including) root '1'.

    Parameters
    ----------
    taxid:
        Starting taxid as a string.
    parent_map:
        Dict mapping taxid → parent taxid, as returned by :func:`load_taxdump`.

    Returns
    -------
    List of taxid strings ``[taxid, parent, grandparent, ..., '1']``.
    """
    lineage: List[str] = []
    while taxid != "1" and taxid in parent_map:
        lineage.append(taxid)
        taxid = parent_map[taxid]
    lineage.append("1")
    return lineage
