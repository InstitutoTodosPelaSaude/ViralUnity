"""Tests for viralunity.scripts.python.taxonomy (shared taxonomy helpers)."""

import os
import tempfile
import unittest

from viralunity.scripts.python.taxonomy import (
    RANKS_OF_INTEREST,
    get_lineage,
    load_taxdump,
)

# Synthetic taxonomy (same structure reused across the test suite):
#   1 (root, no rank)
#   └── 1000 (family)
#         └── 2000 (genus)
#               ├── 3001 (species)
#               │     └── 4001 (no rank / strain)
#               └── 3002 (species)
#   └── 1100 (family)

NODES_RECORDS = [
    ("1", "1", "no rank"),
    ("1000", "1", "family"),
    ("2000", "1000", "genus"),
    ("3001", "2000", "species"),
    ("3002", "2000", "species"),
    ("4001", "3001", "no rank"),
    ("1100", "1", "family"),
]

NAMES_RECORDS = [
    ("1", "root", "scientific name"),
    ("1000", "FamilyA", "scientific name"),
    ("2000", "GenusA", "scientific name"),
    ("3001", "SpeciesA1", "scientific name"),
    ("3002", "SpeciesA2", "scientific name"),
    ("4001", "StrainA1a", "scientific name"),
    ("1100", "FamilyB", "scientific name"),
    # A synonym row that must NOT overwrite the scientific name:
    ("1000", "OldFamilyA", "synonym"),
]


def _write_taxdump(directory: str):
    nodes_path = os.path.join(directory, "nodes.dmp")
    names_path = os.path.join(directory, "names.dmp")
    with open(nodes_path, "w") as f:
        for taxid, parent, rank in NODES_RECORDS:
            f.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t-\t|\n")
    with open(names_path, "w") as f:
        for taxid, name, name_class in NAMES_RECORDS:
            f.write(f"{taxid}\t|\t{name}\t|\t\t|\t{name_class}\t|\n")
    return nodes_path, names_path


class TestRanksOfInterest(unittest.TestCase):
    def test_contains_expected_ranks(self):
        self.assertIn("family", RANKS_OF_INTEREST)
        self.assertIn("genus", RANKS_OF_INTEREST)
        self.assertIn("species", RANKS_OF_INTEREST)

    def test_does_not_contain_strain(self):
        self.assertNotIn("strain", RANKS_OF_INTEREST)
        self.assertNotIn("no rank", RANKS_OF_INTEREST)


class TestLoadTaxdump(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.nodes, self.names = _write_taxdump(self._tmp.name)

    # ---- full (3-tuple) mode ------------------------------------------------

    def test_returns_three_dicts(self):
        result = load_taxdump(self.nodes, self.names)
        self.assertEqual(len(result), 3)
        parent, rank, name = result
        self.assertIsInstance(parent, dict)
        self.assertIsInstance(rank, dict)
        self.assertIsInstance(name, dict)

    def test_parent_map_correct(self):
        parent, _, _ = load_taxdump(self.nodes, self.names)
        self.assertEqual(parent["4001"], "3001")
        self.assertEqual(parent["3001"], "2000")
        self.assertEqual(parent["2000"], "1000")
        self.assertEqual(parent["1000"], "1")

    def test_rank_map_correct(self):
        _, rank, _ = load_taxdump(self.nodes, self.names)
        self.assertEqual(rank["1000"], "family")
        self.assertEqual(rank["2000"], "genus")
        self.assertEqual(rank["3001"], "species")
        self.assertEqual(rank["4001"], "no rank")

    def test_name_map_contains_scientific_names_only(self):
        _, _, name = load_taxdump(self.nodes, self.names)
        self.assertEqual(name["1000"], "FamilyA")
        self.assertEqual(name["3001"], "SpeciesA1")
        # Synonym row must NOT appear as the name for 1000.
        self.assertNotEqual(name.get("1000"), "OldFamilyA")

    def test_name_map_excludes_synonym_rows(self):
        _, _, name = load_taxdump(self.nodes, self.names)
        # No name dict entry should have synonym value.
        self.assertNotIn("OldFamilyA", name.values())

    # ---- nodes-only mode (names_dmp=None) ------------------------------------

    def test_nodes_only_returns_empty_name_map(self):
        parent, rank, name = load_taxdump(self.nodes)
        self.assertEqual(name, {})
        # parent and rank should still be populated.
        self.assertIn("3001", parent)
        self.assertIn("3001", rank)

    def test_nodes_only_parent_map_identical_to_full(self):
        parent_full, _, _ = load_taxdump(self.nodes, self.names)
        parent_nodes, _, _ = load_taxdump(self.nodes)
        self.assertEqual(parent_full, parent_nodes)

    def test_nodes_only_rank_map_identical_to_full(self):
        _, rank_full, _ = load_taxdump(self.nodes, self.names)
        _, rank_nodes, _ = load_taxdump(self.nodes)
        self.assertEqual(rank_full, rank_nodes)

    # ---- empty file edge cases ------------------------------------------------

    def test_empty_nodes_file(self):
        empty_nodes = os.path.join(self._tmp.name, "empty_nodes.dmp")
        open(empty_nodes, "w").close()
        parent, rank, name = load_taxdump(empty_nodes)
        self.assertEqual(parent, {})
        self.assertEqual(rank, {})
        self.assertEqual(name, {})


class TestGetLineage(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        nodes, names = _write_taxdump(self._tmp.name)
        self.parent, _, _ = load_taxdump(nodes, names)

    def test_lineage_includes_self_and_root(self):
        lineage = get_lineage("4001", self.parent)
        self.assertEqual(lineage, ["4001", "3001", "2000", "1000", "1"])

    def test_lineage_from_species(self):
        lineage = get_lineage("3001", self.parent)
        self.assertEqual(lineage, ["3001", "2000", "1000", "1"])

    def test_lineage_from_family(self):
        lineage = get_lineage("1000", self.parent)
        self.assertEqual(lineage, ["1000", "1"])

    def test_lineage_from_root_is_just_root(self):
        lineage = get_lineage("1", self.parent)
        self.assertEqual(lineage, ["1"])

    def test_unknown_taxid_returns_list_with_root(self):
        # If taxid is not in parent_map but is also not "1", the loop exits
        # immediately and "1" is appended.
        lineage = get_lineage("99999", self.parent)
        self.assertEqual(lineage, ["1"])

    def test_lineage_terminates_at_root_not_infinite(self):
        # Ensure no infinite loops for known taxa.
        for tid in ["1000", "2000", "3001", "3002", "4001"]:
            lineage = get_lineage(tid, self.parent)
            self.assertIn("1", lineage)
            self.assertLess(len(lineage), 100)


if __name__ == "__main__":
    unittest.main()
