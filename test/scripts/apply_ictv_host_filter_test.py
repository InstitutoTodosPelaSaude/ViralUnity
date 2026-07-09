"""Tests for viralunity.scripts.python.apply_ictv_host_filter."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.apply_ictv_host_filter import (
    filter_summary,
    lineage_allowed,
    load_allowlist,
    run,
)

# Synthetic taxonomy:
#   1 (root)
#     ├── 1000 (family, vertebrate virus)  <- ALLOWED
#     │     └── 2000 (genus)
#     │           └── 3001 (species)
#     │                 └── 4001 (strain)
#     └── 1100 (family, bacteriophage)      <- NOT allowed
#           └── 2100 (genus)
#                 └── 3100 (species)
NODES_RECORDS = [
    ("1", "1", "no rank"),
    ("1000", "1", "family"),
    ("2000", "1000", "genus"),
    ("3001", "2000", "species"),
    ("4001", "3001", "no rank"),
    ("1100", "1", "family"),
    ("2100", "1100", "genus"),
    ("3100", "2100", "species"),
]
NAMES_RECORDS = [
    ("1", "root"),
    ("1000", "Orthomyxoviridae"),
    ("2000", "Alphainfluenzavirus"),
    ("3001", "Influenza A virus"),
    ("4001", "Influenza A strain X"),
    ("1100", "Steitzviridae"),
    ("2100", "Lambdavirus"),
    ("3100", "Lambdavirus DE3"),
]


def _write_taxdump(directory):
    nodes = os.path.join(directory, "nodes.dmp")
    names = os.path.join(directory, "names.dmp")
    with open(nodes, "w") as f:
        for taxid, parent, rank in NODES_RECORDS:
            f.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t-\t|\n")
    with open(names, "w") as f:
        for taxid, name in NAMES_RECORDS:
            f.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    return nodes, names


def _summary_df():
    rows = [
        ("A", "kraken2", "reads", "species", "3001", "Influenza A virus"),
        ("A", "kraken2", "reads", "no rank", "4001", "Influenza A strain X"),
        ("A", "kraken2", "reads", "family", "1000", "Orthomyxoviridae"),
        ("A", "kraken2", "reads", "species", "3100", "Lambdavirus DE3"),
        ("A", "kraken2", "reads", "family", "1100", "Steitzviridae"),
        ("A", "kraken2", "reads", "species", "0", "unclassified"),
        ("A", "kraken2", "reads", "species", "999999", "Ghost virus"),
    ]
    return pd.DataFrame(rows, columns=["sample", "tool", "mode", "rank", "taxid", "name"])


class TestLoadAllowlist(unittest.TestCase):
    def test_parses_taxids_ignoring_comments_and_blanks(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "allow.txt")
            with open(path, "w") as f:
                f.write("# vertebrate-infecting virus families\n")
                f.write("1000\n")
                f.write("\n")
                f.write("2000  # inline comment kept out\n")
            allow = load_allowlist(path)
            self.assertEqual(allow, {"1000", "2000"})


class TestLineageAllowed(unittest.TestCase):
    def setUp(self):
        from viralunity.scripts.python.taxonomy import load_taxdump

        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        nodes, names = _write_taxdump(self._tmp.name)
        self.parent, _rank, _name = load_taxdump(nodes, names)

    def test_species_under_allowed_family_is_allowed(self):
        self.assertTrue(lineage_allowed("3001", {"1000"}, self.parent))

    def test_strain_under_allowed_family_is_allowed(self):
        self.assertTrue(lineage_allowed("4001", {"1000"}, self.parent))

    def test_allowed_family_itself_is_allowed(self):
        self.assertTrue(lineage_allowed("1000", {"1000"}, self.parent))

    def test_taxon_outside_allowlist_is_rejected(self):
        self.assertFalse(lineage_allowed("3100", {"1000"}, self.parent))

    def test_taxid_zero_rejected(self):
        self.assertFalse(lineage_allowed("0", {"1000"}, self.parent))

    def test_unknown_taxid_rejected(self):
        self.assertFalse(lineage_allowed("999999", {"1000"}, self.parent))


class TestFilterSummary(unittest.TestCase):
    def setUp(self):
        from viralunity.scripts.python.taxonomy import load_taxdump

        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        nodes, names = _write_taxdump(self._tmp.name)
        self.parent, _rank, _name = load_taxdump(nodes, names)
        self.in_path = os.path.join(self._tmp.name, "in.tsv")
        self.out_path = os.path.join(self._tmp.name, "out.tsv")
        self.dropped_path = os.path.join(self._tmp.name, "dropped.tsv")
        _summary_df().to_csv(self.in_path, sep="\t", index=False)

    def test_keeps_only_vertebrate_lineage_rows(self):
        kept, dropped = filter_summary(
            self.in_path, self.out_path, self.dropped_path, {"1000"}, self.parent
        )
        self.assertEqual(kept, 3)
        self.assertEqual(dropped, 4)
        out = pd.read_csv(self.out_path, sep="\t", dtype=str)
        self.assertEqual(set(out["taxid"]), {"3001", "4001", "1000"})

    def test_dropped_sidecar_lists_removed_rows_with_reason(self):
        filter_summary(self.in_path, self.out_path, self.dropped_path, {"1000"}, self.parent)
        dropped = pd.read_csv(self.dropped_path, sep="\t", dtype=str)
        self.assertEqual(set(dropped["taxid"]), {"3100", "1100", "0", "999999"})
        self.assertTrue((dropped["drop_reason"] == "not_vertebrate_virus").all())

    def test_empty_input_yields_empty_outputs(self):
        open(self.in_path, "w").close()
        kept, dropped = filter_summary(
            self.in_path, self.out_path, self.dropped_path, {"1000"}, self.parent
        )
        self.assertEqual((kept, dropped), (0, 0))
        self.assertTrue(os.path.exists(self.out_path))

    def test_empty_allowlist_drops_everything(self):
        kept, dropped = filter_summary(
            self.in_path, self.out_path, self.dropped_path, set(), self.parent
        )
        self.assertEqual(kept, 0)
        self.assertEqual(dropped, 7)


class TestRunEndToEnd(unittest.TestCase):
    def test_run_glues_allowlist_and_taxdump(self):
        with tempfile.TemporaryDirectory() as d:
            nodes, names = _write_taxdump(d)
            in_path = os.path.join(d, "in.tsv")
            out_path = os.path.join(d, "out.tsv")
            dropped_path = os.path.join(d, "dropped.tsv")
            allow_path = os.path.join(d, "allow.txt")
            _summary_df().to_csv(in_path, sep="\t", index=False)
            with open(allow_path, "w") as f:
                f.write("1000\n")
            run(
                summary=in_path,
                output=out_path,
                dropped=dropped_path,
                allowlist_file=allow_path,
                nodes_dmp=nodes,
                names_dmp=names,
            )
            out = pd.read_csv(out_path, sep="\t", dtype=str)
            self.assertEqual(set(out["taxid"]), {"3001", "4001", "1000"})


if __name__ == "__main__":
    unittest.main()
