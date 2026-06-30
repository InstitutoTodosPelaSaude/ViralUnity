"""Tests for viralunity.scripts.python.build_genome_length_table."""

import os
import tempfile
import unittest

from viralunity.scripts.python.build_genome_length_table import (
    _parse_fai,
    _parse_genome2taxid,
    build_genome_length_table,
    write_table,
)

# Synthetic taxonomy (same as taxonomy_test.py):
#   1 (root, no rank)
#   └── 1000 (family)
#         └── 2000 (genus)
#               ├── 3001 (species)
#               └── 3002 (species)
#   └── 1100 (family)
#         └── 2100 (genus)
#               └── 3100 (species)

NODES_RECORDS = [
    ("1", "1", "no rank"),
    ("1000", "1", "family"),
    ("2000", "1000", "genus"),
    ("3001", "2000", "species"),
    ("3002", "2000", "species"),
    ("1100", "1", "family"),
    ("2100", "1100", "genus"),
    ("3100", "2100", "species"),
]

NAMES_RECORDS = [
    ("1", "root", "scientific name"),
    ("1000", "FamilyA", "scientific name"),
    ("2000", "GenusA", "scientific name"),
    ("3001", "SpeciesA1", "scientific name"),
    ("3002", "SpeciesA2", "scientific name"),
    ("1100", "FamilyB", "scientific name"),
    ("2100", "GenusB", "scientific name"),
    ("3100", "SpeciesB1", "scientific name"),
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


def _write_fai(directory: str, entries) -> str:
    """Write a samtools-faidx style .fai file.

    entries: list of (seq_name, length) tuples.
    The additional .fai columns (offset, linewidth, binwidth) are stubbed.
    """
    path = os.path.join(directory, "viral.fasta.fai")
    with open(path, "w") as f:
        for name, length in entries:
            f.write(f"{name}\t{length}\t0\t60\t61\n")
    return path


def _write_genome2taxid(directory: str, entries) -> str:
    """Write a genome2taxid TSV (accession<TAB>taxid, no header)."""
    path = os.path.join(directory, "genome2taxid.tsv")
    with open(path, "w") as f:
        for accession, taxid in entries:
            f.write(f"{accession}\t{taxid}\n")
    return path


class TestParseFai(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)

    def test_basic_parsing(self):
        path = _write_fai(self._tmp.name, [("NC_001461.1", 10700), ("NC_001542.1", 15894)])
        result = _parse_fai(path)
        self.assertEqual(result["NC_001461.1"], 10700)
        self.assertEqual(result["NC_001542.1"], 15894)

    def test_strips_description_from_name(self):
        # Real .fai files have "NC_001461.1 Dengue virus 1, complete genome" as name
        path = os.path.join(self._tmp.name, "test.fai")
        with open(path, "w") as f:
            f.write("NC_001461.1 Dengue virus 1, complete genome\t10700\t0\t60\t61\n")
        result = _parse_fai(path)
        # Key should be first whitespace-delimited token only
        self.assertIn("NC_001461.1", result)
        self.assertEqual(result["NC_001461.1"], 10700)
        self.assertNotIn("NC_001461.1 Dengue virus 1, complete genome", result)

    def test_empty_fai_returns_empty_dict(self):
        path = os.path.join(self._tmp.name, "empty.fai")
        open(path, "w").close()
        self.assertEqual(_parse_fai(path), {})

    def test_skips_malformed_lines(self):
        path = os.path.join(self._tmp.name, "bad.fai")
        with open(path, "w") as f:
            f.write("NC_001.1\n")  # only one column
            f.write("NC_002.1\t5000\t0\t60\t61\n")
        result = _parse_fai(path)
        self.assertNotIn("NC_001.1", result)
        self.assertEqual(result["NC_002.1"], 5000)


class TestParseGenome2Taxid(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)

    def test_basic_parsing(self):
        path = _write_genome2taxid(self._tmp.name, [("NC_001.1", "3001"), ("NC_002.1", "3002")])
        result = _parse_genome2taxid(path)
        self.assertEqual(result["NC_001.1"], "3001")
        self.assertEqual(result["NC_002.1"], "3002")

    def test_empty_file_returns_empty(self):
        path = os.path.join(self._tmp.name, "empty.tsv")
        open(path, "w").close()
        self.assertEqual(_parse_genome2taxid(path), {})

    def test_skips_malformed_lines(self):
        path = os.path.join(self._tmp.name, "bad.tsv")
        with open(path, "w") as f:
            f.write("NC_001.1\n")  # only one column
            f.write("NC_002.1\t3002\n")
        result = _parse_genome2taxid(path)
        self.assertNotIn("NC_001.1", result)
        self.assertEqual(result["NC_002.1"], "3002")


class TestBuildGenomeLengthTable(unittest.TestCase):
    """Integration tests for the full genome-length build pipeline."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.nodes, self.names = _write_taxdump(self._tmp.name)

    # Sentinel to distinguish "caller didn't pass names_dmp" from "caller passed None".
    _USE_DEFAULT_NAMES = object()

    def _run(self, fai_entries, g2t_entries, names_dmp=_USE_DEFAULT_NAMES):
        fai = _write_fai(self._tmp.name, fai_entries)
        g2t = _write_genome2taxid(self._tmp.name, g2t_entries)
        effective_names = self.names if names_dmp is self._USE_DEFAULT_NAMES else names_dmp
        return build_genome_length_table(
            fai_path=fai,
            genome2taxid_path=g2t,
            nodes_dmp=self.nodes,
            names_dmp=effective_names,
        )

    def _as_dict(self, rows):
        """Convert rows to {(rank,taxid): (length, n)} for easy assertion."""
        return {(r, t): (length, n) for r, t, _name, length, n in rows}

    def test_single_species_genome(self):
        rows = self._run(
            fai_entries=[("NC_001.1", 10000)],
            g2t_entries=[("NC_001.1", "3001")],
        )
        d = self._as_dict(rows)
        # Should have species 3001, genus 2000, family 1000
        self.assertIn(("species", "3001"), d)
        self.assertIn(("genus", "2000"), d)
        self.assertIn(("family", "1000"), d)
        # All have length 10000
        self.assertEqual(d[("species", "3001")][0], 10000)
        self.assertEqual(d[("genus", "2000")][0], 10000)
        self.assertEqual(d[("family", "1000")][0], 10000)

    def test_median_across_multiple_genomes_same_species(self):
        # Three genomes under species 3001: lengths 5000, 10000, 15000 → median 10000
        rows = self._run(
            fai_entries=[("A", 5000), ("B", 10000), ("C", 15000)],
            g2t_entries=[("A", "3001"), ("B", "3001"), ("C", "3001")],
        )
        d = self._as_dict(rows)
        self.assertEqual(d[("species", "3001")][0], 10000)
        self.assertEqual(d[("species", "3001")][1], 3)  # n_genomes

    def test_median_with_two_genomes_is_average(self):
        # Two genomes: 5000, 15000 → median 10000
        rows = self._run(
            fai_entries=[("A", 5000), ("B", 15000)],
            g2t_entries=[("A", "3001"), ("B", "3001")],
        )
        d = self._as_dict(rows)
        self.assertEqual(d[("species", "3001")][0], 10000)

    def test_two_species_same_genus_genus_median(self):
        # SpeciesA1 (3001) has one genome of 5000.
        # SpeciesA2 (3002) has one genome of 15000.
        # GenusA (2000) should have median of [5000, 15000] = 10000.
        rows = self._run(
            fai_entries=[("A", 5000), ("B", 15000)],
            g2t_entries=[("A", "3001"), ("B", "3002")],
        )
        d = self._as_dict(rows)
        # Each species has its own single length
        self.assertEqual(d[("species", "3001")][0], 5000)
        self.assertEqual(d[("species", "3002")][0], 15000)
        # Genus-level median across both species' genomes
        self.assertEqual(d[("genus", "2000")][0], 10000)
        self.assertEqual(d[("genus", "2000")][1], 2)

    def test_two_families_distinct(self):
        # One genome under family 1000 and one under 1100 (different families)
        rows = self._run(
            fai_entries=[("A", 10000), ("B", 20000)],
            g2t_entries=[("A", "3001"), ("B", "3100")],  # 3001→1000, 3100→1100
        )
        d = self._as_dict(rows)
        # Families should be distinct
        self.assertEqual(d[("family", "1000")][0], 10000)
        self.assertEqual(d[("family", "1100")][0], 20000)

    def test_accession_missing_in_genome2taxid_is_skipped(self):
        rows = self._run(
            fai_entries=[("NC_KNOWN.1", 10000), ("NC_UNKNOWN.1", 5000)],
            g2t_entries=[("NC_KNOWN.1", "3001")],  # NC_UNKNOWN.1 has no taxid
        )
        d = self._as_dict(rows)
        # Only the known accession contributes
        self.assertEqual(d[("species", "3001")][1], 1)

    def test_taxid_with_no_rank_of_interest_not_in_output(self):
        # If a taxid has no family/genus/species ancestor, it produces no rows.
        # Taxid "99999" is not in our taxonomy, so it has no lineage at any rank.
        rows = self._run(
            fai_entries=[("UNKNOWN.1", 5000)],
            g2t_entries=[("UNKNOWN.1", "99999")],
        )
        # No rows should be emitted (no lineage nodes at ranks of interest)
        self.assertEqual(rows, [])

    def test_without_names_dmp_name_column_is_na(self):
        rows = self._run(
            fai_entries=[("A", 10000)],
            g2t_entries=[("A", "3001")],
            names_dmp=None,  # no names file
        )
        d = {(r, t): name for r, t, name, *_ in rows}
        self.assertEqual(d[("species", "3001")], "NA")

    def test_with_names_dmp_name_column_is_populated(self):
        rows = self._run(
            fai_entries=[("A", 10000)],
            g2t_entries=[("A", "3001")],
        )
        d = {(r, t): name for r, t, name, *_ in rows}
        self.assertEqual(d[("species", "3001")], "SpeciesA1")

    def test_n_genomes_counts_individual_accessions_not_nodes(self):
        # 2 accessions assigned to 3001, both contribute at species, genus, family level
        rows = self._run(
            fai_entries=[("A", 4000), ("B", 6000)],
            g2t_entries=[("A", "3001"), ("B", "3001")],
        )
        d = self._as_dict(rows)
        self.assertEqual(d[("species", "3001")][1], 2)
        self.assertEqual(d[("genus", "2000")][1], 2)
        self.assertEqual(d[("family", "1000")][1], 2)


class TestWriteTable(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)

    def test_write_creates_file_with_header(self):
        path = os.path.join(self._tmp.name, "out.tsv")
        rows = [("species", "3001", "SpeciesA1", 10000, 5)]
        write_table(rows, path)
        with open(path) as f:
            lines = f.readlines()
        self.assertEqual(lines[0].strip(), "rank\ttaxid\tname\tgenome_length_bp\tn_genomes")
        self.assertEqual(lines[1].strip(), "species\t3001\tSpeciesA1\t10000\t5")

    def test_write_empty_table(self):
        path = os.path.join(self._tmp.name, "empty.tsv")
        write_table([], path)
        with open(path) as f:
            lines = f.readlines()
        # Header only
        self.assertEqual(len(lines), 1)
        self.assertIn("rank", lines[0])

    def test_write_creates_parent_dirs(self):
        path = os.path.join(self._tmp.name, "subdir", "out.tsv")
        write_table([("species", "1", "root", 5000, 1)], path)
        self.assertTrue(os.path.exists(path))


if __name__ == "__main__":
    unittest.main()
