"""Tests for viralunity.scripts.python.summarize_krona_taxa."""

import os
import tempfile
import unittest

from viralunity.scripts.python.summarize_krona_taxa import load_diamond_reads, summarize_krona


class TestSummarizeKrona(unittest.TestCase):
    def test_skips_malformed_line(self):
        with tempfile.TemporaryDirectory() as d:
            krona = os.path.join(d, "k.tsv")
            with open(krona, "w") as f:
                f.write("contig1\t111\n")
                f.write("BADLINE_NO_TAB\n")  # only one field -> must not crash

            parent_map = {"111": "1"}
            rank_map = {"111": "species"}

            contig_counts, _read_counts, _totals = summarize_krona(krona, parent_map, rank_map)

            self.assertEqual(contig_counts[("species", "111")], 1)


class TestLoadDiamondReads(unittest.TestCase):
    def test_reads_last_column_and_skips_comments(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "annotated.tsv")
            with open(p, "w") as f:
                f.write("# comment\n")
                f.write("contig1\t111\tVirus A\t7\n")
            reads = load_diamond_reads(p)
            self.assertEqual(reads["contig1"], 7)


if __name__ == "__main__":
    unittest.main()
