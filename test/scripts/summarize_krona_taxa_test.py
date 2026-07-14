"""Tests for viralunity.scripts.python.summarize_krona_taxa."""

import os
import tempfile
import unittest

from viralunity.scripts.python.summarize_krona_taxa import load_diamond_reads, run, summarize_krona


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


class TestRunHeader(unittest.TestCase):
    def _write_taxdump(self, d):
        nodes = os.path.join(d, "nodes.dmp")
        names = os.path.join(d, "names.dmp")
        # 111 (species) -> 1 (root)
        with open(nodes, "w") as f:
            f.write("111\t|\t1\t|\tspecies\t|\n")
            f.write("1\t|\t1\t|\tno rank\t|\n")
        with open(names, "w") as f:
            f.write("111\t|\tVirus A\t|\t\t|\tscientific name\t|\n")
        return nodes, names

    def test_output_header_has_no_source_column(self):
        with tempfile.TemporaryDirectory() as d:
            nodes, names = self._write_taxdump(d)
            krona = os.path.join(d, "k.tsv")
            with open(krona, "w") as f:
                f.write("contig1\t111\n")
            out = os.path.join(d, "summary.tsv")

            run(
                krona=krona,
                diamond_tax=None,
                taxdump_nodes=nodes,
                taxdump_names=names,
                sample="sample-A",
                classifier="kraken2",
                unit="reads",
                output=out,
            )

            with open(out) as f:
                header = f.readline().rstrip("\n").split("\t")

            self.assertNotIn("source", header)
            self.assertEqual(
                header,
                ["sample", "tool", "mode", "rank", "taxid", "name", "count", "percent"],
            )


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
