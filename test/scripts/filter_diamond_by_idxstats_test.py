"""Tests for viralunity.scripts.python.filter_diamond_by_idxstats."""

import os
import tempfile
import unittest
from pathlib import Path

from viralunity.scripts.python.filter_diamond_by_idxstats import (
    diamond_base_ids,
    filter_and_annotate,
    load_mapped_counts,
)


def _write(d, name, content):
    p = os.path.join(d, name)
    with open(p, "w") as f:
        f.write(content)
    return Path(p)


class TestLoadMappedCounts(unittest.TestCase):
    def test_sums_subregions_and_skips_unmapped_star(self):
        with tempfile.TemporaryDirectory() as d:
            idx = _write(
                d,
                "idx.tsv",
                "contigA|r1\t100\t5\t0\ncontigA|r2\t100\t3\t0\ncontigB\t50\t0\t0\n*\t0\t0\t0\n",
            )
            counts = load_mapped_counts(idx)
            self.assertEqual(counts["contigA"], 8)  # 5 + 3 summed across sub-regions
            self.assertEqual(counts["contigB"], 0)
            self.assertNotIn("*", counts)


class TestDiamondBaseIds(unittest.TestCase):
    def test_strips_pipe_suffix_and_skips_comments(self):
        with tempfile.TemporaryDirectory() as d:
            dmd = _write(d, "d.tsv", "# header\ncontigA|orf1\tsseq\t99\ncontigB\tsseq\t80\n")
            self.assertEqual(diamond_base_ids(dmd), {"contigA", "contigB"})


class TestFilterAndAnnotate(unittest.TestCase):
    def test_keeps_supported_contigs_and_appends_count(self):
        with tempfile.TemporaryDirectory() as d:
            dmd = _write(d, "d.tsv", "contigA|orf1\tsseq\t99\ncontigB\tsseq\t80\n")
            out = Path(os.path.join(d, "out.tsv"))
            mapped = {"contigA": 8, "contigB": 0}

            kept_rows, kept_contigs = filter_and_annotate(dmd, out, mapped, min_mapped=1)

            self.assertEqual((kept_rows, kept_contigs), (1, 1))
            lines = [ln for ln in out.read_text().splitlines() if ln.strip()]
            self.assertEqual(
                lines, ["contigA|orf1\tsseq\t99\t8"]
            )  # contigB dropped, count appended


if __name__ == "__main__":
    unittest.main()
