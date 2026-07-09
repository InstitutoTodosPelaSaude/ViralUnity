"""Tests for viralunity.scripts.python.split_search_output."""

import os
import tempfile
import unittest

from viralunity.scripts.python.split_search_output import split_search_output


def _write_tsv(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(row) + "\n")


def _read_lines(path):
    with open(path) as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


class TestSplitSearchOutput(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name
        self.combined = os.path.join(self.d, "combined.tsv")
        self.out_a = os.path.join(self.d, "A.tsv")
        self.out_b = os.path.join(self.d, "B.tsv")

    def test_routes_rows_by_prefix_and_strips_prefix(self):
        _write_tsv(
            self.combined,
            [
                ["sample-A|k141_1", "sseq1", "99.0"],
                ["sample-B|k141_1", "sseq2", "88.0"],
                ["sample-A|k141_2", "sseq3", "77.0"],
            ],
        )
        counts = split_search_output(
            self.combined, {"sample-A": self.out_a, "sample-B": self.out_b}
        )
        self.assertEqual(counts, {"sample-A": 2, "sample-B": 1})
        self.assertEqual(
            _read_lines(self.out_a),
            ["k141_1\tsseq1\t99.0", "k141_2\tsseq3\t77.0"],
        )
        self.assertEqual(_read_lines(self.out_b), ["k141_1\tsseq2\t88.0"])

    def test_sample_with_no_hits_gets_empty_file(self):
        _write_tsv(self.combined, [["sample-A|k141_1", "sseq1", "99.0"]])
        counts = split_search_output(
            self.combined, {"sample-A": self.out_a, "sample-B": self.out_b}
        )
        self.assertEqual(counts["sample-B"], 0)
        self.assertTrue(os.path.exists(self.out_b))
        self.assertEqual(_read_lines(self.out_b), [])

    def test_only_first_delimiter_splits_sample_from_contig(self):
        # contig ids never contain '|', but be defensive: split once.
        _write_tsv(self.combined, [["sample-A|k141_1|region2", "sseq1", "50.0"]])
        split_search_output(self.combined, {"sample-A": self.out_a})
        self.assertEqual(_read_lines(self.out_a), ["k141_1|region2\tsseq1\t50.0"])

    def test_unknown_prefix_is_skipped(self):
        _write_tsv(
            self.combined,
            [
                ["sample-A|k141_1", "sseq1", "99.0"],
                ["sample-Z|k141_9", "sseqZ", "10.0"],
            ],
        )
        counts = split_search_output(self.combined, {"sample-A": self.out_a})
        self.assertEqual(counts, {"sample-A": 1})
        self.assertEqual(_read_lines(self.out_a), ["k141_1\tsseq1\t99.0"])

    def test_empty_combined_input_yields_all_empty_outputs(self):
        open(self.combined, "w").close()
        counts = split_search_output(
            self.combined, {"sample-A": self.out_a, "sample-B": self.out_b}
        )
        self.assertEqual(counts, {"sample-A": 0, "sample-B": 0})
        self.assertEqual(_read_lines(self.out_a), [])
        self.assertEqual(_read_lines(self.out_b), [])

    def test_missing_combined_input_yields_all_empty_outputs(self):
        counts = split_search_output(os.path.join(self.d, "nope.tsv"), {"sample-A": self.out_a})
        self.assertEqual(counts, {"sample-A": 0})
        self.assertEqual(_read_lines(self.out_a), [])


if __name__ == "__main__":
    unittest.main()
