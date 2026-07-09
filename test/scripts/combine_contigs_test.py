"""Tests for viralunity.scripts.python.combine_contigs."""

import os
import tempfile
import unittest

from viralunity.scripts.python.combine_contigs import combine_contigs


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def _read_fasta(path):
    """Return list of (header_without_gt, seq) preserving order."""
    records = []
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            elif line:
                seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


class TestCombineContigs(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name
        self.out = os.path.join(self.d, "combined.fasta")

    def test_prefixes_headers_with_sample_label(self):
        a = os.path.join(self.d, "a.fa")
        b = os.path.join(self.d, "b.fa")
        _write_fasta(a, [("k141_1", "ACGT"), ("k141_2", "TTTT")])
        _write_fasta(b, [("k141_1", "GGGG")])
        n = combine_contigs([("sample-A", a), ("sample-B", b)], self.out)
        self.assertEqual(n, 3)
        self.assertEqual(
            _read_fasta(self.out),
            [
                ("sample-A|k141_1", "ACGT"),
                ("sample-A|k141_2", "TTTT"),
                ("sample-B|k141_1", "GGGG"),
            ],
        )

    def test_strips_description_after_first_whitespace(self):
        a = os.path.join(self.d, "a.fa")
        _write_fasta(a, [("k141_9 flag=1 multi=3.0 len=500", "ACGTACGT")])
        combine_contigs([("sample-A", a)], self.out)
        self.assertEqual(_read_fasta(self.out), [("sample-A|k141_9", "ACGTACGT")])

    def test_missing_file_is_skipped(self):
        a = os.path.join(self.d, "a.fa")
        _write_fasta(a, [("k141_1", "ACGT")])
        missing = os.path.join(self.d, "does_not_exist.fa")
        n = combine_contigs([("sample-A", a), ("sample-B", missing)], self.out)
        self.assertEqual(n, 1)
        self.assertEqual(_read_fasta(self.out), [("sample-A|k141_1", "ACGT")])

    def test_empty_file_is_skipped(self):
        a = os.path.join(self.d, "a.fa")
        empty = os.path.join(self.d, "empty.fa")
        _write_fasta(a, [("k141_1", "ACGT")])
        open(empty, "w").close()
        n = combine_contigs([("sample-A", a), ("sample-B", empty)], self.out)
        self.assertEqual(n, 1)

    def test_all_empty_yields_empty_output(self):
        empty = os.path.join(self.d, "empty.fa")
        open(empty, "w").close()
        n = combine_contigs([("sample-A", empty)], self.out)
        self.assertEqual(n, 0)
        self.assertTrue(os.path.exists(self.out))
        self.assertEqual(_read_fasta(self.out), [])

    def test_bare_gt_header_is_skipped_not_crashing(self):
        # A malformed '>' header with no contig id must not raise IndexError.
        a = os.path.join(self.d, "a.fa")
        with open(a, "w") as fh:
            fh.write(">\n")  # bare header, empty defline
            fh.write("ACGT\n")
            fh.write(">k141_2\nGGGG\n")
        n = combine_contigs([("sample-A", a)], self.out)
        # only the well-formed record is emitted
        self.assertEqual(n, 1)
        self.assertEqual(_read_fasta(self.out), [("sample-A|k141_2", "GGGG")])


if __name__ == "__main__":
    unittest.main()
