"""Tests for viralunity.scripts.python.add_RPM_to_summary."""

import gzip
import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.add_RPM_to_summary import add_rpm, count_fastq_reads


class TestCountFastqReads(unittest.TestCase):
    def test_plain_and_gzip(self):
        with tempfile.TemporaryDirectory() as d:
            plain = os.path.join(d, "r.fastq")
            with open(plain, "w") as f:
                f.write("@r1\nACGT\n+\n!!!!\n@r2\nACGT\n+\n!!!!\n")
            self.assertEqual(count_fastq_reads(plain), 2)

            gz = os.path.join(d, "r.fastq.gz")
            with gzip.open(gz, "wt") as f:
                f.write("@r1\nACGT\n+\n!!!!\n")
            self.assertEqual(count_fastq_reads(gz), 1)

    def test_empty_file_is_zero(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "empty.fastq")
            open(p, "w").close()
            self.assertEqual(count_fastq_reads(p), 0)


class TestAddRpm(unittest.TestCase):
    def _fastq(self, d, n_reads):
        p = os.path.join(d, "reads.fastq")
        with open(p, "w") as f:
            f.write("@r\nACGT\n+\n!!!!\n" * n_reads)
        return p

    def test_rpm_math(self):
        with tempfile.TemporaryDirectory() as d:
            fastq = self._fastq(d, 2)  # 2 total reads
            df = pd.DataFrame({"sample": ["S1"], "mapped_reads": [1]})
            out = add_rpm(df, {"S1": fastq}, reads_col="mapped_reads", rpm_col="rpm")
            self.assertEqual(int(out["total_reads"].iloc[0]), 2)
            self.assertAlmostEqual(out["rpm"].iloc[0], 1 / 2 * 1e6)

    def test_zero_denominator_yields_zero_rpm(self):
        with tempfile.TemporaryDirectory() as d:
            fastq = self._fastq(d, 0)  # empty -> 0 total reads
            df = pd.DataFrame({"sample": ["S1"], "mapped_reads": [5]})
            out = add_rpm(df, {"S1": fastq}, reads_col="mapped_reads", rpm_col="rpm")
            self.assertEqual(out["rpm"].iloc[0], 0.0)

    def test_missing_reads_column_defaults_to_zero(self):
        with tempfile.TemporaryDirectory() as d:
            fastq = self._fastq(d, 3)
            df = pd.DataFrame({"sample": ["S1"]})  # no mapped_reads column
            out = add_rpm(df, {"S1": fastq}, reads_col="mapped_reads", rpm_col="rpm")
            self.assertEqual(out["mapped_reads"].iloc[0], 0)
            self.assertEqual(out["rpm"].iloc[0], 0.0)

    def test_unknown_sample_raises(self):
        with tempfile.TemporaryDirectory() as d:
            fastq = self._fastq(d, 1)
            df = pd.DataFrame({"sample": ["S1", "S2"], "mapped_reads": [1, 1]})
            with self.assertRaises(KeyError):
                add_rpm(df, {"S1": fastq}, reads_col="mapped_reads", rpm_col="rpm")


if __name__ == "__main__":
    unittest.main()
