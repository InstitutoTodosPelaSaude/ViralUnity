"""Tests for viralunity.scripts.python.convert_diamond_output_to_krona_input."""

import os
import subprocess
import sys
import tempfile
import unittest

from viralunity.scripts.python import convert_diamond_output_to_krona_input as mod

SCRIPT = mod.__file__


def _write(d, name, content):
    p = os.path.join(d, name)
    with open(p, "w") as f:
        f.write(content)
    return p


class TestMain(unittest.TestCase):
    def test_main_maps_reads_to_taxids(self):
        with tempfile.TemporaryDirectory() as d:
            diamond = _write(d, "d.tsv", "read1\tACC1|prot\t99\nread2\tACC2|prot\t99\n")
            seqs = _write(d, "s.fasta", ">read1\nAAAA\n>read2\nCCCC\n>read3\nGGGG\n")
            taxids = _write(d, "t.tsv", "ACC1\t111\nACC2\t222\n")
            out = os.path.join(d, "out.tsv")

            mod.main([diamond, seqs, taxids], out, "fasta")

            with open(out) as f:
                rows = dict(line.split() for line in f if line.strip())
            self.assertEqual(rows["read1"], "111")
            self.assertEqual(rows["read2"], "222")
            self.assertEqual(rows["read3"], "0")  # in sequences but no DIAMOND hit


class TestCli(unittest.TestCase):
    """The CLI branch must run (regression: argparse was never imported -> NameError)."""

    def test_cli_runs_without_nameerror(self):
        with tempfile.TemporaryDirectory() as d:
            diamond = _write(d, "d.tsv", "read1\tACC1|prot\t99\n")
            seqs = _write(d, "s.fasta", ">read1\nAAAA\n")
            taxids = _write(d, "t.tsv", "ACC1\t111\n")
            out = os.path.join(d, "out.tsv")

            result = subprocess.run(
                [
                    sys.executable,
                    SCRIPT,
                    "--diamond",
                    diamond,
                    "--sequences",
                    seqs,
                    "--taxids",
                    taxids,
                    "--output",
                    out,
                    "--data-format",
                    "fasta",
                ],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(os.path.exists(out))
            with open(out) as f:
                self.assertIn("read1\t111", f.read())


if __name__ == "__main__":
    unittest.main()
