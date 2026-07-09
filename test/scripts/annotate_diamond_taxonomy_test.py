"""Tests for viralunity.scripts.python.annotate_diamond_taxonomy."""

import os
import tempfile
import unittest

from viralunity.scripts.python.annotate_diamond_taxonomy import run


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)
    return path


def _setup(d, diamond_content):
    diamond = _write(os.path.join(d, "diamond.tsv"), diamond_content)
    taxids = _write(os.path.join(d, "protein2taxid.tsv"), "ACC1\t111\n")
    taxdump = os.path.join(d, "taxdump")
    os.makedirs(taxdump)
    _write(os.path.join(taxdump, "names.dmp"), "111\t|\tVirus A\t|\t\t|\tscientific name\t|\n")
    out = os.path.join(d, "out.tsv")
    return diamond, taxids, taxdump, out


class TestAnnotateDiamondTaxonomy(unittest.TestCase):
    def test_happy_path(self):
        with tempfile.TemporaryDirectory() as d:
            diamond, taxids, taxdump, out = _setup(d, "contig1\tACC1|prot\t10\n")
            run(diamond, taxids, taxdump, out)
            with open(out) as f:
                self.assertEqual(f.read().strip(), "contig1\t111\tVirus A\t10")

    def test_rejects_non_integer_last_column(self):
        # A raw DIAMOND outfmt6 line ends in a float bitscore, not an integer
        # mapped-read count. Writing it verbatim as mapped_reads corrupts counts
        # silently downstream, so the mismatch must fail loudly instead.
        with tempfile.TemporaryDirectory() as d:
            diamond, taxids, taxdump, out = _setup(d, "contig1\tACC1|prot\t99.5\n")
            with self.assertRaises(ValueError):
                run(diamond, taxids, taxdump, out)


if __name__ == "__main__":
    unittest.main()
