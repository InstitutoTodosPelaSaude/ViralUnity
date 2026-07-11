"""Tests for viralunity.scripts.python.extract_viral_contigs."""

import os
import tempfile
import unittest

from viralunity.scripts.python.extract_viral_contigs import (
    viral_contig_ids,
    write_viral_contigs,
)

# cA -> species 3001 -> family 1000 -> Viruses 10239; cB -> Bacteria 2.
PARENT = {"3001": "1000", "1000": "10239", "10239": "1", "5001": "2", "2": "1", "1": "1"}


def _write_krona(path, rows):
    with open(path, "w") as f:
        for contig, taxid in rows:
            f.write(f"{contig}\t{taxid}\n")


class TestViralContigIds(unittest.TestCase):
    def test_keeps_only_viral_lineage(self):
        with tempfile.TemporaryDirectory() as d:
            krona = os.path.join(d, "k.txt")
            _write_krona(krona, [("cA", "3001"), ("cB", "5001"), ("cC", "0")])
            self.assertEqual(viral_contig_ids(krona, PARENT), {"cA"})

    def test_missing_file_is_empty(self):
        self.assertEqual(viral_contig_ids("/no/such/file", PARENT), set())


class TestWriteViralContigs(unittest.TestCase):
    def _taxdump(self, d):
        nodes = os.path.join(d, "nodes.dmp")
        with open(nodes, "w") as f:
            f.write("3001\t|\t1000\t|\tspecies\t|\n")
            f.write("1000\t|\t10239\t|\tfamily\t|\n")
            f.write("10239\t|\t1\t|\tsuperkingdom\t|\n")
            f.write("5001\t|\t2\t|\tspecies\t|\n")
            f.write("2\t|\t1\t|\tsuperkingdom\t|\n")
            f.write("1\t|\t1\t|\tno rank\t|\n")
        return d

    def test_writes_only_viral_records(self):
        with tempfile.TemporaryDirectory() as d:
            self._taxdump(d)
            contigs = os.path.join(d, "contigs.fa")
            with open(contigs, "w") as f:
                f.write(">cA flag=1 len=10\nACGTACGTAC\n")
                f.write(">cB flag=1 len=8\nTTTTTTTT\n")
            krona = os.path.join(d, "k.txt")
            _write_krona(krona, [("cA", "3001"), ("cB", "5001")])
            out_fa = os.path.join(d, "viral.fa")
            out_ids = os.path.join(d, "viral.ids.txt")

            kept = write_viral_contigs(contigs, krona, d, out_fa, out_ids)
            self.assertEqual(kept, 1)
            self.assertEqual(open(out_ids).read().split(), ["cA"])
            fa = open(out_fa).read()
            self.assertIn(">cA", fa)
            self.assertNotIn(">cB", fa)
            self.assertIn("ACGTACGTAC", fa)

    def test_empty_inputs_write_empty_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            self._taxdump(d)
            contigs = os.path.join(d, "contigs.fa")
            open(contigs, "w").close()
            krona = os.path.join(d, "k.txt")
            open(krona, "w").close()
            out_fa = os.path.join(d, "viral.fa")
            out_ids = os.path.join(d, "viral.ids.txt")
            kept = write_viral_contigs(contigs, krona, d, out_fa, out_ids)
            self.assertEqual(kept, 0)
            self.assertEqual(os.path.getsize(out_fa), 0)


if __name__ == "__main__":
    unittest.main()
