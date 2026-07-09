"""Tests for viralunity.scripts.python.annotate_nr_taxonomy."""

import contextlib
import io
import os
import tempfile
import unittest

from viralunity.scripts.python.annotate_nr_taxonomy import (
    RANK_ORDER,
    annotate,
    resolve_lineage_ranks,
)

# Synthetic viral lineage:
#   10239 superkingdom "Viruses"        (-> domain)
#   2559587 realm "Riboviria"           parent 10239
#   11308 family "Orthomyxoviridae"     parent 2559587
#   197911 genus "Alphainfluenzavirus"  parent 11308
#   11320 species "Influenza A virus"   parent 197911
NODES = [
    ("1", "1", "no rank"),
    ("10239", "1", "superkingdom"),
    ("2559587", "10239", "realm"),
    ("11308", "2559587", "family"),
    ("197911", "11308", "genus"),
    ("11320", "197911", "species"),
]
NAMES = [
    ("1", "root"),
    ("10239", "Viruses"),
    ("2559587", "Riboviria"),
    ("11308", "Orthomyxoviridae"),
    ("197911", "Alphainfluenzavirus"),
    ("11320", "Influenza A virus"),
]


def _write_taxdump(d):
    nodes = os.path.join(d, "nodes.dmp")
    names = os.path.join(d, "names.dmp")
    with open(nodes, "w") as f:
        for t, p, r in NODES:
            f.write(f"{t}\t|\t{p}\t|\t{r}\t|\t-\t|\n")
    with open(names, "w") as f:
        for t, n in NAMES:
            f.write(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n")
    return nodes, names


class TestResolveLineageRanks(unittest.TestCase):
    def setUp(self):
        from viralunity.scripts.python.taxonomy import load_taxdump

        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        nodes, names = _write_taxdump(self._tmp.name)
        self.parent, self.rank, self.name = load_taxdump(nodes, names)

    def test_species_resolves_full_available_lineage(self):
        ranks = resolve_lineage_ranks("11320", self.parent, self.rank, self.name)
        d = dict(zip(RANK_ORDER, ranks))
        self.assertEqual(d["domain"], "Viruses")  # superkingdom aliased to domain
        self.assertEqual(d["realm"], "Riboviria")
        self.assertEqual(d["family"], "Orthomyxoviridae")
        self.assertEqual(d["genus"], "Alphainfluenzavirus")
        self.assertEqual(d["species"], "Influenza A virus")
        # ranks absent from the lineage are NA
        self.assertEqual(d["kingdom"], "NA")
        self.assertEqual(d["phylum"], "NA")

    def test_unknown_taxid_is_all_na(self):
        ranks = resolve_lineage_ranks("999999", self.parent, self.rank, self.name)
        self.assertEqual(set(ranks), {"NA"})

    def test_multiple_staxids_uses_first(self):
        # helper takes the first ';'-separated taxid
        ranks = resolve_lineage_ranks("11320;12345", self.parent, self.rank, self.name)
        self.assertEqual(dict(zip(RANK_ORDER, ranks))["species"], "Influenza A virus")


class TestAnnotate(unittest.TestCase):
    def setUp(self):
        from viralunity.scripts.python.taxonomy import load_taxdump

        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        nodes, names = _write_taxdump(self._tmp.name)
        self.parent, self.rank, self.name = load_taxdump(nodes, names)
        self.in_path = os.path.join(self._tmp.name, "nr.tsv")
        self.out_path = os.path.join(self._tmp.name, "annotated.tsv")

    def test_appends_nine_rank_columns_dropping_staxids(self):
        # 12 blast columns + staxids (13th)
        blast = ["q1", "s1", "99", "100", "0", "0", "1", "100", "1", "100", "1e-40", "200"]
        with open(self.in_path, "w") as f:
            f.write("\t".join(blast + ["11320"]) + "\n")
        annotate(self.in_path, self.out_path, self.parent, self.rank, self.name)
        with open(self.out_path) as f:
            fields = f.readline().rstrip("\n").split("\t")
        # 12 blast + 9 ranks = 21 columns, no header
        self.assertEqual(len(fields), 21)
        self.assertEqual(fields[:12], blast)
        self.assertEqual(fields[20], "Influenza A virus")  # species column

    def test_rows_missing_staxids_are_skipped_with_a_warning(self):
        blast = ["q1", "s1", "99", "100", "0", "0", "1", "100", "1", "100", "1e-40", "200"]
        with open(self.in_path, "w") as f:
            f.write("\t".join(blast) + "\n")  # 12 cols: no staxids field
            f.write("\t".join(blast + ["11320"]) + "\n")  # well-formed
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            written = annotate(self.in_path, self.out_path, self.parent, self.rank, self.name)
        self.assertEqual(written, 1)  # only the well-formed row
        msg = err.getvalue()
        self.assertIn("skipped 1", msg)
        self.assertIn("staxids", msg)


if __name__ == "__main__":
    unittest.main()
