"""Tests for viralunity.scripts.python.filter_top_species_hits (ported from the
REVISA nr_validation prototype)."""

import contextlib
import io
import os
import tempfile
import unittest

from viralunity.scripts.python.filter_top_species_hits import (
    choose_representative,
    is_viral,
    lca_consensus,
    process_file,
)

# 21-col row = 12 blast cols + 9 ranks (domain..species). Helper builds one.
BLAST = ["q", "s", "99", "100", "0", "0", "1", "100", "1", "100", "1e-40", "200"]


def _row(qseqid, domain="NA", realm="NA", phylum="NA", family="NA", genus="NA", species="NA"):
    blast = list(BLAST)
    blast[0] = qseqid
    # order: domain kingdom realm phylum class order family genus species
    ranks = [domain, "NA", realm, phylum, "NA", "NA", family, genus, species]
    return blast + ranks


class TestChooseRepresentative(unittest.TestCase):
    def test_picks_first_hit_with_species(self):
        hits = [_row("q", species="NA"), _row("q", species="Influenza A virus")]
        self.assertEqual(choose_representative(hits)[-1], "Influenza A virus")

    def test_falls_back_to_best_hit_when_no_species(self):
        hits = [_row("q", family="Fam1"), _row("q", family="Fam2")]
        self.assertEqual(choose_representative(hits)[-3], "Fam1")  # family col of first


class TestIsViral(unittest.TestCase):
    def test_viral_phylum_suffix(self):
        self.assertTrue(is_viral(_row("q", phylum="Negarnaviricota")))

    def test_non_viral_phylum(self):
        self.assertFalse(is_viral(_row("q", phylum="Pseudomonadota")))

    def test_na_phylum_not_viral(self):
        self.assertFalse(is_viral(_row("q", phylum="NA")))


class TestLcaConsensus(unittest.TestCase):
    def test_majority_species_wins(self):
        hits = [
            _row("q", species="Influenza A virus"),
            _row("q", species="Influenza A virus"),
            _row("q", species="Other virus"),
        ]
        rank, taxon = lca_consensus(hits, 0.5)
        self.assertEqual((rank, taxon), ("species", "Influenza A virus"))

    def test_climbs_to_family_when_species_below_threshold(self):
        hits = [
            _row("q", family="Orthomyxoviridae", species="A"),
            _row("q", family="Orthomyxoviridae", species="B"),
            _row("q", family="Orthomyxoviridae", species="C"),
        ]
        rank, taxon = lca_consensus(hits, 0.5)
        self.assertEqual((rank, taxon), ("family", "Orthomyxoviridae"))

    def test_na_counts_against_consensus(self):
        # 1 of 3 hits has species -> 0.33 < 0.5 -> not species-level
        hits = [_row("q", species="X"), _row("q", species="NA"), _row("q", species="NA")]
        rank, _ = lca_consensus(hits, 0.5)
        self.assertNotEqual(rank, "species")

    def test_no_consensus_is_unclassified(self):
        hits = [_row("q"), _row("q")]  # all NA
        self.assertEqual(lca_consensus(hits, 0.5), ("unclassified", "NA"))


class TestProcessFile(unittest.TestCase):
    def test_one_row_per_query_plus_viral_subset(self):
        with tempfile.TemporaryDirectory() as d:
            in_path = os.path.join(d, "in.tsv")
            out_path = os.path.join(d, "out.tsv")
            viral_path = os.path.join(d, "out.viruses_only.tsv")
            with open(in_path, "w") as f:
                # q1: viral (phylum viricota); q2: bacterial
                f.write("\t".join(_row("q1", phylum="Negarnaviricota", species="Flu")) + "\n")
                f.write("\t".join(_row("q1", phylum="Negarnaviricota", species="Flu")) + "\n")
                f.write("\t".join(_row("q2", phylum="Pseudomonadota", species="E. coli")) + "\n")
            process_file(in_path, out_path, viral_path, threshold=0.5)
            with open(out_path) as f:
                out_lines = [ln for ln in f.read().splitlines() if ln.strip()]
            with open(viral_path) as f:
                viral_lines = [ln for ln in f.read().splitlines() if ln.strip()]
            # header + 2 query rows
            self.assertEqual(len(out_lines), 3)
            # header + 1 viral query row
            self.assertEqual(len(viral_lines), 2)
            self.assertIn("q1", viral_lines[1])

    def test_wrong_column_count_rows_skipped_with_warning(self):
        with tempfile.TemporaryDirectory() as d:
            in_path = os.path.join(d, "in.tsv")
            out_path = os.path.join(d, "out.tsv")
            viral_path = os.path.join(d, "out.viruses_only.tsv")
            with open(in_path, "w") as f:
                f.write("\t".join(_row("q1", species="Flu")) + "\n")  # 21 cols, valid
                f.write("q_bad\tonly\tthree\n")  # malformed: not 21 columns
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                queries, _viral = process_file(in_path, out_path, viral_path, threshold=0.5)
            self.assertEqual(queries, 1)  # only the valid query
            msg = err.getvalue()
            self.assertIn("skipped 1", msg)


if __name__ == "__main__":
    unittest.main()
