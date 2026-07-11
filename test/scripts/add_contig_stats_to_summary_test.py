"""Tests for viralunity.scripts.python.add_contig_stats_to_summary."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.add_contig_stats_to_summary import (
    add_contig_stats,
    build_taxon_best_contig,
    load_contig_depths,
)

# species 3001 under family 1000; species 3002 under the same family.
PARENT = {"3001": "1000", "3002": "1000", "1000": "1", "1": "1"}
RANK = {"3001": "species", "3002": "species", "1000": "family", "1": "no rank"}


class TestLoadContigDepths(unittest.TestCase):
    def test_length_and_median(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "depth.txt")
            with open(p, "w") as f:
                # contigA: 5 positions, depths 1..5 -> median 3, length 5
                for pos, dep in enumerate([1, 2, 3, 4, 5], start=1):
                    f.write(f"contigA\t{pos}\t{dep}\n")
                # contigB: 2 positions, depths 10, 20 -> median 15, length 2
                f.write("contigB\t1\t10\n")
                f.write("contigB\t2\t20\n")
            stats = load_contig_depths(p)
            self.assertEqual(stats["contigA"], (5, 3.0))
            self.assertEqual(stats["contigB"], (2, 15.0))

    def test_missing_file_returns_empty(self):
        self.assertEqual(load_contig_depths("/does/not/exist"), {})


class TestBuildTaxonBestContig(unittest.TestCase):
    def _krona(self, d, rows):
        p = os.path.join(d, "kr.tsv")
        with open(p, "w") as f:
            for contig, taxid in rows:
                f.write(f"{contig}\t{taxid}\n")
        return p

    def test_picks_largest_contig_and_propagates_to_family(self):
        with tempfile.TemporaryDirectory() as d:
            krona = self._krona(d, [("cA", "3001"), ("cB", "3001")])
            depth_by_sample = {"S1": {"cA": (500, 10.0), "cB": (5000, 100.0)}}  # cB is larger
            best = build_taxon_best_contig([("S1", krona)], depth_by_sample, PARENT, RANK)
            # species 3001: largest is cB (5000 bp, depth 100)
            self.assertEqual(best[("S1", "species", "3001")], (5000, 100.0))
            # family 1000 inherits the largest contig of any descendant species
            self.assertEqual(best[("S1", "family", "1000")], (5000, 100.0))

    def test_family_takes_max_across_species(self):
        with tempfile.TemporaryDirectory() as d:
            krona = self._krona(d, [("cA", "3001"), ("cB", "3002")])
            depth_by_sample = {"S1": {"cA": (1000, 5.0), "cB": (8000, 9.0)}}
            best = build_taxon_best_contig([("S1", krona)], depth_by_sample, PARENT, RANK)
            self.assertEqual(best[("S1", "species", "3001")], (1000, 5.0))
            self.assertEqual(best[("S1", "species", "3002")], (8000, 9.0))
            # family 1000: max over both species -> cB
            self.assertEqual(best[("S1", "family", "1000")], (8000, 9.0))

    def test_contig_without_depth_is_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            krona = self._krona(d, [("cA", "3001")])
            best = build_taxon_best_contig([("S1", krona)], {"S1": {}}, PARENT, RANK)
            self.assertNotIn(("S1", "species", "3001"), best)


class TestAddContigStats(unittest.TestCase):
    def test_columns_added_with_na_for_missing(self):
        summary = pd.DataFrame(
            [
                {"sample": "S1", "rank": "species", "taxid": "3001", "name": "Virus A"},
                {"sample": "S1", "rank": "family", "taxid": "1000", "name": "Familyidae"},
                {"sample": "S1", "rank": "species", "taxid": "9999", "name": "Absent"},
            ]
        )
        best = {
            ("S1", "species", "3001"): (5000, 100.0),
            ("S1", "family", "1000"): (5000, 100.0),
        }
        out = add_contig_stats(summary, best)
        self.assertIn("largest_contig_bp", out.columns)
        self.assertIn("largest_contig_median_depth", out.columns)
        r3001 = out[out["taxid"] == "3001"].iloc[0]
        self.assertEqual(r3001["largest_contig_bp"], 5000)
        self.assertEqual(r3001["largest_contig_median_depth"], 100.0)
        # taxon with no assigned contig -> NA
        r9999 = out[out["taxid"] == "9999"].iloc[0]
        self.assertEqual(r9999["largest_contig_bp"], "NA")
        self.assertEqual(r9999["largest_contig_median_depth"], "NA")

    def test_missing_required_column_raises(self):
        summary = pd.DataFrame([{"sample": "S1", "rank": "species"}])  # no taxid
        with self.assertRaises(ValueError):
            add_contig_stats(summary, {})


if __name__ == "__main__":
    unittest.main()
