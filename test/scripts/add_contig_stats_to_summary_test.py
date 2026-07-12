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
        self.assertIn("largest_contig_ref_coverage_pct", out.columns)
        r3001 = out[out["taxid"] == "3001"].iloc[0]
        self.assertEqual(r3001["largest_contig_bp"], 5000)
        self.assertEqual(r3001["largest_contig_median_depth"], 100.0)
        # taxon with no assigned contig -> NA
        r9999 = out[out["taxid"] == "9999"].iloc[0]
        self.assertEqual(r9999["largest_contig_bp"], "NA")
        self.assertEqual(r9999["largest_contig_median_depth"], "NA")
        # no genome_length_bp column -> coverage is all-NA (backward compatible)
        self.assertTrue((out["largest_contig_ref_coverage_pct"] == "NA").all())

    def test_missing_required_column_raises(self):
        summary = pd.DataFrame([{"sample": "S1", "rank": "species"}])  # no taxid
        with self.assertRaises(ValueError):
            add_contig_stats(summary, {})


class TestRefCoveragePct(unittest.TestCase):
    """largest_contig_ref_coverage_pct = largest_contig_bp / genome_length_bp * 100,
    raw and uncapped; NA when there is no contig or no usable reference length."""

    def test_normal_uncapped_and_na_cases(self):
        summary = pd.DataFrame(
            [
                # 7000 / 10000 -> 70.0
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3001",
                    "name": "A",
                    "genome_length_bp": "10000",
                },
                # 11200 / 10000 -> 112.0 (uncapped: contig exceeds the reference)
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3002",
                    "name": "B",
                    "genome_length_bp": "10000",
                },
                # genome_length_bp = NA -> NA
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3003",
                    "name": "C",
                    "genome_length_bp": "NA",
                },
                # genome_length_bp = 0 -> NA (guards divide-by-zero)
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3004",
                    "name": "D",
                    "genome_length_bp": "0",
                },
                # length present but no assigned contig -> NA
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3005",
                    "name": "E",
                    "genome_length_bp": "10000",
                },
            ]
        )
        best = {
            ("S1", "species", "3001"): (7000, 50.0),
            ("S1", "species", "3002"): (11200, 60.0),
            ("S1", "species", "3003"): (5000, 10.0),
            ("S1", "species", "3004"): (5000, 10.0),
            # 3005 intentionally absent from best
        }
        by = add_contig_stats(summary, best).set_index("taxid")["largest_contig_ref_coverage_pct"]
        self.assertEqual(by["3001"], 70.0)
        self.assertEqual(by["3002"], 112.0)  # uncapped
        self.assertEqual(by["3003"], "NA")  # reference length NA
        self.assertEqual(by["3004"], "NA")  # reference length 0
        self.assertEqual(by["3005"], "NA")  # no assigned contig

    def test_rounding_to_two_dp(self):
        summary = pd.DataFrame(
            [
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3001",
                    "name": "A",
                    "genome_length_bp": "29903",
                }
            ]
        )
        best = {("S1", "species", "3001"): (2999, 10.0)}
        out = add_contig_stats(summary, best)
        # 2999 / 29903 * 100 = 10.0291... -> 10.03
        self.assertEqual(out.iloc[0]["largest_contig_ref_coverage_pct"], 10.03)

    def test_column_positioned_between_bp_and_depth(self):
        summary = pd.DataFrame(
            [
                {
                    "sample": "S1",
                    "rank": "species",
                    "taxid": "3001",
                    "name": "A",
                    "genome_length_bp": "10000",
                }
            ]
        )
        best = {("S1", "species", "3001"): (7000, 50.0)}
        cols = list(add_contig_stats(summary, best).columns)
        i = cols.index("largest_contig_bp")
        self.assertEqual(cols[i + 1], "largest_contig_ref_coverage_pct")
        self.assertEqual(cols[i + 2], "largest_contig_median_depth")


if __name__ == "__main__":
    unittest.main()
