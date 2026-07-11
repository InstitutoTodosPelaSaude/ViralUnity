"""Tests for viralunity.scripts.python.split_summary_by_rank."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.split_summary_by_rank import (
    add_higher_rank_names,
    ensure_final_species,
    split_by_rank,
)

# species 3001 -> genus 2001 -> family 1000 -> root
PARENT = {"3001": "2001", "2001": "1000", "1000": "1", "1": "1"}
RANK = {"3001": "species", "2001": "genus", "1000": "family", "1": "no rank"}
NAME = {"3001": "Virus A", "2001": "Genusvirus", "1000": "Familyidae"}


class TestEnsureFinalSpecies(unittest.TestCase):
    def test_added_from_name_when_no_nr(self):
        df = pd.DataFrame([{"taxid": "3001", "name": "Virus A", "rank": "species"}])
        out = ensure_final_species(df)
        self.assertEqual(out["final_species"].iloc[0], "Virus A")
        # positioned right after name
        cols = list(out.columns)
        self.assertEqual(cols[cols.index("name") + 1], "final_species")

    def test_coalesces_with_nr_correction(self):
        df = pd.DataFrame(
            [
                {"taxid": "3001", "name": "Virus A", "nr_correct_species": "Virus B"},
                {"taxid": "3002", "name": "Virus C", "nr_correct_species": "NA"},
            ]
        )
        out = ensure_final_species(df)
        self.assertEqual(list(out["final_species"]), ["Virus B", "Virus C"])

    def test_left_alone_when_present(self):
        df = pd.DataFrame([{"taxid": "3001", "name": "Virus A", "final_species": "kept"}])
        out = ensure_final_species(df)
        self.assertEqual(out["final_species"].iloc[0], "kept")


class TestAddHigherRankNames(unittest.TestCase):
    def test_species_row_gets_family_and_genus(self):
        df = pd.DataFrame([{"taxid": "3001", "name": "Virus A", "rank": "species"}])
        out = add_higher_rank_names(df, PARENT, RANK, NAME)
        self.assertEqual(out["family"].iloc[0], "Familyidae")
        self.assertEqual(out["genus"].iloc[0], "Genusvirus")
        cols = list(out.columns)
        self.assertEqual(cols[cols.index("name") + 1], "family")
        self.assertEqual(cols[cols.index("name") + 2], "genus")


class TestSplitByRank(unittest.TestCase):
    def test_writes_three_files_with_propagated_columns(self):
        df = pd.DataFrame(
            [
                {"sample": "S1", "rank": "species", "taxid": "3001", "name": "Virus A"},
                {"sample": "S1", "rank": "genus", "taxid": "2001", "name": "Genusvirus"},
                {"sample": "S1", "rank": "family", "taxid": "1000", "name": "Familyidae"},
            ]
        )
        with tempfile.TemporaryDirectory() as d:
            paths = {r: os.path.join(d, r, f"{r}.tsv") for r in ("family", "genus", "species")}
            split_by_rank(df, paths, PARENT, RANK, NAME)

            sp = pd.read_csv(paths["species"], sep="\t")
            self.assertIn("family", sp.columns)
            self.assertIn("genus", sp.columns)
            self.assertEqual(sp["family"].iloc[0], "Familyidae")
            self.assertEqual(sp["genus"].iloc[0], "Genusvirus")
            self.assertIn("final_species", sp.columns)

            gn = pd.read_csv(paths["genus"], sep="\t")
            self.assertIn("family", gn.columns)
            self.assertNotIn("genus", gn.columns)  # own-rank column dropped

            fam = pd.read_csv(paths["family"], sep="\t")
            self.assertNotIn("family", fam.columns)
            self.assertNotIn("genus", fam.columns)

    def test_missing_rank_column_raises(self):
        df = pd.DataFrame([{"sample": "S1", "taxid": "3001", "name": "x"}])
        with self.assertRaises(ValueError):
            split_by_rank(df, {"species": "/tmp/x.tsv"}, PARENT, RANK, NAME)


if __name__ == "__main__":
    unittest.main()
