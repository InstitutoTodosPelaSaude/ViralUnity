"""Tests for viralunity.scripts.python.build_ictv_vertebrate_taxids."""

import os
import tempfile
import unittest

from viralunity.scripts.python.build_ictv_vertebrate_taxids import (
    build_name_index,
    extract_vertebrate_taxa,
    is_vertebrate_host,
    resolve_names_to_taxids,
)


class TestIsVertebrateHost(unittest.TestCase):
    def test_plain_vertebrates(self):
        self.assertTrue(is_vertebrate_host("vertebrates"))

    def test_case_insensitive(self):
        self.assertTrue(is_vertebrate_host("Vertebrates"))

    def test_combined_host_including_vertebrates(self):
        self.assertTrue(is_vertebrate_host("invertebrates, vertebrates"))

    def test_invertebrates_only_is_not_vertebrate(self):
        # The substring trap: 'invertebrates' contains 'vertebrates'.
        self.assertFalse(is_vertebrate_host("invertebrates"))

    def test_bacteria_host_is_not_vertebrate(self):
        self.assertFalse(is_vertebrate_host("bacteria"))

    def test_empty_or_nan(self):
        self.assertFalse(is_vertebrate_host(""))
        self.assertFalse(is_vertebrate_host(None))


class TestExtractVertebrateTaxa(unittest.TestCase):
    def test_collects_family_and_genus_names_for_vertebrate_rows(self):
        rows = [
            {
                "Family": "Orthomyxoviridae",
                "Genus": "Alphainfluenzavirus",
                "Host source": "vertebrates",
            },
            {"Family": "Steitzviridae", "Genus": "Lambdavirus", "Host source": "bacteria"},
            {
                "Family": "Flaviviridae",
                "Genus": "Flavivirus",
                "Host source": "invertebrates, vertebrates",
            },
            {"Family": "Baculoviridae", "Genus": "Betabaculovirus", "Host source": "invertebrates"},
        ]
        names = extract_vertebrate_taxa(rows)
        self.assertEqual(
            names,
            {"Orthomyxoviridae", "Alphainfluenzavirus", "Flaviviridae", "Flavivirus"},
        )

    def test_ignores_blank_family_or_genus(self):
        rows = [
            {"Family": "Orthomyxoviridae", "Genus": "", "Host source": "vertebrates"},
            {"Family": "", "Genus": "", "Host source": "vertebrates"},
        ]
        self.assertEqual(extract_vertebrate_taxa(rows), {"Orthomyxoviridae"})


class TestBuildNameIndex(unittest.TestCase):
    def test_maps_lowercased_scientific_names_to_taxids(self):
        with tempfile.TemporaryDirectory() as d:
            names = os.path.join(d, "names.dmp")
            with open(names, "w") as f:
                f.write("1000\t|\tOrthomyxoviridae\t|\t\t|\tscientific name\t|\n")
                f.write("2000\t|\tAlphainfluenzavirus\t|\t\t|\tscientific name\t|\n")
                f.write("2000\t|\tFluvirus\t|\t\t|\tsynonym\t|\n")
            index = build_name_index(names)
            self.assertEqual(index["orthomyxoviridae"], "1000")
            self.assertEqual(index["alphainfluenzavirus"], "2000")
            # synonyms are not indexed (scientific name only)
            self.assertNotIn("fluvirus", index)


class TestResolveNamesToTaxids(unittest.TestCase):
    def test_resolves_known_names_and_reports_unresolved(self):
        index = {"orthomyxoviridae": "1000", "alphainfluenzavirus": "2000"}
        taxids, unresolved = resolve_names_to_taxids(
            {"Orthomyxoviridae", "Alphainfluenzavirus", "Madeupviridae"}, index
        )
        self.assertEqual(taxids, {"1000", "2000"})
        self.assertEqual(unresolved, {"Madeupviridae"})


if __name__ == "__main__":
    unittest.main()
