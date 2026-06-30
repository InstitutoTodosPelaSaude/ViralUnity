"""Tests for viralunity.scripts.python.select_reference_genomes.

Covers the post-filter table selection added so reference assembly reads the
bleed / negative-control-filtered summaries instead of the raw counts table:

  - apply_pass_column_filters: drops only explicit bleed_pass/neg_pass == False
    rows (bool or string), keeps NA / True / missing columns.
  - _fallback_suffixes / resolve_summary_file: prefer the requested suffix, then
    fall back through bleed-only to the raw counts table.
  - collect_summary_files: one file per enabled (method, source) classifier.
"""

import os
import unittest

import numpy as np
import pandas as pd

from viralunity.scripts.python.select_reference_genomes import (
    _fallback_suffixes,
    apply_pass_column_filters,
    collect_summary_files,
    resolve_summary_file,
)


class TestApplyPassColumnFilters(unittest.TestCase):
    def test_drops_explicit_false_bool(self):
        df = pd.DataFrame(
            {
                "name": ["keep", "drop"],
                "bleed_pass": [True, False],
            }
        )
        out = apply_pass_column_filters(df)
        self.assertEqual(list(out["name"]), ["keep"])

    def test_drops_explicit_false_string_and_numpy(self):
        df = pd.DataFrame(
            {
                "name": ["a", "b", "c"],
                "neg_pass": ["True", "False", np.bool_(False)],
            }
        )
        out = apply_pass_column_filters(df)
        self.assertEqual(list(out["name"]), ["a"])

    def test_keeps_na_rows(self):
        # NA neg_pass (no negative controls / taxon absent from controls) is kept.
        df = pd.DataFrame(
            {
                "name": ["x", "y"],
                "bleed_pass": [True, True],
                "neg_pass": [np.nan, "NA"],
            }
        )
        out = apply_pass_column_filters(df)
        # "NA" is a plain string here (not parsed as NaN), but it is not "false".
        self.assertEqual(list(out["name"]), ["x", "y"])

    def test_combined_bleed_and_neg(self):
        df = pd.DataFrame(
            {
                "name": ["both_pass", "bleed_fail", "neg_fail"],
                "bleed_pass": [True, False, True],
                "neg_pass": [True, True, False],
            }
        )
        out = apply_pass_column_filters(df)
        self.assertEqual(list(out["name"]), ["both_pass"])

    def test_missing_columns_passthrough(self):
        df = pd.DataFrame({"name": ["a", "b"], "count": [10, 20]})
        out = apply_pass_column_filters(df)
        self.assertEqual(list(out["name"]), ["a", "b"])


class TestFallbackSuffixes(unittest.TestCase):
    def test_neg_order(self):
        self.assertEqual(
            _fallback_suffixes("_RPM.bleed.neg"),
            ["_RPM.bleed.neg", "_RPM.bleed", ""],
        )

    def test_bleed_dedup(self):
        self.assertEqual(_fallback_suffixes("_RPM.bleed"), ["_RPM.bleed", ""])

    def test_empty_dedup(self):
        self.assertEqual(_fallback_suffixes(""), [""])


class TestResolveSummaryFile(unittest.TestCase):
    def _make(self, root, classifier, suffix):
        d = os.path.join(root, classifier)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, f"{classifier}_taxa_summary{suffix}.tsv")
        with open(path, "w") as fh:
            fh.write("sample\tmode\tcount\n")
        return path

    def test_prefers_requested_suffix(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self._make(root, "kraken2_reads", "")
            self._make(root, "kraken2_reads", "_RPM.bleed")
            neg = self._make(root, "kraken2_reads", "_RPM.bleed.neg")
            got = resolve_summary_file(root, "kraken2_reads", "_RPM.bleed.neg")
            self.assertEqual(got, neg)

    def test_falls_back_to_bleed_then_raw(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            raw = self._make(root, "kraken2_reads", "")
            # neg requested but only raw exists -> raw (bleed missing too)
            got = resolve_summary_file(root, "kraken2_reads", "_RPM.bleed.neg")
            self.assertEqual(got, raw)

            bleed = self._make(root, "kraken2_reads", "_RPM.bleed")
            got = resolve_summary_file(root, "kraken2_reads", "_RPM.bleed.neg")
            self.assertEqual(got, bleed)

    def test_returns_none_when_absent(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self.assertIsNone(resolve_summary_file(root, "kraken2_reads", "_RPM.bleed"))


class TestCollectSummaryFiles(unittest.TestCase):
    def _seed(self, root):
        for classifier in ("kraken2_reads", "kraken2_contigs", "diamond_contigs"):
            d = os.path.join(root, classifier)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, f"{classifier}_taxa_summary_RPM.bleed.tsv"), "w") as fh:
                fh.write("sample\tmode\tcount\n")

    def test_method_source_selection(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self._seed(root)
            files = collect_summary_files(root, "kraken2", "reads", "_RPM.bleed")
            self.assertEqual(
                [os.path.basename(f) for f in files], ["kraken2_reads_taxa_summary_RPM.bleed.tsv"]
            )

    def test_both_collects_all_present(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self._seed(root)
            files = collect_summary_files(root, "both", "both", "_RPM.bleed")
            names = sorted(os.path.basename(f) for f in files)
            # diamond_reads not seeded -> not present
            self.assertEqual(
                names,
                [
                    "diamond_contigs_taxa_summary_RPM.bleed.tsv",
                    "kraken2_contigs_taxa_summary_RPM.bleed.tsv",
                    "kraken2_reads_taxa_summary_RPM.bleed.tsv",
                ],
            )


if __name__ == "__main__":
    unittest.main()
