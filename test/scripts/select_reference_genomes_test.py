"""Tests for viralunity.scripts.python.select_reference_genomes.

Covers the post-filter table selection added so reference assembly reads the
bleed / negative-control-filtered summaries instead of the raw counts table:

  - apply_pass_column_filters: drops only explicit bleed_pass/neg_pass == False
    rows (bool or string), keeps NA / True / missing columns.
  - resolve_summary_file: pick the most-filtered summary table (longest chain),
    degrading to less-filtered tables and finally the raw counts table; ignore
    audit sidecars (*.dropped.tsv / *_nr_flags.tsv).
  - collect_summary_files: one file per enabled (method, source) classifier.
"""

import os
import unittest

import numpy as np
import pandas as pd

from viralunity.scripts.python.select_reference_genomes import (
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


class TestResolveSummaryFile(unittest.TestCase):
    def _make(self, root, classifier, suffix):
        d = os.path.join(root, classifier)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, f"{classifier}_taxa_summary{suffix}.tsv")
        with open(path, "w") as fh:
            fh.write("sample\tmode\tcount\n")
        return path

    def test_prefers_most_filtered_chain(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self._make(root, "kraken2_reads", "")
            self._make(root, "kraken2_reads", "_RPM")
            self._make(root, "kraken2_reads", "_RPM.bleed")
            full = self._make(root, "kraken2_reads", "_RPM.bleed.neg.ictv")
            self.assertEqual(resolve_summary_file(root, "kraken2_reads"), full)

    def test_falls_back_to_less_filtered(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            raw = self._make(root, "kraken2_reads", "")
            self.assertEqual(resolve_summary_file(root, "kraken2_reads"), raw)
            bleed = self._make(root, "kraken2_reads", "_RPM.bleed")
            self.assertEqual(resolve_summary_file(root, "kraken2_reads"), bleed)

    def test_ignores_audit_sidecars(self):
        # A .dropped sidecar has the same chain depth as the real table and a
        # longer name, so it must be excluded explicitly, not out-ranked.
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self._make(root, "kraken2_contigs", "_RPKM.nr.dropped")
            real = self._make(root, "kraken2_contigs", "_RPKM.nr.bleed")
            self.assertEqual(resolve_summary_file(root, "kraken2_contigs"), real)

    def test_returns_none_when_absent(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "kraken2_reads"), exist_ok=True)
            self.assertIsNone(resolve_summary_file(root, "kraken2_reads"))


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
            files = collect_summary_files(root, "kraken2", "reads")
            self.assertEqual(
                [os.path.basename(f) for f in files], ["kraken2_reads_taxa_summary_RPM.bleed.tsv"]
            )

    def test_both_collects_all_present(self):
        import tempfile

        with tempfile.TemporaryDirectory() as root:
            self._seed(root)
            files = collect_summary_files(root, "both", "both")
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
