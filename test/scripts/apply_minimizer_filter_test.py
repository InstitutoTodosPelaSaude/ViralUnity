"""Tests for viralunity.scripts.python.apply_minimizer_filter."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.apply_minimizer_filter import (
    filter_summary,
    parse_report_minimizers,
    passes_minimizers,
    run,
)


def _write_report(path, lines):
    """lines: list of (pct, clade_reads, taxon_reads, n_min, n_distinct, rank, taxid, name)."""
    with open(path, "w") as f:
        for row in lines:
            f.write("\t".join(str(x) for x in row) + "\n")


class TestParseReportMinimizers(unittest.TestCase):
    def test_parses_taxid_to_minimizer_counts(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "s.report.txt")
            _write_report(
                p,
                [
                    ("50.0", "1000", "0", "5000", "4000", "R", "1", "root"),
                    ("40.0", "800", "50", "3000", "2500", "F", "1000", "  Orthomyxoviridae"),
                    ("2.0", "40", "40", "60", "10", "S", "3100", "    Lambdavirus DE3"),
                ],
            )
            counts = parse_report_minimizers(p)
            self.assertEqual(counts["1000"], (3000, 2500))
            self.assertEqual(counts["3100"], (60, 10))

    def test_skips_lines_without_minimizer_columns(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "s.report.txt")
            # A plain 6-column kraken2 line (no --report-minimizer-data) must be skipped.
            with open(p, "w") as f:
                f.write("40.0\t800\t50\tF\t1000\tOrthomyxoviridae\n")
            self.assertEqual(parse_report_minimizers(p), {})


class TestPassesMinimizers(unittest.TestCase):
    def test_distinct_below_threshold_fails(self):
        self.assertFalse(passes_minimizers(10, 60, min_distinct=50, max_dup=None))

    def test_distinct_at_threshold_passes(self):
        self.assertTrue(passes_minimizers(50, 60, min_distinct=50, max_dup=None))

    def test_duplication_above_max_fails(self):
        # total/distinct = 100/10 = 10 > 5
        self.assertFalse(passes_minimizers(10, 100, min_distinct=None, max_dup=5.0))

    def test_duplication_within_max_passes(self):
        self.assertTrue(passes_minimizers(50, 100, min_distinct=None, max_dup=5.0))

    def test_no_thresholds_always_passes(self):
        self.assertTrue(passes_minimizers(1, 999999, min_distinct=None, max_dup=None))


class TestFilterSummary(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.in_path = os.path.join(self._tmp.name, "in.tsv")
        self.out_path = os.path.join(self._tmp.name, "out.tsv")
        self.dropped_path = os.path.join(self._tmp.name, "dropped.tsv")
        df = pd.DataFrame(
            [
                ("A", "kraken2", "reads", "species", "3001"),  # high distinct -> keep
                ("A", "kraken2", "reads", "species", "3100"),  # low distinct -> drop
                ("A", "kraken2", "reads", "species", "5555"),  # no report data -> keep
            ],
            columns=["sample", "tool", "mode", "rank", "taxid"],
        )
        df.to_csv(self.in_path, sep="\t", index=False)
        self.counts = {"A": {"3001": (4000, 3000), "3100": (60, 10)}}

    def test_drops_only_low_distinct_rows(self):
        kept, dropped = filter_summary(
            self.in_path,
            self.out_path,
            self.dropped_path,
            self.counts,
            min_distinct=50,
            max_dup=None,
        )
        self.assertEqual(kept, 2)
        self.assertEqual(dropped, 1)
        out = pd.read_csv(self.out_path, sep="\t", dtype=str)
        self.assertEqual(set(out["taxid"]), {"3001", "5555"})

    def test_missing_report_data_is_kept(self):
        filter_summary(
            self.in_path,
            self.out_path,
            self.dropped_path,
            self.counts,
            min_distinct=50,
            max_dup=None,
        )
        dropped = pd.read_csv(self.dropped_path, sep="\t", dtype=str)
        self.assertEqual(set(dropped["taxid"]), {"3100"})

    def test_empty_input_yields_empty_outputs(self):
        open(self.in_path, "w").close()
        kept, dropped = filter_summary(
            self.in_path, self.out_path, self.dropped_path, self.counts, 50, None
        )
        self.assertEqual((kept, dropped), (0, 0))
        self.assertTrue(os.path.exists(self.out_path))


class TestRunEndToEnd(unittest.TestCase):
    def test_run_maps_reports_to_samples(self):
        with tempfile.TemporaryDirectory() as d:
            report_a = os.path.join(d, "A.report.txt")
            _write_report(
                report_a,
                [
                    ("2.0", "40", "40", "4000", "3000", "S", "3001", "  Keep sp"),
                    ("1.0", "40", "40", "60", "10", "S", "3100", "  Drop sp"),
                ],
            )
            in_path = os.path.join(d, "in.tsv")
            out_path = os.path.join(d, "out.tsv")
            dropped_path = os.path.join(d, "dropped.tsv")
            pd.DataFrame(
                [
                    ("A", "kraken2", "reads", "species", "3001"),
                    ("A", "kraken2", "reads", "species", "3100"),
                ],
                columns=["sample", "tool", "mode", "rank", "taxid"],
            ).to_csv(in_path, sep="\t", index=False)
            run(
                summary=in_path,
                output=out_path,
                dropped=dropped_path,
                reports=[report_a],
                samples=["A"],
                min_distinct=50,
                max_dup=None,
            )
            out = pd.read_csv(out_path, sep="\t", dtype=str)
            self.assertEqual(set(out["taxid"]), {"3001"})


if __name__ == "__main__":
    unittest.main()
