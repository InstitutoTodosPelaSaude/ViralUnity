"""Tests for viralunity.scripts.python.generate_consensus_report.

Covers the pure, I/O-light helpers that back the interactive consensus report:
  - is_segmented (segment-column detection)
  - dedupe_and_sum_reads (per-sample dedupe of read counts, sum of mapped reads)
  - downsample_min_pool (min-pooling that preserves coverage dips)
  - resolve_basewise_path (exact per-base coverage path templates)
  - read_stats_summary / build_stats_table_html (percentage x100 + rounding)
  - load_basewise_table (missing file -> empty frame, not a raise)
  - the figure builders (return a Figure; depth chart carries 20x/100x guides)
"""

import io
import os
import tempfile
import unittest

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from viralunity.scripts.python.generate_consensus_report import (
    build_aggregated_coverage_line_plot,
    build_reads_histogram,
    build_report_metadata,
    build_stats_table_html,
    dedupe_and_sum_reads,
    downsample_min_pool,
    is_segmented,
    load_basewise_table,
    read_stats_summary,
    resolve_basewise_path,
)

UNSEG_CSV = (
    "sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,"
    "average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,"
    "horizontal_coverage\n"
    "sample-A,348018,330568,643219,3245.2445239608064,0.9959535832525165,"
    "0.9959535832525165,0.9475972310470522,0.9959535832525165\n"
    "sample-B,1000,900,800,45.5,0.5,0.25,0.0,0.6\n"
)

SEG_CSV = (
    "sample_name,segment,number_of_reads,number_of_trim_paired_reads,"
    "number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,"
    "percentage_above_1000x,horizontal_coverage\n"
    "sample-A,S1,5000,4000,3000,100.0,0.9,0.8,0.5,0.95\n"
    "sample-A,S2,5000,4000,2000,50.0,0.7,0.4,0.1,0.80\n"
)


class TestIsSegmented(unittest.TestCase):
    def test_true_when_segment_column_present(self):
        df = pd.read_csv(io.StringIO(SEG_CSV))
        self.assertTrue(is_segmented(df))

    def test_false_without_segment_column(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        self.assertFalse(is_segmented(df))


class TestDedupeAndSumReads(unittest.TestCase):
    def test_segmented_dedupes_reads_and_sums_mapped(self):
        df = pd.read_csv(io.StringIO(SEG_CSV))
        out = dedupe_and_sum_reads(df)
        row = out.set_index("sample_name").loc["sample-A"]
        # reads / trimmed are identical across segment rows -> not doubled.
        self.assertEqual(int(row["number_of_reads"]), 5000)
        self.assertEqual(int(row["number_of_trim_paired_reads"]), 4000)
        # mapped is per-segment -> summed (3000 + 2000).
        self.assertEqual(int(row["number_of_mapped_reads"]), 5000)
        self.assertEqual(len(out), 1)

    def test_unsegmented_is_identity_per_sample(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        out = dedupe_and_sum_reads(df)
        self.assertEqual(len(out), 2)
        row = out.set_index("sample_name").loc["sample-A"]
        self.assertEqual(int(row["number_of_mapped_reads"]), 643219)


class TestDownsampleMinPool(unittest.TestCase):
    def test_caps_points_and_preserves_dip(self):
        n = 20000
        positions = np.arange(1, n + 1)
        depths = np.full(n, 1000, dtype=float)
        dip_value = 3.0
        depths[12345] = dip_value  # a narrow, deep dip in the middle
        out_pos, out_depth = downsample_min_pool(positions, depths, max_points=2000)
        self.assertLessEqual(len(out_pos), 2000)
        self.assertEqual(len(out_pos), len(out_depth))
        # the dip's minimum must survive min-pooling (not averaged away).
        self.assertIn(dip_value, list(out_depth))
        self.assertEqual(min(out_depth), dip_value)

    def test_short_series_passthrough(self):
        positions = np.arange(1, 51)
        depths = np.arange(50, 0, -1).astype(float)
        out_pos, out_depth = downsample_min_pool(positions, depths, max_points=2000)
        self.assertEqual(len(out_pos), 50)
        self.assertTrue(np.array_equal(out_pos, positions))
        self.assertTrue(np.array_equal(out_depth, depths))


class TestResolveBasewisePath(unittest.TestCase):
    def test_unsegmented_template(self):
        p = resolve_basewise_path("/out/", "sample-A")
        self.assertEqual(p, "/out/assembly/coverage_stats/sample-A.table_cov_basewise.txt")

    def test_segmented_template(self):
        p = resolve_basewise_path("/out/", "sample-A", segment="S1")
        self.assertEqual(p, "/out/assembly/S1/coverage_stats/sample-A.table_cov_basewise.txt")

    def test_no_trailing_slash(self):
        p = resolve_basewise_path("/out", "sample-A", segment="M")
        self.assertEqual(p, "/out/assembly/M/coverage_stats/sample-A.table_cov_basewise.txt")


class TestReadStatsAndTable(unittest.TestCase):
    def test_read_stats_summary(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "s.csv")
            with open(path, "w") as fh:
                fh.write(UNSEG_CSV)
            df = read_stats_summary(path)
        self.assertIn("sample_name", df.columns)
        self.assertEqual(len(df), 2)

    def test_table_formats_percentages_and_depth(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        html = build_stats_table_html(df)
        # horizontal_coverage is a fraction in [0,1] -> rendered x100 with %.
        self.assertIn("99.6%", html)  # 0.99595... -> 99.6%
        # depth is rounded for display.
        self.assertIn("3245.2", html)
        # raw fraction strings must not leak through.
        self.assertNotIn("0.9959535832525165", html)
        # sample names appear verbatim.
        self.assertIn("sample-A", html)
        self.assertIn("sample-B", html)

    def test_counts_grouped_with_thousands_separators_display_only(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        html = build_stats_table_html(df)
        # 6-7 digit read counts render grouped for the reader...
        self.assertIn("348,018", html)  # number_of_reads
        self.assertIn("330,568", html)  # number_of_trim_paired_reads
        # ...mean depth groups its integer part but keeps one fractional digit.
        self.assertIn("3,245.2", html)
        # ...while the raw ungrouped value is preserved as the numeric sort key.
        self.assertIn('data-sort="348018"', html)

    def test_table_drops_percentage_above_columns_and_humanizes_headers(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        html = build_stats_table_html(df)
        # the percentage_above_*x columns are removed entirely...
        self.assertNotIn("percentage_above", html)
        self.assertNotIn("94.8%", html)  # was percentage_above_1000x for sample-A
        # ...raw snake_case headers are replaced with human-readable labels.
        self.assertNotIn("sample_name", html)
        self.assertNotIn("average_depth", html)
        self.assertIn(">Sample<", html)
        self.assertIn("Total reads", html)
        self.assertIn("Mean depth", html)
        self.assertIn("Genome coverage", html)


class TestReportMetadata(unittest.TestCase):
    def test_illumina_config_is_paired_with_qc(self):
        meta = build_report_metadata({"data": "illumina", "scheme": "NA"}, "/no/dir")
        self.assertEqual(meta["platform"], "illumina")
        self.assertEqual(meta["library_layout"], "paired")
        self.assertTrue(meta["qc_performed"])
        self.assertIsNone(meta["primer_scheme"])

    def test_nanopore_config_is_single_without_qc(self):
        meta = build_report_metadata({"data": "nanopore", "scheme": "NA"}, "/no/dir")
        self.assertEqual(meta["platform"], "nanopore")
        self.assertEqual(meta["library_layout"], "single")
        self.assertFalse(meta["qc_performed"])

    def test_primer_scheme_kept_when_not_na(self):
        meta = build_report_metadata(
            {"data": "illumina", "scheme": "schemes/sarscov2.primers.bed"}, "/no/dir"
        )
        self.assertEqual(meta["primer_scheme"], "schemes/sarscov2.primers.bed")

    def test_infers_illumina_from_qc_dir_when_no_config(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "qc"))
            meta = build_report_metadata(None, d)
        self.assertEqual(meta["platform"], "illumina")
        self.assertEqual(meta["library_layout"], "paired")
        self.assertTrue(meta["qc_performed"])

    def test_infers_nanopore_when_no_config_and_no_qc_dir(self):
        with tempfile.TemporaryDirectory() as d:
            meta = build_report_metadata(None, d)
        self.assertEqual(meta["platform"], "nanopore")
        self.assertEqual(meta["library_layout"], "single")
        self.assertFalse(meta["qc_performed"])


class TestMappingRate(unittest.TestCase):
    def test_paired_mapped_column_is_a_rate_with_raw_count_tooltip(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        html = build_stats_table_html(df, paired=True)
        # sample-A: 643219 / (2 x 330568) = 97.3% (mapped is individual reads,
        # QC-passed is read pairs -> double the denominator when paired).
        self.assertIn("97.3%", html)
        self.assertIn(">Mapped %<", html)
        # the raw count is preserved in the cell tooltip, not as the visible text.
        self.assertIn("643,219 mapped reads", html)
        self.assertNotIn(">643,219<", html)

    def test_single_end_uses_undoubled_denominator(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        html = build_stats_table_html(df, paired=False)
        # single-end (Nanopore): mapped / QC-passed, no doubling. sample-B is
        # 800 / 900 = 88.9%; paired would instead give 800 / 1800 = 44.4%.
        self.assertIn("88.9%", html)
        self.assertNotIn("44.4%", html)

    def test_sortable_headers_have_accessibility_attributes(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        html = build_stats_table_html(df, paired=True)
        self.assertIn('role="button"', html)
        self.assertIn('aria-sort="none"', html)
        self.assertIn('tabindex="0"', html)


class TestLoadBasewiseTable(unittest.TestCase):
    def test_missing_file_returns_empty_frame(self):
        df = load_basewise_table("/no/such/file.txt")
        self.assertTrue(df.empty)
        self.assertEqual(list(df.columns), ["reference_id", "position", "depth"])

    def test_reads_tab_separated_no_header(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "cov.txt")
            with open(path, "w") as fh:
                fh.write("NC_007373.1|Influenza_A|H3N2|segment_1\t1\t0\n")
                fh.write("NC_007373.1|Influenza_A|H3N2|segment_1\t2\t15\n")
            df = load_basewise_table(path)
        self.assertEqual(list(df.columns), ["reference_id", "position", "depth"])
        self.assertEqual(len(df), 2)
        # reference_id with pipes is kept opaque (not split).
        self.assertEqual(df.iloc[0]["reference_id"], "NC_007373.1|Influenza_A|H3N2|segment_1")
        self.assertEqual(int(df.iloc[1]["depth"]), 15)


class TestFigureBuilders(unittest.TestCase):
    def test_reads_histogram_returns_figure(self):
        df = dedupe_and_sum_reads(pd.read_csv(io.StringIO(UNSEG_CSV)))
        fig = build_reads_histogram(df)
        self.assertIsInstance(fig, go.Figure)
        self.assertGreaterEqual(len(fig.data), 1)

    def test_throughput_has_distinct_hues_and_mapped_percent_panel(self):
        df = dedupe_and_sum_reads(pd.read_csv(io.StringIO(UNSEG_CSV)))
        fig = build_reads_histogram(df, paired=True)
        by_name = {t.name: t for t in fig.data}
        self.assertEqual(set(by_name), {"Total reads", "QC-passed reads", "Mapped %"})
        # three distinct hues, none dimmed by opacity.
        colors = {t.marker.color for t in fig.data}
        self.assertEqual(len(colors), 3)
        for t in fig.data:
            self.assertIn(t.opacity, (None, 1.0))
        # lower panel is the mapping rate on a fixed 0-100 axis.
        mapped = by_name["Mapped %"]
        self.assertAlmostEqual(mapped.y[0], 643219 / (2 * 330568) * 100, places=1)
        self.assertEqual(tuple(fig.layout.yaxis2.range), (0, 100))

    def test_aggregated_coverage_is_linear_with_honest_zeros(self):
        series = {"sample-A": (np.array([1, 2, 3, 4]), np.array([0.0, 10.0, 100.0, 0.0]))}
        fig = build_aggregated_coverage_line_plot(series)
        # linear by default (Plotly leaves type unset for linear axes).
        self.assertIn(fig.layout.yaxis.type, (None, "linear"))
        # zeros are plotted as zero, never clamped to 1 for a log axis.
        self.assertEqual(list(fig.data[0].y), [0.0, 10.0, 100.0, 0.0])
        # y-range is fit to the data and anchored at zero.
        self.assertEqual(fig.layout.yaxis.range[0], 0.0)
        self.assertAlmostEqual(fig.layout.yaxis.range[1], 108.0)


if __name__ == "__main__":
    unittest.main()
