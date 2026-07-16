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
import json
import os
import tempfile
import unittest

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from viralunity.scripts.python.generate_consensus_report import (
    EMPHASIS_MUTED,
    _json_for_script,
    _load_coverage_cache,
    _stats_rows_by_sample,
    build_aggregated_coverage_line_plot,
    build_annotation_model,
    build_kpi_summary,
    build_reads_histogram,
    build_report_metadata,
    build_stats_table_html,
    dedupe_and_sum_reads,
    downsample_min_pool,
    is_segmented,
    load_basewise_table,
    parse_gff3,
    parse_primer_bed,
    read_stats_summary,
    resolve_annotation_path,
    resolve_basewise_path,
    ReportParams,
    report_params_from_config,
)

_MODULE = "viralunity.scripts.python.generate_consensus_report"

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


class TestKpiSummary(unittest.TestCase):
    def test_global_unsegmented_counts_and_median(self):
        df = pd.read_csv(io.StringIO(UNSEG_CSV))
        kpi = build_kpi_summary(df, {}, segmented=False)
        g = kpi["global"]
        self.assertEqual(g["samples"], 2)
        self.assertEqual(g["ge90"], 1)  # sample-A 0.996; sample-B 0.60 fails
        # median horizontal coverage = median(0.9959…, 0.60)
        self.assertAlmostEqual(g["median_coverage"], (0.9959535832525165 + 0.6) / 2, places=4)
        self.assertAlmostEqual(g["median_depth"], (3245.2445239608064 + 45.5) / 2, places=2)
        self.assertEqual(kpi["per_segment"], {})

    def test_segmented_global_is_length_weighted(self):
        df = pd.read_csv(io.StringIO(SEG_CSV))  # sample-A: S1 hc0.95/d100, S2 hc0.80/d50
        lengths = {("sample-A", "S1"): 3000, ("sample-A", "S2"): 1000}
        g = build_kpi_summary(df, lengths, segmented=True)["global"]
        self.assertEqual(g["samples"], 1)
        # weighted breadth = (0.95*3000 + 0.80*1000)/4000 = 0.9125 -> clears 90%
        self.assertEqual(g["ge90"], 1)
        self.assertAlmostEqual(g["median_coverage"], 0.9125, places=4)
        # weighted depth = (100*3000 + 50*1000)/4000 = 87.5
        self.assertAlmostEqual(g["median_depth"], 87.5, places=3)

    def test_segmented_weighting_actually_matters(self):
        # equal weights give breadth 0.875 < 0.90, so the sample would NOT pass —
        # confirms the length weighting is doing the work, not a coincidence.
        df = pd.read_csv(io.StringIO(SEG_CSV))
        equal = {("sample-A", "S1"): 1000, ("sample-A", "S2"): 1000}
        self.assertEqual(build_kpi_summary(df, equal, segmented=True)["global"]["ge90"], 0)

    def test_segmented_falls_back_to_equal_weight_when_no_track_lengths(self):
        # Coverage tracks unavailable for every segment (empty lengths) while the
        # stats CSV is valid: the whole-genome KPIs must reflect the CSV via an
        # equal-weight mean, not collapse to 0%/0x (regression guard for the
        # length-weights-vs-CSV-values divergence).
        df = pd.read_csv(io.StringIO(SEG_CSV))
        g = build_kpi_summary(df, {}, segmented=True)["global"]
        self.assertAlmostEqual(g["median_coverage"], (0.95 + 0.80) / 2, places=4)
        self.assertAlmostEqual(g["median_depth"], (100.0 + 50.0) / 2, places=3)
        self.assertNotEqual(g["median_coverage"], 0.0)

    def test_segmented_per_segment_blocks(self):
        df = pd.read_csv(io.StringIO(SEG_CSV))
        lengths = {("sample-A", "S1"): 3000, ("sample-A", "S2"): 1000}
        ps = build_kpi_summary(df, lengths, segmented=True)["per_segment"]
        self.assertEqual(set(ps), {"S1", "S2"})
        self.assertEqual(
            ps["S1"], {"samples": 1, "ge90": 1, "median_coverage": 0.95, "median_depth": 100.0}
        )
        self.assertEqual(
            ps["S2"], {"samples": 1, "ge90": 0, "median_coverage": 0.80, "median_depth": 50.0}
        )


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

    def test_warns_loudly_when_layout_is_inferred(self):
        # Inferring layout silently drives the mapping-rate x2 denominator, so the
        # config-less path must warn (a misidentified platform doubles/halves the
        # reported rate).
        with tempfile.TemporaryDirectory() as d:
            with self.assertLogs(_MODULE, level="WARNING") as cm:
                build_report_metadata(None, d)
        self.assertTrue(any("inferred platform" in m for m in cm.output))

    def test_no_warning_when_config_declares_layout(self):
        with self.assertNoLogs(_MODULE, level="WARNING"):
            build_report_metadata({"data": "illumina", "scheme": "NA"}, "/no/dir")


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


class TestParseGff3(unittest.TestCase):
    def _write(self, text):
        d = tempfile.mkdtemp()
        path = os.path.join(d, "genes.gff3")
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_missing_file_returns_empty_list(self):
        self.assertEqual(parse_gff3("/no/such.gff3", "chrA"), [])

    def test_extracts_genes_on_contig_with_1based_coords(self):
        gff = self._write(
            "##gff-version 3\n"
            "chrA\tRefSeq\tgene\t266\t21555\t.\t+\t.\tID=gene-orf1;Name=ORF1ab\n"
            "chrA\tRefSeq\tgene\t21563\t25384\t.\t+\t.\tID=gene-S;Name=S\n"
        )
        feats = parse_gff3(gff, "chrA")
        self.assertEqual(len(feats), 2)
        self.assertEqual(feats[0]["start"], 266)
        self.assertEqual(feats[0]["end"], 21555)
        self.assertEqual(feats[0]["name"], "ORF1ab")
        self.assertEqual(feats[1]["name"], "S")

    def test_excludes_features_on_other_contigs(self):
        gff = self._write(
            "chrA\tRefSeq\tgene\t1\t99\t.\t+\t.\tName=keep\n"
            "chrB\tRefSeq\tgene\t1\t99\t.\t+\t.\tName=drop\n"
        )
        feats = parse_gff3(gff, "chrA")
        self.assertEqual([f["name"] for f in feats], ["keep"])

    def test_label_precedence_prefers_name_then_gene_then_id(self):
        gff = self._write(
            "chrA\tx\tgene\t1\t9\t.\t+\t.\tID=g1;gene=orf1;Name=ORF1ab\n"
            "chrA\tx\tgene\t10\t19\t.\t+\t.\tID=g2;gene=Spike\n"
            "chrA\tx\tgene\t20\t29\t.\t+\t.\tID=g3\n"
        )
        feats = parse_gff3(gff, "chrA")
        self.assertEqual([f["name"] for f in feats], ["ORF1ab", "Spike", "g3"])

    def test_falls_back_to_cds_when_no_gene_features(self):
        gff = self._write(
            "chrA\tx\tCDS\t5\t50\t.\t+\t0\tID=cds1;Name=nsp1\n"
            "chrA\tx\tregion\t1\t100\t.\t+\t.\tID=r1;Name=whole\n"
        )
        feats = parse_gff3(gff, "chrA")
        self.assertEqual([f["name"] for f in feats], ["nsp1"])


class TestParsePrimerBed(unittest.TestCase):
    def _write(self, text):
        d = tempfile.mkdtemp()
        path = os.path.join(d, "scheme.bed")
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_missing_file_returns_empty_list(self):
        self.assertEqual(parse_primer_bed("/no/such.bed", "chrA"), [])

    def test_artic_pairs_left_right_into_pool_lanes_1based(self):
        bed = self._write(
            "chrA\t30\t54\tscheme_1_LEFT\tpool1\t+\n"
            "chrA\t385\t410\tscheme_1_RIGHT\tpool1\t-\n"
            "chrA\t320\t344\tscheme_2_LEFT\tpool2\t+\n"
            "chrA\t705\t730\tscheme_2_RIGHT\tpool2\t-\n"
        )
        lanes = parse_primer_bed(bed, "chrA")
        self.assertEqual([ln["label"] for ln in lanes], ["Pool A", "Pool B"])
        a = lanes[0]["features"]
        self.assertEqual(len(a), 1)
        # BED 0-based half-open [30,54) -> 1-based inclusive start 31; amplicon
        # spans LEFT.start..RIGHT.end = 31..410.
        self.assertEqual(a[0]["start"], 31)
        self.assertEqual(a[0]["end"], 410)
        self.assertEqual(lanes[1]["features"][0]["start"], 321)
        self.assertEqual(lanes[1]["features"][0]["end"], 730)

    def test_pool_by_amplicon_parity_when_no_pool_column(self):
        bed = self._write(
            "chrA\t0\t20\tscheme_1_LEFT\n"
            "chrA\t100\t120\tscheme_1_RIGHT\n"
            "chrA\t90\t110\tscheme_2_LEFT\n"
            "chrA\t200\t220\tscheme_2_RIGHT\n"
        )
        lanes = parse_primer_bed(bed, "chrA")
        self.assertEqual([ln["label"] for ln in lanes], ["Pool A", "Pool B"])
        self.assertEqual(lanes[0]["features"][0]["name"], "scheme_1")
        self.assertEqual(lanes[1]["features"][0]["name"], "scheme_2")

    def test_falls_back_to_individual_primers_when_unpairable(self):
        bed = self._write("chrA\t0\t20\tprimerX\n" "chrA\t50\t70\tprimerY\n")
        lanes = parse_primer_bed(bed, "chrA")
        self.assertEqual(len(lanes), 1)
        self.assertEqual(lanes[0]["label"], "Primers")
        self.assertEqual({f["name"] for f in lanes[0]["features"]}, {"primerX", "primerY"})
        self.assertEqual(lanes[0]["features"][0]["start"], 1)

    def test_excludes_other_contigs(self):
        bed = self._write(
            "chrA\t0\t20\ts_1_LEFT\n"
            "chrA\t100\t120\ts_1_RIGHT\n"
            "chrB\t0\t20\ts_1_LEFT\n"
            "chrB\t100\t120\ts_1_RIGHT\n"
        )
        lanes = parse_primer_bed(bed, "chrA")
        total = sum(len(ln["features"]) for ln in lanes)
        self.assertEqual(total, 1)

    def test_more_than_two_pools_yield_extra_lanes(self):
        bed = self._write(
            "chrA\t0\t20\ts_1_LEFT\tp1\n"
            "chrA\t100\t120\ts_1_RIGHT\tp1\n"
            "chrA\t90\t110\ts_2_LEFT\tp2\n"
            "chrA\t200\t220\ts_2_RIGHT\tp2\n"
            "chrA\t190\t210\ts_3_LEFT\tp3\n"
            "chrA\t300\t320\ts_3_RIGHT\tp3\n"
        )
        lanes = parse_primer_bed(bed, "chrA")
        self.assertEqual([ln["label"] for ln in lanes], ["Pool A", "Pool B", "Pool C"])


class TestResolveAnnotationPath(unittest.TestCase):
    def test_primer_path(self):
        self.assertEqual(
            resolve_annotation_path("/out/", "primer"),
            "/out/annotation/primer_scheme.bed",
        )

    def test_gene_path_unsegmented(self):
        self.assertEqual(
            resolve_annotation_path("/out/", "gene"),
            "/out/annotation/gene_annotation.gff3",
        )

    def test_gene_path_segmented(self):
        self.assertEqual(
            resolve_annotation_path("/out/", "gene", segment="S1"),
            "/out/annotation/S1.gene_annotation.gff3",
        )


class TestCoverageCacheContigCapture(unittest.TestCase):
    def test_captures_contig_name_per_segment(self):
        with tempfile.TemporaryDirectory() as d:
            cov_dir = os.path.join(d, "assembly", "coverage_stats")
            os.makedirs(cov_dir)
            with open(os.path.join(cov_dir, "sample-A.table_cov_basewise.txt"), "w") as fh:
                fh.write("chrTEST\t1\t10\nchrTEST\t2\t20\n")
            df = pd.DataFrame({"sample_name": ["sample-A"]})
            _cache, _lengths, contigs = _load_coverage_cache(d, df, segmented=False)
        self.assertEqual(contigs[None], "chrTEST")


class TestBuildAnnotationModel(unittest.TestCase):
    def _dir_with_annotation(self):
        d = tempfile.mkdtemp()
        ann = os.path.join(d, "annotation")
        os.makedirs(ann)
        with open(os.path.join(ann, "gene_annotation.gff3"), "w") as fh:
            fh.write("chrTEST\tx\tgene\t1\t100\t.\t+\t.\tName=G1\n")
        with open(os.path.join(ann, "primer_scheme.bed"), "w") as fh:
            fh.write("chrTEST\t0\t20\ts_1_LEFT\tp1\n")
            fh.write("chrTEST\t80\t100\ts_1_RIGHT\tp1\n")
        return d

    def test_builds_gene_and_primer_lanes(self):
        d = self._dir_with_annotation()
        model = build_annotation_model(d, [None], {None: "chrTEST"})
        self.assertTrue(model["has_genes"])
        self.assertTrue(model["has_primers"])
        lanes = model["by_segment"][""]["lanes"]
        kinds = [ln["kind"] for ln in lanes]
        self.assertIn("gene", kinds)
        self.assertIn("primer", kinds)

    def test_no_annotation_dir_is_empty_with_flags_false(self):
        with tempfile.TemporaryDirectory() as d:
            model = build_annotation_model(d, [None], {None: "chrTEST"})
        self.assertFalse(model["has_genes"])
        self.assertFalse(model["has_primers"])
        self.assertEqual(model["by_segment"], {})

    def test_contig_mismatch_draws_no_features(self):
        d = self._dir_with_annotation()
        model = build_annotation_model(d, [None], {None: "chrOTHER"})
        self.assertFalse(model["has_genes"])
        self.assertFalse(model["has_primers"])


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

    def test_aggregated_over_8_samples_falls_back_to_emphasis(self):
        # >8 series exceeds the CVD-safe categorical ceiling, so identity is no
        # longer colour-coded: a single muted hue, dimmed, with the legend off
        # (the stats table carries per-sample identity instead).
        series = {f"s{i}": (np.array([1, 2, 3]), np.array([10.0, 20.0, 30.0])) for i in range(9)}
        fig = build_aggregated_coverage_line_plot(series)
        self.assertFalse(fig.layout.showlegend)
        self.assertTrue(all(t.line.color == EMPHASIS_MUTED for t in fig.data))
        self.assertTrue(all(t.opacity == 0.6 for t in fig.data))

    def test_aggregated_8_samples_keeps_categorical_legend(self):
        # At the 8-series ceiling the categorical palette is still used (legend on,
        # no forced muted hue).
        series = {f"s{i}": (np.array([1, 2, 3]), np.array([10.0, 20.0, 30.0])) for i in range(8)}
        fig = build_aggregated_coverage_line_plot(series)
        self.assertTrue(fig.layout.showlegend)
        self.assertFalse(any(t.line.color == EMPHASIS_MUTED for t in fig.data))

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


def _one_sample_df(name):
    return pd.DataFrame(
        [
            {
                "sample_name": name,
                "number_of_reads": 1000,
                "number_of_trim_paired_reads": 900,
                "number_of_mapped_reads": 800,
                "average_depth": 45.5,
                "horizontal_coverage": 0.6,
            }
        ]
    )


class TestReportParams(unittest.TestCase):
    def test_defaults(self):
        p = ReportParams()
        self.assertEqual(p.pass_threshold, 0.90)
        self.assertEqual(p.warn_threshold, 0.70)
        self.assertEqual(p.chart_color, "#2a78d6")
        self.assertEqual(p.colorbar_thickness, 14)

    def test_pct_labels_are_derived(self):
        p = ReportParams(pass_threshold=0.85, warn_threshold=0.5)
        self.assertEqual(p.pass_pct_label, "85%")
        self.assertEqual(p.warn_pct_label, "50%")

    def test_tier_boundaries(self):
        p = ReportParams(pass_threshold=0.90, warn_threshold=0.70)
        self.assertEqual(p.tier(0.90), "pass")  # inclusive lower bound
        self.assertEqual(p.tier(0.95), "pass")
        self.assertEqual(p.tier(0.70), "warn")  # inclusive lower bound
        self.assertEqual(p.tier(0.89), "warn")
        self.assertEqual(p.tier(0.69), "fail")
        self.assertEqual(p.tier("not-a-number"), "fail")

    def test_client_dict_has_labels_and_hues(self):
        d = ReportParams().as_client_dict()
        self.assertEqual(d["passPctLabel"], "90%")
        self.assertEqual(d["warnPctLabel"], "70%")
        self.assertIn("fail", d["tierColors"])

    def test_from_config_precedence(self):
        # explicit override beats config beats default; None overrides ignored.
        p = report_params_from_config(
            {"report_pass_threshold": 0.8, "report_chart_color": "#123456"},
            warn_threshold=0.5,
            pass_threshold=None,
        )
        self.assertEqual(p.pass_threshold, 0.8)  # from config (override was None)
        self.assertEqual(p.warn_threshold, 0.5)  # from override
        self.assertEqual(p.chart_color, "#123456")  # from config
        self.assertEqual(p.colorbar_thickness, 14)  # default

    def test_from_config_empty(self):
        self.assertEqual(report_params_from_config(None), ReportParams())


class TestHtmlEscaping(unittest.TestCase):
    """Data-controlled strings (sample/segment names) must not inject markup —
    the report is a self-contained file designed to be shared."""

    def test_stats_table_escapes_sample_names(self):
        df = _one_sample_df('"><img src=x onerror=alert(1)>')
        html = build_stats_table_html(df, paired=False)
        self.assertNotIn("<img src=x onerror=alert(1)>", html)
        self.assertIn("&lt;img", html)
        # the attribute-breaking quote is neutralized too.
        self.assertNotIn('data-sort=""><img', html)

    def test_stats_table_escapes_malicious_numeric_columns(self):
        # A numeric column carrying non-numeric markup (a corrupt/crafted stats
        # CSV) must not inject into the cell text or the data-sort attribute; the
        # formatters fall back to str(value), so escaping happens at cell build.
        payload = '"><img src=x onerror=alert(1)>'
        df = _one_sample_df("s1")
        df["number_of_reads"] = payload  # INT_COLS branch
        df["average_depth"] = payload  # average_depth branch
        html = build_stats_table_html(df, paired=False)
        self.assertNotIn("<img src=x onerror=alert(1)>", html)
        self.assertIn("&lt;img", html)
        self.assertNotIn('data-sort=""><img', html)

    def test_mapping_rate_tooltip_escapes_malicious_count(self):
        # The mapped-reads tooltip is built from _fmt_int(raw); a malicious raw
        # count must not break out of the title="" attribute.
        df = _one_sample_df("s1")
        df["number_of_mapped_reads"] = '"><img src=x onerror=alert(1)>'
        html = build_stats_table_html(df, paired=True)
        self.assertNotIn('title=""><img', html)
        self.assertNotIn("<img src=x onerror=alert(1)>", html)

    def test_stats_rows_by_sample_escapes_cell_values(self):
        df = _one_sample_df("<b>x</b>")
        rows = _stats_rows_by_sample(df, segmented=False, paired=False)
        # the key is the raw name (a JS object key, not HTML); the rendered
        # "Sample" cell value (inserted via innerHTML client-side) is escaped.
        self.assertEqual(rows["<b>x</b>"][0]["Sample"], "&lt;b&gt;x&lt;/b&gt;")

    def test_json_for_script_neutralizes_script_breakout(self):
        payload = {"sample": "</script><script>alert(1)</script>"}
        s = _json_for_script(payload)
        self.assertNotIn("</script>", s)
        self.assertIn("\\u003c", s)
        # still valid JSON that round-trips back to the original characters.
        self.assertEqual(json.loads(s), payload)


if __name__ == "__main__":
    unittest.main()
