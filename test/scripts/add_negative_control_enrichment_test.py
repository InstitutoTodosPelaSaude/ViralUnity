"""Tests for viralunity.scripts.python.add_negative_control_enrichment.

Covers:
  - calculate_fold_enrichment / calculate_log2_ratio / calculate_z_score helpers
  - apply_negative_control_enrichment: all three control-count tiers (0/1/≥2)
  - Decision-metric selection (rpkm preferred over rpm when available)
  - neg_pass correctness and NA-as-keep contract
  - Control-SD=0 z-score fallback to log2-ratio gate
  - Taxa absent from controls (zero-background assumption)
  - Output column schema
"""

import math
import unittest

import pandas as pd

from viralunity.scripts.python.add_negative_control_enrichment import (
    apply_negative_control_enrichment,
    calculate_fold_enrichment,
    calculate_log2_ratio,
    calculate_z_score,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helper factories
# ─────────────────────────────────────────────────────────────────────────────

def _row(
    sample: str,
    taxid: str,
    rpm: float,
    rpkm=None,
    rank: str = "species",
    tool: str = "kraken2",
    mode: str = "reads",
    total_reads: int = 1_000_000,
    count: int = 0,
) -> dict:
    r = {
        "sample": sample,
        "tool": tool,
        "mode": mode,
        "rank": rank,
        "taxid": taxid,
        "count": count,
        "total_reads": total_reads,
        "rpm": rpm,
    }
    if rpkm is not None:
        r["rpkm"] = rpkm
    return r


def _df(*rows) -> pd.DataFrame:
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for the three helper functions
# ─────────────────────────────────────────────────────────────────────────────

class TestCalculateFoldEnrichment(unittest.TestCase):
    def test_double_sample_gives_close_to_2(self):
        fe = calculate_fold_enrichment(sample_metric=10.0, control_mean=5.0, pseudocount=0.0)
        self.assertAlmostEqual(fe, 2.0)

    def test_with_pseudocount(self):
        fe = calculate_fold_enrichment(sample_metric=10.0, control_mean=0.0, pseudocount=1.0)
        self.assertAlmostEqual(fe, 11.0)

    def test_equal_sample_and_control_gives_1(self):
        fe = calculate_fold_enrichment(sample_metric=5.0, control_mean=5.0, pseudocount=0.0)
        self.assertAlmostEqual(fe, 1.0)

    def test_zero_both_plus_pseudocount(self):
        fe = calculate_fold_enrichment(sample_metric=0.0, control_mean=0.0, pseudocount=1.0)
        self.assertAlmostEqual(fe, 1.0)


class TestCalculateLog2Ratio(unittest.TestCase):
    def test_double_gives_log2_1_equals_1(self):
        l2r = calculate_log2_ratio(sample_metric=10.0, control_mean=5.0, pseudocount=0.0)
        self.assertAlmostEqual(l2r, 1.0)

    def test_equal_gives_zero(self):
        l2r = calculate_log2_ratio(sample_metric=5.0, control_mean=5.0, pseudocount=0.0)
        self.assertAlmostEqual(l2r, 0.0)

    def test_sample_lower_than_control_negative(self):
        l2r = calculate_log2_ratio(sample_metric=2.0, control_mean=8.0, pseudocount=0.0)
        self.assertLess(l2r, 0.0)

    def test_pseudocount_stabilises_zero_control(self):
        l2r = calculate_log2_ratio(sample_metric=100.0, control_mean=0.0, pseudocount=1.0)
        self.assertGreater(l2r, 0.0)
        self.assertFalse(math.isinf(l2r))

    def test_log2_relationship_to_fold_enrichment(self):
        pc = 1.0
        sample, ctrl = 10.0, 3.0
        fe  = calculate_fold_enrichment(sample, ctrl, pc)
        l2r = calculate_log2_ratio(sample, ctrl, pc)
        self.assertAlmostEqual(l2r, math.log2(fe), places=10)


class TestCalculateZScore(unittest.TestCase):
    def test_basic_z_score(self):
        # controls: [2, 4] → mean=3, sd=√2 ≈ 1.414
        z = calculate_z_score(6.0, [2.0, 4.0])
        self.assertAlmostEqual(z, (6.0 - 3.0) / ((4.0 - 2.0) / (2 ** 0.5)), places=5)

    def test_returns_none_with_one_control(self):
        z = calculate_z_score(6.0, [3.0], min_controls=2)
        self.assertIsNone(z)

    def test_returns_none_with_zero_controls(self):
        z = calculate_z_score(6.0, [], min_controls=2)
        self.assertIsNone(z)

    def test_returns_none_when_sd_zero(self):
        z = calculate_z_score(6.0, [5.0, 5.0])
        self.assertIsNone(z)

    def test_returns_zero_when_sample_equals_mean(self):
        z = calculate_z_score(5.0, [4.0, 6.0])
        self.assertAlmostEqual(z, 0.0)

    def test_negative_z_when_sample_below_mean(self):
        z = calculate_z_score(1.0, [4.0, 6.0])
        self.assertIsNotNone(z)
        self.assertLess(z, 0.0)


# ─────────────────────────────────────────────────────────────────────────────
# Integration tests for apply_negative_control_enrichment
# ─────────────────────────────────────────────────────────────────────────────

class TestZeroControls(unittest.TestCase):
    """n_controls == 0 → all enrichment cols are NA, neg_pass is NA."""

    def setUp(self):
        self.df = _df(
            _row("S1", "3001", rpm=100.0),
            _row("S2", "3001", rpm=50.0),
        )

    def test_neg_pass_is_na_for_all_rows(self):
        out = apply_negative_control_enrichment(self.df, negatives=[])
        self.assertTrue(out["neg_pass"].isna().all())

    def test_neg_decision_is_none_string(self):
        out = apply_negative_control_enrichment(self.df, negatives=[])
        self.assertTrue((out["neg_decision"] == "none").all())

    def test_enrichment_cols_are_na(self):
        out = apply_negative_control_enrichment(self.df, negatives=[])
        for col in ("fold_enrichment", "log2_ratio", "z_score"):
            self.assertTrue(out[col].isna().all(), f"Expected NA for {col}")

    def test_n_negative_controls_is_zero(self):
        out = apply_negative_control_enrichment(self.df, negatives=[])
        self.assertTrue((out["n_negative_controls"] == 0).all())


class TestSingleControl(unittest.TestCase):
    """n_controls == 1 → log2_ratio gate, no z-score."""

    def setUp(self):
        self.df = _df(
            _row("CTRL", "3001", rpm=10.0),   # negative control
            _row("S1",   "3001", rpm=80.0),   # high enrichment
            _row("S2",   "3001", rpm=12.0),   # close to control
        )

    def test_neg_decision_is_log2_ratio(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"], log2_ratio_threshold=1.0
        )
        for _, row in out[out["sample"] != "CTRL"].iterrows():
            self.assertEqual(row["neg_decision"], "log2_ratio")

    def test_high_enrichment_passes(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"], log2_ratio_threshold=1.0, pseudocount=1.0
        )
        s1 = out[out["sample"] == "S1"].iloc[0]
        # log2((80+1)/(10+1)) ≈ 2.88 → passes threshold of 1.0
        self.assertTrue(s1["neg_pass"])

    def test_near_background_fails(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"], log2_ratio_threshold=1.0, pseudocount=1.0
        )
        s2 = out[out["sample"] == "S2"].iloc[0]
        # log2((12+1)/(10+1)) ≈ 0.24 → fails threshold of 1.0
        self.assertFalse(s2["neg_pass"])

    def test_z_score_is_na_for_single_control(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"]
        )
        self.assertTrue(
            out[out["sample"] != "CTRL"]["z_score"].isna().all()
        )

    def test_control_row_neg_pass_is_na(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"], log2_ratio_threshold=1.0
        )
        ctrl = out[out["sample"] == "CTRL"].iloc[0]
        self.assertTrue(pd.isna(ctrl["neg_pass"]))


class TestMultipleControls(unittest.TestCase):
    """n_controls >= 2 → z-score gate with log2-ratio fallback."""

    def setUp(self):
        # Three negative controls around RPM ~5, two biological samples
        self.df = _df(
            _row("CTRL1", "3001", rpm=4.0),
            _row("CTRL2", "3001", rpm=5.0),
            _row("CTRL3", "3001", rpm=6.0),
            _row("S_high", "3001", rpm=200.0),   # clearly above background
            _row("S_low",  "3001", rpm=5.5),     # within background
        )

    def test_neg_decision_is_z_score(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2", "CTRL3"], z_score_threshold=3.0
        )
        sample_rows = out[~out["sample"].isin(["CTRL1", "CTRL2", "CTRL3"])]
        self.assertTrue((sample_rows["neg_decision"] == "z_score").all())

    def test_high_rpm_sample_passes(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2", "CTRL3"], z_score_threshold=3.0
        )
        s_high = out[out["sample"] == "S_high"].iloc[0]
        self.assertTrue(s_high["neg_pass"])

    def test_low_rpm_sample_fails(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2", "CTRL3"], z_score_threshold=3.0
        )
        s_low = out[out["sample"] == "S_low"].iloc[0]
        self.assertFalse(s_low["neg_pass"])

    def test_n_negative_controls_correct(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2", "CTRL3"]
        )
        self.assertTrue((out["n_negative_controls"] == 3).all())

    def test_control_stats_populated(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2", "CTRL3"]
        )
        row = out.iloc[0]  # any row (stats are per taxon group)
        # Mean of [4, 5, 6] = 5
        self.assertAlmostEqual(row["control_mean"], 5.0, places=5)


class TestZeroControlSD(unittest.TestCase):
    """When all controls have identical metric → SD=0, z undefined → fallback to log2-ratio."""

    def setUp(self):
        # Both controls have RPM 10.0 → SD = 0
        self.df = _df(
            _row("CTRL1", "3001", rpm=10.0),
            _row("CTRL2", "3001", rpm=10.0),
            _row("S1",    "3001", rpm=50.0),
        )

    def test_z_score_is_none_when_sd_zero(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2"]
        )
        s1 = out[out["sample"] == "S1"].iloc[0]
        self.assertTrue(pd.isna(s1["z_score"]))

    def test_falls_back_to_log2_ratio_gate(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2"], log2_ratio_threshold=1.0
        )
        s1 = out[out["sample"] == "S1"].iloc[0]
        self.assertEqual(s1["neg_decision"], "log2_ratio_fallback")
        # log2((50+1)/(10+1)) ≈ 2.21 → passes threshold 1.0
        self.assertTrue(s1["neg_pass"])

    def test_falls_back_and_fails_when_below_threshold(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2"], log2_ratio_threshold=10.0
        )
        s1 = out[out["sample"] == "S1"].iloc[0]
        self.assertEqual(s1["neg_decision"], "log2_ratio_fallback")
        self.assertFalse(s1["neg_pass"])


class TestAbsentFromControls(unittest.TestCase):
    """Taxa absent from negative controls → control_mean=0, computed against pseudocount."""

    def setUp(self):
        # CTRL1/CTRL2 have taxid 3001 at rpm 5, but S1 has a taxid 9999 not in controls.
        self.df = _df(
            _row("CTRL1", "3001", rpm=5.0),
            _row("CTRL2", "3001", rpm=5.0),
            _row("CTRL1", "9999", rpm=0.0),  # absent in controls (zero reads)
            _row("CTRL2", "9999", rpm=0.0),
            _row("S1",    "3001", rpm=100.0),
            _row("S1",    "9999", rpm=50.0),
        )

    def test_absent_taxon_uses_zero_background(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2"]
        )
        s1_9999 = out[(out["sample"] == "S1") & (out["taxid"] == "9999")].iloc[0]
        # control_mean for 9999 is 0 (zero rpm in controls)
        self.assertAlmostEqual(s1_9999["control_mean"], 0.0, places=5)

    def test_absent_taxon_still_gets_enrichment_metrics(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2"]
        )
        s1_9999 = out[(out["sample"] == "S1") & (out["taxid"] == "9999")].iloc[0]
        self.assertFalse(pd.isna(s1_9999["fold_enrichment"]))
        self.assertFalse(pd.isna(s1_9999["log2_ratio"]))

    def test_absent_taxon_with_high_sample_passes(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL1", "CTRL2"],
            z_score_threshold=3.0, log2_ratio_threshold=1.0
        )
        s1_9999 = out[(out["sample"] == "S1") & (out["taxid"] == "9999")].iloc[0]
        # z is undefined (sd for 9999 in controls is 0), falls back to log2-ratio
        # log2((50+1)/(0+1)) ≈ 5.67 > 1.0 → passes
        self.assertTrue(s1_9999["neg_pass"])


class TestDecisionMetricSelection(unittest.TestCase):
    """RPKM should be preferred over RPM when available and non-NA."""

    def setUp(self):
        # Control has rpkm 2.0, sample has rpkm 100.0 (very high)
        # rpm values would give different log2-ratios
        self.df = _df(
            _row("CTRL", "3001", rpm=20.0, rpkm=2.0),
            _row("S1",   "3001", rpm=200.0, rpkm=100.0),
        )

    def test_neg_metric_is_rpkm_when_rpkm_column_present(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"]
        )
        self.assertTrue((out["neg_metric"] == "rpkm").all())

    def test_log2_ratio_uses_rpkm_values(self):
        out = apply_negative_control_enrichment(
            self.df, negatives=["CTRL"], pseudocount=1.0
        )
        s1 = out[out["sample"] == "S1"].iloc[0]
        expected_l2r = math.log2((100.0 + 1.0) / (2.0 + 1.0))
        self.assertAlmostEqual(s1["log2_ratio"], expected_l2r, places=5)

    def test_neg_metric_falls_back_to_rpm_when_rpkm_all_na(self):
        df = _df(
            _row("CTRL", "3001", rpm=20.0, rpkm=float("nan")),
            _row("S1",   "3001", rpm=200.0, rpkm=float("nan")),
        )
        out = apply_negative_control_enrichment(df, negatives=["CTRL"])
        self.assertTrue((out["neg_metric"] == "rpm").all())

    def test_neg_metric_is_rpm_when_no_rpkm_column(self):
        df = _df(
            _row("CTRL", "3001", rpm=20.0),  # no rpkm key
            _row("S1",   "3001", rpm=200.0),
        )
        out = apply_negative_control_enrichment(df, negatives=["CTRL"])
        self.assertTrue((out["neg_metric"] == "rpm").all())


class TestNAAsKeepContract(unittest.TestCase):
    """When neg_pass is NA the Krona filter keeps the taxon (conservative).

    This test verifies that zero-control rows have NA neg_pass, not False.
    """

    def test_zero_control_neg_pass_is_na_not_false(self):
        df = _df(
            _row("S1", "3001", rpm=100.0),
            _row("S2", "3001", rpm=200.0),
        )
        out = apply_negative_control_enrichment(df, negatives=[])
        self.assertFalse(out["neg_pass"].eq(False).any(), "Expected NA, not False")
        self.assertTrue(out["neg_pass"].isna().all())


class TestOutputSchema(unittest.TestCase):
    """All expected output columns must be present."""

    EXPECTED_COLS = {
        "is_negative_control",
        "n_negative_controls",
        "neg_metric",
        "control_mean",
        "control_sd",
        "control_median",
        "control_max",
        "fold_enrichment",
        "log2_ratio",
        "z_score",
        "enrichment_pseudocount",
        "z_score_threshold_used",
        "log2_ratio_threshold_used",
        "neg_decision",
        "neg_pass",
    }

    def _out(self, negatives=None):
        df = _df(
            _row("CTRL", "3001", rpm=10.0),
            _row("S1",   "3001", rpm=50.0),
        )
        return apply_negative_control_enrichment(
            df, negatives=negatives or []
        )

    def test_all_columns_present_zero_controls(self):
        out = self._out(negatives=[])
        missing = self.EXPECTED_COLS - set(out.columns)
        self.assertEqual(missing, set(), f"Missing columns: {missing}")

    def test_all_columns_present_with_control(self):
        out = self._out(negatives=["CTRL"])
        missing = self.EXPECTED_COLS - set(out.columns)
        self.assertEqual(missing, set(), f"Missing columns: {missing}")

    def test_original_columns_preserved(self):
        df = _df(_row("S1", "3001", rpm=50.0))
        out = apply_negative_control_enrichment(df, negatives=[])
        for col in df.columns:
            self.assertIn(col, out.columns, f"Original column '{col}' missing from output")

    def test_row_count_unchanged(self):
        df = _df(
            _row("CTRL", "3001", rpm=10.0),
            _row("S1",   "3001", rpm=50.0),
            _row("S2",   "3001", rpm=60.0),
        )
        out = apply_negative_control_enrichment(df, negatives=["CTRL"])
        self.assertEqual(len(out), len(df))


class TestErrorCases(unittest.TestCase):
    def test_missing_rpm_column_raises(self):
        df = pd.DataFrame([{"sample": "S1", "taxid": "3001"}])
        with self.assertRaises(ValueError):
            apply_negative_control_enrichment(df, negatives=[])

    def test_missing_sample_column_raises(self):
        df = pd.DataFrame([{"taxid": "3001", "rpm": 10.0}])
        with self.assertRaises(ValueError):
            apply_negative_control_enrichment(df, negatives=[])

    def test_negatives_not_in_table_raises(self):
        df = _df(_row("S1", "3001", rpm=50.0))
        with self.assertRaises(ValueError):
            apply_negative_control_enrichment(df, negatives=["CTRL_MISSING"])


class TestPseudocountEffect(unittest.TestCase):
    def test_higher_pseudocount_reduces_fold_enrichment(self):
        """With a larger pseudocount both metrics get pulled towards 1."""
        df = _df(
            _row("CTRL", "3001", rpm=0.0),
            _row("S1",   "3001", rpm=100.0),
        )
        out_pc1 = apply_negative_control_enrichment(df, negatives=["CTRL"], pseudocount=1.0)
        out_pc10 = apply_negative_control_enrichment(df, negatives=["CTRL"], pseudocount=10.0)
        fe1  = out_pc1[out_pc1["sample"] == "S1"].iloc[0]["fold_enrichment"]
        fe10 = out_pc10[out_pc10["sample"] == "S1"].iloc[0]["fold_enrichment"]
        # Larger pseudocount → smaller fold-enrichment (shrunk towards 1)
        self.assertGreater(fe1, fe10)

    def test_pseudocount_recorded_in_output(self):
        df = _df(_row("S1", "3001", rpm=50.0))
        out = apply_negative_control_enrichment(df, negatives=[], pseudocount=2.5)
        self.assertAlmostEqual(out.iloc[0]["enrichment_pseudocount"], 2.5)


if __name__ == "__main__":
    unittest.main()
