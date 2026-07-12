"""Tests for viralunity.scripts.python.apply_max_rpm_bleed_filter.

Covers:
  - infer_group_cols (the shared taxon grouping)
  - RPM-only behaviour (no --viral-genomes)
  - RPKM auto-selection and the within-taxon ratio INVARIANCE (rpkm vs rpm give
    the same bleed_pass when the filter is applied under both)
  - the floor gate being the only place rpkm vs rpm actually diverges
  - metric forcing and error cases
"""

import unittest

import pandas as pd

from viralunity.scripts.python.apply_max_rpm_bleed_filter import (
    apply_bleed_filter,
    infer_group_cols,
)


def _row(sample, taxid, rpm, rpkm=None, rank="species", tool="diamond", mode="contigs"):
    r = {
        "sample": sample,
        "tool": tool,
        "mode": mode,
        "rank": rank,
        "taxid": taxid,
        "rpm": rpm,
    }
    if rpkm is not None:
        r["rpkm"] = rpkm
    return r


def _df(*rows):
    return pd.DataFrame(rows)


class TestInferGroupCols(unittest.TestCase):
    def test_includes_tool_mode_rank_taxid(self):
        df = _df(_row("S1", "3001", 10.0))
        self.assertEqual(infer_group_cols(df), ["tool", "mode", "rank", "taxid"])

    def test_requires_taxid(self):
        df = pd.DataFrame([{"sample": "S1", "rank": "species", "rpm": 1.0}])
        with self.assertRaises(ValueError):
            infer_group_cols(df)


class TestRpmOnly(unittest.TestCase):
    """No rpkm column -> metric is rpm, behaviour matches the historical filter."""

    def setUp(self):
        # One taxon, three samples. max_rpm = 1000, threshold = 0.5% = 5.0.
        self.df = _df(
            _row("S1", "3001", rpm=1000.0),
            _row("S2", "3001", rpm=3.0),  # below threshold
            _row("S3", "3001", rpm=10.0),  # above threshold
        )

    def test_metric_is_rpm(self):
        out = apply_bleed_filter(self.df)
        self.assertTrue((out["bleed_metric"] == "rpm").all())

    def test_pass_fail_and_threshold(self):
        out = apply_bleed_filter(self.df, fraction=0.005)
        self.assertTrue((out["bleed_max"] == 1000.0).all())
        self.assertTrue((out["bleed_threshold"] == 5.0).all())
        by = out.set_index("sample")["bleed_pass"]
        self.assertTrue(bool(by["S1"]))
        self.assertFalse(bool(by["S2"]))
        self.assertTrue(bool(by["S3"]))

    def test_below_floor_is_not_applied(self):
        df = _df(_row("S1", "9999", rpm=0.5), _row("S2", "9999", rpm=0.2))
        out = apply_bleed_filter(df, rpm_floor=1.0)  # max 0.5 < 1.0
        self.assertFalse(out["bleed_applied"].any())
        self.assertTrue(out["bleed_pass"].all())  # not applied -> all keep


class TestRpkmInvariance(unittest.TestCase):
    """The core finding: with the filter applied under both metrics, rpkm and
    rpm give the SAME bleed_pass, because it is a within-taxon ratio and genome
    length (constant per taxon) cancels."""

    def setUp(self):
        # rpkm = rpm * 0.1 (a single 10 kb genome), constant within the taxon.
        self.rpm_rows = [
            _row("S1", "3001", rpm=1000.0),
            _row("S2", "3001", rpm=3.0),
            _row("S3", "3001", rpm=10.0),
        ]
        self.rpkm_rows = [
            _row("S1", "3001", rpm=1000.0, rpkm=100.0),
            _row("S2", "3001", rpm=3.0, rpkm=0.3),
            _row("S3", "3001", rpm=10.0, rpkm=1.0),
        ]

    def test_bleed_pass_identical_between_metrics(self):
        rpm_out = apply_bleed_filter(_df(*self.rpm_rows))  # no rpkm col -> rpm
        rpkm_out = apply_bleed_filter(_df(*self.rpkm_rows))  # auto -> rpkm
        self.assertTrue((rpkm_out["bleed_metric"] == "rpkm").all())
        self.assertEqual(
            list(rpm_out.sort_values("sample")["bleed_pass"]),
            list(rpkm_out.sort_values("sample")["bleed_pass"]),
        )

    def test_rpkm_threshold_is_rescaled(self):
        out = apply_bleed_filter(_df(*self.rpkm_rows))
        # bleed_max is now the rpkm max (100), threshold 0.5% of that.
        self.assertTrue((out["bleed_max"] == 100.0).all())
        self.assertTrue((out["bleed_threshold"] == 0.5).all())


class TestFloorDivergence(unittest.TestCase):
    """The ONLY place rpkm vs rpm changes the outcome: the floor gate, because
    the two metrics live on different scales."""

    def setUp(self):
        # max rpm = 5 (>= rpm_floor 1.0 -> applied), but rpkm = rpm*0.01 so
        # max rpkm = 0.05 (< rpkm_floor 0.1 -> NOT applied).
        self.rows = [
            _row("S1", "3001", rpm=5.0, rpkm=0.05),
            _row("S2", "3001", rpm=0.01, rpkm=0.0001),  # would fail under rpm
        ]

    def test_rpm_applies_rpkm_does_not(self):
        forced_rpm = apply_bleed_filter(_df(*self.rows), metric="rpm", rpm_floor=1.0)
        self.assertTrue(forced_rpm["bleed_applied"].all())
        # under rpm, S2 (0.01) is below 0.005*5 = 0.025 -> fails
        self.assertFalse(bool(forced_rpm.set_index("sample")["bleed_pass"]["S2"]))

        auto_rpkm = apply_bleed_filter(_df(*self.rows), rpkm_floor=0.1)
        self.assertTrue((auto_rpkm["bleed_metric"] == "rpkm").all())
        self.assertFalse(auto_rpkm["bleed_applied"].any())  # max rpkm below floor
        self.assertTrue(auto_rpkm["bleed_pass"].all())  # not applied -> all keep


class TestMetricSelectionAndErrors(unittest.TestCase):
    def test_auto_falls_back_to_rpm_when_group_rpkm_all_na(self):
        df = _df(
            _row("S1", "3001", rpm=100.0, rpkm=float("nan")),
            _row("S2", "3001", rpm=50.0, rpkm=float("nan")),
        )
        out = apply_bleed_filter(df)
        self.assertTrue((out["bleed_metric"] == "rpm").all())

    def test_per_group_metric_is_independent(self):
        # taxon 3001 has rpkm, taxon 4000 does not (NA) -> different metrics.
        df = _df(
            _row("S1", "3001", rpm=100.0, rpkm=10.0),
            _row("S2", "3001", rpm=1.0, rpkm=0.1),
            _row("S1", "4000", rpm=100.0, rpkm=float("nan")),
            _row("S2", "4000", rpm=1.0, rpkm=float("nan")),
        )
        out = apply_bleed_filter(df)
        m = out.set_index("taxid")["bleed_metric"]
        self.assertEqual(m["3001"].iloc[0] if hasattr(m["3001"], "iloc") else m["3001"], "rpkm")
        self.assertEqual(m["4000"].iloc[0] if hasattr(m["4000"], "iloc") else m["4000"], "rpm")

    def test_force_rpkm_without_column_raises(self):
        df = _df(_row("S1", "3001", rpm=100.0))
        with self.assertRaises(ValueError):
            apply_bleed_filter(df, metric="rpkm")

    def test_bad_metric_raises(self):
        df = _df(_row("S1", "3001", rpm=100.0))
        with self.assertRaises(ValueError):
            apply_bleed_filter(df, metric="cpm")

    def test_missing_rpm_column_raises(self):
        df = pd.DataFrame(
            [{"sample": "S1", "tool": "x", "mode": "y", "rank": "species", "taxid": "3001"}]
        )
        with self.assertRaises(ValueError):
            apply_bleed_filter(df)

    def test_missing_sample_column_raises(self):
        df = pd.DataFrame(
            [{"tool": "x", "mode": "y", "rank": "species", "taxid": "3001", "rpm": 1.0}]
        )
        with self.assertRaises(ValueError):
            apply_bleed_filter(df)


class TestRpkmInvarianceMultiTaxon(unittest.TestCase):
    """Invariance is per-taxon: two taxa with DIFFERENT genome lengths (hence
    different rpm->rpkm scale factors) coexisting in one table still produce the
    same bleed_pass under rpm and rpkm, because the ratio cancels within each
    taxon independently."""

    def setUp(self):
        # taxon 1000: rpkm = rpm * 0.1 (10 kb). taxon 2000: rpkm = rpm * 0.5 (2 kb).
        self.rpm_rows = [
            _row("S1", "1000", rpm=1000.0),
            _row("S2", "1000", rpm=3.0),
            _row("S3", "1000", rpm=10.0),
            _row("S1", "2000", rpm=500.0),
            _row("S2", "2000", rpm=2.0),
            _row("S3", "2000", rpm=400.0),
        ]
        self.rpkm_rows = [
            _row("S1", "1000", rpm=1000.0, rpkm=100.0),
            _row("S2", "1000", rpm=3.0, rpkm=0.3),
            _row("S3", "1000", rpm=10.0, rpkm=1.0),
            _row("S1", "2000", rpm=500.0, rpkm=250.0),
            _row("S2", "2000", rpm=2.0, rpkm=1.0),
            _row("S3", "2000", rpm=400.0, rpkm=200.0),
        ]

    def test_bleed_pass_identical_across_taxa(self):
        rpm_out = apply_bleed_filter(_df(*self.rpm_rows))  # no rpkm col -> rpm
        rpkm_out = apply_bleed_filter(_df(*self.rpkm_rows))  # auto -> rpkm
        self.assertTrue((rpkm_out["bleed_metric"] == "rpkm").all())
        key = ["taxid", "sample"]
        self.assertEqual(
            list(rpm_out.sort_values(key)["bleed_pass"]),
            list(rpkm_out.sort_values(key)["bleed_pass"]),
        )
        # sanity: the two taxa have genuinely different max metrics/thresholds
        maxes = rpkm_out.groupby("taxid")["bleed_max"].first().to_dict()
        self.assertAlmostEqual(maxes["1000"], 100.0)
        self.assertAlmostEqual(maxes["2000"], 250.0)


if __name__ == "__main__":
    unittest.main()
