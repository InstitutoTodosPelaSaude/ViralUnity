"""Tests for viralunity.scripts.python.add_rpkm_to_summary."""

import unittest

import pandas as pd

from viralunity.scripts.python.add_rpkm_to_summary import add_rpkm


def _summary(rows):
    return pd.DataFrame(rows)


def _genome_lengths(rows):
    return pd.DataFrame(rows)


class TestAddRpkm(unittest.TestCase):
    def _gl(self, rank, taxid, length, n=1):
        return {"rank": rank, "taxid": str(taxid), "genome_length_bp": length, "n_genomes": n}

    def _row(self, rank, taxid, rpm, sample="S1", total_reads=1_000_000):
        return {
            "sample": sample,
            "rank": rank,
            "taxid": str(taxid),
            "rpm": rpm,
            "total_reads": total_reads,
        }

    # ---- basic math --------------------------------------------------------

    def test_rpkm_equals_rpm_times_1000_over_length(self):
        df = _summary([self._row("species", "3001", rpm=100.0)])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        # rpkm = rpm * 1000 / genome_length_bp = 100 * 1000 / 10000 = 10.0
        self.assertAlmostEqual(out.loc[0, "rpkm"], 10.0, places=6)

    def test_rpkm_zero_rpm_gives_zero_rpkm(self):
        df = _summary([self._row("species", "3001", rpm=0.0)])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        self.assertAlmostEqual(out.loc[0, "rpkm"], 0.0)

    def test_rpkm_varies_with_genome_length(self):
        """Longer genome → lower RPKM for same RPM."""
        df = _summary([
            self._row("species", "3001", rpm=1000.0),
            self._row("species", "3002", rpm=1000.0),
        ])
        gl = _genome_lengths([
            self._gl("species", "3001", length=1_000),   # short
            self._gl("species", "3002", length=10_000),  # long
        ])
        out = add_rpkm(df, gl)
        rpkm_short = out.loc[out["taxid"] == "3001", "rpkm"].iloc[0]
        rpkm_long  = out.loc[out["taxid"] == "3002", "rpkm"].iloc[0]
        self.assertGreater(rpkm_short, rpkm_long)

    # ---- rank merging -------------------------------------------------------

    def test_genus_and_family_rows_get_lengths_too(self):
        df = _summary([
            self._row("family", "1000", rpm=50.0),
            self._row("genus",  "2000", rpm=80.0),
            self._row("species","3001", rpm=100.0),
        ])
        gl = _genome_lengths([
            self._gl("family",  "1000", length=10_000),
            self._gl("genus",   "2000", length=9_000),
            self._gl("species", "3001", length=8_000),
        ])
        out = add_rpkm(df, gl)
        # All three should have non-NA rpkm
        self.assertFalse(out["rpkm"].isna().any())

    # ---- NA handling --------------------------------------------------------

    def test_missing_genome_length_gives_na_rpkm(self):
        df = _summary([
            self._row("species", "3001", rpm=100.0),
            self._row("species", "9999", rpm=200.0),  # no entry in genome_lengths
        ])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        rpkm_3001 = out.loc[out["taxid"] == "3001", "rpkm"].iloc[0]
        rpkm_9999 = out.loc[out["taxid"] == "9999", "rpkm"].iloc[0]
        self.assertAlmostEqual(rpkm_3001, 10.0, places=6)
        self.assertTrue(pd.isna(rpkm_9999))

    def test_zero_genome_length_gives_na_rpkm(self):
        df = _summary([self._row("species", "3001", rpm=100.0)])
        gl = _genome_lengths([self._gl("species", "3001", length=0)])
        out = add_rpkm(df, gl)
        self.assertTrue(pd.isna(out.loc[0, "rpkm"]))

    def test_empty_genome_lengths_table_all_na(self):
        df = _summary([self._row("species", "3001", rpm=100.0)])
        gl = _genome_lengths([])
        # Empty genome_lengths df needs columns
        gl = pd.DataFrame(columns=["rank", "taxid", "genome_length_bp", "n_genomes"])
        out = add_rpkm(df, gl)
        self.assertTrue(pd.isna(out.loc[0, "rpkm"]))
        self.assertTrue(pd.isna(out.loc[0, "genome_length_bp"]))

    # ---- output schema ------------------------------------------------------

    def test_adds_genome_length_bp_n_genomes_rpkm_columns(self):
        df = _summary([self._row("species", "3001", rpm=100.0)])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000, n=3)])
        out = add_rpkm(df, gl)
        self.assertIn("genome_length_bp", out.columns)
        self.assertIn("n_genomes", out.columns)
        self.assertIn("rpkm", out.columns)
        self.assertEqual(out.loc[0, "n_genomes"], 3)

    def test_original_columns_preserved(self):
        df = _summary([self._row("species", "3001", rpm=100.0)])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        for col in ["sample", "rank", "taxid", "rpm", "total_reads"]:
            self.assertIn(col, out.columns)

    def test_row_count_unchanged(self):
        df = _summary([
            self._row("species", "3001", rpm=100.0),
            self._row("species", "3002", rpm=50.0),
        ])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        self.assertEqual(len(out), len(df))

    # ---- error cases --------------------------------------------------------

    def test_missing_rank_column_raises(self):
        df = pd.DataFrame([{"taxid": "3001", "rpm": 100.0}])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        with self.assertRaises(ValueError):
            add_rpkm(df, gl)

    def test_missing_taxid_column_raises(self):
        df = pd.DataFrame([{"rank": "species", "rpm": 100.0}])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        with self.assertRaises(ValueError):
            add_rpkm(df, gl)

    def test_missing_rpm_column_raises(self):
        df = pd.DataFrame([{"rank": "species", "taxid": "3001"}])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        with self.assertRaises(ValueError):
            add_rpkm(df, gl)

    # ---- multi-sample correctness -------------------------------------------

    def test_rpkm_correct_across_samples_with_different_depths(self):
        """RPKM should be independent of total_reads when computed via RPM."""
        # Sample A has 1M total reads, RPM 100 → RPM already encodes the depth.
        # RPKM = RPM * 1000 / genome_length, so RPKM doesn't change with total_reads.
        df = _summary([
            self._row("species", "3001", rpm=100.0, sample="A", total_reads=1_000_000),
            self._row("species", "3001", rpm=100.0, sample="B", total_reads=2_000_000),
        ])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        rpkm_A = out.loc[out["sample"] == "A", "rpkm"].iloc[0]
        rpkm_B = out.loc[out["sample"] == "B", "rpkm"].iloc[0]
        self.assertAlmostEqual(rpkm_A, rpkm_B, places=6)

    def test_taxid_as_integer_in_summary_still_merges(self):
        """genome_lengths has str taxid; summary taxids cast to str on merge."""
        df = pd.DataFrame([{
            "sample": "S1",
            "rank": "species",
            "taxid": 3001,  # integer
            "rpm": 100.0,
            "total_reads": 1_000_000,
        }])
        gl = _genome_lengths([self._gl("species", "3001", length=10_000)])
        out = add_rpkm(df, gl)
        self.assertAlmostEqual(out.loc[0, "rpkm"], 10.0, places=6)


if __name__ == "__main__":
    unittest.main()
