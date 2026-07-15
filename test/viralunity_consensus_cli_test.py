"""Tests for viralunity consensus CLI (click-based)."""

import os
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from viralunity.viralunity_consensus_cli import consensus


class Test_ConsensusIlluminaCommand(unittest.TestCase):
    """Tests for `viralunity consensus illumina`."""

    def setUp(self):
        self.runner = CliRunner()
        self._required = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def _invoke(self, extra_args=None):
        args = self._required + (extra_args or [])
        with patch("viralunity.viralunity_consensus_cli.consensus_main", return_value=0):
            return self.runner.invoke(consensus, args, catch_exceptions=False)

    def test_required_args_missing_causes_error(self):
        """Missing required args should exit with non-zero code."""
        result = self.runner.invoke(consensus, ["illumina"])
        self.assertNotEqual(result.exit_code, 0)

    def test_required_args_success_with_reference(self):
        result = self._invoke()
        self.assertEqual(result.exit_code, 0, result.output)

    def test_required_args_success_without_reference(self):
        """Reference flags are optional at parse time; validated by core logic."""
        args = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
        ]
        with patch("viralunity.viralunity_consensus_cli.consensus_main", return_value=0):
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)

    def test_required_args_success_with_segmented_reference(self):
        args = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--segmented-reference",
            "S=/path/to/S.fasta",
            "--segmented-reference",
            "L=/path/to/L.fasta",
        ]
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        called_args = mock_main.call_args[0][0]
        self.assertEqual(
            called_args["segmented_reference"],
            {
                "S": "/path/to/S.fasta",
                "L": "/path/to/L.fasta",
            },
        )

    def test_gene_annotation_threads_into_args(self):
        """--gene-annotation lands in args as a plain path."""
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(
                consensus,
                self._required + ["--gene-annotation", "genes.gff3"],
                catch_exceptions=False,
            )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["gene_annotation"], "genes.gff3")

    def test_segmented_gene_annotation_parses_to_dict(self):
        """Repeated --segmented-gene-annotation SEGMENT=PATH parses to a dict."""
        args = self._required + [
            "--segmented-gene-annotation",
            "S=/path/to/S.gff3",
            "--segmented-gene-annotation",
            "L=/path/to/L.gff3",
        ]
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(
            mock_main.call_args[0][0]["segmented_gene_annotation"],
            {"S": "/path/to/S.gff3", "L": "/path/to/L.gff3"},
        )

    def test_default_values_optional_args(self):
        """Check that all optional args have correct defaults for illumina."""
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "illumina")
        self.assertIsNone(args["gene_annotation"])
        self.assertIsNone(args["segmented_gene_annotation"])
        self.assertEqual(args["run_name"], "undefined")
        self.assertIsNone(args["adapters"])
        self.assertEqual(args["trim_head"], 0)
        self.assertEqual(args["trim_tail"], 0)
        self.assertEqual(args["cut_front_mean_quality"], 10)
        self.assertEqual(args["cut_tail_mean_quality"], 10)
        self.assertEqual(args["cut_right_window_size"], 4)
        self.assertEqual(args["cut_right_mean_quality"], 15)
        self.assertEqual(args["af_threshold"], 0.51)
        self.assertEqual(args["af_isnv_threshold"], 0.0)
        self.assertFalse(args["run_isnv"])
        self.assertEqual(args["minimum_coverage"], 20)
        self.assertEqual(args["minimum_read_length"], 50)
        self.assertEqual(args["threads"], 1)
        self.assertEqual(args["threads_total"], 1)
        self.assertFalse(args["create_config_only"])


class Test_ConsensusNanoporeCommand(unittest.TestCase):
    """Tests for `viralunity consensus nanopore`."""

    def setUp(self):
        self.runner = CliRunner()
        self._required = [
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def test_required_args_success(self):
        with patch("viralunity.viralunity_consensus_cli.consensus_main", return_value=0):
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check nanopore-specific defaults."""
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "nanopore")
        self.assertEqual(args["af_threshold"], 0.51)
        self.assertEqual(args["chunk_size"], 10000)
        self.assertEqual(args["clair3_model"], "r1041_e82_400bps_sup_v500")
        self.assertEqual(args["variant_quality"], 20)
        self.assertEqual(args["variant_depth"], 10)
        self.assertEqual(args["minimum_map_quality"], 30)


class Test_ConsensusCondaPrefix(unittest.TestCase):
    """``--conda-prefix`` (and the ``$VIRALUNITY_CONDA_PREFIX`` env var) must
    land in the args dict passed to ``consensus_main`` so the orchestrator
    can forward it to Snakemake."""

    def setUp(self):
        self.runner = CliRunner()
        self._required_illumina = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]
        self._required_nanopore = [
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def _invoke(self, cli_args, env=None):
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(consensus, cli_args, env=env or {}, catch_exceptions=False)
        return result, mock_main

    def test_explicit_conda_prefix_illumina(self):
        result, mock_main = self._invoke(self._required_illumina + ["--conda-prefix", "/tmp/foo"])
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], "/tmp/foo")

    def test_explicit_conda_prefix_nanopore(self):
        result, mock_main = self._invoke(self._required_nanopore + ["--conda-prefix", "/tmp/foo"])
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], "/tmp/foo")

    def test_default_conda_prefix(self):
        # Ensure env var is unset so we hit the fallback.
        env = {"VIRALUNITY_CONDA_PREFIX": ""}
        env.pop("VIRALUNITY_CONDA_PREFIX")
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("VIRALUNITY_CONDA_PREFIX", None)
            result, mock_main = self._invoke(self._required_illumina)
        self.assertEqual(result.exit_code, 0, result.output)
        expected = str(Path.home() / ".cache" / "viralunity" / "conda-envs")
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], expected)

    def test_env_var_conda_prefix(self):
        with patch.dict(os.environ, {"VIRALUNITY_CONDA_PREFIX": "/srv/shared/envs"}):
            result, mock_main = self._invoke(self._required_illumina)
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], "/srv/shared/envs")


if __name__ == "__main__":
    unittest.main()
