"""Tests for the viralunity get-databases CLI command group.

Focus on the DNA-contamination filter that protects ``diamond makedb`` from
the NCBI Datasets bug where some virus protein FASTA records are actually
nucleotide CDS sequences.
"""

import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

import click
import pytest
from click.testing import CliRunner

from viralunity.viralunity_get_databases_cli import (
    _build_diamond_makedb_cmd,
    _build_diamond_prepdb_cmd,
    _build_update_blastdb_cmd,
    _looks_like_dna,
    _reformat_protein_fasta,
    _safe_extract_zip,
    get_databases,
)


class Test_SafeExtractZip(unittest.TestCase):
    """Zip-slip guard for user-overridable archive downloads."""

    def test_rejects_path_traversal_member(self):
        with tempfile.TemporaryDirectory() as tmp:
            zip_path = os.path.join(tmp, "evil.zip")
            dest = os.path.join(tmp, "dest")
            os.makedirs(dest)
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("../escaped.txt", "pwned")
            with zipfile.ZipFile(zip_path) as zf:
                with self.assertRaises(click.ClickException):
                    _safe_extract_zip(zf, dest)
            self.assertFalse(os.path.exists(os.path.join(tmp, "escaped.txt")))

    def test_extracts_safe_member(self):
        with tempfile.TemporaryDirectory() as tmp:
            zip_path = os.path.join(tmp, "ok.zip")
            dest = os.path.join(tmp, "dest")
            os.makedirs(dest)
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("sub/ok.txt", "fine")
            with zipfile.ZipFile(zip_path) as zf:
                _safe_extract_zip(zf, dest)
            self.assertTrue(os.path.exists(os.path.join(dest, "sub", "ok.txt")))


# ---------------------------------------------------------------------------
# _looks_like_dna unit tests
# ---------------------------------------------------------------------------


class TestLooksLikeDna(unittest.TestCase):
    """Verify the DNA/RNA detection heuristic used to drop contaminated records."""

    def test_pure_dna_above_threshold_is_dna(self):
        seq = "ATGGACAACTCAATTGTTGTAGTCAGAGCTACTAAGGCC"
        self.assertTrue(_looks_like_dna(seq))

    def test_pure_rna_above_threshold_is_dna(self):
        seq = "AUGGACAACUCAAUUGUUGUAGUCAGAGCUACU"
        self.assertTrue(_looks_like_dna(seq))

    def test_dna_with_n_is_dna(self):
        seq = "ATGNNNGACAACTCAATTGTTGTAGTCAGAGCTAC"
        self.assertTrue(_looks_like_dna(seq))

    def test_lowercase_dna_is_dna(self):
        seq = "atggacaactcaattgttgtagtcagagctactaag"
        self.assertTrue(_looks_like_dna(seq))

    def test_real_protein_is_not_dna(self):
        seq = "MPKLPRGLRFGADNEILNDFQELWFPDLFIESSDTHPWYTLKGRVLNAHLDDRLPNVGGRQ"
        self.assertFalse(_looks_like_dna(seq))

    def test_short_acgt_sequence_is_kept(self):
        """Tiny ACGT peptides are biologically possible; we err on the safe
        side and only flag sequences at or above ``min_length``."""
        self.assertFalse(_looks_like_dna("ACGTACG"))

    def test_gaps_and_stops_are_ignored(self):
        seq = "ATG-GAC*AAC TCA-ATT-GTT-GTA-GTC-AGA-GCT-ACT"
        self.assertTrue(_looks_like_dna(seq))

    def test_min_length_parameter_is_respected(self):
        seq = "ATGGAC"
        self.assertFalse(_looks_like_dna(seq, min_length=20))
        self.assertTrue(_looks_like_dna(seq, min_length=5))


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


_MIXED_FASTA = textwrap.dedent("""\
    >NC_139268.1
    ATGGACAACTCAATTGTTGTAGTCAGAGCTACTAAGGCCGCCTTTGTGCCAATCAAACCT
    AAATTGGAAGATGAGGTCAACTATCCTCGAGAGTTCTTTGTAGACGGAAGGATTCCTGCG
    >YP_011242554.1 polyprotein [organism=dengue virus type 2]
    MPKLPRGLRFGADNEILNDFQELWFPDLFIESSDTHPWYTLKGRVLNAHLDDRLPNVGGRQ
    IRRTPHRATVPIASSGLRPVTTVQYDPTALSFLLNARVDIRELRRELLD
    >NC_076867.1
    ATGCTATCTGCAGATGCCAGGACACGGTGGAGTGGAGCTAAACAGGACATAGAGACTCTA
    GCAAGAGGGATTAGTGGAGCCGGGAGATCAGAAGAAATCAGTTTAGATATTGAACCAGAA
    >NP_056776.2:1-3391 polyprotein [organism=zika virus]
    MAKLETVTLSNIGKDGKQTLVLNPRGVNPTNGVAALSQAGAVPALEKRVTVSVSQPSRNR
    KNYKVQVKIQNPTACTANGSCDPSVTRQAYADVTFSFTQYSTRT
    >YP_999999.1 hypothetical protein [organism=mystery virus]
    MSKTASSHNSLSAQLRRAANTRIEVEGNLALSIANDLLLAYGQSPFNSEAECISLSPRFD
    """)


def _write_fasta(tmpdir: Path, name: str, content: str) -> Path:
    path = tmpdir / name
    path.write_text(content)
    return path


# ---------------------------------------------------------------------------
# _reformat_protein_fasta tests
# ---------------------------------------------------------------------------


class TestReformatProteinFasta(unittest.TestCase):
    """Verify DNA records are dropped during the reformat step."""

    def test_drops_dna_records_and_maps_proteins(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "protein.faa", _MIXED_FASTA)
            output_fa = tmpdir / "viral.protein.faa"

            org2taxid = {
                "dengue virus type 2": "11060",
                "zika virus": "64320",
                "mystery virus": "99999",
            }

            taxid_map = _reformat_protein_fasta(
                [input_fa],
                output_fa,
                org2taxid,
            )

            self.assertEqual(
                taxid_map,
                {
                    "YP_011242554.1": "11060",
                    "NP_056776.2": "64320",
                    "YP_999999.1": "99999",
                },
            )

            written = output_fa.read_text()
            self.assertIn(">YP_011242554.1", written)
            self.assertIn(">NP_056776.2", written)
            self.assertIn(">YP_999999.1", written)
            self.assertNotIn("NC_139268.1", written)
            self.assertNotIn("NC_076867.1", written)
            self.assertNotIn("ATGGACAACTCAATTGTTGTAGTCAGAGCTAC", written)

    def test_pure_protein_input_writes_everything(self):
        protein_only = textwrap.dedent("""\
            >YP_1.1 protein A [organism=virus a]
            MPKLPRGLRFGADNEILNDFQ
            >YP_2.1 protein B [organism=virus b]
            MSKTASSHNSLSAQLRRAANT
            """)
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "protein.faa", protein_only)
            output_fa = tmpdir / "out.faa"

            taxid_map = _reformat_protein_fasta(
                [input_fa],
                output_fa,
                {"virus a": "1", "virus b": "2"},
            )

            self.assertEqual(taxid_map, {"YP_1.1": "1", "YP_2.1": "2"})
            written = output_fa.read_text()
            self.assertIn(">YP_1.1", written)
            self.assertIn(">YP_2.1", written)


# ---------------------------------------------------------------------------
# clean-protein-fasta CLI command tests
# ---------------------------------------------------------------------------


class TestCleanProteinFastaCommand(unittest.TestCase):
    """End-to-end tests for ``viralunity get-databases clean-protein-fasta``."""

    def setUp(self):
        self.runner = CliRunner()

    def test_clean_writes_separate_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", _MIXED_FASTA)
            output_fa = tmpdir / "viral.protein.cleaned.faa"

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(output_fa),
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertTrue(output_fa.exists())

            cleaned = output_fa.read_text()
            self.assertIn(">YP_011242554.1", cleaned)
            self.assertNotIn(">NC_139268.1", cleaned)
            self.assertNotIn(">NC_076867.1", cleaned)
            self.assertIn("Kept 3 protein record(s).", result.output)
            self.assertIn("Dropped 2 DNA/RNA record(s)", result.output)

            self.assertEqual(input_fa.read_text(), _MIXED_FASTA)

    def test_clean_in_place_with_backup(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", _MIXED_FASTA)

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(input_fa),
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)

            backup = input_fa.with_suffix(input_fa.suffix + ".with_dna.bak")
            self.assertTrue(backup.exists(), "backup file should be created")
            self.assertEqual(backup.read_text(), _MIXED_FASTA)

            cleaned = input_fa.read_text()
            self.assertIn(">YP_011242554.1", cleaned)
            self.assertNotIn(">NC_139268.1", cleaned)

    def test_clean_in_place_no_backup(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", _MIXED_FASTA)

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(input_fa),
                    "--no-backup",
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)

            backup = input_fa.with_suffix(input_fa.suffix + ".with_dna.bak")
            self.assertFalse(backup.exists(), "backup must not be created with --no-backup")

            cleaned = input_fa.read_text()
            self.assertNotIn(">NC_139268.1", cleaned)
            self.assertIn(">YP_011242554.1", cleaned)

    def test_clean_protein_only_input_is_a_noop(self):
        protein_only = textwrap.dedent("""\
            >YP_1.1 protein A
            MPKLPRGLRFGADNEILNDFQ
            >YP_2.1 protein B
            MSKTASSHNSLSAQLRRAANT
            """)
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", protein_only)
            output_fa = tmpdir / "out.faa"

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(output_fa),
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertIn("Kept 2 protein record(s).", result.output)
            self.assertIn("No DNA/RNA records found.", result.output)

    def test_clean_missing_input_fails(self):
        result = self.runner.invoke(
            get_databases,
            [
                "clean-protein-fasta",
                "--input",
                "/nonexistent/path/protein.faa",
                "--output",
                "/tmp/should_not_be_created.faa",
            ],
        )
        self.assertNotEqual(result.exit_code, 0)


# ---------------------------------------------------------------------------
# get-databases nr: command builders (pure) and CLI wiring
# ---------------------------------------------------------------------------


class TestNrCommandBuilders(unittest.TestCase):
    """The argv builders are pure, so assert exact command lines."""

    def test_update_blastdb_cmd_default_single_thread(self):
        self.assertEqual(
            _build_update_blastdb_cmd("nr", "ncbi", 1),
            ["update_blastdb.pl", "--decompress", "--source", "ncbi", "nr"],
        )

    def test_update_blastdb_cmd_multithread_and_source(self):
        self.assertEqual(
            _build_update_blastdb_cmd("nr_viruses", "gcp", 8),
            [
                "update_blastdb.pl",
                "--decompress",
                "--source",
                "gcp",
                "--num_threads",
                "8",
                "nr_viruses",
            ],
        )

    def test_prepdb_cmd(self):
        self.assertEqual(
            _build_diamond_prepdb_cmd("db/diamond/nr"),
            ["diamond", "prepdb", "--db", "db/diamond/nr"],
        )

    def test_makedb_cmd_minimal(self):
        self.assertEqual(
            _build_diamond_makedb_cmd("nr.faa", "db/diamond/nr", 1),
            ["diamond", "makedb", "--in", "nr.faa", "--db", "db/diamond/nr"],
        )

    def test_makedb_cmd_with_taxonomy_and_threads(self):
        self.assertEqual(
            _build_diamond_makedb_cmd(
                "nr.faa", "db/diamond/nr", 4, "map.gz", "nodes.dmp", "names.dmp"
            ),
            [
                "diamond",
                "makedb",
                "--in",
                "nr.faa",
                "--db",
                "db/diamond/nr",
                "--taxonmap",
                "map.gz",
                "--taxonnodes",
                "nodes.dmp",
                "--taxonnames",
                "names.dmp",
                "--threads",
                "4",
            ],
        )


def _fake_run_command(cmd, cwd=None, timeout=None):
    """Simulate run_command by creating the artifacts each tool would produce.

    Lets the real get_nr control flow (including the _verify_blastdb and
    diamond-makedb output checks) run without invoking any external tool.
    """
    cmd = [str(c) for c in cmd]
    if cmd and cmd[0] == "update_blastdb.pl":
        Path(cwd, f"{cmd[-1]}.pal").touch()  # BLAST+ multi-volume alias file
    elif cmd[:2] == ["diamond", "makedb"]:
        Path(cmd[cmd.index("--db") + 1]).touch()
    # 'diamond prepdb' needs no simulated output.


class TestNrCli(unittest.TestCase):
    """CLI paths for 'get-databases nr' with run_command stubbed (no network)."""

    def setUp(self):
        self.runner = CliRunner()

    @patch(
        "viralunity.viralunity_get_databases_cli.run_command",
        side_effect=_fake_run_command,
    )
    def test_download_then_prepdb(self, mock_run):
        with self.runner.isolated_filesystem():
            result = self.runner.invoke(
                get_databases, ["nr", "--path", "db", "--source", "gcp", "--threads", "4"]
            )
            self.assertEqual(result.exit_code, 0, result.output)
            calls = [[str(c) for c in call.args[0]] for call in mock_run.call_args_list]
            self.assertEqual(calls[0][0], "update_blastdb.pl")
            self.assertIn("--source", calls[0])
            self.assertEqual(mock_run.call_args_list[0].kwargs["cwd"], "db/diamond")
            self.assertEqual(calls[1], ["diamond", "prepdb", "--db", "db/diamond/nr"])
            self.assertIn("--nr-diamond-database  db/diamond/nr", result.output)

    @patch(
        "viralunity.viralunity_get_databases_cli.run_command",
        side_effect=_fake_run_command,
    )
    def test_skip_prepdb_omits_prepdb(self, mock_run):
        with self.runner.isolated_filesystem():
            result = self.runner.invoke(get_databases, ["nr", "--path", "db", "--skip-prepdb"])
            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(mock_run.call_count, 1)
            self.assertEqual(str(mock_run.call_args_list[0].args[0][0]), "update_blastdb.pl")

    @patch(
        "viralunity.viralunity_get_databases_cli.run_command",
        side_effect=_fake_run_command,
    )
    def test_from_blastdb_skips_download_and_preps(self, mock_run):
        with self.runner.isolated_filesystem():
            Path("existing_nr.pal").touch()  # simulate an already-downloaded db
            result = self.runner.invoke(
                get_databases, ["nr", "--path", "db", "--from-blastdb", "existing_nr"]
            )
            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(mock_run.call_count, 1)
            self.assertEqual(
                [str(c) for c in mock_run.call_args_list[0].args[0]],
                ["diamond", "prepdb", "--db", "existing_nr"],
            )
            self.assertIn("--nr-diamond-database  existing_nr", result.output)

    @patch(
        "viralunity.viralunity_get_databases_cli.run_command",
        side_effect=_fake_run_command,
    )
    def test_from_fasta_calls_makedb_with_taxonomy(self, mock_run):
        with self.runner.isolated_filesystem():
            Path("nr.faa").write_text(">p\nMKTLLILAV\n")
            Path("map.tsv").write_text("accession\taccession.version\ttaxid\tgi\n")
            Path("nodes.dmp").touch()
            Path("names.dmp").touch()
            result = self.runner.invoke(
                get_databases,
                [
                    "nr",
                    "--path",
                    "db",
                    "--from-fasta",
                    "nr.faa",
                    "--taxonmap",
                    "map.tsv",
                    "--taxonnodes",
                    "nodes.dmp",
                    "--taxonnames",
                    "names.dmp",
                ],
            )
            self.assertEqual(result.exit_code, 0, result.output)
            cmd = [str(c) for c in mock_run.call_args_list[0].args[0]]
            self.assertEqual(cmd[:2], ["diamond", "makedb"])
            self.assertIn("--taxonmap", cmd)
            self.assertIn("--nr-diamond-database  db/diamond/nr.dmnd", result.output)

    @patch(
        "viralunity.viralunity_get_databases_cli.run_command",
        side_effect=_fake_run_command,
    )
    def test_mutually_exclusive_byo_flags(self, mock_run):
        with self.runner.isolated_filesystem():
            Path("nr.faa").write_text(">p\nMKTLLILAV\n")
            result = self.runner.invoke(
                get_databases,
                ["nr", "--from-blastdb", "some_prefix", "--from-fasta", "nr.faa"],
            )
            self.assertNotEqual(result.exit_code, 0)
            self.assertIn("mutually exclusive", result.output)
            mock_run.assert_not_called()

    @patch(
        "viralunity.viralunity_get_databases_cli.run_command",
        side_effect=_fake_run_command,
    )
    def test_taxon_flags_require_from_fasta(self, mock_run):
        with self.runner.isolated_filesystem():
            Path("map.tsv").touch()
            result = self.runner.invoke(get_databases, ["nr", "--taxonmap", "map.tsv"])
            self.assertNotEqual(result.exit_code, 0)
            self.assertIn("only apply with --from-fasta", result.output)
            mock_run.assert_not_called()


@pytest.mark.empirical
class TestNrFromFastaEmpirical(unittest.TestCase):
    """Exercise the real --from-fasta path (diamond makedb + staxids).

    Proves the makedb/taxonomy machinery end-to-end on a 2-sequence FASTA
    without touching the real 100+ GB nr database. Opt-in (empirical marker)
    and skipped where diamond is unavailable.
    """

    @unittest.skipIf(shutil.which("diamond") is None, "diamond not on PATH")
    def test_build_dmnd_and_recover_staxids(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            Path("nr.faa").write_text(
                ">ABC12345.1 p1\nMKTLLILAVVAAALADQAAAAAAAAAAAA\n"
                ">DEF67890.1 p2\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n"
            )
            Path("map.tsv").write_text(
                "accession\taccession.version\ttaxid\tgi\n"
                "ABC12345\tABC12345.1\t11320\t0\n"
                "DEF67890\tDEF67890.1\t10298\t0\n"
            )
            # NCBI .dmp format: fields joined by '\t|\t', line ends with '\t|'.
            Path("nodes.dmp").write_text(
                "1\t|\t1\t|\tno rank\t|\t\t|\n"
                "11320\t|\t1\t|\tspecies\t|\t\t|\n"
                "10298\t|\t1\t|\tspecies\t|\t\t|\n"
            )
            Path("names.dmp").write_text(
                "1\t|\troot\t|\t\t|\tscientific name\t|\n"
                "11320\t|\tInfluenza A virus\t|\t\t|\tscientific name\t|\n"
                "10298\t|\tHuman herpesvirus 1\t|\t\t|\tscientific name\t|\n"
            )
            result = runner.invoke(
                get_databases,
                [
                    "nr",
                    "--path",
                    "db",
                    "--from-fasta",
                    "nr.faa",
                    "--taxonmap",
                    "map.tsv",
                    "--taxonnodes",
                    "nodes.dmp",
                    "--taxonnames",
                    "names.dmp",
                ],
            )
            self.assertEqual(result.exit_code, 0, result.output)
            dmnd = Path("db/diamond/nr.dmnd")
            self.assertTrue(dmnd.is_file())

            # Back-translate the (non-repetitive) P2 protein so the query hits
            # it, then confirm the hit carries the embedded staxids.
            codon = {
                "M": "ATG",
                "S": "AGC",
                "T": "ACC",
                "N": "AAC",
                "P": "CCG",
                "K": "AAA",
                "Q": "CAG",
                "R": "CGT",
                "D": "GAT",
                "V": "GTG",
                "F": "TTT",
                "G": "GGC",
            }
            peptide = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"
            Path("q.fna").write_text(">q1\n" + "".join(codon[a] for a in peptide) + "\n")
            out = subprocess.run(
                [
                    "diamond",
                    "blastx",
                    "--db",
                    str(dmnd),
                    "--query",
                    "q.fna",
                    "--outfmt",
                    "6",
                    "qseqid",
                    "sseqid",
                    "staxids",
                    "--quiet",
                ],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
            self.assertTrue(out.strip(), "expected at least one diamond hit")
            self.assertIn("10298", out)


if __name__ == "__main__":
    unittest.main()
