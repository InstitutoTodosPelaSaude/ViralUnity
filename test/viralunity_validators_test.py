"""Tests for metagenomics argument validation (viralunity.validators).

Focuses on the cross-dependency checks that are easy to get wrong from the CLI,
in particular the RPKM (``--viral-genomes``) → ``--viral-taxids`` requirement,
which is independent of reference assembly.
"""

import csv
import os
import tempfile
import unittest

from viralunity.exceptions import SampleSheetError, ValidationError, ViralUnityFileNotFoundError
from viralunity.validators import (
    ensure_within_base,
    sanitize_identifier,
    validate_metagenomics_requirements,
    validate_sample_sheet,
)


def _touch(path: str) -> str:
    with open(path, "w") as fh:
        fh.write("")
    return path


class Test_RPKM_Validation(unittest.TestCase):
    """compute_rpkm is derived from --viral-genomes; its deps must be validated."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        # A classification-free config so we reach the viral_genomes check without
        # needing real kraken2/krona/taxdump databases.
        self.base_args = {
            "run_denovo_assembly": False,
            "run_kraken2_reads": False,
            "run_kraken2_contigs": False,
            "run_diamond_reads": False,
            "run_diamond_contigs": False,
            "run_reference_assembly": False,
        }
        self.genomes = _touch(os.path.join(self.tmp, "viral.genomes.fasta"))
        self.g2t = _touch(os.path.join(self.tmp, "genome2taxid.tsv"))

    def tearDown(self):
        self._tmp.cleanup()

    def test_no_viral_genomes_is_ok(self):
        """No --viral-genomes → RPKM off → no extra requirements."""
        validate_metagenomics_requirements({**self.base_args, "viral_genomes": "NA"})

    def test_viral_genomes_without_taxids_raises(self):
        """--viral-genomes without --viral-taxids must fail up-front, not at runtime."""
        args = {**self.base_args, "viral_genomes": self.genomes, "viral_taxids": "NA"}
        with self.assertRaises(ValidationError):
            validate_metagenomics_requirements(args)

    def test_viral_genomes_missing_file_raises(self):
        args = {
            **self.base_args,
            "viral_genomes": os.path.join(self.tmp, "does_not_exist.fasta"),
            "viral_taxids": self.g2t,
        }
        with self.assertRaises(ViralUnityFileNotFoundError):
            validate_metagenomics_requirements(args)

    def test_viral_taxids_missing_file_raises(self):
        args = {
            **self.base_args,
            "viral_genomes": self.genomes,
            "viral_taxids": os.path.join(self.tmp, "missing_g2t.tsv"),
        }
        with self.assertRaises(ViralUnityFileNotFoundError):
            validate_metagenomics_requirements(args)

    def test_viral_genomes_with_taxids_is_ok(self):
        args = {**self.base_args, "viral_genomes": self.genomes, "viral_taxids": self.g2t}
        validate_metagenomics_requirements(args)


class Test_SampleSheetIntegrity(unittest.TestCase):
    """Guardrails against silent sample-sheet data loss / corruption."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.r1 = _touch(os.path.join(self.tmp, "s_R1.fastq.gz"))
        self.r2 = _touch(os.path.join(self.tmp, "s_R2.fastq.gz"))
        self.np = _touch(os.path.join(self.tmp, "s.fastq.gz"))

    def tearDown(self):
        self._tmp.cleanup()

    def _sheet(self, rows):
        path = os.path.join(self.tmp, "sheet.csv")
        with open(path, "w", newline="") as fh:
            csv.writer(fh).writerows(rows)
        return path

    def test_duplicate_sample_ids_rejected(self):
        """Two rows with the same id must error, not silently collapse to one."""
        sheet = self._sheet([["dup", self.r1, self.r2], ["dup", self.r1, self.r2]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "illumina")

    def test_ragged_row_rejected_not_crash(self):
        """A short/ragged row raises SampleSheetError (not a TypeError from NaN)."""
        sheet = self._sheet([["s1", self.r1, self.r2], ["s2", self.r1]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "illumina")

    def test_extra_column_rejected(self):
        """A nanopore sheet with an Illumina-style 3rd column is rejected."""
        sheet = self._sheet([["s1", self.np, "extra"]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "nanopore")

    def test_trailing_comma_tolerated(self):
        """A stray trailing empty field is tolerated, not treated as a real column."""
        sheet = self._sheet([["s1", self.r1, self.r2, ""]])
        samples = validate_sample_sheet(sheet, "illumina")
        self.assertEqual(samples, {"s1": [self.r1, self.r2]})

    def test_sample_named_NA_kept_as_string(self):
        """A sample literally named 'NA' stays a string (pandas would coerce to NaN)."""
        sheet = self._sheet([["NA", self.np]])
        samples = validate_sample_sheet(sheet, "nanopore")
        self.assertEqual(list(samples.keys()), ["NA"])

    def test_numeric_sample_name_kept_as_string(self):
        """A numeric-looking id like '001' is not coerced to an int."""
        sheet = self._sheet([["001", self.np]])
        samples = validate_sample_sheet(sheet, "nanopore")
        self.assertEqual(list(samples.keys()), ["001"])

    def test_empty_sample_name_rejected(self):
        sheet = self._sheet([["", self.np]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "nanopore")

    def test_unsafe_sample_name_rejected(self):
        """A sample id with a path separator must be rejected (injection surface)."""
        sheet = self._sheet([["../evil", self.np]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "nanopore")


class Test_Sanitization(unittest.TestCase):
    """Untrusted-input guards for run names and output paths."""

    def test_sanitize_identifier_accepts_safe(self):
        self.assertEqual(sanitize_identifier("run_2026.01-A"), "run_2026.01-A")

    def test_sanitize_identifier_rejects_path_separator(self):
        for bad in ["../etc", "a/b", "a\\b", "..", ".", "", "  ", "a b", "a;rm -rf"]:
            with self.assertRaises(ValidationError):
                sanitize_identifier(bad, field="run name")

    def test_ensure_within_base_allows_child(self):
        with tempfile.TemporaryDirectory() as base:
            target = ensure_within_base("runs/x", base)
            self.assertTrue(target.startswith(os.path.abspath(base)))

    def test_ensure_within_base_rejects_traversal(self):
        with tempfile.TemporaryDirectory() as base:
            with self.assertRaises(ValidationError):
                ensure_within_base("../../etc/passwd", base)


if __name__ == "__main__":
    unittest.main()
