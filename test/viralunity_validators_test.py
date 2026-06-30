"""Tests for metagenomics argument validation (viralunity.validators).

Focuses on the cross-dependency checks that are easy to get wrong from the CLI,
in particular the RPKM (``--viral-genomes``) → ``--viral-taxids`` requirement,
which is independent of reference assembly.
"""

import os
import tempfile
import unittest

from viralunity.exceptions import ValidationError, ViralUnityFileNotFoundError
from viralunity.validators import validate_metagenomics_requirements


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


if __name__ == "__main__":
    unittest.main()
