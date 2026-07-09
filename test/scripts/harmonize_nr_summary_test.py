"""Tests for viralunity.scripts.python.harmonize_nr_summary."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.harmonize_nr_summary import (
    harmonize,
    load_nr_verdicts,
    majority_bool,
)

# Taxonomy: species 3001 (viral, under family 1000) and 3100 (bacterial phage sp).
NODES = [
    ("1", "1", "no rank"),
    ("1000", "1", "family"),
    ("3001", "1000", "species"),
    ("1100", "1", "family"),
    ("3100", "1100", "species"),
]
NAMES = [
    ("1", "root"),
    ("1000", "Orthomyxoviridae"),
    ("3001", "Influenza A virus"),
    ("1100", "Steitzviridae"),
    ("3100", "Lambdavirus DE3"),
]


def _write_taxdump(d):
    nodes = os.path.join(d, "nodes.dmp")
    names = os.path.join(d, "names.dmp")
    with open(nodes, "w") as f:
        for t, p, r in NODES:
            f.write(f"{t}\t|\t{p}\t|\t{r}\t|\t-\t|\n")
    with open(names, "w") as f:
        for t, n in NAMES:
            f.write(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n")
    return d


class TestMajorityBool(unittest.TestCase):
    def test_true_false_ambiguous_none(self):
        self.assertIs(majority_bool(2, 3), True)
        self.assertIs(majority_bool(1, 3), False)
        self.assertEqual(majority_bool(1, 2), "ambiguous")
        self.assertIsNone(majority_bool(0, 0))


class TestLoadNrVerdicts(unittest.TestCase):
    def test_two_part_key_and_viral_flag(self):
        with tempfile.TemporaryDirectory() as d:
            nr = os.path.join(d, "nr.tsv")
            header = [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "domain",
                "kingdom",
                "realm",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "consensus_rank",
                "consensus_taxon",
            ]
            with open(nr, "w") as f:
                f.write("\t".join(header) + "\n")
                row = (
                    ["sample-A|k141_1"]
                    + ["x"] * 11
                    + [
                        "Viruses",
                        "NA",
                        "Riboviria",
                        "Negarnaviricota",
                        "NA",
                        "NA",
                        "Orthomyxoviridae",
                        "NA",
                        "Influenza A virus",
                        "species",
                        "Influenza A virus",
                    ]
                )
                f.write("\t".join(row) + "\n")
            verdicts = load_nr_verdicts(nr)
            self.assertIn(("sample-A", "k141_1"), verdicts)
            v = verdicts[("sample-A", "k141_1")]
            self.assertTrue(v["is_viral"])
            self.assertEqual(v["species"], "Influenza A virus")


class TestHarmonize(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name
        _write_taxdump(self.d)
        # summary: one viral species row (3001) and one bacterial species row (3100)
        pd.DataFrame(
            [
                ("sample-A", "diamond", "contigs", "species", "3001", "Influenza A virus"),
                ("sample-A", "diamond", "contigs", "species", "3100", "Lambdavirus DE3"),
                ("sample-A", "diamond", "contigs", "family", "1000", "Orthomyxoviridae"),
            ],
            columns=["sample", "tool", "mode", "rank", "taxid", "name"],
        ).to_csv(os.path.join(self.d, "summary.tsv"), sep="\t", index=False)
        # krona map: contig -> refseq taxid (climbs to the species)
        with open(os.path.join(self.d, "sample-A.krona.txt"), "w") as f:
            f.write("k141_1\t3001\n")  # -> Influenza (viral per NR)
            f.write("k141_2\t3100\n")  # -> Lambdavirus (NOT viral per NR)
        # NR verdict table: k141_1 viral, k141_2 non-viral (bacterial phylum)
        header = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "domain",
            "kingdom",
            "realm",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "consensus_rank",
            "consensus_taxon",
        ]
        with open(os.path.join(self.d, "nr.tsv"), "w") as f:
            f.write("\t".join(header) + "\n")
            f.write(
                "\t".join(
                    ["sample-A|k141_1"]
                    + ["x"] * 11
                    + [
                        "Viruses",
                        "NA",
                        "Riboviria",
                        "Negarnaviricota",
                        "NA",
                        "NA",
                        "Orthomyxoviridae",
                        "NA",
                        "Influenza A virus",
                        "species",
                        "Influenza A virus",
                    ]
                )
                + "\n"
            )
            f.write(
                "\t".join(
                    ["sample-A|k141_2"]
                    + ["x"] * 11
                    + [
                        "Bacteria",
                        "NA",
                        "NA",
                        "Pseudomonadota",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "Escherichia coli",
                        "species",
                        "Escherichia coli",
                    ]
                )
                + "\n"
            )

    def test_drops_nonviral_species_appends_nr_pass_writes_flags(self):
        out = os.path.join(self.d, "out.tsv")
        dropped = os.path.join(self.d, "dropped.tsv")
        flags = os.path.join(self.d, "flags.tsv")
        harmonize(
            summary_path=os.path.join(self.d, "summary.tsv"),
            nr_table_path=os.path.join(self.d, "nr.tsv"),
            krona_files=[("sample-A", os.path.join(self.d, "sample-A.krona.txt"))],
            output_path=out,
            dropped_path=dropped,
            flags_path=flags,
            taxdump_dir=self.d,
        )
        kept = pd.read_csv(out, sep="\t", dtype=str)
        # bacterial species 3100 removed; viral 3001 and family row kept
        self.assertEqual(set(kept["taxid"]), {"3001", "1000"})
        self.assertIn("nr_pass", kept.columns)
        # viral species row nr_pass True; family row NA (kept)
        flu = kept[kept["taxid"] == "3001"].iloc[0]
        self.assertEqual(flu["nr_is_virus"], "True")
        # dropped sidecar holds the bacterial species
        drp = pd.read_csv(dropped, sep="\t", dtype=str)
        self.assertEqual(set(drp["taxid"]), {"3100"})
        # flags file exists with a header
        self.assertTrue(os.path.exists(flags))


if __name__ == "__main__":
    unittest.main()
