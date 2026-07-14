"""Tests for viralunity.scripts.python.harmonize_nr_summary."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.harmonize_nr_summary import (
    NA,
    _add_final_species,
    _build_flags,
    aggregate_species_row,
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
        # final_species present, immediately right of nr_correct_species
        cols = list(kept.columns)
        self.assertIn("final_species", cols)
        self.assertEqual(cols[cols.index("nr_correct_species") + 1], "final_species")
        # NR agreed on the flu species -> final_species falls back to name
        self.assertEqual(flu["final_species"], "Influenza A virus")
        # family row (no NR correction) -> final_species echoes the family name
        fam = kept[kept["taxid"] == "1000"].iloc[0]
        self.assertEqual(fam["final_species"], "Orthomyxoviridae")


class TestFinalSpecies(unittest.TestCase):
    """_add_final_species: coalesce(nr_correct_species, name), positioned + all rows."""

    def _df(self):
        return pd.DataFrame(
            [
                # NR disagreed -> nr_correct_species carries the correction
                ("species", "Influenza B virus", "Influenza A virus"),
                # NR agreed (nr_correct_species NA) -> keep original name
                ("species", "Enterovirus C", NA),
                # family row -> echoes its own name
                ("family", "Picornaviridae", NA),
                # empty-string correction is treated as absent
                ("species", "Mayaro virus", ""),
            ],
            columns=["rank", "name", "nr_correct_species"],
        )

    def test_coalesce_and_positioning(self):
        out = _add_final_species(self._df())
        cols = list(out.columns)
        self.assertEqual(cols[cols.index("nr_correct_species") + 1], "final_species")
        vals = list(out["final_species"])
        self.assertEqual(
            vals,
            ["Influenza A virus", "Enterovirus C", "Picornaviridae", "Mayaro virus"],
        )

    def test_does_not_mutate_input(self):
        df = self._df()
        _add_final_species(df)
        self.assertNotIn("final_species", df.columns)


class TestAggregateSpeciesRow(unittest.TestCase):
    """The species-agreement path inside aggregate_species_row."""

    def _viral(self, species):
        return {"is_viral": True, "species": species}

    def test_species_mismatch_reports_correct_species(self):
        # NR agrees the taxon is viral, but the majority of hits name a
        # different species than the RefSeq row -> nr_species_correct False and
        # the most-common NR species surfaced as nr_correct_species.
        verdicts = {
            ("s", "c1"): self._viral("Influenza B virus"),
            ("s", "c2"): self._viral("Influenza B virus"),
        }
        agg = aggregate_species_row("s", "3001", "Influenza A virus", ["c1", "c2"], verdicts)
        self.assertEqual(agg["nr_is_virus"], "True")
        self.assertEqual(agg["nr_species_correct"], "False")
        self.assertEqual(agg["nr_correct_species"], "Influenza B virus")

    def test_species_match_is_correct(self):
        verdicts = {("s", "c1"): self._viral("Influenza A virus")}
        agg = aggregate_species_row("s", "3001", "Influenza A virus", ["c1"], verdicts)
        self.assertEqual(agg["nr_species_correct"], "True")
        self.assertEqual(agg["nr_correct_species"], NA)

    def test_no_nr_data_is_na(self):
        agg = aggregate_species_row("s", "3001", "Influenza A virus", [], {})
        self.assertEqual(agg["nr_is_virus"], NA)
        self.assertEqual(agg["nr_species_correct"], NA)


class TestBuildFlags(unittest.TestCase):
    """RefSeq-vs-NR species disagreement flags (misid_novel / misid_known)."""

    def _rec(self, species):
        return {
            "is_viral": True,
            "species": species,
            "phylum": "Negarnaviricota",
            "pident": "90.0",
            "evalue": "1e-9",
            "bitscore": "200",
        }

    def test_misid_novel_when_nr_species_not_in_sample(self):
        verdicts = {("s", "c1"): self._rec("Influenza B virus")}
        contig_info = {("s", "c1"): {"hit_taxid": "3001", "sp_name": "Influenza A virus"}}
        rows = _build_flags(verdicts, contig_info, {"s": {"influenza a virus"}})
        self.assertEqual(len(rows), 1)
        # columns: sample, contig, reason, refseq_species, nr_species, in_sample, ...
        self.assertEqual(rows[0][2], "misid_novel")
        self.assertEqual(rows[0][5], "False")

    def test_misid_known_when_nr_species_in_sample(self):
        verdicts = {("s", "c1"): self._rec("Influenza B virus")}
        contig_info = {("s", "c1"): {"hit_taxid": "3001", "sp_name": "Influenza A virus"}}
        rows = _build_flags(verdicts, contig_info, {"s": {"influenza b virus"}})
        self.assertEqual(rows[0][2], "misid_known")
        self.assertEqual(rows[0][5], "True")

    def test_no_flag_when_nr_agrees_with_refseq(self):
        verdicts = {("s", "c1"): self._rec("Influenza A virus")}
        contig_info = {("s", "c1"): {"hit_taxid": "3001", "sp_name": "Influenza A virus"}}
        self.assertEqual(_build_flags(verdicts, contig_info, {"s": set()}), [])


class TestHarmonizeInputContracts(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name
        _write_taxdump(self.d)

    def _run(self, summary_path):
        out = os.path.join(self.d, "out.tsv")
        dropped = os.path.join(self.d, "dropped.tsv")
        flags = os.path.join(self.d, "flags.tsv")
        result = harmonize(
            summary_path=summary_path,
            nr_table_path=os.path.join(self.d, "does_not_exist.tsv"),
            krona_files=[],
            output_path=out,
            dropped_path=dropped,
            flags_path=flags,
            taxdump_dir=self.d,
        )
        return result, out, dropped, flags

    def test_empty_input_writes_empty_outputs_and_flags_header(self):
        summary = os.path.join(self.d, "summary.tsv")
        open(summary, "w").close()  # 0-byte input
        (kept, dropped_n), out, dropped, flags = self._run(summary)
        self.assertEqual((kept, dropped_n), (0, 0))
        self.assertEqual(os.path.getsize(out), 0)
        self.assertEqual(os.path.getsize(dropped), 0)
        # flags always carries its fixed schema header
        with open(flags) as fh:
            self.assertTrue(fh.readline().startswith("sample\tcontig\treason"))

    def test_missing_required_column_raises(self):
        summary = os.path.join(self.d, "summary.tsv")
        # no 'taxid' column
        pd.DataFrame(
            [("sample-A", "species", "Influenza A virus")],
            columns=["sample", "rank", "name"],
        ).to_csv(summary, sep="\t", index=False)
        with self.assertRaises(ValueError) as ctx:
            self._run(summary)
        self.assertIn("taxid", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
