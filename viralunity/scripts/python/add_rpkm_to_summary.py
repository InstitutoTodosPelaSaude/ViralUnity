#!/usr/bin/env python3
"""Add RPKM (Reads Per Kilobase per Million reads) to a ViralUnity taxa summary.

RPKM is computed from the already-present ``rpm`` column by dividing by the
representative genome length in kilobases:

    rpkm = rpm * 1000 / genome_length_bp

where ``genome_length_bp`` is the median genome length at each (rank, taxid)
node from the RefSeq viral genome database, as produced by
``build_genome_length_table.py``.

Rows for which no genome length is available receive ``genome_length_bp = NA``
and ``rpkm = NA``.  Higher-rank (family, genus) RPKM values are approximate
because they use the median genome length across all RefSeq sequences under
that node.

Input  : *_taxa_summary_RPM.tsv   (must contain rpm, rank, taxid columns)
         genome_lengths.tsv       (rank, taxid, name, genome_length_bp, n_genomes)
Output : *_taxa_summary_RPKM.tsv  (adds genome_length_bp, n_genomes, rpkm columns)
"""

import argparse
import sys

import pandas as pd


def add_rpkm(
    df: pd.DataFrame,
    genome_lengths: pd.DataFrame,
    rpm_col: str = "rpm",
    rpkm_col: str = "rpkm",
) -> pd.DataFrame:
    """Merge genome lengths and compute RPKM.

    Parameters
    ----------
    df:
        Input summary TSV loaded as a DataFrame.  Must contain ``rank``,
        ``taxid``, and *rpm_col* columns.
    genome_lengths:
        Genome-length table with columns ``rank``, ``taxid``,
        ``genome_length_bp``, ``n_genomes``.
    rpm_col:
        Name of the RPM column to use as the normalisation base.
    rpkm_col:
        Name of the RPKM column to create.

    Returns
    -------
    DataFrame with ``genome_length_bp``, ``n_genomes``, and *rpkm_col* added.
    Rows with no matching genome length receive NA in all three columns.
    """
    if "rank" not in df.columns or "taxid" not in df.columns:
        raise ValueError("Input must contain 'rank' and 'taxid' columns.")
    if rpm_col not in df.columns:
        raise ValueError(f"Input missing RPM column '{rpm_col}'. Run add_RPM_to_summary.py first.")

    gl = genome_lengths[["rank", "taxid", "genome_length_bp", "n_genomes"]].copy()
    gl["taxid"] = gl["taxid"].astype(str)

    out = df.copy()
    out["taxid"] = out["taxid"].astype(str)

    out = out.merge(gl, on=["rank", "taxid"], how="left")

    # rpkm = rpm * 1000 / genome_length_bp
    out[rpkm_col] = pd.NA
    valid = out["genome_length_bp"].notna() & (out["genome_length_bp"].astype(float) > 0)
    out.loc[valid, rpkm_col] = (
        out.loc[valid, rpm_col].astype(float)
        * 1000.0
        / out.loc[valid, "genome_length_bp"].astype(float)
    )

    n_with_length = valid.sum()
    n_total = len(out)
    sys.stderr.write(
        f"[add_rpkm] {n_with_length}/{n_total} rows have RPKM "
        f"({n_total - n_with_length} without genome length -> NA).\n"
    )

    return out


def run_cli() -> None:
    ap = argparse.ArgumentParser(
        description="Add RPKM to a ViralUnity taxa summary TSV using per-taxon genome lengths."
    )
    ap.add_argument("--summary", required=True, help="Input *_RPM.tsv taxa summary.")
    ap.add_argument(
        "--genome-lengths",
        required=True,
        help="Genome-length table produced by build_genome_length_table.py.",
    )
    ap.add_argument("--out", required=True, help="Output *_RPKM.tsv path.")
    ap.add_argument("--rpm-col", default="rpm", help="Name of RPM column (default: rpm).")
    ap.add_argument(
        "--rpkm-col", default="rpkm", help="Name of RPKM column to create (default: rpkm)."
    )
    args = ap.parse_args()

    df = pd.read_csv(args.summary, sep="\t", dtype={"taxid": str})
    gl = pd.read_csv(args.genome_lengths, sep="\t", dtype={"taxid": str})
    out = add_rpkm(df, gl, rpm_col=args.rpm_col, rpkm_col=args.rpkm_col)
    out.to_csv(args.out, sep="\t", index=False, na_rep="NA")


def run_snakemake() -> None:
    df = pd.read_csv(str(snakemake.input.summary), sep="\t", dtype={"taxid": str})
    gl = pd.read_csv(str(snakemake.input.genome_lengths), sep="\t", dtype={"taxid": str})

    rpm_col = getattr(snakemake.params, "rpm_col", "rpm")
    rpkm_col = getattr(snakemake.params, "rpkm_col", "rpkm")

    out = add_rpkm(df, gl, rpm_col=rpm_col, rpkm_col=rpkm_col)
    out.to_csv(str(snakemake.output[0]), sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
