"""Validation functions for ViralUnity pipeline arguments and data."""

import csv
import logging
import os
import re
from typing import Any, Dict, List, Optional, cast

from viralunity.constants import DataType
from viralunity.exceptions import (
    AdaptersNotFoundError,
    DiamondDatabaseNotFoundError,
    GeneAnnotationNotFoundError,
    Kraken2DatabaseNotFoundError,
    KronaDatabaseNotFoundError,
    PrimerSchemeNotFoundError,
    ReferenceNotFoundError,
    SampleConfigurationNotFoundError,
    SampleSheetError,
    TaxdumpNotFoundError,
    ValidationError,
    ViralUnityFileNotFoundError,
)
from viralunity.reference_splitter import (
    count_records,
    split_annotation_by_segment,
    split_multifasta,
)

logger = logging.getLogger(__name__)


def validate_file_exists(file_path: str, description: str = "File") -> None:
    """Validate that a file exists.

    Args:
        file_path: Path to the file
        description: Description of the file for error messages

    Raises:
        ViralUnityFileNotFoundError: If the file does not exist
    """
    if not os.path.isfile(file_path):
        raise ViralUnityFileNotFoundError(f"{description} does not exist: {file_path}")


def validate_directory_exists(directory_path: str, description: str = "Directory") -> None:
    """Validate that a directory exists.

    Args:
        directory_path: Path to the directory
        description: Description of the directory for error messages

    Raises:
        ViralUnityFileNotFoundError: If the directory does not exist
    """
    if not os.path.isdir(directory_path):
        raise ViralUnityFileNotFoundError(f"{description} does not exist: {directory_path}")


# ---------------------------------------------------------------------------
# Untrusted-input sanitization
# ---------------------------------------------------------------------------
#
# Sample ids and run names become shell tokens and filesystem path components
# inside the Snakemake rules (``sample-<id>``, ``<output>/<run_name>/...``).
# When ViralUnity is embedded in a service that ingests uploaded data, these
# must be constrained so they cannot inject path traversal or shell
# metacharacters. The helpers here are opt-in checks the CLI/service layer can
# call; they intentionally allow only a conservative identifier charset.

_SAFE_IDENTIFIER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def sanitize_identifier(value: Any, field: str = "identifier") -> str:
    """Validate that *value* is a safe filesystem/shell identifier.

    Allows letters, digits, ``.``, ``_`` and ``-`` (and must start with a
    letter/digit). Rejects empty values, path separators, ``..``, NUL, and any
    other metacharacter.

    Returns:
        The stripped value.

    Raises:
        ValidationError: If the value is not a safe identifier.
    """
    if value is None or not str(value).strip():
        raise ValidationError(f"{field} must be a non-empty string.")
    stripped = str(value).strip()
    if stripped in (".", "..") or not _SAFE_IDENTIFIER_RE.match(stripped):
        raise ValidationError(
            f"{field} may only contain letters, digits, '.', '_', '-' and must not "
            f"contain path separators or spaces: {value!r}"
        )
    return stripped


def ensure_within_base(path: str, base: str) -> str:
    """Resolve *path* against *base* and ensure it does not escape *base*.

    Absolute paths are resolved as-is; relative paths are joined to *base*.

    Returns:
        The absolute, normalized target path.

    Raises:
        ValidationError: If the resolved path is outside *base*.
    """
    base_abs = os.path.abspath(base)
    target = os.path.abspath(path if os.path.isabs(path) else os.path.join(base_abs, path))
    if os.path.commonpath([base_abs, target]) != base_abs:
        raise ValidationError(f"Path escapes base directory {base!r}: {path!r}")
    return target


# Objective validity bounds for numeric parameters. These are correctness
# constraints (a thread count must be >= 1; an allele frequency is a fraction in
# [0, 1]), NOT scientific tuning choices — analysis knobs such as
# z_score_threshold / log10_ratio_threshold are intentionally left unbounded.
_NUMERIC_BOUNDS = {
    "threads": (1, None),
    "threads_total": (1, None),
    "minimum_coverage": (1, None),
    "minimum_depth": (1, None),
    "minimum_length": (0, None),
    "minimum_read_length": (0, None),
    "af_threshold": (0.0, 1.0),
    "af_isnv_threshold": (0.0, 1.0),
    "bleed_fraction": (0.0, 1.0),
}
# Parameters that must be strictly positive.
_POSITIVE_NUMERIC_KEYS = ("evalue", "enrichment_pseudocount")


def validate_numeric_parameters(args: Dict[str, Any]) -> None:
    """Range-check numeric parameters that have objective validity bounds.

    Only keys present in ``args`` are checked, so the same function is safe for
    both the consensus and metagenomics pipelines. Raises on the first
    out-of-range value.

    Raises:
        ValidationError: If a present numeric parameter is non-numeric or falls
            outside its valid range.
    """

    def _as_number(key: str) -> Optional[float]:
        if key not in args or args[key] is None:
            return None
        val = args[key]
        # bool is a subclass of int; reject it explicitly.
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValidationError(f"{key} must be a number, got {val!r}.")
        return val

    for key, (lo, hi) in _NUMERIC_BOUNDS.items():
        val = _as_number(key)
        if val is None:
            continue
        if lo is not None and val < lo:
            raise ValidationError(f"{key} must be >= {lo}, got {val}.")
        if hi is not None and val > hi:
            raise ValidationError(f"{key} must be <= {hi}, got {val}.")

    for key in _POSITIVE_NUMERIC_KEYS:
        val = _as_number(key)
        if val is None:
            continue
        if val <= 0:
            raise ValidationError(f"{key} must be > 0, got {val}.")


def validate_sample_sheet(sample_sheet_path: str, data_type: str) -> Dict[str, List[str]]:
    """Validate and parse sample sheet file.

    Args:
        sample_sheet_path: Path to the sample sheet CSV file
        data_type: Type of sequencing data (illumina or nanopore)

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        SampleSheetError: If the sample sheet is invalid
        ViralUnityFileNotFoundError: If the sample sheet file does not exist
    """
    validate_file_exists(sample_sheet_path, "Sample sheet file")

    expected_columns = 3 if data_type == DataType.ILLUMINA else 2

    try:
        with open(sample_sheet_path, newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
    except (OSError, UnicodeDecodeError) as e:
        raise SampleSheetError(f"Failed to read sample sheet: {e}") from e

    samples: Dict[str, List[str]] = {}

    for line_number, row in enumerate(rows, start=1):
        # Skip blank lines (rows with no cells, or only empty cells).
        if not row or all(not cell.strip() for cell in row):
            continue

        # Tolerate trailing empty fields (e.g. a stray trailing comma) but
        # reject rows whose real column count does not match the data type,
        # instead of silently NaN-padding them the way pandas did.
        trimmed = list(row)
        while trimmed and not trimmed[-1].strip():
            trimmed.pop()

        if len(trimmed) != expected_columns:
            raise SampleSheetError(
                f"{str(data_type).capitalize()} sample sheet requires exactly "
                f"{expected_columns} columns; row {line_number} has {len(trimmed)}: {row!r}"
            )

        sample_name = trimmed[0].strip()
        if not sample_name:
            raise SampleSheetError(f"Empty sample name on row {line_number}.")

        # Sample ids become shell tokens and path components in the workflow
        # (``sample-<id>``); reject anything that could inject a path or
        # metacharacter when the sheet comes from untrusted input.
        if sample_name in (".", "..") or not _SAFE_IDENTIFIER_RE.match(sample_name):
            raise SampleSheetError(
                f"Unsafe sample id '{sample_name}' on row {line_number}: use only letters, "
                f"digits, '.', '_', '-' (no path separators or spaces)."
            )

        if sample_name in samples:
            raise SampleSheetError(
                f"Duplicate sample id '{sample_name}' on row {line_number}. "
                f"Sample ids must be unique so no sample is silently dropped."
            )

        file_paths = [cell.strip() for cell in trimmed[1:]]
        for offset, file_path in enumerate(file_paths, start=2):
            if not file_path:
                raise SampleSheetError(
                    f"Missing file path in column {offset} for sample "
                    f"'{sample_name}' on row {line_number}."
                )

        if data_type == DataType.ILLUMINA:
            validate_file_exists(file_paths[0], f"R1 file for sample {sample_name}")
            validate_file_exists(file_paths[1], f"R2 file for sample {sample_name}")
        else:  # nanopore
            validate_file_exists(file_paths[0], f"File for sample {sample_name}")

        samples[sample_name] = file_paths

    if not samples:
        raise SampleSheetError("No valid samples found in sample sheet")

    return samples


def validate_illumina_requirements(args: Dict[str, Any]) -> None:
    """Validate Illumina-specific requirements (fastp: adapters optional)."""
    if args.get("data_type") != DataType.ILLUMINA:
        return

    adapters = args.get("adapters")
    if adapters and str(adapters).strip() and str(adapters).strip() != "NA":
        try:
            validate_file_exists(adapters, "Illumina adapter sequences file")
        except ViralUnityFileNotFoundError as e:
            raise AdaptersNotFoundError(
                f"Illumina adapter sequences file does not exist: {e}"
            ) from e


def validate_nanopore_requirements(args: Dict[str, Any]) -> None:
    """Validate Nanopore-specific requirements (e.g. polishing options)."""
    if args.get("data_type") != DataType.NANOPORE:
        return
    # Optional: validate medaka_model if provided (Medaka accepts known model names)
    # For now no strict validation; extend as needed.


def _segment_input_dir(args: Dict[str, Any]) -> str:
    """Directory where auto-split per-segment inputs are written.

    Mirrors ``ConfigGenerator.add_output``'s ``output/run_name/`` layout so the
    generated per-segment references live inside the run's directory.
    """
    return os.path.join(args["output"], args.get("run_name") or "", "input_references")


def _maybe_split_multifasta_reference(args: Dict[str, Any]) -> bool:
    """Auto-split a single multi-record ``--reference`` into per-segment files.

    A ``--reference`` FASTA with more than one record is treated as a segmented
    virus: it is split into one file per record and ``args["reference"]`` is
    replaced with the ``{segment: path}`` mapping the segmented workflows
    consume. ``--single-reference`` forces the historical single-reference
    behaviour (all records aligned together in one pass).

    Returns:
        True if the reference was split (segmented mode engaged), else False.
    """
    reference = args.get("reference")
    if args.get("single_reference"):
        return False
    if not isinstance(reference, str):
        return False
    if count_records(reference) <= 1:
        return False

    segment_map = split_multifasta(reference, _segment_input_dir(args))
    logger.info(
        "Reference '%s' contains %d records; running in segmented mode with " "segments: %s",
        reference,
        len(segment_map),
        ", ".join(segment_map),
    )
    args["reference"] = segment_map
    return True


def _maybe_split_multifasta_annotation(args: Dict[str, Any]) -> None:
    """Split a single combined gene annotation to match an auto-split reference.

    Only called when the reference was auto-split from a multi-FASTA, so the
    segment keys are sanitised header tokens that annotation seqids can match.
    A ``--segmented-gene-annotation`` (already per-segment) is left untouched.
    """
    gene_annotation = args.get("gene_annotation")
    if not gene_annotation or isinstance(gene_annotation, dict):
        return
    if args.get("segmented_gene_annotation"):
        return

    try:
        validate_file_exists(cast(str, gene_annotation), "Gene annotation file")
    except ViralUnityFileNotFoundError as e:
        raise GeneAnnotationNotFoundError(f"Gene annotation file does not exist: {e}") from e

    out_dir = os.path.join(_segment_input_dir(args), "annotation")
    try:
        annotation_map = split_annotation_by_segment(
            cast(str, gene_annotation), out_dir, list(args["reference"].keys())
        )
    except ValueError as e:
        raise ValidationError(str(e)) from e
    logger.info(
        "Gene annotation '%s' split into %d per-segment files.",
        gene_annotation,
        len(annotation_map),
    )
    args["gene_annotation"] = annotation_map


def validate_consensus_requirements(args: Dict[str, Any]) -> None:
    """Validate consensus pipeline requirements.

    Exactly one of --reference or --segmented-reference must be provided.
    When --segmented-reference is used, the values are parsed from
    SEGMENT=PATH format and stored as a dict in args["reference"].

    Args:
        args: Dictionary of pipeline arguments

    Raises:
        ValidationError: If consensus requirements are not met
    """
    reference = args.get("reference")
    segmented_reference = args.get("segmented_reference")

    if reference and segmented_reference:
        raise ValidationError(
            "--reference and --segmented-reference are mutually exclusive. "
            "Please provide only one."
        )

    if not reference and not segmented_reference:
        raise ValidationError(
            "A reference is required. Provide --reference for a single reference "
            "or --segmented-reference for segmented viruses."
        )

    if segmented_reference:
        if isinstance(segmented_reference, dict):
            parsed_segments = segmented_reference
        else:
            parsed_segments = {}
            for entry in segmented_reference:
                if "=" not in entry:
                    raise ValidationError(
                        f"Invalid segmented reference format: '{entry}'. "
                        f"Expected SEGMENT=PATH (e.g. S=/path/to/S.fasta)"
                    )
                segment_name, segment_path = entry.split("=", 1)
                segment_name = segment_name.strip()
                segment_path = segment_path.strip()
                if not segment_name or not segment_path:
                    raise ValidationError(
                        f"Invalid segmented reference format: '{entry}'. "
                        f"Both segment name and path are required."
                    )
                parsed_segments[segment_name] = segment_path

        args["reference"] = parsed_segments
        args["segmented_reference"] = None
        reference = parsed_segments

    autosplit = False
    if isinstance(reference, dict):
        for segment_name, segment_path in reference.items():
            try:
                validate_file_exists(segment_path, f"Reference file for segment '{segment_name}'")
            except ViralUnityFileNotFoundError as e:
                raise ReferenceNotFoundError(str(e)) from e
    else:
        try:
            validate_file_exists(cast(str, reference), "Reference sequence file")
        except ViralUnityFileNotFoundError as e:
            raise ReferenceNotFoundError(f"Reference sequence file does not exist: {e}") from e
        autosplit = _maybe_split_multifasta_reference(args)

    primer_scheme = args.get("primer_scheme")
    if primer_scheme:
        try:
            validate_file_exists(primer_scheme, "Primer scheme file")
        except ViralUnityFileNotFoundError as e:
            raise PrimerSchemeNotFoundError(f"Primer scheme file does not exist: {e}") from e

    # When the reference was auto-split from a multi-FASTA, its segment keys are
    # the sanitised header tokens, so a single combined gene annotation can be
    # split by seqid to match. (Explicit --segmented-reference keys are user
    # chosen and may not match annotation seqids, so we don't auto-split there.)
    if autosplit:
        _maybe_split_multifasta_annotation(args)

    _validate_gene_annotation(args)


def _validate_gene_annotation(args: Dict[str, Any]) -> None:
    """Validate the optional gene-annotation GFF3 input.

    Both ``--gene-annotation`` (a single path) and ``--segmented-gene-annotation``
    (SEGMENT=PATH pairs, collapsed into ``args["gene_annotation"]`` as a dict)
    are optional and mutually exclusive. Providing neither is legal.
    """
    gene_annotation = args.get("gene_annotation")
    segmented = args.get("segmented_gene_annotation")

    if gene_annotation and segmented:
        raise ValidationError(
            "--gene-annotation and --segmented-gene-annotation are mutually "
            "exclusive. Please provide only one."
        )

    if segmented:
        if isinstance(segmented, dict):
            parsed_segments = segmented
        else:
            parsed_segments = {}
            for entry in segmented:
                if "=" not in entry:
                    raise ValidationError(
                        f"Invalid segmented gene annotation format: '{entry}'. "
                        f"Expected SEGMENT=PATH (e.g. S=/path/to/S.gff3)"
                    )
                segment_name, segment_path = entry.split("=", 1)
                segment_name = segment_name.strip()
                segment_path = segment_path.strip()
                if not segment_name or not segment_path:
                    raise ValidationError(
                        f"Invalid segmented gene annotation format: '{entry}'. "
                        f"Both segment name and path are required."
                    )
                parsed_segments[segment_name] = segment_path

        args["gene_annotation"] = parsed_segments
        args["segmented_gene_annotation"] = None
        gene_annotation = parsed_segments

    if isinstance(gene_annotation, dict):
        for segment_name, segment_path in gene_annotation.items():
            try:
                validate_file_exists(
                    segment_path, f"Gene annotation file for segment '{segment_name}'"
                )
            except ViralUnityFileNotFoundError as e:
                raise GeneAnnotationNotFoundError(str(e)) from e
    elif gene_annotation:
        try:
            validate_file_exists(cast(str, gene_annotation), "Gene annotation file")
        except ViralUnityFileNotFoundError as e:
            raise GeneAnnotationNotFoundError(f"Gene annotation file does not exist: {e}") from e


def validate_metagenomics_requirements(args: Dict[str, Any]) -> None:
    """Validate metagenomics pipeline requirements.

    Kraken2 and Diamond are optional; validate only the resources for the tools
    the user has enabled (run_kraken2_reads, run_kraken2_contigs, run_diamond_reads,
    run_diamond_contigs).
    """
    run_denovo = args.get("run_denovo_assembly", True)
    run_k2_reads = args.get("run_kraken2_reads", True)
    run_k2_contigs = args.get("run_kraken2_contigs", True)
    run_diamond_reads = args.get("run_diamond_reads", False)
    run_diamond_contigs = args.get("run_diamond_contigs", False)

    if not run_denovo:
        if run_k2_contigs:
            raise ValidationError(
                "Cannot run kraken2 on contigs (--run-kraken2-contigs) when denovo assembly is disabled (--no-denovo-assembly)."
            )
        if run_diamond_contigs:
            raise ValidationError(
                "Cannot run diamond on contigs (--run-diamond-contigs) when denovo assembly is disabled (--no-denovo-assembly)."
            )

    any_kraken2 = run_k2_reads or run_k2_contigs
    any_diamond = run_diamond_reads or run_diamond_contigs
    any_classification = any_kraken2 or any_diamond

    # Krona database: required whenever any classification is run (Krona plots for both tools)
    if any_classification:
        krona_db = args.get("krona_database")
        if not krona_db or krona_db == "NA":
            raise KronaDatabaseNotFoundError(
                "Krona database directory is required when running any classification "
                "(Kraken2 and/or Diamond). Set --krona-database."
            )
        try:
            validate_directory_exists(krona_db, "Krona database directory")
        except ViralUnityFileNotFoundError as e:
            raise KronaDatabaseNotFoundError(f"Krona database directory does not exist: {e}") from e

    # Kraken2 database: required only when Kraken2 is enabled
    if any_kraken2:
        kraken2_db = args.get("kraken2_database")
        if not kraken2_db or kraken2_db == "NA":
            raise Kraken2DatabaseNotFoundError(
                "Kraken2 database directory is required when running Kraken2. "
                "Set --kraken2-database or disable Kraken2 with --no-kraken2-reads / --no-kraken2-contigs."
            )
        try:
            validate_directory_exists(kraken2_db, "Kraken2 database directory")
        except ViralUnityFileNotFoundError as e:
            raise Kraken2DatabaseNotFoundError(
                f"Kraken2 database directory does not exist: {e}"
            ) from e

    # Taxdump: required for taxonomic summaries whenever any classification is run
    if any_classification:
        taxdump = args.get("taxdump", "NA")
        if not taxdump or taxdump == "NA":
            raise TaxdumpNotFoundError(
                "taxdump directory (NCBI nodes.dmp, names.dmp) is required for taxonomic summaries "
                "when running any classification. Set --taxdump."
            )
        try:
            validate_directory_exists(taxdump, "Taxdump directory")
        except ViralUnityFileNotFoundError as e:
            raise TaxdumpNotFoundError(f"Taxdump directory does not exist: {e}") from e
        nodes = os.path.join(taxdump, "nodes.dmp")
        names = os.path.join(taxdump, "names.dmp")
        if not os.path.isfile(nodes) or not os.path.isfile(names):
            raise TaxdumpNotFoundError(
                f"Taxdump directory must contain nodes.dmp and names.dmp: {taxdump}"
            )

    # Diamond: require database and assembly summary only when Diamond is enabled
    if any_diamond:
        diamond_db = args.get("diamond_database", "NA")
        if not diamond_db or diamond_db == "NA":
            raise DiamondDatabaseNotFoundError(
                "diamond_database is required when running Diamond. "
                "Set --diamond-database or do not use --run-diamond-reads / --run-diamond-contigs."
            )
        taxids = args.get("taxids", "NA")
        if not taxids or taxids == "NA":
            raise DiamondDatabaseNotFoundError(
                "taxids mapping file is required when running Diamond. "
                "Set --taxids or do not use --run-diamond-reads / --run-diamond-contigs."
            )
        if not os.path.isfile(taxids):
            raise DiamondDatabaseNotFoundError(f"Taxid mapping file not found: {taxids}")

    # Deacon index: when provided for host depletion, must exist
    deacon_idx = args.get("deacon_index", "NA")
    if deacon_idx and str(deacon_idx).strip() not in ("", "NA"):
        try:
            validate_file_exists(deacon_idx, "Deacon index file")
        except ViralUnityFileNotFoundError as e:
            raise ViralUnityFileNotFoundError(f"Deacon index file does not exist: {e}") from e

    # RPKM normalisation (compute_rpkm) is enabled whenever a viral genomes FASTA
    # is provided, independent of reference assembly. The genome-length table it
    # needs requires both the FASTA and the genome2taxid (--viral-taxids) mapping,
    # so validate them up-front rather than failing deep inside Snakemake.
    viral_genomes = args.get("viral_genomes", "NA")
    if viral_genomes and str(viral_genomes).strip() not in ("", "NA"):
        if not os.path.isfile(viral_genomes):
            raise ViralUnityFileNotFoundError(f"Viral genomes file does not exist: {viral_genomes}")
        viral_taxids = args.get("viral_taxids", "NA")
        if not viral_taxids or str(viral_taxids).strip() in ("", "NA"):
            raise ValidationError(
                "--viral-genomes enables RPKM normalisation, which requires a "
                "genome2taxid mapping. Provide --viral-taxids, or omit "
                "--viral-genomes to skip RPKM."
            )
        if not os.path.isfile(viral_taxids):
            raise ViralUnityFileNotFoundError(f"Viral taxids file does not exist: {viral_taxids}")

    # ICTV vertebrate-virus host filter: allowlist taxid file must exist when on.
    if args.get("run_ictv_host_filter", False):
        ictv_file = args.get("ictv_vertebrate_taxids_file", "NA")
        if not ictv_file or str(ictv_file).strip() in ("", "NA"):
            raise ValidationError(
                "--run-ictv-host-filter requires --ictv-vertebrate-taxids-file "
                "(the vertebrate-virus taxid allowlist built by "
                "build_ictv_vertebrate_taxids.py)."
            )
        if not os.path.isfile(ictv_file):
            raise ViralUnityFileNotFoundError(
                f"ICTV vertebrate taxids file does not exist: {ictv_file}"
            )

    # NR validation: contig tracks only; needs denovo + diamond_contigs and an
    # nr database that resolves as a BLAST+ db or a native .dmnd.
    if args.get("run_nr_validation", False):
        if not args.get("run_denovo_assembly", False) or not args.get("run_diamond_contigs", False):
            raise ValidationError(
                "--run-nr-validation requires --run-denovo-assembly and "
                "--run-diamond-contigs (the viral-contig set is diamond-defined)."
            )
        nr_db = args.get("nr_diamond_database", "NA")
        if not nr_db or str(nr_db).strip() in ("", "NA"):
            raise ValidationError(
                "--run-nr-validation requires --nr-diamond-database (a BLAST+ nr "
                "database or a .dmnd)."
            )
        is_dmnd = str(nr_db).endswith(".dmnd") and os.path.isfile(nr_db)
        is_blastdb = any(
            os.path.isfile(f"{nr_db}{suffix}") for suffix in (".pal", ".pin", ".000.pin", ".dmnd")
        )
        if not (is_dmnd or is_blastdb):
            raise ViralUnityFileNotFoundError(
                f"nr database not found (expected a .dmnd or BLAST+ nr db): {nr_db}"
            )

    validate_reference_assembly_requirements(args)


def validate_reference_assembly_requirements(args: Dict[str, Any]) -> None:
    """Validate cross-dependencies for Reference Assembly in metagenomics workflows."""
    if not args.get("run_reference_assembly"):
        return

    method = args.get("method")
    source = args.get("source")
    strategy = args.get("reference_selection_strategy")

    # If --similarity is used, it must be on contigs, and the respective classification tool must run on contigs
    if strategy == "similarity":
        run_denovo = args.get("run_denovo_assembly", False)
        if not run_denovo:
            raise ValidationError(
                "Reference selection strategy 'similarity' requires --run-denovo-assembly."
            )
        if source == "reads":
            raise ValidationError(
                "Reference selection strategy 'similarity' can only be used if --source includes 'contigs'."
            )

        # Check if the chosen method actually runs on contigs
        if method == "kraken2" and not args.get("run_kraken2_contigs", True):
            raise ValidationError(
                "Strategy 'similarity' with method 'kraken2' requires --run-kraken2-contigs."
            )
        if method == "diamond" and not args.get("run_diamond_contigs", False):
            raise ValidationError(
                "Strategy 'similarity' with method 'diamond' requires --run-diamond-contigs."
            )
        if method == "both" and not (
            args.get("run_kraken2_contigs", True) or args.get("run_diamond_contigs", False)
        ):
            raise ValidationError(
                "Strategy 'similarity' with method 'both' requires at least one mapping tool on contigs."
            )

    if method in ["kraken2", "both"]:
        if source in ["reads", "both"] and not args.get("run_kraken2_reads", True):
            raise ValidationError(
                "Method includes 'kraken2' on 'reads' but --no-kraken2-reads was passed."
            )
        if source in ["contigs", "both"] and not args.get("run_kraken2_contigs", True):
            raise ValidationError(
                "Method includes 'kraken2' on 'contigs' but --no-kraken2-contigs was passed."
            )

    if method in ["diamond", "both"]:
        if source in ["reads", "both"] and not args.get("run_diamond_reads", False):
            raise ValidationError(
                "Method includes 'diamond' on 'reads' but --run-diamond-reads is not enabled."
            )
        if source in ["contigs", "both"] and not args.get("run_diamond_contigs", False):
            raise ValidationError(
                "Method includes 'diamond' on 'contigs' but --run-diamond-contigs is not enabled."
            )

    viral_genomes = args.get("viral_genomes")
    if viral_genomes and not os.path.isfile(viral_genomes):
        raise ViralUnityFileNotFoundError(f"Viral genomes file does not exist: {viral_genomes}")

    if strategy == "taxid":
        viral_taxids = args.get("viral_taxids")
        if viral_taxids and not os.path.isfile(viral_taxids):
            raise ViralUnityFileNotFoundError(f"Viral taxids file does not exist: {viral_taxids}")


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
#
# CLI arguments that point to a filesystem location. The lists are kept here
# (next to the validators) so they can be reused both at validation time and
# inside ``viralunity_meta.main`` / ``viralunity_consensus.main`` to make the
# paths absolute before anything else sees them. Resolving paths up-front
# means a user who runs ``viralunity meta nanopore ... --host-reference
# databases/host/host.fasta --config-file scratch/run.yml`` does NOT silently
# get the host reference looked up under ``scratch/databases/host/...`` just
# because the generated config happens to live in ``scratch/``.

META_PATH_ARG_KEYS = (
    "sample_sheet",
    "config_file",
    "output",
    "kraken2_database",
    "krona_database",
    "taxdump",
    "host_reference",
    "deacon_index",
    "taxids",
    "diamond_database",
    "ictv_vertebrate_taxids_file",
    "nr_diamond_database",
    "viral_genomes",
    "viral_taxids",
    "adapters",
)

CONSENSUS_PATH_ARG_KEYS = (
    "sample_sheet",
    "config_file",
    "output",
    "reference",
    "primer_scheme",
    "gene_annotation",
    "adapters",
)


def _is_path_sentinel(value: Any) -> bool:
    """Return True if a path-typed argument value should be left untouched.

    Sentinels are: ``None``, non-string scalars (ints, bools), the empty
    string, and the literal placeholder ``"NA"`` (after stripping
    whitespace). This matches the conventions used by the existing
    validators in this module.
    """
    if value is None:
        return True
    if not isinstance(value, str):
        return True
    stripped = value.strip()
    if not stripped:
        return True
    if stripped == "NA":
        return True
    return False


def resolve_path_args(
    args: Dict[str, Any],
    keys,
    base_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Rewrite path-typed argument values to absolute paths in place.

    Each value listed in ``keys`` that is a non-sentinel relative path string
    is replaced with ``os.path.abspath(os.path.join(base_dir, value))``.
    Absolute paths, sentinels (``None``, ``""``, ``"NA"``), and non-string
    values are left unchanged. Missing keys are ignored.

    The ``reference`` argument of the consensus pipeline can be a dict
    (segmented reference, ``segment -> path``); each value of the dict is
    resolved while preserving the keys.

    Args:
        args: Mutable dict of CLI arguments.
        keys: Iterable of argument keys whose values are filesystem paths.
        base_dir: Base directory to resolve relative paths against.
            Defaults to the current working directory at call time.

    Returns:
        The same ``args`` dict (modified in place). Returned for chaining.
    """
    base_dir = base_dir or os.getcwd()

    for key in keys:
        if key not in args:
            continue
        value = args[key]

        if isinstance(value, dict):
            args[key] = {
                seg: (
                    v
                    if _is_path_sentinel(v) or os.path.isabs(v)
                    else os.path.abspath(os.path.join(base_dir, v))
                )
                for seg, v in value.items()
            }
            continue

        if _is_path_sentinel(value):
            continue
        if os.path.isabs(value):
            continue

        args[key] = os.path.abspath(os.path.join(base_dir, value))

    return args


def get_samples_from_args(args: Dict[str, Any]) -> Dict[str, List[str]]:
    """Extract and validate samples from arguments.

    Args:
        args: Dictionary of pipeline arguments

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If samples cannot be determined from arguments
    """
    sample_sheet = args.get("sample_sheet")
    samples = args.get("samples")
    data_type = args.get("data_type")

    if sample_sheet:
        # A path was provided: it must exist. Do not silently fall back to
        # `samples`, which would report a misleading "nothing provided" error
        # when the real problem is a mistyped or missing sample-sheet path.
        validate_file_exists(sample_sheet, "Sample sheet file")
        return validate_sample_sheet(sample_sheet, cast(str, data_type))
    if samples:
        return samples
    raise SampleConfigurationNotFoundError(
        "Either 'sample_sheet' or 'samples' must be provided. " f"Sample sheet path: {sample_sheet}"
    )
