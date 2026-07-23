"""Split a single multi-record reference (and annotation) into per-segment files.

Segmented consensus runs historically required one ``--segmented-reference
SEGMENT=PATH`` per segment. To make the common case easier, a user can instead
pass a single multi-FASTA to ``--reference``; when it holds more than one
record the pipeline treats it as segmented. This module turns such a file into
the exact ``{segment: path}`` mapping the segmented workflows already consume,
so nothing under ``viralunity/scripts/`` has to change.

Segment names are derived from the FASTA headers: the first whitespace/pipe
token (usually the accession) with ``/ \\ | , ~`` and spaces mapped to ``_``.
This mirrors the sanitisation in
``viralunity/scripts/python/generate_consensus_report.py`` so segment keys,
report contig ids, and nanopore-sanitised headers all agree.
"""

import gzip
import logging
import os
import re
from typing import Dict, Iterable, List, Tuple

logger = logging.getLogger(__name__)

# Same character class the nanopore workflow and consensus report use.
_SANITIZE_RE = re.compile(r"[/\\|,~ ]")


def _open_text(path: str):
    """Open a plain or gzipped text file for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def sanitize_segment_name(header: str) -> str:
    """Derive a filesystem-safe segment key from a FASTA header.

    Strips a leading ``>``, takes the first whitespace/pipe-delimited token,
    then replaces ``/ \\ | , ~`` and spaces with ``_``.

    Raises:
        ValueError: if the header yields an empty token.
    """
    token = header.lstrip(">").strip()
    token = re.split(r"[\s|]", token)[0] if token else ""
    key = _SANITIZE_RE.sub("_", token)
    if not key:
        raise ValueError(f"Could not derive a segment name from header: {header!r}")
    return key


def parse_fasta(path: str) -> List[Tuple[str, List[str]]]:
    """Parse a (optionally gzipped) FASTA into ``(header, body_lines)`` records.

    ``header`` is the header line without the leading ``>`` and without the
    trailing newline. ``body_lines`` are the raw sequence lines (newlines kept)
    so each record can be re-emitted verbatim.
    """
    records: List[Tuple[str, List[str]]] = []
    current_header = None
    current_body: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, current_body))
                current_header = line[1:].strip()
                current_body = []
            elif current_header is not None:
                current_body.append(line)
    if current_header is not None:
        records.append((current_header, current_body))
    return records


def count_records(path: str) -> int:
    """Count the number of records (``>`` headers) in a FASTA file."""
    count = 0
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def _dedup_key(key: str, seen: Dict[str, int]) -> str:
    """Return ``key`` (or ``key_2``, ``key_3`` ...) so every key is unique."""
    if key not in seen:
        seen[key] = 1
        return key
    seen[key] += 1
    new_key = f"{key}_{seen[key]}"
    logger.warning(
        "Duplicate segment name %r derived from reference headers; using %r instead.",
        key,
        new_key,
    )
    return new_key


def split_multifasta(path: str, out_dir: str) -> Dict[str, str]:
    """Split a multi-record FASTA into one single-record file per segment.

    Each record's original header is preserved verbatim inside its file; the
    segment key (dict key, and later the wildcard/directory name) is the
    sanitised header token. Returns ``{segment_key: absolute_path}`` in the
    original record order.

    Raises:
        ValueError: if ``path`` contains no records.
    """
    records = parse_fasta(path)
    if not records:
        raise ValueError(f"Reference FASTA contains no sequences: {path}")

    os.makedirs(out_dir, exist_ok=True)
    mapping: Dict[str, str] = {}
    seen: Dict[str, int] = {}
    for header, body in records:
        key = _dedup_key(sanitize_segment_name(header), seen)
        seg_path = os.path.abspath(os.path.join(out_dir, f"{key}.fasta"))
        text = ">" + header + "\n" + "".join(body)
        if not text.endswith("\n"):
            text += "\n"
        with open(seg_path, "w") as fh:
            fh.write(text)
        mapping[key] = seg_path
    return mapping


def split_annotation_by_segment(
    path: str, out_dir: str, segment_keys: Iterable[str]
) -> Dict[str, str]:
    """Split a single multi-record GFF3/BED by seqid into per-segment files.

    Column-1 seqids are sanitised the same way as segment keys and matched
    against ``segment_keys``. Leading comment/directive lines (``#``/``##``,
    ``track``/``browser``) are preserved at the top of every emitted file.
    Returns ``{segment_key: absolute_path}`` for segments that received at
    least one feature.

    Raises:
        ValueError: if no feature line matches any segment.
    """
    key_by_sanitized = {sanitize_segment_name(k): k for k in segment_keys}

    header_lines: List[str] = []
    grouped: Dict[str, List[str]] = {}
    unmatched: set = set()
    seen_feature = False

    with _open_text(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if line.startswith(("#", "track", "browser")):
                if not seen_feature:
                    header_lines.append(line)
                continue
            seen_feature = True
            seqid = line.split("\t", 1)[0].strip()
            sanitized = sanitize_segment_name(seqid)
            matched = key_by_sanitized.get(sanitized)
            if matched is None:
                unmatched.add(seqid)
                continue
            grouped.setdefault(matched, []).append(line)

    if unmatched:
        logger.warning(
            "Gene-annotation seqids without a matching reference segment (skipped): %s",
            ", ".join(sorted(unmatched)),
        )
    if not grouped:
        raise ValueError(
            f"No annotation features in {path} matched any reference segment "
            f"({', '.join(key_by_sanitized.values())})."
        )

    ext = os.path.splitext(path)[1] or ".gff3"
    os.makedirs(out_dir, exist_ok=True)
    mapping: Dict[str, str] = {}
    for key, lines in grouped.items():
        ann_path = os.path.abspath(os.path.join(out_dir, f"{key}{ext}"))
        with open(ann_path, "w") as fh:
            fh.write("".join(header_lines) + "".join(lines))
        mapping[key] = ann_path
    return mapping
