#!/usr/bin/env python3
"""Stream the GSE239676 files and update the predownload structure audit.

This resource-preparation audit never modifies the downloaded source files and
does not inspect expression values beyond the MatrixMarket header/dimensions.
The three GEO files with a ``.gz`` suffix are opened as gzip first and fall
back to plain UTF-8 text when no gzip stream is present.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import re
from pathlib import Path
from typing import Dict, Iterator, List, TextIO, Tuple

from f0_utils import normalize_sha256, read_tsv, rel, sha256_file, write_tsv


AUDIT_FIELDS = [
    "dataset_id",
    "file_name",
    "relative_path",
    "source_url",
    "file_role",
    "download_status",
    "file_size_bytes",
    "sha256",
    "read_status",
    "n_rows",
    "n_columns",
    "n_sample_or_cell_columns",
    "gene_or_feature_id_type",
    "metadata_or_sample_id_examples",
    "value_type_or_key_fields",
    "first_feature_examples",
    "major_limitations",
    "method_implication",
    "audit_decision",
]

FILE_ROLES = {
    "GSE239676_barcodes.tsv.gz": "geo_cell_barcodes_plain_text_with_gz_suffix",
    "GSE239676_features.tsv.gz": "geo_gene_features_plain_text_with_gz_suffix",
    "GSE239676_meta.tsv.gz": "geo_cell_metadata_plain_text_with_gz_suffix",
    "GSE239676_count_matrix.mtx.gz": "geo_sparse_integer_count_matrix_market",
    "GSE239676_series_matrix.txt.gz": "geo_series_matrix_metadata",
}


@contextlib.contextmanager
def open_text_auto(path: Path) -> Iterator[Tuple[TextIO, str]]:
    """Open gzip text first, falling back to plain text and report the mode."""
    gzip_handle: TextIO | None = None
    try:
        gzip_handle = gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
        gzip_handle.read(1)
        gzip_handle.seek(0)
    except (gzip.BadGzipFile, OSError, EOFError):
        if gzip_handle is not None:
            gzip_handle.close()
    else:
        try:
            yield gzip_handle, "gzip_read_ok"
        finally:
            gzip_handle.close()
        return
    plain_handle = path.open("rt", encoding="utf-8", errors="replace", newline="")
    try:
        yield plain_handle, "plain_text_read_ok_gz_suffix"
    finally:
        plain_handle.close()


def pipe_examples(values: List[str], limit: int = 8) -> str:
    return "|".join(values[:limit])


def audit_one_column(path: Path) -> Tuple[int, List[str], str]:
    count = 0
    examples: List[str] = []
    with open_text_auto(path) as (handle, read_status):
        for raw in handle:
            value = raw.rstrip("\r\n")
            if not value:
                continue
            count += 1
            if len(examples) < 8:
                examples.append(value)
    return count, examples, read_status


def audit_metadata(path: Path) -> Dict[str, object]:
    with open_text_auto(path) as (handle, read_status):
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        unique: Dict[str, set[str]] = {field: set() for field in fields}
        examples: List[str] = []
        n_rows = 0
        for row in reader:
            n_rows += 1
            if len(examples) < 5:
                examples.append("|".join(str(row.get(field, "")) for field in fields))
            for field in fields:
                unique[field].add(str(row.get(field, "")))
    return {
        "n_rows": n_rows,
        "fields": fields,
        "unique": unique,
        "examples": examples,
        "read_status": read_status,
    }


def audit_matrix_market(path: Path) -> Dict[str, object]:
    with open_text_auto(path) as (handle, read_status):
        banner = handle.readline().strip()
        comments: List[str] = []
        dimension_line = ""
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("%"):
                if len(comments) < 3:
                    comments.append(line)
                continue
            dimension_line = line
            break
    parts = dimension_line.split()
    if len(parts) != 3 or not all(part.isdigit() for part in parts):
        raise RuntimeError(f"Invalid MatrixMarket dimension line in {path}: {dimension_line!r}")
    n_features, n_cells, n_nonzero = (int(part) for part in parts)
    return {
        "banner": banner,
        "comments": comments,
        "n_features": n_features,
        "n_cells": n_cells,
        "n_nonzero": n_nonzero,
        "read_status": read_status,
    }


def audit_series_matrix(path: Path) -> Dict[str, object]:
    sample_ids: List[str] = []
    sample_titles: List[str] = []
    sample_sources: List[str] = []
    summaries: List[str] = []
    designs: List[str] = []
    with open_text_auto(path) as (handle, read_status):
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("!Series_sample_id"):
                sample_ids.extend(re.findall(r"GSM\d+", line))
            elif line.startswith("!Sample_geo_accession"):
                sample_ids.extend(re.findall(r"GSM\d+", line))
            elif line.startswith("!Sample_title"):
                sample_titles.extend(
                    value.strip().strip('"')
                    for value in next(csv.reader([line], delimiter="\t"))[1:]
                )
            elif line.startswith("!Sample_source_name_ch1"):
                sample_sources.extend(
                    value.strip().strip('"')
                    for value in next(csv.reader([line], delimiter="\t"))[1:]
                )
            elif line.startswith("!Series_summary"):
                summaries.append(line.split("\t", 1)[-1].strip('"'))
            elif line.startswith("!Series_overall_design"):
                designs.append(line.split("\t", 1)[-1].strip('"'))
    return {
        "sample_ids": sorted(set(sample_ids)),
        "sample_titles": sample_titles,
        "sample_sources": sample_sources,
        "summaries": summaries,
        "designs": designs,
        "read_status": read_status,
    }


def source_urls(root: Path) -> Dict[str, str]:
    urls: Dict[str, str] = {}
    for manifest in [
        root / "reports/download_resources/GSE239676_download_manifest.tsv",
        root / "data/metadata/download_manifest.tsv",
    ]:
        if not manifest.exists():
            continue
        for row in read_tsv(manifest):
            if row.get("dataset_id") == "GSE239676":
                urls[row.get("file_name", "")] = row.get("source_url", "")
    return urls


def build_rows(root: Path) -> List[Dict[str, object]]:
    base = root / "data/public_downloads/GSE239676"
    paths = {name: base / name for name in FILE_ROLES}
    missing = [name for name, path in paths.items() if not path.exists()]
    if missing:
        raise RuntimeError("Missing GSE239676 files: " + ", ".join(missing))

    barcode_count, barcode_examples, barcode_status = audit_one_column(paths["GSE239676_barcodes.tsv.gz"])
    feature_count, feature_examples, feature_status = audit_one_column(paths["GSE239676_features.tsv.gz"])
    metadata = audit_metadata(paths["GSE239676_meta.tsv.gz"])
    matrix = audit_matrix_market(paths["GSE239676_count_matrix.mtx.gz"])
    series = audit_series_matrix(paths["GSE239676_series_matrix.txt.gz"])

    patients = sorted(metadata["unique"].get("Patient", set()))
    samples = sorted(metadata["unique"].get("Sample", set()))
    tissues = sorted(metadata["unique"].get("Tissue", set()))
    has_peritoneal = any(
        title.startswith("PC ") or "Peritoneal" in source
        for title, source in zip(series["sample_titles"], series["sample_sources"])
    )
    has_liver = any(
        title.startswith("LM ") or "Liver" in source
        for title, source in zip(series["sample_titles"], series["sample_sources"])
    )
    structure_ok = (
        matrix["banner"].lower() == "%%matrixmarket matrix coordinate integer general"
        and matrix["n_features"] == feature_count
        and matrix["n_cells"] == barcode_count == metadata["n_rows"]
        and len(patients) == 20
        and has_peritoneal
        and has_liver
    )
    consistency = (
        f"matrix={matrix['n_features']}x{matrix['n_cells']}; features={feature_count}; "
        f"barcodes={barcode_count}; metadata_rows={metadata['n_rows']}; "
        f"patients={len(patients)}; metadata_samples={len(samples)}; tissues={','.join(tissues)}; "
        f"GEO_samples={len(series['sample_ids'])}; peritoneal_in_GEO={has_peritoneal}; liver_in_GEO={has_liver}"
    )
    decision = "downloaded_structure_preaudit_pass" if structure_ok else "pause_structure_mismatch"
    limitation = (
        f"The matrix contains only {feature_count} feature rows, substantially fewer than the 26571 rows in the "
        "GSE183904 main dataset, which indicates a restricted or prefiltered feature space; the exact author "
        "filtering rule remains to be confirmed in the F2.4 cohort audit. MLMOD human_scoring_signature coverage "
        "must be audited before scoring, and fixed UCell maxRank=1500 represents a different fraction of the "
        "available feature space across cohorts, so absolute UCell scores must not be treated as directly "
        "comparable. "
        "External validation cohort only; do not use expression values, labels or outcomes to tune the F2 signature, "
        "UCell parameters, thresholds or candidate ranking. Full F2.4 cohort audit remains required."
    )
    implication = (
        "Matrix dimensions and metadata composition are eligible for later F2.4 external-cohort audit after signature "
        "and scoring rules are frozen. F2.4 must report signature-gene coverage and evaluate fixed-maxRank "
        "within-cohort score behavior; validation should use prespecified within-cohort contrasts rather than "
        "cross-cohort absolute-score equality. This preaudit does not authorize biological analysis."
    )
    urls = source_urls(root)

    common: Dict[str, object] = {
        "dataset_id": "GSE239676",
        "download_status": "complete",
        "major_limitations": limitation,
        "method_implication": implication,
        "audit_decision": decision,
    }

    def base_row(name: str, read_status: str) -> Dict[str, object]:
        path = paths[name]
        return {
            **common,
            "file_name": name,
            "relative_path": rel(path, root),
            "source_url": urls.get(name, ""),
            "file_role": FILE_ROLES[name],
            "file_size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "read_status": read_status,
        }

    rows: List[Dict[str, object]] = []
    rows.append(
        {
            **base_row("GSE239676_barcodes.tsv.gz", barcode_status),
            "n_rows": barcode_count,
            "n_columns": 1,
            "n_sample_or_cell_columns": barcode_count,
            "gene_or_feature_id_type": "not_applicable_cell_barcode",
            "metadata_or_sample_id_examples": pipe_examples(barcode_examples),
            "value_type_or_key_fields": "cell_barcode_one_per_line; " + consistency,
            "first_feature_examples": "",
        }
    )
    rows.append(
        {
            **base_row("GSE239676_features.tsv.gz", feature_status),
            "n_rows": feature_count,
            "n_columns": 1,
            "n_sample_or_cell_columns": "",
            "gene_or_feature_id_type": "human_gene_symbol",
            "metadata_or_sample_id_examples": "",
            "value_type_or_key_fields": "gene_symbol_one_per_line; " + consistency,
            "first_feature_examples": pipe_examples(feature_examples),
        }
    )
    rows.append(
        {
            **base_row("GSE239676_meta.tsv.gz", metadata["read_status"]),
            "n_rows": metadata["n_rows"],
            "n_columns": len(metadata["fields"]),
            "n_sample_or_cell_columns": metadata["n_rows"],
            "gene_or_feature_id_type": "not_applicable_cell_metadata",
            "metadata_or_sample_id_examples": pipe_examples(metadata["examples"], 5),
            "value_type_or_key_fields": (
                f"fields={','.join(metadata['fields'])}; patients={len(patients)}; samples={len(samples)}; "
                f"tissues={','.join(tissues)}; metadata/GEO mapping=As_or_PB_to_PC,Li_to_LM,Ov_to_OM; {consistency}"
            ),
            "first_feature_examples": "",
        }
    )
    rows.append(
        {
            **base_row("GSE239676_count_matrix.mtx.gz", matrix["read_status"]),
            "n_rows": matrix["n_features"],
            "n_columns": matrix["n_cells"],
            "n_sample_or_cell_columns": matrix["n_cells"],
            "gene_or_feature_id_type": "feature_rows_linked_to_features_file",
            "metadata_or_sample_id_examples": pipe_examples(barcode_examples),
            "value_type_or_key_fields": (
                f"{matrix['banner']}; coordinate_integer; nonzero_entries={matrix['n_nonzero']}; {consistency}"
            ),
            "first_feature_examples": pipe_examples(feature_examples),
        }
    )
    rows.append(
        {
            **base_row("GSE239676_series_matrix.txt.gz", series["read_status"]),
            "n_rows": len(series["sample_ids"]),
            "n_columns": "",
            "n_sample_or_cell_columns": len(series["sample_ids"]),
            "gene_or_feature_id_type": "not_applicable",
            "metadata_or_sample_id_examples": pipe_examples(series["sample_ids"]),
            "value_type_or_key_fields": (
                "GEO metadata confirms 20 treatment-naive stage-IV patients with primary, liver-metastasis and "
                "peritoneal-carcinomatosis specimens; raw data withheld for privacy; " + consistency
            ),
            "first_feature_examples": "",
        }
    )
    return rows


def update_audit_table(root: Path, new_rows: List[Dict[str, object]]) -> Path:
    output = root / "results/F0_audit/predownloaded_resource_structure_audit.tsv"
    existing = read_tsv(output)
    replacement_names = {str(row["file_name"]) for row in new_rows}
    kept = [
        row
        for row in existing
        if not (row.get("dataset_id") == "GSE239676" and row.get("file_name") in replacement_names)
    ]
    for row in kept:
        row["sha256"] = normalize_sha256(row.get("sha256", ""))
    write_tsv(output, [*kept, *new_rows], AUDIT_FIELDS)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".", help="Project root; defaults to current directory.")
    return parser.parse_args()


def main() -> int:
    root = Path(parse_args().project_root).resolve()
    rows = build_rows(root)
    output = update_audit_table(root, rows)
    print(f"GSE239676 preaudit rows written: {len(rows)}")
    print(f"Output: {output}")
    for row in rows:
        print(f"  {row['file_name']}: {row['read_status']} | {row['audit_decision']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
