#!/usr/bin/env python3
"""Run the F0 data reconnaissance and reproducibility audit.

This script is intentionally standard-library only.  Formal execution requires
the explicit ``--execute`` flag; without it the script performs a dry run that
checks required inputs and prints the planned output paths.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import os
import platform
import re
import shutil
import sys
import tarfile
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


RANDOM_SEED = 42
RUN_ID = dt.datetime.now(dt.timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S")
F0_OUTPUTS = [
    "data/metadata/project_structure_ready.txt",
    "logs/F0_setup/analysis_log.md",
    "data/metadata/sample_info.tsv",
    "data/metadata/data_audit.tsv",
    "data/metadata/processed_input_manifest.tsv",
    "data/metadata/F0_dataset_inventory.tsv",
    "data/metadata/F0_file_manifest.tsv",
    "data/metadata/F0_metadata_field_inventory.tsv",
    "data/metadata/F0_author_processing_audit.tsv",
    "data/metadata/F0_data_readiness_by_F_section.tsv",
    "data/metadata/F0_external_resource_inventory.tsv",
    "data/metadata/F0_method_prior_decision.tsv",
    "data/metadata/decision_evidence_log.tsv",
    "data/metadata/excluded_samples.tsv",
    "results/F0_audit/F0_gate_checklist.tsv",
    "results/F0_audit/F0_global_data_reconnaissance_report.md",
    "results/F0_audit/F0_execution_report.md",
]

REQUIRED_DIRS = [
    "docs",
    "scripts/environment_setup",
    "logs/environment_setup",
    "reports/environment_setup",
    "environment",
    "data/public_downloads",
    "data/public_downloads/SCENIC_resources",
    "data/processed_input/GSE183904",
    "data/metadata",
    "scripts/F0_setup",
    "scripts/F1_single_cell",
    "scripts/F2_scoring",
    "scripts/F3_function",
    "scripts/F4_communication",
    "scripts/F5_bulk_marker",
    "scripts/F6_immunity",
    "scripts/F7_genomics_pathway",
    "scripts/F8_single_gene",
    "results/F0_audit",
    "results/F1_qc",
    "objects/F1_single_cell",
    "results/F1_annotation",
    "results/F1_malignancy",
    "results/F2_scoring",
    "results/F3_function",
    "results/F4_communication",
    "results/F5_bulk_marker",
    "results/F6_immunity",
    "results/F7_genomics_pathway",
    "results/F8_single_gene",
    "logs/F0_setup",
    "logs/F1_single_cell",
    "logs/F2_scoring",
    "logs/F3_function",
    "logs/F4_communication",
    "logs/F5_bulk_marker",
    "logs/F6_immunity",
    "logs/F7_genomics_pathway",
    "logs/F8_single_gene",
    "tmp",
]

REQUIRED_INPUTS = [
    "data/public_downloads/GSE183904_RAW.tar",
    "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
    "data/metadata/download_manifest.tsv",
    "data/metadata/preupload_resources_manifest.tsv",
    "data/metadata/preupload_pending_resources.tsv",
    "data/metadata/cell_type_marker_panel.tsv",
    "data/metadata/pipeline_parameters.yaml",
    "data/metadata/software_versions.tsv",
    "environment/execution_environment_inventory.tsv",
    "environment/environment_lock_manifest.tsv",
    "environment/random_seed_registry.tsv",
    "results/F0_audit/gse183904_csv_structure_precheck.tsv",
    "results/F0_audit/non_GSE183904_data_structure_precheck.tsv",
    "results/F0_audit/predownloaded_resource_structure_audit.tsv",
    "results/F0_audit/bulk_GEO_series_matrix_deep_precheck.tsv",
    "results/F0_audit/bulk_GEO_sample_overlap_precheck.tsv",
    "results/F0_audit/bulk_GEO_processing_metadata_precheck.tsv",
    "results/F0_audit/TCGA_STAD_survival_table_precheck.tsv",
    "results/F0_audit/TCGA_STAD_star_counts_inverse_validation.tsv",
    "results/F0_audit/GSE206785_dataset_structure_precheck.tsv",
    "results/F0_audit/GSE206785_metadata_precheck.tsv",
]


def rel(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: str(row.get(field, "")) for field in fields})
    tmp.replace(path)


def clean_csv_token(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1].replace('""', '"')
    return value


def append_log(root: Path, message: str) -> None:
    log_path = root / "logs/F0_setup/analysis_log.md"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8", newline="\n") as fh:
        fh.write(f"- {now_iso()} | {message}\n")


def parse_geo_line(line: str) -> Tuple[str, List[str]]:
    parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
    key = parts[0]
    values = [value.strip().strip('"') for value in parts[1:]]
    return key, values


def parse_gse183904_series(series_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    sample_fields: Dict[str, List[str]] = {}
    patient_by_sample: Dict[str, str] = {}
    mapping_header: List[str] | None = None
    in_mapping = False

    with gzip.open(series_path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for raw_line in fh:
            if not raw_line.startswith("!"):
                continue
            key, values = parse_geo_line(raw_line)
            if key.startswith("!Sample_"):
                sample_fields[key] = values
            if key == "!Series_summary" and values:
                text = values[0].strip()
                if text.startswith("Manuscript_ID"):
                    mapping_header = re.split(r"\s+", text)
                    in_mapping = True
                    continue
                if in_mapping and mapping_header and text:
                    tokens = re.split(r"\s+", text)
                    if len(tokens) >= 5 and tokens[0].startswith("NGC"):
                        patient = tokens[0].strip()
                        for sample_file in tokens[1:5]:
                            if sample_file != "-":
                                sample_id = sample_file.replace(".csv", "")
                                patient_by_sample[sample_id] = patient
                    elif text.startswith("!"):
                        in_mapping = False

    required = ["!Sample_geo_accession", "!Sample_title", "!Sample_source_name_ch1"]
    missing = [key for key in required if key not in sample_fields]
    if missing:
        raise RuntimeError(f"GSE183904 series matrix missing fields: {', '.join(missing)}")
    return sample_fields, patient_by_sample


def group_from_title(title: str) -> Tuple[str, str, str]:
    text = title.lower()
    if "primary gastric tissue" in text and "normal" in text:
        return "Normal_Gastric", "stomach", "normal"
    if "primary gastric tissue" in text and "tumor" in text:
        return "Primary_Tumor", "stomach", "tumor"
    if "periton" in text and "normal" in text:
        return "Normal_Peritoneum", "peritoneum", "normal"
    if "periton" in text and "tumor" in text:
        return "Peritoneal_Metastasis", "peritoneum", "tumor"
    return "Unclear", "unknown", "unknown"


def build_sample_info(
    sample_fields: Dict[str, List[str]],
    patient_by_sample: Dict[str, str],
    tar_members: Sequence[str],
) -> List[Dict[str, object]]:
    accessions = sample_fields["!Sample_geo_accession"]
    titles = sample_fields["!Sample_title"]
    sources = sample_fields.get("!Sample_source_name_ch1", [""] * len(accessions))
    chars = sample_fields.get("!Sample_characteristics_ch1", [""] * len(accessions))
    member_by_accession = {
        member.split("_", 1)[0]: member for member in tar_members if member.endswith(".csv.gz")
    }

    rows: List[Dict[str, object]] = []
    for idx, geo_accession in enumerate(accessions):
        title = titles[idx]
        sample_match = re.search(r"(sample\d+)", title)
        sample_id = sample_match.group(1) if sample_match else f"sample{idx + 1}"
        group, tissue_site, tumor_status = group_from_title(title)
        member_name = member_by_accession.get(geo_accession, f"{geo_accession}_{sample_id}.csv.gz")
        include_in_group = "false" if group == "Normal_Peritoneum" else "true"
        rows.append(
            {
                "geo_accession": geo_accession,
                "sample_id": sample_id,
                "sample_file": member_name,
                "patient_id": patient_by_sample.get(sample_id, "unknown"),
                "sample_title": title,
                "source_name_ch1": sources[idx] if idx < len(sources) else "",
                "sample_characteristics_ch1": chars[idx] if idx < len(chars) else "",
                "tissue_site": tissue_site,
                "tumor_status": tumor_status,
                "group_analysis": group,
                "source_of_group": "GEO_series_matrix_sample_title_and_characteristics",
                "metadata_confidence": "high" if group != "Unclear" else "low",
                "include_in_f1": "true" if group != "Unclear" else "pending",
                "include_in_group_comparison": include_in_group,
                "is_paired": "unknown",
                "pairing_scope": "patient_id_from_GEO_title_mapping_available_but_F1_pairing_not_preapproved",
                "note": (
                    "Normal_Peritoneum is reference/display only"
                    if group == "Normal_Peritoneum"
                    else "PM sample-level comparisons are directional only when group is Peritoneal_Metastasis"
                    if group == "Peritoneal_Metastasis"
                    else ""
                ),
            }
        )
    return rows


def tar_csv_members(tar_path: Path) -> List[tarfile.TarInfo]:
    with tarfile.open(tar_path, "r") as tf:
        members = [
            member
            for member in tf.getmembers()
            if member.isfile() and Path(member.name).name.endswith(".csv.gz")
        ]
    return sorted(members, key=lambda m: Path(m.name).name)


def extract_members(tar_path: Path, members: Sequence[tarfile.TarInfo], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        for member in members:
            target = out_dir / Path(member.name).name
            if not target.resolve().parent.samefile(out_dir.resolve()):
                raise RuntimeError(f"Unsafe tar member target: {member.name}")
            src = tf.extractfile(member)
            if src is None:
                raise RuntimeError(f"Could not read tar member: {member.name}")
            tmp = target.with_suffix(target.suffix + ".tmp")
            with src, tmp.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            tmp.replace(target)


def audit_csv_gz(path: Path) -> Dict[str, object]:
    start = time.perf_counter()
    rows = 0
    bad_column_count_rows = 0
    invalid_numeric_value_rows = 0
    missing_value_rows = 0
    mt_gene_count = 0
    hb_gene_count = 0
    mt_genes: List[str] = []
    hb_genes: List[str] = []
    first_genes: List[str] = []
    gene_hash = hashlib.sha256()

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        header = fh.readline()
        if not header:
            raise RuntimeError(f"Empty matrix file: {path}")
        header = header.rstrip("\n\r")
        header_fields = [clean_csv_token(value) for value in header.split(",")]
        expected_commas = len(header_fields) - 1
        first_header_fields = "|".join(header_fields[:6])
        barcode_suffixes = []
        for barcode in header_fields[1:6]:
            suffix = barcode.rsplit("_", 1)[-1] if "_" in barcode else ""
            barcode_suffixes.append(suffix)

        for line in fh:
            rows += 1
            stripped = line.rstrip("\n\r")
            comma_count = stripped.count(",")
            if comma_count != expected_commas:
                bad_column_count_rows += 1
            first_comma = stripped.find(",")
            if first_comma < 0:
                gene = clean_csv_token(stripped)
                values = ""
            else:
                gene = clean_csv_token(stripped[:first_comma])
                values = stripped[first_comma + 1 :]
            if len(first_genes) < 5:
                first_genes.append(gene)
            gene_hash.update(gene.encode("utf-8"))
            gene_hash.update(b"\n")
            upper_gene = gene.upper()
            if upper_gene.startswith("MT-"):
                mt_gene_count += 1
                if len(mt_genes) < 8:
                    mt_genes.append(gene)
            if re.match(r"^HB[A-Z0-9]", upper_gene) and not upper_gene.startswith(("HBS", "HBP")):
                hb_gene_count += 1
                if len(hb_genes) < 8:
                    hb_genes.append(gene)
            if values.startswith(",") or values.endswith(",") or ",," in values:
                missing_value_rows += 1
            if re.search(r"[^0-9,]", values):
                invalid_numeric_value_rows += 1

    suspected = "nonnegative_integer_count_like"
    if bad_column_count_rows or invalid_numeric_value_rows or missing_value_rows:
        suspected = "format_or_numeric_issue_detected"
    return {
        "matrix_rows_genes": rows,
        "matrix_cols_cells": expected_commas,
        "header_total_columns_including_gene_col": len(header_fields),
        "first_header_fields": first_header_fields,
        "first_genes": "|".join(first_genes),
        "barcode_suffix_examples": "|".join(barcode_suffixes),
        "bad_column_count_rows": bad_column_count_rows,
        "invalid_numeric_value_rows": invalid_numeric_value_rows,
        "missing_value_rows": missing_value_rows,
        "mt_gene_count": mt_gene_count,
        "mt_gene_examples": "|".join(mt_genes),
        "hb_gene_count": hb_gene_count,
        "hb_gene_examples": "|".join(hb_genes),
        "gene_order_sha256": gene_hash.hexdigest(),
        "suspected_matrix_type": suspected,
        "scan_seconds": round(time.perf_counter() - start, 2),
    }


def build_processed_manifest(
    root: Path,
    members: Sequence[tarfile.TarInfo],
    sample_info: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    by_file = {str(row["sample_file"]): row for row in sample_info}
    rows = []
    for member in members:
        name = Path(member.name).name
        path = root / "data/processed_input/GSE183904" / name
        sample = by_file.get(name, {})
        rows.append(
            {
                "source_archive": "data/public_downloads/GSE183904_RAW.tar",
                "archive_member_name": name,
                "extracted_path": rel(path, root),
                "geo_accession": sample.get("geo_accession", name.split("_", 1)[0]),
                "sample_id": sample.get("sample_id", ""),
                "file_size": path.stat().st_size if path.exists() else "",
                "sha256": sha256_file(path) if path.exists() else "",
                "extraction_date": now_iso(),
                "file_role": "author_filtered_raw_gene_count_csv_gz_for_F1",
                "note": "Compressed csv.gz retained; no permanent plain CSV extraction.",
            }
        )
    return rows


def build_data_audit(
    root: Path,
    processed_manifest: Sequence[Dict[str, object]],
    sample_info: Sequence[Dict[str, object]],
    precheck_rows: Sequence[Dict[str, str]],
) -> Tuple[List[Dict[str, object]], List[str]]:
    sample_by_file = {str(row["sample_file"]): row for row in sample_info}
    precheck_by_member = {row["member_name"]: row for row in precheck_rows}
    mismatch_notes: List[str] = []
    out_rows: List[Dict[str, object]] = []

    for manifest_row in processed_manifest:
        path = root / str(manifest_row["extracted_path"])
        stats = audit_csv_gz(path)
        member = str(manifest_row["archive_member_name"])
        sample = sample_by_file.get(member, {})
        pre = precheck_by_member.get(member, {})
        compare_fields = [
            "matrix_rows_genes",
            "matrix_cols_cells",
            "mt_gene_count",
            "hb_gene_count",
            "gene_order_sha256",
            "suspected_matrix_type",
        ]
        row_mismatches = []
        for field in compare_fields:
            observed = str(stats.get(field, ""))
            expected = str(pre.get(field, ""))
            if expected and observed != expected:
                row_mismatches.append(f"{field}: observed={observed}; precheck={expected}")
        pre_decision = pre.get("audit_decision_precheck", "")
        ok = (
            stats["suspected_matrix_type"] == "nonnegative_integer_count_like"
            and not row_mismatches
            and pre_decision == "enter_full_F1_candidate"
            and sample.get("include_in_f1") == "true"
        )
        if row_mismatches:
            mismatch_notes.append(f"{member} | " + " | ".join(row_mismatches))
        out_rows.append(
            {
                "geo_accession": manifest_row["geo_accession"],
                "sample_id": manifest_row["sample_id"],
                "sample_file": member,
                "extracted_path": manifest_row["extracted_path"],
                "matrix_orientation": "gene_rows_by_cell_columns",
                **stats,
                "precheck_audit_decision": pre_decision,
                "audit_decision": "enter_full_F1" if ok else "pause_for_review",
                "decision_scope": "file_format_only",
                "include_in_f1": sample.get("include_in_f1", "pending"),
                "raw_droplet_available": "false",
                "empty_droplet_background_available": "false",
                "normalization_artifact_flag": "false"
                if stats["suspected_matrix_type"] == "nonnegative_integer_count_like"
                else "true",
                "precheck_comparison_status": "match" if not row_mismatches else "mismatch",
                "note": "; ".join(row_mismatches),
            }
        )
    return out_rows, mismatch_notes


def summarize_unique(values: Iterable[str], max_items: int = 8) -> str:
    counts: Dict[str, int] = {}
    total = 0
    for value in values:
        total += 1
        key = value if value != "" else "<empty>"
        counts[key] = counts.get(key, 0) + 1
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    shown = [f"{key}:{count}" for key, count in ordered[:max_items]]
    if len(ordered) > max_items:
        shown.append(f"...plus_{len(ordered) - max_items}_more")
    return "|".join(shown) if shown else "none"


def build_metadata_inventory(
    sample_fields: Dict[str, List[str]],
    precheck_files: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for key, values in sample_fields.items():
        field = key.replace("!Sample_", "")
        missing = sum(1 for value in values if value == "")
        rows.append(
            {
                "dataset_id": "GSE183904",
                "source_file": "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
                "field_name": field,
                "field_category": "sample_metadata",
                "n_records_checked": len(values),
                "unique_values_summary": summarize_unique(values),
                "missing_count": missing,
                "missing_rate": round(missing / len(values), 6) if values else "",
                "example_values": "|".join(values[:5]),
                "usable_for_F_section": "F0,F1",
                "interpretation_risk": "patient_pairing_requires_title_mapping_review"
                if field in {"title", "characteristics_ch1"}
                else "low",
                "note": "Current public matrix files do not provide cell-level metadata.",
            }
        )
    for row in precheck_files.get("GSE206785_metadata_precheck", []):
        rows.append(
            {
                "dataset_id": row.get("dataset_id", "GSE206785"),
                "source_file": row.get("file_name", ""),
                "field_name": row.get("column", ""),
                "field_category": "cell_metadata",
                "n_records_checked": row.get("n_rows", ""),
                "unique_values_summary": row.get("top_values", ""),
                "missing_count": row.get("missing", ""),
                "missing_rate": "",
                "example_values": row.get("top_values", ""),
                "usable_for_F_section": "F2_external_sc_validation_after_approval",
                "interpretation_risk": "processed_external_scRNA_no_PM_group"
                if row.get("column") in {"Group", "Sample"}
                else "section_specific_review_required",
                "note": row.get("precheck_note", ""),
            }
        )
    for row in precheck_files.get("bulk_GEO_series_matrix_deep_precheck", []):
        rows.append(
            {
                "dataset_id": row.get("dataset_id", ""),
                "source_file": row.get("file_name", ""),
                "field_name": "characteristics_keys",
                "field_category": "clinical",
                "n_records_checked": row.get("n_samples", ""),
                "unique_values_summary": row.get("characteristics_keys", ""),
                "missing_count": "",
                "missing_rate": "",
                "example_values": row.get("clinical_like_characteristics_examples", ""),
                "usable_for_F_section": "F2.4,F5_after_dataset_specific_audit",
                "interpretation_risk": "processed_microarray_scale_and_endpoint_audit_required",
                "note": row.get("method_implication", ""),
            }
        )
    return rows


def build_author_processing_audit() -> List[Dict[str, object]]:
    return [
        {
            "dataset_id": "GSE183904",
            "source_reference_or_file": "GSE183904_series_matrix.txt.gz !Sample_data_processing",
            "processing_step": "QC_filtering",
            "author_reported_status": "reported",
            "method_or_threshold_if_reported": "Seurat v3; genes/features shared by >=3 cells; 500<=nFeature<6000; percent.mt<=20",
            "evidence_location": "GEO sample data processing lines",
            "confidence_level": "high",
            "implication_for_downstream_plan": "F1 should verify author-filtered raw count matrices and perform conservative residual QC rather than pretending raw droplets are available.",
            "requires_special_handling": "record baseline thresholds in F1 qc_thresholds_by_sample.tsv",
        },
        {
            "dataset_id": "GSE183904",
            "source_reference_or_file": "GSE183904_series_matrix.txt.gz !Series_overall_design",
            "processing_step": "raw_FASTQ_or_empty_droplet_availability",
            "author_reported_status": "raw_files_not_submitted",
            "method_or_threshold_if_reported": "Raw files were not submitted due to patient privacy concerns.",
            "evidence_location": "GEO series overall design",
            "confidence_level": "high",
            "implication_for_downstream_plan": "Do not claim FASTQ/BCL/Cell Ranger raw output reprocessing for GSE183904; SoupX is not default without background droplets.",
            "requires_special_handling": "raw_droplet_available=false; empty_droplet_background_available=false",
        },
        {
            "dataset_id": "GSE183904",
            "source_reference_or_file": "GSE183904_series_matrix.txt.gz !Sample_data_processing",
            "processing_step": "count_generation",
            "author_reported_status": "reported",
            "method_or_threshold_if_reported": "Cell Ranger v3.0 aligned FASTQ to hg38 and generated single-cell feature counts.",
            "evidence_location": "GEO sample data processing lines",
            "confidence_level": "high",
            "implication_for_downstream_plan": "Public CSV.gz files are author-filtered raw gene count matrices, not normalized expression.",
            "requires_special_handling": "F0 data_audit.tsv must verify nonnegative integer matrix structure.",
        },
        {
            "dataset_id": "GSE183904",
            "source_reference_or_file": "GSE183904_series_matrix.txt.gz and available public files",
            "processing_step": "doublet_detection",
            "author_reported_status": "not_reported_in_public_metadata",
            "method_or_threshold_if_reported": "not_available",
            "evidence_location": "No public doublet field found in current F0 inputs.",
            "confidence_level": "medium",
            "implication_for_downstream_plan": "F1 Gate1 must independently assess doublet risk; scDblFinder by sample is the preregistered main method if installed/approved.",
            "requires_special_handling": "Do not treat absence of public doublet metadata as evidence that no doublets exist.",
        },
        {
            "dataset_id": "GSE183904",
            "source_reference_or_file": "GSE183904_series_matrix.txt.gz and available public files",
            "processing_step": "ambient_RNA_correction",
            "author_reported_status": "not_reported_in_public_metadata",
            "method_or_threshold_if_reported": "not_available",
            "evidence_location": "No public ambient correction field found in current F0 inputs.",
            "confidence_level": "medium",
            "implication_for_downstream_plan": "F1 Gate1 should evaluate cross-expression and decontX/celda feasibility; SoupX is not default without empty droplets.",
            "requires_special_handling": "Report contamination estimates separately from corrected-count usage.",
        },
        {
            "dataset_id": "GSE206785",
            "source_reference_or_file": "GSE206785_dataset_structure_precheck.tsv",
            "processing_step": "external_scRNA_processing",
            "author_reported_status": "processed_matrix_local_candidate",
            "method_or_threshold_if_reported": "log1p(count)-like cell-by-gene matrix",
            "evidence_location": "results/F0_audit/GSE206785_dataset_structure_precheck.tsv",
            "confidence_level": "medium",
            "implication_for_downstream_plan": "Candidate for processed-expression external scoring only after approval; not for raw QC/doublet/ambient reprocessing.",
            "requires_special_handling": "F2 external validation audit required before use.",
        },
    ]


def build_dataset_inventory(
    sample_info: Sequence[Dict[str, object]],
    precheck_files: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, object]]:
    group_counts: Dict[str, int] = {}
    for row in sample_info:
        group = str(row["group_analysis"])
        group_counts[group] = group_counts.get(group, 0) + 1
    group_summary = "|".join(f"{k}:{v}" for k, v in sorted(group_counts.items()))
    rows: List[Dict[str, object]] = [
        {
            "dataset_id": "GSE183904",
            "dataset_name_or_accession": "GSE183904",
            "data_domain": "scRNA",
            "intended_F_sections": "F1,F2,F3,F4,F8",
            "current_availability_status": "project_local_available",
            "primary_or_supplementary_role": "primary_scRNA_discovery",
            "sample_or_cell_count_if_known": f"40 samples; groups {group_summary}",
            "patient_count_if_known": "31 title-mapped manuscript IDs if all mappings pass review",
            "event_count_if_known": "not_applicable",
            "platform_or_assay": "10x Genomics single-cell RNA-seq; GPL24676; hg38",
            "gene_id_type_if_known": "gene_symbol_or_feature_name_from_author_CSV",
            "raw_or_processed_status_if_known": "author-filtered raw gene count CSV.gz; not FASTQ/raw droplets",
            "metadata_available": "GEO sample metadata available; no public cell-level metadata in current inputs",
            "clinical_endpoint_available": "not_for_survival",
            "major_limitations": "PM n=3; Normal_Peritoneum n=1; raw droplets and empty droplets unavailable",
            "next_required_audit_step": "F0 gate review, then F1 Gate1 QC/doublet/ambient assessment after approval",
        }
    ]
    for row in precheck_files.get("non_GSE183904", []):
        rows.append(
            {
                "dataset_id": row.get("dataset_or_resource_id", ""),
                "dataset_name_or_accession": row.get("dataset_or_resource_id", ""),
                "data_domain": row.get("data_domain", ""),
                "intended_F_sections": row.get("intended_F_sections", ""),
                "current_availability_status": "project_local_available"
                if row.get("local_relative_path")
                else "unknown",
                "primary_or_supplementary_role": "candidate_or_supporting_resource",
                "sample_or_cell_count_if_known": row.get("n_sample_or_value_columns", ""),
                "patient_count_if_known": "",
                "event_count_if_known": "",
                "platform_or_assay": row.get("platform_or_release", ""),
                "gene_id_type_if_known": row.get("metadata_or_annotation_columns", ""),
                "raw_or_processed_status_if_known": row.get("value_type_or_expression_unit", ""),
                "metadata_available": row.get("data_structure", ""),
                "clinical_endpoint_available": "section_specific_audit_required",
                "major_limitations": row.get("known_limitations", ""),
                "next_required_audit_step": row.get("server_executor_required_followup", ""),
            }
        )
    for row in precheck_files.get("predownloaded", []):
        rows.append(
            {
                "dataset_id": row.get("dataset_id", ""),
                "dataset_name_or_accession": row.get("dataset_id", ""),
                "data_domain": "preloaded_resource",
                "intended_F_sections": "section_specific",
                "current_availability_status": "project_local_available"
                if row.get("download_status") == "complete"
                else row.get("download_status", "unknown"),
                "primary_or_supplementary_role": row.get("file_role", ""),
                "sample_or_cell_count_if_known": row.get("n_sample_or_cell_columns", ""),
                "patient_count_if_known": "",
                "event_count_if_known": "",
                "platform_or_assay": "",
                "gene_id_type_if_known": row.get("gene_or_feature_id_type", ""),
                "raw_or_processed_status_if_known": row.get("value_type_or_key_fields", ""),
                "metadata_available": row.get("metadata_or_sample_id_examples", ""),
                "clinical_endpoint_available": "section_specific_audit_required",
                "major_limitations": row.get("major_limitations", ""),
                "next_required_audit_step": row.get("method_implication", ""),
            }
        )
    return rows


def build_file_manifest(
    root: Path,
    processed_manifest: Sequence[Dict[str, object]],
    generated_paths: Sequence[str],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for source_path, source_name in [
        ("data/metadata/download_manifest.tsv", "download_manifest"),
        ("data/metadata/preupload_resources_manifest.tsv", "preupload_resources"),
    ]:
        for row in read_tsv(root / source_path):
            relative = row.get("relative_path_under_data_public_downloads", "")
            if relative:
                path = root / "data/public_downloads" / relative
            else:
                file_name = row.get("file_name", "")
                path = root / "data/public_downloads" / file_name
            rows.append(
                {
                    "file_name": row.get("file_name", ""),
                    "relative_path_if_available": rel(path, root) if path.exists() else relative,
                    "dataset_id": row.get("dataset_id", ""),
                    "file_role": row.get("file_role", ""),
                    "data_domain": "mixed",
                    "source_type": "researcher_provided_local"
                    if source_name == "download_manifest"
                    else "locally_precached",
                    "source_url_or_local_manifest": row.get("source_url", row.get("source_url_or_derivation", "")),
                    "file_size_bytes": row.get("file_size", row.get("file_size_bytes", "")),
                    "sha256": row.get("sha256", ""),
                    "compression_format": Path(row.get("file_name", "")).suffix.lstrip("."),
                    "read_status": row.get("read_status", "not_checked_in_F0_file_manifest"),
                    "availability_status": "project_local_available" if path.exists() else "not_found_at_expected_path",
                    "used_in_F0": "true"
                    if row.get("dataset_id") in {"GSE183904"} or row.get("file_name") == "GSE183904_series_matrix.txt.gz"
                    else "inventory_only",
                    "planned_F_section_use": "F1-F8",
                    "audit_status": "recorded_from_existing_manifest",
                    "artifact_class": "audit_trail",
                    "publication_destination": "not_planned",
                    "review_priority": "standard",
                    "note": row.get("note", ""),
                }
            )
    for row in read_tsv(root / "data/metadata/preupload_pending_resources.tsv"):
        rows.append(
            {
                "file_name": row.get("resource_name", ""),
                "relative_path_if_available": "",
                "dataset_id": row.get("resource_group", ""),
                "file_role": "pending_resource",
                "data_domain": "mixed",
                "source_type": "pending",
                "source_url_or_local_manifest": "data/metadata/preupload_pending_resources.tsv",
                "file_size_bytes": "",
                "sha256": "",
                "compression_format": "",
                "read_status": "not_available",
                "availability_status": "pending",
                "used_in_F0": "inventory_only",
                "planned_F_section_use": "section_specific",
                "audit_status": "pending_not_used",
                "artifact_class": "audit_trail",
                "publication_destination": "not_planned",
                "review_priority": "spot_check",
                "note": row.get("recommended_next_step", ""),
            }
        )
    for row in processed_manifest:
        rows.append(
            {
                "file_name": row.get("archive_member_name", ""),
                "relative_path_if_available": row.get("extracted_path", ""),
                "dataset_id": "GSE183904",
                "file_role": row.get("file_role", ""),
                "data_domain": "scRNA",
                "source_type": "generated_by_F0",
                "source_url_or_local_manifest": row.get("source_archive", ""),
                "file_size_bytes": row.get("file_size", ""),
                "sha256": row.get("sha256", ""),
                "compression_format": "gz",
                "read_status": "gzip_stream_audited",
                "availability_status": "project_local_available",
                "used_in_F0": "true",
                "planned_F_section_use": "F1",
                "audit_status": "generated_and_audited_by_F0",
                "artifact_class": "core",
                "publication_destination": "not_planned",
                "review_priority": "full",
                "note": row.get("note", ""),
            }
        )
    for out in generated_paths:
        path = root / out
        rows.append(
            {
                "file_name": Path(out).name,
                "relative_path_if_available": out,
                "dataset_id": "F0",
                "file_role": "F0_generated_output",
                "data_domain": "audit",
                "source_type": "generated_by_F0",
                "source_url_or_local_manifest": "scripts/F0_setup/run_F0_full_audit.py",
                "file_size_bytes": path.stat().st_size if path.exists() else "",
                "sha256": sha256_file(path) if path.exists() and Path(out).name != "F0_file_manifest.tsv" else "",
                "compression_format": Path(out).suffix.lstrip("."),
                "read_status": "written",
                "availability_status": "project_local_available" if path.exists() else "planned",
                "used_in_F0": "true",
                "planned_F_section_use": "F0,F1-F8_as_relevant",
                "audit_status": "generated_by_current_script",
                "artifact_class": "gate_decision"
                if Path(out).name in {"F0_gate_checklist.tsv", "F0_execution_report.md"}
                else "audit_trail",
                "publication_destination": "undecided"
                if Path(out).name in {"F0_global_data_reconnaissance_report.md", "F0_execution_report.md"}
                else "not_planned",
                "review_priority": "standard",
                "note": "Self-manifest sha256 intentionally blank for F0_file_manifest.tsv"
                if Path(out).name == "F0_file_manifest.tsv"
                else "",
            }
        )
    return rows


def build_data_readiness() -> List[Dict[str, object]]:
    return [
        {
            "F_section": "F1",
            "required_data_domains": "GSE183904 scRNA count matrices; sample metadata; marker panel",
            "available_datasets_or_files": "GSE183904_RAW.tar; GSE183904_series_matrix.txt.gz; cell_type_marker_panel.tsv",
            "missing_or_pending_items": "F1 R dependency fill still required for scDblFinder/SoupX_or_celda/CopyKAT as selected later",
            "minimum_data_needed_to_start": "F0_scRNA_F1_gate PASS or PASS_WITH_NOTED_LIMITATIONS",
            "specialized_audit_required_in_section": "Gate1 QC/doublet/ambient by sample",
            "current_readiness_status": "ready_with_limitations",
            "blocking_for_next_section": "yes_until_F0_gate_passes",
            "recommended_first_action": "Run F1 Gate1 plan only after F0 review and user approval",
            "note": "F1 must not use MLMOD score/signature/prognosis.",
        },
        {
            "F_section": "F2",
            "required_data_domains": "F1 objects; GSE235046/SRP444325 signature source; bulk/external sc validation candidates",
            "available_datasets_or_files": "GSE235046 metadata/count-like table; SRP444325 download manifest; GSE206785/GSE239676 candidates; TCGA/bulk GEO candidates",
            "missing_or_pending_items": "Signature must be frozen before validation; SRP444325 raw reprocessing requires separate approved plan",
            "minimum_data_needed_to_start": "F1 approved malignant/candidate cell object plus F2.1 plan",
            "specialized_audit_required_in_section": "SRA reprocessing or approved fallback; bulk/external sc cohort audit",
            "current_readiness_status": "pending_local_copy_or_download",
            "blocking_for_next_section": "no_for_F1; yes_for_F2_start",
            "recommended_first_action": "After F1, audit signature-source route and external validation isolation",
            "note": "Do not round fractional GSE235046 table for primary DESeq2.",
        },
        {
            "F_section": "F3",
            "required_data_domains": "F2 frozen states; F1 malignant epithelial object; SCENIC resources",
            "available_datasets_or_files": "SCENIC resource manifest and candidate local files",
            "missing_or_pending_items": "SCENIC compatibility and gene overlap audit",
            "minimum_data_needed_to_start": "F2 candidate/high-low definitions approved",
            "specialized_audit_required_in_section": "pySCENIC/ctxcore compatibility and resource version audit",
            "current_readiness_status": "pending_local_copy_or_download",
            "blocking_for_next_section": "no_for_F1; yes_for_F3_start",
            "recommended_first_action": "Run resource compatibility smoke test before GRN plan",
            "note": "SCENIC smoke test is technical only, not a biological screen.",
        },
        {
            "F_section": "F4",
            "required_data_domains": "F1 annotated object; F2 candidate state; ligand-receptor resources",
            "available_datasets_or_files": "No formal LR output yet",
            "missing_or_pending_items": "CellChat/LIANA resources and sample-aware pseudobulk plan",
            "minimum_data_needed_to_start": "F2 state approved and F1 cell type labels approved",
            "specialized_audit_required_in_section": "LR database/version and sample-aware expression support audit",
            "current_readiness_status": "unknown",
            "blocking_for_next_section": "no_for_F1; yes_for_F4_start",
            "recommended_first_action": "Prepare F4 LR method plan after F2 gate approval",
            "note": "CellChat network is background; sample-aware LR expression is primary.",
        },
        {
            "F_section": "F5",
            "required_data_domains": "Bulk expression and clinical/survival endpoints",
            "available_datasets_or_files": "TCGA-STAD Xena/GDC/cBioPortal candidates; GEO/ACRG candidate prechecks",
            "missing_or_pending_items": "Dataset-specific clinical endpoint and expression-scale audit",
            "minimum_data_needed_to_start": "F2 clinical/signature evidence state and approved bulk cohort audit",
            "specialized_audit_required_in_section": "Probe mapping, expression scale, survival endpoint, overlap audit",
            "current_readiness_status": "ready_with_limitations",
            "blocking_for_next_section": "no_for_F1; yes_for_F5_start",
            "recommended_first_action": "Freeze F5 cohort ranking and endpoint availability before marker modeling",
            "note": "SuperSeries overlap must be handled; processed arrays are not count-model inputs.",
        },
        {
            "F_section": "F6",
            "required_data_domains": "Bulk expression; F1 reference; F5 genes as applicable",
            "available_datasets_or_files": "TCGA/bulk candidates and GSE183904 reference candidate",
            "missing_or_pending_items": "BayesPrism benchmark plan; TIDE/IPS input/output approvals",
            "minimum_data_needed_to_start": "F5 core genes and approved deconvolution method plan",
            "specialized_audit_required_in_section": "Simulation/LOPO benchmark and immune method legality/version audit",
            "current_readiness_status": "unknown",
            "blocking_for_next_section": "no_for_F1; yes_for_F6_start",
            "recommended_first_action": "Write F6_deconvolution_method_plan.tsv and pause",
            "note": "No web/API immune outputs should be generated in F0.",
        },
        {
            "F_section": "F7",
            "required_data_domains": "TCGA expression, mutation, CNV, methylation; F5 allowed genes",
            "available_datasets_or_files": "TCGA star_counts/star_tpm, survival, clinical, CNV, MAF, 450K candidates",
            "missing_or_pending_items": "MAF merge/TMB denominator; methylation probe filtering; expression overlap",
            "minimum_data_needed_to_start": "F5 allowed_for_F7 genes and F2/F5 interpretation boundaries",
            "specialized_audit_required_in_section": "Multiomics sample intersection and data-type-specific QC",
            "current_readiness_status": "ready_with_limitations",
            "blocking_for_next_section": "no_for_F1; yes_for_F7_start",
            "recommended_first_action": "Run F7 specialized TCGA audit before modeling",
            "note": "TCGA star_counts reverse validation permits VST/WGCNA after F7 audit, not automatic DE.",
        },
        {
            "F_section": "F8",
            "required_data_domains": "F5 selected gene(s); F1/F2 objects; optional perturbation resources",
            "available_datasets_or_files": "Depends on F5/F2 outputs",
            "missing_or_pending_items": "No final gene selected yet",
            "minimum_data_needed_to_start": "F5 approved final gene list and F2/F3/F4 context",
            "specialized_audit_required_in_section": "Single-gene evidence and perturbation model reliability audit",
            "current_readiness_status": "not_ready",
            "blocking_for_next_section": "no_for_F1; yes_for_F8_start",
            "recommended_first_action": "Wait for upstream outputs",
            "note": "F8.3 evidence ceiling remains model_supported_hypothesis.",
        },
    ]


def build_external_resource_inventory(root: Path) -> List[Dict[str, object]]:
    rows = []
    for row in read_tsv(root / "data/metadata/preupload_resources_manifest.tsv"):
        dataset = row.get("dataset_id", "")
        if dataset.startswith("SCENIC") or dataset in {"inferCNV_gene_order"} or "SCENIC" in row.get("relative_path_under_data_public_downloads", ""):
            rows.append(
                {
                    "resource_name": row.get("file_name", ""),
                    "resource_type": dataset,
                    "planned_F_sections": "F3" if dataset.startswith("SCENIC") else "F1",
                    "current_availability_status": "project_local_available"
                    if row.get("exists") == "TRUE"
                    else "not_available",
                    "source_url_or_database": row.get("source_url_or_derivation", ""),
                    "version_or_release_if_known": "hg38_v10" if "SCENIC" in dataset else "",
                    "license_or_login_requirement": "public_resource_verify_terms_before_publication",
                    "local_manifest_if_any": "data/metadata/preupload_resources_manifest.tsv",
                    "sha256_if_file": row.get("sha256", ""),
                    "large_or_unstable_download_risk": "yes"
                    if int(row.get("file_size_bytes", "0") or "0") > 100_000_000
                    else "no",
                    "next_required_audit_step": "section_specific_compatibility_check",
                    "note": row.get("note", ""),
                }
            )
    for resource_name, resource_type, sections, note in [
        ("MSigDB_GOBP_FERROPTOSIS", "external_database", "F2,F3", "Version and license audit before enrichment/scoring."),
        ("Reactome", "external_database", "F3,F7", "Version audit before pathway interpretation."),
        ("KEGG", "external_database", "F3,F7", "Access and license constraints must be checked before use."),
        ("JASPAR/DoRothEA", "external_database", "F3,F8", "Use only after GRN method plan approval."),
        ("TIDE/IPS/SubMap", "web_or_API_resource", "F6", "Prepare only after final gene/signature inputs are frozen."),
    ]:
        rows.append(
            {
                "resource_name": resource_name,
                "resource_type": resource_type,
                "planned_F_sections": sections,
                "current_availability_status": "not_yet_audited",
                "source_url_or_database": resource_name,
                "version_or_release_if_known": "",
                "license_or_login_requirement": "requires_section_specific_audit",
                "local_manifest_if_any": "",
                "sha256_if_file": "",
                "large_or_unstable_download_risk": "unknown",
                "next_required_audit_step": "audit when section starts; do not use in F0",
                "note": note,
            }
        )
    return rows


def build_method_prior_decision() -> List[Dict[str, object]]:
    return [
        {
            "dataset_or_resource_id": "GSE183904",
            "intended_F_sections": "F1-F4,F8",
            "known_data_structure": "author-filtered raw gene count CSV.gz matrices; gene rows by cell columns",
            "known_limitation_from_precheck": "No FASTQ/raw droplets/empty droplets; PM n=3; Normal_Peritoneum n=1",
            "default_or_recommended_route": "F1 conservative QC reanalysis using author-filtered counts; verify QC, doublet and ambient risk",
            "method_not_allowed_by_default": "Do not claim FASTQ/Cell Ranger raw reprocessing; do not force SoupX without background droplets; do not use MLMOD in F1",
            "required_pre_execution_audit": "F0 data_audit and F1 Gate1 QC/doublet/ambient plan",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "If count structure conflicts with precheck, pause and revise F1 input strategy",
            "interpretation_boundary": "F1 produces a reusable annotated object, not biological MLMOD conclusions",
            "note": "Normal_Peritoneum reference-only; PM sample-level statistics directional only.",
        },
        {
            "dataset_or_resource_id": "GSE235046/SRP444325",
            "intended_F_sections": "F2.1",
            "known_data_structure": "Mouse BMDM RNA-seq; public GEO count-like table has decimal values; ENA paired FASTQ route prepared",
            "known_limitation_from_precheck": "GEO count-like table is not all integer raw counts",
            "default_or_recommended_route": "Approved SRA reprocessing STAR/RSEM/tximport/DESeq2 or equivalent auditable count model",
            "method_not_allowed_by_default": "Do not round public decimal table for primary DESeq2",
            "required_pre_execution_audit": "F2.1 raw reprocessing resource and design audit",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "TMM + limma-voom/limma conditional fallback if raw route not approved or fails",
            "interpretation_boundary": "Cross-species macrophage-to-epithelial transfer remains a mechanistic hypothesis requiring controls",
            "note": "Interaction/Torin annotations do not filter main signature membership.",
        },
        {
            "dataset_or_resource_id": "TCGA-STAD_Xena_GDC_cBioPortal",
            "intended_F_sections": "F2.4,F5,F7",
            "known_data_structure": "Xena star_counts stored as log2(count+1) but inverse-validatable; TPM, clinical, CNV, MAF, 450K candidate files exist",
            "known_limitation_from_precheck": "Stored star_counts are not directly raw integer values; multiomics files need sample intersection and QC",
            "default_or_recommended_route": "Use processed expression for scoring/Cox/correlation; reconstructed counts for VST/WGCNA only after F7 audit",
            "method_not_allowed_by_default": "Do not use stored log2 values as raw-count DE input; DESeq2 differential expression requires separate approval",
            "required_pre_execution_audit": "F2.4/F5/F7 specialized expression, clinical, mutation, CNV and methylation audits",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Download direct GDC raw_count only if later approved DE requires it",
            "interpretation_boundary": "Bulk association supports prognosis/clinical association, not single-cell causal mechanism",
            "note": "TCGA_STAD_star_counts_inverse_validation.tsv must be rechecked in F0/F7.",
        },
        {
            "dataset_or_resource_id": "bulk_GEO_ACRG_candidates",
            "intended_F_sections": "F2.4,F5",
            "known_data_structure": "Processed microarray series matrices with heterogeneous scale and overlapping SuperSeries risk",
            "known_limitation_from_precheck": "GSE15459/GSE84426 linear-like scales; GSE66229/GSE26253 log-like scales; sample overlap risks",
            "default_or_recommended_route": "Probe-to-gene audit, within-cohort scoring/limma/Cox/KM as appropriate",
            "method_not_allowed_by_default": "Do not use count models; do not combine cohorts before scale/probe/overlap audit",
            "required_pre_execution_audit": "Expression scale, platform annotation, clinical endpoint and sample overlap audit",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Use only cohorts with sufficient endpoints and non-overlap after audit",
            "interpretation_boundary": "Validation strength depends on independent patients and endpoint completeness",
            "note": "GSE62254 subset of GSE66229; GSE84426 subset/same-release risk with GSE84437.",
        },
        {
            "dataset_or_resource_id": "GSE206785",
            "intended_F_sections": "F2_external_sc_validation",
            "known_data_structure": "Cell-by-gene log1p(count)-like processed matrix with metadata",
            "known_limitation_from_precheck": "No PM group; not raw counts; no raw QC/doublet/ambient reprocessing",
            "default_or_recommended_route": "Processed-expression MLMOD scoring and broad localization validation after F2 freeze",
            "method_not_allowed_by_default": "Do not use integer-count pseudobulk models or PM enrichment validation",
            "required_pre_execution_audit": "External sc cohort audit after F2 signature/high-low rules are frozen",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Use as not_evaluable if gene coverage or annotation mapping fails",
            "interpretation_boundary": "Can support cross-dataset scoring transfer, not raw-processing reproducibility",
            "note": "Must not tune F2 signature using GSE206785.",
        },
        {
            "dataset_or_resource_id": "SCENIC_hg38_v10",
            "intended_F_sections": "F3.3",
            "known_data_structure": "TF list, motif-to-TF table, rankings.feather resources",
            "known_limitation_from_precheck": "Compatibility with pyarrow/ctxcore/pySCENIC not yet audited",
            "default_or_recommended_route": "Run technical compatibility smoke test before full GRNBoost2 plan",
            "method_not_allowed_by_default": "Do not use smoke test as biological filter or downsample based on MLMOD",
            "required_pre_execution_audit": "Resource version, gene overlap and multi-seed GRN stability plan",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Rent approved server if full SCENIC resources exceed local capacity",
            "interpretation_boundary": "GRN perturbation supports model-supported hypotheses only",
            "note": "Formal GRNBoost2 requires at least 10 seeds if run.",
        },
    ]


def build_decision_evidence_log() -> List[Dict[str, object]]:
    return [
        {
            "decision_id": "F0_DECISION_001",
            "stage": "F0",
            "decision_topic": "primary_scRNA_dataset",
            "decision_value": "GSE183904 is the primary discovery scRNA dataset",
            "decision_reason": "Project hypothesis and main plan are built around gastric cancer scRNA reanalysis with 40 local sample matrices.",
            "evidence_type": "local_file_and_GEO_metadata",
            "evidence_source": "GSE183904_RAW.tar; GSE183904_series_matrix.txt.gz",
            "source_url_or_file": "data/public_downloads/GSE183904_RAW.tar",
            "evidence_strength": "high_for_file_availability",
            "date": now_iso(),
            "requires_sensitivity_analysis": "no",
            "note": "Scientific conclusions require downstream gates.",
        },
        {
            "decision_id": "F0_DECISION_002",
            "stage": "F0/F1",
            "decision_topic": "GSE183904_matrix_boundary",
            "decision_value": "Use author-filtered raw gene count matrices, not FASTQ/raw droplets",
            "decision_reason": "GEO metadata reports raw files not submitted and supplementary matrix table contains raw gene counts after Cell Ranger and Seurat filtering.",
            "evidence_type": "GEO_metadata_and_stream_audit",
            "evidence_source": "GSE183904_series_matrix.txt.gz; data_audit.tsv",
            "source_url_or_file": "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
            "evidence_strength": "high",
            "date": now_iso(),
            "requires_sensitivity_analysis": "yes",
            "note": "F1 must verify residual QC/doublet/ambient risk.",
        },
        {
            "decision_id": "F0_DECISION_003",
            "stage": "F0/F1",
            "decision_topic": "group_handling",
            "decision_value": "Normal_Peritoneum is reference-only; Peritoneal_Metastasis is tumor but PM n=3 is directional only",
            "decision_reason": "GSE183904 sample titles show one normal peritoneum and three peritoneal tumor samples.",
            "evidence_type": "GEO_sample_metadata",
            "evidence_source": "sample_info.tsv",
            "source_url_or_file": "data/metadata/sample_info.tsv",
            "evidence_strength": "high_for_group_labels",
            "date": now_iso(),
            "requires_sensitivity_analysis": "yes",
            "note": "PM cannot be treated as independent validation source.",
        },
        {
            "decision_id": "F0_DECISION_004",
            "stage": "F0/F1",
            "decision_topic": "marker_panel_lock",
            "decision_value": "cell_type_marker_panel.tsv is read-only preregistered annotation input",
            "decision_reason": "Project instructions forbid modifying marker panel during F0/F1 execution; issues must be reported separately.",
            "evidence_type": "project_contract",
            "evidence_source": "AGENTS.md and main plan",
            "source_url_or_file": "data/metadata/cell_type_marker_panel.tsv",
            "evidence_strength": "project_rule",
            "date": now_iso(),
            "requires_sensitivity_analysis": "no",
            "note": "Script validates basic structure only.",
        },
    ]


def build_excluded_samples(sample_info: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    excluded = [row for row in sample_info if row.get("include_in_f1") != "true"]
    if not excluded:
        return [
            {
                "sample_id": "none",
                "geo_accession": "none",
                "exclusion_scope": "F1_object_construction",
                "exclusion_status": "no_exclusion",
                "reason": "All 40 GSE183904 samples are eligible for F1 object construction if data_audit.tsv remains PASS.",
                "evidence_file": "data/metadata/sample_info.tsv; data/metadata/data_audit.tsv",
                "note": "Normal_Peritoneum is excluded from main group comparison but not from F1 object construction by default.",
            }
        ]
    rows = []
    for row in excluded:
        rows.append(
            {
                "sample_id": row.get("sample_id", ""),
                "geo_accession": row.get("geo_accession", ""),
                "exclusion_scope": "F1_object_construction",
                "exclusion_status": "pending_or_excluded",
                "reason": "F0 group or metadata unresolved",
                "evidence_file": "data/metadata/sample_info.tsv",
                "note": row.get("note", ""),
            }
        )
    return rows


def build_gate_checklist(
    root: Path,
    sample_info: Sequence[Dict[str, object]],
    data_audit: Sequence[Dict[str, object]],
    processed_manifest: Sequence[Dict[str, object]],
    mismatches: Sequence[str],
) -> List[Dict[str, object]]:
    tar_path = root / "data/public_downloads/GSE183904_RAW.tar"
    rows = []

    def add(item: str, required: str, observed: str, ok: bool, level: str, evidence: str, note: str = "") -> None:
        rows.append(
            {
                "gate_item": item,
                "required_status": required,
                "observed_status": observed,
                "pass_fail": "PASS" if ok else "FAIL",
                "blocking_level": level,
                "evidence_file": evidence,
                "note": note,
            }
        )

    add(
        "project_structure_ready",
        "project_structure_ready.txt exists",
        "exists" if (root / "data/metadata/project_structure_ready.txt").exists() else "missing",
        (root / "data/metadata/project_structure_ready.txt").exists(),
        "blocking",
        "data/metadata/project_structure_ready.txt",
    )
    add(
        "analysis_log",
        "logs/F0_setup/analysis_log.md exists",
        "exists" if (root / "logs/F0_setup/analysis_log.md").exists() else "missing",
        (root / "logs/F0_setup/analysis_log.md").exists(),
        "blocking",
        "logs/F0_setup/analysis_log.md",
    )
    add(
        "archive_readability",
        "GSE183904_RAW.tar readable and has 40 csv.gz members",
        f"{len(processed_manifest)} processed members",
        tar_path.exists() and len(processed_manifest) == 40,
        "blocking",
        "data/public_downloads/GSE183904_RAW.tar; data/metadata/processed_input_manifest.tsv",
    )
    add(
        "processed_input_manifest",
        "40 rows with file_size and sha256",
        f"{len(processed_manifest)} rows; empty_sha={sum(1 for r in processed_manifest if not r.get('sha256'))}",
        len(processed_manifest) == 40 and all(r.get("sha256") and r.get("file_size") for r in processed_manifest),
        "blocking",
        "data/metadata/processed_input_manifest.tsv",
    )
    add(
        "sample_info",
        "40 rows; no Unclear group",
        f"{len(sample_info)} rows; unclear={sum(1 for r in sample_info if r.get('group_analysis') == 'Unclear')}",
        len(sample_info) == 40 and all(r.get("group_analysis") != "Unclear" for r in sample_info),
        "blocking",
        "data/metadata/sample_info.tsv",
    )
    add(
        "data_audit",
        "40 rows; include_in_f1=true samples enter_full_F1 and no normalization artifact",
        f"{len(data_audit)} rows; pause={sum(1 for r in data_audit if r.get('audit_decision') != 'enter_full_F1')}",
        len(data_audit) == 40
        and all(
            r.get("audit_decision") == "enter_full_F1" and r.get("normalization_artifact_flag") == "false"
            for r in data_audit
            if r.get("include_in_f1") == "true"
        ),
        "blocking",
        "data/metadata/data_audit.tsv",
    )
    add(
        "precheck_comparison",
        "formal data_audit agrees with gse183904_csv_structure_precheck.tsv",
        f"{len(mismatches)} mismatching sample(s)",
        len(mismatches) == 0,
        "blocking",
        "data/metadata/data_audit.tsv; results/F0_audit/gse183904_csv_structure_precheck.tsv",
        "Must pause if any mismatch exists.",
    )
    add(
        "marker_panel",
        "marker panel exists and is not modified by F0",
        "exists" if (root / "data/metadata/cell_type_marker_panel.tsv").exists() else "missing",
        (root / "data/metadata/cell_type_marker_panel.tsv").exists(),
        "blocking",
        "data/metadata/cell_type_marker_panel.tsv",
    )
    for required_file in [
        "data/metadata/F0_dataset_inventory.tsv",
        "data/metadata/F0_file_manifest.tsv",
        "data/metadata/F0_metadata_field_inventory.tsv",
        "data/metadata/F0_author_processing_audit.tsv",
        "data/metadata/F0_data_readiness_by_F_section.tsv",
        "data/metadata/F0_external_resource_inventory.tsv",
        "data/metadata/F0_method_prior_decision.tsv",
        "data/metadata/decision_evidence_log.tsv",
        "data/metadata/excluded_samples.tsv",
    ]:
        add(
            Path(required_file).name,
            f"{required_file} exists",
            "exists" if (root / required_file).exists() else "missing",
            (root / required_file).exists(),
            "blocking",
            required_file,
        )
    return rows


def write_project_structure(root: Path) -> None:
    for directory in REQUIRED_DIRS:
        (root / directory).mkdir(parents=True, exist_ok=True)
    path = root / "data/metadata/project_structure_ready.txt"
    lines = [
        f"project_root={root.resolve()}",
        f"created_or_verified_at={now_iso()}",
        "operator=Codex",
        "directory_list=" + "|".join(REQUIRED_DIRS),
        "note=Directories were created if missing; existing inputs were not overwritten.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_reports(
    root: Path,
    sample_info: Sequence[Dict[str, object]],
    data_audit: Sequence[Dict[str, object]],
    gate_rows: Sequence[Dict[str, object]],
    mismatches: Sequence[str],
) -> None:
    group_counts: Dict[str, int] = {}
    for row in sample_info:
        group = str(row["group_analysis"])
        group_counts[group] = group_counts.get(group, 0) + 1
    n_enter = sum(1 for row in data_audit if row.get("audit_decision") == "enter_full_F1")
    blocking_failed = [row for row in gate_rows if row["blocking_level"] == "blocking" and row["pass_fail"] != "PASS"]
    sc_gate = "PASS" if not blocking_failed and n_enter == 40 else "FAIL"
    inventory_status = "partial_with_pending_local_inputs"

    recon = [
        "# F0 Global Data Reconnaissance Report",
        "",
        f"Run ID: {RUN_ID}",
        f"Generated at: {now_iso()}",
        "",
        "## Current Usable Data",
        "",
        f"- GSE183904: 40 author-filtered raw gene-count CSV.gz matrices; group counts: {group_counts}.",
        "- GSE183904 is the only dataset eligible to start F1 after F0 gate approval.",
        "- TCGA-STAD, bulk GEO/ACRG, GSE206785, GSE239676, GSE235046/SRP444325 and SCENIC resources are inventory candidates only until their section-specific audits approve use.",
        "",
        "## Boundaries",
        "",
        "- F0 does not produce biological conclusions.",
        "- GSE183904 does not provide FASTQ/raw droplets/empty droplets in current public inputs.",
        "- Normal_Peritoneum is reference/display only for group comparisons; Peritoneal_Metastasis sample-level inference is directional only because n=3.",
        "- If any formal stream audit value conflicts with the precheck table, F0 must pause before F1.",
        "",
        "## Method Decisions Deferred To Later Gates",
        "",
        "- F1 decides residual QC, doublet and ambient risk handling from F0 evidence plus Gate1 pilot.",
        "- F2 signature construction must not use external validation cohorts for tuning.",
        "- F5-F7 bulk/multiomics routes require expression-scale, clinical endpoint and sample-overlap audits before modeling.",
        "",
        "## F0 Gate Summary",
        "",
        f"- F0_scRNA_F1_gate: {sc_gate}",
        f"- F0_project_data_inventory_status: {inventory_status}",
        f"- Samples entering F1 object construction if approved: {n_enter}",
        f"- Blocking checklist failures: {len(blocking_failed)}",
    ]
    if mismatches:
        recon.extend(["", "## Precheck Mismatches", ""])
        recon.extend([f"- {note}" for note in mismatches])
    (root / "results/F0_audit/F0_global_data_reconnaissance_report.md").write_text(
        "\n".join(recon) + "\n", encoding="utf-8"
    )

    report = [
        "# F0 Execution Report",
        "",
        f"Run ID: {RUN_ID}",
        f"Generated at: {now_iso()}",
        "",
        "## Inputs Checked",
        "",
        "- data/public_downloads/GSE183904_RAW.tar",
        "- data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
        "- Existing manifests and F0 precheck tables under data/metadata and results/F0_audit",
        "",
        "## Main Observations",
        "",
        f"- GSE183904 sample_info rows: {len(sample_info)}",
        f"- GSE183904 data_audit rows: {len(data_audit)}",
        f"- Enter_full_F1 samples: {n_enter}",
        f"- Group counts: {group_counts}",
        "",
        "## Gate Decision",
        "",
        f"F0_scRNA_F1_gate: {sc_gate}",
        f"F0_project_data_inventory_status: {inventory_status}",
        "",
        "F1 may start only after Claude Code review and user approval.",
    ]
    (root / "results/F0_audit/F0_execution_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )


def dry_run(root: Path) -> int:
    missing = [item for item in REQUIRED_INPUTS if not (root / item).exists()]
    print("F0 dry run: no outputs will be written.")
    print(f"Project root: {root.resolve()}")
    print(f"Required inputs: {len(REQUIRED_INPUTS)}")
    if missing:
        print("Missing inputs:")
        for item in missing:
            print(f"  - {item}")
    else:
        print("All required inputs are present.")
    print("Planned formal outputs:")
    for item in F0_OUTPUTS:
        print(f"  - {item}")
    return 1 if missing else 0


def execute(root: Path) -> int:
    missing = [item for item in REQUIRED_INPUTS if not (root / item).exists()]
    if missing:
        raise RuntimeError("Missing required inputs: " + ", ".join(missing))

    write_project_structure(root)
    append_log(root, f"F0 execution started; run_id={RUN_ID}; python={sys.version.split()[0]}; os={platform.platform()}")

    tar_path = root / "data/public_downloads/GSE183904_RAW.tar"
    tar_sha = sha256_file(tar_path)
    download_rows = read_tsv(root / "data/metadata/download_manifest.tsv")
    manifest_tar = next((row for row in download_rows if row.get("file_name") == "GSE183904_RAW.tar"), {})
    expected_sha = manifest_tar.get("sha256", "").lower()
    if expected_sha and tar_sha.lower() != expected_sha:
        append_log(root, f"WARNING archive sha mismatch: observed={tar_sha}; manifest={expected_sha}")

    members = tar_csv_members(tar_path)
    member_names = [Path(member.name).name for member in members]
    if len(members) != 40:
        append_log(root, f"BLOCKING expected 40 csv.gz members but observed {len(members)}")
    extract_members(tar_path, members, root / "data/processed_input/GSE183904")

    sample_fields, patient_by_sample = parse_gse183904_series(
        root / "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz"
    )
    sample_info = build_sample_info(sample_fields, patient_by_sample, member_names)
    sample_fields_out = [
        "geo_accession",
        "sample_id",
        "sample_file",
        "patient_id",
        "sample_title",
        "source_name_ch1",
        "sample_characteristics_ch1",
        "tissue_site",
        "tumor_status",
        "group_analysis",
        "source_of_group",
        "metadata_confidence",
        "include_in_f1",
        "include_in_group_comparison",
        "is_paired",
        "pairing_scope",
        "note",
    ]
    write_tsv(root / "data/metadata/sample_info.tsv", list(sample_info), sample_fields_out)

    processed_manifest = build_processed_manifest(root, members, sample_info)
    processed_fields = [
        "source_archive",
        "archive_member_name",
        "extracted_path",
        "geo_accession",
        "sample_id",
        "file_size",
        "sha256",
        "extraction_date",
        "file_role",
        "note",
    ]
    write_tsv(root / "data/metadata/processed_input_manifest.tsv", list(processed_manifest), processed_fields)

    precheck_rows = read_tsv(root / "results/F0_audit/gse183904_csv_structure_precheck.tsv")
    data_audit, mismatches = build_data_audit(root, processed_manifest, sample_info, precheck_rows)
    data_audit_fields = [
        "geo_accession",
        "sample_id",
        "sample_file",
        "extracted_path",
        "matrix_orientation",
        "matrix_rows_genes",
        "matrix_cols_cells",
        "header_total_columns_including_gene_col",
        "first_header_fields",
        "first_genes",
        "barcode_suffix_examples",
        "bad_column_count_rows",
        "invalid_numeric_value_rows",
        "missing_value_rows",
        "mt_gene_count",
        "mt_gene_examples",
        "hb_gene_count",
        "hb_gene_examples",
        "gene_order_sha256",
        "suspected_matrix_type",
        "scan_seconds",
        "precheck_audit_decision",
        "audit_decision",
        "decision_scope",
        "include_in_f1",
        "raw_droplet_available",
        "empty_droplet_background_available",
        "normalization_artifact_flag",
        "precheck_comparison_status",
        "note",
    ]
    write_tsv(root / "data/metadata/data_audit.tsv", list(data_audit), data_audit_fields)
    for mismatch in mismatches:
        append_log(root, "PRECHECK_MISMATCH " + mismatch)

    precheck_files = {
        "non_GSE183904": read_tsv(root / "results/F0_audit/non_GSE183904_data_structure_precheck.tsv"),
        "predownloaded": read_tsv(root / "results/F0_audit/predownloaded_resource_structure_audit.tsv"),
        "bulk_GEO_series_matrix_deep_precheck": read_tsv(root / "results/F0_audit/bulk_GEO_series_matrix_deep_precheck.tsv"),
        "GSE206785_metadata_precheck": read_tsv(root / "results/F0_audit/GSE206785_metadata_precheck.tsv"),
    }

    dataset_inventory = build_dataset_inventory(sample_info, precheck_files)
    dataset_fields = [
        "dataset_id",
        "dataset_name_or_accession",
        "data_domain",
        "intended_F_sections",
        "current_availability_status",
        "primary_or_supplementary_role",
        "sample_or_cell_count_if_known",
        "patient_count_if_known",
        "event_count_if_known",
        "platform_or_assay",
        "gene_id_type_if_known",
        "raw_or_processed_status_if_known",
        "metadata_available",
        "clinical_endpoint_available",
        "major_limitations",
        "next_required_audit_step",
    ]
    write_tsv(root / "data/metadata/F0_dataset_inventory.tsv", dataset_inventory, dataset_fields)

    metadata_inventory = build_metadata_inventory(sample_fields, precheck_files)
    metadata_fields = [
        "dataset_id",
        "source_file",
        "field_name",
        "field_category",
        "n_records_checked",
        "unique_values_summary",
        "missing_count",
        "missing_rate",
        "example_values",
        "usable_for_F_section",
        "interpretation_risk",
        "note",
    ]
    write_tsv(root / "data/metadata/F0_metadata_field_inventory.tsv", metadata_inventory, metadata_fields)

    author_rows = build_author_processing_audit()
    author_fields = [
        "dataset_id",
        "source_reference_or_file",
        "processing_step",
        "author_reported_status",
        "method_or_threshold_if_reported",
        "evidence_location",
        "confidence_level",
        "implication_for_downstream_plan",
        "requires_special_handling",
    ]
    write_tsv(root / "data/metadata/F0_author_processing_audit.tsv", author_rows, author_fields)

    readiness_rows = build_data_readiness()
    readiness_fields = [
        "F_section",
        "required_data_domains",
        "available_datasets_or_files",
        "missing_or_pending_items",
        "minimum_data_needed_to_start",
        "specialized_audit_required_in_section",
        "current_readiness_status",
        "blocking_for_next_section",
        "recommended_first_action",
        "note",
    ]
    write_tsv(root / "data/metadata/F0_data_readiness_by_F_section.tsv", readiness_rows, readiness_fields)

    resource_rows = build_external_resource_inventory(root)
    resource_fields = [
        "resource_name",
        "resource_type",
        "planned_F_sections",
        "current_availability_status",
        "source_url_or_database",
        "version_or_release_if_known",
        "license_or_login_requirement",
        "local_manifest_if_any",
        "sha256_if_file",
        "large_or_unstable_download_risk",
        "next_required_audit_step",
        "note",
    ]
    write_tsv(root / "data/metadata/F0_external_resource_inventory.tsv", resource_rows, resource_fields)

    method_rows = build_method_prior_decision()
    method_fields = [
        "dataset_or_resource_id",
        "intended_F_sections",
        "known_data_structure",
        "known_limitation_from_precheck",
        "default_or_recommended_route",
        "method_not_allowed_by_default",
        "required_pre_execution_audit",
        "approval_required_before_change",
        "sensitivity_or_fallback_route",
        "interpretation_boundary",
        "note",
    ]
    write_tsv(root / "data/metadata/F0_method_prior_decision.tsv", method_rows, method_fields)

    decision_rows = build_decision_evidence_log()
    decision_fields = [
        "decision_id",
        "stage",
        "decision_topic",
        "decision_value",
        "decision_reason",
        "evidence_type",
        "evidence_source",
        "source_url_or_file",
        "evidence_strength",
        "date",
        "requires_sensitivity_analysis",
        "note",
    ]
    write_tsv(root / "data/metadata/decision_evidence_log.tsv", decision_rows, decision_fields)

    excluded_rows = build_excluded_samples(sample_info)
    excluded_fields = [
        "sample_id",
        "geo_accession",
        "exclusion_scope",
        "exclusion_status",
        "reason",
        "evidence_file",
        "note",
    ]
    write_tsv(root / "data/metadata/excluded_samples.tsv", excluded_rows, excluded_fields)

    gate_rows = build_gate_checklist(root, sample_info, data_audit, processed_manifest, mismatches)
    gate_fields = [
        "gate_item",
        "required_status",
        "observed_status",
        "pass_fail",
        "blocking_level",
        "evidence_file",
        "note",
    ]
    write_tsv(root / "results/F0_audit/F0_gate_checklist.tsv", gate_rows, gate_fields)
    write_reports(root, sample_info, data_audit, gate_rows, mismatches)

    file_manifest = build_file_manifest(root, processed_manifest, F0_OUTPUTS)
    file_fields = [
        "file_name",
        "relative_path_if_available",
        "dataset_id",
        "file_role",
        "data_domain",
        "source_type",
        "source_url_or_local_manifest",
        "file_size_bytes",
        "sha256",
        "compression_format",
        "read_status",
        "availability_status",
        "used_in_F0",
        "planned_F_section_use",
        "audit_status",
        "artifact_class",
        "publication_destination",
        "review_priority",
        "note",
    ]
    write_tsv(root / "data/metadata/F0_file_manifest.tsv", file_manifest, file_fields)

    append_log(root, f"F0 execution completed; run_id={RUN_ID}; gate_rows={len(gate_rows)}; mismatches={len(mismatches)}")
    return 0


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root. Defaults to current working directory.",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually write F0 outputs. Omit for dry-run input/output review.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run(root)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
