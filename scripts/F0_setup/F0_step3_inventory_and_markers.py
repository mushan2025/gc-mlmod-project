#!/usr/bin/env python3
"""F0 step 3: build inventories and audit the read-only marker panel.

Dependencies: successful step 1/2 outputs plus existing project manifests and
precheck tables. Outputs: dataset/file/metadata/author/resource inventories and
the conditional marker_panel_issue_report.tsv. The marker evidence table is
used only for reference-ID integrity, never to alter or expand the panel.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from F0_step2_sample_info_and_audit import parse_gse183904_series
from f0_utils import (
    F0_OUTPUTS,
    append_log,
    current_run_id,
    dry_run_report,
    normalize_sha256,
    now_iso,
    parse_stage_args,
    read_tsv,
    rel,
    require_paths,
    sha256_file,
    write_tsv,
)


STAGE_NAME = "F0 step 3 inventory and marker audit"
STAGE_REQUIRED = [
    "data/metadata/project_structure_ready.txt",
    "logs/F0_setup/analysis_log.md",
    "data/metadata/processed_input_manifest.tsv",
    "data/metadata/sample_info.tsv",
    "data/metadata/data_audit.tsv",
    "data/metadata/download_manifest.tsv",
    "data/metadata/preupload_resources_manifest.tsv",
    "data/metadata/preupload_pending_resources.tsv",
    "data/metadata/cell_type_marker_panel.tsv",
    "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
    "results/F0_audit/non_GSE183904_data_structure_precheck.tsv",
    "results/F0_audit/predownloaded_resource_structure_audit.tsv",
    "results/F0_audit/bulk_GEO_series_matrix_deep_precheck.tsv",
    "results/F0_audit/GSE206785_metadata_precheck.tsv",
]
STAGE_OUTPUTS = [
    "data/metadata/F0_dataset_inventory.tsv",
    "data/metadata/F0_file_manifest.tsv",
    "data/metadata/F0_metadata_field_inventory.tsv",
    "data/metadata/F0_author_processing_audit.tsv",
    "data/metadata/F0_external_resource_inventory.tsv",
    "results/F0_audit/marker_panel_issue_report.tsv (conditional)",
]

MARKER_GENE_PATTERN = re.compile(r"^[A-Z][A-Z0-9-]*$")
MARKER_REQUIRED_FIELDS = [
    "cell_type",
    "positive_markers",
    "minimum_rule",
    "confidence",
    "evidence_ids",
]


def summarize_unique(values: Iterable[str], max_items: int = 8) -> str:
    counts: Dict[str, int] = {}
    for value in values:
        key = value if value != "" else "<empty>"
        counts[key] = counts.get(key, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    shown = [f"{key}:{count}" for key, count in ordered[:max_items]]
    if len(ordered) > max_items:
        shown.append(f"...plus_{len(ordered) - max_items}_more")
    return "|".join(shown) if shown else "none"


def build_metadata_inventory(
    sample_fields: Dict[str, List[str]],
    prechecks: Dict[str, List[Dict[str, str]]],
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
                "interpretation_risk": (
                    "patient_pairing_requires_title_mapping_review"
                    if field in {"title", "characteristics_ch1"}
                    else "low"
                ),
                "note": "Current public matrix files do not provide cell-level metadata.",
            }
        )
    for row in prechecks.get("GSE206785_metadata_precheck", []):
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
                "interpretation_risk": (
                    "processed_external_scRNA_no_PM_group"
                    if row.get("column") in {"Group", "Sample"}
                    else "section_specific_review_required"
                ),
                "note": row.get("precheck_note", ""),
            }
        )
    for row in prechecks.get("bulk_GEO_series_matrix_deep_precheck", []):
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
    gse239_meta = next(
        (
            row
            for row in prechecks.get("predownloaded", [])
            if row.get("dataset_id") == "GSE239676" and row.get("file_name") == "GSE239676_meta.tsv.gz"
        ),
        None,
    )
    if gse239_meta:
        rows.append(
            {
                "dataset_id": "GSE239676",
                "source_file": gse239_meta.get("relative_path", ""),
                "field_name": "Sample,Patient,Tissue,HPinfec",
                "field_category": "cell_metadata_structure_preaudit",
                "n_records_checked": gse239_meta.get("n_rows", ""),
                "unique_values_summary": gse239_meta.get("value_type_or_key_fields", ""),
                "missing_count": "not_scanned_in_lightweight_preaudit",
                "missing_rate": "",
                "example_values": gse239_meta.get("metadata_or_sample_id_examples", ""),
                "usable_for_F_section": "F2.4_after_signature_freeze_and_cohort_approval",
                "interpretation_risk": "external_validation_isolation_required",
                "note": "Structure only; no expression values or biological outcomes were inspected.",
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
            "implication_for_downstream_plan": "F1 should verify author-filtered raw counts and perform conservative residual QC.",
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
            "implication_for_downstream_plan": "Do not claim FASTQ/BCL/Cell Ranger raw-output reprocessing; SoupX is not default without background droplets.",
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
            "implication_for_downstream_plan": "Public CSV.gz files are author-filtered raw gene-count matrices, not normalized expression.",
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
            "implication_for_downstream_plan": "F1 Gate1 must independently assess doublet risk.",
            "requires_special_handling": "Do not treat absent reporting as evidence that no doublets exist.",
        },
        {
            "dataset_id": "GSE183904",
            "source_reference_or_file": "GSE183904_series_matrix.txt.gz and available public files",
            "processing_step": "ambient_RNA_correction",
            "author_reported_status": "not_reported_in_public_metadata",
            "method_or_threshold_if_reported": "not_available",
            "evidence_location": "No public ambient-correction field found in current F0 inputs.",
            "confidence_level": "medium",
            "implication_for_downstream_plan": "F1 Gate1 should evaluate cross-expression and DecontX/celda feasibility.",
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
            "implication_for_downstream_plan": "Processed-expression external scoring only after approval; not raw QC reprocessing.",
            "requires_special_handling": "F2 external-validation audit required before use.",
        },
    ]


def gse239_inventory_row(predownloaded: Sequence[Dict[str, str]]) -> Dict[str, object] | None:
    rows = [row for row in predownloaded if row.get("dataset_id") == "GSE239676"]
    if not rows:
        return None
    matrix = next((row for row in rows if row.get("file_name") == "GSE239676_count_matrix.mtx.gz"), {})
    metadata = next((row for row in rows if row.get("file_name") == "GSE239676_meta.tsv.gz"), {})
    statuses = {row.get("audit_decision", "") for row in rows}
    ready = statuses == {"downloaded_structure_preaudit_pass"} and len(rows) >= 5
    return {
        "dataset_id": "GSE239676",
        "dataset_name_or_accession": "GSE239676",
        "data_domain": "scRNA_external_validation_candidate",
        "intended_F_sections": "F2.4,F8_external_validation",
        "current_availability_status": "project_local_structure_preaudit_pass" if ready else "pause_structure_mismatch",
        "primary_or_supplementary_role": "external_scRNA_validation_only",
        "sample_or_cell_count_if_known": f"{matrix.get('n_columns', '')} cells; GEO samples documented in series metadata",
        "patient_count_if_known": "20",
        "event_count_if_known": "not_applicable",
        "platform_or_assay": "single-cell RNA-seq; MatrixMarket integer coordinate matrix",
        "gene_id_type_if_known": "human_gene_symbol",
        "raw_or_processed_status_if_known": "sparse integer count matrix candidate; full F2.4 audit pending",
        "metadata_available": metadata.get("value_type_or_key_fields", ""),
        "clinical_endpoint_available": "not_for_survival_validation",
        "major_limitations": metadata.get("major_limitations", ""),
        "next_required_audit_step": "F2.4 external cohort audit after signature/scoring/high-low rules are frozen",
    }


def build_dataset_inventory(
    sample_info: Sequence[Dict[str, str]],
    prechecks: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, object]]:
    group_counts: Dict[str, int] = {}
    for row in sample_info:
        group = row.get("group_analysis", "")
        group_counts[group] = group_counts.get(group, 0) + 1
    group_summary = "|".join(f"{key}:{value}" for key, value in sorted(group_counts.items()))
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
            "platform_or_assay": "10x Genomics scRNA-seq; GPL24676; hg38",
            "gene_id_type_if_known": "gene_symbol_or_feature_name_from_author_CSV",
            "raw_or_processed_status_if_known": "author-filtered raw gene count CSV.gz; not FASTQ/raw droplets",
            "metadata_available": "GEO sample metadata; no public cell-level metadata in current inputs",
            "clinical_endpoint_available": "not_for_survival",
            "major_limitations": "PM n=3; Normal_Peritoneum n=1; raw/empty droplets unavailable",
            "next_required_audit_step": "F0 gate review, then approved F1 Gate1 QC/doublet/ambient assessment",
        }
    ]
    for row in prechecks.get("non_GSE183904", []):
        rows.append(
            {
                "dataset_id": row.get("dataset_or_resource_id", ""),
                "dataset_name_or_accession": row.get("dataset_or_resource_id", ""),
                "data_domain": row.get("data_domain", ""),
                "intended_F_sections": row.get("intended_F_sections", ""),
                "current_availability_status": "project_local_available" if row.get("local_relative_path") else "unknown",
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
    gse239 = gse239_inventory_row(prechecks.get("predownloaded", []))
    if gse239:
        rows.append(gse239)
    for row in prechecks.get("predownloaded", []):
        if row.get("dataset_id") == "GSE239676":
            continue
        rows.append(
            {
                "dataset_id": row.get("dataset_id", ""),
                "dataset_name_or_accession": row.get("dataset_id", ""),
                "data_domain": "preloaded_resource",
                "intended_F_sections": "section_specific",
                "current_availability_status": (
                    "project_local_available" if row.get("download_status") == "complete" else row.get("download_status", "unknown")
                ),
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
    processed_manifest: Sequence[Dict[str, str]],
    generated_paths: Sequence[str],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for source_path, source_name in [
        ("data/metadata/download_manifest.tsv", "download_manifest"),
        ("data/metadata/preupload_resources_manifest.tsv", "preupload_resources"),
    ]:
        for row in read_tsv(root / source_path):
            relative = row.get("relative_path_under_data_public_downloads", "")
            path = (
                root / "data/public_downloads" / relative
                if relative
                else root / "data/public_downloads" / row.get("file_name", "")
            )
            rows.append(
                {
                    "file_name": row.get("file_name", ""),
                    "relative_path_if_available": rel(path, root) if path.exists() else relative,
                    "dataset_id": row.get("dataset_id", ""),
                    "file_role": row.get("file_role", ""),
                    "data_domain": "mixed",
                    "source_type": "researcher_provided_local" if source_name == "download_manifest" else "locally_precached",
                    "source_url_or_local_manifest": row.get("source_url", row.get("source_url_or_derivation", "")),
                    "file_size_bytes": row.get("file_size", row.get("file_size_bytes", "")),
                    "sha256": normalize_sha256(row.get("sha256", "")),
                    "compression_format": Path(row.get("file_name", "")).suffix.lstrip("."),
                    "read_status": row.get("read_status", "not_checked_in_F0_file_manifest"),
                    "availability_status": "project_local_available" if path.exists() else "not_found_at_expected_path",
                    "used_in_F0": "true" if row.get("dataset_id") == "GSE183904" else "inventory_only",
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
                "source_type": "generated_by_F0_step1",
                "source_url_or_local_manifest": row.get("source_archive", ""),
                "file_size_bytes": row.get("file_size", ""),
                "sha256": normalize_sha256(row.get("sha256", "")),
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
    for output in generated_paths:
        path = root / output
        is_self = Path(output).name == "F0_file_manifest.tsv"
        rows.append(
            {
                "file_name": Path(output).name,
                "relative_path_if_available": output,
                "dataset_id": "F0",
                "file_role": "F0_generated_output",
                "data_domain": "audit",
                "source_type": "generated_by_staged_F0_scripts",
                "source_url_or_local_manifest": "scripts/F0_setup/run_F0_full_audit.py",
                "file_size_bytes": path.stat().st_size if path.exists() else "",
                "sha256": sha256_file(path) if path.exists() and not is_self else "",
                "compression_format": Path(output).suffix.lstrip("."),
                "read_status": "written" if path.exists() else "not_yet_written",
                "availability_status": "project_local_available" if path.exists() else "planned",
                "used_in_F0": "true",
                "planned_F_section_use": "F0,F1-F8_as_relevant",
                "audit_status": "generated_by_current_staged_scripts" if path.exists() else "planned",
                "artifact_class": "gate_decision" if Path(output).name in {"F0_gate_checklist.tsv", "F0_execution_report.md"} else "audit_trail",
                "publication_destination": "undecided" if Path(output).name in {"F0_global_data_reconnaissance_report.md", "F0_execution_report.md"} else "not_planned",
                "review_priority": "standard",
                "note": "Self-manifest SHA256 intentionally blank" if is_self else "",
            }
        )
    return rows


def build_external_resource_inventory(root: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for row in read_tsv(root / "data/metadata/preupload_resources_manifest.tsv"):
        dataset = row.get("dataset_id", "")
        relative = row.get("relative_path_under_data_public_downloads", "")
        if dataset.startswith("SCENIC") or dataset == "inferCNV_gene_order" or "SCENIC" in relative:
            size = int(row.get("file_size_bytes", "0") or "0")
            rows.append(
                {
                    "resource_name": row.get("file_name", ""),
                    "resource_type": dataset,
                    "planned_F_sections": "F3" if dataset.startswith("SCENIC") else "F1",
                    "current_availability_status": "project_local_available" if row.get("exists") == "TRUE" else "not_available",
                    "source_url_or_database": row.get("source_url_or_derivation", ""),
                    "version_or_release_if_known": "hg38_v10" if "SCENIC" in dataset else "",
                    "license_or_login_requirement": "public_resource_verify_terms_before_publication",
                    "local_manifest_if_any": "data/metadata/preupload_resources_manifest.tsv",
                    "sha256_if_file": normalize_sha256(row.get("sha256", "")),
                    "large_or_unstable_download_risk": "yes" if size > 100_000_000 else "no",
                    "next_required_audit_step": "section_specific_compatibility_check",
                    "note": row.get("note", ""),
                }
            )
    for name, resource_type, sections, note in [
        ("MSigDB_GOBP_FERROPTOSIS", "external_database", "F2,F3", "Version and license audit before enrichment/scoring."),
        ("Reactome", "external_database", "F3,F7", "Version audit before pathway interpretation."),
        ("KEGG", "external_database", "F3,F7", "Access and license constraints must be checked before use."),
        ("JASPAR/DoRothEA", "external_database", "F3,F8", "Use only after GRN method-plan approval."),
        ("TIDE/IPS/SubMap", "web_or_API_resource", "F6", "Prepare only after final gene/signature inputs are frozen."),
    ]:
        rows.append(
            {
                "resource_name": name,
                "resource_type": resource_type,
                "planned_F_sections": sections,
                "current_availability_status": "not_yet_audited",
                "source_url_or_database": name,
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


def audit_marker_panel(root: Path) -> List[Dict[str, object]]:
    panel_path = root / "data/metadata/cell_type_marker_panel.tsv"
    evidence_path = root / "data/metadata/cell_marker_reference_evidence.tsv"
    panel = read_tsv(panel_path)
    evidence = read_tsv(evidence_path) if evidence_path.exists() else []
    panel_fields = set(panel[0].keys()) if panel else set()
    evidence_ids = {row.get("evidence_id", "") for row in evidence}
    issues: List[Dict[str, object]] = []

    def add(cell_type: str, issue_type: str, field: str, observed: str, action: str) -> None:
        issues.append(
            {
                "issue_id": f"MARKER_ISSUE_{len(issues) + 1:03d}",
                "cell_type": cell_type or "table_level",
                "issue_type": issue_type,
                "field_name": field,
                "observed_value": observed,
                "severity": "warning",
                "blocking_for_F0": "false",
                "recommended_action": action,
                "evidence_source": "cell_type_marker_panel.tsv; cell_marker_reference_evidence.tsv (integrity-only)",
            }
        )

    for field in MARKER_REQUIRED_FIELDS:
        if field not in panel_fields:
            add("", "missing_required_column", field, "missing", "Researcher review required; do not modify panel automatically.")
    if not evidence_path.exists():
        add(
            "",
            "missing_internal_evidence_reference",
            "evidence_ids",
            rel(evidence_path, root),
            "Restore the read-only internal reference before F1 annotation review; do not alter the panel.",
        )

    seen: set[str] = set()
    for row in panel:
        cell_type = row.get("cell_type", "")
        if cell_type in seen:
            add(cell_type, "duplicate_cell_type", "cell_type", cell_type, "Resolve duplicate definition after researcher approval.")
        seen.add(cell_type)
        positive = row.get("positive_markers", "").strip()
        if positive in {"", "NA"}:
            add(cell_type, "empty_positive_markers", "positive_markers", positive, "Define positive markers only after researcher approval.")
        else:
            for gene in [item.strip() for item in positive.split(",") if item.strip()]:
                if not MARKER_GENE_PATTERN.fullmatch(gene):
                    add(cell_type, "invalid_marker_format", "positive_markers", gene, "Verify HGNC-style symbol and approve any correction.")
        for evidence_id in [item.strip() for item in row.get("evidence_ids", "").split(",") if item.strip() and item.strip() != "NA"]:
            if evidence_id not in evidence_ids:
                add(cell_type, "missing_evidence_reference", "evidence_ids", evidence_id, "Verify evidence ID; do not alter panel automatically.")
    return issues


def execute(root: Path) -> int:
    require_paths(root, STAGE_REQUIRED, STAGE_NAME)
    append_log(root, f"F0 step3 started; run_id={current_run_id()}")
    sample_info = read_tsv(root / "data/metadata/sample_info.tsv")
    processed_manifest = read_tsv(root / "data/metadata/processed_input_manifest.tsv")
    data_audit = read_tsv(root / "data/metadata/data_audit.tsv")
    if len(processed_manifest) != 40 or len(sample_info) != 40 or len(data_audit) != 40:
        raise RuntimeError(
            "F0 step3 requires 40 processed-manifest, sample-info and data-audit rows"
        )
    if any(row.get("audit_decision") != "enter_full_F1" for row in data_audit):
        raise RuntimeError("F0 step3 cannot continue while data_audit.tsv contains paused samples")
    sample_fields, _ = parse_gse183904_series(
        root / "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz"
    )
    prechecks = {
        "non_GSE183904": read_tsv(root / "results/F0_audit/non_GSE183904_data_structure_precheck.tsv"),
        "predownloaded": read_tsv(root / "results/F0_audit/predownloaded_resource_structure_audit.tsv"),
        "bulk_GEO_series_matrix_deep_precheck": read_tsv(root / "results/F0_audit/bulk_GEO_series_matrix_deep_precheck.tsv"),
        "GSE206785_metadata_precheck": read_tsv(root / "results/F0_audit/GSE206785_metadata_precheck.tsv"),
    }

    dataset_rows = build_dataset_inventory(sample_info, prechecks)
    write_tsv(
        root / "data/metadata/F0_dataset_inventory.tsv",
        dataset_rows,
        [
            "dataset_id", "dataset_name_or_accession", "data_domain", "intended_F_sections",
            "current_availability_status", "primary_or_supplementary_role", "sample_or_cell_count_if_known",
            "patient_count_if_known", "event_count_if_known", "platform_or_assay", "gene_id_type_if_known",
            "raw_or_processed_status_if_known", "metadata_available", "clinical_endpoint_available",
            "major_limitations", "next_required_audit_step",
        ],
    )
    metadata_rows = build_metadata_inventory(sample_fields, prechecks)
    write_tsv(
        root / "data/metadata/F0_metadata_field_inventory.tsv",
        metadata_rows,
        [
            "dataset_id", "source_file", "field_name", "field_category", "n_records_checked",
            "unique_values_summary", "missing_count", "missing_rate", "example_values",
            "usable_for_F_section", "interpretation_risk", "note",
        ],
    )
    author_rows = build_author_processing_audit()
    write_tsv(
        root / "data/metadata/F0_author_processing_audit.tsv",
        author_rows,
        [
            "dataset_id", "source_reference_or_file", "processing_step", "author_reported_status",
            "method_or_threshold_if_reported", "evidence_location", "confidence_level",
            "implication_for_downstream_plan", "requires_special_handling",
        ],
    )
    resource_rows = build_external_resource_inventory(root)
    write_tsv(
        root / "data/metadata/F0_external_resource_inventory.tsv",
        resource_rows,
        [
            "resource_name", "resource_type", "planned_F_sections", "current_availability_status",
            "source_url_or_database", "version_or_release_if_known", "license_or_login_requirement",
            "local_manifest_if_any", "sha256_if_file", "large_or_unstable_download_risk",
            "next_required_audit_step", "note",
        ],
    )

    issues = audit_marker_panel(root)
    issue_path = root / "results/F0_audit/marker_panel_issue_report.tsv"
    if issues:
        write_tsv(
            issue_path,
            issues,
            [
                "issue_id", "cell_type", "issue_type", "field_name", "observed_value",
                "severity", "blocking_for_F0", "recommended_action", "evidence_source",
            ],
        )
    elif issue_path.exists():
        issue_path.unlink()

    file_rows = build_file_manifest(root, processed_manifest, F0_OUTPUTS)
    write_tsv(
        root / "data/metadata/F0_file_manifest.tsv",
        file_rows,
        [
            "file_name", "relative_path_if_available", "dataset_id", "file_role", "data_domain",
            "source_type", "source_url_or_local_manifest", "file_size_bytes", "sha256",
            "compression_format", "read_status", "availability_status", "used_in_F0",
            "planned_F_section_use", "audit_status", "artifact_class", "publication_destination",
            "review_priority", "note",
        ],
    )
    marker_status = "no_issues" if not issues else f"{len(issues)}_nonblocking_issue(s)"
    append_log(
        root,
        f"F0 step3 completed; dataset_rows={len(dataset_rows)}; metadata_rows={len(metadata_rows)}; "
        f"resource_rows={len(resource_rows)}; marker_panel={marker_status}; F1 reminder: add marker-panel and evidence-reference method priors.",
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, STAGE_REQUIRED, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
