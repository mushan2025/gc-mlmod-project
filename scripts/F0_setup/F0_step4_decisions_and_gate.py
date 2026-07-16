#!/usr/bin/env python3
"""F0 step 4: write method decisions, gate checklist, and final reports.

Dependencies: all successful step 1-3 outputs. Outputs: readiness/method/
decision/exclusion tables, the final gate checklist and reports, and a refreshed
F0 file manifest. This step finalizes F0 only; it never starts F1.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence

from F0_step3_inventory_and_markers import build_file_manifest
from f0_utils import (
    F0_OUTPUTS,
    append_log,
    current_run_id,
    dry_run_report,
    now_iso,
    parse_stage_args,
    read_tsv,
    require_paths,
    write_tsv,
)


STAGE_NAME = "F0 step 4 decisions and gate"
STAGE_REQUIRED = [
    "data/metadata/project_structure_ready.txt",
    "logs/F0_setup/analysis_log.md",
    "data/metadata/processed_input_manifest.tsv",
    "data/metadata/sample_info.tsv",
    "data/metadata/data_audit.tsv",
    "data/metadata/F0_dataset_inventory.tsv",
    "data/metadata/F0_file_manifest.tsv",
    "data/metadata/F0_metadata_field_inventory.tsv",
    "data/metadata/F0_author_processing_audit.tsv",
    "data/metadata/F0_external_resource_inventory.tsv",
]
STAGE_OUTPUTS = [
    "data/metadata/F0_data_readiness_by_F_section.tsv",
    "data/metadata/F0_method_prior_decision.tsv",
    "data/metadata/decision_evidence_log.tsv",
    "data/metadata/excluded_samples.tsv",
    "results/F0_audit/F0_gate_checklist.tsv",
    "results/F0_audit/F0_global_data_reconnaissance_report.md",
    "results/F0_audit/F0_execution_report.md",
]


def build_data_readiness() -> List[Dict[str, object]]:
    return [
        {
            "F_section": "F1",
            "required_data_domains": "GSE183904 scRNA count matrices; sample metadata; marker panel",
            "available_datasets_or_files": "GSE183904_RAW.tar; GSE183904_series_matrix.txt.gz; cell_type_marker_panel.tsv",
            "missing_or_pending_items": "F1 R dependencies remain required for selected doublet/ambient/CNV methods",
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
            "available_datasets_or_files": "GSE235046 metadata/count-like table; SRP444325 manifest; GSE206785 candidate; GSE239676 structure preaudit passed; TCGA/bulk GEO candidates",
            "missing_or_pending_items": "Signature must be frozen before validation; SRP444325 raw reprocessing needs its approved plan; external cohorts need F2.4 audit",
            "minimum_data_needed_to_start": "F1 approved malignant/candidate-cell object plus F2.1 plan",
            "specialized_audit_required_in_section": "SRA reprocessing or approved fallback; bulk/external sc cohort audit",
            "current_readiness_status": "ready_with_limitations",
            "blocking_for_next_section": "no_for_F1; yes_for_F2_start",
            "recommended_first_action": "After F1, audit signature-source route and preserve external-validation isolation",
            "note": "Do not round fractional GSE235046 table for primary DESeq2; do not tune from GSE239676.",
        },
        {
            "F_section": "F3",
            "required_data_domains": "F2 frozen states; F1 malignant epithelial object; SCENIC resources",
            "available_datasets_or_files": "SCENIC resource manifest and candidate local files",
            "missing_or_pending_items": "SCENIC compatibility and gene-overlap audit",
            "minimum_data_needed_to_start": "F2 candidate/high-low definitions approved",
            "specialized_audit_required_in_section": "pySCENIC/ctxcore compatibility and resource-version audit",
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
            "minimum_data_needed_to_start": "F2 state and F1 cell-type labels approved",
            "specialized_audit_required_in_section": "LR database/version and sample-aware expression-support audit",
            "current_readiness_status": "unknown",
            "blocking_for_next_section": "no_for_F1; yes_for_F4_start",
            "recommended_first_action": "Prepare F4 LR method plan after F2 gate approval",
            "note": "CellChat network is background; sample-aware LR expression is primary.",
        },
        {
            "F_section": "F5",
            "required_data_domains": "Bulk expression and clinical/survival endpoints",
            "available_datasets_or_files": "TCGA-STAD Xena/GDC/cBioPortal candidates; GEO/ACRG prechecks",
            "missing_or_pending_items": "Dataset-specific clinical endpoint and expression-scale audit",
            "minimum_data_needed_to_start": "F2 clinical/signature evidence state and approved bulk-cohort audit",
            "specialized_audit_required_in_section": "Probe mapping, scale, endpoint and overlap audit",
            "current_readiness_status": "ready_with_limitations",
            "blocking_for_next_section": "no_for_F1; yes_for_F5_start",
            "recommended_first_action": "Freeze F5 cohort ranking and endpoint availability before modeling",
            "note": "SuperSeries overlap must be handled; processed arrays are not count-model inputs.",
        },
        {
            "F_section": "F6",
            "required_data_domains": "Bulk expression; F1 reference; F5 genes as applicable",
            "available_datasets_or_files": "TCGA/bulk candidates and GSE183904 reference candidate",
            "missing_or_pending_items": "BayesPrism benchmark plan; TIDE/IPS input/output approvals",
            "minimum_data_needed_to_start": "F5 core genes and approved deconvolution method plan",
            "specialized_audit_required_in_section": "Simulation/LOPO benchmark and immune-method legality/version audit",
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
            "note": "Reverse validation permits VST/WGCNA after F7 audit, not automatic DE.",
        },
        {
            "F_section": "F8",
            "required_data_domains": "F5 selected gene(s); F1/F2 objects; optional perturbation resources",
            "available_datasets_or_files": "Depends on F5/F2 outputs",
            "missing_or_pending_items": "No final gene selected yet",
            "minimum_data_needed_to_start": "F5 approved final-gene list and F2/F3/F4 context",
            "specialized_audit_required_in_section": "Single-gene evidence and perturbation-model reliability audit",
            "current_readiness_status": "not_ready",
            "blocking_for_next_section": "no_for_F1; yes_for_F8_start",
            "recommended_first_action": "Wait for upstream outputs",
            "note": "F8.3 evidence ceiling remains model_supported_hypothesis.",
        },
    ]


def build_method_prior_decision() -> List[Dict[str, object]]:
    return [
        {
            "dataset_or_resource_id": "GSE183904",
            "intended_F_sections": "F1-F4,F8",
            "known_data_structure": "author-filtered raw gene-count CSV.gz; gene rows by cell columns",
            "known_limitation_from_precheck": "No FASTQ/raw/empty droplets; PM n=3; Normal_Peritoneum n=1",
            "default_or_recommended_route": "F1 conservative QC reanalysis; verify QC, doublet and ambient risk",
            "method_not_allowed_by_default": "No FASTQ/Cell Ranger raw claim; no forced SoupX without background; no MLMOD in F1",
            "required_pre_execution_audit": "F0 data_audit and F1 Gate1 QC/doublet/ambient plan",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "If count structure conflicts, pause and revise F1 input strategy",
            "interpretation_boundary": "F1 produces a reusable annotated object, not MLMOD conclusions",
            "note": "Normal_Peritoneum reference-only; PM sample-level statistics directional only.",
        },
        {
            "dataset_or_resource_id": "GSE235046/SRP444325",
            "intended_F_sections": "F2.1",
            "known_data_structure": "Mouse BMDM RNA-seq; public table has decimals; ENA paired FASTQ route prepared",
            "known_limitation_from_precheck": "GEO count-like table is not all integer raw counts",
            "default_or_recommended_route": "Approved SRA reprocessing STAR/RSEM/tximport/DESeq2 or equivalent auditable count model",
            "method_not_allowed_by_default": "Do not round decimal table for primary DESeq2",
            "required_pre_execution_audit": "F2.1 raw-reprocessing resource and design audit",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "TMM plus limma-voom/limma conditional fallback",
            "interpretation_boundary": "Cross-species macrophage-to-epithelial transfer remains a hypothesis requiring controls",
            "note": "Interaction/Torin annotations do not filter signature membership.",
        },
        {
            "dataset_or_resource_id": "TCGA-STAD_Xena_GDC_cBioPortal",
            "intended_F_sections": "F2.4,F5,F7",
            "known_data_structure": "Xena star_counts log2(count+1) but inverse-validatable; TPM/clinical/CNV/MAF/450K candidates exist",
            "known_limitation_from_precheck": "Stored values are not directly raw integers; multiomics needs sample intersection/QC",
            "default_or_recommended_route": "Processed expression for scoring/Cox; reconstructed counts for VST/WGCNA only after F7 audit",
            "method_not_allowed_by_default": "Do not use stored log2 values as raw-count DE input",
            "required_pre_execution_audit": "F2.4/F5/F7 expression, clinical, mutation, CNV and methylation audits",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Download direct GDC raw_count only if later approved DE requires it",
            "interpretation_boundary": "Bulk association supports prognosis association, not single-cell causal mechanism",
            "note": "TCGA_STAD_star_counts_inverse_validation.tsv must be rechecked in F7.",
        },
        {
            "dataset_or_resource_id": "bulk_GEO_ACRG_candidates",
            "intended_F_sections": "F2.4,F5",
            "known_data_structure": "Processed microarray series matrices with heterogeneous scale and overlap risk",
            "known_limitation_from_precheck": "Linear/log-like scales differ; SuperSeries sample overlap exists",
            "default_or_recommended_route": "Probe-to-gene audit and within-cohort scoring/limma/Cox/KM as appropriate",
            "method_not_allowed_by_default": "No count models; no combining before scale/probe/overlap audit",
            "required_pre_execution_audit": "Scale, annotation, clinical endpoint and overlap audit",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Use only cohorts with sufficient endpoints and non-overlap",
            "interpretation_boundary": "Validation strength depends on independent patients and endpoint completeness",
            "note": "GSE62254 subset of GSE66229; GSE84426 subset/same-release risk with GSE84437.",
        },
        {
            "dataset_or_resource_id": "GSE206785",
            "intended_F_sections": "F2_external_sc_validation",
            "known_data_structure": "Cell-by-gene log1p(count)-like processed matrix with metadata",
            "known_limitation_from_precheck": "No PM group; not raw counts",
            "default_or_recommended_route": "Processed-expression scoring/localization after F2 freeze",
            "method_not_allowed_by_default": "No integer-count pseudobulk or PM validation",
            "required_pre_execution_audit": "External cohort audit after F2 rules freeze",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "not_evaluable if coverage or mapping fails",
            "interpretation_boundary": "Cross-dataset scoring transfer only, not raw-processing reproducibility",
            "note": "Must not tune F2 signature using GSE206785.",
        },
        {
            "dataset_or_resource_id": "SCENIC_hg38_v10",
            "intended_F_sections": "F3.3",
            "known_data_structure": "TF list, motif-to-TF table, rankings.feather resources",
            "known_limitation_from_precheck": "pyarrow/ctxcore/pySCENIC compatibility not yet audited",
            "default_or_recommended_route": "Technical compatibility smoke test before full GRNBoost2",
            "method_not_allowed_by_default": "Do not use smoke test as biological filter or MLMOD-balanced downsample",
            "required_pre_execution_audit": "Version, overlap and multi-seed stability plan",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "Approved server if full SCENIC exceeds local capacity",
            "interpretation_boundary": "GRN perturbation supports model-supported hypotheses only",
            "note": "Formal GRNBoost2 requires at least 10 seeds if run.",
        },
    ]


def build_decision_evidence_log() -> List[Dict[str, object]]:
    date = now_iso()
    return [
        {
            "decision_id": "F0_DECISION_001",
            "stage": "F0",
            "decision_topic": "primary_scRNA_dataset",
            "decision_value": "GSE183904 is the primary discovery scRNA dataset",
            "decision_reason": "The project hypothesis and plan use 40 local GSE183904 sample matrices.",
            "evidence_type": "local_file_and_GEO_metadata",
            "evidence_source": "GSE183904_RAW.tar; GSE183904_series_matrix.txt.gz",
            "source_url_or_file": "data/public_downloads/GSE183904_RAW.tar",
            "evidence_strength": "high_for_file_availability",
            "date": date,
            "requires_sensitivity_analysis": "no",
            "note": "Scientific conclusions require downstream gates.",
        },
        {
            "decision_id": "F0_DECISION_002",
            "stage": "F0/F1",
            "decision_topic": "GSE183904_matrix_boundary",
            "decision_value": "Use author-filtered raw gene-count matrices, not FASTQ/raw droplets",
            "decision_reason": "GEO reports raw files not submitted; supplementary matrices are post-Cell-Ranger counts after author filtering.",
            "evidence_type": "GEO_metadata_and_stream_audit",
            "evidence_source": "GSE183904_series_matrix.txt.gz; data_audit.tsv",
            "source_url_or_file": "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
            "evidence_strength": "high",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "F1 must verify residual QC/doublet/ambient risk.",
        },
        {
            "decision_id": "F0_DECISION_003",
            "stage": "F0/F1",
            "decision_topic": "group_handling",
            "decision_value": "Normal_Peritoneum is reference-only; PM n=3 is directional only",
            "decision_reason": "GEO sample titles show one normal peritoneum and three peritoneal-tumor samples.",
            "evidence_type": "GEO_sample_metadata",
            "evidence_source": "sample_info.tsv",
            "source_url_or_file": "data/metadata/sample_info.tsv",
            "evidence_strength": "high_for_group_labels",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "PM cannot be treated as an independent validation source.",
        },
        {
            "decision_id": "F0_DECISION_004",
            "stage": "F0/F1",
            "decision_topic": "marker_panel_lock",
            "decision_value": "cell_type_marker_panel.tsv is read-only preregistered annotation input",
            "decision_reason": "Project rules require reporting issues without modifying the panel.",
            "evidence_type": "project_contract",
            "evidence_source": "AGENTS.md and main plan",
            "source_url_or_file": "data/metadata/cell_type_marker_panel.tsv",
            "evidence_strength": "project_rule",
            "date": date,
            "requires_sensitivity_analysis": "no",
            "note": "F0 validates structure and evidence-ID integrity only.",
        },
    ]


def build_excluded_samples(sample_info: Sequence[Dict[str, str]]) -> List[Dict[str, object]]:
    excluded = [row for row in sample_info if row.get("include_in_f1") != "true"]
    if not excluded:
        return [
            {
                "sample_id": "none",
                "geo_accession": "none",
                "exclusion_scope": "F1_object_construction",
                "exclusion_status": "no_exclusion",
                "reason": "All 40 samples are eligible if data_audit.tsv remains PASS.",
                "evidence_file": "data/metadata/sample_info.tsv; data/metadata/data_audit.tsv",
                "note": "Normal_Peritoneum is excluded from main comparison, not object construction by default.",
            }
        ]
    return [
        {
            "sample_id": row.get("sample_id", ""),
            "geo_accession": row.get("geo_accession", ""),
            "exclusion_scope": "F1_object_construction",
            "exclusion_status": "pending_or_excluded",
            "reason": "F0 group or metadata unresolved",
            "evidence_file": "data/metadata/sample_info.tsv",
            "note": row.get("note", ""),
        }
        for row in excluded
    ]


def build_gate_checklist(
    root: Path,
    sample_info: Sequence[Dict[str, str]],
    data_audit: Sequence[Dict[str, str]],
    processed_manifest: Sequence[Dict[str, str]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []

    def add(
        item: str,
        required: str,
        observed: str,
        status: str,
        level: str,
        evidence: str,
        note: str = "",
    ) -> None:
        rows.append(
            {
                "gate_item": item,
                "required_status": required,
                "observed_status": observed,
                "pass_fail": status,
                "blocking_level": level,
                "evidence_file": evidence,
                "note": note,
            }
        )

    def pass_fail(ok: bool) -> str:
        return "PASS" if ok else "FAIL"

    add(
        "project_structure_ready",
        "project_structure_ready.txt exists",
        "exists" if (root / "data/metadata/project_structure_ready.txt").exists() else "missing",
        pass_fail((root / "data/metadata/project_structure_ready.txt").exists()),
        "blocking",
        "data/metadata/project_structure_ready.txt",
    )
    add(
        "analysis_log",
        "analysis_log.md exists",
        "exists" if (root / "logs/F0_setup/analysis_log.md").exists() else "missing",
        pass_fail((root / "logs/F0_setup/analysis_log.md").exists()),
        "blocking",
        "logs/F0_setup/analysis_log.md",
    )
    add(
        "archive_readability",
        "40 processed csv.gz members",
        f"{len(processed_manifest)} processed members",
        pass_fail(len(processed_manifest) == 40),
        "blocking",
        "data/public_downloads/GSE183904_RAW.tar; data/metadata/processed_input_manifest.tsv",
    )
    manifest_ok = len(processed_manifest) == 40 and all(row.get("sha256") and row.get("file_size") for row in processed_manifest)
    add(
        "processed_input_manifest",
        "40 rows with file_size and SHA256",
        f"{len(processed_manifest)} rows; empty_sha={sum(1 for row in processed_manifest if not row.get('sha256'))}",
        pass_fail(manifest_ok),
        "blocking",
        "data/metadata/processed_input_manifest.tsv",
    )
    sample_ok = len(sample_info) == 40 and all(
        row.get("group_analysis") != "Unclear"
        and row.get("sample_id_match_status") == "match"
        and row.get("source_of_group")
        and row.get("metadata_confidence")
        for row in sample_info
    )
    add(
        "sample_info",
        "40 rows; no Unclear group; manifest/GEO-title sample IDs match",
        f"{len(sample_info)} rows; unclear={sum(1 for row in sample_info if row.get('group_analysis') == 'Unclear')}; "
        f"sample_id_mismatch={sum(1 for row in sample_info if row.get('sample_id_match_status') != 'match')}",
        pass_fail(sample_ok),
        "blocking",
        "data/metadata/sample_info.tsv",
    )
    audit_ok = len(data_audit) == 40 and all(
        row.get("audit_decision") == "enter_full_F1"
        and row.get("normalization_artifact_flag") == "false"
        and row.get("observed_numeric_type") == "nonnegative_integer_count_like"
        and row.get("suspected_matrix_type") == "author_filtered_raw_gene_count_matrix"
        and row.get("processed_input_manifest_match") == "true"
        and row.get("per_gene_mean_distribution_status") == "consistent_with_sparse_right_skew"
        and row.get("format_decision_scope") == "file_format_only"
        and row.get("decision_scope") == "file_format_and_public_processing_boundary"
        for row in data_audit
        if row.get("include_in_f1") == "true"
    )
    add(
        "data_audit",
        "40 rows; included samples pass numeric/distribution/mapping checks and both decision layers",
        f"{len(data_audit)} rows; pause={sum(1 for row in data_audit if row.get('audit_decision') != 'enter_full_F1')}",
        pass_fail(audit_ok),
        "blocking",
        "data/metadata/data_audit.tsv",
    )
    mismatches = [row for row in data_audit if row.get("precheck_comparison_status") != "match"]
    add(
        "precheck_comparison",
        "formal data_audit agrees with preregistered key fields",
        f"{len(mismatches)} mismatching sample(s)",
        pass_fail(not mismatches),
        "blocking",
        "data/metadata/data_audit.tsv; results/F0_audit/gse183904_csv_structure_precheck.tsv",
        "Must pause on any mismatch.",
    )
    issue_path = root / "results/F0_audit/marker_panel_issue_report.tsv"
    issue_count = len(read_tsv(issue_path)) if issue_path.exists() else 0
    marker_status = "PASS_WITH_NOTED_ISSUES" if issue_count else "PASS"
    add(
        "marker_panel",
        "panel exists; content/evidence references audited; issues only reported",
        f"exists; issues={issue_count}",
        marker_status,
        "nonblocking",
        "data/metadata/cell_type_marker_panel.tsv; results/F0_audit/marker_panel_issue_report.tsv (conditional)",
        "F0 never modifies the marker panel.",
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
        exists = (root / required_file).exists()
        add(Path(required_file).name, f"{required_file} exists", "exists" if exists else "missing", pass_fail(exists), "blocking", required_file)
    return rows


def write_reports(
    root: Path,
    sample_info: Sequence[Dict[str, str]],
    data_audit: Sequence[Dict[str, str]],
    gate_rows: Sequence[Dict[str, object]],
) -> None:
    group_counts: Dict[str, int] = {}
    for row in sample_info:
        group = row.get("group_analysis", "")
        group_counts[group] = group_counts.get(group, 0) + 1
    n_enter = sum(1 for row in data_audit if row.get("audit_decision") == "enter_full_F1")
    blocking_failed = [
        row for row in gate_rows if row["blocking_level"] == "blocking" and row["pass_fail"] == "FAIL"
    ]
    limited = any(row["pass_fail"] == "PASS_WITH_NOTED_ISSUES" for row in gate_rows)
    gate = "FAIL" if blocking_failed or n_enter != 40 else "PASS_WITH_NOTED_LIMITATIONS" if limited else "PASS"
    run_id = current_run_id()
    recon = [
        "# F0 Global Data Reconnaissance Report", "", f"Run ID: {run_id}", f"Generated at: {now_iso()}", "",
        "## Current Usable Data", "",
        f"- GSE183904: 40 author-filtered raw gene-count CSV.gz matrices; group counts: {group_counts}.",
        "- GSE183904 is the only dataset eligible to start F1 after F0 gate approval.",
        "- GSE239676 passed a structure-only preaudit (8,630 features, 222,240 cells, 20 patients, PC/LM labels) but remains isolated until F2.4 approval.",
        "- Other bulk, multiomics and external resources remain candidates until section-specific audits.", "",
        "## Boundaries", "",
        "- F0 does not produce biological conclusions.",
        "- GSE183904 lacks FASTQ/raw/empty droplets in current public inputs.",
        "- Normal_Peritoneum is reference/display only; PM sample-level inference is directional because n=3.",
        "- Any formal/precheck mismatch blocks F1.", "",
        "## F0 Gate Summary", "", f"- F0_scRNA_F1_gate: {gate}",
        "- F0_project_data_inventory_status: partial_with_section_specific_audits_pending",
        f"- Samples entering F1 object construction if approved: {n_enter}",
        f"- Blocking checklist failures: {len(blocking_failed)}",
    ]
    (root / "results/F0_audit/F0_global_data_reconnaissance_report.md").write_text("\n".join(recon) + "\n", encoding="utf-8")
    report = [
        "# F0 Execution Report", "", f"Run ID: {run_id}", f"Generated at: {now_iso()}", "",
        "## Inputs Checked", "", "- GSE183904 archive and GEO metadata", "- Existing manifests/prechecks",
        "- GSE239676 structure-preaudit records", "- Read-only marker panel plus integrity-only evidence references", "",
        "## Main Observations", "", f"- sample_info rows: {len(sample_info)}", f"- data_audit rows: {len(data_audit)}",
        f"- enter_full_F1 samples: {n_enter}", f"- group counts: {group_counts}", "",
        "## Gate Decision", "", f"F0_scRNA_F1_gate: {gate}",
        "F0_project_data_inventory_status: partial_with_section_specific_audits_pending", "",
        "F1 may start only after Claude Code reviews the executed outputs and the user approves the gate.",
    ]
    (root / "results/F0_audit/F0_execution_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def execute(root: Path) -> int:
    require_paths(root, STAGE_REQUIRED, STAGE_NAME)
    append_log(root, f"F0 step4 started; run_id={current_run_id()}")
    sample_info = read_tsv(root / "data/metadata/sample_info.tsv")
    data_audit = read_tsv(root / "data/metadata/data_audit.tsv")
    processed_manifest = read_tsv(root / "data/metadata/processed_input_manifest.tsv")

    readiness = build_data_readiness()
    write_tsv(root / "data/metadata/F0_data_readiness_by_F_section.tsv", readiness, [
        "F_section", "required_data_domains", "available_datasets_or_files", "missing_or_pending_items",
        "minimum_data_needed_to_start", "specialized_audit_required_in_section", "current_readiness_status",
        "blocking_for_next_section", "recommended_first_action", "note",
    ])
    methods = build_method_prior_decision()
    write_tsv(root / "data/metadata/F0_method_prior_decision.tsv", methods, [
        "dataset_or_resource_id", "intended_F_sections", "known_data_structure", "known_limitation_from_precheck",
        "default_or_recommended_route", "method_not_allowed_by_default", "required_pre_execution_audit",
        "approval_required_before_change", "sensitivity_or_fallback_route", "interpretation_boundary", "note",
    ])
    decisions = build_decision_evidence_log()
    write_tsv(root / "data/metadata/decision_evidence_log.tsv", decisions, [
        "decision_id", "stage", "decision_topic", "decision_value", "decision_reason", "evidence_type",
        "evidence_source", "source_url_or_file", "evidence_strength", "date", "requires_sensitivity_analysis", "note",
    ])
    excluded = build_excluded_samples(sample_info)
    write_tsv(root / "data/metadata/excluded_samples.tsv", excluded, [
        "sample_id", "geo_accession", "exclusion_scope", "exclusion_status", "reason", "evidence_file", "note",
    ])

    gate_rows = build_gate_checklist(root, sample_info, data_audit, processed_manifest)
    write_tsv(root / "results/F0_audit/F0_gate_checklist.tsv", gate_rows, [
        "gate_item", "required_status", "observed_status", "pass_fail", "blocking_level", "evidence_file", "note",
    ])
    write_reports(root, sample_info, data_audit, gate_rows)
    failures = sum(1 for row in gate_rows if row["pass_fail"] == "FAIL" and row["blocking_level"] == "blocking")
    append_log(
        root,
        f"F0 step4 finalized; gate_rows={len(gate_rows)}; blocking_failures={failures}; "
        "F1 marker-panel/evidence method-prior reminder retained for F1 execution planning.",
    )

    file_rows = build_file_manifest(root, processed_manifest, F0_OUTPUTS)
    write_tsv(root / "data/metadata/F0_file_manifest.tsv", file_rows, [
        "file_name", "relative_path_if_available", "dataset_id", "file_role", "data_domain", "source_type",
        "source_url_or_local_manifest", "file_size_bytes", "sha256", "compression_format", "read_status",
        "availability_status", "used_in_F0", "planned_F_section_use", "audit_status", "artifact_class",
        "publication_destination", "review_priority", "note",
    ])
    if failures:
        raise RuntimeError(f"F0 gate contains {failures} blocking failure(s); pause before F1")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, STAGE_REQUIRED, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
