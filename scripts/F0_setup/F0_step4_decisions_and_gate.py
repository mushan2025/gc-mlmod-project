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
            "available_datasets_or_files": "GSE183904_RAW.tar; GSE183904_series_matrix.txt.gz; processing-history source audit; cell_type_marker_panel.tsv",
            "missing_or_pending_items": "Exact Cell Ranger cell-calling settings, DoubletFinder parameters/barcodes and public export timing remain unknown; F1 R dependencies remain required",
            "minimum_data_needed_to_start": "F0_scRNA_F1_gate PASS or PASS_WITH_NOTED_LIMITATIONS",
            "specialized_audit_required_in_section": "Gate1 executes the preregistered six-layer per-sample QC, freezes nmads 4/3/2 masks, and performs per-capture doublet plus retained-cell ambient review",
            "current_readiness_status": "ready_with_limitations",
            "blocking_for_next_section": "yes_until_F0_gate_passes",
            "recommended_first_action": "Run F1 Gate1 plan only after F0 review and user approval",
            "note": "Author SCTransform/integration/clustering are post-export context, not transformations embedded in the public raw counts; F1 must not use MLMOD score/signature/prognosis.",
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
            "known_data_structure": "Cell Ranger-derived public called/retained-cell raw gene-count CSV.gz; gene rows by cell columns; F0 measures public-full and author-like min.cells=3 QC spaces separately",
            "known_limitation_from_precheck": "No FASTQ/raw/empty droplets; public feature rows do not necessarily reproduce the author's per-sample min.cells=3 analysis space; cell metrics can be feature-space dependent near thresholds; author used DoubletFinder but public-matrix doublet status, parameters and excluded barcodes are unresolved; PM n=3; Normal_Peritoneum n=1",
            "default_or_recommended_route": "F1 uses the preregistered six-layer per-sample QC framework, public_full_feature_space for independent QC, author-like space for provenance only, per-sample scDblFinder, retained-cell DecontX diagnostics and frozen nmads 4/3/2 masks",
            "method_not_allowed_by_default": "No FASTQ/Cell Ranger reproduction claim; no true knee/emptyDrops/SoupX/CellBender without raw droplets; no MLMOD or outcome-guided QC; no permanent gene deletion from author min.cells=3",
            "required_pre_execution_audit": "F0 data_audit plus processing-history/export-boundary audit and F1 Gate1 QC/doublet/ambient plan",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "nmads 4/3/2 cell-inclusion masks; scDblFinder dbr 0.5x/1x/1.5x; DoubletFinder same-matrix method sensitivity; if count/QC provenance is not evaluable, pause and revise F1 input strategy",
            "interpretation_boundary": "F1 produces a reusable annotated object, not MLMOD conclusions",
            "note": "Do not force F1 to reproduce the paper's 152423-cell object; all local exclusions require new barcode-level reasons and at least two preregistered, nonredundant evidence families unless hard-anomaly. This rule prevents duplicate counting and does not claim statistical independence. Normal_Peritoneum is reference-only and PM statistics are directional only.",
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
            "decision_value": "Use Cell Ranger-derived public raw gene-count matrices with separately verified cell-QC and feature-filter boundaries, not FASTQ/raw droplets",
            "decision_reason": "GEO reports raw files not submitted; supplementary matrices are post-Cell-Ranger counts, and F0 independently tests both the reported retained-cell thresholds and per-sample min.cells=3 feature rule.",
            "evidence_type": "original_paper_GEO_metadata_and_stream_audit",
            "evidence_source": "GSE183904 processing-history source audit; GSE183904_series_matrix.txt.gz; data_audit.tsv",
            "source_url_or_file": "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
            "evidence_strength": "high",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "Author SCTransform, integration and clustering are not encoded in the public integer counts; F1 must verify residual QC/doublet/ambient risk and must not assume public feature rows equal the author Seurat analysis space.",
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
        {
            "decision_id": "F0_DECISION_005",
            "stage": "F0/F1",
            "decision_topic": "GSE183904_doublet_and_export_boundary",
            "decision_value": "Author DoubletFinder use is documented, but public-matrix doublet status remains unresolved",
            "decision_reason": "The original paper reports DoubletFinder removal, while public matrices and the paper final 40-tissue-sample object have different cell totals without a barcode-level exclusion trace.",
            "evidence_type": "original_paper_local_matrix_cross_source_reconciliation",
            "evidence_source": "F0_author_processing_audit.tsv; data_audit.tsv; Kumar et al. Cancer Discovery 2022",
            "source_url_or_file": "data/metadata/F0_author_processing_audit.tsv",
            "evidence_strength": "high_for_counts_and_reported_method_medium_for_export_timing",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "F1 uses per-sample scDblFinder as the preregistered primary call and records every local exclusion; the paper final cell count is not a target.",
        },
        {
            "decision_id": "F0_DECISION_006",
            "stage": "F0/F1",
            "decision_topic": "GSE183904_feature_filter_boundary",
            "decision_value": "Public feature rows are audited independently of retained-cell QC",
            "decision_reason": "The original paper reports a per-sample three-cell feature rule, which can be tested directly in each public matrix and need not be inferred from the cell thresholds.",
            "evidence_type": "original_paper_and_full_stream_sample_by_gene_audit",
            "evidence_source": "F0_author_processing_audit.tsv; data_audit.tsv; Kumar et al. Cancer Discovery 2022",
            "source_url_or_file": "data/metadata/F0_author_processing_audit.tsv",
            "evidence_strength": "high_after_formal_full_stream_audit",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "If low-detection public feature rows are present, F1 preserves archived counts, uses the author-like space only for provenance and applies separate feature-coverage rules for each downstream method.",
        },
        {
            "decision_id": "F0_DECISION_007",
            "stage": "F0/F1",
            "decision_topic": "F1_six_layer_QC_framework",
            "decision_value": "Use per-sample diagnostics, hard-anomaly safety net, multi-family QC, provisional diagnostic clustering, per-capture ambient/doublet assessment and frozen lenient/main/strict masks",
            "decision_reason": "Fixed global thresholds are not sufficient for heterogeneous gastric samples, while raw-droplet-dependent steps cannot be reconstructed from the public called-cell matrices.",
            "evidence_type": "method_papers_official_documentation_and_project_input_audit",
            "evidence_source": "scuttle/OSCA; miQC; SampleQC; scDblFinder; emptyDrops; SoupX; DecontX; GSE183904 data_audit",
            "source_url_or_file": "胃癌MLMOD亚群主线研究方案.txt",
            "evidence_strength": "high_for_framework_project_specific_for_operational_thresholds",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "Barcode knee, Cell Ranger cell calling, emptyDrops, SoupX and CellBender are not_evaluable_input_limited; same-matrix doublet methods cannot supply a second nonredundant evidence family, and no same-matrix family is claimed to be statistically independent evidence.",
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
    author_processing: Sequence[Dict[str, str]],
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
        row.get("audit_decision") == "enter_full_F1_independent_reQC"
        and row.get("normalization_artifact_flag") == "false"
        and row.get("observed_numeric_type") == "nonnegative_integer_count_like"
        and row.get("suspected_matrix_type") == "public_called_cell_raw_gene_count_matrix"
        and row.get("processed_input_manifest_match") == "true"
        and row.get("per_gene_mean_distribution_status") == "consistent_with_sparse_right_skew"
        and row.get("author_cell_qc_reproduction_status") in {"pass", "measured_mismatch"}
        and row.get("author_cell_qc_reproduction_status_public_space") in {"pass", "measured_mismatch"}
        and row.get("author_cell_qc_reproduction_status_author_like_space") in {"pass", "measured_mismatch"}
        and row.get("author_feature_filter_reproduction_status") in {"pass", "measured_mismatch"}
        and row.get("format_decision_scope") == "file_format_only"
        and row.get("decision_scope") == "file_format_and_public_processing_boundary"
        for row in data_audit
        if row.get("include_in_f1") == "true"
    )
    add(
        "data_audit",
        "40 rows; included samples pass numeric/distribution/mapping checks; both QC spaces and processing boundaries are evaluable",
        f"{len(data_audit)} rows; pause={sum(1 for row in data_audit if row.get('audit_decision') != 'enter_full_F1_independent_reQC')}",
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
    history_rows = [
        row for row in author_processing if row.get("dataset_id") == "GSE183904"
    ]
    required_history_steps = {
        "tissue_dissociation",
        "library_preparation",
        "sequencing",
        "count_generation",
        "cell_calling",
        "QC_filtering",
        "cell_QC_filtering_verification",
        "feature_filtering_verification",
        "public_matrix_export",
        "doublet_detection",
        "export_boundary_reconciliation",
        "ambient_RNA_correction",
        "normalization",
        "batch_correction",
        "clustering_and_annotation",
        "raw_FASTQ_or_empty_droplet_availability",
    }
    history_steps = {row.get("processing_step", "") for row in history_rows}
    history_core_fields = [
        "history_record_id",
        "processing_scope",
        "relation_to_public_matrix",
        "author_reported_status",
        "source_reference_or_file",
        "evidence_basis",
        "confidence_level",
        "unresolved_detail",
        "implication_for_downstream_plan",
        "requires_special_handling",
    ]
    mapped_rows_ok = bool(history_rows) and all(
        all(row.get(field, "").strip() for field in history_core_fields)
        for row in history_rows
    )
    doublet_row = next(
        (row for row in history_rows if row.get("processing_step") == "doublet_detection"),
        {},
    )
    qc_verification_row = next(
        (
            row
            for row in history_rows
            if row.get("processing_step") == "cell_QC_filtering_verification"
        ),
        {},
    )
    feature_filter_verification_row = next(
        (
            row
            for row in history_rows
            if row.get("processing_step") == "feature_filtering_verification"
        ),
        {},
    )
    reconciliation_row = next(
        (row for row in history_rows if row.get("processing_step") == "export_boundary_reconciliation"),
        {},
    )
    public_cell_count = sum(int(row.get("matrix_cols_cells", "0") or 0) for row in data_audit)
    cell_qc_mismatch_count = sum(
        row.get("author_cell_qc_reproduction_status") == "measured_mismatch"
        for row in data_audit
    )
    cell_qc_not_evaluable_count = sum(
        row.get("author_cell_qc_reproduction_status") not in {"pass", "measured_mismatch"}
        for row in data_audit
    )
    public_space_qc_mismatch_cells = sum(
        int(row.get("author_cell_threshold_mismatch_count_public_space", "0") or 0)
        for row in data_audit
    )
    author_like_qc_mismatch_cells = sum(
        int(row.get("author_cell_threshold_mismatch_count_author_like_space", "0") or 0)
        for row in data_audit
    )
    feature_filter_mismatch_count = sum(
        row.get("author_feature_filter_reproduction_status") == "measured_mismatch"
        for row in data_audit
    )
    feature_filter_not_evaluable_count = sum(
        row.get("author_feature_filter_reproduction_status") not in {"pass", "measured_mismatch"}
        for row in data_audit
    )
    sample_gene_rows_below_3 = sum(
        int(row.get("feature_rows_detected_lt_3_count", "0") or 0)
        for row in data_audit
    )
    if cell_qc_not_evaluable_count:
        expected_cell_qc_verification_status = "not_evaluable"
    elif cell_qc_mismatch_count:
        expected_cell_qc_verification_status = "measured_mismatch_in_one_or_both_feature_spaces"
    else:
        expected_cell_qc_verification_status = "verified_in_both_feature_spaces"
    expected_feature_filter_verification_status = (
        "public_feature_rows_include_genes_below_reported_per_sample_threshold"
        if feature_filter_mismatch_count
        else "verified_against_all_public_sample_feature_rows"
    )
    history_ok = (
        required_history_steps.issubset(history_steps)
        and mapped_rows_ok
        and doublet_row.get("author_reported_status") == "reported_in_original_paper"
        and doublet_row.get("relation_to_public_matrix") == "export_boundary_unresolved"
        and qc_verification_row.get("author_reported_status")
        == expected_cell_qc_verification_status
        and qc_verification_row.get("observed_author_qc_mismatch_cell_count_public_space")
        == str(public_space_qc_mismatch_cells)
        and qc_verification_row.get("observed_author_qc_mismatch_cell_count_author_like_space")
        == str(author_like_qc_mismatch_cells)
        and qc_verification_row.get("observed_samples_with_author_cell_qc_mismatch")
        == str(cell_qc_mismatch_count)
        and qc_verification_row.get("observed_samples_author_cell_qc_not_evaluable") == "0"
        and feature_filter_not_evaluable_count == 0
        and feature_filter_verification_row.get("author_reported_status")
        == expected_feature_filter_verification_status
        and feature_filter_verification_row.get("observed_samples_with_feature_filter_mismatch")
        == str(feature_filter_mismatch_count)
        and feature_filter_verification_row.get("observed_sample_gene_rows_below_3")
        == str(sample_gene_rows_below_3)
        and feature_filter_verification_row.get("observed_samples_feature_filter_not_evaluable")
        == "0"
        and reconciliation_row.get("observed_public_cell_count") == str(public_cell_count)
        and reconciliation_row.get("paper_reported_final_tissue_cell_count") == "152423"
        and reconciliation_row.get("count_difference") == str(public_cell_count - 152423)
    )
    add(
        "author_processing_provenance",
        "all matrix-relevant stages sourced; author report/F0 verification/inference/unknown separated; unknowns mapped to downstream action",
        f"{len(history_rows)} GSE183904 rows; missing_stages={','.join(sorted(required_history_steps - history_steps)) or 'none'}; "
        f"doublet_status={doublet_row.get('author_reported_status', 'missing')}; "
        f"QC_verification={qc_verification_row.get('author_reported_status', 'missing')}; "
        f"feature_filter_verification={feature_filter_verification_row.get('author_reported_status', 'missing')}; "
        f"public_vs_paper_cells={reconciliation_row.get('observed_public_cell_count', 'missing')}_vs_"
        f"{reconciliation_row.get('paper_reported_final_tissue_cell_count', 'missing')}",
        pass_fail(history_ok),
        "blocking",
        "docs/source_verification/GSE183904_processing_history_source_audit.tsv; data/metadata/F0_author_processing_audit.tsv",
        "Unknown history is allowed only when explicit and linked to a conservative F1 action.",
    )
    cell_qc_boundary_status = qc_verification_row.get("author_reported_status", "missing")
    add(
        "public_cell_qc_boundary",
        "author cell thresholds are tested in public-full and author-like min.cells=3 feature spaces",
        f"status={cell_qc_boundary_status}; mismatch_samples={cell_qc_mismatch_count}; "
        f"public_space_mismatch_cells={public_space_qc_mismatch_cells}; "
        f"author_like_space_mismatch_cells={author_like_qc_mismatch_cells}",
        (
            "PASS_WITH_NOTED_ISSUES"
            if cell_qc_mismatch_count and cell_qc_not_evaluable_count == 0
            else pass_fail(cell_qc_not_evaluable_count == 0)
        ),
        "nonblocking",
        "data/metadata/data_audit.tsv; data/metadata/F0_author_processing_audit.tsv",
        "A measured near-boundary or space-dependent mismatch limits exact author-object claims but does not invalidate independent F1 re-QC.",
    )
    feature_filter_boundary_status = feature_filter_verification_row.get(
        "author_reported_status", "missing"
    )
    add(
        "public_feature_filter_boundary",
        "author per-sample min.cells=3 rule is tested against every public sample-by-gene row",
        f"status={feature_filter_boundary_status}; mismatch_samples={feature_filter_mismatch_count}; "
        f"sample_gene_rows_below_3={sample_gene_rows_below_3}",
        (
            "PASS_WITH_NOTED_ISSUES"
            if feature_filter_mismatch_count and feature_filter_not_evaluable_count == 0
            else pass_fail(feature_filter_not_evaluable_count == 0)
        ),
        "nonblocking",
        "data/metadata/data_audit.tsv; data/metadata/F0_author_processing_audit.tsv",
        "A measured feature-space export mismatch limits provenance claims but does not invalidate raw counts.",
    )
    unresolved_rows = [
        row
        for row in history_rows
        if row.get("relation_to_public_matrix") in {"export_boundary_unresolved", "availability_boundary"}
        or row.get("author_reported_status", "").startswith("not_")
        or row.get("author_reported_status") == "raw_files_not_submitted"
        or row.get("author_reported_status")
        == "measured_mismatch_in_one_or_both_feature_spaces"
        or row.get("author_reported_status")
        == "public_feature_rows_include_genes_below_reported_per_sample_threshold"
    ]
    add(
        "author_processing_residual_unknowns",
        "unresolved upstream details are explicit and do not masquerade as completed processing",
        f"{len(unresolved_rows)} provenance row(s) retain explicit unknown or unavailable details",
        "PASS_WITH_NOTED_ISSUES" if unresolved_rows else "PASS",
        "nonblocking",
        "data/metadata/F0_author_processing_audit.tsv",
        "These limitations constrain F1 methods and make a clean F0 PASS inappropriate.",
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
    author_processing: Sequence[Dict[str, str]],
) -> None:
    group_counts: Dict[str, int] = {}
    for row in sample_info:
        group = row.get("group_analysis", "")
        group_counts[group] = group_counts.get(group, 0) + 1
    n_enter = sum(
        1
        for row in data_audit
        if row.get("audit_decision") == "enter_full_F1_independent_reQC"
    )
    blocking_failed = [
        row for row in gate_rows if row["blocking_level"] == "blocking" and row["pass_fail"] == "FAIL"
    ]
    limited = any(row["pass_fail"] == "PASS_WITH_NOTED_ISSUES" for row in gate_rows)
    gate = "FAIL" if blocking_failed or n_enter != 40 else "PASS_WITH_NOTED_LIMITATIONS" if limited else "PASS"
    run_id = current_run_id()
    gse_history = [row for row in author_processing if row.get("dataset_id") == "GSE183904"]
    doublet_row = next((row for row in gse_history if row.get("processing_step") == "doublet_detection"), {})
    qc_row = next(
        (
            row
            for row in gse_history
            if row.get("processing_step") == "cell_QC_filtering_verification"
        ),
        {},
    )
    feature_filter_row = next(
        (
            row
            for row in gse_history
            if row.get("processing_step") == "feature_filtering_verification"
        ),
        {},
    )
    reconciliation_row = next((row for row in gse_history if row.get("processing_step") == "export_boundary_reconciliation"), {})
    recon = [
        "# F0 Global Data Reconnaissance Report", "", f"Run ID: {run_id}", f"Generated at: {now_iso()}", "",
        "## Current Usable Data", "",
        f"- GSE183904: 40 Cell Ranger-derived public called/retained-cell raw gene-count CSV.gz matrices with two-space cell-QC and feature-filter boundaries audited separately; group counts: {group_counts}.",
        "- GSE183904 is the only dataset eligible to start F1 after F0 gate approval.",
        "- GSE239676 passed a structure-only preaudit (8,630 features, 222,240 cells, 20 patients, PC/LM labels) but remains isolated until F2.4 approval.",
        "- Other bulk, multiomics and external resources remain candidates until section-specific audits.", "",
        "## GSE183904 Processing Provenance", "",
        "- Wet-lab acquisition, dissociation, 10x 5-prime library preparation, HiSeq4000 sequencing, Cell Ranger v3.0/hg38 count generation and author Seurat thresholds are source-audited.",
        f"- Public-cell threshold verification: {qc_row.get('author_reported_status', 'missing')}; public-full-space mismatch cells={qc_row.get('observed_author_qc_mismatch_cell_count_public_space', 'missing')}; author-like-space mismatch cells={qc_row.get('observed_author_qc_mismatch_cell_count_author_like_space', 'missing')}.",
        f"- Public-feature threshold verification: {feature_filter_row.get('author_reported_status', 'missing')}; mismatch samples={feature_filter_row.get('observed_samples_with_feature_filter_mismatch', 'missing')}; sample-by-gene rows below 3 detected cells={feature_filter_row.get('observed_sample_gene_rows_below_3', 'missing')}.",
        f"- Author doublet method: {doublet_row.get('method_or_threshold_if_reported', 'missing')}",
        f"- Public versus paper final tissue-cell count: {reconciliation_row.get('observed_public_cell_count', 'missing')} versus {reconciliation_row.get('paper_reported_final_tissue_cell_count', 'missing')}; difference={reconciliation_row.get('count_difference', 'missing')}.",
        "- Cell Ranger cell-calling details, exact DoubletFinder parameters/barcodes, ambient correction status and the public export timing relative to doublet removal remain unresolved; true knee/cell calling, emptyDrops, SoupX and CellBender are not_evaluable_input_limited.",
        "- Author SCTransform, integration, clustering and annotation describe the author's downstream object and are not embedded in the public integer-count CSV values.", "",
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
        "## Inputs Checked", "", "- GSE183904 archive, GEO metadata, original-paper processing-history source audit", "- Existing manifests/prechecks",
        "- GSE239676 structure-preaudit records", "- Read-only marker panel plus integrity-only evidence references", "",
        "## Main Observations", "", f"- sample_info rows: {len(sample_info)}", f"- data_audit rows: {len(data_audit)}",
        f"- enter_full_F1_independent_reQC samples: {n_enter}", f"- group counts: {group_counts}",
        f"- author processing-history rows for GSE183904: {len(gse_history)}",
        f"- cell-QC boundary: {qc_row.get('author_reported_status', 'missing')}; public-full-space mismatch cells={qc_row.get('observed_author_qc_mismatch_cell_count_public_space', 'missing')}; author-like-space mismatch cells={qc_row.get('observed_author_qc_mismatch_cell_count_author_like_space', 'missing')}",
        f"- feature-filter boundary: {feature_filter_row.get('author_reported_status', 'missing')}; mismatch samples={feature_filter_row.get('observed_samples_with_feature_filter_mismatch', 'missing')}",
        f"- public/paper tissue-cell reconciliation: {reconciliation_row.get('observed_public_cell_count', 'missing')} / {reconciliation_row.get('paper_reported_final_tissue_cell_count', 'missing')}",
        "- unresolved upstream details are carried forward as explicit F1 constraints", "",
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
    author_processing = read_tsv(root / "data/metadata/F0_author_processing_audit.tsv")

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

    gate_rows = build_gate_checklist(
        root,
        sample_info,
        data_audit,
        processed_manifest,
        author_processing,
    )
    write_tsv(root / "results/F0_audit/F0_gate_checklist.tsv", gate_rows, [
        "gate_item", "required_status", "observed_status", "pass_fail", "blocking_level", "evidence_file", "note",
    ])
    write_reports(root, sample_info, data_audit, gate_rows, author_processing)
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
