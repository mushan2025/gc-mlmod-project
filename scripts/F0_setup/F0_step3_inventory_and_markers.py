#!/usr/bin/env python3
"""F0 Step3：整理数据清单、处理历史证据，并审计 marker panel。

为什么要做：Step2 证明了矩阵能否计算，但还需要回答“项目有哪些数据和资源、
每份数据能支持哪一步、作者公开说明了哪些处理、哪些环节仍未知”。Step3 把
这些信息整理成结构化表，供 Step4 作 gate 判断，也供后续 F1-F8 追溯。

主要输入：成功完成的 Step1/Step2 输出、项目已有 manifest、各数据集预检查表、
处理历史来源审计，以及只读 marker panel。

主要输出：数据集/文件/metadata/处理历史/外部资源清单；仅在 marker panel 发现
问题时生成 ``marker_panel_issue_report.tsv``。

重要边界：marker 证据表只用于检查 evidence ID 是否存在。脚本不会自动修改、
补充或扩展 marker panel，避免程序未经研究者批准改变细胞注释依据。
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
    "docs/source_verification/GSE183904_processing_history_source_audit.tsv",
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

# marker 基因仅做格式审计，不在本步骤评价其生物学正确性。
MARKER_GENE_PATTERN = re.compile(r"^[A-Z][A-Z0-9-]*$")
MARKER_REQUIRED_FIELDS = [
    "cell_type",
    "positive_markers",
    "minimum_rule",
    "confidence",
    "evidence_ids",
]

# 处理史表同时容纳来源声明和本项目观察，但用两个状态字段严格区分：
# author_reported_status = 原论文/GEO/工具说明实际报告了什么；
# record_status = F0 重算或跨来源核对得到什么。二者不能互相替代。
PROCESSING_HISTORY_FIELDS = [
    "history_record_id",
    "stage_order",
    "dataset_id",
    "processing_step",
    "processing_scope",
    "relation_to_public_matrix",
    "author_reported_status",
    "method_or_threshold_if_reported",
    "source_reference_or_file",
    "evidence_location",
    "evidence_basis",
    "source_accessed_date",
    "confidence_level",
    "matrix_content_effect",
    "unresolved_detail",
    "implication_for_downstream_plan",
    "requires_special_handling",
    "record_status",
    "observed_public_cell_count",
    "paper_reported_final_tissue_cell_count",
    "count_difference",
    "count_difference_fraction",
    "observed_source_reported_qc_pass_cell_count",
    "observed_fixed_qc_pass_cell_count",
    "observed_fixed_qc_fail_cell_count",
    "observed_fail_nFeature_low_count",
    "observed_fail_nFeature_high_count",
    "observed_fail_nCount_count",
    "observed_fail_percent_mt_count",
    "observed_fail_percent_hb_count",
    "observed_samples_fixed_qc_not_evaluable",
    "observed_sample1_pilot_validation_status",
    "observed_sample_feature_rows_below_min_cells3",
    "observed_working_feature_count_min",
    "observed_working_feature_count_max",
    "observed_samples_working_feature_not_evaluable",
]

# 这 7 个文件构成正式 F0 代码范围；其大写 SHA256 会写入文件 manifest。
F0_SCRIPT_PATHS = [
    "scripts/F0_setup/f0_utils.py",
    "scripts/F0_setup/F0_step1_structure_and_extract.py",
    "scripts/F0_setup/F0_step2_sample_info_and_audit.py",
    "scripts/F0_setup/F0_step3_inventory_and_markers.py",
    "scripts/F0_setup/F0_step4_decisions_and_gate.py",
    "scripts/F0_setup/run_F0_full_audit.py",
    "scripts/F0_setup/validate_F0_readonly.py",
]

# F0 依赖环境与后续 F 节分开管理。这两个文件是 pip/Conda 两条可替代路线，
# 都进入 manifest 锁定，但正式机器只需采用其中一条获批路线。
F0_ENVIRONMENT_LOCK_PATHS = [
    "environment/F0/requirements.txt",
    "environment/F0/environment.yml",
]


def summarize_unique(values: Iterable[str], max_items: int = 8) -> str:
    """压缩显示字段的不同取值及频数，避免清单单元格无限变长。"""

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
    """汇总各数据集 metadata 字段的记录数、缺失情况和可用边界。

    这张表帮助后续区分“字段真实不存在”和“分析者忘记读取”。它只登记字段
    可用性，不把尚未审核的数据集自动批准为验证队列。
    """

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


def build_author_processing_audit(
    source_rows: Sequence[Dict[str, str]],
    data_audit: Sequence[Dict[str, str]],
) -> List[Dict[str, object]]:
    """把来源处理史与 F0 实测结果合并成一条可审计时间线。

    来源行必须覆盖组织解离、建库、测序、count 生成、cell calling、QC、导出、
    doublet、ambient、标准化、整合/聚类及 raw droplets 可用性等关键阶段。

    F0 另外加入三类记录：固定 QC 重算、min.cells=3 工作空间重算，以及
    158,641 个公开细胞与论文 152,423 个最终细胞的跨来源核对。细胞数差异只
    记录事实，不擅自归因于 DoubletFinder 或其他单一步骤。
    """

    if not source_rows:
        raise RuntimeError("GSE183904 processing-history source audit is empty")
    missing_fields = [field for field in PROCESSING_HISTORY_FIELDS[:17] if field not in source_rows[0]]
    if missing_fields:
        raise RuntimeError(
            "GSE183904 processing-history source audit missing fields: "
            + ", ".join(missing_fields)
        )
    record_ids = [row.get("history_record_id", "") for row in source_rows]
    if len(record_ids) != len(set(record_ids)) or any(not value for value in record_ids):
        raise RuntimeError("GSE183904 processing-history source audit has blank or duplicate record IDs")
    required_steps = {
        "tissue_dissociation",
        "library_preparation",
        "sequencing",
        "count_generation",
        "cell_calling",
        "QC_filtering",
        "public_matrix_export",
        "doublet_detection",
        "ambient_RNA_correction",
        "normalization",
        "batch_correction",
        "clustering_and_annotation",
        "raw_FASTQ_or_empty_droplet_availability",
    }
    observed_steps = {
        row.get("processing_step", "")
        for row in source_rows
        if row.get("dataset_id") == "GSE183904"
    }
    missing_steps = sorted(required_steps - observed_steps)
    if missing_steps:
        raise RuntimeError(
            "GSE183904 processing-history source audit missing stages: "
            + ", ".join(missing_steps)
        )

    # 先原样保留来源审计，再追加 F0 自己产生的观察记录。
    rows: List[Dict[str, object]] = [dict(row) for row in source_rows]
    public_cells = sum(int(row.get("matrix_cols_cells", "0") or 0) for row in data_audit)
    source_qc_pass_cells = sum(
        int(row.get("source_reported_qc_pass_count", "0") or 0) for row in data_audit
    )
    fixed_qc_pass_cells = sum(
        int(row.get("final_fixed_qc_pass_count", "0") or 0) for row in data_audit
    )
    fixed_qc_fail_cells = sum(
        int(row.get("final_fixed_qc_fail_count", "0") or 0) for row in data_audit
    )
    rule_fail_counts = {
        field: sum(int(row.get(field, "0") or 0) for row in data_audit)
        for field in (
            "fail_nFeature_low_count",
            "fail_nFeature_high_count",
            "fail_nCount_count",
            "fail_percent_mt_count",
            "fail_percent_hb_count",
        )
    }
    fixed_qc_not_evaluable_rows = [
        row for row in data_audit if row.get("fixed_qc_rule_recalculation_status") != "pass"
    ]
    pilot_rows = [
        row for row in data_audit if row.get("pilot_validation_applicable") == "true"
    ]
    pilot_status = (
        pilot_rows[0].get("pilot_validation_status", "missing")
        if len(pilot_rows) == 1
        else f"invalid_pilot_row_count_{len(pilot_rows)}"
    )
    fixed_qc_evaluable = (
        bool(data_audit)
        and not fixed_qc_not_evaluable_rows
        and fixed_qc_pass_cells + fixed_qc_fail_cells == public_cells
        and pilot_status == "pass"
    )
    sample_feature_rows_below_3 = sum(
        int(row.get("feature_rows_detected_lt_3_count", "0") or 0)
        for row in data_audit
    )
    working_feature_not_evaluable_rows = [
        row
        for row in data_audit
        if row.get("working_feature_space_recalculation_status") != "pass"
    ]
    working_feature_counts = [
        int(row.get("qc_retained_feature_count", "0") or 0)
        for row in data_audit
        if row.get("working_feature_space_recalculation_status") == "pass"
    ]
    working_feature_evaluable = (
        bool(data_audit)
        and not working_feature_not_evaluable_rows
        and len(working_feature_counts) == len(data_audit)
        and all(value > 0 for value in working_feature_counts)
    )
    fixed_qc_status = (
        "F0_fixed_QC_recalculation_pass_all_samples"
        if fixed_qc_evaluable
        else "not_evaluable"
    )
    rows.append(
        {
            "history_record_id": "H018",
            "stage_order": "92",
            "dataset_id": "GSE183904",
            "processing_step": "fixed_QC_rule_recalculation",
            "processing_scope": "public_matrix_verification",
            "relation_to_public_matrix": "public_matrix_observed",
            "author_reported_status": "not_applicable_F0_observation",
            "record_status": fixed_qc_status,
            "method_or_threshold_if_reported": (
                f"For {public_cells} public cells, F0 first retained per-sample features detected in at "
                "least 3 cells, then recomputed nCount, nFeature, percent.mt and percent.HB. Source-reported "
                "rules were 500<=nFeature<6000 and percent.mt<=20; project rules were nCount>1000 and "
                f"percent.HB<5. Source-rule pass={source_qc_pass_cells}; final fixed-QC pass={fixed_qc_pass_cells}."
            ),
            "source_reference_or_file": "data/metadata/data_audit.tsv",
            "evidence_location": "F0 full-stream min.cells=3 working-space QC fields and frozen sample1 regression",
            "evidence_basis": "F0_full_stream_recomputation",
            "source_accessed_date": now_iso()[:10],
            "confidence_level": "high" if fixed_qc_evaluable else "not_evaluable",
            "matrix_content_effect": (
                "Quantifies which currently public cells would be retained by the approved project QC rule; "
                "it does not alter the archived 26571-row matrices in F0."
            ),
            "unresolved_detail": (
                "Cell Ranger-called barcodes and author per-step excluded-barcode lists are unavailable, so "
                "the public cells' complete pre-export filtering history cannot be reconstructed."
            ),
            "implication_for_downstream_plan": (
                "F1 must independently apply the same min.cells=3 working space and fixed QC inequalities, "
                "then record every cell-level decision before doublet and ambient-RNA assessment."
            ),
            "requires_special_handling": (
                "Do not claim that this recalculation reproduces the author's final object, and do not target "
                "the paper's final cell count."
            ),
            "observed_public_cell_count": public_cells,
            "observed_source_reported_qc_pass_cell_count": source_qc_pass_cells,
            "observed_fixed_qc_pass_cell_count": fixed_qc_pass_cells,
            "observed_fixed_qc_fail_cell_count": fixed_qc_fail_cells,
            "observed_fail_nFeature_low_count": rule_fail_counts["fail_nFeature_low_count"],
            "observed_fail_nFeature_high_count": rule_fail_counts["fail_nFeature_high_count"],
            "observed_fail_nCount_count": rule_fail_counts["fail_nCount_count"],
            "observed_fail_percent_mt_count": rule_fail_counts["fail_percent_mt_count"],
            "observed_fail_percent_hb_count": rule_fail_counts["fail_percent_hb_count"],
            "observed_samples_fixed_qc_not_evaluable": len(fixed_qc_not_evaluable_rows),
            "observed_sample1_pilot_validation_status": pilot_status,
        }
    )
    working_feature_status = (
        "F0_working_feature_space_recalculation_pass_all_samples"
        if working_feature_evaluable
        else "not_evaluable"
    )
    rows.append(
        {
            "history_record_id": "H020",
            "stage_order": "94",
            "dataset_id": "GSE183904",
            "processing_step": "working_feature_space_recalculation",
            "processing_scope": "public_matrix_verification",
            "relation_to_public_matrix": "public_matrix_observed",
            "author_reported_status": "not_applicable_F0_observation",
            "record_status": working_feature_status,
            "method_or_threshold_if_reported": (
                "The original paper reports considering per-sample features detected in at least 3 cells. "
                f"F0 found {sample_feature_rows_below_3} sample-by-gene rows below that threshold; the "
                f"resulting working feature counts range from {min(working_feature_counts) if working_feature_counts else 'NA'} "
                f"to {max(working_feature_counts) if working_feature_counts else 'NA'} across {len(data_audit)} samples."
            ),
            "source_reference_or_file": (
                "data/metadata/data_audit.tsv; "
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC9394383/"
            ),
            "evidence_location": (
                "F0 full-stream per-gene detected-cell counts and original-paper Bioinformatic QC methods"
            ),
            "evidence_basis": "F0_full_stream_recomputation_plus_original_paper",
            "source_accessed_date": now_iso()[:10],
            "confidence_level": "high" if working_feature_evaluable else "not_evaluable",
            "matrix_content_effect": (
                "Defines the only per-sample feature space used to compute F0/F1 cell-QC metrics while "
                "leaving the archived raw count rows intact."
            ),
            "unresolved_detail": (
                "Public CSV feature rows do not identify the exact Seurat object export step; low-detection "
                "rows in the archive therefore cannot establish whether the author's working object was filtered."
            ),
            "implication_for_downstream_plan": (
                "Preserve every archived CSV feature row unchanged; F1 recomputes min.cells=3 per sample for "
                "QC/object construction, while DE, pseudobulk and scoring retain method-specific coverage rules."
            ),
            "requires_special_handling": (
                "Do not label low-detection archived rows as an author-processing mismatch and do not use "
                "min.cells=3 as a permanent downstream gene filter."
            ),
            "observed_public_cell_count": public_cells,
            "observed_sample_feature_rows_below_min_cells3": sample_feature_rows_below_3,
            "observed_working_feature_count_min": (
                min(working_feature_counts) if working_feature_counts else ""
            ),
            "observed_working_feature_count_max": (
                max(working_feature_counts) if working_feature_counts else ""
            ),
            "observed_samples_working_feature_not_evaluable": len(
                working_feature_not_evaluable_rows
            ),
        }
    )

    paper_final_cells = 152_423
    count_difference = public_cells - paper_final_cells
    difference_fraction = count_difference / public_cells if public_cells else 0.0
    rows.append(
        {
            "history_record_id": "H019",
            "stage_order": "125",
            "dataset_id": "GSE183904",
            "processing_step": "export_boundary_reconciliation",
            "processing_scope": "cross_source_reconciliation",
            "relation_to_public_matrix": "export_boundary_unresolved",
            "author_reported_status": "not_applicable_cross_source_reconciliation",
            "record_status": "cross_source_count_difference_observed",
            "method_or_threshold_if_reported": (
                f"Public matrices contain {public_cells} cells; the original-paper 40-tissue-sample "
                f"analysis reports {paper_final_cells}; difference={count_difference} "
                f"({difference_fraction:.4%} of public cells)."
            ),
            "source_reference_or_file": (
                "data/metadata/data_audit.tsv; "
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC9394383/"
            ),
            "evidence_location": "F0 matrix column counts and original-paper Figure 1/results",
            "evidence_basis": "F0_local_matrix_plus_original_paper",
            "source_accessed_date": now_iso()[:10],
            "confidence_level": "high_for_counts_medium_for_cause",
            "matrix_content_effect": "Shows that the public matrices are not numerically identical to the author's final 40-tissue-sample analysis object.",
            "unresolved_detail": (
                "The difference is consistent with later doublet or other post-export filtering but "
                "cannot identify the cause or affected barcodes without author code and exclusion lists."
            ),
            "implication_for_downstream_plan": (
                "Treat public-matrix doublet status as unresolved and reconstruct F1 exclusions de novo."
            ),
            "requires_special_handling": (
                "Do not force F1 to reproduce 152423 cells and do not attribute all 6218 cells to DoubletFinder."
            ),
            "observed_public_cell_count": public_cells,
            "paper_reported_final_tissue_cell_count": paper_final_cells,
            "count_difference": count_difference,
            "count_difference_fraction": format(difference_fraction, ".12g"),
            "observed_source_reported_qc_pass_cell_count": source_qc_pass_cells,
            "observed_fixed_qc_pass_cell_count": fixed_qc_pass_cells,
            "observed_fixed_qc_fail_cell_count": fixed_qc_fail_cells,
            "observed_samples_fixed_qc_not_evaluable": len(fixed_qc_not_evaluable_rows),
            "observed_sample1_pilot_validation_status": pilot_status,
        }
    )
    rows.append(
        {
            "history_record_id": "H900",
            "stage_order": "900",
            "dataset_id": "GSE206785",
            "processing_step": "external_scRNA_processing",
            "processing_scope": "external_validation_candidate",
            "relation_to_public_matrix": "separate_dataset",
            "author_reported_status": "not_applicable_project_inventory",
            "record_status": "processed_matrix_local_candidate",
            "method_or_threshold_if_reported": "log1p(count)-like cell-by-gene matrix",
            "source_reference_or_file": "results/F0_audit/GSE206785_dataset_structure_precheck.tsv",
            "evidence_location": "results/F0_audit/GSE206785_dataset_structure_precheck.tsv",
            "evidence_basis": "project_local_structure_precheck",
            "source_accessed_date": now_iso()[:10],
            "confidence_level": "medium",
            "matrix_content_effect": "Processed expression only; no raw-QC reconstruction.",
            "unresolved_detail": "Full external-cohort audit remains pending.",
            "implication_for_downstream_plan": "Processed-expression external scoring only after approval.",
            "requires_special_handling": "F2 external-validation audit required before use.",
        }
    )
    return sorted(rows, key=lambda row: (row.get("dataset_id", ""), int(row.get("stage_order", 0))))


def gse239_inventory_row(predownloaded: Sequence[Dict[str, str]]) -> Dict[str, object] | None:
    """把 GSE239676 的多文件预检查压缩为一条外部验证候选记录。"""

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
    """建立全项目数据集总表，并注明用途、现状、限制和下一次专项审计。

    “本地已有”不等于“已经批准使用”。外部验证集仍必须遵守签名冻结和队列
    隔离规则，避免验证数据反向参与参数选择。
    """

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
            "raw_or_processed_status_if_known": "Cell Ranger-derived public raw gene count CSV.gz; F0 separately audits retained-cell and feature-filter boundaries; not FASTQ/raw droplets",
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
    """建立文件级来源、大小、SHA256、用途和审核优先级清单。

    正式 F0 的 7 个脚本和 2 个环境锁文件也进入此表，从而能够证明某次结果
    由哪一版代码和依赖规格生成。
    manifest 自身的 SHA256 故意留空，因为把自己的校验和写入自身会形成循环。
    """

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
    processing_source = root / "docs/source_verification/GSE183904_processing_history_source_audit.tsv"
    rows.append(
        {
            "file_name": processing_source.name,
            "relative_path_if_available": rel(processing_source, root),
            "dataset_id": "GSE183904",
            "file_role": "processing_history_primary_source_audit",
            "data_domain": "scientific_provenance",
            "source_type": "project_documented_primary_source_audit",
            "source_url_or_local_manifest": rel(processing_source, root),
            "file_size_bytes": processing_source.stat().st_size if processing_source.exists() else "",
            "sha256": sha256_file(processing_source) if processing_source.exists() else "",
            "compression_format": "tsv",
            "read_status": "read_by_F0_step3" if processing_source.exists() else "missing",
            "availability_status": "project_local_available" if processing_source.exists() else "not_found_at_expected_path",
            "used_in_F0": "true",
            "planned_F_section_use": "F0,F1",
            "audit_status": "primary_sources_reviewed_and_structured",
            "artifact_class": "audit_trail",
            "publication_destination": "not_planned",
            "review_priority": "full",
            "note": "Separates author report, F0 verification, cross-source inference and unresolved history.",
        }
    )
    for script_relative_path in F0_SCRIPT_PATHS:
        script_path = root / script_relative_path
        rows.append(
            {
                "file_name": script_path.name,
                "relative_path_if_available": script_relative_path,
                "dataset_id": "F0",
                "file_role": (
                    "F0_readonly_validation_script"
                    if script_path.name == "validate_F0_readonly.py"
                    else "F0_formal_execution_script"
                ),
                "data_domain": "execution_code",
                "source_type": "git_tracked_project_code",
                "source_url_or_local_manifest": "git repository and current F0_file_manifest.tsv",
                "file_size_bytes": script_path.stat().st_size if script_path.exists() else "",
                "sha256": sha256_file(script_path) if script_path.exists() else "",
                "compression_format": "py",
                "read_status": "checksum_recorded" if script_path.exists() else "missing",
                "availability_status": (
                    "project_local_available" if script_path.exists() else "not_found_at_expected_path"
                ),
                "used_in_F0": (
                    "pre_execution_validation"
                    if script_path.name == "validate_F0_readonly.py"
                    else "true"
                ),
                "planned_F_section_use": "F0",
                "audit_status": (
                    "code_checksum_locked_at_F0_execution" if script_path.exists() else "missing"
                ),
                "artifact_class": "audit_trail",
                "publication_destination": "not_planned",
                "review_priority": "full",
                "note": "SHA256 is written in uppercase; git commit records the reviewed source state.",
            }
        )
    for lock_relative_path in F0_ENVIRONMENT_LOCK_PATHS:
        lock_path = root / lock_relative_path
        is_conda = lock_path.name == "environment.yml"
        rows.append(
            {
                "file_name": lock_path.name,
                "relative_path_if_available": lock_relative_path,
                "dataset_id": "F0",
                "file_role": (
                    "F0_conda_environment_specification"
                    if is_conda
                    else "F0_pip_dependency_lock"
                ),
                "data_domain": "execution_environment",
                "source_type": "git_tracked_F0_environment_lock",
                "source_url_or_local_manifest": "environment/environment_lock_manifest.tsv",
                "file_size_bytes": lock_path.stat().st_size if lock_path.exists() else "",
                "sha256": sha256_file(lock_path) if lock_path.exists() else "",
                "compression_format": lock_path.suffix.lstrip("."),
                "read_status": "checksum_recorded" if lock_path.exists() else "missing",
                "availability_status": (
                    "project_local_available" if lock_path.exists() else "not_found_at_expected_path"
                ),
                "used_in_F0": "environment_reproduction_contract",
                "planned_F_section_use": "F0_only",
                "audit_status": (
                    "environment_checksum_locked_at_F0_execution"
                    if lock_path.exists()
                    else "missing"
                ),
                "artifact_class": "audit_trail",
                "publication_destination": "not_planned",
                "review_priority": "full",
                "note": (
                    "F0-only dependency specification; pip and Conda are alternative setup routes, "
                    "not two environments to combine-install."
                ),
            }
        )
    for output in generated_paths:
        path = root / output
        output_name = Path(output).name
        is_self = output_name == "F0_file_manifest.tsv"
        gate_decision_names = {
            "F0_author_processing_audit.tsv",
            "F0_gate_checklist.tsv",
            "F0_execution_report.md",
        }
        rows.append(
            {
                "file_name": output_name,
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
                "artifact_class": "gate_decision" if output_name in gate_decision_names else "audit_trail",
                "publication_destination": "undecided" if output_name in {"F0_global_data_reconnaissance_report.md", "F0_execution_report.md"} else "not_planned",
                "review_priority": "full" if output_name == "F0_author_processing_audit.tsv" else "standard",
                "note": "Self-manifest SHA256 intentionally blank" if is_self else "",
            }
        )
    return rows


def build_external_resource_inventory(root: Path) -> List[Dict[str, object]]:
    """登记 SCENIC、inferCNV 和后续数据库资源的当前可用状态。

    尚未做版本、许可或兼容性审计的数据库只标记为候选，F0 不会下载或使用它们。
    """

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
    """只读检查 marker panel 的必需列、格式、重复项和证据 ID 完整性。

    检查结果均为 warning，不在 F0 中自动改 panel。细胞类型 marker 的生物学
    合理性将在 F1 注释阶段结合表达结果和来源证据审核。
    """

    panel_path = root / "data/metadata/cell_type_marker_panel.tsv"
    evidence_path = root / "data/metadata/cell_marker_reference_evidence.tsv"
    panel = read_tsv(panel_path)
    evidence = read_tsv(evidence_path) if evidence_path.exists() else []
    panel_fields = set(panel[0].keys()) if panel else set()
    evidence_ids = {row.get("evidence_id", "") for row in evidence}
    issues: List[Dict[str, object]] = []

    def add(cell_type: str, issue_type: str, field: str, observed: str, action: str) -> None:
        """以统一字段记录一条 marker panel 审计问题。"""

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
    """正式执行 Step3，生成清单和处理史；Step2 有暂停样本时拒绝继续。"""

    # 1. 先确认 Step1/Step2 三张关键表均有 40 行且全部允许独立重新 QC。
    require_paths(root, STAGE_REQUIRED, STAGE_NAME)
    append_log(root, f"F0 step3 started; run_id={current_run_id()}")
    sample_info = read_tsv(root / "data/metadata/sample_info.tsv")
    processed_manifest = read_tsv(root / "data/metadata/processed_input_manifest.tsv")
    data_audit = read_tsv(root / "data/metadata/data_audit.tsv")
    if len(processed_manifest) != 40 or len(sample_info) != 40 or len(data_audit) != 40:
        raise RuntimeError(
            "F0 step3 requires 40 processed-manifest, sample-info and data-audit rows"
        )
    if any(
        row.get("audit_decision") != "enter_full_F1_independent_reQC"
        for row in data_audit
    ):
        raise RuntimeError("F0 step3 cannot continue while data_audit.tsv contains paused samples")
    # 2. 收集主队列和候选数据集的 metadata/结构预检查信息。
    sample_fields, _ = parse_gse183904_series(
        root / "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz"
    )
    prechecks = {
        "non_GSE183904": read_tsv(root / "results/F0_audit/non_GSE183904_data_structure_precheck.tsv"),
        "predownloaded": read_tsv(root / "results/F0_audit/predownloaded_resource_structure_audit.tsv"),
        "bulk_GEO_series_matrix_deep_precheck": read_tsv(root / "results/F0_audit/bulk_GEO_series_matrix_deep_precheck.tsv"),
        "GSE206785_metadata_precheck": read_tsv(root / "results/F0_audit/GSE206785_metadata_precheck.tsv"),
    }

    # 3. 分别生成数据集、metadata、处理历史和外部资源清单。
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
    processing_source_rows = read_tsv(
        root / "docs/source_verification/GSE183904_processing_history_source_audit.tsv"
    )
    author_rows = build_author_processing_audit(processing_source_rows, data_audit)
    write_tsv(
        root / "data/metadata/F0_author_processing_audit.tsv",
        author_rows,
        PROCESSING_HISTORY_FIELDS,
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

    # 4. marker panel 保持只读；只有发现问题时才输出单独的 warning 表。
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

    # 5. 最后登记输入、输出和脚本 SHA256，锁定本次 F0 所用文件版本。
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
        f"author_processing_rows={len(author_rows)}; resource_rows={len(resource_rows)}; "
        f"marker_panel={marker_status}; F1 reminder: use processing-history constraints and "
        "add marker-panel/evidence-reference method priors.",
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """默认 dry run；显式提供 ``--execute`` 后才写正式清单。"""

    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, STAGE_REQUIRED, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
