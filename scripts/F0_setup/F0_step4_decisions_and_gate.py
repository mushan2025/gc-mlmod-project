#!/usr/bin/env python3
"""F0 Step4：形成方法决策、十项 gate checklist 和最终审计报告。

为什么要做：前 3 步产生了大量事实表，但不能仅凭“脚本跑完”就进入 F1。
Step4 把事实转换成明确的就绪状态、方法限制和通过/阻断判断，并为每个判断
保留证据文件。

主要输入：全部成功的 Step1-Step3 输出。
主要输出：各 F 节数据就绪表、方法前提表、决策证据日志、样本排除表、十项
F0 gate checklist、全局数据摸底报告和正式执行报告。

重要边界：本步骤只完成 F0。即使 gate 通过，也必须再经过 Claude Code 审核
和用户批准；脚本不会自动开始 F1。
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from F0_step3_inventory_and_markers import (
    F0_ENVIRONMENT_LOCK_PATHS,
    F0_SCRIPT_PATHS,
    build_file_manifest,
)
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

# 再次冻结关键契约，用于检查 Step2 输出没有在阶段之间被意外改变。
EXPECTED_GLOBIN_PANEL = "HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ"
PILOT_FILE_NAME = "GSM5573466_sample1.csv.gz"


def build_data_readiness() -> List[Dict[str, object]]:
    """登记 F1-F8 当前具备什么数据、还缺什么以及何时必须暂停。

    这是一张路线导航表，不代表 F2-F8 已获准执行。每个 F 节仍需在开始前做
    自己的数据专项审计和资源评估。
    """

    return [
        {
            "F_section": "F1",
            "required_data_domains": "GSE183904 scRNA count matrices; sample metadata; marker panel",
            "available_datasets_or_files": "GSE183904_RAW.tar; GSE183904_series_matrix.txt.gz; processing-history source audit; cell_type_marker_panel.tsv",
            "missing_or_pending_items": "Exact Cell Ranger cell-calling settings, DoubletFinder parameters/barcodes and public export timing remain unknown; F1 R dependencies remain required",
            "minimum_data_needed_to_start": "F0_scRNA_F1_gate PASS or PASS_WITH_NOTED_LIMITATIONS",
            "specialized_audit_required_in_section": "Gate1 preserves the full raw matrices, builds one per-sample min.cells=3 QC working space, applies the fixed five-inequality cell mask, then performs per-sample doublet and retained-cell ambient assessment",
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
    """根据已知数据形态预先限制可用方法，防止把错误输入交给算法。

    例如，没有 raw droplets 就不能声称重做真实 cell calling；小数或标准化
    矩阵也不能伪装成 DESeq2 所需的原始整数 counts。
    """

    return [
        {
            "dataset_or_resource_id": "GSE183904",
            "intended_F_sections": "F1-F4,F8",
            "known_data_structure": "Cell Ranger-derived public called/retained-cell raw gene-count CSV.gz; 26571 archived gene rows by cell columns; F0 builds one per-sample min.cells=3 working feature space solely for the approved QC metrics",
            "known_limitation_from_precheck": "No FASTQ/raw/empty droplets; the exact author export timing, Cell Ranger cell-calling details, DoubletFinder parameters and excluded barcodes are unresolved; PM n=3; Normal_Peritoneum n=1",
            "default_or_recommended_route": "F1 preserves all archived raw rows, recomputes nCount/nFeature/percent.mt/percent.HB after per-sample min.cells=3, applies 500<=nFeature<6000, nCount>1000, percent.mt<=20 and percent.HB<5, then runs per-sample scDblFinder and retained-cell DecontX diagnostics",
            "method_not_allowed_by_default": "No FASTQ/Cell Ranger reproduction claim; no true knee/emptyDrops/SoupX/CellBender without raw droplets; no adaptive MAD/evidence-voting QC or MLMOD/outcome-guided QC; no permanent downstream gene deletion from min.cells=3",
            "required_pre_execution_audit": "F0 data_audit plus processing-history/export-boundary audit and F1 Gate1 QC/doublet/ambient plan",
            "approval_required_before_change": "yes",
            "sensitivity_or_fallback_route": "DoubletFinder same-matrix sensitivity; DecontX corrected-count sensitivity only after its own gate; if fixed QC cannot be recalculated or sample1 regression fails, pause and revise the implementation",
            "interpretation_boundary": "F1 produces a reusable annotated object, not MLMOD conclusions",
            "note": "Do not force F1 to reproduce the paper's 152423-cell object. Every local exclusion receives a barcode-level reason from the fixed QC mask or primary scDblFinder call; ambient scores do not delete cells. Normal_Peritoneum is reference-only and PM statistics are directional only.",
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
    """记录关键研究决策、理由、证据来源、证据强度和敏感性分析要求。

    这样后续审稿或复查时可以区分：哪些是原始数据事实，哪些是文献报告，
    哪些是本项目明确作出的选择。
    """

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
            "decision_value": "Use Cell Ranger-derived public called/retained-cell raw gene-count matrices, not FASTQ/raw droplets; preserve all archived rows and construct one min.cells=3 QC working space per sample",
            "decision_reason": "GEO reports raw files not submitted; supplementary matrices contain nonnegative integer counts. The approved project rule applies min.cells=3 only when computing the four QC metrics and does not treat low-detection archived rows as a processing mismatch.",
            "evidence_type": "original_paper_GEO_metadata_and_stream_audit",
            "evidence_source": "GSE183904 processing-history source audit; GSE183904_series_matrix.txt.gz; data_audit.tsv",
            "source_url_or_file": "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
            "evidence_strength": "high",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "Author SCTransform, integration and clustering are not encoded in the public integer counts; F1 must independently record QC/doublet/ambient decisions and cannot reconstruct author-excluded barcodes.",
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
            "decision_topic": "GSE183904_QC_working_feature_space",
            "decision_value": "Per sample, retain features detected in at least 3 cells before calculating the only nCount/nFeature/percent.mt/percent.HB QC metric set",
            "decision_reason": "The original paper reports considering features detected in at least three cells. Applying that rule once before QC gives an explicit, source-aligned working space while preserving all archived count rows for method-specific downstream filtering.",
            "evidence_type": "original_paper_and_full_stream_sample_by_gene_audit",
            "evidence_source": "F0_author_processing_audit.tsv; data_audit.tsv; Kumar et al. Cancer Discovery 2022",
            "source_url_or_file": "data/metadata/F0_author_processing_audit.tsv",
            "evidence_strength": "high_after_formal_full_stream_audit",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "Low-detection rows remain in the immutable raw input. DE, pseudobulk and scoring later apply their own coverage rules rather than inheriting min.cells=3 permanently.",
        },
        {
            "decision_id": "F0_DECISION_007",
            "stage": "F0/F1",
            "decision_topic": "F1_fixed_QC_framework",
            "decision_value": "Apply 500<=nFeature<6000, nCount>1000, percent.mt<=20 and percent.HB<5 after per-sample min.cells=3; use exact frozen human globin genes",
            "decision_reason": "The nFeature and mitochondrial boundaries are source-reported for this dataset; nCount and globin limits are explicit project additions. A full sample1 pilot confirms that the rule is computable on the actual public matrix, and retaining one deterministic mask keeps every exclusion auditable.",
            "evidence_type": "original_paper_GEO_processing_record_project_decision_and_local_pilot",
            "evidence_source": "Kumar et al. Cancer Discovery 2022; GSE183904 GEO; GSE183904 sample1 full-stream pilot; main plan",
            "source_url_or_file": "胃癌MLMOD亚群主线研究方案.txt",
            "evidence_strength": "high_for_source_thresholds_and_local_computability_project_specific_for_added_thresholds",
            "date": date,
            "requires_sensitivity_analysis": "yes",
            "note": "The fixed mask is not claimed as a universal QC law. Barcode knee, Cell Ranger cell calling, emptyDrops, SoupX and CellBender remain not_evaluable_input_limited; scDblFinder is primary, DoubletFinder is sensitivity, and DecontX scores do not delete cells.",
        },
    ]


def build_excluded_samples(sample_info: Sequence[Dict[str, str]]) -> List[Dict[str, object]]:
    """生成样本层排除记录；即使没有排除，也显式写出 ``no_exclusion``。

    注意：Normal_Peritoneum 默认仍进入对象构建，只是不进入主要组间比较。
    这里登记的是样本去留，不是 Step2 重算得到的细胞去留。
    """

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
    """把 F0 证据归并为恰好十项 gate，并判断能否提交审核。

    每项都包含：要求、实测状态、PASS/FAIL、是否阻断和证据文件。
    ``PASS_WITH_NOTED_ISSUES`` 表示核心契约成立，但仍有公开数据本身无法提供的
    上游信息；这些未知项必须已映射到保守的 F1 处理，不能静默忽略。
    """

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
        """以统一格式向 gate checklist 添加一项判断。"""

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
        """把布尔检查转换成审核表使用的 PASS/FAIL 文本。"""

        return "PASS" if ok else "FAIL"

    # Gate 1：项目目录和追加式日志是否存在。
    structure_ready = (root / "data/metadata/project_structure_ready.txt").exists()
    analysis_log_ready = (root / "logs/F0_setup/analysis_log.md").exists()
    add(
        "project_audit_scaffold",
        "project structure marker and append-only F0 analysis log exist",
        f"structure={'exists' if structure_ready else 'missing'}; log={'exists' if analysis_log_ready else 'missing'}",
        pass_fail(structure_ready and analysis_log_ready),
        "blocking",
        "data/metadata/project_structure_ready.txt; logs/F0_setup/analysis_log.md",
    )
    # Gate 2：40 个提取矩阵是否都有文件大小和合法的大写 SHA256。
    manifest_ok = len(processed_manifest) == 40 and all(
        row.get("sha256")
        and row.get("sha256") == row.get("sha256", "").upper()
        and len(row.get("sha256", "")) == 64
        and all(character in "0123456789ABCDEF" for character in row.get("sha256", ""))
        and row.get("file_size")
        for row in processed_manifest
    )
    add(
        "archive_and_processed_input_manifest",
        "archive is readable as exactly 40 processed csv.gz members, each with file_size and uppercase SHA256",
        f"{len(processed_manifest)} rows; empty_sha={sum(1 for row in processed_manifest if not row.get('sha256'))}",
        pass_fail(manifest_ok),
        "blocking",
        "data/metadata/processed_input_manifest.tsv",
    )
    # Gate 3：40 个样本的文件名、GEO title 和分组是否完全闭合。
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
    # Gate 4：矩阵审计、固定 QC 可计算性和 sample1 冻结回归是否通过。
    included_audits = [row for row in data_audit if row.get("include_in_f1") == "true"]
    pilot_rows = [row for row in data_audit if row.get("file_name") == PILOT_FILE_NAME]
    nonpilot_rows = [row for row in data_audit if row.get("file_name") != PILOT_FILE_NAME]
    pilot_ok = (
        len(pilot_rows) == 1
        and pilot_rows[0].get("pilot_validation_applicable") == "true"
        and pilot_rows[0].get("pilot_validation_status") == "pass"
        and pilot_rows[0].get("qc_retained_feature_count") == "19294"
        and pilot_rows[0].get("source_reported_qc_pass_count") == "2684"
        and pilot_rows[0].get("additional_fail_nCount_after_source_count") == "53"
        and pilot_rows[0].get("additional_fail_percent_hb_after_source_nCount_count") == "0"
        and pilot_rows[0].get("final_fixed_qc_pass_count") == "2631"
        and set(pilot_rows[0].get("globin_panel_used_for_qc", "").split("|"))
        == {"HBA1", "HBA2", "HBB", "HBD"}
    )
    audit_ok = (
        len(data_audit) == 40
        and len(included_audits) == 40
        and pilot_ok
        and all(
            row.get("audit_decision") == "enter_full_F1_independent_reQC"
            and row.get("normalization_artifact_flag") == "false"
            and row.get("observed_numeric_type") == "nonnegative_integer_count_like"
            and row.get("suspected_matrix_type") == "public_called_cell_raw_gene_count_matrix"
            and row.get("matrix_orientation") == "gene_by_cell"
            and row.get("matrix_orientation_validation_status") == "pass_gene_by_cell"
            and row.get("row_label_cell_barcode_like_count") == "0"
            and row.get("cell_barcode_pattern") == "10x_16nt_barcode_with_numeric_suffix"
            and row.get("processed_input_manifest_match") == "true"
            and row.get("per_gene_mean_distribution_status") == "consistent_with_sparse_right_skew"
            and row.get("public_processing_evidence_status")
            == "public_input_shape_verified_fixed_QC_recalculated_processing_history_pending_F0_step3"
            and row.get("working_feature_space_recalculation_status") == "pass"
            and row.get("fixed_qc_rule_recalculation_status") == "pass"
            and row.get("globin_panel_expected") == EXPECTED_GLOBIN_PANEL
            and row.get("globin_panel_present", "") != ""
            and row.get("raw_droplet_available") == "false"
            and row.get("empty_droplet_background_available") == "false"
            and row.get("format_decision_scope") == "file_format_only"
            and row.get("decision_scope") == "file_format_and_public_processing_boundary"
            for row in included_audits
        )
        and all(
            row.get("pilot_validation_applicable") == "false"
            and row.get("pilot_validation_status") == "not_applicable"
            for row in nonpilot_rows
        )
    )
    add(
        "data_audit",
        "40 rows; gene-by-cell orientation is confirmed; included samples pass numeric/distribution/mapping checks; one min.cells=3 QC space and the fixed rule are evaluable; sample1 frozen regression passes",
        f"{len(data_audit)} rows; included={len(included_audits)}; pause={sum(1 for row in data_audit if row.get('audit_decision') != 'enter_full_F1_independent_reQC')}; pilot={'pass' if pilot_ok else 'fail'}",
        pass_fail(audit_ok),
        "blocking",
        "data/metadata/data_audit.tsv",
    )
    # Gate 5：正式完整审计是否与预登记结构事实一致。
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
    # Gate 6：作者声明、F0 实测和未知处理史是否区分清楚且覆盖关键阶段。
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
        "fixed_QC_rule_recalculation",
        "working_feature_space_recalculation",
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
    generated_record_status_ok = all(
        row.get("record_status", "").strip()
        for row in history_rows
        if row.get("history_record_id") in {"H018", "H019", "H020"}
    )
    doublet_row = next(
        (row for row in history_rows if row.get("processing_step") == "doublet_detection"),
        {},
    )
    fixed_qc_row = next(
        (
            row
            for row in history_rows
            if row.get("processing_step") == "fixed_QC_rule_recalculation"
        ),
        {},
    )
    working_feature_row = next(
        (
            row
            for row in history_rows
            if row.get("processing_step") == "working_feature_space_recalculation"
        ),
        {},
    )
    reconciliation_row = next(
        (row for row in history_rows if row.get("processing_step") == "export_boundary_reconciliation"),
        {},
    )
    public_cell_count = sum(int(row.get("matrix_cols_cells", "0") or 0) for row in data_audit)
    source_qc_pass_cells = sum(
        int(row.get("source_reported_qc_pass_count", "0") or 0) for row in data_audit
    )
    fixed_qc_pass_cells = sum(
        int(row.get("final_fixed_qc_pass_count", "0") or 0) for row in data_audit
    )
    fixed_qc_fail_cells = sum(
        int(row.get("final_fixed_qc_fail_count", "0") or 0) for row in data_audit
    )
    fixed_qc_not_evaluable_count = sum(
        row.get("fixed_qc_rule_recalculation_status") != "pass" for row in data_audit
    )
    working_feature_not_evaluable_count = sum(
        row.get("working_feature_space_recalculation_status") != "pass" for row in data_audit
    )
    sample_feature_rows_below_3 = sum(
        int(row.get("feature_rows_detected_lt_3_count", "0") or 0)
        for row in data_audit
    )
    working_feature_counts = [
        int(row.get("qc_retained_feature_count", "0") or 0)
        for row in data_audit
        if row.get("working_feature_space_recalculation_status") == "pass"
    ]
    fixed_status_expected = (
        "F0_fixed_QC_recalculation_pass_all_samples"
        if fixed_qc_not_evaluable_count == 0 and pilot_ok
        else "not_evaluable"
    )
    working_status_expected = (
        "F0_working_feature_space_recalculation_pass_all_samples"
        if working_feature_not_evaluable_count == 0 and len(working_feature_counts) == len(data_audit)
        else "not_evaluable"
    )
    history_ok = (
        required_history_steps.issubset(history_steps)
        and mapped_rows_ok
        and generated_record_status_ok
        and doublet_row.get("author_reported_status") == "reported_in_original_paper"
        and doublet_row.get("relation_to_public_matrix") == "export_boundary_unresolved"
        and fixed_qc_row.get("record_status") == fixed_status_expected
        and str(fixed_qc_row.get("observed_source_reported_qc_pass_cell_count", ""))
        == str(source_qc_pass_cells)
        and str(fixed_qc_row.get("observed_fixed_qc_pass_cell_count", ""))
        == str(fixed_qc_pass_cells)
        and str(fixed_qc_row.get("observed_fixed_qc_fail_cell_count", ""))
        == str(fixed_qc_fail_cells)
        and str(fixed_qc_row.get("observed_samples_fixed_qc_not_evaluable", ""))
        == str(fixed_qc_not_evaluable_count)
        and fixed_qc_row.get("observed_sample1_pilot_validation_status")
        == ("pass" if pilot_ok else "fail")
        and working_feature_row.get("record_status") == working_status_expected
        and str(working_feature_row.get("observed_sample_feature_rows_below_min_cells3", ""))
        == str(sample_feature_rows_below_3)
        and str(working_feature_row.get("observed_working_feature_count_min", ""))
        == str(min(working_feature_counts) if working_feature_counts else "")
        and str(working_feature_row.get("observed_working_feature_count_max", ""))
        == str(max(working_feature_counts) if working_feature_counts else "")
        and str(working_feature_row.get("observed_samples_working_feature_not_evaluable", ""))
        == str(working_feature_not_evaluable_count)
        and str(reconciliation_row.get("observed_public_cell_count", "")) == str(public_cell_count)
        and str(reconciliation_row.get("paper_reported_final_tissue_cell_count", "")) == "152423"
        and str(reconciliation_row.get("count_difference", "")) == str(public_cell_count - 152423)
        and reconciliation_row.get("record_status") == "cross_source_count_difference_observed"
    )
    unresolved_rows = [
        row
        for row in history_rows
        if row.get("relation_to_public_matrix") in {"export_boundary_unresolved", "availability_boundary"}
        or (
            not row.get("record_status", "")
            and row.get("author_reported_status", "").startswith("not_")
        )
        or row.get("author_reported_status") == "raw_files_not_submitted"
    ]
    add(
        "author_processing_provenance",
        "all matrix-relevant stages sourced; author report/F0 verification/inference/unknown separated; unknowns mapped to downstream action",
        f"{len(history_rows)} GSE183904 rows; missing_stages={','.join(sorted(required_history_steps - history_steps)) or 'none'}; "
        f"doublet_status={doublet_row.get('author_reported_status', 'missing')}; "
        f"fixed_QC_recalculation={fixed_qc_row.get('record_status', 'missing')}; "
        f"working_feature_recalculation={working_feature_row.get('record_status', 'missing')}; "
        f"public_vs_paper_cells={reconciliation_row.get('observed_public_cell_count', 'missing')}_vs_"
        f"{reconciliation_row.get('paper_reported_final_tissue_cell_count', 'missing')}",
        (
            "FAIL"
            if not history_ok
            else "PASS_WITH_NOTED_ISSUES"
            if unresolved_rows
            else "PASS"
        ),
        "blocking",
        "docs/source_verification/GSE183904_processing_history_source_audit.tsv; data/metadata/F0_author_processing_audit.tsv",
        f"{len(unresolved_rows)} unresolved or unavailable provenance row(s); unknown history is allowed only when explicit and linked to a conservative F1 action.",
    )
    # Gate 7：固定 QC 在全部样本可计算；不同样本保留数不同本身不算失败。
    add(
        "fixed_QC_rule_recalculation",
        "all 40 samples have an evaluable fixed QC mask; sample1 reproduces the frozen regression",
        f"status={fixed_qc_row.get('record_status', 'missing')}; source_pass={source_qc_pass_cells}; "
        f"final_pass={fixed_qc_pass_cells}; final_fail={fixed_qc_fail_cells}; pilot={'pass' if pilot_ok else 'fail'}",
        pass_fail(fixed_qc_not_evaluable_count == 0 and pilot_ok),
        "blocking",
        "data/metadata/data_audit.tsv; data/metadata/F0_author_processing_audit.tsv",
        "Cell counts are expected to differ by sample; only non-evaluable calculation or pilot disagreement blocks F0.",
    )
    # Gate 8：每个样本都能建立非空的 min.cells=3 QC 工作 feature 空间。
    add(
        "working_feature_space_recalculation",
        "all 40 samples have a nonempty per-sample min.cells=3 QC working feature space",
        f"status={working_feature_row.get('record_status', 'missing')}; "
        f"working_features={min(working_feature_counts) if working_feature_counts else 'NA'}-"
        f"{max(working_feature_counts) if working_feature_counts else 'NA'}; "
        f"archived_sample_gene_rows_below_3={sample_feature_rows_below_3}",
        pass_fail(working_feature_not_evaluable_count == 0 and len(working_feature_counts) == 40),
        "blocking",
        "data/metadata/data_audit.tsv; data/metadata/F0_author_processing_audit.tsv",
        "Low-detection archived rows are preserved and are not classified as an author-processing mismatch.",
    )
    # Gate 9：marker panel 问题只作非阻断 warning，且 F0 不修改 panel。
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
    # Gate 10：全部契约表存在，7 个正式脚本和 2 个 F0 环境锁均有合法大写 SHA256。
    required_contract_files = [
        "data/metadata/F0_dataset_inventory.tsv",
        "data/metadata/F0_file_manifest.tsv",
        "data/metadata/F0_metadata_field_inventory.tsv",
        "data/metadata/F0_author_processing_audit.tsv",
        "data/metadata/F0_data_readiness_by_F_section.tsv",
        "data/metadata/F0_external_resource_inventory.tsv",
        "data/metadata/F0_method_prior_decision.tsv",
        "data/metadata/decision_evidence_log.tsv",
        "data/metadata/excluded_samples.tsv",
    ]
    missing_contract_files = [
        required_file
        for required_file in required_contract_files
        if not (root / required_file).exists()
    ]
    file_manifest_rows = (
        read_tsv(root / "data/metadata/F0_file_manifest.tsv")
        if (root / "data/metadata/F0_file_manifest.tsv").exists()
        else []
    )
    file_manifest_by_path = {
        row.get("relative_path_if_available", ""): row for row in file_manifest_rows
    }
    script_checksum_failures = []
    for script_path in F0_SCRIPT_PATHS:
        manifest_row = file_manifest_by_path.get(script_path, {})
        checksum = manifest_row.get("sha256", "")
        if (
            len(checksum) != 64
            or checksum != checksum.upper()
            or any(character not in "0123456789ABCDEF" for character in checksum)
        ):
            script_checksum_failures.append(script_path)
    environment_lock_checksum_failures = []
    for lock_path in F0_ENVIRONMENT_LOCK_PATHS:
        manifest_row = file_manifest_by_path.get(lock_path, {})
        checksum = manifest_row.get("sha256", "")
        if (
            len(checksum) != 64
            or checksum != checksum.upper()
            or any(character not in "0123456789ABCDEF" for character in checksum)
        ):
            environment_lock_checksum_failures.append(lock_path)
    add(
        "required_F0_contract_tables",
        "all inventory, readiness, method, decision and exclusion tables exist; every F0 script and F0 environment lock has an uppercase SHA256 in F0_file_manifest.tsv",
        f"present={len(required_contract_files) - len(missing_contract_files)}/{len(required_contract_files)}; "
        f"missing={','.join(missing_contract_files) or 'none'}; "
        f"script_checksum_failures={','.join(script_checksum_failures) or 'none'}; "
        f"environment_lock_checksum_failures={','.join(environment_lock_checksum_failures) or 'none'}",
        pass_fail(
            not missing_contract_files
            and not script_checksum_failures
            and not environment_lock_checksum_failures
        ),
        "blocking",
        "; ".join(required_contract_files),
    )
    if len(rows) != 10:
        raise RuntimeError(f"F0 gate checklist contract requires 10 grouped items; observed {len(rows)}")
    return rows


def write_reports(
    root: Path,
    sample_info: Sequence[Dict[str, str]],
    data_audit: Sequence[Dict[str, str]],
    gate_rows: Sequence[Dict[str, object]],
    author_processing: Sequence[Dict[str, str]],
) -> None:
    """把结构化表中的关键结果整理为两份便于人工审核的 Markdown 报告。

    报告只概括表格，不替代原始 TSV。整体 gate 若有阻断失败则为 FAIL；若无
    阻断但存在已登记限制，则为 PASS_WITH_NOTED_LIMITATIONS。
    """

    group_counts: Dict[str, int] = {}
    for row in sample_info:
        group = row.get("group_analysis", "")
        group_counts[group] = group_counts.get(group, 0) + 1
    n_enter = sum(
        1
        for row in data_audit
        if row.get("audit_decision") == "enter_full_F1_independent_reQC"
    )
    enter_files = {
        row.get("file_name", "")
        for row in data_audit
        if row.get("audit_decision") == "enter_full_F1_independent_reQC"
    }
    n_object_eligible = sum(
        row.get("include_in_f1") == "true" and row.get("sample_file", "") in enter_files
        for row in sample_info
    )
    n_excluded_or_pending = len(sample_info) - n_object_eligible
    n_main_group = sum(
        row.get("include_in_f1") == "true"
        and row.get("include_in_group_comparison") == "true"
        and row.get("sample_file", "") in enter_files
        for row in sample_info
    )
    blocking_failed = [
        row for row in gate_rows if row["blocking_level"] == "blocking" and row["pass_fail"] == "FAIL"
    ]
    limited = any(row["pass_fail"] == "PASS_WITH_NOTED_ISSUES" for row in gate_rows)
    gate = "FAIL" if blocking_failed or n_enter != 40 else "PASS_WITH_NOTED_LIMITATIONS" if limited else "PASS"
    run_id = current_run_id()
    gse_history = [row for row in author_processing if row.get("dataset_id") == "GSE183904"]
    doublet_row = next((row for row in gse_history if row.get("processing_step") == "doublet_detection"), {})
    fixed_qc_row = next(
        (
            row
            for row in gse_history
            if row.get("processing_step") == "fixed_QC_rule_recalculation"
        ),
        {},
    )
    working_feature_row = next(
        (
            row
            for row in gse_history
            if row.get("processing_step") == "working_feature_space_recalculation"
        ),
        {},
    )
    reconciliation_row = next((row for row in gse_history if row.get("processing_step") == "export_boundary_reconciliation"), {})
    recon = [
        "# F0 Global Data Reconnaissance Report", "", f"Run ID: {run_id}", f"Generated at: {now_iso()}", "",
        "## Current Usable Data", "",
        f"- GSE183904: 40 Cell Ranger-derived public called/retained-cell raw gene-count CSV.gz matrices; all 26571 archived rows are preserved, and one per-sample min.cells=3 working feature space is used for fixed-QC recalculation; group counts: {group_counts}.",
        "- GSE183904 is the only dataset eligible to start F1 after F0 gate approval.",
        "- GSE239676 passed a structure-only preaudit (8,630 features, 222,240 cells, 20 patients, PC/LM labels) but remains isolated until F2.4 approval.",
        "- Other bulk, multiomics and external resources remain candidates until section-specific audits.", "",
        "## GSE183904 Processing Provenance", "",
        "- Wet-lab acquisition, dissociation, 10x 5-prime library preparation, HiSeq4000 sequencing, Cell Ranger v3.0/hg38 count generation and author Seurat thresholds are source-audited.",
        f"- Fixed-QC recalculation: {fixed_qc_row.get('record_status', 'missing')}; source-reported nFeature/percent.mt rules pass={fixed_qc_row.get('observed_source_reported_qc_pass_cell_count', 'missing')}; final source-plus-project fixed-QC pass={fixed_qc_row.get('observed_fixed_qc_pass_cell_count', 'missing')}; fail={fixed_qc_row.get('observed_fixed_qc_fail_cell_count', 'missing')}; sample1 pilot={fixed_qc_row.get('observed_sample1_pilot_validation_status', 'missing')}.",
        f"- Working feature-space recalculation: {working_feature_row.get('record_status', 'missing')}; retained features per sample={working_feature_row.get('observed_working_feature_count_min', 'missing')}-{working_feature_row.get('observed_working_feature_count_max', 'missing')}; archived sample-by-gene rows below 3 detected cells={working_feature_row.get('observed_sample_feature_rows_below_min_cells3', 'missing')}.",
        f"- Author doublet method: {doublet_row.get('method_or_threshold_if_reported', 'missing')}",
        f"- Public versus paper final tissue-cell count: {reconciliation_row.get('observed_public_cell_count', 'missing')} versus {reconciliation_row.get('paper_reported_final_tissue_cell_count', 'missing')}; difference={reconciliation_row.get('count_difference', 'missing')}.",
        "- Cell Ranger cell-calling details, exact DoubletFinder parameters/barcodes, ambient correction status and the public export timing relative to doublet removal remain unresolved; true knee/cell calling, emptyDrops, SoupX and CellBender are not_evaluable_input_limited.",
        "- Author SCTransform, integration, clustering and annotation describe the author's downstream object and are not embedded in the public integer-count CSV values.", "",
        "## Boundaries", "",
        "- F0 does not produce biological conclusions.",
        "- GSE183904 lacks FASTQ/raw/empty droplets in current public inputs.",
        "- Normal_Peritoneum is reference/display only; PM sample-level inference is directional because n=3.",
        "- The fixed-QC counts describe an independent recalculation on currently public cells; they do not prove when the author filtered cells or reproduce excluded barcodes.",
        "- Any formal/precheck mismatch, non-evaluable fixed rule or sample1 regression failure blocks F1.", "",
        "## F0 Gate Summary", "", f"- F0_scRNA_F1_gate: {gate}",
        "- F0_project_data_inventory_status: partial_with_pending_local_inputs",
        f"- GSE183904 contains {len(sample_info)} samples; {n_object_eligible} are allowed into F1 object construction if approved, {n_excluded_or_pending} are excluded or pending, and {n_main_group} are allowed into the main group comparison.",
        "- Main comparison groups are Normal_Gastric, Primary_Tumor and Peritoneal_Metastasis; Normal_Peritoneum is reference-only.",
        f"- Blocking checklist failures: {len(blocking_failed)}",
    ]
    (root / "results/F0_audit/F0_global_data_reconnaissance_report.md").write_text("\n".join(recon) + "\n", encoding="utf-8")
    report = [
        "# F0 Execution Report", "", f"Run ID: {run_id}", f"Generated at: {now_iso()}", "",
        "## Inputs Checked", "", "- GSE183904 archive, GEO metadata, original-paper processing-history source audit", "- Existing manifests/prechecks",
        "- GSE239676 structure-preaudit records", "- Read-only marker panel plus integrity-only evidence references", "",
        "## Runtime", "", f"- Python {sys.version.split()[0]}; NumPy {np.__version__}; Windows-native proposed execution path.", "",
        "## Main Observations", "", f"- sample_info rows: {len(sample_info)}", f"- data_audit rows: {len(data_audit)}",
        f"- enter_full_F1_independent_reQC samples: {n_enter}", f"- group counts: {group_counts}",
        f"- object-eligible / excluded-or-pending / main-group samples: {n_object_eligible} / {n_excluded_or_pending} / {n_main_group}",
        f"- author processing-history rows for GSE183904: {len(gse_history)}",
        f"- fixed-QC recalculation: {fixed_qc_row.get('record_status', 'missing')}; source-rule pass={fixed_qc_row.get('observed_source_reported_qc_pass_cell_count', 'missing')}; final pass/fail={fixed_qc_row.get('observed_fixed_qc_pass_cell_count', 'missing')}/{fixed_qc_row.get('observed_fixed_qc_fail_cell_count', 'missing')}; sample1 pilot={fixed_qc_row.get('observed_sample1_pilot_validation_status', 'missing')}",
        f"- min.cells=3 working feature space: {working_feature_row.get('record_status', 'missing')}; retained feature range={working_feature_row.get('observed_working_feature_count_min', 'missing')}-{working_feature_row.get('observed_working_feature_count_max', 'missing')}",
        f"- public/paper tissue-cell reconciliation: {reconciliation_row.get('observed_public_cell_count', 'missing')} / {reconciliation_row.get('paper_reported_final_tissue_cell_count', 'missing')}",
        "- unresolved upstream details are carried forward as explicit F1 constraints", "",
        "## Gate Decision", "", f"F0_scRNA_F1_gate: {gate}",
        "F0_project_data_inventory_status: partial_with_pending_local_inputs", "",
        "F1 may start only after Claude Code reviews the executed outputs and the user approves the gate.",
    ]
    (root / "results/F0_audit/F0_execution_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def execute(root: Path) -> int:
    """正式执行 Step4，写出决策表、gate、报告并检查 17 项输出契约。"""

    # 1. 读取前三步已经冻结的输入事实。
    require_paths(root, STAGE_REQUIRED, STAGE_NAME)
    append_log(root, f"F0 step4 started; run_id={current_run_id()}")
    sample_info = read_tsv(root / "data/metadata/sample_info.tsv")
    data_audit = read_tsv(root / "data/metadata/data_audit.tsv")
    processed_manifest = read_tsv(root / "data/metadata/processed_input_manifest.tsv")
    author_processing = read_tsv(root / "data/metadata/F0_author_processing_audit.tsv")

    # 2. 写出“后续能做什么、不能做什么、依据是什么”的四类决策表。
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

    # 3. 生成十项 gate 和两份人工可读报告；这里仍不会启动 F1。
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

    # 4. 报告生成后刷新 manifest，使最终文件状态和脚本 SHA256 都进入清单。
    file_rows = build_file_manifest(root, processed_manifest, F0_OUTPUTS)
    write_tsv(root / "data/metadata/F0_file_manifest.tsv", file_rows, [
        "file_name", "relative_path_if_available", "dataset_id", "file_role", "data_domain", "source_type",
        "source_url_or_local_manifest", "file_size_bytes", "sha256", "compression_format", "read_status",
        "availability_status", "used_in_F0", "planned_F_section_use", "audit_status", "artifact_class",
        "publication_destination", "review_priority", "note",
    ])
    # 5. 最终完整性检查：17 个正式输出必须存在且非空，阻断失败必须报错。
    missing_or_empty_outputs = [
        relative_path
        for relative_path in F0_OUTPUTS
        if not (root / relative_path).exists() or (root / relative_path).stat().st_size == 0
    ]
    if missing_or_empty_outputs:
        append_log(
            root,
            "BLOCKING missing_or_empty_F0_outputs=" + ",".join(missing_or_empty_outputs),
        )
        raise RuntimeError(
            "F0 output contract incomplete: " + ", ".join(missing_or_empty_outputs)
        )
    if failures:
        raise RuntimeError(f"F0 gate contains {failures} blocking failure(s); pause before F1")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """默认 dry run；显式提供 ``--execute`` 后才写决策表和报告。"""

    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, STAGE_REQUIRED, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
