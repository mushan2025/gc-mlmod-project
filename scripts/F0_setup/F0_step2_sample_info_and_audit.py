#!/usr/bin/env python3
"""F0 Step2：建立样本信息并完整审计 40 个 GSE183904 表达矩阵。

本步骤回答两个问题：
1. 每个表达矩阵究竟对应哪个 GEO 样本、患者和组织分组？
2. 公开矩阵是否真的是可用于独立 QC 的非负整数 count 矩阵？

主要输入：Step1 生成的 ``processed_input_manifest.tsv`` 和提取矩阵、GEO
series matrix，以及此前登记的结构预检查表。

主要输出：
- ``sample_info.tsv``：样本、患者、组织来源和是否允许进入 F1；
- ``data_audit.tsv``：每个矩阵的结构、数值、QC 重算和 sample1 回归结果。

实现特点：逐个样本、逐个基因流式读取，不合并 40 个大矩阵。任何未知分组、
样本映射冲突、关键预检查不一致或数值/格式问题都会阻断后续阶段。
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from f0_utils import (
    append_log,
    clean_csv_token,
    current_run_id,
    dry_run_report,
    normalize_sha256,
    parse_stage_args,
    read_tsv,
    require_paths,
    sha256_equal,
    sha256_file,
    validate_f0_python_environment,
    write_tsv,
)


STAGE_NAME = "F0 step 2 sample info and matrix audit"
STAGE_REQUIRED = [
    "data/metadata/processed_input_manifest.tsv",
    "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
    "results/F0_audit/gse183904_csv_structure_precheck.tsv",
    "environment/F0/requirements.txt",
    "environment/F0/environment.yml",
]
STAGE_OUTPUTS = [
    "data/metadata/sample_info.tsv",
    "data/metadata/data_audit.tsv",
]

# 这些正则表达式用于确认矩阵内容和命名格式，而不是进行生物学筛选。
NONNEGATIVE_INTEGER_ROW = re.compile(r"\d+(?:,\d+)*\Z")
GENE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
CELL_BARCODE = re.compile(r"[ACGTN]{16}_\d+\Z", re.IGNORECASE)
SAMPLE_ID = re.compile(r"sample\d+\Z", re.IGNORECASE)

# 标准化伪影哨兵在查看新算出的细胞 nCount 前已经预登记。
# 常见 CP1k/CP10k/CP100k/CPM 会把每个细胞的总量强行缩放到相近固定值；
# 因此若绝大多数细胞 nCount 异常集中在这些值附近，就要怀疑输入并非原始 counts。
# 这些阈值设置得很严格：触发只代表“必须复核”，不能单独证明作者使用了哪种方法。
FIXED_LIBRARY_SIZE_TARGETS = (1_000, 10_000, 100_000, 1_000_000)
FIXED_TARGET_RELATIVE_TOLERANCE = 0.01
FIXED_TARGET_MIN_CELL_FRACTION = 0.90
FIXED_TARGET_MAX_RELATIVE_IQR = 0.02
DOMINANT_ROUND_MULTIPLE = 100
DOMINANT_ROUND_TOP_N = 10
DOMINANT_ROUND_MIN_CELL_FRACTION = 0.80
NEAR_CONSTANT_MAX_RELATIVE_IQR = 0.01
NEAR_CONSTANT_MAX_RELATIVE_RANGE = 0.05
CONCENTRATION_RULE_MIN_CELLS = 200
LOW_NCOUNT_REVIEW_BOUNDARY = 500
HIGH_NCOUNT_MEDIAN_REVIEW_BOUNDARY = 100_000
HIGH_NCOUNT_MAX_REVIEW_BOUNDARY = 1_000_000
SPARSE_MATRIX_MIN_ZERO_RATE = 0.50
# 冻结 QC 阈值。变量名中的 EXCLUSIVE 表示该边界本身不通过。
SOURCE_MIN_NFEATURE = 500
SOURCE_MAX_NFEATURE_EXCLUSIVE = 6000
SOURCE_MAX_PERCENT_MT = 20.0
PROJECT_MIN_NCOUNT_EXCLUSIVE = 1000
PROJECT_MAX_PERCENT_HB_EXCLUSIVE = 5.0
QC_MIN_CELLS_PER_FEATURE = 3

# 获批 F0/F1 方案冻结的人 globin 转录本 panel。
# 基因名先转成大写再精确匹配；宽泛的 ``^HB`` 会把 HBEGF、HBS1L、HBP1
# 等非 globin 基因误计入红细胞比例，因此明确禁止。
GLOBIN_PANEL = (
    "HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "HBE1",
    "HBG1",
    "HBG2",
    "HBM",
    "HBQ1",
    "HBZ",
)
GLOBIN_PANEL_SET = frozenset(GLOBIN_PANEL)

# sample1 的结果已在正式执行前冻结，用作“同一代码、同一数据是否仍给出同一结果”的回归检查。
PILOT_FILE_NAME = "GSM5573466_sample1.csv.gz"
PILOT_EXPECTED = {
    "matrix_cols_cells": 2685,
    "row_label_cell_barcode_like_count": 0,
    "qc_retained_feature_count": 19294,
    "source_reported_qc_pass_count": 2684,
    "additional_fail_nCount_after_source_count": 53,
    "additional_fail_percent_hb_after_source_nCount_count": 0,
    "final_fixed_qc_pass_count": 2631,
}
PILOT_EXPECTED_GLOBIN_GENES_USED = frozenset({"HBA1", "HBA2", "HBB", "HBD"})


def parse_geo_line(line: str) -> Tuple[str, List[str]]:
    """解析 GEO series matrix 的一行，去除字段外围引号。"""

    parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
    return parts[0], [value.strip().strip('"') for value in parts[1:]]


def parse_gse183904_series(series_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """读取 GEO 元数据，并提取样本字段和 patient_id 对应关系。

    GEO 文件中患者与 sampleN 的映射藏在 Series_summary 文本中，因此这里
    单独解析；无法确认的患者不会被猜测，而会在 sample_info 中记为 unknown。
    """

    sample_fields: Dict[str, List[str]] = {}
    patient_by_sample: Dict[str, str] = {}
    mapping_header: List[str] | None = None
    in_mapping = False

    with gzip.open(series_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
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
                                patient_by_sample[sample_file.replace(".csv", "")] = patient
                    elif text.startswith("!"):
                        in_mapping = False

    required = ["!Sample_geo_accession", "!Sample_title", "!Sample_source_name_ch1"]
    missing = [key for key in required if key not in sample_fields]
    if missing:
        raise RuntimeError(f"GSE183904 series matrix missing fields: {', '.join(missing)}")
    return sample_fields, patient_by_sample


def group_from_title(title: str) -> Tuple[str, str, str, str, str]:
    """按预登记关键词把 GEO title 映射为四种组织分组。

    没有命中规则时返回 ``Unclear``，由 gate 阻断，而不是根据文件顺序或主观
    判断补分组。
    """

    text = title.lower()
    if "primary gastric tissue" in text and "normal" in text:
        return (
            "Normal_Gastric",
            "gastric",
            "non_tumor",
            "none",
            'GEO title contains "Primary Gastric Tissue" and "Normal" -> Normal_Gastric',
        )
    if "primary gastric tissue" in text and "tumor" in text:
        return (
            "Primary_Tumor",
            "gastric",
            "tumor",
            "none",
            'GEO title contains "Primary Gastric Tissue" and "Tumor" -> Primary_Tumor',
        )
    if "periton" in text and "normal" in text:
        return (
            "Normal_Peritoneum",
            "peritoneum",
            "non_tumor",
            "none",
            'GEO title contains "Peritonium tissue" and "Normal" -> Normal_Peritoneum',
        )
    if "periton" in text and "tumor" in text:
        return (
            "Peritoneal_Metastasis",
            "peritoneum",
            "tumor",
            "peritoneum",
            'GEO title contains "Peritonium tissue" and "Tumor" -> Peritoneal_Metastasis',
        )
    return "Unclear", "unclear", "unclear", "unclear", "no preregistered title rule matched"


def build_sample_info(
    sample_fields: Dict[str, List[str]],
    patient_by_sample: Dict[str, str],
    manifest_rows: Sequence[Dict[str, str]],
) -> Tuple[List[Dict[str, object]], List[str]]:
    """交叉核对 GEO accession、title 中的 sampleN 和压缩包文件名。

    输出每个样本的分组、patient_id、纳入状态和所有已知限制。Normal_Peritoneum
    被保留为参考，但不进入主要组间比较；PM 因样本数仅 3，后续只能方向性解释。
    """

    accessions = sample_fields["!Sample_geo_accession"]
    titles = sample_fields["!Sample_title"]
    sources = sample_fields.get("!Sample_source_name_ch1", [""] * len(accessions))
    chars = sample_fields.get("!Sample_characteristics_ch1", [""] * len(accessions))
    manifest_by_accession: Dict[str, Dict[str, str]] = {}
    for row in manifest_rows:
        accession = row.get("geo_accession", "")
        if not accession or accession in manifest_by_accession:
            raise RuntimeError(f"Duplicate or empty GEO accession in processed manifest: {accession!r}")
        manifest_by_accession[accession] = row

    if len(set(accessions)) != len(accessions):
        raise RuntimeError("GSE183904 series matrix contains duplicate GEO accessions")

    rows: List[Dict[str, object]] = []
    mismatch_notes: List[str] = []
    for idx, accession in enumerate(accessions):
        title = titles[idx]
        sample_match = re.search(r"\b(sample\d+)\b", title, flags=re.IGNORECASE)
        geo_title_sample_id = sample_match.group(1).lower() if sample_match else ""
        manifest_row = manifest_by_accession.get(accession, {})
        sample_id = manifest_row.get("sample_id", "")
        member_name = manifest_row.get("archive_member_name", "")
        if not manifest_row:
            sample_id_match_status = "manifest_accession_missing"
            mismatch_notes.append(f"{accession}: GEO accession is absent from processed manifest")
        elif not geo_title_sample_id:
            sample_id_match_status = "geo_title_sample_id_missing"
            mismatch_notes.append(f"{accession}: GEO title lacks sampleN token; title={title!r}")
        elif not SAMPLE_ID.fullmatch(sample_id):
            sample_id_match_status = "manifest_sample_id_invalid"
            mismatch_notes.append(
                f"{accession}: manifest sample_id has unexpected format; manifest={sample_id!r}"
            )
        elif geo_title_sample_id != sample_id.lower():
            sample_id_match_status = "mismatch"
            mismatch_notes.append(
                f"{accession}: sample_id mismatch; manifest={sample_id}; GEO_title={geo_title_sample_id}"
            )
        else:
            sample_id_match_status = "match"

        group, tissue_site, tumor_status, metastasis_site, group_rule = group_from_title(title)
        group_original = title.split(":", 1)[1].strip() if ":" in title else title
        patient_id = patient_by_sample.get(sample_id, "unknown")
        metadata_issues: List[str] = []
        if patient_id == "unknown":
            metadata_issues.append("patient_id_missing")
        metadata_issues.append("pairing_unknown")
        if "peritonium" in group_original.lower():
            metadata_issues.append("geo_typo_peritonium")
        if "  " in group_original:
            metadata_issues.append("double_space_in_group_original")
        if sample_id_match_status != "match":
            metadata_issues.append(f"sample_id_{sample_id_match_status}")

        include_in_f1 = "true" if group != "Unclear" and sample_id_match_status == "match" else "pending"
        include_in_group = (
            "false"
            if group == "Normal_Peritoneum"
            else "true"
            if include_in_f1 == "true"
            else "pending"
        )
        if group == "Normal_Peritoneum":
            include_reason = "single_normal_peritoneum_reference_only"
        elif include_in_f1 == "pending":
            include_reason = "metadata_or_sample_file_mapping_requires_review"
        else:
            include_reason = "preregistered_GSE183904_group_and_file_mapping_confirmed"
        rows.append(
            {
                "dataset_id": "GSE183904",
                "sample_file": member_name,
                "geo_accession": accession,
                "sample_id": sample_id,
                "geo_title_sample_id": geo_title_sample_id,
                "sample_id_match_status": sample_id_match_status,
                "patient_id": patient_id,
                "group_original": group_original,
                "group_analysis": group,
                "group_analysis_rule": group_rule,
                "sample_title": title,
                "source_name_ch1": sources[idx] if idx < len(sources) else "",
                "sample_characteristics_ch1": chars[idx] if idx < len(chars) else "",
                "tissue_site": tissue_site,
                "tumor_status": tumor_status,
                "metastasis_site": metastasis_site,
                "is_paired": "unknown",
                "paired_normal_id": "unknown",
                "include_in_f1": include_in_f1,
                "include_in_group_comparison": include_in_group,
                "include_reason": include_reason,
                "source_of_group": "GEO_Sample_title",
                "metadata_confidence": (
                    "high" if group != "Unclear" and sample_id_match_status == "match" else "low"
                ),
                "metadata_issue": ";".join(metadata_issues) if metadata_issues else "none",
                "pairing_scope": "patient_id_from_GEO_title_mapping_available_but_F1_pairing_not_preapproved",
                "note": (
                    "Normal_Peritoneum is reference/display only"
                    if group == "Normal_Peritoneum"
                    else "PM sample-level comparisons are directional only because n=3"
                    if group == "Peritoneal_Metastasis"
                    else ""
                ),
            }
        )
    missing_from_geo = sorted(set(manifest_by_accession) - set(accessions))
    for accession in missing_from_geo:
        mismatch_notes.append(f"{accession}: processed manifest accession is absent from GEO series matrix")
    return rows, mismatch_notes


def format_number(value: float | int | None) -> str:
    """把数值转换为稳定、紧凑的文本格式，供 TSV 输出使用。"""

    if value is None or not math.isfinite(float(value)):
        return ""
    numeric = float(value)
    if numeric.is_integer():
        return str(int(numeric))
    return format(numeric, ".12g")


def classify_anomalous_values(values: str, expected_cells: int) -> Dict[str, object]:
    """在快速整数解析失败时，逐值区分缺失、小数、非数值和负整数。

    公开 count 矩阵应由非负整数构成。异常分类的目的不是修补数据，而是提供
    可定位的失败原因并暂停执行。
    """

    tokens = values.split(",")
    counts: Dict[str, object] = {
        "observed_values": len(tokens),
        "integer_count": 0,
        "missing_value_count": 0,
        "noninteger_float_count": 0,
        "invalid_nonnumeric_count": 0,
        "negative_integer_count": 0,
        "zero_value_count": 0,
        "numeric_min": None,
        "numeric_max": None,
    }
    integer_values: List[int] = []
    all_integer = len(tokens) == expected_cells
    for raw_value in tokens:
        value = clean_csv_token(raw_value)
        if value == "":
            counts["missing_value_count"] = int(counts["missing_value_count"]) + 1
            all_integer = False
            continue
        try:
            integer = int(value)
        except ValueError:
            all_integer = False
            try:
                numeric = float(value)
            except ValueError:
                counts["invalid_nonnumeric_count"] = int(counts["invalid_nonnumeric_count"]) + 1
                continue
            if not math.isfinite(numeric):
                counts["invalid_nonnumeric_count"] = int(counts["invalid_nonnumeric_count"]) + 1
                continue
            counts["noninteger_float_count"] = int(counts["noninteger_float_count"]) + 1
            if numeric == 0:
                counts["zero_value_count"] = int(counts["zero_value_count"]) + 1
        else:
            integer_values.append(integer)
            counts["integer_count"] = int(counts["integer_count"]) + 1
            numeric = float(integer)
            if integer < 0:
                counts["negative_integer_count"] = int(counts["negative_integer_count"]) + 1
            if integer == 0:
                counts["zero_value_count"] = int(counts["zero_value_count"]) + 1

        current_min = counts["numeric_min"]
        current_max = counts["numeric_max"]
        counts["numeric_min"] = numeric if current_min is None else min(float(current_min), numeric)
        counts["numeric_max"] = numeric if current_max is None else max(float(current_max), numeric)

    counts["integer_values"] = integer_values if all_integer else None
    return counts


def assess_ncount_distribution(ncount: np.ndarray, complete: bool) -> Dict[str, object]:
    """检查基于全部存档行的细胞总 count 分布是否像原始 counts。

    这里的 ``raw_full_nCount`` 只用于识别固定库大小等标准化伪影，不参与细胞
    去留。若大量细胞总量几乎相同或集中在 1,000/10,000 等固定值，脚本会
    标记复核，而不会擅自断言具体标准化方法。
    """

    blank = {
        "raw_full_nCount_min": "",
        "raw_full_nCount_Q1": "",
        "raw_full_nCount_median": "",
        "raw_full_nCount_Q3": "",
        "raw_full_nCount_max": "",
        "ncount_distinct_count": "",
        "ncount_relative_iqr": "",
        "ncount_relative_range": "",
        "dominant_round_ncount_fraction": "",
        "fixed_target_near_value": "",
        "fixed_target_near_fraction": "",
        "ncount_range_status": "not_evaluable",
        "normalization_artifact_flag": "not_evaluable",
        "normalization_artifact_reason": "nCount_not_evaluable_due_to_column_or_numeric_anomaly",
    }
    if not complete or ncount.size == 0:
        return blank

    q1, median, q3 = np.quantile(ncount, [0.25, 0.5, 0.75], method="linear")
    minimum = int(np.min(ncount))
    maximum = int(np.max(ncount))
    distinct_values, frequencies = np.unique(ncount, return_counts=True)
    distinct_count = int(distinct_values.size)
    relative_iqr = float((q3 - q1) / median) if median > 0 else math.inf
    relative_range = float((maximum - minimum) / median) if median > 0 else math.inf

    frequency_order = np.argsort(frequencies)[::-1][:DOMINANT_ROUND_TOP_N]
    dominant_round_cells = sum(
        int(frequencies[idx])
        for idx in frequency_order
        if int(distinct_values[idx]) % DOMINANT_ROUND_MULTIPLE == 0
    )
    dominant_round_fraction = dominant_round_cells / int(ncount.size)

    target_fractions = {
        target: float(
            np.mean(
                np.abs(ncount.astype(np.float64) - target)
                <= max(1.0, target * FIXED_TARGET_RELATIVE_TOLERANCE)
            )
        )
        for target in FIXED_LIBRARY_SIZE_TARGETS
    }
    best_target, best_target_fraction = max(target_fractions.items(), key=lambda item: item[1])

    triggers: List[str] = []
    if distinct_count == 1:
        triggers.append("all_cell_nCount_identical")
    if (
        ncount.size >= CONCENTRATION_RULE_MIN_CELLS
        and dominant_round_fraction >= DOMINANT_ROUND_MIN_CELL_FRACTION
    ):
        triggers.append("at_least_80pct_cells_in_top10_round100_totals")
    if (
        ncount.size >= CONCENTRATION_RULE_MIN_CELLS
        and best_target_fraction >= FIXED_TARGET_MIN_CELL_FRACTION
        and relative_iqr <= FIXED_TARGET_MAX_RELATIVE_IQR
    ):
        triggers.append(f"fixed_library_size_target_{best_target}_concentration")
    if (
        ncount.size >= CONCENTRATION_RULE_MIN_CELLS
        and relative_iqr <= NEAR_CONSTANT_MAX_RELATIVE_IQR
        and relative_range <= NEAR_CONSTANT_MAX_RELATIVE_RANGE
    ):
        triggers.append("near_constant_nonround_library_size")

    range_warnings: List[str] = []
    if minimum < LOW_NCOUNT_REVIEW_BOUNDARY:
        range_warnings.append(f"min_below_{LOW_NCOUNT_REVIEW_BOUNDARY}")
    if median > HIGH_NCOUNT_MEDIAN_REVIEW_BOUNDARY:
        range_warnings.append(f"median_above_{HIGH_NCOUNT_MEDIAN_REVIEW_BOUNDARY}")
    if maximum > HIGH_NCOUNT_MAX_REVIEW_BOUNDARY:
        range_warnings.append(f"max_above_{HIGH_NCOUNT_MAX_REVIEW_BOUNDARY}")

    trigger_text = ",".join(triggers) if triggers else "none"
    reason = (
        f"trigger={trigger_text}; distinct={distinct_count}; relative_iqr={relative_iqr:.6g}; "
        f"relative_range={relative_range:.6g}; dominant_round_fraction={dominant_round_fraction:.6g}; "
        f"best_fixed_target={best_target}; best_fixed_target_fraction={best_target_fraction:.6g}"
    )
    return {
        "raw_full_nCount_min": format_number(minimum),
        "raw_full_nCount_Q1": format_number(q1),
        "raw_full_nCount_median": format_number(median),
        "raw_full_nCount_Q3": format_number(q3),
        "raw_full_nCount_max": format_number(maximum),
        "ncount_distinct_count": distinct_count,
        "ncount_relative_iqr": format_number(relative_iqr),
        "ncount_relative_range": format_number(relative_range),
        "dominant_round_ncount_fraction": format(dominant_round_fraction, ".12g"),
        "fixed_target_near_value": best_target,
        "fixed_target_near_fraction": format(best_target_fraction, ".12g"),
        "ncount_range_status": (
            "within_preregistered_review_boundaries"
            if not range_warnings
            else "review_required:" + ",".join(range_warnings)
        ),
        "normalization_artifact_flag": "true" if triggers else "false",
        "normalization_artifact_reason": reason,
    }


def assess_fixed_qc(
    ncount: np.ndarray,
    nfeature: np.ndarray,
    mt_ncount: np.ndarray,
    hb_ncount: np.ndarray,
    complete: bool,
) -> Dict[str, object]:
    """在唯一的 ``min.cells=3`` 工作空间中重算固定细胞 QC 规则。

    四项指标含义：
    - nCount：细胞总 UMI/count，过低提示信息量不足；
    - nFeature：检出的基因数，过低常见于低质量细胞，过高可能提示混合细胞；
    - percent.mt：线粒体转录本比例，过高常提示细胞受损；
    - percent.HB：globin 转录本比例，过高提示红细胞/血液成分干扰。

    本函数只返回各规则和最终联合规则的计数，不在 F0 中真正删除细胞。
    """

    metric_fields = {
        "qc_nCount_min": "",
        "qc_nCount_Q1": "",
        "qc_nCount_median": "",
        "qc_nCount_Q3": "",
        "qc_nCount_max": "",
        "qc_nFeature_min": "",
        "qc_nFeature_Q1": "",
        "qc_nFeature_median": "",
        "qc_nFeature_Q3": "",
        "qc_nFeature_max": "",
        "qc_percent_mt_min": "",
        "qc_percent_mt_Q1": "",
        "qc_percent_mt_median": "",
        "qc_percent_mt_Q3": "",
        "qc_percent_mt_max": "",
        "qc_percent_hb_min": "",
        "qc_percent_hb_Q1": "",
        "qc_percent_hb_median": "",
        "qc_percent_hb_Q3": "",
        "qc_percent_hb_max": "",
        "fail_nFeature_low_count": "",
        "fail_nFeature_high_count": "",
        "fail_nCount_count": "",
        "fail_percent_mt_count": "",
        "fail_percent_hb_count": "",
        "source_reported_qc_fail_count": "",
        "source_reported_qc_pass_count": "",
        "additional_fail_nCount_after_source_count": "",
        "additional_fail_percent_hb_after_source_nCount_count": "",
        "final_fixed_qc_fail_count": "",
        "final_fixed_qc_pass_count": "",
        "fixed_qc_rule_recalculation_status": "not_evaluable",
        "fixed_qc_rule_recalculation_note": (
            "fixed_QC_not_evaluable_due_to_column_numeric_or_denominator_anomaly"
        ),
    }
    if (
        not complete
        or ncount.size == 0
        or nfeature.size != ncount.size
        or mt_ncount.size != ncount.size
        or hb_ncount.size != ncount.size
        or np.any(ncount <= 0)
    ):
        return metric_fields

    # 百分比的分母和分子都来自同一个 min.cells=3 工作 feature 空间。
    percent_mt = np.divide(
        mt_ncount.astype(np.float64) * 100.0,
        ncount.astype(np.float64),
    )
    percent_hb = np.divide(
        hb_ncount.astype(np.float64) * 100.0,
        ncount.astype(np.float64),
    )
    # 明确写出五个失败布尔数组，便于核对等号究竟属于保留还是排除。
    feature_low = nfeature < SOURCE_MIN_NFEATURE
    feature_high = nfeature >= SOURCE_MAX_NFEATURE_EXCLUSIVE
    count_low = ncount <= PROJECT_MIN_NCOUNT_EXCLUSIVE
    mt_high = percent_mt > SOURCE_MAX_PERCENT_MT
    hb_high = percent_hb >= PROJECT_MAX_PERCENT_HB_EXCLUSIVE
    source_fail = feature_low | feature_high | mt_high
    source_pass = ~source_fail
    after_source_and_ncount = source_pass & ~count_low
    # “或”表示违反任一规则即不通过；联合计数不会重复计算同一个细胞。
    final_fail = source_fail | count_low | hb_high
    ncount_q1, ncount_median, ncount_q3 = np.quantile(
        ncount, [0.25, 0.5, 0.75], method="linear"
    )
    nfeature_q1, nfeature_median, nfeature_q3 = np.quantile(
        nfeature, [0.25, 0.5, 0.75], method="linear"
    )
    mt_q1, mt_median, mt_q3 = np.quantile(
        percent_mt, [0.25, 0.5, 0.75], method="linear"
    )
    hb_q1, hb_median, hb_q3 = np.quantile(
        percent_hb, [0.25, 0.5, 0.75], method="linear"
    )
    low_count = int(np.count_nonzero(feature_low))
    high_count = int(np.count_nonzero(feature_high))
    ncount_low_count = int(np.count_nonzero(count_low))
    mt_high_count = int(np.count_nonzero(mt_high))
    hb_high_count = int(np.count_nonzero(hb_high))
    source_fail_count = int(np.count_nonzero(source_fail))
    source_pass_count = int(ncount.size) - source_fail_count
    additional_ncount = int(np.count_nonzero(source_pass & count_low))
    additional_hb = int(np.count_nonzero(after_source_and_ncount & hb_high))
    final_fail_count = int(np.count_nonzero(final_fail))
    final_pass_count = int(ncount.size) - final_fail_count
    note = (
        "fixed_rule=500<=nFeature<6000,nCount>1000,percent.mt<=20,percent.HB<5; "
        "source_rules=nFeature_and_percent.mt; project_rules=nCount_and_percent.HB; "
        f"source_pass={source_pass_count}; additional_nCount_fail={additional_ncount}; "
        f"additional_percent.HB_fail={additional_hb}; final_pass={final_pass_count}"
    )
    return {
        "qc_nCount_min": format_number(int(np.min(ncount))),
        "qc_nCount_Q1": format_number(ncount_q1),
        "qc_nCount_median": format_number(ncount_median),
        "qc_nCount_Q3": format_number(ncount_q3),
        "qc_nCount_max": format_number(int(np.max(ncount))),
        "qc_nFeature_min": format_number(int(np.min(nfeature))),
        "qc_nFeature_Q1": format_number(nfeature_q1),
        "qc_nFeature_median": format_number(nfeature_median),
        "qc_nFeature_Q3": format_number(nfeature_q3),
        "qc_nFeature_max": format_number(int(np.max(nfeature))),
        "qc_percent_mt_min": format_number(float(np.min(percent_mt))),
        "qc_percent_mt_Q1": format_number(mt_q1),
        "qc_percent_mt_median": format_number(mt_median),
        "qc_percent_mt_Q3": format_number(mt_q3),
        "qc_percent_mt_max": format_number(float(np.max(percent_mt))),
        "qc_percent_hb_min": format_number(float(np.min(percent_hb))),
        "qc_percent_hb_Q1": format_number(hb_q1),
        "qc_percent_hb_median": format_number(hb_median),
        "qc_percent_hb_Q3": format_number(hb_q3),
        "qc_percent_hb_max": format_number(float(np.max(percent_hb))),
        "fail_nFeature_low_count": low_count,
        "fail_nFeature_high_count": high_count,
        "fail_nCount_count": ncount_low_count,
        "fail_percent_mt_count": mt_high_count,
        "fail_percent_hb_count": hb_high_count,
        "source_reported_qc_fail_count": source_fail_count,
        "source_reported_qc_pass_count": source_pass_count,
        "additional_fail_nCount_after_source_count": additional_ncount,
        "additional_fail_percent_hb_after_source_nCount_count": additional_hb,
        "final_fixed_qc_fail_count": final_fail_count,
        "final_fixed_qc_pass_count": final_pass_count,
        "fixed_qc_rule_recalculation_status": "pass",
        "fixed_qc_rule_recalculation_note": note,
    }


def assess_working_feature_space(
    detected_in_0_cells: int,
    detected_in_1_cell: int,
    detected_in_2_cells: int,
    examples: Sequence[str],
    total_features: int,
    complete: bool,
) -> Dict[str, object]:
    """汇总每个样本中低于 ``min.cells=3`` 的 feature 数和最终工作维度。

    这些低检出 feature 只是不参与细胞 QC 指标计算；原始矩阵行仍完整保留，
    也不能据此声称作者导出数据时没有做过基因过滤。
    """

    below_threshold = detected_in_0_cells + detected_in_1_cell + detected_in_2_cells
    if not complete:
        return {
            "qc_min_cells_per_feature": QC_MIN_CELLS_PER_FEATURE,
            "feature_rows_detected_in_0_cells": "",
            "feature_rows_detected_in_1_cell": "",
            "feature_rows_detected_in_2_cells": "",
            "feature_rows_detected_lt_3_count": "",
            "feature_rows_detected_lt_3_examples": "",
            "qc_retained_feature_count": "",
            "working_feature_space_recalculation_status": "not_evaluable",
            "working_feature_space_recalculation_note": (
                "per_sample_feature_detection_not_evaluable_due_to_column_or_numeric_anomaly"
            ),
        }
    return {
        "qc_min_cells_per_feature": QC_MIN_CELLS_PER_FEATURE,
        "feature_rows_detected_in_0_cells": detected_in_0_cells,
        "feature_rows_detected_in_1_cell": detected_in_1_cell,
        "feature_rows_detected_in_2_cells": detected_in_2_cells,
        "feature_rows_detected_lt_3_count": below_threshold,
        "feature_rows_detected_lt_3_examples": "|".join(examples),
        "qc_retained_feature_count": total_features - below_threshold,
        "working_feature_space_recalculation_status": "pass",
        "working_feature_space_recalculation_note": (
            "working_rule=feature_detected_in_at_least_3_cells; "
            f"detected_in_0={detected_in_0_cells}; detected_in_1={detected_in_1_cell}; "
            f"detected_in_2={detected_in_2_cells}; below_3={below_threshold}"
        ),
    }


def summarize_gene_means(
    gene_means: Sequence[float],
    matrix_rows: int,
    zero_value_count: int,
    total_values_checked: int,
) -> Dict[str, object]:
    """用零值率和基因均值分布确认矩阵具有单细胞 count 的稀疏右偏特征。"""

    zero_rate = zero_value_count / total_values_checked if total_values_checked else math.nan
    if len(gene_means) != matrix_rows or not gene_means:
        return {
            "zero_value_count": zero_value_count,
            "zero_value_rate": format(zero_rate, ".12g") if math.isfinite(zero_rate) else "",
            "per_gene_mean_distribution_status": "not_evaluable",
            "per_gene_mean_distribution_note": (
                f"not_evaluable; complete_gene_means={len(gene_means)}; matrix_rows={matrix_rows}; "
                f"zero_value_rate={format(zero_rate, '.12g') if math.isfinite(zero_rate) else 'NA'}"
            ),
        }

    means = np.asarray(gene_means, dtype=np.float64)
    q1, median, q3 = np.quantile(means, [0.25, 0.5, 0.75], method="linear")
    mean_of_means = float(np.mean(means))
    maximum = float(np.max(means))
    low_001 = float(np.mean(means <= 0.01))
    low_01 = float(np.mean(means <= 0.1))
    tail_ratio = maximum / median if median > 0 else math.inf
    status = (
        "consistent_with_sparse_right_skew"
        if zero_rate >= SPARSE_MATRIX_MIN_ZERO_RATE and maximum > q3 and mean_of_means > median
        else "review_required"
    )
    note = (
        f"status={status}; zero_value_rate={zero_rate:.6g}; gene_mean_le_0.01_fraction={low_001:.6g}; "
        f"gene_mean_le_0.1_fraction={low_01:.6g}; gene_mean_Q1={q1:.6g}; "
        f"gene_mean_median={median:.6g}; gene_mean_Q3={q3:.6g}; gene_mean_mean={mean_of_means:.6g}; "
        f"gene_mean_max={maximum:.6g}; max_to_median={tail_ratio:.6g}"
    )
    return {
        "zero_value_count": zero_value_count,
        "zero_value_rate": format(zero_rate, ".12g"),
        "per_gene_mean_distribution_status": status,
        "per_gene_mean_distribution_note": note,
    }


def validate_frozen_pilot(
    file_name: str,
    observed: Dict[str, object],
    globin_genes_used: Sequence[str],
) -> Dict[str, object]:
    """仅对指定 sample1 检查冻结数字；其他样本标记为不适用。

    该回归检查用于发现代码或输入发生意外变化，不要求其他样本复制 sample1
    的细胞数或基因数。
    """

    if file_name != PILOT_FILE_NAME:
        return {
            "pilot_validation_applicable": "false",
            "pilot_validation_status": "not_applicable",
            "pilot_validation_note": "Frozen regression target applies only to GSM5573466_sample1.csv.gz.",
        }

    checks: List[str] = []
    passed = True
    for field, expected in PILOT_EXPECTED.items():
        raw_observed = observed.get(field, "")
        try:
            numeric_observed = int(raw_observed)
        except (TypeError, ValueError):
            numeric_observed = None
        field_passed = numeric_observed == expected
        passed = passed and field_passed
        checks.append(
            f"{field}={raw_observed}(expected={expected},status={'pass' if field_passed else 'fail'})"
        )

    observed_globin = frozenset(gene.upper() for gene in globin_genes_used)
    globin_passed = observed_globin == PILOT_EXPECTED_GLOBIN_GENES_USED
    passed = passed and globin_passed
    checks.append(
        "globin_panel_used_for_qc="
        + "|".join(sorted(observed_globin))
        + "(expected="
        + "|".join(sorted(PILOT_EXPECTED_GLOBIN_GENES_USED))
        + f",status={'pass' if globin_passed else 'fail'})"
    )
    orientation_passed = (
        observed.get("matrix_orientation") == "gene_by_cell"
        and observed.get("matrix_orientation_validation_status") == "pass_gene_by_cell"
    )
    passed = passed and orientation_passed
    checks.append(
        "matrix_orientation="
        f"{observed.get('matrix_orientation', '')}; "
        "validation_status="
        f"{observed.get('matrix_orientation_validation_status', '')}"
        f"(expected=gene_by_cell/pass_gene_by_cell,status={'pass' if orientation_passed else 'fail'})"
    )
    return {
        "pilot_validation_applicable": "true",
        "pilot_validation_status": "pass" if passed else "fail",
        "pilot_validation_note": "; ".join(checks),
    }


def audit_csv_gz(path: Path) -> Dict[str, object]:
    """完整流式审计一个“基因 × 细胞”的压缩 CSV 矩阵。

    每次只把一个基因的整行 count 读入内存，因此内存主要随单个样本的细胞数
    增长，而不是随 40 个样本的总矩阵增长。函数同时完成结构检查、数值检查、
    工作 feature 空间构建、QC 重算和可复现性摘要。
    """

    start = time.perf_counter()
    rows = 0
    bad_column_count_rows = 0
    missing_value_count = 0
    noninteger_float_count = 0
    invalid_nonnumeric_count = 0
    negative_integer_count = 0
    integer_count = 0
    zero_value_count = 0
    total_values_checked = 0
    mt_gene_count = 0
    mt_genes: List[str] = []
    qc_mt_genes: List[str] = []
    globin_panel_present: set[str] = set()
    globin_panel_used_for_qc: set[str] = set()
    false_hb_prefix_genes_excluded: set[str] = set()
    first_genes: List[str] = []
    gene_names: set[str] = set()
    gene_duplicate_count = 0
    gene_name_issue_count = 0
    row_label_cell_barcode_like_count = 0
    gene_means: List[float] = []
    gene_detected_in_0_cells = 0
    gene_detected_in_1_cell = 0
    gene_detected_in_2_cells = 0
    gene_detected_lt_3_examples: List[str] = []
    gene_hash = hashlib.sha256()
    global_min: float | None = None
    global_max: float | None = None
    ncount_complete = True

    # 第一行是细胞 barcode；后续每行是一个基因及其在所有细胞中的 counts。
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline()
        if not header:
            raise RuntimeError(f"Empty matrix file: {path}")
        header_fields = [clean_csv_token(value) for value in next(csv.reader([header.rstrip("\r\n")]))]
        expected_cells = len(header_fields) - 1
        if expected_cells <= 0:
            raise RuntimeError(f"Matrix header has no cell columns: {path}")
        first_header_blank = header_fields[0] == ""
        cell_barcodes = header_fields[1:]
        cell_barcode_duplicate_count = len(cell_barcodes) - len(set(cell_barcodes))
        cell_barcode_pattern = (
            "10x_16nt_barcode_with_numeric_suffix"
            if cell_barcodes and all(CELL_BARCODE.fullmatch(value) for value in cell_barcodes)
            else "mixed_or_unrecognized_cell_barcode"
        )
        # raw_full_ncount 使用全部存档行，只审计输入形态。
        # qc_* 数组只累加“至少在 3 个细胞检出”的基因，用于唯一一套细胞 QC。
        raw_full_ncount = np.zeros(expected_cells, dtype=np.int64)
        qc_ncount = np.zeros(expected_cells, dtype=np.int64)
        qc_nfeature = np.zeros(expected_cells, dtype=np.int32)
        qc_mt_ncount = np.zeros(expected_cells, dtype=np.int64)
        qc_hb_ncount = np.zeros(expected_cells, dtype=np.int64)
        first_header_fields = "|".join(header_fields[:6])
        barcode_suffixes = [
            barcode.rsplit("_", 1)[-1] if "_" in barcode else "" for barcode in header_fields[1:6]
        ]

        # 流式逐基因扫描：不创建完整矩阵，也不修改原文件。
        for line in handle:
            rows += 1
            stripped = line.rstrip("\r\n")
            first_comma = stripped.find(",")
            if first_comma < 0:
                gene = clean_csv_token(stripped)
                values = ""
                observed_cells = 0
            else:
                gene = clean_csv_token(stripped[:first_comma])
                values = stripped[first_comma + 1 :]
                observed_cells = values.count(",") + 1
            upper_gene = gene.upper()
            is_mt_gene = upper_gene.startswith("MT-")
            is_globin_panel_gene = upper_gene in GLOBIN_PANEL_SET
            if CELL_BARCODE.fullmatch(gene):
                row_label_cell_barcode_like_count += 1
            if observed_cells != expected_cells:
                bad_column_count_rows += 1
                ncount_complete = False
            total_values_checked += observed_cells

            row_array: np.ndarray | None = None
            if NONNEGATIVE_INTEGER_ROW.fullmatch(values):
                parsed = np.fromstring(values, dtype=np.int64, sep=",")
                if parsed.size == expected_cells:
                    row_array = parsed
                    integer_count += int(parsed.size)
                    zero_value_count += int(np.count_nonzero(parsed == 0))
                    row_min = float(np.min(parsed))
                    row_max = float(np.max(parsed))
                else:
                    ncount_complete = False
                    anomalies = classify_anomalous_values(values, expected_cells)
            else:
                anomalies = classify_anomalous_values(values, expected_cells)

            if row_array is None:
                missing_value_count += int(anomalies["missing_value_count"])
                noninteger_float_count += int(anomalies["noninteger_float_count"])
                invalid_nonnumeric_count += int(anomalies["invalid_nonnumeric_count"])
                negative_integer_count += int(anomalies["negative_integer_count"])
                integer_count += int(anomalies["integer_count"])
                zero_value_count += int(anomalies["zero_value_count"])
                row_min = anomalies["numeric_min"]
                row_max = anomalies["numeric_max"]
                if any(
                    int(anomalies[field]) > 0
                    for field in (
                        "missing_value_count",
                        "noninteger_float_count",
                        "invalid_nonnumeric_count",
                        "negative_integer_count",
                    )
                ):
                    ncount_complete = False
                integer_values = anomalies["integer_values"]
                if integer_values is not None:
                    row_array = np.asarray(integer_values, dtype=np.int64)
                else:
                    ncount_complete = False

            if row_min is not None:
                global_min = float(row_min) if global_min is None else min(global_min, float(row_min))
            if row_max is not None:
                global_max = float(row_max) if global_max is None else max(global_max, float(row_max))
            if row_array is not None and row_array.size == expected_cells:
                raw_full_ncount += row_array
                detected_cells = int(np.count_nonzero(row_array > 0))
                # 只有达到 min.cells=3 的基因才进入 QC 工作空间。
                if detected_cells >= QC_MIN_CELLS_PER_FEATURE:
                    qc_ncount += row_array
                    qc_nfeature += row_array > 0
                    if is_mt_gene:
                        qc_mt_ncount += row_array
                        qc_mt_genes.append(gene)
                    if is_globin_panel_gene:
                        qc_hb_ncount += row_array
                        globin_panel_used_for_qc.add(upper_gene)
                if detected_cells == 0:
                    gene_detected_in_0_cells += 1
                elif detected_cells == 1:
                    gene_detected_in_1_cell += 1
                elif detected_cells == 2:
                    gene_detected_in_2_cells += 1
                if detected_cells < QC_MIN_CELLS_PER_FEATURE and len(gene_detected_lt_3_examples) < 10:
                    gene_detected_lt_3_examples.append(f"{gene}:{detected_cells}")
                gene_means.append(float(np.sum(row_array, dtype=np.int64)) / expected_cells)
            else:
                ncount_complete = False

            if len(first_genes) < 5:
                first_genes.append(gene)
            if gene in gene_names:
                gene_duplicate_count += 1
            else:
                gene_names.add(gene)
            if not GENE_NAME.fullmatch(gene):
                gene_name_issue_count += 1
            gene_hash.update(gene.encode("utf-8"))
            gene_hash.update(b"\n")
            if is_mt_gene:
                mt_gene_count += 1
                if len(mt_genes) < 8:
                    mt_genes.append(gene)
            if is_globin_panel_gene:
                globin_panel_present.add(upper_gene)
            elif upper_gene.startswith("HB"):
                false_hb_prefix_genes_excluded.add(upper_gene)

    # 文件扫描结束后，把累计量转换成每个样本的一行结构化审计结果。
    # 矩阵方向不能由 .csv.gz 扩展名推断。只有“列名像 10x barcode、行名像基因、
    # 行名不整体像 barcode、左上角为空”同时成立，才确认 gene × cell。
    row_gene_identity_confirmed = (
        rows > 0
        and gene_name_issue_count == 0
        and row_label_cell_barcode_like_count == 0
    )
    orientation_confirmed = (
        first_header_blank
        and cell_barcode_pattern == "10x_16nt_barcode_with_numeric_suffix"
        and row_gene_identity_confirmed
    )
    orientation_evidence = (
        f"first_header_blank={'true' if first_header_blank else 'false'}; "
        f"column_barcode_pattern={cell_barcode_pattern}; "
        f"row_gene_name_issue_count={gene_name_issue_count}; "
        f"row_label_cell_barcode_like_count={row_label_cell_barcode_like_count}; "
        f"rows={rows}; columns={expected_cells}"
    )

    numeric_anomaly_count = (
        missing_value_count
        + noninteger_float_count
        + invalid_nonnumeric_count
        + negative_integer_count
    )
    observed_numeric_type = (
        "nonnegative_integer_count_like"
        if bad_column_count_rows == 0 and numeric_anomaly_count == 0
        else "format_or_numeric_issue_detected"
    )
    integer_value_rate = integer_count / total_values_checked if total_values_checked else math.nan
    ncount_summary = assess_ncount_distribution(raw_full_ncount, ncount_complete)
    fixed_qc_summary = assess_fixed_qc(
        qc_ncount,
        qc_nfeature,
        qc_mt_ncount,
        qc_hb_ncount,
        ncount_complete,
    )
    working_feature_summary = assess_working_feature_space(
        gene_detected_in_0_cells,
        gene_detected_in_1_cell,
        gene_detected_in_2_cells,
        gene_detected_lt_3_examples,
        rows,
        ncount_complete,
    )
    gene_mean_summary = summarize_gene_means(
        gene_means,
        rows,
        zero_value_count,
        total_values_checked,
    )
    globin_present_ordered = [gene for gene in GLOBIN_PANEL if gene in globin_panel_present]
    globin_missing_ordered = [gene for gene in GLOBIN_PANEL if gene not in globin_panel_present]
    globin_used_ordered = [gene for gene in GLOBIN_PANEL if gene in globin_panel_used_for_qc]
    globin_excluded_by_min_cells = [
        gene
        for gene in GLOBIN_PANEL
        if gene in globin_panel_present and gene not in globin_panel_used_for_qc
    ]
    result: Dict[str, object] = {
        "matrix_orientation": "gene_by_cell" if orientation_confirmed else "not_confirmed",
        "matrix_orientation_validation_status": (
            "pass_gene_by_cell" if orientation_confirmed else "fail_not_gene_by_cell"
        ),
        "matrix_orientation_validation_evidence": orientation_evidence,
        "matrix_rows_genes": rows,
        "matrix_cols_cells": expected_cells,
        "header_total_columns_including_gene_col": len(header_fields),
        "first_header_fields": first_header_fields,
        "first_genes": "|".join(first_genes),
        "barcode_suffix_examples": "|".join(barcode_suffixes),
        "row_identity": (
            "gene_symbol_or_ensembl_like_name"
            if row_gene_identity_confirmed
            else "cell_barcode_like_or_invalid_row_labels"
        ),
        "column_identity": (
            "cell_barcode"
            if cell_barcode_pattern == "10x_16nt_barcode_with_numeric_suffix"
            else "not_confirmed_cell_barcode"
        ),
        "first_header_blank": "true" if first_header_blank else "false",
        "gene_name_type": (
            "gene_symbol_or_ensembl_like_name"
            if row_gene_identity_confirmed
            else "cell_barcode_like_or_invalid_row_labels"
        ),
        "gene_name_issue_count": gene_name_issue_count,
        "row_label_cell_barcode_like_count": row_label_cell_barcode_like_count,
        "cell_barcode_pattern": cell_barcode_pattern,
        "gene_duplicate_count": gene_duplicate_count,
        "cell_barcode_duplicate_count": cell_barcode_duplicate_count,
        "bad_column_count_rows": bad_column_count_rows,
        "integer_check_method": "full_stream",
        "integer_parser": "numpy_fromstring_with_exact_anomaly_fallback",
        "total_values_checked": total_values_checked,
        "integer_value_rate": format(integer_value_rate, ".12g") if math.isfinite(integer_value_rate) else "",
        "missing_value_count": missing_value_count,
        "noninteger_float_count": noninteger_float_count,
        "invalid_nonnumeric_count": invalid_nonnumeric_count,
        "negative_integer_count": negative_integer_count,
        "numeric_anomaly_count": numeric_anomaly_count,
        "min_value": format_number(global_min),
        "max_value": format_number(global_max),
        "has_negative_value": "true" if global_min is not None and global_min < 0 else "false",
        **ncount_summary,
        **working_feature_summary,
        **fixed_qc_summary,
        **gene_mean_summary,
        "mt_gene_count": mt_gene_count,
        "mt_gene_examples": "|".join(mt_genes),
        "qc_mt_gene_count": len(qc_mt_genes),
        "qc_mt_gene_examples": "|".join(qc_mt_genes[:8]),
        "globin_panel_expected": "|".join(GLOBIN_PANEL),
        "globin_panel_present": "|".join(globin_present_ordered),
        "globin_panel_missing": "|".join(globin_missing_ordered) or "none",
        "globin_panel_used_for_qc": "|".join(globin_used_ordered),
        "globin_panel_excluded_by_min_cells3": "|".join(globin_excluded_by_min_cells) or "none",
        "false_hb_prefix_genes_excluded": (
            "|".join(sorted(false_hb_prefix_genes_excluded)) or "none"
        ),
        "gene_order_sha256": gene_hash.hexdigest().upper(),
        "observed_numeric_type": observed_numeric_type,
        "scan_seconds": round(time.perf_counter() - start, 2),
    }
    result.update(validate_frozen_pilot(path.name, result, globin_used_ordered))
    return result


def legacy_numeric_precheck(stats: Dict[str, object], precheck: Dict[str, str]) -> Tuple[str, List[str]]:
    """把旧的行级异常预检查与本次更严格的逐值审计作兼容核对。

    旧 HB 计数没有在这里比较，因为旧定义使用宽泛前缀，新方案已冻结为精确
    globin panel；把二者强行比较会制造没有方法学意义的“不一致”。
    """

    mismatches: List[str] = []
    missing_rows = precheck.get("missing_value_rows", "")
    invalid_rows = precheck.get("invalid_numeric_value_rows", "")
    if missing_rows == "0" and int(stats["missing_value_count"]) != 0:
        mismatches.append(
            f"legacy missing_value_rows=0 but missing_value_count={stats['missing_value_count']}"
        )
    if invalid_rows == "0":
        observed = (
            int(stats["noninteger_float_count"])
            + int(stats["invalid_nonnumeric_count"])
            + int(stats["negative_integer_count"])
        )
        if observed != 0:
            mismatches.append(f"legacy invalid_numeric_value_rows=0 but numeric anomalies={observed}")
    if missing_rows not in {"", "0"} or invalid_rows not in {"", "0"}:
        mismatches.append("nonzero legacy row-level anomaly counts cannot be equated to new value-level counts")
    return ("match" if not mismatches else "mismatch"), mismatches


def build_data_audit(
    root: Path,
    processed_manifest: Sequence[Dict[str, str]],
    sample_info: Sequence[Dict[str, object]],
    precheck_rows: Sequence[Dict[str, str]],
) -> Tuple[List[Dict[str, object]], List[str]]:
    """对 40 个矩阵运行完整审计，并与 Step1、GEO 和预检查逐项核对。

    只有文件结构、数值类型、样本映射、工作空间、固定 QC 和 pilot 均可评估
    时，样本才得到 ``enter_full_F1_independent_reQC``。该状态只表示允许 F1
    重新 QC，不表示作者完整处理史已经还原。
    """

    sample_by_file = {str(row["sample_file"]): row for row in sample_info}
    precheck_by_member = {row["member_name"]: row for row in precheck_rows}
    mismatch_notes: List[str] = []
    output: List[Dict[str, object]] = []

    # 每个样本独立审计，避免把不同测序深度和样本来源混为一个总体。
    for manifest_row in processed_manifest:
        path = root / manifest_row["extracted_path"]
        stats = audit_csv_gz(path)
        member = manifest_row["archive_member_name"]
        sample = sample_by_file.get(member, {})
        precheck = precheck_by_member.get(member, {})
        row_mismatches: List[str] = []
        if not precheck:
            row_mismatches.append("precheck row missing for archive member")
        # 复核此前已登记的关键结构事实；新旧结果不一致时必须暂停解释。
        for field in [
            "matrix_rows_genes",
            "matrix_cols_cells",
            "mt_gene_count",
            "gene_order_sha256",
        ]:
            observed = str(stats.get(field, ""))
            expected = str(precheck.get(field, ""))
            equal = (
                sha256_equal(observed, expected)
                if field == "gene_order_sha256"
                else observed == expected
            )
            if expected and not equal:
                row_mismatches.append(f"{field}: observed={observed}; precheck={expected}")
        observed_numeric_type = str(stats.get("observed_numeric_type", ""))
        expected_numeric_type = str(precheck.get("suspected_matrix_type", ""))
        if expected_numeric_type and observed_numeric_type != expected_numeric_type:
            row_mismatches.append(
                "observed_numeric_type: "
                f"observed={observed_numeric_type}; precheck.suspected_matrix_type={expected_numeric_type}"
            )

        legacy_status, legacy_mismatches = legacy_numeric_precheck(stats, precheck)
        row_mismatches.extend(legacy_mismatches)
        pre_decision = precheck.get("audit_decision_precheck", "")
        processed_match = (
            sample.get("sample_file") == member
            and sample.get("geo_accession") == manifest_row.get("geo_accession")
            and sample.get("sample_id") == manifest_row.get("sample_id")
            and sample.get("sample_id_match_status") == "match"
        )
        if not processed_match:
            row_mismatches.append("processed manifest, GEO accession/title and sample_info mapping do not match")
        precheck_match = not row_mismatches
        # 将每类失败原因单独保存，便于用户和审核者判断问题出在哪里。
        structural_failures: List[str] = []
        if stats["matrix_orientation_validation_status"] != "pass_gene_by_cell":
            structural_failures.append("matrix_orientation_not_confirmed_as_gene_by_cell")
        if stats["observed_numeric_type"] != "nonnegative_integer_count_like":
            structural_failures.append("numeric_or_column_structure_issue")
        if stats["first_header_blank"] != "true":
            structural_failures.append("first_header_not_blank")
        if stats["row_identity"] != "gene_symbol_or_ensembl_like_name":
            structural_failures.append("gene_name_identity_issue")
        if int(stats["gene_duplicate_count"]) != 0:
            structural_failures.append("duplicate_gene_names")
        if int(stats["cell_barcode_duplicate_count"]) != 0:
            structural_failures.append("duplicate_cell_barcodes")
        if stats["per_gene_mean_distribution_status"] != "consistent_with_sparse_right_skew":
            structural_failures.append("per_gene_mean_distribution_not_confirmed")
        if stats["normalization_artifact_flag"] != "false":
            structural_failures.append("normalization_artifact_or_not_evaluable")
        fixed_qc_status = str(stats["fixed_qc_rule_recalculation_status"])
        if fixed_qc_status != "pass":
            structural_failures.append("fixed_QC_rule_not_evaluable")
        working_feature_status = str(stats["working_feature_space_recalculation_status"])
        if working_feature_status != "pass":
            structural_failures.append("min_cells3_working_feature_space_not_evaluable")
        pilot_applicable = str(stats["pilot_validation_applicable"])
        pilot_status = str(stats["pilot_validation_status"])
        if pilot_applicable == "true" and pilot_status != "pass":
            structural_failures.append("sample1_frozen_pilot_regression_failed")
        if pilot_applicable == "false" and pilot_status != "not_applicable":
            structural_failures.append("nonpilot_sample_has_invalid_pilot_status")
        if str(stats.get("globin_panel_expected", "")) != "|".join(GLOBIN_PANEL):
            structural_failures.append("frozen_globin_panel_definition_mismatch")
        if not str(stats.get("globin_panel_present", "")):
            structural_failures.append("no_frozen_globin_panel_gene_present_in_matrix")
        structural_ok = not structural_failures
        processing_boundary_ok = (
            structural_ok
            and processed_match
            and not row_mismatches
            and pre_decision == "enter_full_F1_candidate"
        )
        suspected_matrix_type = (
            "public_called_cell_raw_gene_count_matrix"
            if processing_boundary_ok
            else "processing_boundary_not_confirmed"
        )
        ok = processing_boundary_ok
        if row_mismatches:
            mismatch_notes.append(f"{member} | " + " | ".join(row_mismatches))
        failure_details = [*structural_failures, *row_mismatches]
        if pre_decision != "enter_full_F1_candidate":
            failure_details.append(f"unexpected_precheck_decision={pre_decision or 'missing'}")
        if ok:
            decision_reason = (
                "full-stream numeric/structure audit passed; per-cell nCount showed no preregistered "
                "normalization artifact; row/column identities confirmed gene-by-cell orientation; all 26571 "
                "archived feature rows remain unchanged; one per-sample "
                "min.cells=3 working feature space was built and nCount, nFeature, percent.mt and percent.HB "
                "were recalculated there; the fixed source-aligned plus project QC rule was fully evaluable; "
                "the frozen globin panel was applied by exact gene matching; sparse right-skew evidence was "
                "present; processed manifest and GEO sample mapping matched; "
                f"working_features={stats['qc_retained_feature_count']}; "
                f"source_rule_pass={stats['source_reported_qc_pass_count']}; "
                f"final_fixed_QC_pass={stats['final_fixed_qc_pass_count']}; pilot={pilot_status}"
            )
        else:
            decision_reason = (
                "pause because one or more format, distribution, duplicate, mapping or precheck conditions failed: "
                + "; ".join(failure_details)
                + f"; artifact_evidence={stats['normalization_artifact_reason']}"
            )
        boundary_notes = [
            "Rows detected in fewer than three cells remain in the immutable archived matrix; they are "
            "excluded only from this sample's QC working feature space and are not labeled as a provenance mismatch."
        ]
        output.append(
            {
                "dataset_id": "GSE183904",
                "file_name": member,
                "geo_accession": manifest_row["geo_accession"],
                "sample_id": manifest_row["sample_id"],
                "extracted_path": manifest_row["extracted_path"],
                "source_archive": manifest_row.get("source_archive", ""),
                "processed_input_manifest_match": "true" if processed_match else "false",
                "file_format": "csv",
                "compressed": "true",
                **stats,
                "legacy_numeric_precheck_status": legacy_status,
                "legacy_hb_precheck_status": "not_compared_definition_changed_to_frozen_exact_globin_panel",
                "precheck_audit_decision": pre_decision,
                "public_processing_evidence_status": (
                    "public_input_shape_verified_fixed_QC_recalculated_processing_history_pending_F0_step3"
                    if ok
                    else "not_confirmed_due_to_F0_step2_pause_condition"
                ),
                "suspected_matrix_type": suspected_matrix_type,
                "audit_decision": "enter_full_F1_independent_reQC" if ok else "pause_for_review",
                "format_decision_scope": "file_format_only",
                "decision_scope": "file_format_and_public_processing_boundary",
                "decision_reason": decision_reason,
                "include_in_f1": sample.get("include_in_f1", "pending"),
                "raw_droplet_available": "false",
                "empty_droplet_background_available": "false",
                "precheck_comparison_status": "match" if precheck_match else "mismatch",
                "note": "; ".join([*row_mismatches, *boundary_notes]),
            }
        )
    return output, mismatch_notes


def validate_step1_manifest(root: Path, rows: Sequence[Dict[str, str]]) -> None:
    """再次核对 Step1 的 40 个提取文件及 SHA256，防止阶段间文件被替换。"""

    if len(rows) != 40:
        raise RuntimeError(f"processed_input_manifest.tsv must contain 40 rows; observed {len(rows)}")
    seen: set[str] = set()
    for row in rows:
        member = row.get("archive_member_name", "")
        if not member or member in seen:
            raise RuntimeError(f"Duplicate or empty archive member in processed manifest: {member!r}")
        seen.add(member)
        if (
            row.get("file_role") != "expected_gene_by_cell_matrix_pending_validation"
            or row.get("expected_matrix_orientation") != "gene_by_cell"
            or row.get("orientation_validation_status")
            != "pending_F0_step2_full_stream_validation"
        ):
            raise RuntimeError(
                f"Step1 orientation contract is invalid for {member}: "
                f"file_role={row.get('file_role', '')}; "
                f"expected={row.get('expected_matrix_orientation', '')}; "
                f"status={row.get('orientation_validation_status', '')}"
            )
        path = root / row.get("extracted_path", "")
        if not path.exists():
            raise RuntimeError(f"Extracted matrix missing: {path}")
        observed_sha = sha256_file(path)
        expected_sha = normalize_sha256(row.get("sha256", ""))
        if not expected_sha or not sha256_equal(observed_sha, expected_sha):
            raise RuntimeError(
                f"Extracted matrix SHA256 mismatch for {member}: observed={observed_sha}; expected={expected_sha}"
            )


def execute(root: Path) -> int:
    """正式执行 Step2，写出两个核心表，并在任何阻断条件出现时停止。"""

    # 1. 输入与 Step1 文件完整性核验。
    require_paths(root, STAGE_REQUIRED, STAGE_NAME)
    validate_f0_python_environment(
        root,
        actual_python_version=sys.version.split()[0],
        actual_numpy_version=np.__version__,
    )
    manifest = read_tsv(root / "data/metadata/processed_input_manifest.tsv")
    validate_step1_manifest(root, manifest)
    append_log(
        root,
        f"F0 step2 started; run_id={current_run_id()}; validated_step1_files=40; "
        f"python={sys.version.split()[0]}; numpy={np.__version__}",
    )

    # 2. 建立样本映射。文件名和 GEO title 必须独立一致，不能按顺序猜测。
    sample_fields, patient_by_sample = parse_gse183904_series(
        root / "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz"
    )
    sample_info, sample_id_mismatches = build_sample_info(sample_fields, patient_by_sample, manifest)
    sample_fields_out = [
        "dataset_id",
        "sample_file",
        "geo_accession",
        "sample_id",
        "geo_title_sample_id",
        "sample_id_match_status",
        "patient_id",
        "group_original",
        "group_analysis",
        "group_analysis_rule",
        "sample_title",
        "source_name_ch1",
        "sample_characteristics_ch1",
        "tissue_site",
        "tumor_status",
        "metastasis_site",
        "is_paired",
        "paired_normal_id",
        "include_in_f1",
        "include_in_group_comparison",
        "include_reason",
        "source_of_group",
        "metadata_confidence",
        "metadata_issue",
        "pairing_scope",
        "note",
    ]
    write_tsv(root / "data/metadata/sample_info.tsv", sample_info, sample_fields_out)

    # 3. 对 40 个矩阵进行完整流式审计并与预检查对照。
    precheck = read_tsv(root / "results/F0_audit/gse183904_csv_structure_precheck.tsv")
    data_audit, mismatches = build_data_audit(root, manifest, sample_info, precheck)
    append_log(
        root,
        "SCHEMA_MIGRATION precheck.suspected_matrix_type -> observed_numeric_type; "
        "precheck.audit_decision_precheck=enter_full_F1_candidate -> "
        "audit_decision=enter_full_F1_independent_reQC only after formal audit; "
        "legacy broad HB count is not compared because HB_percent now uses the frozen exact globin panel; "
        "Step1 orientation is an expectation only and Step2 validates row/column identities; "
        "suspected_matrix_type records public input shape only; raw_full_nCount and the single min.cells=3 "
        "QC metric space are stored separately",
    )
    data_fields = [
        "dataset_id",
        "file_name",
        "geo_accession",
        "sample_id",
        "extracted_path",
        "source_archive",
        "processed_input_manifest_match",
        "file_format",
        "compressed",
        "matrix_orientation",
        "matrix_orientation_validation_status",
        "matrix_orientation_validation_evidence",
        "matrix_rows_genes",
        "matrix_cols_cells",
        "header_total_columns_including_gene_col",
        "first_header_fields",
        "first_genes",
        "row_identity",
        "column_identity",
        "first_header_blank",
        "gene_name_type",
        "gene_name_issue_count",
        "row_label_cell_barcode_like_count",
        "cell_barcode_pattern",
        "barcode_suffix_examples",
        "gene_duplicate_count",
        "cell_barcode_duplicate_count",
        "bad_column_count_rows",
        "integer_check_method",
        "integer_parser",
        "total_values_checked",
        "integer_value_rate",
        "missing_value_count",
        "noninteger_float_count",
        "invalid_nonnumeric_count",
        "negative_integer_count",
        "numeric_anomaly_count",
        "min_value",
        "max_value",
        "has_negative_value",
        "zero_value_count",
        "zero_value_rate",
        "raw_full_nCount_min",
        "raw_full_nCount_Q1",
        "raw_full_nCount_median",
        "raw_full_nCount_Q3",
        "raw_full_nCount_max",
        "qc_min_cells_per_feature",
        "feature_rows_detected_in_0_cells",
        "feature_rows_detected_in_1_cell",
        "feature_rows_detected_in_2_cells",
        "feature_rows_detected_lt_3_count",
        "feature_rows_detected_lt_3_examples",
        "qc_retained_feature_count",
        "working_feature_space_recalculation_status",
        "working_feature_space_recalculation_note",
        "qc_nCount_min",
        "qc_nCount_Q1",
        "qc_nCount_median",
        "qc_nCount_Q3",
        "qc_nCount_max",
        "qc_nFeature_min",
        "qc_nFeature_Q1",
        "qc_nFeature_median",
        "qc_nFeature_Q3",
        "qc_nFeature_max",
        "qc_percent_mt_min",
        "qc_percent_mt_Q1",
        "qc_percent_mt_median",
        "qc_percent_mt_Q3",
        "qc_percent_mt_max",
        "qc_percent_hb_min",
        "qc_percent_hb_Q1",
        "qc_percent_hb_median",
        "qc_percent_hb_Q3",
        "qc_percent_hb_max",
        "fail_nFeature_low_count",
        "fail_nFeature_high_count",
        "fail_nCount_count",
        "fail_percent_mt_count",
        "fail_percent_hb_count",
        "source_reported_qc_fail_count",
        "source_reported_qc_pass_count",
        "additional_fail_nCount_after_source_count",
        "additional_fail_percent_hb_after_source_nCount_count",
        "final_fixed_qc_fail_count",
        "final_fixed_qc_pass_count",
        "fixed_qc_rule_recalculation_status",
        "fixed_qc_rule_recalculation_note",
        "pilot_validation_applicable",
        "pilot_validation_status",
        "pilot_validation_note",
        "ncount_distinct_count",
        "ncount_relative_iqr",
        "ncount_relative_range",
        "dominant_round_ncount_fraction",
        "fixed_target_near_value",
        "fixed_target_near_fraction",
        "ncount_range_status",
        "per_gene_mean_distribution_status",
        "per_gene_mean_distribution_note",
        "normalization_artifact_flag",
        "normalization_artifact_reason",
        "mt_gene_count",
        "mt_gene_examples",
        "qc_mt_gene_count",
        "qc_mt_gene_examples",
        "globin_panel_expected",
        "globin_panel_present",
        "globin_panel_missing",
        "globin_panel_used_for_qc",
        "globin_panel_excluded_by_min_cells3",
        "false_hb_prefix_genes_excluded",
        "gene_order_sha256",
        "observed_numeric_type",
        "public_processing_evidence_status",
        "suspected_matrix_type",
        "scan_seconds",
        "legacy_numeric_precheck_status",
        "legacy_hb_precheck_status",
        "precheck_audit_decision",
        "audit_decision",
        "format_decision_scope",
        "decision_scope",
        "decision_reason",
        "include_in_f1",
        "raw_droplet_available",
        "empty_droplet_background_available",
        "precheck_comparison_status",
        "note",
    ]
    write_tsv(root / "data/metadata/data_audit.tsv", data_audit, data_fields)
    for mismatch in mismatches:
        append_log(root, "PRECHECK_MISMATCH " + mismatch)
    for mismatch in sample_id_mismatches:
        append_log(root, "SAMPLE_ID_MISMATCH " + mismatch)

    # 4. 汇总所有阻断条件。任一问题存在时，Step3 不应被自动执行。
    unclear = [row for row in sample_info if row.get("group_analysis") == "Unclear"]
    artifacts = [row for row in data_audit if row.get("normalization_artifact_flag") != "false"]
    working_feature_failures = [
        row
        for row in data_audit
        if row.get("working_feature_space_recalculation_status") != "pass"
    ]
    fixed_qc_failures = [
        row
        for row in data_audit
        if row.get("fixed_qc_rule_recalculation_status") != "pass"
    ]
    pilot_failures = [
        row
        for row in data_audit
        if row.get("pilot_validation_applicable") == "true"
        and row.get("pilot_validation_status") != "pass"
    ]
    paused = [
        row
        for row in data_audit
        if row.get("audit_decision") != "enter_full_F1_independent_reQC"
    ]
    append_log(
        root,
        f"F0 step2 completed; sample_rows={len(sample_info)}; audit_rows={len(data_audit)}; "
        f"precheck_mismatches={len(mismatches)}; sample_id_mismatches={len(sample_id_mismatches)}; "
        f"unclear_groups={len(unclear)}; artifacts_or_not_evaluable={len(artifacts)}; "
        f"working_feature_failures={len(working_feature_failures)}; "
        f"fixed_qc_failures={len(fixed_qc_failures)}; pilot_failures={len(pilot_failures)}; "
        f"paused={len(paused)}",
    )
    if len(sample_info) != 40 or unclear or sample_id_mismatches or mismatches or artifacts or paused:
        raise RuntimeError(
            "F0 step2 pause condition reached: "
            f"sample_rows={len(sample_info)}, unclear={len(unclear)}, sample_id_mismatches={len(sample_id_mismatches)}, "
            f"precheck_mismatches={len(mismatches)}, normalization_artifacts_or_not_evaluable={len(artifacts)}, "
            f"paused={len(paused)}"
        )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """默认仅列出输入输出；显式提供 ``--execute`` 后才运行完整审计。"""

    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, STAGE_REQUIRED, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
