#!/usr/bin/env python3
"""F0 step 2: build sample_info and full-stream audit 40 GSE183904 matrices.

Dependencies: step 1 processed_input_manifest.tsv and extracted CSV.gz files,
the GSE183904 GEO series matrix, and the preregistered structure precheck.
Outputs: sample_info.tsv and data_audit.tsv. Any unresolved group, key-field
precheck mismatch, or numeric/format issue stops the staged wrapper.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
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
    write_tsv,
)


STAGE_NAME = "F0 step 2 sample info and matrix audit"
STAGE_REQUIRED = [
    "data/metadata/processed_input_manifest.tsv",
    "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
    "results/F0_audit/gse183904_csv_structure_precheck.tsv",
]
STAGE_OUTPUTS = [
    "data/metadata/sample_info.tsv",
    "data/metadata/data_audit.tsv",
]

NONNEGATIVE_INTEGER_ROW = re.compile(r"\d+(?:,\d+)*\Z")
GENE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
CELL_BARCODE = re.compile(r"[ACGTN]{16}_\d+\Z", re.IGNORECASE)
SAMPLE_ID = re.compile(r"sample\d+\Z", re.IGNORECASE)

# Project-specific, high-specificity sentinels preregistered before inspecting
# the new per-cell nCount results. Fixed targets reflect common CP1k/CP10k/
# CP100k/CPM library-size normalization scales; the concentration thresholds
# are deliberately stringent and trigger review rather than proving a method.
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
AUTHOR_MIN_NFEATURE = 500
AUTHOR_MAX_NFEATURE_EXCLUSIVE = 6000
AUTHOR_MAX_PERCENT_MT = 20.0
AUTHOR_MIN_CELLS_PER_FEATURE = 3


def parse_geo_line(line: str) -> Tuple[str, List[str]]:
    parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
    return parts[0], [value.strip().strip('"') for value in parts[1:]]


def parse_gse183904_series(series_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
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
    if value is None or not math.isfinite(float(value)):
        return ""
    numeric = float(value)
    if numeric.is_integer():
        return str(int(numeric))
    return format(numeric, ".12g")


def classify_anomalous_values(values: str, expected_cells: int) -> Dict[str, object]:
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
    blank = {
        "per_cell_nCount_min": "",
        "per_cell_nCount_Q1": "",
        "per_cell_nCount_median": "",
        "per_cell_nCount_Q3": "",
        "per_cell_nCount_max": "",
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
        "per_cell_nCount_min": format_number(minimum),
        "per_cell_nCount_Q1": format_number(q1),
        "per_cell_nCount_median": format_number(median),
        "per_cell_nCount_Q3": format_number(q3),
        "per_cell_nCount_max": format_number(maximum),
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


def assess_cell_qc_space(
    ncount: np.ndarray,
    nfeature: np.ndarray,
    mt_ncount: np.ndarray,
    complete: bool,
    metric_prefix: str,
    threshold_suffix: str,
) -> Dict[str, object]:
    metric_fields = {
        f"{metric_prefix}_nCount_min": "",
        f"{metric_prefix}_nCount_Q1": "",
        f"{metric_prefix}_nCount_median": "",
        f"{metric_prefix}_nCount_Q3": "",
        f"{metric_prefix}_nCount_max": "",
        f"{metric_prefix}_nFeature_min": "",
        f"{metric_prefix}_nFeature_Q1": "",
        f"{metric_prefix}_nFeature_median": "",
        f"{metric_prefix}_nFeature_Q3": "",
        f"{metric_prefix}_nFeature_max": "",
        f"{metric_prefix}_percent_mt_min": "",
        f"{metric_prefix}_percent_mt_Q1": "",
        f"{metric_prefix}_percent_mt_median": "",
        f"{metric_prefix}_percent_mt_Q3": "",
        f"{metric_prefix}_percent_mt_max": "",
    }
    threshold_fields = {
        f"author_nFeature_lt_500_count_{threshold_suffix}": "",
        f"author_nFeature_ge_6000_count_{threshold_suffix}": "",
        f"author_percent_mt_gt_20_count_{threshold_suffix}": "",
        f"author_cell_threshold_mismatch_count_{threshold_suffix}": "",
        f"author_cell_qc_reproduction_status_{threshold_suffix}": "not_evaluable",
        f"author_cell_qc_reproduction_note_{threshold_suffix}": (
            "cell_QC_not_evaluable_due_to_column_numeric_or_denominator_anomaly"
        ),
    }
    if (
        not complete
        or ncount.size == 0
        or nfeature.size != ncount.size
        or mt_ncount.size != ncount.size
        or np.any(ncount <= 0)
    ):
        return {**metric_fields, **threshold_fields}

    percent_mt = np.divide(
        mt_ncount.astype(np.float64) * 100.0,
        ncount.astype(np.float64),
    )
    feature_low = nfeature < AUTHOR_MIN_NFEATURE
    feature_high = nfeature >= AUTHOR_MAX_NFEATURE_EXCLUSIVE
    mt_high = percent_mt > AUTHOR_MAX_PERCENT_MT
    any_mismatch = feature_low | feature_high | mt_high
    ncount_q1, ncount_median, ncount_q3 = np.quantile(
        ncount, [0.25, 0.5, 0.75], method="linear"
    )
    nfeature_q1, nfeature_median, nfeature_q3 = np.quantile(
        nfeature, [0.25, 0.5, 0.75], method="linear"
    )
    mt_q1, mt_median, mt_q3 = np.quantile(
        percent_mt, [0.25, 0.5, 0.75], method="linear"
    )
    low_count = int(np.count_nonzero(feature_low))
    high_count = int(np.count_nonzero(feature_high))
    mt_high_count = int(np.count_nonzero(mt_high))
    mismatch_count = int(np.count_nonzero(any_mismatch))
    status = "pass" if mismatch_count == 0 else "measured_mismatch"
    note = (
        f"reported_baseline=500<=nFeature<6000_and_percent.mt<=20; "
        f"nFeature_lt_500={low_count}; nFeature_ge_6000={high_count}; "
        f"percent.mt_gt_20={mt_high_count}; unique_mismatching_cells={mismatch_count}"
    )
    return {
        f"{metric_prefix}_nCount_min": format_number(int(np.min(ncount))),
        f"{metric_prefix}_nCount_Q1": format_number(ncount_q1),
        f"{metric_prefix}_nCount_median": format_number(ncount_median),
        f"{metric_prefix}_nCount_Q3": format_number(ncount_q3),
        f"{metric_prefix}_nCount_max": format_number(int(np.max(ncount))),
        f"{metric_prefix}_nFeature_min": format_number(int(np.min(nfeature))),
        f"{metric_prefix}_nFeature_Q1": format_number(nfeature_q1),
        f"{metric_prefix}_nFeature_median": format_number(nfeature_median),
        f"{metric_prefix}_nFeature_Q3": format_number(nfeature_q3),
        f"{metric_prefix}_nFeature_max": format_number(int(np.max(nfeature))),
        f"{metric_prefix}_percent_mt_min": format_number(float(np.min(percent_mt))),
        f"{metric_prefix}_percent_mt_Q1": format_number(mt_q1),
        f"{metric_prefix}_percent_mt_median": format_number(mt_median),
        f"{metric_prefix}_percent_mt_Q3": format_number(mt_q3),
        f"{metric_prefix}_percent_mt_max": format_number(float(np.max(percent_mt))),
        f"author_nFeature_lt_500_count_{threshold_suffix}": low_count,
        f"author_nFeature_ge_6000_count_{threshold_suffix}": high_count,
        f"author_percent_mt_gt_20_count_{threshold_suffix}": mt_high_count,
        f"author_cell_threshold_mismatch_count_{threshold_suffix}": mismatch_count,
        f"author_cell_qc_reproduction_status_{threshold_suffix}": status,
        f"author_cell_qc_reproduction_note_{threshold_suffix}": note,
    }


def assess_author_feature_filter_boundary(
    detected_in_0_cells: int,
    detected_in_1_cell: int,
    detected_in_2_cells: int,
    examples: Sequence[str],
    total_features: int,
    complete: bool,
) -> Dict[str, object]:
    below_threshold = detected_in_0_cells + detected_in_1_cell + detected_in_2_cells
    if not complete:
        return {
            "author_min_cells_per_feature_reported": AUTHOR_MIN_CELLS_PER_FEATURE,
            "feature_rows_detected_in_0_cells": "",
            "feature_rows_detected_in_1_cell": "",
            "feature_rows_detected_in_2_cells": "",
            "feature_rows_detected_lt_3_count": "",
            "feature_rows_detected_lt_3_examples": "",
            "author_like_retained_feature_count": "",
            "author_feature_filter_reproduction_status": "not_evaluable",
            "author_feature_filter_reproduction_note": (
                "per_sample_feature_detection_not_evaluable_due_to_column_or_numeric_anomaly"
            ),
        }
    status = "pass" if below_threshold == 0 else "measured_mismatch"
    return {
        "author_min_cells_per_feature_reported": AUTHOR_MIN_CELLS_PER_FEATURE,
        "feature_rows_detected_in_0_cells": detected_in_0_cells,
        "feature_rows_detected_in_1_cell": detected_in_1_cell,
        "feature_rows_detected_in_2_cells": detected_in_2_cells,
        "feature_rows_detected_lt_3_count": below_threshold,
        "feature_rows_detected_lt_3_examples": "|".join(examples),
        "author_like_retained_feature_count": total_features - below_threshold,
        "author_feature_filter_reproduction_status": status,
        "author_feature_filter_reproduction_note": (
            f"reported_per_sample_baseline=feature_detected_in_at_least_3_cells; "
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


def audit_csv_gz(path: Path) -> Dict[str, object]:
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
    hb_gene_count = 0
    mt_genes: List[str] = []
    hb_genes: List[str] = []
    first_genes: List[str] = []
    gene_names: set[str] = set()
    gene_duplicate_count = 0
    gene_name_issue_count = 0
    gene_means: List[float] = []
    gene_detected_in_0_cells = 0
    gene_detected_in_1_cell = 0
    gene_detected_in_2_cells = 0
    gene_detected_lt_3_examples: List[str] = []
    gene_hash = hashlib.sha256()
    global_min: float | None = None
    global_max: float | None = None
    ncount_complete = True

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
        ncount = np.zeros(expected_cells, dtype=np.int64)
        nfeature = np.zeros(expected_cells, dtype=np.int32)
        mt_ncount = np.zeros(expected_cells, dtype=np.int64)
        author_like_ncount = np.zeros(expected_cells, dtype=np.int64)
        author_like_nfeature = np.zeros(expected_cells, dtype=np.int32)
        author_like_mt_ncount = np.zeros(expected_cells, dtype=np.int64)
        first_header_fields = "|".join(header_fields[:6])
        barcode_suffixes = [
            barcode.rsplit("_", 1)[-1] if "_" in barcode else "" for barcode in header_fields[1:6]
        ]

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
                ncount += row_array
                nfeature += row_array > 0
                if is_mt_gene:
                    mt_ncount += row_array
                detected_cells = int(np.count_nonzero(row_array > 0))
                if detected_cells >= AUTHOR_MIN_CELLS_PER_FEATURE:
                    author_like_ncount += row_array
                    author_like_nfeature += row_array > 0
                    if is_mt_gene:
                        author_like_mt_ncount += row_array
                if detected_cells == 0:
                    gene_detected_in_0_cells += 1
                elif detected_cells == 1:
                    gene_detected_in_1_cell += 1
                elif detected_cells == 2:
                    gene_detected_in_2_cells += 1
                if detected_cells < AUTHOR_MIN_CELLS_PER_FEATURE and len(gene_detected_lt_3_examples) < 10:
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
            if re.match(r"^HB[A-Z0-9]", upper_gene) and not upper_gene.startswith(("HBS", "HBP")):
                hb_gene_count += 1
                if len(hb_genes) < 8:
                    hb_genes.append(gene)

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
    ncount_summary = assess_ncount_distribution(ncount, ncount_complete)
    public_space_qc_summary = assess_cell_qc_space(
        ncount,
        nfeature,
        mt_ncount,
        ncount_complete,
        "public_full_feature",
        "public_space",
    )
    author_like_qc_summary = assess_cell_qc_space(
        author_like_ncount,
        author_like_nfeature,
        author_like_mt_ncount,
        ncount_complete,
        "author_like",
        "author_like_space",
    )
    author_feature_filter_summary = assess_author_feature_filter_boundary(
        gene_detected_in_0_cells,
        gene_detected_in_1_cell,
        gene_detected_in_2_cells,
        gene_detected_lt_3_examples,
        rows,
        ncount_complete,
    )
    public_space_status = str(
        public_space_qc_summary["author_cell_qc_reproduction_status_public_space"]
    )
    author_like_status = str(
        author_like_qc_summary["author_cell_qc_reproduction_status_author_like_space"]
    )
    if "not_evaluable" in {public_space_status, author_like_status}:
        author_cell_qc_reproduction_status = "not_evaluable"
    elif public_space_status == "pass" and author_like_status == "pass":
        author_cell_qc_reproduction_status = "pass"
    else:
        author_cell_qc_reproduction_status = "measured_mismatch"
    author_cell_qc_reproduction_note = (
        f"public_full_feature_space={public_space_status}; "
        f"author_like_min_cells3_feature_space={author_like_status}; "
        "the two spaces are reported separately and are not interchangeable"
    )
    if ncount_complete and np.all(ncount > 0) and np.all(author_like_ncount > 0):
        nfeature_decrease = nfeature.astype(np.int64) - author_like_nfeature.astype(np.int64)
        ncount_decrease = ncount - author_like_ncount
        public_percent_mt = mt_ncount.astype(np.float64) * 100.0 / ncount.astype(np.float64)
        author_like_percent_mt = (
            author_like_mt_ncount.astype(np.float64)
            * 100.0
            / author_like_ncount.astype(np.float64)
        )
        percent_mt_change = np.abs(public_percent_mt - author_like_percent_mt)
        feature_space_difference_summary: Dict[str, object] = {
            "feature_space_nFeature_changed_cell_count": int(
                np.count_nonzero(nfeature_decrease != 0)
            ),
            "feature_space_nFeature_max_decrease": int(np.max(nfeature_decrease)),
            "feature_space_nCount_changed_cell_count": int(np.count_nonzero(ncount_decrease != 0)),
            "feature_space_nCount_max_decrease": int(np.max(ncount_decrease)),
            "feature_space_percent_mt_changed_cell_count": int(
                np.count_nonzero(percent_mt_change > 1e-12)
            ),
            "feature_space_percent_mt_max_absolute_change": format_number(
                float(np.max(percent_mt_change))
            ),
        }
    else:
        feature_space_difference_summary = {
            "feature_space_nFeature_changed_cell_count": "",
            "feature_space_nFeature_max_decrease": "",
            "feature_space_nCount_changed_cell_count": "",
            "feature_space_nCount_max_decrease": "",
            "feature_space_percent_mt_changed_cell_count": "",
            "feature_space_percent_mt_max_absolute_change": "",
        }
    gene_mean_summary = summarize_gene_means(
        gene_means,
        rows,
        zero_value_count,
        total_values_checked,
    )
    return {
        "matrix_rows_genes": rows,
        "matrix_cols_cells": expected_cells,
        "header_total_columns_including_gene_col": len(header_fields),
        "first_header_fields": first_header_fields,
        "first_genes": "|".join(first_genes),
        "barcode_suffix_examples": "|".join(barcode_suffixes),
        "row_identity": (
            "gene_symbol_or_ensembl_like_name" if gene_name_issue_count == 0 else "mixed_or_invalid_gene_name"
        ),
        "column_identity": "cell_barcode",
        "first_header_blank": "true" if first_header_blank else "false",
        "gene_name_type": (
            "gene_symbol_or_ensembl_like_name" if gene_name_issue_count == 0 else "mixed_or_invalid_gene_name"
        ),
        "gene_name_issue_count": gene_name_issue_count,
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
        **public_space_qc_summary,
        **author_like_qc_summary,
        "author_cell_qc_reproduction_status": author_cell_qc_reproduction_status,
        "author_cell_qc_reproduction_note": author_cell_qc_reproduction_note,
        **feature_space_difference_summary,
        **author_feature_filter_summary,
        **gene_mean_summary,
        "mt_gene_count": mt_gene_count,
        "mt_gene_examples": "|".join(mt_genes),
        "hb_gene_count": hb_gene_count,
        "hb_gene_examples": "|".join(hb_genes),
        "gene_order_sha256": gene_hash.hexdigest().upper(),
        "observed_numeric_type": observed_numeric_type,
        "scan_seconds": round(time.perf_counter() - start, 2),
    }


def legacy_numeric_precheck(stats: Dict[str, object], precheck: Dict[str, str]) -> Tuple[str, List[str]]:
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
    sample_by_file = {str(row["sample_file"]): row for row in sample_info}
    precheck_by_member = {row["member_name"]: row for row in precheck_rows}
    mismatch_notes: List[str] = []
    output: List[Dict[str, object]] = []

    for manifest_row in processed_manifest:
        path = root / manifest_row["extracted_path"]
        stats = audit_csv_gz(path)
        member = manifest_row["archive_member_name"]
        sample = sample_by_file.get(member, {})
        precheck = precheck_by_member.get(member, {})
        row_mismatches: List[str] = []
        if not precheck:
            row_mismatches.append("precheck row missing for archive member")
        for field in [
            "matrix_rows_genes",
            "matrix_cols_cells",
            "mt_gene_count",
            "hb_gene_count",
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
        precheck_match = not row_mismatches
        pre_decision = precheck.get("audit_decision_precheck", "")
        processed_match = (
            sample.get("sample_file") == member
            and sample.get("geo_accession") == manifest_row.get("geo_accession")
            and sample.get("sample_id") == manifest_row.get("sample_id")
            and sample.get("sample_id_match_status") == "match"
        )
        if not processed_match:
            row_mismatches.append("processed manifest, GEO accession/title and sample_info mapping do not match")
        structural_failures: List[str] = []
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
        author_cell_qc_status = str(stats["author_cell_qc_reproduction_status"])
        if author_cell_qc_status == "not_evaluable":
            structural_failures.append("reported_author_cell_QC_thresholds_not_evaluable")
        feature_filter_status = str(stats["author_feature_filter_reproduction_status"])
        if feature_filter_status not in {"pass", "measured_mismatch"}:
            structural_failures.append("reported_author_feature_filter_not_evaluable")
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
            if feature_filter_status == "pass":
                feature_boundary_reason = (
                    "public feature rows satisfy the reported per-sample min.cells=3 boundary"
                )
            else:
                feature_boundary_reason = (
                    f"public feature rows include {stats['feature_rows_detected_lt_3_count']} genes detected "
                    "in fewer than three cells, so the author's per-sample feature filtering is not "
                    "embedded in this public file"
                )
            decision_reason = (
                "full-stream numeric/structure audit passed; per-cell nCount showed no preregistered "
                "normalization artifact; nFeature and percent.mt were recomputed in both the public full "
                "feature space and an author-like per-sample min.cells=3 feature space; sparse right-skew "
                "evidence was present; processed manifest and GEO sample mapping matched; "
                f"author_cell_qc_status={author_cell_qc_status}; {feature_boundary_reason}"
            )
        else:
            decision_reason = (
                "pause because one or more format, distribution, duplicate, mapping or precheck conditions failed: "
                + "; ".join(failure_details)
                + f"; artifact_evidence={stats['normalization_artifact_reason']}"
            )
        boundary_notes: List[str] = []
        if author_cell_qc_status == "measured_mismatch":
            boundary_notes.append(
                "nonblocking processing-boundary limitation: at least one retained public cell differs "
                "from an author-reported cell threshold in one or both explicitly measured feature spaces"
            )
        if feature_filter_status == "measured_mismatch":
            boundary_notes.append(
                "nonblocking processing-boundary limitation: public feature rows do not reflect the "
                "author-reported per-sample min.cells=3 filter"
            )
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
                "matrix_orientation": "gene_by_cell",
                **stats,
                "legacy_numeric_precheck_status": legacy_status,
                "precheck_audit_decision": pre_decision,
                "public_processing_evidence_status": (
                    "public_input_shape_verified_author_processing_boundaries_measured"
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
    if len(rows) != 40:
        raise RuntimeError(f"processed_input_manifest.tsv must contain 40 rows; observed {len(rows)}")
    seen: set[str] = set()
    for row in rows:
        member = row.get("archive_member_name", "")
        if not member or member in seen:
            raise RuntimeError(f"Duplicate or empty archive member in processed manifest: {member!r}")
        seen.add(member)
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
    require_paths(root, STAGE_REQUIRED, STAGE_NAME)
    manifest = read_tsv(root / "data/metadata/processed_input_manifest.tsv")
    validate_step1_manifest(root, manifest)
    append_log(root, f"F0 step2 started; run_id={current_run_id()}; validated_step1_files=40")

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

    precheck = read_tsv(root / "results/F0_audit/gse183904_csv_structure_precheck.tsv")
    data_audit, mismatches = build_data_audit(root, manifest, sample_info, precheck)
    append_log(
        root,
        "SCHEMA_MIGRATION precheck.suspected_matrix_type -> observed_numeric_type; "
        "precheck.audit_decision_precheck=enter_full_F1_candidate -> "
        "audit_decision=enter_full_F1_independent_reQC only after formal audit; "
        "suspected_matrix_type now records public input shape only",
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
        "per_cell_nCount_min",
        "per_cell_nCount_Q1",
        "per_cell_nCount_median",
        "per_cell_nCount_Q3",
        "per_cell_nCount_max",
        "public_full_feature_nCount_min",
        "public_full_feature_nCount_Q1",
        "public_full_feature_nCount_median",
        "public_full_feature_nCount_Q3",
        "public_full_feature_nCount_max",
        "public_full_feature_nFeature_min",
        "public_full_feature_nFeature_Q1",
        "public_full_feature_nFeature_median",
        "public_full_feature_nFeature_Q3",
        "public_full_feature_nFeature_max",
        "public_full_feature_percent_mt_min",
        "public_full_feature_percent_mt_Q1",
        "public_full_feature_percent_mt_median",
        "public_full_feature_percent_mt_Q3",
        "public_full_feature_percent_mt_max",
        "author_nFeature_lt_500_count_public_space",
        "author_nFeature_ge_6000_count_public_space",
        "author_percent_mt_gt_20_count_public_space",
        "author_cell_threshold_mismatch_count_public_space",
        "author_cell_qc_reproduction_status_public_space",
        "author_cell_qc_reproduction_note_public_space",
        "author_like_nCount_min",
        "author_like_nCount_Q1",
        "author_like_nCount_median",
        "author_like_nCount_Q3",
        "author_like_nCount_max",
        "author_like_nFeature_min",
        "author_like_nFeature_Q1",
        "author_like_nFeature_median",
        "author_like_nFeature_Q3",
        "author_like_nFeature_max",
        "author_like_percent_mt_min",
        "author_like_percent_mt_Q1",
        "author_like_percent_mt_median",
        "author_like_percent_mt_Q3",
        "author_like_percent_mt_max",
        "author_nFeature_lt_500_count_author_like_space",
        "author_nFeature_ge_6000_count_author_like_space",
        "author_percent_mt_gt_20_count_author_like_space",
        "author_cell_threshold_mismatch_count_author_like_space",
        "author_cell_qc_reproduction_status_author_like_space",
        "author_cell_qc_reproduction_note_author_like_space",
        "author_cell_qc_reproduction_status",
        "author_cell_qc_reproduction_note",
        "feature_space_nFeature_changed_cell_count",
        "feature_space_nFeature_max_decrease",
        "feature_space_nCount_changed_cell_count",
        "feature_space_nCount_max_decrease",
        "feature_space_percent_mt_changed_cell_count",
        "feature_space_percent_mt_max_absolute_change",
        "author_min_cells_per_feature_reported",
        "feature_rows_detected_in_0_cells",
        "feature_rows_detected_in_1_cell",
        "feature_rows_detected_in_2_cells",
        "feature_rows_detected_lt_3_count",
        "feature_rows_detected_lt_3_examples",
        "author_like_retained_feature_count",
        "author_feature_filter_reproduction_status",
        "author_feature_filter_reproduction_note",
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
        "hb_gene_count",
        "hb_gene_examples",
        "gene_order_sha256",
        "observed_numeric_type",
        "public_processing_evidence_status",
        "suspected_matrix_type",
        "scan_seconds",
        "legacy_numeric_precheck_status",
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

    unclear = [row for row in sample_info if row.get("group_analysis") == "Unclear"]
    artifacts = [row for row in data_audit if row.get("normalization_artifact_flag") != "false"]
    feature_filter_boundary_mismatches = [
        row
        for row in data_audit
        if row.get("author_feature_filter_reproduction_status") == "measured_mismatch"
    ]
    cell_qc_boundary_mismatches = [
        row
        for row in data_audit
        if row.get("author_cell_qc_reproduction_status") == "measured_mismatch"
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
        f"cell_qc_boundary_mismatches={len(cell_qc_boundary_mismatches)}; "
        f"feature_filter_boundary_mismatches={len(feature_filter_boundary_mismatches)}; "
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
    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, STAGE_REQUIRED, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
