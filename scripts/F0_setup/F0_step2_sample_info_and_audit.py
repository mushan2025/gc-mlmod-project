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
import re
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

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
    manifest_rows: Sequence[Dict[str, str]],
) -> List[Dict[str, object]]:
    accessions = sample_fields["!Sample_geo_accession"]
    titles = sample_fields["!Sample_title"]
    sources = sample_fields.get("!Sample_source_name_ch1", [""] * len(accessions))
    chars = sample_fields.get("!Sample_characteristics_ch1", [""] * len(accessions))
    member_by_accession = {
        row.get("geo_accession", ""): row.get("archive_member_name", "") for row in manifest_rows
    }

    rows: List[Dict[str, object]] = []
    for idx, accession in enumerate(accessions):
        title = titles[idx]
        sample_match = re.search(r"(sample\d+)", title)
        sample_id = sample_match.group(1) if sample_match else f"sample{idx + 1}"
        group, tissue_site, tumor_status = group_from_title(title)
        member_name = member_by_accession.get(accession, f"{accession}_{sample_id}.csv.gz")
        rows.append(
            {
                "geo_accession": accession,
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
                "include_in_group_comparison": "false" if group == "Normal_Peritoneum" else "true",
                "is_paired": "unknown",
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
    return rows


def classify_anomalous_values(values: str) -> Dict[str, int]:
    counts = {
        "missing_value_count": 0,
        "noninteger_float_count": 0,
        "invalid_nonnumeric_count": 0,
        "negative_integer_count": 0,
    }
    for raw_value in values.split(","):
        value = clean_csv_token(raw_value)
        if value == "":
            counts["missing_value_count"] += 1
            continue
        try:
            integer = int(value)
        except ValueError:
            try:
                float(value)
            except ValueError:
                counts["invalid_nonnumeric_count"] += 1
            else:
                counts["noninteger_float_count"] += 1
        else:
            if integer < 0:
                counts["negative_integer_count"] += 1
    return counts


def audit_csv_gz(path: Path) -> Dict[str, object]:
    start = time.perf_counter()
    rows = 0
    bad_column_count_rows = 0
    missing_value_count = 0
    noninteger_float_count = 0
    invalid_nonnumeric_count = 0
    negative_integer_count = 0
    total_values_checked = 0
    mt_gene_count = 0
    hb_gene_count = 0
    mt_genes: List[str] = []
    hb_genes: List[str] = []
    first_genes: List[str] = []
    gene_hash = hashlib.sha256()

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline()
        if not header:
            raise RuntimeError(f"Empty matrix file: {path}")
        header_fields = [clean_csv_token(value) for value in header.rstrip("\r\n").split(",")]
        expected_cells = len(header_fields) - 1
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
            if observed_cells != expected_cells:
                bad_column_count_rows += 1
            total_values_checked += observed_cells

            if not NONNEGATIVE_INTEGER_ROW.fullmatch(values):
                anomalies = classify_anomalous_values(values)
                missing_value_count += anomalies["missing_value_count"]
                noninteger_float_count += anomalies["noninteger_float_count"]
                invalid_nonnumeric_count += anomalies["invalid_nonnumeric_count"]
                negative_integer_count += anomalies["negative_integer_count"]

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

    numeric_anomaly_count = (
        missing_value_count
        + noninteger_float_count
        + invalid_nonnumeric_count
        + negative_integer_count
    )
    suspected = (
        "nonnegative_integer_count_like"
        if bad_column_count_rows == 0 and numeric_anomaly_count == 0
        else "format_or_numeric_issue_detected"
    )
    return {
        "matrix_rows_genes": rows,
        "matrix_cols_cells": expected_cells,
        "header_total_columns_including_gene_col": len(header_fields),
        "first_header_fields": first_header_fields,
        "first_genes": "|".join(first_genes),
        "barcode_suffix_examples": "|".join(barcode_suffixes),
        "bad_column_count_rows": bad_column_count_rows,
        "integer_check_method": "full_stream_two_level_exact",
        "total_values_checked": total_values_checked,
        "missing_value_count": missing_value_count,
        "noninteger_float_count": noninteger_float_count,
        "invalid_nonnumeric_count": invalid_nonnumeric_count,
        "negative_integer_count": negative_integer_count,
        "numeric_anomaly_count": numeric_anomaly_count,
        "mt_gene_count": mt_gene_count,
        "mt_gene_examples": "|".join(mt_genes),
        "hb_gene_count": hb_gene_count,
        "hb_gene_examples": "|".join(hb_genes),
        "gene_order_sha256": gene_hash.hexdigest().upper(),
        "suspected_matrix_type": suspected,
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
            "suspected_matrix_type",
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

        legacy_status, legacy_mismatches = legacy_numeric_precheck(stats, precheck)
        row_mismatches.extend(legacy_mismatches)
        pre_decision = precheck.get("audit_decision_precheck", "")
        ok = (
            stats["suspected_matrix_type"] == "nonnegative_integer_count_like"
            and not row_mismatches
            and pre_decision == "enter_full_F1_candidate"
            and sample.get("include_in_f1") == "true"
        )
        if row_mismatches:
            mismatch_notes.append(f"{member} | " + " | ".join(row_mismatches))
        output.append(
            {
                "geo_accession": manifest_row["geo_accession"],
                "sample_id": manifest_row["sample_id"],
                "sample_file": member,
                "extracted_path": manifest_row["extracted_path"],
                "matrix_orientation": "gene_rows_by_cell_columns",
                **stats,
                "legacy_numeric_precheck_status": legacy_status,
                "precheck_audit_decision": pre_decision,
                "audit_decision": "enter_full_F1" if ok else "pause_for_review",
                "decision_scope": "file_format_only",
                "include_in_f1": sample.get("include_in_f1", "pending"),
                "raw_droplet_available": "false",
                "empty_droplet_background_available": "false",
                "normalization_artifact_flag": (
                    "false" if stats["suspected_matrix_type"] == "nonnegative_integer_count_like" else "true"
                ),
                "precheck_comparison_status": "match" if not row_mismatches else "mismatch",
                "note": "; ".join(row_mismatches),
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
    sample_info = build_sample_info(sample_fields, patient_by_sample, manifest)
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
    write_tsv(root / "data/metadata/sample_info.tsv", sample_info, sample_fields_out)

    precheck = read_tsv(root / "results/F0_audit/gse183904_csv_structure_precheck.tsv")
    data_audit, mismatches = build_data_audit(root, manifest, sample_info, precheck)
    data_fields = [
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
        "integer_check_method",
        "total_values_checked",
        "missing_value_count",
        "noninteger_float_count",
        "invalid_nonnumeric_count",
        "negative_integer_count",
        "numeric_anomaly_count",
        "mt_gene_count",
        "mt_gene_examples",
        "hb_gene_count",
        "hb_gene_examples",
        "gene_order_sha256",
        "suspected_matrix_type",
        "scan_seconds",
        "legacy_numeric_precheck_status",
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
    write_tsv(root / "data/metadata/data_audit.tsv", data_audit, data_fields)
    for mismatch in mismatches:
        append_log(root, "PRECHECK_MISMATCH " + mismatch)

    unclear = [row for row in sample_info if row.get("group_analysis") == "Unclear"]
    artifacts = [row for row in data_audit if row.get("normalization_artifact_flag") == "true"]
    paused = [row for row in data_audit if row.get("audit_decision") != "enter_full_F1"]
    append_log(
        root,
        f"F0 step2 completed; sample_rows={len(sample_info)}; audit_rows={len(data_audit)}; "
        f"mismatches={len(mismatches)}; unclear_groups={len(unclear)}; artifacts={len(artifacts)}; paused={len(paused)}",
    )
    if len(sample_info) != 40 or unclear or mismatches or artifacts or paused:
        raise RuntimeError(
            "F0 step2 pause condition reached: "
            f"sample_rows={len(sample_info)}, unclear={len(unclear)}, "
            f"mismatches={len(mismatches)}, normalization_artifacts={len(artifacts)}, paused={len(paused)}"
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
