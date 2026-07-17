#!/usr/bin/env python3
"""Read-only validation for the F0 fixed-QC implementation.

This script writes only inside a temporary directory. It checks the frozen
inequality boundaries, exact globin-panel matching, and the real sample1
regression target without creating any formal F0 output.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import tarfile
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np

from F0_step2_sample_info_and_audit import assess_fixed_qc, audit_csv_gz
from F0_step3_inventory_and_markers import build_author_processing_audit
from F0_step4_decisions_and_gate import (
    EXPECTED_GLOBIN_PANEL,
    PILOT_FILE_NAME,
    build_gate_checklist,
)
from f0_utils import read_tsv


def require_equal(observed: object, expected: object, label: str) -> None:
    if observed != expected:
        raise AssertionError(f"{label}: observed={observed!r}; expected={expected!r}")


def validate_fixed_boundaries() -> None:
    ncount = np.asarray([2000, 2000, 2000, 2000, 1000, 1001, 2000], dtype=np.int64)
    nfeature = np.asarray([499, 500, 5999, 6000, 550, 550, 550], dtype=np.int32)
    mt_ncount = np.asarray([0, 400, 0, 0, 0, 201, 0], dtype=np.int64)
    hb_ncount = np.asarray([0, 0, 0, 0, 0, 0, 100], dtype=np.int64)
    result = assess_fixed_qc(ncount, nfeature, mt_ncount, hb_ncount, complete=True)

    expected = {
        "fail_nFeature_low_count": 1,
        "fail_nFeature_high_count": 1,
        "fail_nCount_count": 1,
        "fail_percent_mt_count": 1,
        "fail_percent_hb_count": 1,
        "source_reported_qc_pass_count": 4,
        "additional_fail_nCount_after_source_count": 1,
        "additional_fail_percent_hb_after_source_nCount_count": 1,
        "final_fixed_qc_fail_count": 5,
        "final_fixed_qc_pass_count": 2,
        "fixed_qc_rule_recalculation_status": "pass",
    }
    for field, value in expected.items():
        require_equal(result[field], value, f"fixed-boundary {field}")


def validate_exact_globin_matching(temp_dir: Path) -> None:
    path = temp_dir / "synthetic_nonpilot.csv.gz"
    barcodes = [
        "AAAAAAAAAAAAAAAA_1",
        "CCCCCCCCCCCCCCCC_1",
        "GGGGGGGGGGGGGGGG_1",
    ]
    rows = [
        ("HBA1", [1, 1, 1]),
        ("HBEGF", [2, 2, 2]),
        ("HBP1", [3, 3, 3]),
        ("HBS1L", [4, 4, 4]),
        ("MT-ND1", [1, 1, 1]),
        ("GENE1", [1000, 1000, 1000]),
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("," + ",".join(barcodes) + "\n")
        for gene, values in rows:
            handle.write(gene + "," + ",".join(str(value) for value in values) + "\n")

    result = audit_csv_gz(path)
    require_equal(result["globin_panel_used_for_qc"], "HBA1", "exact globin inclusion")
    require_equal(
        result["false_hb_prefix_genes_excluded"],
        "HBEGF|HBP1|HBS1L",
        "false HB-prefix exclusion",
    )
    require_equal(result["pilot_validation_status"], "not_applicable", "synthetic pilot state")


def validate_real_sample1(project_root: Path, temp_dir: Path) -> None:
    archive = project_root / "data/public_downloads/GSE183904_RAW.tar"
    if not archive.exists():
        raise FileNotFoundError(f"Required pilot archive is missing: {archive}")

    target_name = "GSM5573466_sample1.csv.gz"
    with tarfile.open(archive, "r") as tar:
        members = [member for member in tar.getmembers() if Path(member.name).name == target_name]
        require_equal(len(members), 1, "sample1 archive-member count")
        source = tar.extractfile(members[0])
        if source is None:
            raise RuntimeError("sample1 archive member is not a readable file")
        target = temp_dir / target_name
        with target.open("wb") as handle:
            shutil.copyfileobj(source, handle)

    result = audit_csv_gz(target)
    require_equal(result["pilot_validation_status"], "pass", "real sample1 regression")
    require_equal(result["qc_retained_feature_count"], 19294, "sample1 working features")
    require_equal(result["source_reported_qc_pass_count"], 2684, "sample1 source-rule pass")
    require_equal(result["final_fixed_qc_pass_count"], 2631, "sample1 final fixed-QC pass")


def validate_downstream_gate_contract(project_root: Path) -> None:
    audits = []
    samples = []
    processed_manifest = []
    for index in range(40):
        pilot = index == 0
        file_name = PILOT_FILE_NAME if pilot else f"synthetic_sample{index + 1}.csv.gz"
        audits.append(
            {
                "file_name": file_name,
                "include_in_f1": "true",
                "matrix_cols_cells": "2685",
                "source_reported_qc_pass_count": "2684",
                "final_fixed_qc_pass_count": "2631",
                "final_fixed_qc_fail_count": "54",
                "fail_nFeature_low_count": "0",
                "fail_nFeature_high_count": "0",
                "fail_nCount_count": "53",
                "fail_percent_mt_count": "1",
                "fail_percent_hb_count": "0",
                "fixed_qc_rule_recalculation_status": "pass",
                "working_feature_space_recalculation_status": "pass",
                "pilot_validation_applicable": "true" if pilot else "false",
                "pilot_validation_status": "pass" if pilot else "not_applicable",
                "additional_fail_nCount_after_source_count": "53",
                "additional_fail_percent_hb_after_source_nCount_count": "0",
                "qc_retained_feature_count": "19294",
                "feature_rows_detected_lt_3_count": "7277",
                "globin_panel_expected": EXPECTED_GLOBIN_PANEL,
                "globin_panel_present": "HBA1",
                "globin_panel_used_for_qc": "HBA1|HBA2|HBB|HBD" if pilot else "HBA1",
                "audit_decision": "enter_full_F1_independent_reQC",
                "normalization_artifact_flag": "false",
                "observed_numeric_type": "nonnegative_integer_count_like",
                "suspected_matrix_type": "public_called_cell_raw_gene_count_matrix",
                "processed_input_manifest_match": "true",
                "per_gene_mean_distribution_status": "consistent_with_sparse_right_skew",
                "public_processing_evidence_status": (
                    "public_input_shape_verified_fixed_QC_recalculated_"
                    "processing_history_pending_F0_step3"
                ),
                "raw_droplet_available": "false",
                "empty_droplet_background_available": "false",
                "format_decision_scope": "file_format_only",
                "decision_scope": "file_format_and_public_processing_boundary",
                "precheck_comparison_status": "match",
            }
        )
        samples.append(
            {
                "sample_file": file_name,
                "group_analysis": "Normal_Gastric",
                "sample_id_match_status": "match",
                "source_of_group": "GEO",
                "metadata_confidence": "high",
                "include_in_f1": "true",
                "include_in_group_comparison": "true",
            }
        )
        processed_manifest.append({"sha256": "A" * 64, "file_size": "1"})

    source_rows = read_tsv(
        project_root / "docs/source_verification/GSE183904_processing_history_source_audit.tsv"
    )
    history_rows = build_author_processing_audit(source_rows, audits)
    gate_rows = build_gate_checklist(
        project_root,
        samples,
        audits,
        processed_manifest,
        history_rows,
    )
    require_equal(len(gate_rows), 10, "grouped gate-item count")
    status_by_item = {row["gate_item"]: row["pass_fail"] for row in gate_rows}
    for item in (
        "sample_info",
        "data_audit",
        "precheck_comparison",
        "fixed_QC_rule_recalculation",
        "working_feature_space_recalculation",
    ):
        require_equal(status_by_item[item], "PASS", f"synthetic gate {item}")
    require_equal(
        status_by_item["author_processing_provenance"],
        "PASS_WITH_NOTED_ISSUES",
        "synthetic provenance limitation state",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    args = parser.parse_args(argv)
    project_root = Path(args.project_root).resolve()

    validate_fixed_boundaries()
    with tempfile.TemporaryDirectory(prefix="f0_readonly_validation_") as temp_name:
        temp_dir = Path(temp_name)
        validate_exact_globin_matching(temp_dir)
        validate_real_sample1(project_root, temp_dir)
    validate_downstream_gate_contract(project_root)

    print("F0 read-only validation: PASS")
    print("- fixed inequality boundaries: PASS")
    print("- exact globin-panel matching: PASS")
    print("- real sample1 frozen regression: PASS")
    print("- Step2-to-Step4 field and gate contract: PASS")
    print("No formal F0 outputs were written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
