#!/usr/bin/env python3
"""只读验证 F0 固定 QC 的代码实现是否忠实于方案。

这个脚本相当于正式执行前的“单元考试”，重点检查：
1. Step1 是否只登记“预期方向、待 Step2 验证”；
2. Step2 能否接受 gene × cell 并拒绝 cell × gene；
3. 每个 QC 不等式在等号边界处是否写对；
4. HB 百分比是否只使用冻结 globin panel，而没有误收 HBEGF 等基因；
5. 真实 sample1 是否复现已冻结的基准数字；
6. Step2 的字段能否顺利传到 Step3 和 Step4 的十项 gate；
7. F0 独立 Python 环境文件是否与当前获批版本一致。

脚本只在系统临时目录中写合成测试文件，结束后自动删除，不会生成任何正式
F0 输出，也不会删除真实细胞。
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np

from F0_step1_structure_and_extract import build_processed_manifest
from F0_step2_sample_info_and_audit import (
    assess_fixed_qc,
    audit_csv_gz,
    validate_step1_manifest,
)
from F0_step3_inventory_and_markers import (
    F0_ENVIRONMENT_LOCK_PATHS,
    build_author_processing_audit,
    build_file_manifest,
)
from F0_step4_decisions_and_gate import (
    EXPECTED_GLOBIN_PANEL,
    PILOT_FILE_NAME,
    build_gate_checklist,
)
from f0_utils import read_tsv, validate_f0_python_environment


def require_equal(observed: object, expected: object, label: str) -> None:
    """比较实测值与预期值；不一致时给出容易定位的错误标签。"""

    if observed != expected:
        raise AssertionError(f"{label}: observed={observed!r}; expected={expected!r}")


def validate_f0_environment_lock(project_root: Path) -> None:
    """确认 F0 的 pip/Conda 两套替代规格与实际获批运行时一致。"""

    validate_f0_python_environment(
        project_root,
        actual_python_version=sys.version.split()[0],
        actual_numpy_version=np.__version__,
    )

    lock_registry = {
        row.get("lock_item", ""): row
        for row in read_tsv(project_root / "environment/environment_lock_manifest.tsv")
    }
    for lock_item in ("F0_pip_requirements", "F0_conda_environment"):
        require_equal(lock_registry.get(lock_item, {}).get("lock_status"), "recorded", lock_item)

    # 在内存中生成 file manifest，确认两个环境文件都会得到大写 SHA256。
    manifest_rows = build_file_manifest(project_root, [], [])
    manifest_by_path = {
        row.get("relative_path_if_available", ""): row for row in manifest_rows
    }
    for lock_path in F0_ENVIRONMENT_LOCK_PATHS:
        checksum = manifest_by_path.get(lock_path, {}).get("sha256", "")
        require_equal(len(checksum), 64, f"{lock_path} SHA256 length")
        require_equal(checksum, checksum.upper(), f"{lock_path} SHA256 uppercase")


def validate_step1_pending_orientation_contract(temp_dir: Path) -> None:
    """确认 Step1 manifest 只写预期方向，不提前声称已经验证。"""

    root = temp_dir / "step1_contract"
    matrix_dir = root / "data/processed_input/GSE183904"
    matrix_dir.mkdir(parents=True)
    members = []
    for index in range(1, 41):
        name = f"GSM{index:07d}_sample{index}.csv.gz"
        (matrix_dir / name).write_bytes(
            f"temporary orientation-contract test {index}".encode("ascii")
        )
        members.append(tarfile.TarInfo(name=name))
    manifest_rows = build_processed_manifest(root, members)
    validate_step1_manifest(root, manifest_rows)
    row = manifest_rows[0]
    require_equal(
        row["file_role"],
        "expected_gene_by_cell_matrix_pending_validation",
        "Step1 pending-orientation file role",
    )
    require_equal(row["expected_matrix_orientation"], "gene_by_cell", "Step1 expected orientation")
    require_equal(
        row["orientation_validation_status"],
        "pending_F0_step2_full_stream_validation",
        "Step1 pending orientation status",
    )


def validate_fixed_boundaries() -> None:
    """用人工构造的 7 个细胞检查所有 QC 等号边界。

    例如 nFeature=500 应保留、nFeature=6000 应排除；percent.mt=20 应
    保留，而 percent.HB=5 应排除。这样可防止 ``<``、``<=`` 写反。
    """

    # 每个位置代表一个合成细胞，四个数组分别提供该细胞的 QC 原始量。
    ncount = np.asarray([2000, 2000, 2000, 2000, 1000, 1001, 2000], dtype=np.int64)
    nfeature = np.asarray([499, 500, 5999, 6000, 550, 550, 550], dtype=np.int32)
    mt_ncount = np.asarray([0, 400, 0, 0, 0, 201, 0], dtype=np.int64)
    hb_ncount = np.asarray([0, 0, 0, 0, 0, 0, 100], dtype=np.int64)
    result = assess_fixed_qc(ncount, nfeature, mt_ncount, hb_ncount, complete=True)

    # 这些预期计数由冻结规则人工推导，不从正式脚本动态生成。
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
    """构造小矩阵，确认只有 HBA1 被计入，三个假 ``HB`` 前缀基因被排除。"""

    path = temp_dir / "synthetic_nonpilot.csv.gz"
    barcodes = [
        "AAAAAAAAAAAAAAAA_1",
        "CCCCCCCCCCCCCCCC.1_1",
        "GGGGGGGGGGGGGGGG-1_1",
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
    require_equal(result["matrix_orientation"], "gene_by_cell", "gene-by-cell orientation")
    require_equal(
        result["matrix_orientation_validation_status"],
        "pass_gene_by_cell",
        "gene-by-cell validation state",
    )
    require_equal(result["globin_panel_used_for_qc"], "HBA1", "exact globin inclusion")
    require_equal(
        result["false_hb_prefix_genes_excluded"],
        "HBEGF|HBP1|HBS1L",
        "false HB-prefix exclusion",
    )
    require_equal(result["pilot_validation_status"], "not_applicable", "synthetic pilot state")


def validate_reversed_orientation_rejected(temp_dir: Path) -> None:
    """构造 cell × gene 矩阵，确认 Step2 不会把它误报为 gene × cell。"""

    path = temp_dir / "synthetic_cell_by_gene.csv.gz"
    genes = ["HBA1", "GENE1", "MT-ND1"]
    barcode_rows = [
        ("AAAAAAAAAAAAAAAA_1", [1, 2, 3]),
        ("CCCCCCCCCCCCCCCC_1", [4, 5, 6]),
        ("GGGGGGGGGGGGGGGG_1", [7, 8, 9]),
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("," + ",".join(genes) + "\n")
        for barcode, values in barcode_rows:
            handle.write(barcode + "," + ",".join(str(value) for value in values) + "\n")

    result = audit_csv_gz(path)
    require_equal(result["matrix_orientation"], "not_confirmed", "cell-by-gene rejection")
    require_equal(
        result["matrix_orientation_validation_status"],
        "fail_not_gene_by_cell",
        "cell-by-gene validation state",
    )
    require_equal(result["row_label_cell_barcode_like_count"], 3, "barcode-like row labels")


def validate_real_sample1(project_root: Path, temp_dir: Path) -> None:
    """从真实压缩包临时提取 sample1，并复现冻结的 pilot 数字。"""

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

    # 调用与正式 Step2 完全相同的审计函数，避免测试和正式代码各算一套。
    result = audit_csv_gz(target)
    require_equal(result["pilot_validation_status"], "pass", "real sample1 regression")
    require_equal(result["matrix_orientation"], "gene_by_cell", "sample1 matrix orientation")
    require_equal(
        result["matrix_orientation_validation_status"],
        "pass_gene_by_cell",
        "sample1 orientation validation",
    )
    require_equal(result["qc_retained_feature_count"], 19294, "sample1 working features")
    require_equal(result["source_reported_qc_pass_count"], 2684, "sample1 source-rule pass")
    require_equal(result["final_fixed_qc_pass_count"], 2631, "sample1 final fixed-QC pass")


def validate_downstream_gate_contract(project_root: Path) -> None:
    """用 40 个合成样本验证 Step2→Step3→Step4 的字段和 gate 契约。

    这里不冒充真实 F0 结果，只测试字段名、状态值和 gate 判断能否闭合。
    如果上游字段被改名而下游忘记同步，本测试会立即失败。
    """

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
                "matrix_orientation": "gene_by_cell",
                "matrix_orientation_validation_status": "pass_gene_by_cell",
                "row_label_cell_barcode_like_count": "0",
                "cell_barcode_pattern": "10x_16nt_barcode_with_numeric_suffix",
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

    # 处理史来源表使用真实只读文件；表达审计部分使用上面构造的合成结果。
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
    """依次运行四类只读验证；任何断言失败都会以非 0 状态退出。"""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    args = parser.parse_args(argv)
    project_root = Path(args.project_root).resolve()

    validate_f0_environment_lock(project_root)
    validate_fixed_boundaries()
    with tempfile.TemporaryDirectory(prefix="f0_readonly_validation_") as temp_name:
        temp_dir = Path(temp_name)
        validate_step1_pending_orientation_contract(temp_dir)
        validate_exact_globin_matching(temp_dir)
        validate_reversed_orientation_rejected(temp_dir)
        validate_real_sample1(project_root, temp_dir)
    validate_downstream_gate_contract(project_root)

    print("F0 read-only validation: PASS")
    print("- F0 Python environment locks: PASS")
    print("- Step1 pending-orientation contract: PASS")
    print("- gene-by-cell accepted / cell-by-gene rejected: PASS")
    print("- fixed inequality boundaries: PASS")
    print("- exact globin-panel matching: PASS")
    print("- real sample1 frozen regression: PASS")
    print("- Step2-to-Step4 field and gate contract: PASS")
    print("No formal F0 outputs were written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
