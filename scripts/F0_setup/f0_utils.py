#!/usr/bin/env python3
"""F0 分阶段审计脚本共用的基础工具和文件契约。

这个脚本本身不做生物学分析，主要负责三类公共工作：
1. 集中登记 F0 必需输入、正式输出和目录，避免不同阶段使用不同文件名；
2. 安全读写 TSV、计算 SHA256、追加运行日志；
3. 为每个阶段提供统一的命令行参数和只读 dry run。

可复现性约定：F0 写出的 SHA256 全部使用大写。比较校验和时也会先统一
转换为大写，因此已有清单中的大小写差异不会造成假的校验失败。
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


# F0 本身没有随机算法，但仍登记统一随机种子，便于后续阶段沿用和审计。
RANDOM_SEED = 42
RUN_ID = dt.datetime.now(dt.timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S")

# F0 正式运行必须生成的完整输出契约。Step4 会检查这些文件是否齐全。
F0_OUTPUTS = [
    "data/metadata/project_structure_ready.txt",
    "logs/F0_setup/analysis_log.md",
    "data/metadata/sample_info.tsv",
    "data/metadata/data_audit.tsv",
    "data/metadata/processed_input_manifest.tsv",
    "data/metadata/F0_dataset_inventory.tsv",
    "data/metadata/F0_file_manifest.tsv",
    "data/metadata/F0_metadata_field_inventory.tsv",
    "data/metadata/F0_author_processing_audit.tsv",
    "data/metadata/F0_data_readiness_by_F_section.tsv",
    "data/metadata/F0_external_resource_inventory.tsv",
    "data/metadata/F0_method_prior_decision.tsv",
    "data/metadata/decision_evidence_log.tsv",
    "data/metadata/excluded_samples.tsv",
    "results/F0_audit/F0_gate_checklist.tsv",
    "results/F0_audit/F0_global_data_reconnaissance_report.md",
    "results/F0_audit/F0_execution_report.md",
]

# 项目标准目录。Step1 只会补建缺失目录，不会删除已有文件。
REQUIRED_DIRS = [
    "docs",
    "scripts/environment_setup",
    "logs/environment_setup",
    "reports/environment_setup",
    "environment",
    "data/public_downloads",
    "data/public_downloads/SCENIC_resources",
    "data/processed_input/GSE183904",
    "data/metadata",
    "scripts/F0_setup",
    "scripts/F1_single_cell",
    "scripts/F2_scoring",
    "scripts/F3_function",
    "scripts/F4_communication",
    "scripts/F5_bulk_marker",
    "scripts/F6_immunity",
    "scripts/F7_genomics_pathway",
    "scripts/F8_single_gene",
    "results/F0_audit",
    "results/F1_qc",
    "objects/F1_single_cell",
    "results/F1_annotation",
    "results/F1_malignancy",
    "results/F2_scoring",
    "results/F3_function",
    "results/F4_communication",
    "results/F5_bulk_marker",
    "results/F6_immunity",
    "results/F7_genomics_pathway",
    "results/F8_single_gene",
    "logs/F0_setup",
    "logs/F1_single_cell",
    "logs/F2_scoring",
    "logs/F3_function",
    "logs/F4_communication",
    "logs/F5_bulk_marker",
    "logs/F6_immunity",
    "logs/F7_genomics_pathway",
    "logs/F8_single_gene",
    "tmp",
]

# 启动 F0 前必须存在的输入。缺少任一文件时，dry run 和正式执行都会停止。
REQUIRED_INPUTS = [
    "data/public_downloads/GSE183904_RAW.tar",
    "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
    "data/metadata/download_manifest.tsv",
    "data/metadata/preupload_resources_manifest.tsv",
    "data/metadata/preupload_pending_resources.tsv",
    "data/metadata/cell_type_marker_panel.tsv",
    "data/metadata/pipeline_parameters.yaml",
    "data/metadata/software_versions.tsv",
    "environment/execution_environment_inventory.tsv",
    "environment/environment_lock_manifest.tsv",
    "environment/random_seed_registry.tsv",
    "results/F0_audit/gse183904_csv_structure_precheck.tsv",
    "results/F0_audit/non_GSE183904_data_structure_precheck.tsv",
    "results/F0_audit/predownloaded_resource_structure_audit.tsv",
    "results/F0_audit/bulk_GEO_series_matrix_deep_precheck.tsv",
    "results/F0_audit/bulk_GEO_sample_overlap_precheck.tsv",
    "results/F0_audit/bulk_GEO_processing_metadata_precheck.tsv",
    "results/F0_audit/TCGA_STAD_survival_table_precheck.tsv",
    "results/F0_audit/TCGA_STAD_star_counts_inverse_validation.tsv",
    "results/F0_audit/GSE206785_dataset_structure_precheck.tsv",
    "results/F0_audit/GSE206785_metadata_precheck.tsv",
    "docs/source_verification/GSE183904_processing_history_source_audit.tsv",
]


def rel(path: Path, root: Path) -> str:
    """把绝对路径转换为相对项目根目录的路径，便于项目整体迁移。"""

    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def now_iso() -> str:
    """返回带时区的当前时间，用于日志和 manifest。"""

    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def normalize_sha256(value: object) -> str:
    """把任意 SHA256 文本去除首尾空白并转为大写。"""
    return str(value or "").strip().upper()


def sha256_equal(left: object, right: object) -> bool:
    """比较两个 SHA256；只忽略字母大小写，不忽略真正的内容差异。"""

    return normalize_sha256(left) == normalize_sha256(right)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """分块读取文件并计算大写 SHA256，避免大文件一次性占满内存。"""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    """读取制表符分隔表，每一行返回为“列名: 内容”的字典。"""

    # utf-8-sig 既能读取普通 UTF-8，也能自动去除 Windows 导出文件开头的 BOM。
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fields: Sequence[str]) -> None:
    """按固定字段顺序安全写出 TSV。

    先写入同目录临时文件，再替换正式文件。这样即使写入中途失败，也不容易
    留下一个看似存在但内容残缺的正式结果。
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: str(row.get(field, "")) for field in fields})
    tmp.replace(path)


def clean_csv_token(value: str) -> str:
    """清理 CSV 单元格外围引号，并还原 CSV 中转义的双引号。"""

    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1].replace('""', '"')
    return value


def append_log(root: Path, message: str) -> None:
    """向 F0 分析日志追加一条带时间戳的记录。"""

    log_path = root / "logs/F0_setup/analysis_log.md"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"- {now_iso()} | {message}\n")


def require_paths(root: Path, paths: Sequence[str], stage_name: str) -> None:
    """检查某阶段的必需文件；缺失时立即报错，禁止带病继续运行。"""

    missing = [path for path in paths if not (root / path).exists()]
    if missing:
        raise RuntimeError(f"{stage_name} missing required paths: " + ", ".join(missing))


def dry_run_report(
    root: Path,
    stage_name: str,
    required_paths: Sequence[str],
    planned_outputs: Sequence[str],
) -> int:
    """只检查输入并展示预期输出，不写任何正式分析文件。

    返回值 0 表示输入齐全；返回值 1 表示存在缺失输入。这个返回值可以被
    PowerShell 或其他调度工具识别。
    """

    missing = [path for path in required_paths if not (root / path).exists()]
    print(f"{stage_name} dry run: no outputs will be written.")
    print(f"Project root: {root.resolve()}")
    print(f"Required inputs: {len(required_paths)}")
    if missing:
        print("Missing inputs:")
        for path in missing:
            print(f"  - {path}")
    else:
        print("All required inputs are present.")
    print("Planned formal outputs:")
    for path in planned_outputs:
        print(f"  - {path}")
    return 1 if missing else 0


def parse_stage_args(description: str, argv: Sequence[str] | None = None) -> argparse.Namespace:
    """为所有 F0 脚本提供统一的 ``--project-root`` 和 ``--execute`` 参数。"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--project-root", default=".", help="Project root; defaults to current directory.")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Write this stage's formal outputs; omit for a read-only dry run.",
    )
    return parser.parse_args(argv)


def current_run_id() -> str:
    """取得本次运行编号；外部指定 F0_RUN_ID 时优先使用该编号。"""

    return os.environ.get("F0_RUN_ID", RUN_ID)
