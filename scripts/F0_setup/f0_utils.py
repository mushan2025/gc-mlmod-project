#!/usr/bin/env python3
"""Shared utilities and contracts for the staged F0 audit scripts.

All SHA256 values written by F0 are uppercase. Comparisons normalize both
sides to uppercase so manifests remain case-insensitive as required.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


RANDOM_SEED = 42
RUN_ID = dt.datetime.now(dt.timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S")

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

REQUIRED_INPUTS = [
    "data/public_downloads/GSE183904_RAW.tar",
    "data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz",
    "docs/source_verification/GSE183904_processing_history_source_audit.tsv",
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
]


def rel(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def normalize_sha256(value: object) -> str:
    """Return a normalized uppercase SHA256 string for storage/comparison."""
    return str(value or "").strip().upper()


def sha256_equal(left: object, right: object) -> bool:
    return normalize_sha256(left) == normalize_sha256(right)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    # utf-8-sig accepts ordinary UTF-8 and strips a BOM from Windows exports.
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fields: Sequence[str]) -> None:
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
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1].replace('""', '"')
    return value


def append_log(root: Path, message: str) -> None:
    log_path = root / "logs/F0_setup/analysis_log.md"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"- {now_iso()} | {message}\n")


def require_paths(root: Path, paths: Sequence[str], stage_name: str) -> None:
    missing = [path for path in paths if not (root / path).exists()]
    if missing:
        raise RuntimeError(f"{stage_name} missing required paths: " + ", ".join(missing))


def dry_run_report(
    root: Path,
    stage_name: str,
    required_paths: Sequence[str],
    planned_outputs: Sequence[str],
) -> int:
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
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--project-root", default=".", help="Project root; defaults to current directory.")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Write this stage's formal outputs; omit for a read-only dry run.",
    )
    return parser.parse_args(argv)


def current_run_id() -> str:
    return os.environ.get("F0_RUN_ID", RUN_ID)
