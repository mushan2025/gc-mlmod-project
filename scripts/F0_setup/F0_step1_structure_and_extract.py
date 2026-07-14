#!/usr/bin/env python3
"""F0 step 1: verify structure/archive, extract 40 matrices, and manifest them.

Dependencies: the 21 preregistered F0 inputs, especially GSE183904_RAW.tar and
download_manifest.tsv. Outputs: project_structure_ready.txt,
processed_input_manifest.tsv, and the initialized analysis_log.md.
"""

from __future__ import annotations

import platform
import shutil
import sys
import tarfile
from pathlib import Path
from typing import Dict, List, Sequence

from f0_utils import (
    REQUIRED_DIRS,
    REQUIRED_INPUTS,
    append_log,
    current_run_id,
    dry_run_report,
    normalize_sha256,
    now_iso,
    parse_stage_args,
    read_tsv,
    rel,
    require_paths,
    sha256_equal,
    sha256_file,
    write_tsv,
)


STAGE_NAME = "F0 step 1 structure and extract"
STAGE_OUTPUTS = [
    "data/metadata/project_structure_ready.txt",
    "data/metadata/processed_input_manifest.tsv",
    "logs/F0_setup/analysis_log.md",
]


def write_project_structure(root: Path) -> None:
    for directory in REQUIRED_DIRS:
        (root / directory).mkdir(parents=True, exist_ok=True)
    path = root / "data/metadata/project_structure_ready.txt"
    lines = [
        f"project_root={root.resolve()}",
        f"created_or_verified_at={now_iso()}",
        "operator=Codex",
        "directory_list=" + "|".join(REQUIRED_DIRS),
        "note=Directories were created if missing; existing inputs were not overwritten.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def initialize_log(root: Path) -> None:
    path = root / "logs/F0_setup/analysis_log.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.stat().st_size == 0:
        path.write_text("# F0 Analysis Log\n\n", encoding="utf-8")


def tar_csv_members(tar_path: Path) -> List[tarfile.TarInfo]:
    with tarfile.open(tar_path, "r") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and Path(member.name).name.endswith(".csv.gz")
        ]
    return sorted(members, key=lambda member: Path(member.name).name)


def extract_members(tar_path: Path, members: Sequence[tarfile.TarInfo], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    resolved_out = out_dir.resolve()
    with tarfile.open(tar_path, "r") as archive:
        for member in members:
            target = out_dir / Path(member.name).name
            if not target.resolve().parent.samefile(resolved_out):
                raise RuntimeError(f"Unsafe tar member target: {member.name}")
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError(f"Could not read tar member: {member.name}")
            tmp = target.with_suffix(target.suffix + ".tmp")
            with source, tmp.open("wb") as destination:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
            tmp.replace(target)


def archive_expected_sha(root: Path) -> str:
    rows = read_tsv(root / "data/metadata/download_manifest.tsv")
    row = next((item for item in rows if item.get("file_name") == "GSE183904_RAW.tar"), None)
    if row is None or not row.get("sha256"):
        raise RuntimeError("download_manifest.tsv lacks the expected SHA256 for GSE183904_RAW.tar")
    return normalize_sha256(row["sha256"])


def build_processed_manifest(root: Path, members: Sequence[tarfile.TarInfo]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for member in members:
        name = Path(member.name).name
        match_parts = name.removesuffix(".csv.gz").split("_", 1)
        accession = match_parts[0]
        sample_id = match_parts[1] if len(match_parts) == 2 else ""
        path = root / "data/processed_input/GSE183904" / name
        rows.append(
            {
                "source_archive": "data/public_downloads/GSE183904_RAW.tar",
                "archive_member_name": name,
                "extracted_path": rel(path, root),
                "geo_accession": accession,
                "sample_id": sample_id,
                "file_size": path.stat().st_size,
                "sha256": sha256_file(path),
                "extraction_date": now_iso(),
                "file_role": "author_filtered_raw_gene_count_csv_gz_for_F1",
                "note": "Compressed csv.gz retained; no permanent plain CSV extraction.",
            }
        )
    return rows


def execute(root: Path) -> int:
    require_paths(root, REQUIRED_INPUTS, STAGE_NAME)
    write_project_structure(root)
    initialize_log(root)
    append_log(
        root,
        f"F0 step1 started; run_id={current_run_id()}; python={sys.version.split()[0]}; os={platform.platform()}",
    )

    tar_path = root / "data/public_downloads/GSE183904_RAW.tar"
    observed_sha = sha256_file(tar_path)
    expected_sha = archive_expected_sha(root)
    if not sha256_equal(observed_sha, expected_sha):
        append_log(root, f"BLOCKING archive SHA256 mismatch: observed={observed_sha}; expected={expected_sha}")
        raise RuntimeError(
            f"GSE183904_RAW.tar SHA256 mismatch: observed={observed_sha}; expected={expected_sha}"
        )

    members = tar_csv_members(tar_path)
    if len(members) != 40:
        append_log(root, f"BLOCKING expected 40 csv.gz members but observed {len(members)}")
        raise RuntimeError(f"Expected 40 csv.gz archive members; observed {len(members)}")
    member_names = [Path(member.name).name for member in members]
    if len(set(member_names)) != 40:
        append_log(root, "BLOCKING archive csv.gz member basenames are not unique")
        raise RuntimeError("GSE183904 archive contains duplicate csv.gz member basenames")

    extract_members(tar_path, members, root / "data/processed_input/GSE183904")
    manifest = build_processed_manifest(root, members)
    fields = [
        "source_archive",
        "archive_member_name",
        "extracted_path",
        "geo_accession",
        "sample_id",
        "file_size",
        "sha256",
        "extraction_date",
        "file_role",
        "note",
    ]
    write_tsv(root / "data/metadata/processed_input_manifest.tsv", manifest, fields)
    append_log(
        root,
        f"F0 step1 completed; archive_sha256={observed_sha}; extracted_members={len(manifest)}",
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, REQUIRED_INPUTS, STAGE_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
