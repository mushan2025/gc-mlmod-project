#!/usr/bin/env python3
"""Validate whether Xena star_counts values invert to integer counts."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import platform
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def validate(input_path: Path, tolerance: float) -> dict[str, str | int | float]:
    numeric_values = 0
    missing_values = 0
    nonfinite_values = 0
    noninteger_after_inverse = 0
    max_inverse_integer_error = 0.0
    min_stored_value = math.inf
    max_stored_value = -math.inf
    feature_rows = 0

    with gzip.open(input_path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        expected_columns = len(header)

        for feature_rows, row in enumerate(reader, start=1):
            if len(row) != expected_columns:
                raise ValueError(
                    f"Row {feature_rows} has {len(row)} columns; "
                    f"expected {expected_columns}."
                )

            for value_text in row[1:]:
                if value_text == "" or value_text.upper() == "NA":
                    missing_values += 1
                    continue

                value = float(value_text)
                if not math.isfinite(value):
                    nonfinite_values += 1
                    continue

                reconstructed = math.pow(2.0, value) - 1.0
                error = abs(reconstructed - round(reconstructed))
                numeric_values += 1
                max_inverse_integer_error = max(max_inverse_integer_error, error)
                min_stored_value = min(min_stored_value, value)
                max_stored_value = max(max_stored_value, value)
                if error > tolerance:
                    noninteger_after_inverse += 1

    validation_status = (
        "PASS"
        if nonfinite_values == 0 and noninteger_after_inverse == 0
        else "FAIL"
    )
    return {
        "source_file": input_path.as_posix(),
        "source_sha256": sha256_file(input_path),
        "n_feature_rows": feature_rows,
        "n_total_columns": expected_columns,
        "n_sample_columns": expected_columns - 1,
        "n_numeric_values": numeric_values,
        "n_missing_values": missing_values,
        "n_nonfinite_values": nonfinite_values,
        "transform_tested": "round(2^x - 1)",
        "integer_tolerance": tolerance,
        "n_noninteger_after_inverse": noninteger_after_inverse,
        "max_inverse_integer_error": max_inverse_integer_error,
        "min_stored_value": min_stored_value,
        "max_stored_value": max_stored_value,
        "validation_status": validation_status,
        "allowed_use": "reconstructed_integer_counts_for_DESeq2_VST_after_F7_input_audit",
        "forbidden_use": "automatic_authorization_for_DESeq2_differential_expression",
        "python_version": platform.python_version(),
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "data/public_downloads/TCGA_STAD/TCGA-STAD.star_counts.tsv.gz"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "results/F0_audit/TCGA_STAD_star_counts_inverse_validation.tsv"
        ),
    )
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    result = validate(args.input, args.tolerance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result), delimiter="\t")
        writer.writeheader()
        writer.writerow(result)

    print(
        f"{result['validation_status']}: {result['n_numeric_values']} values; "
        f"max error={result['max_inverse_integer_error']}"
    )


if __name__ == "__main__":
    main()
