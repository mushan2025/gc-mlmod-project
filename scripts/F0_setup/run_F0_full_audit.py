#!/usr/bin/env python3
"""Run the four approved F0 audit stages in order.

Without ``--execute`` this wrapper checks the current preregistered inputs and
planned output contract without writing files. Formal execution requires
``--execute`` and stops immediately if any stage reaches a pause or blocking
condition. It never starts F1.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import F0_step1_structure_and_extract as step1
import F0_step2_sample_info_and_audit as step2
import F0_step3_inventory_and_markers as step3
import F0_step4_decisions_and_gate as step4
from f0_utils import F0_OUTPUTS, REQUIRED_INPUTS, dry_run_report, parse_stage_args


STAGE_NAME = "F0 full staged audit"


def execute(root: Path) -> int:
    for stage in [step1, step2, step3, step4]:
        result = stage.execute(root)
        if result != 0:
            return result
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, REQUIRED_INPUTS, F0_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
