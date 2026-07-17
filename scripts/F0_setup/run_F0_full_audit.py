#!/usr/bin/env python3
"""按固定顺序运行四个已批准的 F0 审计阶段。

这是 F0 的总入口脚本，负责依次调用 Step1、Step2、Step3 和 Step4。

输入：各阶段在 ``f0_utils.REQUIRED_INPUTS`` 中登记的文件。
输出：各阶段共同生成 ``f0_utils.F0_OUTPUTS`` 中登记的 17 个正式产物。

安全规则：
- 不加 ``--execute`` 时只做 dry run，不写任何正式文件；
- 正式执行必须显式加入 ``--execute``；
- 任一阶段返回非 0 状态时立即停止；
- 本脚本只完成 F0，绝不会自动启动 F1。
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
    """依次执行四个阶段，并把第一个失败状态返回给调用者。"""

    # 顺序不能调换：后一步会读取前一步生成的表格和审计状态。
    for stage in [step1, step2, step3, step4]:
        result = stage.execute(root)
        if result != 0:
            return result
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """解析命令行，并在 dry run 与正式执行之间作明确分流。"""

    args = parse_stage_args(__doc__ or STAGE_NAME, argv)
    root = Path(args.project_root).resolve()
    # 默认路径是只读检查。只有用户明确写出 --execute 才会调用 execute()。
    if not args.execute:
        return dry_run_report(root, STAGE_NAME, REQUIRED_INPUTS, F0_OUTPUTS)
    return execute(root)


if __name__ == "__main__":
    raise SystemExit(main())
