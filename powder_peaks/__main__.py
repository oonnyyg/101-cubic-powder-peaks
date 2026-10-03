"""命令行入口: python -m powder_peaks 输入.json --out 输出目录"""

from __future__ import annotations

import argparse
import sys

from .core import compute_powder_pattern
from .report import write_reports
from .validate import InputError, load_config


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="powder_peaks",
        description="计算理想立方晶胞的粉末衍射峰位置和相对强度",
    )
    parser.add_argument("input", help="输入 JSON 文件路径")
    parser.add_argument("--out", required=True, help="输出目录")
    args = parser.parse_args(argv)

    try:
        config = load_config(args.input)
        result = compute_powder_pattern(config)
    except InputError as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    csv_path, json_path = write_reports(result, args.out)
    print(f"反射表: {csv_path}")
    print(f"峰表: {json_path}")
    print(
        f"共 {len(result['reflections'])} 条反射, {len(result['peaks'])} 个非零强度峰"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
