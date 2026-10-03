"""命令行入口：python -m powder_peaks 输入.json --out 输出目录"""

import argparse
import csv
import json
import sys
from pathlib import Path

from .core import compute_powder_pattern
from .validate import validate_input


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="powder_peaks",
        description="理想立方晶胞粉末衍射峰位置与相对强度计算",
    )
    parser.add_argument("input", help="输入 JSON 文件路径")
    parser.add_argument("--out", required=True, help="输出目录")
    args = parser.parse_args(argv)

    try:
        with open(args.input, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        a, wavelength, atoms, q_max = validate_input(data)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    reflections, peaks = compute_powder_pattern(a, wavelength, atoms, q_max)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "reflections.csv"
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["h", "k", "l", "q", "d", "two_theta_deg", "intensity"])
        for r in reflections:
            writer.writerow([r["h"], r["k"], r["l"], r["q"],
                             f"{r['d']:.10g}", f"{r['two_theta']:.10g}",
                             f"{r['intensity']:.10g}"])

    json_path = out_dir / "peaks.json"
    peak_table = [
        {
            "q": p["q"],
            "d": p["d"],
            "two_theta_deg": p["two_theta"],
            "intensity": p["intensity"],
            "relative_intensity": p["relative_intensity"],
            "n_reflections": p["n_reflections"],
        }
        for p in peaks
    ]
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump({"peaks": peak_table}, fh, ensure_ascii=False, indent=2)

    print(f"已写出 {csv_path} 与 {json_path}："
          f"{len(reflections)} 条反射，{len(peaks)} 个非零峰")
    return 0


if __name__ == "__main__":
    sys.exit(main())
