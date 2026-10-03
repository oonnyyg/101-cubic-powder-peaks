"""反射 CSV 与峰表 JSON 的输出。"""

from __future__ import annotations

import csv
import json
from pathlib import Path


def write_reports(result: dict, out_dir: str | Path) -> tuple[Path, Path]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    csv_path = out_path / "reflections.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["h", "k", "l", "q", "d", "two_theta_deg", "intensity"])
        for refl in result["reflections"]:
            writer.writerow(
                [
                    refl["h"],
                    refl["k"],
                    refl["l"],
                    refl["q"],
                    f"{refl['d']:.10g}",
                    f"{refl['two_theta']:.10g}",
                    f"{refl['intensity']:.10g}",
                ]
            )
    json_path = out_path / "peaks.json"
    peaks_payload = [
        {
            "q": peak["q"],
            "d": peak["d"],
            "two_theta_deg": peak["two_theta"],
            "intensity": peak["intensity"],
            "relative_intensity": peak["relative_intensity"],
            "reflection_count": peak["reflection_count"],
            "extinct_count": peak["extinct_count"],
        }
        for peak in result["peaks"]
    ]
    json_path.write_text(
        json.dumps({"peaks": peaks_payload}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return csv_path, json_path
