"""输入数据的读取与校验。"""

from __future__ import annotations

import json
import math
from pathlib import Path


class InputError(ValueError):
    """输入数据非法时抛出。"""


def _finite_number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InputError(f"{name} 必须是数值，得到 {value!r}")
    result = float(value)
    if not math.isfinite(result):
        raise InputError(f"{name} 必须是有限数值，得到 {value!r}")
    return result


def _positive_number(value, name: str) -> float:
    result = _finite_number(value, name)
    if result <= 0.0:
        raise InputError(f"{name} 必须为正数，得到 {value!r}")
    return result


def _fraction_coordinate(value, name: str) -> float:
    result = _finite_number(value, name)
    if not 0.0 <= result < 1.0:
        raise InputError(f"{name} 必须在 [0,1) 内，得到 {value!r}")
    return result


def validate_config(data) -> dict:
    if not isinstance(data, dict):
        raise InputError("输入 JSON 顶层必须是对象")

    try:
        a = _positive_number(data["a"], "晶格边长 a")
        wavelength = _positive_number(data["wavelength"], "波长 wavelength")
        q_max_raw = data["Q"]
        atoms_raw = data["atoms"]
    except KeyError as exc:
        raise InputError(f"缺少必需字段: {exc.args[0]}") from None

    if isinstance(q_max_raw, bool) or not isinstance(q_max_raw, int):
        raise InputError(f"Q 必须是正整数，得到 {q_max_raw!r}")
    if not 1 <= q_max_raw <= 400:
        raise InputError(f"Q 必须在 1..400 内，得到 {q_max_raw!r}")

    if not isinstance(atoms_raw, list) or not atoms_raw:
        raise InputError("atoms 必须是非空列表")

    atoms = []
    for index, atom in enumerate(atoms_raw):
        prefix = f"atoms[{index}]"
        if not isinstance(atom, dict):
            raise InputError(f"{prefix} 必须是对象")
        try:
            x = _fraction_coordinate(atom["x"], f"{prefix}.x")
            y = _fraction_coordinate(atom["y"], f"{prefix}.y")
            z = _fraction_coordinate(atom["z"], f"{prefix}.z")
            amplitude = _positive_number(atom["amplitude"], f"{prefix}.amplitude")
            occupancy = _finite_number(atom["occupancy"], f"{prefix}.occupancy")
        except KeyError as exc:
            raise InputError(f"{prefix} 缺少字段: {exc.args[0]}") from None
        if not 0.0 < occupancy <= 1.0:
            raise InputError(f"{prefix}.occupancy 必须在 (0,1] 内，得到 {occupancy!r}")
        atoms.append(
            {"x": x, "y": y, "z": z, "amplitude": amplitude, "occupancy": occupancy}
        )

    return {"a": a, "wavelength": wavelength, "Q": q_max_raw, "atoms": atoms}


def load_config(path: str | Path) -> dict:
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except OSError as exc:
        raise InputError(f"无法读取输入文件 {path}: {exc}") from None
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise InputError(f"输入文件不是合法 JSON: {exc}") from None
    return validate_config(data)
