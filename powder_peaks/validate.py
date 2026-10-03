"""输入校验：所有非法输入抛出 ValueError，不写任何报告。"""

import math


def _finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} 必须是数值")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} 必须是有限数值")
    return value


def validate_input(data):
    """校验并规范化输入字典，返回 (a, wavelength, atoms, q_max)。"""
    if not isinstance(data, dict):
        raise ValueError("输入必须是 JSON 对象")

    a = _finite_number(data.get("a"), "晶格边长 a")
    if a <= 0:
        raise ValueError("晶格边长 a 必须为正")

    wavelength = _finite_number(data.get("wavelength"), "波长 wavelength")
    if wavelength <= 0:
        raise ValueError("波长 wavelength 必须为正")

    q_max = data.get("Q")
    if isinstance(q_max, bool) or not isinstance(q_max, int):
        raise ValueError("Q 必须是正整数")
    if not 1 <= q_max <= 400:
        raise ValueError("Q 必须满足 1 <= Q <= 400")

    atoms = data.get("atoms")
    if not isinstance(atoms, list) or not atoms:
        raise ValueError("atoms 必须是非空列表")

    norm_atoms = []
    for i, atom in enumerate(atoms):
        if not isinstance(atom, dict):
            raise ValueError(f"原子 {i} 必须是对象")
        coords = atom.get("position")
        if not isinstance(coords, list) or len(coords) != 3:
            raise ValueError(f"原子 {i} 的 position 必须是三个分数坐标")
        coords = tuple(_finite_number(c, f"原子 {i} 坐标") for c in coords)
        for c in coords:
            if not 0.0 <= c < 1.0:
                raise ValueError(f"原子 {i} 分数坐标必须在 [0,1) 内")
        amplitude = _finite_number(atom.get("amplitude"), f"原子 {i} 散射振幅")
        if amplitude <= 0:
            raise ValueError(f"原子 {i} 散射振幅必须为正")
        occupancy = _finite_number(atom.get("occupancy"), f"原子 {i} 占位率")
        if not 0.0 < occupancy <= 1.0:
            raise ValueError(f"原子 {i} 占位率必须在 (0,1] 内")
        norm_atoms.append(
            {"position": coords, "amplitude": amplitude, "occupancy": occupancy}
        )

    return a, wavelength, norm_atoms, q_max
