"""衍射峰位置与强度计算核心。"""

from __future__ import annotations

import cmath
import math

INTENSITY_EPSILON = 1e-12


def enumerate_reflections(q_max: int):
    """枚举 0 < h^2+k^2+l^2 <= q_max 的全部有符号整数三元组。"""
    limit = math.isqrt(q_max)
    reflections = []
    for h in range(-limit, limit + 1):
        for k in range(-limit, limit + 1):
            for l in range(-limit, limit + 1):
                q = h * h + k * k + l * l
                if 0 < q <= q_max:
                    reflections.append((h, k, l, q))
    reflections.sort(key=lambda item: (item[3], item[0], item[1], item[2]))
    return reflections


def structure_factor(h: int, k: int, l: int, atoms) -> complex:
    """F(hkl) = sum_j occ_j * f_j * exp(2*pi*i*(h*x_j + k*y_j + l*z_j))"""
    total = 0j
    for atom in atoms:
        phase = 2.0 * math.pi * (h * atom["x"] + k * atom["y"] + l * atom["z"])
        total += atom["occupancy"] * atom["amplitude"] * cmath.exp(1j * phase)
    return total


def compute_powder_pattern(config: dict) -> dict:
    a = config["a"]
    wavelength = config["wavelength"]
    atoms = config["atoms"]

    # 消光阈值: 全部占位率乘振幅之和的平方的 1e-12 倍
    max_amplitude = sum(atom["occupancy"] * atom["amplitude"] for atom in atoms)
    threshold = INTENSITY_EPSILON * max_amplitude * max_amplitude

    reflections = []
    for h, k, l, q in enumerate_reflections(config["Q"]):
        d = a / math.sqrt(q)
        sin_theta = wavelength / (2.0 * d)
        if sin_theta > 1.0:
            continue
        two_theta = math.degrees(2.0 * math.asin(sin_theta))
        intensity = abs(structure_factor(h, k, l, atoms)) ** 2
        if intensity <= threshold:
            intensity = 0.0
        reflections.append(
            {
                "h": h,
                "k": k,
                "l": l,
                "q": q,
                "d": d,
                "two_theta": two_theta,
                "intensity": intensity,
            }
        )

    # 按 q 合并峰: 累加各反射强度
    peaks_by_q: dict[int, dict] = {}
    for refl in reflections:
        q = refl["q"]
        peak = peaks_by_q.get(q)
        if peak is None:
            peak = {
                "q": q,
                "d": refl["d"],
                "two_theta": refl["two_theta"],
                "intensity": 0.0,
                "reflection_count": 0,
                "extinct_count": 0,
            }
            peaks_by_q[q] = peak
        peak["intensity"] += refl["intensity"]
        peak["reflection_count"] += 1
        if refl["intensity"] == 0.0:
            peak["extinct_count"] += 1

    peaks = [p for p in (peaks_by_q[q] for q in sorted(peaks_by_q)) if p["intensity"] > 0.0]
    if peaks:
        strongest = max(p["intensity"] for p in peaks)
        for peak in peaks:
            peak["relative_intensity"] = 100.0 * peak["intensity"] / strongest

    return {"reflections": reflections, "peaks": peaks}
