"""核心计算：枚举反射、结构因子、按 q 合并峰。"""

import cmath
import math


def _structure_factor(h, k, l, atoms):
    """F(hkl) = sum_j occ_j * f_j * exp(2*pi*i*(h*x_j + k*y_j + l*z_j))"""
    total = 0j
    for atom in atoms:
        x, y, z = atom["position"]
        phase = 2.0 * math.pi * (h * x + k * y + l * z)
        total += atom["occupancy"] * atom["amplitude"] * cmath.exp(1j * phase)
    return total


def compute_powder_pattern(a, wavelength, atoms, q_max):
    """返回 (reflections, peaks)，峰按 q 合并，最强峰归一为 100。"""
    limit = int(math.isqrt(q_max))
    f_scale = sum(atom["occupancy"] * atom["amplitude"] for atom in atoms)
    threshold = 1e-12 * f_scale * f_scale

    reflections = []
    for h in range(-limit, limit + 1):
        for k in range(-limit, limit + 1):
            for m in range(-limit, limit + 1):
                q = h * h + k * k + m * m
                if q == 0 or q > q_max:
                    continue
                d = a / math.sqrt(q)
                ratio = wavelength / (2.0 * d)
                if ratio > 1.0:
                    continue
                two_theta = math.degrees(2.0 * math.asin(ratio))
                intensity = abs(_structure_factor(h, k, m, atoms)) ** 2
                if intensity <= threshold:
                    intensity = 0.0
                reflections.append({
                    "h": h, "k": k, "l": m, "q": q, "d": d,
                    "two_theta": two_theta, "intensity": intensity,
                })

    reflections.sort(key=lambda r: (r["q"], r["h"], r["k"], r["l"]))

    by_q = {}
    for r in reflections:
        group = by_q.setdefault(r["q"], {
            "q": r["q"], "d": r["d"], "two_theta": r["two_theta"],
            "intensity": 0.0, "n_reflections": 0,
        })
        group["intensity"] += r["intensity"]
        group["n_reflections"] += 1

    peaks = [by_q[q] for q in sorted(by_q) if by_q[q]["intensity"] > 0.0]
    if peaks:
        i_max = max(p["intensity"] for p in peaks)
        for p in peaks:
            p["relative_intensity"] = 100.0 * p["intensity"] / i_max

    return reflections, peaks
