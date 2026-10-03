import math
import unittest

from powder_peaks.core import (
    compute_powder_pattern,
    enumerate_reflections,
    structure_factor,
)


def simple_atom(**overrides):
    atom = {"x": 0.0, "y": 0.0, "z": 0.0, "amplitude": 1.0, "occupancy": 1.0}
    atom.update(overrides)
    return atom


def config(atoms, q=20, a=4.0, wavelength=1.5406):
    return {"a": a, "wavelength": wavelength, "Q": q, "atoms": atoms}


class EnumerateTests(unittest.TestCase):
    def test_counts_and_signs(self):
        refs = enumerate_reflections(3)
        triples = {(h, k, l) for h, k, l, _ in refs}
        self.assertEqual(len(refs), 6 + 12 + 8)  # q=1,2,3
        self.assertIn((-1, 0, 0), triples)
        self.assertIn((1, -1, 0), triples)
        self.assertIn((1, 1, -1), triples)

    def test_no_reduction_of_multiples(self):
        # (100) 与 (200) 均保留, 不约分
        triples = {(h, k, l) for h, k, l, _ in enumerate_reflections(4)}
        self.assertIn((1, 0, 0), triples)
        self.assertIn((2, 0, 0), triples)

    def test_sorted_by_q(self):
        refs = enumerate_reflections(10)
        qs = [q for *_, q in refs]
        self.assertEqual(qs, sorted(qs))


class StructureFactorTests(unittest.TestCase):
    def test_single_atom_at_origin(self):
        f = structure_factor(2, 3, 4, [simple_atom(amplitude=2.5, occupancy=0.5)])
        self.assertAlmostEqual(f.real, 1.25)
        self.assertAlmostEqual(f.imag, 0.0)

    def test_bcc_extinction(self):
        # 体心: h+k+l 为奇数时消光
        atoms = [simple_atom(), simple_atom(x=0.5, y=0.5, z=0.5)]
        self.assertAlmostEqual(abs(structure_factor(1, 0, 0, atoms)), 0.0, places=12)
        self.assertAlmostEqual(abs(structure_factor(2, 0, 0, atoms)), 2.0)

    def test_fcc_extinction(self):
        # 面心: h,k,l 奇偶混杂时消光
        atoms = [
            simple_atom(),
            simple_atom(x=0.0, y=0.5, z=0.5),
            simple_atom(x=0.5, y=0.0, z=0.5),
            simple_atom(x=0.5, y=0.5, z=0.0),
        ]
        self.assertAlmostEqual(abs(structure_factor(1, 0, 0, atoms)), 0.0, places=12)
        self.assertAlmostEqual(abs(structure_factor(1, 1, 1, atoms)), 4.0)


class PatternTests(unittest.TestCase):
    def test_bragg_filter_and_angles(self):
        result = compute_powder_pattern(config([simple_atom()], q=4))
        for refl in result["reflections"]:
            d = refl["d"]
            self.assertLessEqual(1.5406 / (2.0 * d), 1.0)
            expected = math.degrees(2.0 * math.asin(1.5406 / (2.0 * d)))
            self.assertAlmostEqual(refl["two_theta"], expected)

    def test_peak_merging_and_normalization(self):
        result = compute_powder_pattern(config([simple_atom()], q=12))
        peaks = result["peaks"]
        self.assertTrue(peaks)
        self.assertAlmostEqual(max(p["relative_intensity"] for p in peaks), 100.0)
        by_q = {p["q"]: p for p in peaks}
        self.assertEqual(by_q[1]["reflection_count"], 6)
        self.assertEqual(by_q[2]["reflection_count"], 12)
        self.assertEqual(by_q[3]["reflection_count"], 8)
        # 单原子原点: 每个反射强度为 1, 峰强度 = 反射数
        self.assertAlmostEqual(by_q[1]["intensity"], 6.0)
        self.assertAlmostEqual(by_q[2]["intensity"], 12.0)

    def test_zero_intensity_reflections_kept_but_peaks_dropped(self):
        atoms = [simple_atom(), simple_atom(x=0.5, y=0.5, z=0.5)]
        result = compute_powder_pattern(config(atoms, q=3))
        zero = [r for r in result["reflections"] if r["intensity"] == 0.0]
        self.assertTrue(zero)  # 消光反射保留在反射表中
        qs_in_peaks = {p["q"] for p in result["peaks"]}
        self.assertNotIn(1, qs_in_peaks)  # bcc: q=1 (h+k+l 奇) 整峰消光
        self.assertIn(2, qs_in_peaks)

    def test_extinct_count_reported(self):
        atoms = [simple_atom(), simple_atom(x=0.5, y=0.5, z=0.5)]
        result = compute_powder_pattern(config(atoms, q=6))
        by_q = {p["q"]: p for p in result["peaks"]}
        # q=6 的 (211) 类: h+k+l 偶数保留, 无消光
        self.assertEqual(by_q[6]["extinct_count"], 0)
        # q=4 的 (200): 6 条反射全部保留
        self.assertEqual(by_q[4]["reflection_count"], 6)

    def test_empty_peak_table_when_all_extinct(self):
        # 波长过长: 无反射满足布拉格条件 -> 空峰表
        result = compute_powder_pattern(
            config([simple_atom()], q=1, a=1.0, wavelength=10.0)
        )
        self.assertEqual(result["reflections"], [])
        self.assertEqual(result["peaks"], [])


if __name__ == "__main__":
    unittest.main()
