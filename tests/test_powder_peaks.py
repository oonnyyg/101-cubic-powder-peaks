import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from powder_peaks.core import compute_powder_pattern
from powder_peaks.validate import validate_input


def atom(x, y, z, f=1.0, occ=1.0):
    return {"position": [x, y, z], "amplitude": f, "occupancy": occ}


class TestCore(unittest.TestCase):
    def test_simple_cubic_all_present(self):
        refs, peaks = compute_powder_pattern(4.0, 1.0, [atom(0, 0, 0)], 10)
        qs = {r["q"] for r in refs}
        self.assertEqual(qs, {1, 2, 3, 4, 5, 6, 8, 9, 10})
        for r in refs:
            self.assertAlmostEqual(r["intensity"], 1.0)
        self.assertEqual(len(peaks), len(qs))
        p1 = next(p for p in peaks if p["q"] == 1)
        p2 = next(p for p in peaks if p["q"] == 2)
        self.assertAlmostEqual(p2["relative_intensity"],
                               2.0 * p1["relative_intensity"])
        self.assertAlmostEqual(max(p["relative_intensity"] for p in peaks), 100.0)

    def test_no_index_reduction(self):
        refs, _ = compute_powder_pattern(4.0, 0.5, [atom(0, 0, 0)], 10)
        hkl = {(r["h"], r["k"], r["l"]) for r in refs}
        self.assertIn((1, 0, 0), hkl)
        self.assertIn((2, 0, 0), hkl)
        self.assertIn((-2, 0, 0), hkl)

    def test_bcc_systematic_absence(self):
        atoms = [atom(0, 0, 0), atom(0.5, 0.5, 0.5)]
        refs, peaks = compute_powder_pattern(3.0, 1.0, atoms, 20)
        for r in refs:
            if (r["h"] + r["k"] + r["l"]) % 2 == 1:
                self.assertEqual(r["intensity"], 0.0)
            else:
                self.assertAlmostEqual(r["intensity"], 4.0)
        self.assertEqual({p["q"] for p in peaks},
                         {2, 4, 6, 8, 10, 12, 14, 16, 18, 20})

    def test_fcc_mixed_parity_extinct(self):
        atoms = [atom(0, 0, 0), atom(0, 0.5, 0.5),
                 atom(0.5, 0, 0.5), atom(0.5, 0.5, 0)]
        refs, _ = compute_powder_pattern(5.0, 1.0, atoms, 20)
        for r in refs:
            parities = {r["h"] % 2, r["k"] % 2, r["l"] % 2}
            if len(parities) == 2:
                self.assertEqual(r["intensity"], 0.0)
            else:
                self.assertAlmostEqual(r["intensity"], 16.0)

    def test_zero_intensity_reflections_kept(self):
        atoms = [atom(0, 0, 0), atom(0.5, 0.5, 0.5)]
        refs, peaks = compute_powder_pattern(3.0, 1.0, atoms, 5)
        self.assertTrue(any(r["intensity"] == 0.0 for r in refs))
        self.assertTrue(all(p["intensity"] > 0.0 for p in peaks))

    def test_peak_merges_multiplicity(self):
        refs, peaks = compute_powder_pattern(4.0, 0.5, [atom(0, 0, 0)], 2)
        p1 = next(p for p in peaks if p["q"] == 1)
        self.assertEqual(p1["n_reflections"], 6)
        self.assertAlmostEqual(p1["intensity"], 6.0)

    def test_bragg_geometry(self):
        refs, _ = compute_powder_pattern(4.0, 1.0, [atom(0, 0, 0)], 3)
        r1 = next(r for r in refs if (r["h"], r["k"], r["l"]) == (1, 0, 0))
        self.assertAlmostEqual(r1["d"], 4.0)
        self.assertAlmostEqual(r1["two_theta"],
                               math.degrees(2 * math.asin(1.0 / 8.0)))

    def test_resolution_cutoff_empty(self):
        refs, peaks = compute_powder_pattern(1.0, 10.0, [atom(0, 0, 0)], 10)
        self.assertEqual(refs, [])
        self.assertEqual(peaks, [])

    def test_relative_normalization(self):
        atoms = [atom(0, 0, 0, f=2.0), atom(0.5, 0.5, 0.5, f=1.0)]
        _, peaks = compute_powder_pattern(4.0, 0.5, atoms, 10)
        self.assertAlmostEqual(max(p["relative_intensity"] for p in peaks), 100.0)


class TestValidate(unittest.TestCase):
    def base(self):
        return {"a": 4.0, "wavelength": 1.54, "Q": 10,
                "atoms": [atom(0, 0, 0)]}

    def test_ok(self):
        a, wl, atoms, q = validate_input(self.base())
        self.assertEqual((a, wl, q), (4.0, 1.54, 10))

    def check_bad(self, **kw):
        data = self.base()
        data.update(kw)
        with self.assertRaises(ValueError):
            validate_input(data)

    def test_bad_values(self):
        self.check_bad(a=0)
        self.check_bad(a=-1)
        self.check_bad(a=float("nan"))
        self.check_bad(wavelength=-1.0)
        self.check_bad(wavelength=float("inf"))
        self.check_bad(Q=0)
        self.check_bad(Q=401)
        self.check_bad(Q=2.5)
        self.check_bad(atoms=[])
        self.check_bad(atoms=[atom(1.0, 0, 0)])
        self.check_bad(atoms=[atom(-0.1, 0, 0)])
        self.check_bad(atoms=[atom(0, 0, 0, f=0)])
        self.check_bad(atoms=[atom(0, 0, 0, occ=0)])
        self.check_bad(atoms=[atom(0, 0, 0, occ=1.5)])


class TestCli(unittest.TestCase):
    def run_cli(self, tmp, payload):
        inp = Path(tmp) / "in.json"
        inp.write_text(json.dumps(payload), encoding="utf-8")
        out = Path(tmp) / "out"
        proc = subprocess.run(
            [sys.executable, "-m", "powder_peaks", str(inp), "--out", str(out)],
            capture_output=True, text=True,
            cwd=Path(__file__).resolve().parent.parent)
        return proc, out

    def test_cli_success(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            payload = {"a": 4.0, "wavelength": 1.0, "Q": 5,
                       "atoms": [atom(0, 0, 0), atom(0.5, 0.5, 0.5)]}
            proc, out = self.run_cli(tmp, payload)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            table = json.loads((out / "peaks.json").read_text(encoding="utf-8"))
            self.assertEqual([p["q"] for p in table["peaks"]], [2, 4])
            lines = (out / "reflections.csv").read_text(
                encoding="utf-8").strip().splitlines()
            self.assertEqual(lines[0], "h,k,l,q,d,two_theta_deg,intensity")
            self.assertGreater(len(lines), 1)

    def test_cli_invalid_writes_nothing(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            payload = {"a": -1.0, "wavelength": 1.0, "Q": 5,
                       "atoms": [atom(0, 0, 0)]}
            proc, out = self.run_cli(tmp, payload)
            self.assertEqual(proc.returncode, 1)
            self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
