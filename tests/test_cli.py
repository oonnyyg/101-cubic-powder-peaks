import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from powder_peaks.validate import InputError, validate_config

VALID = {
    "a": 4.0,
    "wavelength": 1.5406,
    "Q": 8,
    "atoms": [{"x": 0.0, "y": 0.0, "z": 0.0, "amplitude": 1.0, "occupancy": 1.0}],
}


class ValidateTests(unittest.TestCase):
    def check_invalid(self, mutate):
        data = json.loads(json.dumps(VALID))
        mutate(data)
        with self.assertRaises(InputError):
            validate_config(data)

    def test_valid(self):
        config = validate_config(json.loads(json.dumps(VALID)))
        self.assertEqual(config["Q"], 8)
        self.assertEqual(len(config["atoms"]), 1)

    def test_bad_lengths(self):
        self.check_invalid(lambda d: d.update(a=0))
        self.check_invalid(lambda d: d.update(a=-1))
        self.check_invalid(lambda d: d.update(wavelength=float("nan")))
        self.check_invalid(lambda d: d.update(wavelength="1.5"))

    def test_bad_q(self):
        self.check_invalid(lambda d: d.update(Q=0))
        self.check_invalid(lambda d: d.update(Q=401))
        self.check_invalid(lambda d: d.update(Q=2.5))
        self.check_invalid(lambda d: d.update(Q=True))

    def test_bad_atoms(self):
        self.check_invalid(lambda d: d.update(atoms=[]))
        self.check_invalid(lambda d: d["atoms"][0].update(x=1.0))
        self.check_invalid(lambda d: d["atoms"][0].update(y=-0.1))
        self.check_invalid(lambda d: d["atoms"][0].update(amplitude=0))
        self.check_invalid(lambda d: d["atoms"][0].update(amplitude=-2))
        self.check_invalid(lambda d: d["atoms"][0].update(occupancy=0))
        self.check_invalid(lambda d: d["atoms"][0].update(occupancy=1.5))
        self.check_invalid(lambda d: d["atoms"][0].update(occupancy=float("inf")))
        self.check_invalid(lambda d: d["atoms"][0].pop("z"))

    def test_missing_fields(self):
        self.check_invalid(lambda d: d.pop("a"))
        self.check_invalid(lambda d: d.pop("atoms"))


class CliTests(unittest.TestCase):
    def run_cli(self, args):
        return subprocess.run(
            [sys.executable, "-m", "powder_peaks", *args],
            capture_output=True,
            text=True,
            cwd=Path.cwd(),
        )

    def test_success_writes_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            input_path = tmp_path / "in.json"
            input_path.write_text(json.dumps(VALID), encoding="utf-8")
            out_dir = tmp_path / "out"
            proc = self.run_cli([str(input_path), "--out", str(out_dir)])
            self.assertEqual(proc.returncode, 0, proc.stderr)
            csv_path = out_dir / "reflections.csv"
            json_path = out_dir / "peaks.json"
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())
            lines = csv_path.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(lines[0], "h,k,l,q,d,two_theta_deg,intensity")
            self.assertGreater(len(lines), 1)
            peaks = json.loads(json_path.read_text(encoding="utf-8"))["peaks"]
            self.assertTrue(peaks)
            self.assertEqual(peaks[0]["q"], 1)
            self.assertAlmostEqual(max(p["relative_intensity"] for p in peaks), 100.0)

    def test_invalid_input_fails_without_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            bad = json.loads(json.dumps(VALID))
            bad["a"] = -1
            input_path = tmp_path / "bad.json"
            input_path.write_text(json.dumps(bad), encoding="utf-8")
            out_dir = tmp_path / "out"
            proc = self.run_cli([str(input_path), "--out", str(out_dir)])
            self.assertEqual(proc.returncode, 1)
            self.assertIn("错误", proc.stderr)
            self.assertFalse((out_dir / "reflections.csv").exists())
            self.assertFalse((out_dir / "peaks.json").exists())

    def test_missing_file_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "out"
            proc = self.run_cli([str(Path(tmp) / "nope.json"), "--out", str(out_dir)])
            self.assertEqual(proc.returncode, 1)
            self.assertFalse(out_dir.exists())


if __name__ == "__main__":
    unittest.main()
