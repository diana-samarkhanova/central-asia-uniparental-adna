"""Guards against false reproduction claims and accidental release overwrite."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "analysis/reproduce.py"
SPEC = importlib.util.spec_from_file_location("reproduce", SCRIPT)
reproduce = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reproduce)


class ReproductionRunnerTests(unittest.TestCase):
    def test_changed_input_is_rejected_before_statistical_run(self):
        with tempfile.TemporaryDirectory() as folder:
            changed = Path(folder) / "input.csv"
            changed.write_text("latitude_0_1deg\n0.0\n")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                reproduce.verify_analytical_input(changed)
        self.assertEqual(reproduce.verify_analytical_input(reproduce.PUBLIC_INPUT),
                         reproduce.FROZEN_INPUT_SHA256)

    def test_cannot_write_inside_release_or_existing_results(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / "repo"
            root.mkdir()
            with self.assertRaises(ValueError):
                reproduce.reserve_output(root / "results", root)
            old = Path(folder) / "old"
            old.mkdir()
            (old / "important.txt").write_text("preserve")
            with self.assertRaises(ValueError):
                reproduce.reserve_output(old, root)
            self.assertEqual((old / "important.txt").read_text(), "preserve")

    def test_compare_requires_files_and_matching_schema(self):
        with tempfile.TemporaryDirectory() as folder:
            a, b = Path(folder) / "a.csv", Path(folder) / "b.csv"
            a.write_text("marker,estimate\nmtDNA,0.5\n")
            self.assertEqual(reproduce.compare_csv(a, b)["status"], "missing")
            b.write_text("marker,value\nmtDNA,0.5\n")
            self.assertEqual(reproduce.compare_csv(a, b)["status"], "schema_mismatch")

    def test_roundoff_allowed_but_scientific_or_label_changes_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            a, b = Path(folder) / "a.csv", Path(folder) / "b.csv"
            a.write_text("marker,estimate\nmtDNA,0.5\n")
            b.write_text("marker,estimate\nmtDNA,0.500000000000001\n")
            self.assertEqual(reproduce.compare_csv(a, b)["status"], "match")
            b.write_text("marker,estimate\nmtDNA,0.6\n")
            self.assertEqual(reproduce.compare_csv(a, b)["status"], "different")
            b.write_text("marker,estimate\nY,0.5\n")
            self.assertEqual(reproduce.compare_csv(a, b)["different_columns"], ["marker"])

    def test_empty_and_nonfinite_values_do_not_mask_difference(self):
        with tempfile.TemporaryDirectory() as folder:
            a, b = Path(folder) / "a.csv", Path(folder) / "b.csv"
            a.write_text("marker,estimate\nmtDNA,\nY,inf\n")
            b.write_text("marker,estimate\nmtDNA,\nY,inf\n")
            self.assertEqual(reproduce.compare_csv(a, b)["status"], "match")
            b.write_text("marker,estimate\nmtDNA,0\nY,inf\n")
            self.assertEqual(reproduce.compare_csv(a, b)["status"], "different")


if __name__ == "__main__":
    unittest.main()
