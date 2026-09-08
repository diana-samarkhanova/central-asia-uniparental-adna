"""Regression checks for independent raw-source reproduction."""

from __future__ import annotations

import copy
import os
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))

from run_analysis import (  # noqa: E402
    BIN_LABELS,
    EXPECTED_INPUT_SHA256,
    bootstrap_site_profiles,
    frozen_catalogue_order,
    named_rng,
)
from test_full_outputs import check_analysis_provenance, sha256  # noqa: E402
from extract_aadr_input import (  # noqa: E402
    EXPECTED_OUTPUT_SHA256,
    extract_aadr_input,
)


class RawReproductionTests(unittest.TestCase):
    def test_extractor_preserves_sources_existing_outputs_and_work(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.anno"
            source.write_text("preserve source")
            with self.assertRaisesRegex(ValueError, "new file"):
                extract_aadr_input(source, root / "new-work", source)
            output = root / "output.csv"
            output.write_text("preserve output")
            with self.assertRaisesRegex(ValueError, "new file"):
                extract_aadr_input(source, root / "new-work", output)
            work = root / "old-work"
            work.mkdir()
            (work / "important.txt").write_text("preserve work")
            with self.assertRaisesRegex(ValueError, "absent or empty"):
                extract_aadr_input(source, work, root / "new-output.csv")
            self.assertEqual(source.read_text(), "preserve source")
            self.assertEqual(output.read_text(), "preserve output")
            self.assertEqual((work / "important.txt").read_text(), "preserve work")

    def test_extractor_rejects_changed_source_before_creating_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "incorrect_aadr.anno"
            source.write_text("changed source\n", encoding="utf-8")
            work = root / "work"
            output = root / "result.csv"
            with self.assertRaisesRegex(ValueError, "AADR SHA-256 mismatch"):
                extract_aadr_input(source, work, output)
            self.assertFalse(work.exists())
            self.assertFalse(output.exists())

    @unittest.skipUnless(
        os.environ.get("ADNA_AADR_FIXTURE"),
        "Set ADNA_AADR_FIXTURE to the verified frozen upstream AADR file",
    )
    def test_verified_aadr_rebuilds_exact_released_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "input.csv"
            report = extract_aadr_input(
                Path(os.environ["ADNA_AADR_FIXTURE"]), root / "work", output
            )
            self.assertEqual(sha256(output), EXPECTED_OUTPUT_SHA256)
            self.assertEqual(
                tuple(report[key] for key in ("primary_records", "primary_sites", "mt_calls", "y_calls")),
                (489, 136, 438, 229),
            )
            self.assertTrue(report["matches_frozen_input_byte_for_byte"])
            self.assertFalse(report["cross_database_audits_recomputed"])

    def test_country_date_ties_have_deterministic_record_order(self) -> None:
        source = pd.DataFrame(
            {
                "individual_id": ["example_c", "example_b", "example_d", "example_a"],
                "country": ["Uzbekistan", "Kazakhstan", "Kazakhstan", "Kazakhstan"],
                "date_bp": [6000, 4000, 3000, 4000],
            }
        )
        expected = ["example_a", "example_b", "example_d", "example_c"]
        self.assertEqual(frozen_catalogue_order(source).individual_id.tolist(), expected)
        self.assertEqual(
            frozen_catalogue_order(source.iloc[::-1]).individual_id.tolist(), expected
        )
        self.assertEqual(source.individual_id.iloc[0], "example_c")

    def test_raw_order_reproduces_saved_catalogue_bootstrap_draws(self) -> None:
        # Every cluster contributes to all periods, with unequal call counts.
        # Reversing records changes first-occurrence cluster order, hence which
        # cluster receives each multinomial weight at a fixed seed.
        rows = []
        for period_index, period in enumerate(BIN_LABELS):
            for site_index in range(3):
                for call_index in range(site_index + 1):
                    rows.append(
                        {
                            "individual_id": f"example_{len(rows):04d}",
                            "country": "Kazakhstan",
                            "locality": f"synthetic_site_{site_index}",
                            "date_bp": 5400 - period_index * 500,
                            "analysis_bin": period,
                            "marker": "A" if (site_index + period_index + call_index) % 3 else "B",
                        }
                    )
        released = pd.DataFrame(rows)
        source_order = released.iloc[::-1].reset_index(drop=True)
        _, expected_tv, _ = bootstrap_site_profiles(
            released, "marker", ["A", "B"], 8, named_rng(741, "order-regression")
        )
        _, actual_tv, _ = bootstrap_site_profiles(
            frozen_catalogue_order(source_order),
            "marker", ["A", "B"], 8, named_rng(741, "order-regression"),
        )
        _, wrong_tv, _ = bootstrap_site_profiles(
            source_order, "marker", ["A", "B"], 8, named_rng(741, "order-regression")
        )
        pd.testing.assert_frame_equal(actual_tv, expected_tv)
        self.assertFalse(wrong_tv.equals(expected_tv))

    def test_fresh_raw_manifest_uses_frozen_inputs_without_recompute_field(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "raw.py"
            recompute = Path(temporary) / "recompute.py"
            script.write_text("raw extraction\n", encoding="utf-8")
            recompute.write_text("catalogue rerun\n", encoding="utf-8")
            manifest = {
                "source_code_sha256": sha256(script),
                "inputs": {f"raw/{name}": digest for name, digest in EXPECTED_INPUT_SHA256.items()},
            }
            check_analysis_provenance(manifest, script, recompute)
            altered = copy.deepcopy(manifest)
            altered["inputs"]["raw/amtdb"] = "0" * 64
            with self.assertRaises(AssertionError):
                check_analysis_provenance(altered, script, recompute)
            rerun = dict(manifest, recomputed_from={"script_sha256": sha256(recompute)})
            check_analysis_provenance(rerun, script, recompute)
            rerun["recomputed_from"]["script_sha256"] = "0" * 64
            with self.assertRaises(AssertionError):
                check_analysis_provenance(rerun, script, recompute)


if __name__ == "__main__":
    unittest.main()
