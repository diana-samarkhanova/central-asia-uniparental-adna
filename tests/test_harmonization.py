from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import pandas as pd
import numpy as np

from run_analysis import (
    major_haplogroup,
    mask_coordinates,
    named_rng,
    valid_call,
    y_isogg_family_prefix,
)


class HaplogroupHarmonizationTests(unittest.TestCase):
    def test_mtdna_hv_is_not_collapsed_into_h(self) -> None:
        self.assertEqual(major_haplogroup("HV6", "mt"), "HV")
        self.assertEqual(major_haplogroup("haplogroup HV1a", "mt"), "HV")
        self.assertEqual(major_haplogroup("H1a", "mt"), "H")

    def test_basal_y_compounds_are_not_assigned_by_first_letter(self) -> None:
        for call in (
            "BT",
            "CT",
            "CF",
            "F",
            "GHIJK",
            "HIJK",
            "IJK",
            "IJ",
            "K",
            "K2",
            "K2B",
        ):
            with self.subTest(call=call):
                self.assertEqual(
                    major_haplogroup(call, "Y"), "Basal/unresolved"
                )
        self.assertEqual(major_haplogroup("R1a-Z93", "Y"), "R")
        self.assertEqual(major_haplogroup("K2a", "Y"), "K")

    def test_missing_calls(self) -> None:
        for value in (
            "",
            ".",
            "..",
            "?",
            "n/a",
            "N/A (female)",
            "NA",
            "unknown",
            None,
            pd.NA,
            pd.NaT,
            np.nan,
            np.array([np.nan]),
            pd.Series([pd.NA]),
        ):
            with self.subTest(value=value):
                self.assertFalse(valid_call(value))
                self.assertEqual(major_haplogroup(value, "mt"), "")

    def test_y_isogg_prefix_family_encoding(self) -> None:
        self.assertEqual(y_isogg_family_prefix("R1a1a1"), "R1a")
        self.assertEqual(y_isogg_family_prefix("J2a1a4b"), "J2a")
        self.assertEqual(y_isogg_family_prefix("Q1b2b1b2~"), "Q1b")
        self.assertEqual(y_isogg_family_prefix("CF"), "Basal/unresolved")
        self.assertEqual(y_isogg_family_prefix(pd.NA), "")

    def test_coordinate_masking_is_non_mutating_and_one_decimal(self) -> None:
        source = pd.DataFrame(
            {
                "latitude": [43.1234, "51.26", ""],
                "longitude": [76.9876, "71.44", ""],
                "individual_id": ["A", "B", "C"],
            }
        )
        masked = mask_coordinates(source, digits=1)
        self.assertEqual(masked["latitude"].iloc[:2].tolist(), [43.1, 51.3])
        self.assertEqual(masked["longitude"].iloc[:2].tolist(), [77.0, 71.4])
        self.assertEqual(source.loc[0, "latitude"], 43.1234)
        self.assertEqual(masked["individual_id"].tolist(), ["A", "B", "C"])

    def test_named_rng_streams_are_repeatable_and_independent(self) -> None:
        first = named_rng(20260725, "mtDNA:bootstrap").integers(0, 1000, 8)
        repeated = named_rng(20260725, "mtDNA:bootstrap").integers(0, 1000, 8)
        other = named_rng(20260725, "Y:bootstrap").integers(0, 1000, 8)
        self.assertEqual(first.tolist(), repeated.tolist())
        self.assertNotEqual(first.tolist(), other.tolist())


if __name__ == "__main__":
    unittest.main()
