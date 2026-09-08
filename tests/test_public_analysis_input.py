import csv
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "derived" / "central_asia_analysis_input_v1.csv"
EXPECTED_COLUMNS = [
    "project_record_id",
    "country",
    "site_key",
    "analysis_bin",
    "study_key",
    "molecular_sex",
    "mt_category",
    "y_category",
    "y_prefix_category",
    "mt_called",
    "y_called",
    "strict_qc",
    "population_outlier",
    "direct_date",
    "date_bp",
    "date_sd_bp",
    "kin_representative_mt",
    "kin_representative_y",
    "latitude_0_1deg",
    "longitude_0_1deg",
]


class PublicAnalysisInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with INPUT.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            cls.header = reader.fieldnames
            cls.rows = list(reader)

    def test_schema_and_fixed_counts(self):
        self.assertEqual(self.header, EXPECTED_COLUMNS)
        self.assertEqual(len(self.rows), 489)
        self.assertEqual(
            len({(row["country"], row["site_key"]) for row in self.rows}),
            136,
        )
        self.assertEqual(
            sum(row["mt_called"] == "True" for row in self.rows), 438
        )
        self.assertEqual(
            sum(row["y_called"] == "True" for row in self.rows), 229
        )

    def test_only_project_codes_are_exposed(self):
        record_ids = [row["project_record_id"] for row in self.rows]
        self.assertEqual(len(record_ids), len(set(record_ids)))
        self.assertTrue(all(re.fullmatch(r"CAU[0-9]{6}", value) for value in record_ids))
        self.assertTrue(
            all(re.fullmatch(r"SITE[0-9]{4}", row["site_key"]) for row in self.rows)
        )
        self.assertTrue(
            all(re.fullmatch(r"STUDY[0-9]{3}", row["study_key"]) for row in self.rows)
        )
        payload = INPUT.read_text(encoding="utf-8")
        self.assertNotRegex(payload, r"CKZ00[1-4]")
        self.assertNotIn("Bestamak", payload)
        self.assertNotIn("Biestamak", payload)

    def test_marker_and_coordinate_guards(self):
        self.assertTrue(
            all(
                row["molecular_sex"] == "M"
                for row in self.rows
                if row["y_called"] == "True"
            )
        )
        for row in self.rows:
            coordinate_values = (
                row["latitude_0_1deg"], row["longitude_0_1deg"]
            )
            self.assertEqual(
                coordinate_values[0] == "",
                coordinate_values[1] == "",
                "latitude and longitude must be jointly present or absent",
            )
            for column in ("latitude_0_1deg", "longitude_0_1deg"):
                if row[column] == "":
                    continue
                value = float(row[column])
                self.assertAlmostEqual(value, round(value, 1), places=12)


if __name__ == "__main__":
    unittest.main()
