from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

from run_analysis import (
    BIN_LABELS,
    callability_table,
    observed_profile_statistics,
    summarize_bootstrap,
)


class StatisticalReportingTests(unittest.TestCase):
    def test_bootstrap_median_is_not_substituted_for_observed_estimate(self) -> None:
        draws = pd.DataFrame(
            {"group": ["A", "A", "A", "A"], "value": [0.8, 0.9, 1.0, 1.1]}
        )
        summary = summarize_bootstrap(
            draws, value="value", group="group", observed={"A": 0.25}
        )
        self.assertEqual(float(summary.loc[0, "estimate"]), 0.25)
        self.assertAlmostEqual(float(summary.loc[0, "bootstrap_median"]), 0.95)

    def test_original_site_profile_statistics_are_reported(self) -> None:
        rows = []
        for index, period in enumerate(BIN_LABELS):
            p_a = index / (len(BIN_LABELS) - 1)
            rows.append(
                {"analysis_bin": period, "A": p_a, "B": 1.0 - p_a}
            )
        profile = pd.DataFrame(rows)
        diversity, turnover = observed_profile_statistics(profile, ["A", "B"])
        self.assertAlmostEqual(diversity[BIN_LABELS[0]], 1.0)
        first_transition = f"{BIN_LABELS[0]} -> {BIN_LABELS[1]}"
        self.assertAlmostEqual(
            turnover[first_transition], 1.0 / (len(BIN_LABELS) - 1)
        )

    def test_sparse_callability_simulation_is_reproducible(self) -> None:
        rows = []
        for period_index, period in enumerate(BIN_LABELS):
            for person in range(6):
                rows.append(
                    {
                        "analysis_bin": period,
                        "molecular_sex": "M" if person < 3 else "F",
                        "mt_call": "H" if person < (period_index % 6) else "",
                        "y_call": "R" if person < (period_index % 3) else "",
                    }
                )
        data = pd.DataFrame(rows)
        _, first = callability_table(data, monte_carlo_resamples=199, seed=17)
        _, second = callability_table(data, monte_carlo_resamples=199, seed=17)
        self.assertTrue(first["sparse_expected_cells"].all())
        self.assertEqual(set(first["p_value_method"]), {"fixed_margin_monte_carlo"})
        np.testing.assert_array_equal(
            first["monte_carlo_p"].to_numpy(), second["monte_carlo_p"].to_numpy()
        )


if __name__ == "__main__":
    unittest.main()
