#!/usr/bin/env python3
"""Integrity checks for the released Central Asia synthesis outputs."""

from __future__ import annotations

import hashlib
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


FORBIDDEN_GRANULAR_TABLES = (
    "aadr_central_asia_unique_individual_catalogue.csv",
    "aadr_deduplication_audit.csv",
    "aadr_primary_analysis_catalogue.csv",
    "cross_database_exact_id_audit.csv",
    "extended_legacy_mtdna_not_in_aadr.csv",
    "site_profiles_mtdna.csv",
    "site_profiles_y.csv",
    "cluster_model_residual_diagnostics.csv",
    "dispersion_profiles_mtdna.csv",
    "dispersion_profiles_y.csv",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", nargs="?", type=Path,
                        default=Path("outputs/adna_central_asia_final"))
    parser.add_argument("--input", type=Path,
                        help="Require this input path in the sensitivity provenance")
    args = parser.parse_args()
    out = args.out
    tables = out / "tables"
    summary = json.loads((out / "results_summary.json").read_text(encoding="utf-8"))

    assert summary["aadr_catalogue_unique_archaeological_individuals"] == 501
    assert summary["primary_analysis_unique_individuals_3500BCE_to_1500CE"] == 489
    assert summary["primary_mt_calls"] == 438
    assert summary["primary_y_calls_in_molecular_males"] == 229
    assert summary["primary_sites"] == 136
    for filename in FORBIDDEN_GRANULAR_TABLES:
        assert not (tables / filename).exists(), (
            f"granular table must not be committed: {filename}"
        )

    for marker in ("mtdna", "y"):
        for weighting in ("individual_weighted", "site_balanced"):
            df = pd.read_csv(
                tables / f"composition_{marker}_{weighting}.csv"
            )
            category_cols = [
                c for c in df.columns if c not in {"analysis_bin", "n_calls", "n_sites"}
            ]
            sums = df.loc[df["n_calls"] > 0, category_cols].sum(axis=1)
            assert np.allclose(sums, 1.0, atol=1e-10)

    expected_diversity = {
        "mtDNA": np.array([8.09304820062489, 7.37924895726579, 7.05598715030903,
                           10.6233042205366, 10.214034086894, 7.28135866692849,
                           6.54759028203901, 5.27766344032807]),
        "Y": np.array([3.28516928489741, 6.59103477659705, 2.143236483694,
                       3.46596156321367, 4.27668004338971, 4.77593764148247,
                       2.60647903043411, 4.33939062525723]),
    }
    expected_turnover = {
        "mtDNA": np.array([0.475190434565435, 0.312396284271284,
                           0.428058053058053, 0.296958348556814,
                           0.579528985507246, 0.553125, 0.332142857142857]),
        "Y": np.array([0.442063492063492, 0.594907407407407,
                       0.24537037037037, 0.196623093681917,
                       0.291316526610644, 0.648809523809524,
                       0.638888888888889]),
    }
    diversity_summary = pd.read_csv(tables / "diversity_site_bootstrap_summary.csv")
    turnover_summary = pd.read_csv(tables / "turnover_site_bootstrap_summary.csv")
    for frame, expected in (
        (diversity_summary, expected_diversity),
        (turnover_summary, expected_turnover),
    ):
        assert "bootstrap_median" in frame
        assert not np.allclose(frame["estimate"], frame["bootstrap_median"])
        for marker in ("mtDNA", "Y"):
            observed = frame.loc[frame["marker"] == marker, "estimate"].to_numpy(float)
            assert np.allclose(observed, expected[marker], atol=1e-12)

    paired = pd.read_csv(tables / "paired_male_turnover_bootstrap_summary.csv")
    assert int(paired.loc[0, "n_individuals"]) == 216
    assert int(paired.loc[0, "bootstrap_replicates"]) == 50000
    assert np.isclose(
        paired.loc[0, "delta_y_minus_mt"], -0.10638088568760834
    )
    assert paired.loc[0, "delta_ci_low"] < 0 < paired.loc[0, "delta_ci_high"]
    assert "two_sided_tail_probability" not in paired
    assert "bootstrap_two_sided_sign_tail_probability" not in paired
    assert "bootstrap_sign_tail_interpretation" not in paired

    resolution = pd.read_csv(tables / "paired_male_y_resolution_sensitivity.csv")
    assert list(resolution["y_encoding"]) == [
        "broad_L1", "AADR_ISOGG_prefix_family"
    ]
    assert np.isclose(resolution.loc[0, "delta_y_minus_mt"], -0.10638088568760834)
    assert np.isclose(resolution.loc[1, "delta_y_minus_mt"], -0.0112922946956559)
    assert np.isclose(resolution.loc[1, "delta_change_vs_broad"], 0.0950885909919523)
    assert np.isclose(resolution.loc[1, "delta_change_bootstrap_median"], 0.0928571428571429)
    assert np.isclose(resolution.loc[1, "delta_change_ci_low"], 0.0233704723536657)
    assert np.isclose(resolution.loc[1, "delta_change_ci_high"], 0.1760522471849909)
    assert resolution.loc[1, "delta_change_ci_low"] > 0
    assert "not a phylogenetic re-call" in resolution.loc[1, "mapping_caveat"]

    callability = pd.read_csv(tables / "marker_callability_tests.csv")
    assert set(callability["p_value_method"]) == {"fixed_margin_monte_carlo"}
    assert set(callability["monte_carlo_resamples"]) == {99999}
    observed_callability = dict(zip(callability["marker"], callability["monte_carlo_p"]))
    assert np.isclose(observed_callability["mtDNA"], 0.00341)
    assert np.isclose(observed_callability["Y"], 0.00530)
    assert set(callability["n_expected_lt_5"]) == {3, 8}

    date_summary = pd.read_csv(tables / "date_uncertainty_summary.csv")
    assert set(date_summary["scenario_draws"]) == {5000}
    assert set(date_summary["analysis_type"]) == {
        "chronological_bin_assignment_scenario_sensitivity"
    }
    assert date_summary["shared_individual_date_scenarios_across_markers"].all()
    assert (~date_summary["is_calibrated_date_posterior"]).all()

    tests = pd.read_csv(tables / "global_composition_cluster_tests.csv")
    assert set(tests["marker"]) == {"mtDNA", "Y"}
    for column in (
        "cluster_wild_p",
        "holm_cluster_wild_p",
        "repeated_site_permutation_p",
        "holm_repeated_site_p",
    ):
        assert ((tests[column] > 0) & (tests[column] <= 1)).all()
    assert (tests["n_site_clusters"] < tests["n_site_period_profiles"]).all()
    assert (tests["repeated_sites"] > 0).all()
    assert set(tests["leverage_adjustment"]) == {"HC2"}
    assert not (tables / "global_composition_permutation_tests.csv").exists()
    mt = tests.loc[tests["marker"] == "mtDNA"].iloc[0]
    assert int(mt["n_site_clusters"]) == 132
    assert int(mt["repeated_site_profiles"]) == 69
    assert int(mt["repeated_sites"]) == 31
    assert np.isclose(mt["repeated_site_permutation_p"], 0.1194)
    assert np.isclose(mt["holm_repeated_site_p"], 0.2388)

    dispersion = pd.read_csv(
        tables / "composition_dispersion_cluster_tests.csv"
    )
    assert set(dispersion["marker"]) == {"mtDNA", "Y"}
    mt_dispersion = dispersion.loc[dispersion["marker"] == "mtDNA"].iloc[0]
    assert int(mt_dispersion["n_site_clusters"]) == 127
    assert np.isclose(mt_dispersion["holm_cluster_wild_p"], 0.0275)
    diagnostics = pd.read_csv(
        tables / "cluster_model_diagnostic_summary.csv"
    )
    assert set(diagnostics["marker"]) == {"mtDNA", "Y"}
    assert (diagnostics["max_full_leverage"] < 1).all()

    sensitivity = pd.read_csv(
        tables / "global_composition_sensitivity_tests.csv"
    )
    assert len(sensitivity) == 57
    assert set(sensitivity["permutations"]) == {1999}
    invalid = sensitivity[~sensitivity["inference_valid"]]
    assert len(invalid) == 1
    assert invalid.iloc[0]["marker"] == "Y"
    assert invalid.iloc[0]["analysis"] == (
        "Site-period profiles with >= 5 calls"
    )
    assert pd.isna(invalid.iloc[0]["cluster_wild_p"])

    analysis_manifest = json.loads(
        (out / "analysis_manifest.json").read_text(encoding="utf-8")
    )
    sensitivity_manifest = json.loads(
        (out / "global_sensitivity_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    run_analysis = Path(__file__).with_name("run_analysis.py")
    run_sensitivity = Path(__file__).with_name(
        "run_global_sensitivities.py"
    )
    assert analysis_manifest["source_code_sha256"] == sha256(run_analysis)
    assert sensitivity_manifest["source_code_sha256"] == sha256(
        run_sensitivity
    )
    assert sensitivity_manifest["primary_analysis_code_sha256"] == sha256(
        run_analysis
    )
    sensitivity_input = sensitivity_manifest["analytical_input"]
    recorded_input = Path(sensitivity_input["file"])
    if not recorded_input.is_absolute():
        recorded_input = Path(__file__).resolve().parents[1] / recorded_input
    if args.input is not None:
        assert recorded_input.resolve() == args.input.resolve()
    assert sensitivity_input["public_project_coded_input"] is True
    frozen_input_sha256 = "692a69cf38cc736f96ea5aa6b3b15024a9a49d50c03ae3e1470a8dc475504cc3"
    assert sensitivity_input["sha256"] == frozen_input_sha256
    assert sha256(recorded_input) == frozen_input_sha256
    assert analysis_manifest["paired_bootstrap_replicates"] == 50000
    assert analysis_manifest["callability_fixed_margin_monte_carlo_resamples"] == 99999
    assert analysis_manifest["date_scenario_draws"] == 5000
    assert "date_draws" not in analysis_manifest
    assert (
        analysis_manifest["public_coordinate_policy"][
            "rounding_decimal_degrees"
        ]
        == 1
    )
    assert analysis_manifest["public_project_coded_input"]["used_for_this_run"]
    assert analysis_manifest["paired_surrogate_probability_policy"][
        "removed_from_user_facing_outputs"
    ]

    turnover_sensitivity = pd.read_csv(tables / "turnover_sensitivity.csv")
    mt_all = turnover_sensitivity.loc[
        (turnover_sensitivity["marker"] == "mtDNA")
        & (turnover_sensitivity["analysis"] == "All marker-qualified calls")
    ].iloc[0]
    assert int(mt_all["n_sites"]) == 132

    mt_site = pd.read_csv(tables / "composition_mtdna_site_balanced.csv")
    b7 = mt_site.loc[mt_site["analysis_bin"] == "B7 651-1000 CE"].iloc[0]
    assert np.isclose(b7["A"], 0.25)
    assert np.isclose(b7["D"], 0.25)

    print("All aggregate-release integrity checks passed.")


if __name__ == "__main__":
    main()
