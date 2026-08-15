#!/usr/bin/env python3
"""Integrity checks for the released Central Asia synthesis outputs."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from run_analysis import (
    callability_table,
    major_haplogroup,
    valid_call,
    y_isogg_family_prefix,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "outputs/adna_central_asia_final")
    tables = out / "tables"
    summary = json.loads((out / "results_summary.json").read_text(encoding="utf-8"))
    catalogue = pd.read_csv(
        tables / "aadr_central_asia_unique_individual_catalogue.csv",
        keep_default_na=False,
    )
    primary = pd.read_csv(
        tables / "aadr_primary_analysis_catalogue.csv", keep_default_na=False
    )

    assert catalogue["individual_id"].is_unique
    assert len(catalogue) == 501
    assert len(primary) == 489
    assert (primary.loc[primary["y_call"] != "", "molecular_sex"] == "M").all()
    assert int((primary["mt_call"] != "").sum()) == 438
    assert int((primary["y_call"] != "").sum()) == 229
    assert summary["primary_sites"] == primary["locality"].nunique()
    hv = primary["mt_call"].str.upper().str.startswith("HV")
    assert hv.sum() == 36
    assert (primary.loc[hv, "mt_l1"] == "HV").all()
    assert (
        primary.loc[primary["y_call"].eq("IJK"), "y_l1"]
        == "Basal/unresolved"
    ).all()

    # Regression guard for the original missing-value parsing defect: a
    # non-call such as "n/a (<2x)" must not become mitochondrial lineage N.
    for missing in (
        "",
        ".",
        "NA",
        "n/a (<2x)",
        "unknown",
        None,
        pd.NA,
        pd.NaT,
        np.nan,
    ):
        assert not valid_call(missing)
        assert major_haplogroup(missing, "mt") == ""
    # Non-scalar inputs must be rejected cleanly rather than reaching an
    # ambiguous NumPy/Pandas truth-value test.
    for non_scalar in (
        np.array([np.nan]),
        np.array(["R1a", "R1b"], dtype=object),
        pd.Series([pd.NA]),
    ):
        assert not valid_call(non_scalar)
        assert major_haplogroup(non_scalar, "Y") == ""
    assert valid_call("N")
    assert major_haplogroup("N", "mt") == "N"
    assert y_isogg_family_prefix("R1a1a1") == "R1a"
    assert y_isogg_family_prefix("J2a1a4b") == "J2a"
    assert y_isogg_family_prefix("Q1b2b1b2~") == "Q1b"
    assert y_isogg_family_prefix("CF") == "Basal/unresolved"
    assert y_isogg_family_prefix("n/a (female)") == ""

    # Metadata and marker calls are selected independently across duplicate
    # AADR representations. I4773 is a frozen real-data regression case.
    i4773 = catalogue.loc[catalogue["individual_id"] == "I4773"].iloc[0]
    assert i4773["genetic_id"] == "I4773.SG"
    assert i4773["mt_call_genetic_id"] == "I4773.AG"
    assert i4773["mt_call"] == "U5a1a2a"
    marker_specific_mt = catalogue[
        (catalogue["mt_call_genetic_id"] != "")
        & (catalogue["mt_call_genetic_id"] != catalogue["genetic_id"])
    ]
    assert len(marker_specific_mt) == 10

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
        "mtDNA": np.array(
            [
                8.09304820062489,
                7.37924895726579,
                7.05598715030903,
                10.6233042205366,
                10.214034086894,
                7.28135866692849,
                6.54759028203901,
                5.27766344032807,
            ]
        ),
        "Y": np.array(
            [
                3.28516928489741,
                6.59103477659705,
                2.143236483694,
                3.46596156321367,
                4.27668004338971,
                4.77593764148247,
                2.60647903043411,
                4.33939062525723,
            ]
        ),
    }
    expected_turnover = {
        "mtDNA": np.array(
            [
                0.475190434565435,
                0.312396284271284,
                0.428058053058053,
                0.296958348556814,
                0.579528985507246,
                0.553125,
                0.332142857142857,
            ]
        ),
        "Y": np.array(
            [
                0.442063492063492,
                0.594907407407407,
                0.24537037037037,
                0.196623093681917,
                0.291316526610644,
                0.648809523809524,
                0.638888888888889,
            ]
        ),
    }
    diversity_summary = pd.read_csv(
        tables / "diversity_site_bootstrap_summary.csv"
    )
    turnover_summary = pd.read_csv(
        tables / "turnover_site_bootstrap_summary.csv"
    )
    for frame in (diversity_summary, turnover_summary):
        assert "bootstrap_median" in frame
        assert np.isfinite(frame["bootstrap_median"]).all()
        assert np.isfinite(frame["estimate"]).all()
    for marker in ("mtDNA", "Y"):
        observed_diversity = diversity_summary.loc[
            diversity_summary["marker"] == marker, "estimate"
        ].to_numpy(float)
        observed_turnover = turnover_summary.loc[
            turnover_summary["marker"] == marker, "estimate"
        ].to_numpy(float)
        assert np.allclose(
            observed_diversity, expected_diversity[marker], atol=1e-12
        )
        assert np.allclose(
            observed_turnover, expected_turnover[marker], atol=1e-12
        )
    assert not np.allclose(
        diversity_summary["estimate"],
        diversity_summary["bootstrap_median"],
    )
    assert not np.allclose(
        turnover_summary["estimate"],
        turnover_summary["bootstrap_median"],
    )

    crosswalk = pd.read_csv(tables / "cross_database_exact_id_audit.csv")
    assert len(crosswalk) == 119
    assert (
        (crosswalk["source"] == "aYChr-DB v5")
        & (crosswalk["match_status"] == "not_exactly_matched")
    ).sum() == 0
    legacy = pd.read_csv(tables / "extended_legacy_mtdna_not_in_aadr.csv")
    assert len(legacy) == 19
    assert legacy["mt_hg"].notna().all()

    paired = pd.read_csv(tables / "paired_male_turnover_bootstrap_summary.csv")
    paired_draws = pd.read_csv(
        tables / "paired_male_turnover_bootstrap_draws.csv"
    )
    assert int(paired.loc[0, "n_individuals"]) == 216
    assert int(paired.loc[0, "bootstrap_replicates"]) == 50000
    assert np.isclose(
        paired.loc[0, "delta_y_minus_mt"], -0.10638088568760834
    )
    assert paired.loc[0, "delta_ci_low"] < 0 < paired.loc[0, "delta_ci_high"]
    assert "two_sided_tail_probability" not in paired
    sign_tail_column = "bootstrap_two_sided_sign_tail_probability"
    assert sign_tail_column in paired
    assert "not a null-hypothesis p-value" in paired.loc[
        0, "bootstrap_sign_tail_interpretation"
    ]
    expected_sign_tail = min(
        1.0,
        2
        * min(
            float((paired_draws["delta_y_minus_mt"] <= 0).mean()),
            float((paired_draws["delta_y_minus_mt"] >= 0).mean()),
        ),
    )
    assert np.isclose(paired.loc[0, sign_tail_column], expected_sign_tail)
    assert np.isclose(
        paired.loc[0, "bootstrap_median_delta"],
        paired_draws["delta_y_minus_mt"].median(),
    )
    assert np.isclose(
        paired.loc[0, "delta_ci_low"],
        paired_draws["delta_y_minus_mt"].quantile(0.025),
    )
    assert np.isclose(
        paired.loc[0, "delta_ci_high"],
        paired_draws["delta_y_minus_mt"].quantile(0.975),
    )

    resolution = pd.read_csv(
        tables / "paired_male_y_resolution_sensitivity.csv"
    )
    assert list(resolution["y_encoding"]) == [
        "broad_L1",
        "AADR_ISOGG_prefix_family",
    ]
    assert (resolution["n_individuals"] == 216).all()
    assert (resolution["bootstrap_replicates"] == 50000).all()
    assert resolution.loc[0, "category_count"] < resolution.loc[
        1, "category_count"
    ]
    assert "not a phylogenetic re-call" in resolution.loc[
        1, "mapping_caveat"
    ]
    assert np.isclose(resolution.loc[0, "delta_change_vs_broad"], 0.0)
    assert not np.isclose(
        resolution.loc[1, "delta_y_minus_mt"],
        resolution.loc[0, "delta_y_minus_mt"],
    )
    assert np.isclose(
        resolution.loc[1, "y_mean_adjacent_tv"],
        0.5061118892001245,
    )
    assert np.isclose(
        resolution.loc[1, "delta_y_minus_mt"],
        -0.0112922946956559,
    )
    assert np.isclose(
        resolution.loc[1, "delta_change_vs_broad"],
        0.0950885909919523,
    )
    family_draws = pd.read_csv(
        tables / "paired_male_y_resolution_family_bootstrap_draws.csv"
    )
    assert resolution.loc[0, "n_individuals"] == resolution.loc[
        1, "n_individuals"
    ]
    assert resolution.loc[0, "n_sites"] == resolution.loc[1, "n_sites"]
    assert np.array_equal(
        paired_draws["replicate"].to_numpy(),
        family_draws["replicate"].to_numpy(),
    )
    assert np.array_equal(
        paired_draws["mt_mean_adjacent_tv"].to_numpy(),
        family_draws["mt_mean_adjacent_tv"].to_numpy(),
    )
    assert np.isclose(
        resolution.loc[1, "bootstrap_median_delta"],
        family_draws["delta_y_minus_mt"].median(),
    )
    assert np.isclose(
        resolution.loc[1, "delta_ci_low"],
        family_draws["delta_y_minus_mt"].quantile(0.025),
    )
    assert np.isclose(
        resolution.loc[1, "delta_ci_high"],
        family_draws["delta_y_minus_mt"].quantile(0.975),
    )
    for column in (
        "delta_change_bootstrap_median",
        "delta_change_ci_low",
        "delta_change_ci_high",
    ):
        assert np.isclose(resolution.loc[0, column], 0.0)
    delta_change_draws = (
        family_draws["delta_y_minus_mt"]
        - paired_draws["delta_y_minus_mt"]
    )
    assert np.isclose(
        resolution.loc[1, "delta_change_bootstrap_median"],
        delta_change_draws.median(),
    )
    assert np.isclose(
        resolution.loc[1, "delta_change_ci_low"],
        delta_change_draws.quantile(0.025),
    )
    assert np.isclose(
        resolution.loc[1, "delta_change_ci_high"],
        delta_change_draws.quantile(0.975),
    )

    callability_tests = pd.read_csv(
        tables / "marker_callability_tests.csv"
    )
    assert set(callability_tests["marker"]) == {"mtDNA", "Y"}
    assert callability_tests["sparse_expected_cells"].all()
    assert set(callability_tests["p_value_method"]) == {
        "fixed_margin_monte_carlo"
    }
    assert (callability_tests["min_expected_count"] < 5).all()
    assert set(callability_tests["monte_carlo_resamples"]) == {99999}
    assert (
        (callability_tests["monte_carlo_p"] > 0)
        & (callability_tests["monte_carlo_p"] <= 1)
    ).all()
    assert np.allclose(
        callability_tests["p_value"], callability_tests["monte_carlo_p"]
    )
    # Named RNG streams make the Monte Carlo result reproducible.
    _, callability_a = callability_table(
        primary, monte_carlo_resamples=999, seed=8102026
    )
    _, callability_b = callability_table(
        primary, monte_carlo_resamples=999, seed=8102026
    )
    assert np.array_equal(
        callability_a["monte_carlo_p"].to_numpy(),
        callability_b["monte_carlo_p"].to_numpy(),
    )

    date_summary = pd.read_csv(tables / "date_uncertainty_summary.csv")
    assert set(date_summary["scenario_draws"]) == {5000}
    assert set(date_summary["analysis_type"]) == {
        "chronological_bin_assignment_scenario_sensitivity"
    }
    assert date_summary["shared_individual_date_scenarios_across_markers"].all()
    assert (~date_summary["is_calibrated_date_posterior"]).all()
    expected_date_draw_columns = [
        "scenario_draw",
        "analysis_type",
        "date_scenario_model",
        "is_calibrated_date_posterior",
        "mean_adjacent_tv",
        "marker",
    ]
    for marker, filename in (
        ("mtDNA", "date_uncertainty_mtdna_draws.csv"),
        ("Y", "date_uncertainty_y_draws.csv"),
    ):
        date_draws = pd.read_csv(tables / filename)
        assert list(date_draws.columns) == expected_date_draw_columns
        assert len(date_draws) == 5000
        assert set(date_draws["marker"]) == {marker}

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

    dispersion = pd.read_csv(
        tables / "composition_dispersion_cluster_tests.csv"
    )
    assert set(dispersion["marker"]) == {"mtDNA", "Y"}
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
    assert analysis_manifest["paired_bootstrap_replicates"] == 50000
    assert (
        analysis_manifest[
            "callability_fixed_margin_monte_carlo_resamples"
        ]
        == 99999
    )
    assert analysis_manifest["date_scenario_draws"] == 5000
    assert "date_draws" not in analysis_manifest
    assert not analysis_manifest["date_assignment_scenario_sensitivity"].get(
        "is_calibrated_date_posterior", False
    )
    assert (
        analysis_manifest["public_coordinate_policy"][
            "rounding_decimal_degrees"
        ]
        == 1
    )

    for frame in (catalogue, primary, legacy):
        for column in ("latitude", "longitude"):
            values = pd.to_numeric(frame[column], errors="coerce").dropna()
            assert np.allclose(values, values.round(1), atol=1e-12)

    print("All integrity checks passed.")


if __name__ == "__main__":
    main()
