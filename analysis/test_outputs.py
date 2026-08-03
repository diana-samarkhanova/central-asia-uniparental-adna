#!/usr/bin/env python3
"""Integrity checks for the released Central Asia synthesis outputs."""

from __future__ import annotations

import hashlib
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
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "outputs/adna_central_asia_final")
    tables = out / "tables"
    summary = json.loads((out / "results_summary.json").read_text(encoding="utf-8"))

    assert summary["aadr_catalogue_unique_archaeological_individuals"] == 501
    assert summary["primary_analysis_unique_individuals_3500BCE_to_1500CE"] == 489
    assert summary["primary_mt_calls"] == 438
    assert summary["primary_y_calls_in_molecular_males"] == 229
    assert summary["primary_sites"] == 137
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

    paired = pd.read_csv(tables / "paired_male_turnover_bootstrap_summary.csv")
    assert int(paired.loc[0, "n_individuals"]) == 216
    assert int(paired.loc[0, "bootstrap_replicates"]) == 50000
    assert np.isclose(
        paired.loc[0, "delta_y_minus_mt"], -0.10638088568760834
    )
    assert paired.loc[0, "delta_ci_low"] < 0 < paired.loc[0, "delta_ci_high"]

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
        analysis_manifest["public_coordinate_policy"][
            "rounding_decimal_degrees"
        ]
        == 1
    )

    print("All aggregate-release integrity checks passed.")


if __name__ == "__main__":
    main()
