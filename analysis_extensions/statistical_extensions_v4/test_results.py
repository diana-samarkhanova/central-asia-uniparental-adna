#!/usr/bin/env python3
"""Integrity and semantic checks for statistical extensions v4."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
    input_path = ROOT / manifest["input"]
    assert input_path.exists()
    assert sha256(input_path) == manifest["input_sha256"]
    assert sha256(HERE / "run_stat_extensions.py") == manifest["analysis_script_sha256"]
    assert sha256(HERE / "README.md") == manifest["readme_sha256"]
    assert sha256(HERE / "test_results.py") == manifest["test_script_sha256"]
    for name, expected in manifest["outputs"].items():
        path = RESULTS / name
        assert path.exists(), name
        assert sha256(path) == expected, name

    stable = pd.read_csv(RESULTS / "stable_profile_cluster_tests.csv")
    assert set(stable.minimum_calls) == {1, 2, 3}
    assert set(stable.marker) == {"mtDNA", "Y"}
    assert len(stable) == 6
    assert stable.valid.all()
    assert stable.raw_cluster_wild_p.between(1 / 10_000, 1).all()
    assert np.isfinite(stable[["pseudo_f", "partial_r2"]]).all().all()
    primary = stable[stable.minimum_calls.eq(2)].set_index("marker")
    assert primary.loc["mtDNA", "n_analyzed_profiles"] == 71
    assert primary.loc["Y", "n_analyzed_profiles"] == 48

    audit = pd.read_csv(RESULTS / "hierarchical_model_identifiability_audit.csv")
    assert (audit.rank_increment_period_after_site_publication <= audit.nominal_period_df).all()
    assert (
        audit.loc[audit.minimum_calls.eq(2), "rank_increment_period_after_site_publication"]
        < audit.loc[audit.minimum_calls.eq(2), "nominal_period_df"]
    ).all()

    predictive = pd.read_csv(RESULTS / "multinomial_predictive_sensitivity.csv")
    assert len(predictive) == 4 and predictive.valid.all()
    assert (
        predictive.site_bootstrap_delta_ci_low
        <= predictive.delta_log_loss_full_minus_reduced
    ).all()
    assert (
        predictive.delta_log_loss_full_minus_reduced
        <= predictive.site_bootstrap_delta_ci_high
    ).all()

    adjusted = pd.read_csv(RESULTS / "paired_male_tv_finite_sample_adjustments.csv")
    assert len(adjusted) == 28
    assert adjusted.observed_site_balanced_tv.between(0, 1).all()
    assert adjusted.noise_floor_subtracted_tv.between(0, 1).all()
    assert adjusted.parametric_bias_corrected_tv.between(0, 1).all()
    mt = adjusted[adjusted.marker.eq("mtDNA")].pivot(
        index="transition", columns="y_encoding", values="null_expected_tv_noise_floor"
    )
    assert np.allclose(mt.iloc[:, 0], mt.iloc[:, 1], atol=0, rtol=0)

    summary = pd.read_csv(RESULTS / "paired_male_tv_adjusted_cluster_bootstrap_summary.csv")
    assert len(summary) == 6
    broad_raw = summary[
        summary.y_encoding.eq("broad_L1") & summary.adjustment.eq("raw")
    ].iloc[0]
    broad_adjustment = adjusted[adjusted.y_encoding.eq("broad_L1")]
    assert np.isclose(
        broad_raw.mt_point_mean_tv,
        broad_adjustment.loc[
            broad_adjustment.marker.eq("mtDNA"), "observed_site_balanced_tv"
        ].mean(),
    )
    assert np.isclose(
        broad_raw.y_point_mean_tv,
        broad_adjustment.loc[
            broad_adjustment.marker.eq("Y"), "observed_site_balanced_tv"
        ].mean(),
    )
    assert np.isclose(broad_raw.mt_point_mean_tv, 0.517404, atol=1e-6)
    assert np.isclose(broad_raw.y_point_mean_tv, 0.411023, atol=1e-6)
    assert np.isclose(
        broad_raw.delta_point_y_minus_mt,
        broad_raw.y_point_mean_tv - broad_raw.mt_point_mean_tv,
    )
    for prefix in ["percentile", "basic", "median_centered"]:
        assert (summary[f"{prefix}_delta_ci_low"] <= summary[f"{prefix}_delta_ci_high"]).all()
    assert (summary.rejected_empty_period_fraction < 0.01).all()

    equal_comp = pd.read_csv(RESULTS / "equal_500y_composition.csv")
    equal_tv = pd.read_csv(RESULTS / "equal_500y_turnover.csv")
    assert len(equal_comp) == 40
    assert len(equal_tv) == 36
    assert equal_tv.tv.between(0, 1).all()
    assert (equal_tv.interval_years == 500).all()

    continuous = pd.read_csv(RESULTS / "continuous_time_cluster_tests.csv")
    assert len(continuous) == 12 and continuous.valid.all()
    assert np.isfinite(
        continuous[["pseudo_f", "partial_r2", "raw_cluster_wild_p"]]
    ).all().all()

    print("PASS: input/script/README and all declared output checksums")
    print("PASS: 6 stable-profile tests, 4 predictive sensitivities")
    print("PASS: 28 TV adjustments, 6 paired summaries, 12 continuous-time tests")
    print("PASS: equal-width 500-y composition and turnover tables")


if __name__ == "__main__":
    main()
