#!/usr/bin/env python3
"""Recompute statistical outputs from the frozen analytical catalogue.

This is the license-safe second reproducibility tier. It is also used to
regenerate resampling outputs after changes that do not affect AADR
deduplication, haplogroup harmonization or database crosswalking.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from run_analysis import (
    COUNTRIES,
    SITE_LOCALITY_ALIASES,
    add_y_resolution_sensitivity_encoding,
    analysis_cell_adequacy,
    apply_site_locality_aliases,
    callability_table,
    count_matrix,
    date_uncertainty,
    dispersion_distance_table,
    figure_composition,
    figure_diversity_turnover,
    figure_sampling,
    holm_adjust,
    mean_adjacent_tv,
    model_residual_diagnostics,
    named_rng,
    observed_profile_statistics,
    paired_marker_turnover_bootstrap,
    profile,
    repeated_site_period_test,
    site_dominance,
    site_cluster_effect_jackknife,
    site_cluster_wild_period_test,
    site_profile_table,
    summarize_bootstrap,
    bootstrap_site_profiles,
    unrelated_subset,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-output", required=True, type=Path)
    parser.add_argument(
        "--catalogue",
        type=Path,
        help=(
            "Analytical input CSV. Defaults to the extended-package primary "
            "catalogue; the public route passes data/derived/"
            "central_asia_analysis_input_v1.csv."
        ),
    )
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--paired-bootstrap", type=int, default=50000)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--callability-resamples", type=int, default=99999)
    parser.add_argument("--date-draws", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260725)
    parser.add_argument("--save-draws", action="store_true")
    parser.add_argument(
        "--aggregate-only",
        action="store_true",
        help="Do not write coded site-profile/residual/dispersion row tables.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out = args.analysis_output
    tables = out / "tables"
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    catalogue_path = (
        args.catalogue
        if args.catalogue is not None
        else tables / "aadr_primary_analysis_catalogue.csv"
    )
    data = pd.read_csv(catalogue_path, keep_default_na=False)
    public_input = "project_record_id" in data.columns
    if public_input:
        expected_public_columns = {
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
        }
        if set(data.columns) != expected_public_columns:
            raise ValueError(
                "Public analytical input schema mismatch: "
                f"{sorted(set(data.columns) ^ expected_public_columns)}"
            )
        data = data.rename(
            columns={
                "project_record_id": "individual_id",
                "site_key": "locality",
                "study_key": "publication",
                "mt_category": "mt_l1_pooled",
                "y_category": "y_l1_pooled",
                "y_prefix_category": "y_isogg_prefix_family_pooled",
                "latitude_0_1deg": "latitude",
                "longitude_0_1deg": "longitude",
            }
        )
        for flag in (
            "mt_called",
            "y_called",
            "strict_qc",
            "population_outlier",
            "direct_date",
            "kin_representative_mt",
            "kin_representative_y",
        ):
            data[flag] = data[flag].astype(str).str.lower().eq("true")
        data["mt_call"] = data["mt_l1_pooled"].where(data["mt_called"], "")
        data["y_call"] = data["y_l1_pooled"].where(data["y_called"], "")
    else:
        data = apply_site_locality_aliases(data)
        data.to_csv(catalogue_path, index=False)
        complete_catalogue_path = (
            tables / "aadr_central_asia_unique_individual_catalogue.csv"
        )
        if complete_catalogue_path.exists():
            complete_catalogue = pd.read_csv(
                complete_catalogue_path, keep_default_na=False
            )
            apply_site_locality_aliases(complete_catalogue).to_csv(
                complete_catalogue_path, index=False
            )
    for column in ("strict_qc", "population_outlier", "direct_date"):
        data[column] = data[column].astype(str).str.lower().eq("true")
    for column in ("date_bp", "date_sd_bp", "latitude", "longitude"):
        data[column] = pd.to_numeric(data[column], errors="coerce")
    summary_path = out / "results_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    categories = {
        "mtDNA": summary["mt_l1_categories"],
        "Y": summary["y_l1_categories"],
    }
    columns = {"mtDNA": "mt_l1_pooled", "Y": "y_l1_pooled"}
    if public_input:
        resolution_data = data.copy()
        observed_y_family_categories = (
            resolution_data.loc[
                resolution_data["y_isogg_prefix_family_pooled"] != "",
                "y_isogg_prefix_family_pooled",
            ]
            .value_counts()
            .sort_values(ascending=False)
            .index.tolist()
        )
        previous_resolution = summary.get(
            "paired_male_y_resolution_sensitivity", []
        )
        previous_family = next(
            (
                row
                for row in previous_resolution
                if row.get("y_encoding") == "AADR_ISOGG_prefix_family"
            ),
            None,
        )
        y_family_categories = (
            str(previous_family["categories"]).split(";")
            if previous_family and previous_family.get("categories")
            else observed_y_family_categories
        )
        if set(y_family_categories) != set(observed_y_family_categories):
            raise ValueError("Public Y-prefix category set does not match summary")
        y_family_source_column = "project-coded y_prefix_category"
    else:
        resolution_data, y_family_categories, y_family_source_column = (
            add_y_resolution_sensitivity_encoding(data, 5)
        )

    # The derived-data route must refresh every table whose site counts or
    # site-balanced values depend on the normalized locality key.  Earlier
    # versions refreshed the inferential tables but accidentally left these
    # descriptive tables at their pre-alias values.
    count_matrix(data).to_csv(tables / "counts_country_by_bin.csv", index=False)
    analysis_cell_adequacy(data).to_csv(
        tables / "country_period_marker_adequacy.csv", index=False
    )
    site_dominance(data).to_csv(
        tables / "site_and_country_dominance.csv", index=False
    )

    sensitivity_rows = []
    filters = {
        "All marker-qualified calls": data,
        "AADR assessment-positive subset": data[data["strict_qc"]],
        "Exclude population-outlier labels": data[~data["population_outlier"]],
        "Direct dates only": data[data["direct_date"]],
    }
    paired = data[
        data["molecular_sex"].astype(str).str.startswith("M")
        & (data["mt_call"] != "")
        & (data["y_call"] != "")
    ]
    for marker, call_col, pooled_col in (
        ("mtDNA", "mt_call", "mt_l1_pooled"),
        ("Y", "y_call", "y_l1_pooled"),
    ):
        marker_filters = dict(filters)
        representative_flag = (
            "kin_representative_mt"
            if marker == "mtDNA"
            else "kin_representative_y"
        )
        if representative_flag in data.columns:
            marker_filters["One representative per <=2d kin component"] = (
                data[data[representative_flag] & (data[pooled_col] != "")]
            )
        else:
            marker_filters["One representative per <=2d kin component"] = (
                unrelated_subset(data, call_col)
            )
        marker_filters["Male-paired subset"] = paired
        for name, subset in marker_filters.items():
            called = subset[subset[pooled_col] != ""]
            sensitivity_rows.append(
                {
                    "analysis": name,
                    "marker": marker,
                    "n_calls": len(called),
                    "n_sites": len(
                        called.drop_duplicates(["country", "locality"])
                    ),
                    "mean_adjacent_tv": mean_adjacent_tv(
                        subset, pooled_col, categories[marker]
                    ),
                }
            )
        for country in COUNTRIES:
            subset = data[data["country"] != country]
            called = subset[subset[pooled_col] != ""]
            sensitivity_rows.append(
                {
                    "analysis": f"Leave out country: {country}",
                    "marker": marker,
                    "n_calls": len(called),
                    "n_sites": len(
                        called.drop_duplicates(["country", "locality"])
                    ),
                    "mean_adjacent_tv": mean_adjacent_tv(
                        subset, pooled_col, categories[marker]
                    ),
                }
            )
    pd.DataFrame(sensitivity_rows).to_csv(
        tables / "turnover_sensitivity.csv", index=False
    )

    diversity_parts = []
    turnover_parts = []
    global_rows = []
    dispersion_rows = []
    dispersion_summary_parts = []
    residual_parts = []
    residual_summaries = []
    bootstrap_diagnostics = {}
    observed_site_profiles = {}

    for marker in ("mtDNA", "Y"):
        column = columns[marker]
        marker_categories = categories[marker]
        diversity_draws, turnover_draws, diagnostics = bootstrap_site_profiles(
            data,
            column,
            marker_categories,
            args.bootstrap,
            named_rng(args.seed, f"{marker}:site-cluster-bootstrap"),
        )
        bootstrap_diagnostics[marker] = diagnostics
        if args.save_draws:
            diversity_draws.to_csv(
                tables / f"bootstrap_diversity_{marker.lower()}_draws.csv",
                index=False,
            )
            turnover_draws.to_csv(
                tables / f"bootstrap_turnover_{marker.lower()}_draws.csv",
                index=False,
            )
        observed_site_profile = profile(
            data,
            column,
            marker_categories,
            site_balanced=True,
        )
        observed_site_profiles[marker] = observed_site_profile
        profile(data, column, marker_categories, site_balanced=False).to_csv(
            tables
            / f"composition_{marker.lower()}_individual_weighted.csv",
            index=False,
        )
        observed_site_profile.to_csv(
            tables / f"composition_{marker.lower()}_site_balanced.csv",
            index=False,
        )
        observed_diversity, observed_turnover = observed_profile_statistics(
            observed_site_profile, marker_categories
        )
        diversity = summarize_bootstrap(
            diversity_draws,
            "q1",
            "analysis_bin",
            observed_diversity,
        )
        diversity["marker"] = marker
        diversity_parts.append(diversity)
        turnover = summarize_bootstrap(
            turnover_draws,
            "tv",
            "transition",
            observed_turnover,
        )
        turnover["marker"] = marker
        turnover_parts.append(turnover)

        profiles = site_profile_table(
            data, column, marker_categories, min_calls=1
        )
        if not args.aggregate_only:
            profiles.to_csv(
                tables / f"site_profiles_{marker.lower()}.csv", index=False
            )
        result = site_cluster_wild_period_test(
            profiles,
            args.permutations,
            named_rng(args.seed, f"{marker}:primary-cluster-wild"),
        )
        result.update(site_cluster_effect_jackknife(profiles))
        result.update(
            repeated_site_period_test(
                profiles,
                args.permutations,
                named_rng(
                    args.seed, f"{marker}:repeated-site-permutation"
                ),
            )
        )
        result["marker"] = marker
        global_rows.append(result)

        residual_table, residual_summary = model_residual_diagnostics(
            profiles, marker
        )
        residual_parts.append(residual_table)
        residual_summaries.append(residual_summary)

        dispersion = dispersion_distance_table(profiles)
        if not args.aggregate_only:
            dispersion.to_csv(
                tables / f"dispersion_profiles_{marker.lower()}.csv", index=False
            )
        dispersion_result = site_cluster_wild_period_test(
            dispersion,
            args.permutations,
            named_rng(args.seed, f"{marker}:dispersion-cluster-wild"),
        )
        dispersion_result["marker"] = marker
        dispersion_rows.append(dispersion_result)
        dispersion_summary = (
            dispersion.groupby("analysis_bin", observed=True)["h_distance"]
            .agg(
                n_profiles="size",
                median="median",
                q1=lambda values: values.quantile(0.25),
                q3=lambda values: values.quantile(0.75),
            )
            .reset_index()
        )
        dispersion_summary["marker"] = marker
        dispersion_summary_parts.append(dispersion_summary)

    diversity_table = pd.concat(diversity_parts, ignore_index=True)
    turnover_table = pd.concat(turnover_parts, ignore_index=True)
    global_table = pd.DataFrame(global_rows)
    global_table["holm_cluster_wild_p"] = holm_adjust(
        global_table["cluster_wild_p"]
    )
    global_table["holm_repeated_site_p"] = holm_adjust(
        global_table["repeated_site_permutation_p"]
    )
    global_table = global_table.sort_values("marker").reset_index(drop=True)
    dispersion_table = pd.DataFrame(dispersion_rows)
    dispersion_table["holm_cluster_wild_p"] = holm_adjust(
        dispersion_table["cluster_wild_p"]
    )
    dispersion_table = dispersion_table.sort_values("marker").reset_index(
        drop=True
    )
    dispersion_summary_table = pd.concat(
        dispersion_summary_parts, ignore_index=True
    )
    residual_table = pd.concat(residual_parts, ignore_index=True)
    residual_summary_table = pd.DataFrame(residual_summaries)

    diversity_table.to_csv(
        tables / "diversity_site_bootstrap_summary.csv", index=False
    )
    turnover_table.to_csv(
        tables / "turnover_site_bootstrap_summary.csv", index=False
    )
    global_table.to_csv(
        tables / "global_composition_cluster_tests.csv", index=False
    )
    dispersion_table.to_csv(
        tables / "composition_dispersion_cluster_tests.csv", index=False
    )
    dispersion_summary_table.to_csv(
        tables / "composition_dispersion_by_period.csv", index=False
    )
    if not args.aggregate_only:
        residual_table.to_csv(
            tables / "cluster_model_residual_diagnostics.csv", index=False
        )
    residual_summary_table.to_csv(
        tables / "cluster_model_diagnostic_summary.csv", index=False
    )

    paired_draws, paired_summary, paired_diagnostics = (
        paired_marker_turnover_bootstrap(
            data,
            categories["mtDNA"],
            categories["Y"],
            args.paired_bootstrap,
            named_rng(args.seed, "paired-marker:site-cluster-bootstrap"),
            y_marker_col="y_l1_pooled",
            y_encoding_label="broad_L1",
        )
    )
    family_draws, family_summary, family_paired_diagnostics = (
        paired_marker_turnover_bootstrap(
            resolution_data,
            categories["mtDNA"],
            y_family_categories,
            args.paired_bootstrap,
            named_rng(args.seed, "paired-marker:site-cluster-bootstrap"),
            y_marker_col="y_isogg_prefix_family_pooled",
            y_encoding_label="AADR_ISOGG_prefix_family",
        )
    )
    paired_resolution_sensitivity = pd.concat(
        [paired_summary, family_summary], ignore_index=True
    )
    paired_resolution_sensitivity["category_count"] = [
        len(categories["Y"]),
        len(y_family_categories),
    ]
    paired_resolution_sensitivity["categories"] = [
        ";".join(categories["Y"]),
        ";".join(y_family_categories),
    ]
    paired_resolution_sensitivity["mapping_source"] = [
        "marker-specific AADR Y call, first-letter L1 encoding",
        f"{y_family_source_column}, first letter + integer + branch letter",
    ]
    paired_resolution_sensitivity["mapping_caveat"] = [
        "Broad descriptive encoding",
        (
            "Nomenclature-prefix sensitivity only; not a phylogenetic re-call "
            "and not uniform evolutionary depth across haplogroups"
        ),
    ]
    broad_delta = float(paired_summary.loc[0, "delta_y_minus_mt"])
    paired_resolution_sensitivity["delta_change_vs_broad"] = (
        paired_resolution_sensitivity["delta_y_minus_mt"] - broad_delta
    )
    shared_resolution_design = (
        int(paired_summary.loc[0, "n_individuals"])
        == int(family_summary.loc[0, "n_individuals"])
        and int(paired_summary.loc[0, "n_sites"])
        == int(family_summary.loc[0, "n_sites"])
        and paired_draws["replicate"].equals(family_draws["replicate"])
        and paired_draws["mt_mean_adjacent_tv"].equals(
            family_draws["mt_mean_adjacent_tv"]
        )
    )
    if not shared_resolution_design:
        raise RuntimeError(
            "Y-resolution sensitivity requires the same paired individuals, "
            "sites and cluster-bootstrap draws"
        )
    delta_change_draws = (
        family_draws["delta_y_minus_mt"]
        - paired_draws["delta_y_minus_mt"]
    )
    paired_resolution_sensitivity["delta_change_bootstrap_median"] = [
        0.0,
        delta_change_draws.median(),
    ]
    paired_resolution_sensitivity["delta_change_ci_low"] = [
        0.0,
        delta_change_draws.quantile(0.025),
    ]
    paired_resolution_sensitivity["delta_change_ci_high"] = [
        0.0,
        delta_change_draws.quantile(0.975),
    ]
    if args.save_draws:
        paired_draws.to_csv(
            tables / "paired_male_turnover_bootstrap_draws.csv",
            index=False,
        )
        family_draws.to_csv(
            tables / "paired_male_y_resolution_family_bootstrap_draws.csv",
            index=False,
        )
    surrogate_probability_columns = [
        "bootstrap_two_sided_sign_tail_probability",
        "bootstrap_sign_tail_interpretation",
    ]
    paired_summary = paired_summary.drop(
        columns=surrogate_probability_columns, errors="ignore"
    )
    paired_resolution_sensitivity = paired_resolution_sensitivity.drop(
        columns=surrogate_probability_columns, errors="ignore"
    )
    paired_summary.to_csv(
        tables / "paired_male_turnover_bootstrap_summary.csv", index=False
    )
    paired_resolution_sensitivity.to_csv(
        tables / "paired_male_y_resolution_sensitivity.csv", index=False
    )

    callability, callability_tests = callability_table(
        data,
        monte_carlo_resamples=args.callability_resamples,
        seed=args.seed,
    )
    callability.to_csv(
        tables / "marker_callability_by_bin.csv", index=False
    )
    callability_tests.to_csv(
        tables / "marker_callability_tests.csv", index=False
    )

    date_rows = []
    for marker in ("mtDNA", "Y"):
        draws = date_uncertainty(
            data,
            columns[marker],
            categories[marker],
            args.date_draws,
            named_rng(args.seed, "shared:date-assignment-scenarios"),
        )
        # Match the full run_analysis.py draw schema exactly.  Keeping the
        # marker in each file also prevents provenance from depending on its
        # filename alone.
        draws["marker"] = marker
        if args.save_draws:
            draws.to_csv(
                tables / f"date_uncertainty_{marker.lower()}_draws.csv",
                index=False,
            )
        date_rows.append(
            {
                "marker": marker,
                "analysis_type": (
                    "chronological_bin_assignment_scenario_sensitivity"
                ),
                "scenario_draws": args.date_draws,
                "observed_mean_adjacent_tv": mean_adjacent_tv(
                    data, columns[marker], categories[marker]
                ),
                "scenario_median": draws["mean_adjacent_tv"].median(),
                "scenario_interval_low": draws["mean_adjacent_tv"].quantile(
                    0.025
                ),
                "scenario_interval_high": draws["mean_adjacent_tv"].quantile(
                    0.975
                ),
                "shared_individual_date_scenarios_across_markers": True,
                "is_calibrated_date_posterior": False,
                "interpretation": (
                    "Assumption-based boundary sensitivity; not a calibrated-"
                    "date posterior interval"
                ),
            }
        )
    date_summary_table = pd.DataFrame(date_rows)
    date_summary_table.to_csv(
        tables / "date_uncertainty_summary.csv", index=False
    )
    if {"latitude", "longitude"}.issubset(data.columns):
        figure_sampling(
            data,
            count_matrix(data),
            figures / "figure_1_sampling.png",
        )
    figure_composition(
        observed_site_profiles["mtDNA"],
        categories["mtDNA"],
        "mtDNA",
        figures / "figure_2_mtdna_composition.png",
    )
    figure_composition(
        observed_site_profiles["Y"],
        categories["Y"],
        "Y chromosome",
        figures / "figure_3_y_composition.png",
    )
    figure_diversity_turnover(
        diversity_table,
        turnover_table,
        figures / "figure_4_diversity_turnover.png",
    )

    summary["global_composition_tests"] = global_table.to_dict("records")
    summary["composition_dispersion_diagnostics"] = dispersion_table.to_dict(
        "records"
    )
    summary["cluster_model_diagnostics"] = residual_summary_table.to_dict(
        "records"
    )
    summary["paired_male_turnover_comparison"] = paired_summary.to_dict(
        "records"
    )
    summary["paired_male_y_resolution_sensitivity"] = (
        paired_resolution_sensitivity.to_dict("records")
    )
    summary["callability_tests"] = callability_tests.to_dict("records")
    summary["date_assignment_scenario_sensitivity"] = (
        date_summary_table.to_dict("records")
    )
    summary["primary_sites"] = int(
        data[["country", "locality"]].drop_duplicates().shape[0]
    )
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )

    manifest_path = out / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("date_draws", None)
    manifest.setdefault("haplogroup_harmonization", {})[
        "Y_resolution_sensitivity"
    ] = (
        "AADR ISOGG call prefix: first letter + first integer + immediately "
        "following branch letter; nomenclature sensitivity only, not a "
        "phylogenetic re-call"
    )
    manifest.update(
        {
            "code_revision_date": "2026-08-21",
            "random_seed": args.seed,
            "rng_streams": (
                "Stable SHA-256-named NumPy SeedSequence streams; each "
                "procedure is independent of unrelated random calls"
            ),
            "bootstrap_replicates": args.bootstrap,
            "paired_bootstrap_replicates": args.paired_bootstrap,
            "cluster_wild_resamples": args.permutations,
            "repeated_site_permutations": args.permutations,
            "callability_fixed_margin_monte_carlo_resamples": (
                args.callability_resamples
            ),
            "date_scenario_draws": args.date_draws,
            "date_assignment_scenario_sensitivity": {
                "shared_across_markers": True,
                "direct_dates": "Normal in BP using catalogue mean and SD",
                "indirect_dates": (
                    "Uniform in BP over mean +/- sqrt(3) times catalogue SD"
                ),
                "interpretation": (
                    "Assumption-based chronological-bin boundary sensitivity; "
                    "not calibrated radiocarbon posterior uncertainty"
                ),
            },
            "paired_surrogate_probability_policy": {
                "removed_from_user_facing_outputs": True,
                "removed_fields": surrogate_probability_columns,
                "reason": (
                    "The ordinary uncentred bootstrap sign fraction is not a "
                    "null-imposed hypothesis-test P value."
                ),
                "recoverable_from_saved_draws_in_extended_package": True,
            },
            "public_project_coded_input": {
                "used_for_this_run": public_input,
                "contains_source_person_ids": False if public_input else None,
                "contains_terminal_haplogroup_calls": False if public_input else None,
                "site_identifier": "project code" if public_input else "normalized locality",
                "linkage_warning": (
                    "Project-coded, not anonymous; public attribute combinations "
                    "may be linkable to the AADR source."
                    if public_input
                    else None
                ),
            },
            "site_locality_normalization": {
                "cluster_key": ["country", "locality"],
                "raw_text_column": "locality_raw",
                "aliases": [
                    {
                        "country": country,
                        "source_locality": source,
                        "normalized_locality": normalized,
                    }
                    for (country, source), normalized in (
                        SITE_LOCALITY_ALIASES.items()
                    )
                ],
            },
            "site_cluster_bootstrap": {
                "cluster": "country + locality",
                "draw": (
                    "one multinomial multiplicity per cluster, shared across "
                    "all periods and, in the paired analysis, both markers"
                ),
                "empty_period_handling": "reject complete draw and redraw",
                "marker_diagnostics": bootstrap_diagnostics,
                "paired_diagnostics": paired_diagnostics,
                "paired_y_resolution_diagnostics": {
                    "broad_L1": paired_diagnostics,
                    "AADR_ISOGG_prefix_family": family_paired_diagnostics,
                    "shared_cluster_draws_across_encodings": True,
                    "warning": (
                        "Resolution comparison is an encoding sensitivity, "
                        "not a demographic or phylogenetic test"
                    ),
                },
            },
            "recomputed_from": {
                "file": str(catalogue_path),
                "sha256": sha256(catalogue_path),
                "script": "analysis/recompute_from_catalogue.py",
                "script_sha256": sha256(Path(__file__)),
                "reason": (
                    (
                        "Recomputed from the schema-locked project-coded AADR-"
                        "derived analytical input; source identifiers, locality "
                        "strings and terminal calls are absent"
                        if public_input
                        else "Normalized one verified Bestamak spelling alias "
                        "while preserving locality_raw"
                    )
                    + "; Figure 4 uses the observed site-balanced statistic as "
                    "its point and cluster-bootstrap percentiles as its interval; "
                    "date sensitivity uses 5,000 shared assignment scenarios"
                ),
            },
            "source_code_sha256": sha256(
                Path(__file__).with_name("run_analysis.py")
            ),
        }
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "global_tests": global_table.to_dict("records"),
                "paired": paired_summary.to_dict("records"),
                "bootstrap_diagnostics": bootstrap_diagnostics,
                "paired_diagnostics": paired_diagnostics,
            },
            indent=2,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
