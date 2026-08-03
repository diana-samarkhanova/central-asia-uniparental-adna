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
    date_uncertainty,
    dispersion_distance_table,
    figure_diversity_turnover,
    holm_adjust,
    model_residual_diagnostics,
    named_rng,
    paired_marker_turnover_bootstrap,
    repeated_site_period_test,
    site_cluster_effect_jackknife,
    site_cluster_wild_period_test,
    site_profile_table,
    summarize_bootstrap,
    bootstrap_site_profiles,
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
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--paired-bootstrap", type=int, default=50000)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--date-draws", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260725)
    parser.add_argument("--save-draws", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out = args.analysis_output
    tables = out / "tables"
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    catalogue_path = tables / "aadr_primary_analysis_catalogue.csv"
    data = pd.read_csv(catalogue_path, keep_default_na=False)
    for column in ("strict_qc", "population_outlier", "direct_date"):
        data[column] = data[column].astype(str).str.lower().eq("true")
    summary_path = out / "results_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    categories = {
        "mtDNA": summary["mt_l1_categories"],
        "Y": summary["y_l1_categories"],
    }
    columns = {"mtDNA": "mt_l1_pooled", "Y": "y_l1_pooled"}

    diversity_parts = []
    turnover_parts = []
    global_rows = []
    dispersion_rows = []
    dispersion_summary_parts = []
    residual_parts = []
    residual_summaries = []
    bootstrap_diagnostics = {}

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
        diversity = summarize_bootstrap(
            diversity_draws, "q1", "analysis_bin"
        )
        diversity["marker"] = marker
        diversity_parts.append(diversity)
        turnover = summarize_bootstrap(
            turnover_draws, "tv", "transition"
        )
        turnover["marker"] = marker
        turnover_parts.append(turnover)

        profiles = site_profile_table(
            data, column, marker_categories, min_calls=1
        )
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
        )
    )
    if args.save_draws:
        paired_draws.to_csv(
            tables / "paired_male_turnover_bootstrap_draws.csv",
            index=False,
        )
    paired_summary.to_csv(
        tables / "paired_male_turnover_bootstrap_summary.csv", index=False
    )

    date_rows = []
    for marker in ("mtDNA", "Y"):
        draws = date_uncertainty(
            data,
            columns[marker],
            categories[marker],
            args.date_draws,
            named_rng(args.seed, f"{marker}:date-uncertainty"),
        )
        if args.save_draws:
            draws.to_csv(
                tables / f"date_uncertainty_{marker.lower()}_draws.csv",
                index=False,
            )
        date_rows.append(
            {
                "marker": marker,
                "draws": args.date_draws,
                "mean_adjacent_tv_median": draws[
                    "mean_adjacent_tv"
                ].median(),
                "ci_low": draws["mean_adjacent_tv"].quantile(0.025),
                "ci_high": draws["mean_adjacent_tv"].quantile(0.975),
            }
        )
    pd.DataFrame(date_rows).to_csv(
        tables / "date_uncertainty_summary.csv", index=False
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
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )

    manifest_path = out / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(
        {
            "code_revision_date": "2026-08-03",
            "random_seed": args.seed,
            "rng_streams": (
                "Stable SHA-256-named NumPy SeedSequence streams; each "
                "procedure is independent of unrelated random calls"
            ),
            "bootstrap_replicates": args.bootstrap,
            "paired_bootstrap_replicates": args.paired_bootstrap,
            "cluster_wild_resamples": args.permutations,
            "repeated_site_permutations": args.permutations,
            "date_draws": args.date_draws,
            "site_cluster_bootstrap": {
                "cluster": "country + locality",
                "draw": (
                    "one multinomial multiplicity per cluster, shared across "
                    "all periods and, in the paired analysis, both markers"
                ),
                "empty_period_handling": "reject complete draw and redraw",
                "marker_diagnostics": bootstrap_diagnostics,
                "paired_diagnostics": paired_diagnostics,
            },
            "recomputed_from": {
                "file": str(catalogue_path),
                "sha256": sha256(catalogue_path),
                "script": "analysis/recompute_from_catalogue.py",
                "script_sha256": sha256(Path(__file__)),
                "reason": (
                    "Corrected cross-period locality clustering and separated "
                    "named random streams; upstream catalogue unchanged"
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
