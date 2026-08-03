#!/usr/bin/env python3
"""Marker-wide temporal composition tests across pre-defined sensitivities."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_analysis import (
    named_rng,
    repeated_site_period_test,
    site_cluster_wild_period_test,
    site_profile_table,
    unrelated_subset,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--analysis-output", required=True, type=Path)
    p.add_argument("--permutations", type=int, default=1999)
    p.add_argument("--seed", type=int, default=20260726)
    return p.parse_args()


def figure_global_sensitivities(table: pd.DataFrame, path: Path) -> None:
    key_order = [
        "All marker-qualified calls",
        "AADR assessment-positive subset",
        "Exclude population-outlier labels",
        "Direct dates only",
        "One representative per <=2d kin component",
        "Kazakhstan only",
        "Leave out country: Kazakhstan",
        "Exclude sparse late bins B7-B8",
        "Site-period profiles with >= 2 calls",
        "Site-period profiles with >= 3 calls",
        "Site-period profiles with >= 5 calls",
        "HC3 leverage adjustment",
    ]
    short = {
        "All marker-qualified calls": "All profiles",
        "AADR assessment-positive subset": "AADR assessment-positive",
        "Exclude population-outlier labels": "Exclude outlier labels",
        "Direct dates only": "Direct dates only",
        "One representative per <=2d kin component": "One per ≤2d kin component",
        "Kazakhstan only": "Kazakhstan only",
        "Leave out country: Kazakhstan": "Leave out Kazakhstan",
        "Exclude sparse late bins B7-B8": "Exclude B7–B8",
        "Site-period profiles with >= 2 calls": "Profiles with ≥2 calls",
        "Site-period profiles with >= 3 calls": "Profiles with ≥3 calls",
        "Site-period profiles with >= 5 calls": "Profiles with ≥5 calls",
        "HC3 leverage adjustment": "HC3 adjustment",
    }
    colors = {"mtDNA": "#D95F02", "Y": "#1B9E77"}
    offsets = {"mtDNA": -0.16, "Y": 0.16}
    fig, axes = plt.subplots(
        1, 2, figsize=(13.5, 7.2), gridspec_kw={"width_ratios": [1.45, 1]}
    )
    ax = axes[0]
    key = table[table["analysis"].isin(key_order)].copy()
    for marker in ("mtDNA", "Y"):
        group = key[key["marker"] == marker].set_index("analysis")
        values = [
            group.loc[name, "cluster_wild_p"]
            if name in group.index
            else np.nan
            for name in key_order
        ]
        y = np.arange(len(key_order)) + offsets[marker]
        ax.scatter(values, y, s=48, color=colors[marker], label=marker, zorder=3)
        for row_index, name in enumerate(key_order):
            if name not in group.index:
                continue
            if not bool(group.loc[name, "inference_valid"]):
                ax.text(
                    0.92,
                    row_index + offsets[marker],
                    f"{marker}: not estimable",
                    ha="right",
                    va="center",
                    fontsize=7.5,
                    color="#666666",
                )
    ax.axvline(0.05, color="#8C2D24", linestyle="--", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_xlim(0.00035, 1.0)
    ax.set_yticks(np.arange(len(key_order)))
    ax.set_yticklabels([short[name] for name in key_order])
    ax.invert_yaxis()
    ax.set_xlabel("Exploratory cluster-wild P")
    ax.set_title("A  Design sensitivity")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1]
    lopo = table[
        table["analysis"].str.startswith("Leave out publication:")
    ].copy()
    for marker, y_center in (("mtDNA", 0), ("Y", 1)):
        values = np.sort(
            lopo.loc[lopo["marker"] == marker, "cluster_wild_p"].to_numpy()
        )
        jitter = np.linspace(-0.16, 0.16, len(values))
        ax.scatter(
            values,
            y_center + jitter,
            s=45,
            color=colors[marker],
            alpha=0.82,
            label=marker,
        )
    ax.axvline(0.05, color="#8C2D24", linestyle="--", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_xlim(0.00035, 1.0)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["mtDNA", "Y"])
    ax.set_xlabel("Leave-one-publication-out P")
    ax.set_title("B  Publication influence")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
    fig.suptitle(
        "Temporal association is sensitive to geography, profile size and source study",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    tables = args.analysis_output / "tables"
    df = pd.read_csv(
        tables / "aadr_primary_analysis_catalogue.csv", keep_default_na=False
    )
    for col in ("strict_qc", "population_outlier", "direct_date"):
        df[col] = df[col].astype(str).str.lower().eq("true")
    summary = json.loads(
        (args.analysis_output / "results_summary.json").read_text(encoding="utf-8")
    )
    rows = []
    base_filters = {
        "All marker-qualified calls": df,
        "AADR assessment-positive subset": df[df["strict_qc"]],
        "Exclude population-outlier labels": df[~df["population_outlier"]],
        "Direct dates only": df[df["direct_date"]],
    }
    for marker, call_col, pooled_col, categories in [
        ("mtDNA", "mt_call", "mt_l1_pooled", summary["mt_l1_categories"]),
        ("Y", "y_call", "y_l1_pooled", summary["y_l1_categories"]),
    ]:
        filters = dict(base_filters)
        filters["One representative per <=2d kin component"] = unrelated_subset(
            df, call_col
        )
        filters["Kazakhstan only"] = df[df["country"] == "Kazakhstan"]
        filters["Exclude sparse late bins B7-B8"] = df[
            ~df["analysis_bin"].isin(
                ["B7 651-1000 CE", "B8 1001-1500 CE"]
            )
        ]
        for country in sorted(df["country"].unique()):
            filters[f"Leave out country: {country}"] = df[
                df["country"] != country
            ]
        for publication in sorted(df["publication"].unique()):
            has_marker_call = (
                df.loc[df["publication"] == publication, pooled_col] != ""
            ).any()
            if not has_marker_call:
                continue
            filters[f"Leave out publication: {publication}"] = df[
                df["publication"] != publication
            ]
        repeated_sensitivity_names = set(base_filters) | {
            "One representative per <=2d kin component"
        }
        for name, subset in filters.items():
            table = site_profile_table(subset, pooled_col, categories, min_calls=1)
            result = site_cluster_wild_period_test(
                table,
                args.permutations,
                named_rng(args.seed, f"{marker}:{name}:cluster-wild"),
            )
            if name in repeated_sensitivity_names:
                result.update(
                    repeated_site_period_test(
                        table,
                        args.permutations,
                        named_rng(
                            args.seed,
                            f"{marker}:{name}:repeated-site-permutation",
                        ),
                    )
                )
            rows.append(
                {
                    "analysis": name,
                    "marker": marker,
                    "n_calls": int((subset[pooled_col] != "").sum()),
                    "n_sites": subset.loc[
                        subset[pooled_col] != "", "locality"
                    ].nunique(),
                    **result,
                    "permutations": args.permutations,
                }
            )
        for minimum_calls in (2, 3, 5):
            name = f"Site-period profiles with >= {minimum_calls} calls"
            table = site_profile_table(
                df, pooled_col, categories, min_calls=minimum_calls
            )
            result = site_cluster_wild_period_test(
                table,
                args.permutations,
                named_rng(args.seed, f"{marker}:{name}:cluster-wild"),
            )
            rows.append(
                {
                    "analysis": name,
                    "marker": marker,
                    "n_calls": int(
                        table["n_calls"].sum()
                    ),
                    "n_sites": int(table["site"].nunique()),
                    **result,
                    "permutations": args.permutations,
                }
            )
        table = site_profile_table(df, pooled_col, categories, min_calls=1)
        hc3_result = site_cluster_wild_period_test(
            table,
            args.permutations,
            named_rng(args.seed, f"{marker}:HC3:cluster-wild"),
            leverage_adjustment="HC3",
        )
        rows.append(
            {
                "analysis": "HC3 leverage adjustment",
                "marker": marker,
                "n_calls": int((df[pooled_col] != "").sum()),
                "n_sites": int(
                    df.loc[df[pooled_col] != "", "locality"].nunique()
                ),
                **hc3_result,
                "permutations": args.permutations,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(tables / "global_composition_sensitivity_tests.csv", index=False)
    manifest = {
        "seed": args.seed,
        "rng_streams": (
            "Stable SHA-256-named NumPy SeedSequence streams for each marker, "
            "sensitivity and resampling procedure"
        ),
        "cluster_wild_resamples_per_test": args.permutations,
        "repeated_site_permutations_per_test": args.permutations,
        "number_of_tests": len(out),
        "source_code_sha256": sha256(Path(__file__)),
        "primary_analysis_code_sha256": sha256(
            Path(__file__).with_name("run_analysis.py")
        ),
        "interpretation": (
            "Exploratory sensitivity analyses. Cluster-wild tests preserve "
            "country-locality dependence; repeated-site tests use site fixed "
            "effects and within-site period-label permutations."
        ),
    }
    (args.analysis_output / "global_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    figure_global_sensitivities(
        out, args.analysis_output / "figures" / "figure_5_sensitivity.png"
    )
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
