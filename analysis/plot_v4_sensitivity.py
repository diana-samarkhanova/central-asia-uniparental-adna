#!/usr/bin/env python3
"""Render Figure 5 v4 from saved sensitivity tables without recalculating tests.

The all-profile, >=2-call and >=3-call rows use the statistical extension's
raw cluster-wild P values. Other design rows and leave-one-publication rows
use exploratory raw P values from the main analysis output. Input hashes,
plotted values and rendering provenance accompany the PNG.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COLORS = {"mtDNA": "#D55E00", "Y": "#009E73"}
OFFSETS = {"mtDNA": -0.12, "Y": 0.12}
DESIGNS = (
    ("All marker-qualified calls", "All profiles (archive)", 1),
    ("AADR assessment-positive subset", "AADR assessment-positive", None),
    ("Exclude population-outlier labels", "Exclude outlier labels", None),
    ("Direct dates only", "Direct dates only", None),
    ("One representative per <=2d kin component", "One per ≤2d kin component", None),
    ("Kazakhstan only", "Kazakhstan only", None),
    ("Leave out country: Kazakhstan", "Leave out Kazakhstan", None),
    ("Exclude sparse late bins B7-B8", "Exclude B7–B8", None),
    ("Site-period profiles with >= 2 calls", "Profiles with ≥2 calls\n(revised primary)", 2),
    ("Site-period profiles with >= 3 calls", "Profiles with ≥3 calls\n(stringent sensitivity)", 3),
    ("Site-period profiles with >= 5 calls", "Profiles with ≥5 calls", None),
    ("HC3 leverage adjustment", "HC3 adjustment", None),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recorded_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def is_valid(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if str(value).strip().lower() in {"true", "false"}:
        return str(value).strip().lower() == "true"
    raise ValueError(f"Invalid inference-valid flag: {value!r}")


def raw_p(row: pd.Series, column: str, valid_column: str) -> float | None:
    if not is_valid(row[valid_column]):
        return None
    value = float(row[column])
    if not np.isfinite(value) or not 0 < value <= 1:
        raise ValueError(f"Valid test must have a finite raw P in (0, 1]: {value}")
    return value


def one_row(table: pd.DataFrame, mask: pd.Series, description: str) -> pd.Series:
    subset = table.loc[mask]
    if len(subset) != 1:
        raise ValueError(f"Expected one {description} row; found {len(subset)}")
    return subset.iloc[0]


def count_text(values: pd.Series) -> str:
    counts = sorted({int(value) for value in values})
    return "/".join(f"{value:,}" for value in counts)


def plot_figure(global_tests: pd.DataFrame, stable_tests: pd.DataFrame,
                output: Path) -> dict[str, object]:
    design_points: list[dict[str, object]] = []
    for index, (name, label, minimum) in enumerate(DESIGNS):
        for marker in COLORS:
            if minimum is None:
                row = one_row(global_tests,
                              global_tests.analysis.eq(name) & global_tests.marker.eq(marker),
                              f"global {name} / {marker}")
                value = raw_p(row, "cluster_wild_p", "inference_valid")
                source = "global_composition_sensitivity_tests.csv"
                resamples = int(row["cluster_wild_resamples"])
            else:
                row = one_row(stable_tests,
                              stable_tests.minimum_calls.eq(minimum) & stable_tests.marker.eq(marker),
                              f"stable minimum {minimum} / {marker}")
                value = raw_p(row, "raw_cluster_wild_p", "valid")
                source = "stable_profile_cluster_tests.csv"
                resamples = int(row["wild_resamples"])
            design_points.append({"row": index, "analysis": name, "label": label,
                                  "marker": marker, "raw_p": value,
                                  "black_outline": minimum is not None,
                                  "resamples": resamples, "source": source})

    publication = global_tests.loc[
        global_tests.analysis.str.startswith("Leave out publication:")
    ].copy()
    if publication.empty:
        raise ValueError("No leave-one-publication sensitivity rows were found")
    if publication.duplicated(["analysis", "marker"]).any():
        raise ValueError("Duplicate publication/marker sensitivity rows")
    if not set(publication.marker).issubset(COLORS):
        raise ValueError("Unexpected marker in publication sensitivity rows")
    # Match the v4 figure's reverse publication order, with paired markers
    # aligned on a shared row. A missing marker row is not an inferred P value.
    study_order = sorted(publication.analysis.unique(), reverse=True)
    study_position = {name: index for index, name in enumerate(study_order)}
    publication_points: list[dict[str, object]] = []
    for _, row in publication.iterrows():
        publication_points.append({
            "row": study_position[row.analysis], "analysis": row.analysis,
            "marker": row.marker,
            "raw_p": raw_p(row, "cluster_wild_p", "inference_valid"),
            "resamples": int(row.cluster_wild_resamples),
            "source": "global_composition_sensitivity_tests.csv",
        })

    positive = [point["raw_p"] for point in design_points + publication_points
                if point["raw_p"] is not None]
    minimum_x = min(0.00002, min(positive) / 2)
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 14, "axes.labelsize": 11,
                         "savefig.facecolor": "white"}):
        fig, axes = plt.subplots(1, 2, figsize=(16.8, 10.2),
                                 gridspec_kw={"width_ratios": [1.45, 1]})
        # Separate footer space keeps explanatory notes outside both data axes.
        fig.subplots_adjust(left=0.225, right=0.985, bottom=0.205,
                            top=0.86, wspace=0.11)
        left, right = axes
        primary_index = next(index for index, (_, _, minimum) in enumerate(DESIGNS)
                             if minimum == 2)
        left.axhspan(primary_index - 0.46, primary_index + 0.46,
                     color="#E8F1F8", zorder=0)
        for point in design_points:
            marker = point["marker"]
            y = point["row"] + OFFSETS[marker]
            if point["raw_p"] is None:
                left.text(0.96, y, f"{marker}: not estimable", ha="right", va="center",
                          fontsize=9, color="#666666")
                continue
            left.scatter(point["raw_p"], y,
                         s=77 if point["black_outline"] else 51,
                         facecolor=COLORS[marker],
                         edgecolor="black" if point["black_outline"] else "none",
                         linewidth=1.2, zorder=3)
        left.set_ylim(len(DESIGNS) - 0.35, -0.7)
        left.set_yticks(np.arange(len(DESIGNS)))
        left.set_yticklabels([label for _, label, _ in DESIGNS])
        left.set_title("A  Design and profile-size sensitivity", loc="left", pad=10)
        left.set_xlabel("Raw cluster-wild P")
        for marker in COLORS:
            left.scatter([], [], s=77, color=COLORS[marker], edgecolor="black",
                         linewidth=1.2, label=marker)
        left.legend(frameon=False, loc="upper right", borderaxespad=0.3)

        for point in publication_points:
            if point["raw_p"] is None:
                continue
            right.scatter(point["raw_p"], point["row"] + OFFSETS[point["marker"]],
                          s=45, color=COLORS[point["marker"]], alpha=0.78,
                          edgecolor="none", zorder=3)
        right.set_ylim(len(study_order) - 0.35, -0.7)
        right.set_yticks([])
        right.set_title("B  Publication influence", loc="left", pad=10)
        right.set_xlabel("Leave-one-publication raw P (exploratory)")
        for axis in axes:
            axis.set_xscale("log")
            axis.set_xlim(minimum_x, 1)
            axis.axvline(0.05, color="#8C2D21", linestyle="--", linewidth=1.25)
            axis.grid(axis="x", color="#D9D9D9", linewidth=0.75)
            axis.spines[["top", "right"]].set_visible(False)
        right.spines["left"].set_visible(False)
        fig.suptitle("Temporal association depends on profile size, geography and source study",
                     x=0.52, y=0.967, fontsize=19, fontweight="bold")
        left_note = (
            f"Black-edged points: extension analysis, {count_text(stable_tests.wild_resamples)} resamples.\n"
            f"Other points: exploratory, {count_text(global_tests.cluster_wild_resamples)} resamples.\n"
            "Blue band: revised primary analysis (profiles with ≥2 calls)."
        )
        right_note = (
            "Each row omits one AADR publication label.\n"
            f"{count_text(publication.cluster_wild_resamples)} resamples; raw P values, not multiplicity-adjusted.\n"
            "Dashed line in both panels: P = 0.05."
        )
        fig.text(left.get_position().x0, 0.111, left_note, va="top", fontsize=9,
                 color="#555555", linespacing=1.55)
        fig.text(right.get_position().x0, 0.111, right_note, va="top", fontsize=9,
                 color="#555555", linespacing=1.55)
        fig.savefig(output, dpi=300)
        plt.close(fig)
    return {"design_points": design_points, "publication_points": publication_points,
            "p_values": "raw; no multiplicity adjustment applied in this figure",
            "reference_line_p": 0.05,
            "revised_primary_minimum_calls": 2,
            "colors": COLORS, "publication_order": study_order}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-output", type=Path, required=True,
                        help="Main output directory containing tables/")
    parser.add_argument("--extension-results", type=Path, required=True,
                        help="Directory containing stable_profile_cluster_tests.csv")
    parser.add_argument("--output", type=Path, required=True, help="Figure PNG path")
    args = parser.parse_args()
    global_path = (args.analysis_output.expanduser().resolve()
                   / "tables/global_composition_sensitivity_tests.csv")
    stable_path = (args.extension_results.expanduser().resolve()
                   / "stable_profile_cluster_tests.csv")
    output = args.output.expanduser().resolve()
    provenance_path = output.with_name(output.stem + ".provenance.json")
    if output.suffix.lower() != ".png":
        parser.error("--output must end in .png")
    for source in [global_path, stable_path]:
        if not source.is_file():
            parser.error(f"Missing input table: {source}")
    for destination in [output, provenance_path]:
        for source in [global_path, stable_path, Path(__file__).resolve()]:
            if destination == source or (destination.exists() and destination.samefile(source)):
                parser.error(f"Output would overwrite an input or generator: {destination}")
    global_tests, stable_tests = pd.read_csv(global_path), pd.read_csv(stable_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs = [{"path": recorded_path(path), "sha256": sha256(path)}
              for path in [global_path, stable_path]]
    plot_description = plot_figure(global_tests, stable_tests, output)
    provenance = {
        "figure": recorded_path(output), "figure_sha256": sha256(output),
        "generator": recorded_path(Path(__file__)), "generator_sha256": sha256(Path(__file__)),
        "inputs": inputs,
        "environment": {"python": platform.python_version(),
                        "matplotlib": matplotlib.__version__,
                        "numpy": np.__version__, "pandas": pd.__version__},
        "rendering": {"format": "PNG", "dpi": 300, "width_inches": 16.8,
                      "height_inches": 10.2, "notes_outside_data_axes": True},
        **plot_description,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False,
                                          allow_nan=False) + "\n", encoding="utf-8")
    print(f"Rendered {output}")
    print(f"Provenance {provenance_path}")


if __name__ == "__main__":
    main()
