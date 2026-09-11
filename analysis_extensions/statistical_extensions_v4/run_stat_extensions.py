#!/usr/bin/env python3
"""Statistical extensions for the Central Asia uniparental aDNA review.

This script is deliberately self-contained and reads the schema-locked,
project-coded AADR-derived analytical input included with the release. It does
not require source person identifiers, locality strings, terminal calls or
exact coordinates.

Extensions
----------
1. Stable-profile Hellinger analyses (>=2 calls as the revised primary
   estimand; >=3 calls as a stringent sparse-profile sensitivity).
2. An identifiability audit and a ridge-regularized multinomial predictive
   sensitivity with country, publication, and site intercept terms.
3. Two assumption-explicit finite-sample adjustments for paired-male total
   variation (TV), evaluated for broad-L1 and AADR ISOGG-prefix Y encodings.
4. Equal-width 500-y bins, descriptive TV per 100 years, and continuous-time
   Hellinger tests that do not interpret binned TV as a biological velocity.

The multinomial analysis is not presented as a fitted Dirichlet-multinomial
mixed model.  Site and period are nearly nested in this archive, especially
after excluding singleton profiles.  A nominal site-random-effect model would
therefore derive its period contrast mainly from distributional assumptions.
The penalized model below is a predictive sensitivity, not confirmatory
inference; the cluster-wild Hellinger analysis remains the inferential model.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import re
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler, SplineTransformer


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data/derived/central_asia_analysis_input_v1.csv"
OUT = Path(__file__).resolve().parent / "results"
SEED = 20260821
N_WILD = 9_999
N_PARAMETRIC = 10_000
N_CLUSTER_BOOT = 10_000
DEFAULT_INPUT = INPUT
DEFAULT_OUT = OUT
DEFAULT_SEED = SEED
DEFAULT_N_WILD = N_WILD
DEFAULT_N_PARAMETRIC = N_PARAMETRIC
DEFAULT_N_CLUSTER_BOOT = N_CLUSTER_BOOT
OUTPUT_NAMES = (
    "stable_profile_cluster_tests.csv",
    "stable_profile_vectors.csv",
    "hierarchical_model_identifiability_audit.csv",
    "multinomial_predictive_sensitivity.csv",
    "paired_male_tv_adjusted_cluster_bootstrap_draws_broad.csv",
    "paired_male_tv_adjusted_cluster_bootstrap_draws_isogg_prefix.csv",
    "paired_male_tv_finite_sample_adjustments.csv",
    "paired_male_tv_adjusted_cluster_bootstrap_summary.csv",
    "equal_500y_composition.csv",
    "equal_500y_turnover.csv",
    "original_bin_linear_scaled_tv.csv",
    "continuous_time_cluster_tests.csv",
)


def positive_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def manifest_path(path: Path) -> str:
    """Store repository paths portably, retaining absolute external paths."""
    path = path.resolve()
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def same_file(first: Path, second: Path) -> bool:
    """Detect lexical aliases, symlinks and existing hard links."""
    return first.resolve() == second.resolve() or (
        first.exists() and second.exists() and first.samefile(second)
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n", 1)[0],
        epilog=("With no arguments, outputs and README are regenerated at the "
                "archived locations. Use --outdir for an isolated run."),
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                        help="Schema-locked analytical CSV (default: release input)")
    parser.add_argument("--outdir", type=Path,
                        help="Results directory (default: archived results directory)")
    parser.add_argument("--readme-output", type=Path,
                        help="Generated README; defaults to OUTDIR/README.md when --outdir is given")
    parser.add_argument("--wild-resamples", type=positive_int, default=DEFAULT_N_WILD)
    parser.add_argument("--parametric-draws", type=positive_int, default=DEFAULT_N_PARAMETRIC)
    parser.add_argument("--paired-bootstrap", type=positive_int, default=DEFAULT_N_CLUSTER_BOOT)
    parser.add_argument("--seed", type=positive_int, default=DEFAULT_SEED)
    args = parser.parse_args(argv)
    explicit_outdir = args.outdir is not None
    args.input = args.input.expanduser().resolve()
    args.outdir = (args.outdir or DEFAULT_OUT).expanduser().resolve()
    args.readme_output = (
        args.readme_output
        or (args.outdir / "README.md" if explicit_outdir
            else Path(__file__).resolve().parent / "README.md")
    ).expanduser().resolve()
    if not args.input.is_file():
        parser.error(f"input is not a file: {args.input}")
    if args.outdir.exists() and not args.outdir.is_dir():
        parser.error(f"outdir is not a directory: {args.outdir}")
    if args.readme_output.exists() and not args.readme_output.is_file():
        parser.error(f"readme-output is not a file: {args.readme_output}")
    if args.readme_output == args.outdir or args.readme_output in args.outdir.parents:
        parser.error("readme-output must not be the results directory or its parent")
    result_paths = [args.outdir / name for name in (*OUTPUT_NAMES, "run_manifest.json")]
    destinations = [*result_paths, args.readme_output]
    protected = [args.input, Path(__file__).resolve(),
                 Path(__file__).resolve().parent / "test_results.py"]
    for destination in destinations:
        if any(same_file(destination, path) for path in protected):
            parser.error(f"output would overwrite input or a script: {destination}")
    if any(same_file(args.readme_output, path) for path in result_paths):
        parser.error("readme-output must not collide with a generated result")
    return args

BIN_LABELS = [
    "B1 3500-2501 BCE",
    "B2 2500-1801 BCE",
    "B3 1800-901 BCE",
    "B4 900-201 BCE",
    "B5 200 BCE-300 CE",
    "B6 301-650 CE",
    "B7 651-1000 CE",
    "B8 1001-1500 CE",
]

# These midpoints describe the manuscript's bins on the catalogue year_ce
# axis. They are used only for the explicitly labelled linear-scaling
# sensitivity and not as an estimate of a biological rate.
BIN_MIDPOINT = {
    "B1 3500-2501 BCE": (-3500 - 2501) / 2,
    "B2 2500-1801 BCE": (-2500 - 1801) / 2,
    "B3 1800-901 BCE": (-1800 - 901) / 2,
    "B4 900-201 BCE": (-900 - 201) / 2,
    "B5 200 BCE-300 CE": (-200 + 300) / 2,
    "B6 301-650 CE": (301 + 650) / 2,
    "B7 651-1000 CE": (651 + 1000) / 2,
    "B8 1001-1500 CE": (1001 + 1500) / 2,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def named_rng(name: str) -> np.random.Generator:
    payload = f"{SEED}:{name}".encode()
    child = int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")
    return np.random.default_rng(child)


def valid_call(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() not in {
        "",
        "..",
        "nan",
        "n/a",
        "n/a (female)",
        "n/a (sex unknown)",
    }


def y_isogg_family_prefix(call: object) -> str:
    if not valid_call(call):
        return ""
    normalized = re.sub(r"\s+", "", str(call).upper()).rstrip("~*")
    if normalized in {
        "BT",
        "CT",
        "CF",
        "F",
        "GHIJK",
        "HIJK",
        "IJK",
        "IJ",
        "K",
        "K2",
        "K2B",
    }:
        return "Basal/unresolved"
    match = re.match(r"^([A-Z])(\d+)([A-Z]?)", normalized)
    if match:
        branch = match.group(3).lower() if match.group(3) else ""
        return f"{match.group(1)}{match.group(2)}{branch}"
    match = re.match(r"^([A-Z])", normalized)
    return match.group(1) if match else "Other/unresolved"


def pool_rare(series: pd.Series, minimum: int = 5) -> tuple[pd.Series, list[str]]:
    counts = series[series != ""].value_counts()
    keep = counts[counts >= minimum].index.tolist()
    pooled = series.where(series.isin(keep), np.where(series.eq(""), "", "Other"))
    categories = (
        pooled[pooled != ""].value_counts().sort_values(ascending=False).index.tolist()
    )
    return pooled, categories


def total_variation(a: np.ndarray, b: np.ndarray) -> float:
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return math.nan
    return float(0.5 * np.abs(a - b).sum())


def holm(values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(values), dtype=float)
    out = np.full(values.shape, np.nan)
    finite_idx = np.flatnonzero(np.isfinite(values))
    if finite_idx.size == 0:
        return out
    order = finite_idx[np.argsort(values[finite_idx])]
    adjusted = np.maximum.accumulate(
        np.minimum(1.0, values[order] * np.arange(len(order), 0, -1))
    )
    out[order] = adjusted
    return out


def markdown_table(frame: pd.DataFrame, digits: int = 6) -> str:
    """Render a small DataFrame without the optional ``tabulate`` package."""
    def format_value(value: object) -> str:
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{digits}g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    columns = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(format_value(value) for value in row) + " |")
    return "\n".join(lines)


def design_country(table: pd.DataFrame) -> np.ndarray:
    country = pd.get_dummies(table["country"], drop_first=True).astype(float)
    return np.column_stack([np.ones(len(table)), country.to_numpy()])


def residual_maker(x: np.ndarray) -> np.ndarray:
    return np.eye(len(x)) - x @ np.linalg.pinv(x)


def sse(m: np.ndarray, y: np.ndarray) -> float:
    residual = m @ y
    return float(np.sum(residual**2))


def site_profile_table(
    data: pd.DataFrame,
    marker_col: str,
    categories: list[str],
    min_calls: int,
    time_col: str = "analysis_bin",
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    use = data[(data[time_col].astype(str) != "") & (data[marker_col] != "")].copy()
    for (country, locality, time_value), group in use.groupby(
        ["country", "locality", time_col], observed=True, dropna=False
    ):
        if len(group) < min_calls:
            continue
        counts = group[marker_col].value_counts()
        p = np.array([counts.get(category, 0) for category in categories], float)
        p /= p.sum()
        row: dict[str, object] = {
            "country": country,
            "site": locality,
            "time_group": str(time_value),
            "n_calls": len(group),
            "mean_year_ce": float(pd.to_numeric(group["year_ce"]).mean()),
            "publication": "|".join(sorted(group["publication"].astype(str).unique())),
        }
        row.update({f"h_{cat}": val for cat, val in zip(categories, np.sqrt(p))})
        row.update({f"p_{cat}": val for cat, val in zip(categories, p)})
        rows.append(row)
    return pd.DataFrame(rows)


def cluster_wild_test(
    table: pd.DataFrame,
    effect_matrix: np.ndarray,
    n_resamples: int,
    rng: np.random.Generator,
    label: str,
) -> dict[str, object]:
    """Null-imposed HC2 Rademacher wild bootstrap by country-locality."""
    source_n = len(table)
    country_counts = table["country"].value_counts()
    singleton = country_counts[country_counts < 2].index.tolist()
    keep = ~table["country"].isin(singleton)
    table = table.loc[keep].reset_index(drop=True)
    effect_matrix = np.asarray(effect_matrix)[np.asarray(keep)]
    # Concatenated marker-profile tables contain all-NA columns belonging to
    # the other marker. Retain only complete response dimensions for the
    # current marker subset.
    y_cols = [
        column
        for column in table
        if column.startswith("h_") and table[column].notna().all()
    ]
    y = table[y_cols].to_numpy(float)
    x0 = design_country(table)
    x1 = np.column_stack([x0, effect_matrix])
    m0 = residual_maker(x0)
    m1 = residual_maker(x1)
    s0 = sse(m0, y)
    s1 = sse(m1, y)
    rank0, rank1 = np.linalg.matrix_rank(x0), np.linalg.matrix_rank(x1)
    df_effect = int(rank1 - rank0)
    df_resid = int(len(table) - rank1)
    leverage = 1.0 - np.diag(m0)
    invalid: list[str] = []
    if df_effect <= 0:
        invalid.append("non-positive effect degrees of freedom")
    if df_resid < 5:
        invalid.append("fewer than five residual degrees of freedom")
    if s1 <= np.finfo(float).eps:
        invalid.append("zero full-model residual SSE")
    if leverage.max(initial=0.0) >= 1 - 1e-8:
        invalid.append("unit reduced-model leverage")
    observed = math.nan
    p_value = math.nan
    partial_r2 = (s0 - s1) / s0 if s0 > 0 else math.nan
    if not invalid:
        observed = ((s0 - s1) / df_effect) / (s1 / df_resid)
        fitted0 = y - m0 @ y
        residual0 = (m0 @ y) / np.sqrt(np.clip(1.0 - leverage, 1e-10, None))[:, None]
        cluster = (table["country"].astype(str) + "|||" + table["site"].astype(str)).to_numpy()
        _, cluster_index = np.unique(cluster, return_inverse=True)
        ge = 0
        for _ in range(n_resamples):
            signs = rng.choice((-1.0, 1.0), size=cluster_index.max() + 1)
            y_star = fitted0 + residual0 * signs[cluster_index, None]
            s0_star = sse(m0, y_star)
            s1_star = sse(m1, y_star)
            f_star = ((s0_star - s1_star) / df_effect) / (s1_star / df_resid)
            ge += f_star >= observed
        p_value = (ge + 1) / (n_resamples + 1)
    cluster_counts = (
        table["country"].astype(str) + "|||" + table["site"].astype(str)
    ).value_counts()
    return {
        "effect": label,
        "pseudo_f": observed,
        "partial_r2": partial_r2,
        "raw_cluster_wild_p": p_value,
        "df_effect": df_effect,
        "df_resid": df_resid,
        "n_input_profiles": source_n,
        "n_analyzed_profiles": len(table),
        "n_calls_in_analyzed_profiles": int(table["n_calls"].sum()),
        "n_sites": int(cluster_counts.size),
        "n_multitime_sites": int((cluster_counts > 1).sum()),
        "dropped_singleton_countries": "|".join(singleton),
        "wild_resamples": n_resamples,
        "valid": not invalid,
        "invalid_reason": "; ".join(invalid),
    }


def period_effect_matrix(table: pd.DataFrame, levels: list[str]) -> np.ndarray:
    cat = pd.Categorical(table["time_group"], categories=levels)
    return pd.get_dummies(cat, drop_first=True).astype(float).to_numpy()


def stable_profile_tests(
    data: pd.DataFrame, marker_spec: dict[str, tuple[str, list[str]]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tests: list[dict[str, object]] = []
    profiles: list[pd.DataFrame] = []
    for marker, (column, categories) in marker_spec.items():
        for minimum, role in [(1, "archival-scope context"), (2, "revised primary"), (3, "stringent sensitivity")]:
            table = site_profile_table(data, column, categories, minimum)
            table.insert(0, "marker", marker)
            table.insert(1, "minimum_calls", minimum)
            table.insert(2, "analysis_role", role)
            profiles.append(table)
            effect = period_effect_matrix(table, BIN_LABELS)
            result = cluster_wild_test(
                table,
                effect,
                N_WILD,
                named_rng(f"stable:{marker}:{minimum}"),
                "categorical period",
            )
            result.update(
                {
                    "marker": marker,
                    "minimum_calls": minimum,
                    "analysis_role": role,
                    "estimand": (
                        "country-adjusted difference among equally weighted "
                        f"site-period lineage profiles with at least {minimum} "
                        "marker calls; does not estimate population frequencies"
                    ),
                }
            )
            tests.append(result)
    out = pd.DataFrame(tests)
    # Family-wise control is shown within each threshold across the two markers.
    out["holm_within_threshold_p"] = np.nan
    for minimum, idx in out.groupby("minimum_calls").groups.items():
        out.loc[idx, "holm_within_threshold_p"] = holm(
            out.loc[idx, "raw_cluster_wild_p"]
        )
    return out, pd.concat(profiles, ignore_index=True)


def identifiability_audit(
    profiles: pd.DataFrame, marker: str, minimum: int
) -> dict[str, object]:
    table = profiles[
        (profiles["marker"] == marker) & (profiles["minimum_calls"] == minimum)
    ].copy()
    period = pd.get_dummies(
        pd.Categorical(table["time_group"], categories=BIN_LABELS), drop_first=True
    ).astype(float)
    country = pd.get_dummies(table["country"], drop_first=True).astype(float)
    publication = pd.get_dummies(table["publication"], drop_first=True).astype(float)
    site_id = table["country"].astype(str) + "|||" + table["site"].astype(str)
    site = pd.get_dummies(site_id, drop_first=True).astype(float)
    intercept = np.ones((len(table), 1))
    nuisance = np.column_stack([intercept, country, publication, site])
    full = np.column_stack([nuisance, period])
    site_periods = table.assign(site_id=site_id).groupby("site_id")["time_group"].nunique()
    publications = data_for_audit_publications(table)
    return {
        "marker": marker,
        "minimum_calls": minimum,
        "n_profiles": len(table),
        "n_sites": site_id.nunique(),
        "n_multiperiod_sites": int((site_periods > 1).sum()),
        "proportion_single_period_sites": float((site_periods == 1).mean()),
        "n_profile_publication_combinations": table["publication"].nunique(),
        "n_publications_after_splitting_combinations": publications,
        "rank_nuisance_country_publication_site": int(np.linalg.matrix_rank(nuisance)),
        "rank_increment_period_after_site_publication": int(
            np.linalg.matrix_rank(full) - np.linalg.matrix_rank(nuisance)
        ),
        "nominal_period_df": len(BIN_LABELS) - 1,
        "interpretation": (
            "A full site-intercept model identifies period chiefly from the "
            "small set of multi-period sites; site and period are otherwise nested."
        ),
    }


def data_for_audit_publications(table: pd.DataFrame) -> int:
    values: set[str] = set()
    for value in table["publication"].astype(str):
        values.update(value.split("|"))
    return len(values)


def weighted_multinomial_predictive_sensitivity(
    data: pd.DataFrame,
    marker: str,
    marker_col: str,
    minimum: int,
) -> dict[str, object]:
    """Nested site-group CV for a regularized multinomial sensitivity.

    Site dummy terms are L2-shrunk nuisance intercept proxies.  The model is
    evaluated only on held-out sites; therefore it cannot obtain an artificial
    advantage by memorizing a held-out site's outcome.  This is a predictive
    robustness check and does not furnish a mixed-model hypothesis test.
    """
    use = data[data[marker_col] != ""].copy()
    profile_size = use.groupby(["country", "locality", "analysis_bin"])[marker_col].transform("size")
    use = use[profile_size >= minimum].copy()
    use["site_id"] = use["country"].astype(str) + "|||" + use["locality"].astype(str)
    use["profile_id"] = use["site_id"] + "|||" + use["analysis_bin"].astype(str)
    sizes = use["profile_id"].value_counts()
    use["weight"] = use["profile_id"].map(lambda x: 1.0 / sizes[x])
    y = use[marker_col].astype(str).to_numpy()
    groups = use["site_id"].to_numpy()
    classes = np.unique(y)
    outer_splits = min(4, len(np.unique(groups)))
    if outer_splits < 3:
        return {
            "marker": marker,
            "minimum_calls": minimum,
            "valid": False,
            "invalid_reason": "fewer than three site clusters",
        }
    # A compact prespecified grid keeps this sensitivity reproducible and
    # avoids over-tuning a small, sparse archive.
    grid = [0.1, 1.0, 10.0]
    rows: list[dict[str, object]] = []

    def make_model(features: list[str], c_value: float) -> Pipeline:
        encoder = ColumnTransformer(
            [("categorical", OneHotEncoder(handle_unknown="ignore"), features)],
            remainder="drop",
        )
        classifier = LogisticRegression(
            C=c_value,
            l1_ratio=0,
            solver="lbfgs",
            max_iter=1_000,
            tol=1e-6,
        )
        return Pipeline([("encode", encoder), ("model", classifier)])

    def align_probability(model: Pipeline, probability: np.ndarray) -> np.ndarray:
        """Align fold-specific model classes to the global marker categories."""
        model_classes = model.named_steps["model"].classes_
        aligned = np.full((len(probability), len(classes)), 1e-15, dtype=float)
        global_index = {value: i for i, value in enumerate(classes)}
        for local_index, value in enumerate(model_classes):
            aligned[:, global_index[value]] = probability[:, local_index]
        aligned /= aligned.sum(axis=1, keepdims=True)
        return aligned

    feature_sets = {
        "reduced_country_publication_site": ["country", "publication", "site_id"],
        "full_plus_period": ["country", "publication", "site_id", "analysis_bin"],
    }
    outer = GroupKFold(n_splits=outer_splits)
    for fold, (train_idx, test_idx) in enumerate(outer.split(use, y, groups)):
        train = use.iloc[train_idx]
        test = use.iloc[test_idx]
        train_y = y[train_idx]
        test_y = y[test_idx]
        inner_groups = train["site_id"].to_numpy()
        inner_splits = min(3, len(np.unique(inner_groups)))
        for model_name, features in feature_sets.items():
            candidate_scores: list[tuple[float, float]] = []
            inner = GroupKFold(n_splits=inner_splits)
            for c_value in grid:
                fold_losses: list[float] = []
                for inner_train, inner_test in inner.split(train, train_y, inner_groups):
                    model = make_model(features, c_value)
                    model.fit(
                        train.iloc[inner_train][features],
                        train_y[inner_train],
                        model__sample_weight=train.iloc[inner_train]["weight"].to_numpy(),
                    )
                    probability = align_probability(
                        model, model.predict_proba(train.iloc[inner_test][features])
                    )
                    fold_losses.append(
                        log_loss(
                            train_y[inner_test],
                            probability,
                            labels=classes,
                            sample_weight=train.iloc[inner_test]["weight"].to_numpy(),
                        )
                    )
                candidate_scores.append((float(np.mean(fold_losses)), c_value))
            best_loss, best_c = min(candidate_scores)
            model = make_model(features, best_c)
            model.fit(
                train[features],
                train_y,
                model__sample_weight=train["weight"].to_numpy(),
            )
            probability = align_probability(model, model.predict_proba(test[features]))
            # Record weighted numerator/denominator so aggregation is exact.
            probability = np.clip(probability, 1e-15, 1.0)
            class_index = {value: i for i, value in enumerate(classes)}
            individual_loss = np.array(
                [-math.log(probability[i, class_index[value]]) for i, value in enumerate(test_y)]
            )
            for row_i, loss_i in zip(test_idx, individual_loss):
                rows.append(
                    {
                        "row_index": int(row_i),
                        "profile_id": use.iloc[row_i]["profile_id"],
                        "site_id": use.iloc[row_i]["site_id"],
                        "weight": float(use.iloc[row_i]["weight"]),
                        "model": model_name,
                        "outer_fold": fold,
                        "selected_c": best_c,
                        "inner_cv_loss": best_loss,
                        "negative_log_likelihood": float(loss_i),
                    }
                )
    losses = pd.DataFrame(rows)
    model_loss = (
        losses.assign(weighted=lambda x: x["weight"] * x["negative_log_likelihood"])
        .groupby("model")
        .apply(lambda x: x["weighted"].sum() / x["weight"].sum(), include_groups=False)
    )
    wide = losses.pivot(index="row_index", columns="model", values="negative_log_likelihood")
    weight = losses.drop_duplicates("row_index").set_index("row_index")["weight"]
    site = losses.drop_duplicates("row_index").set_index("row_index")["site_id"]
    delta = wide["full_plus_period"] - wide["reduced_country_publication_site"]
    site_delta = (
        pd.DataFrame({"delta": delta, "weight": weight, "site": site})
        .assign(weighted=lambda x: x.delta * x.weight)
        .groupby("site")
        .agg(weighted=("weighted", "sum"), weight=("weight", "sum"))
    )
    site_delta["mean"] = site_delta["weighted"] / site_delta["weight"]
    rng = named_rng(f"multinomial-bootstrap:{marker}:{minimum}")
    numerator = site_delta["weighted"].to_numpy()
    denominator = site_delta["weight"].to_numpy()
    sampled = rng.integers(0, len(site_delta), size=(20_000, len(site_delta)))
    # Ratio-of-sums preserves the profile-weighted estimand when resampling
    # whole sites, including sites observed in more than one period.
    boot = numerator[sampled].sum(axis=1) / denominator[sampled].sum(axis=1)
    observed_delta = float((delta * weight).sum() / weight.sum())
    return {
        "marker": marker,
        "minimum_calls": minimum,
        "valid": True,
        "n_individual_calls": len(use),
        "n_site_period_profiles": use["profile_id"].nunique(),
        "n_sites": use["site_id"].nunique(),
        "n_lineage_categories": len(classes),
        "reduced_heldout_site_log_loss": float(model_loss["reduced_country_publication_site"]),
        "full_period_heldout_site_log_loss": float(model_loss["full_plus_period"]),
        "delta_log_loss_full_minus_reduced": observed_delta,
        "site_bootstrap_delta_ci_low": float(np.quantile(boot, 0.025)),
        "site_bootstrap_delta_ci_high": float(np.quantile(boot, 0.975)),
        "outer_group_folds": outer_splits,
        "model_interpretation": (
            "Negative delta favors period for held-out-site prediction. Ridge-penalized "
            "site/publication terms are nuisance intercept proxies, not a fitted random-effects "
            "sampling model; interval is descriptive site bootstrap, not a P value."
        ),
    }


def site_counts_by_period(
    data: pd.DataFrame, marker_col: str, categories: list[str]
) -> dict[str, list[tuple[tuple[str, str], int, np.ndarray]]]:
    output: dict[str, list[tuple[tuple[str, str], int, np.ndarray]]] = {
        period: [] for period in BIN_LABELS
    }
    called = data[data[marker_col] != ""].copy()
    for (country, site, period), group in called.groupby(
        ["country", "locality", "analysis_bin"], observed=True
    ):
        count = group[marker_col].value_counts()
        vector = np.array([count.get(category, 0) for category in categories], int)
        output[str(period)].append(((str(country), str(site)), len(group), vector))
    return output


def mean_site_profile(entries: list[tuple[tuple[str, str], int, np.ndarray]]) -> np.ndarray:
    return np.mean([counts / n for _, n, counts in entries], axis=0)


def expected_plugin_tv(
    sizes_a: list[int],
    sizes_b: list[int],
    q_a: np.ndarray,
    q_b: np.ndarray,
    n_draws: int,
    rng: np.random.Generator,
) -> np.ndarray:
    result = np.empty(n_draws)
    for replicate in range(n_draws):
        pa = np.mean([rng.multinomial(n, q_a) / n for n in sizes_a], axis=0)
        pb = np.mean([rng.multinomial(n, q_b) / n for n in sizes_b], axis=0)
        result[replicate] = total_variation(pa, pb)
    return result


def finite_sample_tv_adjustments(
    paired: pd.DataFrame,
    marker: str,
    marker_col: str,
    categories: list[str],
    encoding: str,
) -> pd.DataFrame:
    by_period = site_counts_by_period(paired, marker_col, categories)
    rows: list[dict[str, object]] = []
    seed_encoding = "common_paired_mt" if marker == "mtDNA" else encoding
    for first, second in zip(BIN_LABELS[:-1], BIN_LABELS[1:]):
        a, b = by_period[first], by_period[second]
        observed_a, observed_b = mean_site_profile(a), mean_site_profile(b)
        observed = total_variation(observed_a, observed_b)
        sizes_a, sizes_b = [x[1] for x in a], [x[1] for x in b]
        aggregate_a = np.sum([x[2] for x in a], axis=0).astype(float)
        aggregate_b = np.sum([x[2] for x in b], axis=0).astype(float)
        # Jeffreys smoothing makes every multinomial simulation well-defined.
        q_a = (aggregate_a + 0.5) / (aggregate_a.sum() + 0.5 * len(categories))
        q_b = (aggregate_b + 0.5) / (aggregate_b.sum() + 0.5 * len(categories))
        q_null = (aggregate_a + aggregate_b + 0.5) / (
            aggregate_a.sum() + aggregate_b.sum() + 0.5 * len(categories)
        )
        null_draws = expected_plugin_tv(
            sizes_a,
            sizes_b,
            q_null,
            q_null,
            N_PARAMETRIC,
            named_rng(f"tv-null:{seed_encoding}:{marker}:{first}"),
        )
        alt_draws = expected_plugin_tv(
            sizes_a,
            sizes_b,
            q_a,
            q_b,
            N_PARAMETRIC,
            named_rng(f"tv-alt:{seed_encoding}:{marker}:{first}"),
        )
        true_working_tv = total_variation(q_a, q_b)
        alt_bias = float(alt_draws.mean() - true_working_tv)
        rows.append(
            {
                "y_encoding": encoding,
                "marker": marker,
                "transition": f"{first} -> {second}",
                "first_bin": first,
                "second_bin": second,
                "n_sites_first": len(a),
                "n_sites_second": len(b),
                "n_calls_first": sum(sizes_a),
                "n_calls_second": sum(sizes_b),
                "observed_site_balanced_tv": observed,
                "null_expected_tv_noise_floor": float(null_draws.mean()),
                "null_noise_ci_low": float(np.quantile(null_draws, 0.025)),
                "null_noise_ci_high": float(np.quantile(null_draws, 0.975)),
                "noise_floor_subtracted_tv": max(0.0, observed - float(null_draws.mean())),
                "working_model_true_tv": true_working_tv,
                "working_model_expected_plugin_tv": float(alt_draws.mean()),
                "working_model_estimated_bias": alt_bias,
                "parametric_bias_corrected_tv": max(0.0, observed - alt_bias),
                "simulation_draws": N_PARAMETRIC,
                "assumption": (
                    "Fixed observed site sample sizes; multinomial sampling; Jeffreys-smoothed "
                    "individual-count compositions; no extra site heterogeneity. Noise-floor "
                    "subtraction is a sensitivity, not an unbiased estimator."
                ),
            }
        )
    return pd.DataFrame(rows)


def paired_cluster_bootstrap_adjusted_delta(
    paired: pd.DataFrame,
    y_col: str,
    mt_categories: list[str],
    y_categories: list[str],
    adjustment: pd.DataFrame,
    encoding: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    paired = paired.copy()
    paired["cluster"] = paired["country"].astype(str) + "|||" + paired["locality"].astype(str)
    clusters = paired["cluster"].drop_duplicates().tolist()
    cluster_index = {cluster: i for i, cluster in enumerate(clusters)}
    profiles: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {
        period: {} for period in BIN_LABELS
    }
    for (cluster, period), group in paired.groupby(["cluster", "analysis_bin"], observed=True):
        mt_count = group["mt_l1_pooled"].value_counts()
        y_count = group[y_col].value_counts()
        mt = np.array([mt_count.get(c, 0) for c in mt_categories], float)
        yy = np.array([y_count.get(c, 0) for c in y_categories], float)
        profiles[str(period)][str(cluster)] = (mt / mt.sum(), yy / yy.sum())
    noise = {
        (row.marker, row.transition): row.null_expected_tv_noise_floor
        for row in adjustment.itertuples()
    }
    bias = {
        (row.marker, row.transition): row.working_model_estimated_bias
        for row in adjustment.itertuples()
    }
    rng = named_rng(f"paired-adjusted-cluster:{encoding}")
    p_cluster = np.full(len(clusters), 1 / len(clusters))
    rows: list[dict[str, float | int]] = []
    rejected = 0
    for replicate in range(N_CLUSTER_BOOT):
        while True:
            weight = rng.multinomial(len(clusters), p_cluster)
            mt_period: dict[str, np.ndarray] = {}
            y_period: dict[str, np.ndarray] = {}
            complete = True
            for period in BIN_LABELS:
                present = list(profiles[period])
                w = np.array([weight[cluster_index[c]] for c in present])
                if w.sum() == 0:
                    complete = False
                    break
                mt_period[period] = np.average(
                    np.stack([profiles[period][c][0] for c in present]), axis=0, weights=w
                )
                y_period[period] = np.average(
                    np.stack([profiles[period][c][1] for c in present]), axis=0, weights=w
                )
            if complete:
                break
            rejected += 1
        mt_raw: list[float] = []
        y_raw: list[float] = []
        mt_noise: list[float] = []
        y_noise: list[float] = []
        mt_bias: list[float] = []
        y_bias: list[float] = []
        for first, second in zip(BIN_LABELS[:-1], BIN_LABELS[1:]):
            transition = f"{first} -> {second}"
            mt_value = total_variation(mt_period[first], mt_period[second])
            y_value = total_variation(y_period[first], y_period[second])
            mt_raw.append(mt_value)
            y_raw.append(y_value)
            mt_noise.append(max(0.0, mt_value - noise[("mtDNA", transition)]))
            y_noise.append(max(0.0, y_value - noise[("Y", transition)]))
            mt_bias.append(max(0.0, mt_value - bias[("mtDNA", transition)]))
            y_bias.append(max(0.0, y_value - bias[("Y", transition)]))
        row: dict[str, float | int] = {"replicate": replicate}
        for label, mt_values, y_values in [
            ("raw", mt_raw, y_raw),
            ("noise_floor_subtracted", mt_noise, y_noise),
            ("parametric_bias_corrected", mt_bias, y_bias),
        ]:
            mt_mean, y_mean = float(np.mean(mt_values)), float(np.mean(y_values))
            row[f"mt_{label}_mean_tv"] = mt_mean
            row[f"y_{label}_mean_tv"] = y_mean
            row[f"delta_y_minus_mt_{label}"] = y_mean - mt_mean
        rows.append(row)
    draws = pd.DataFrame(rows)
    summary_rows: list[dict[str, object]] = []
    point_column = {
        "raw": "observed_site_balanced_tv",
        "noise_floor_subtracted": "noise_floor_subtracted_tv",
        "parametric_bias_corrected": "parametric_bias_corrected_tv",
    }
    for label in ["raw", "noise_floor_subtracted", "parametric_bias_corrected"]:
        delta = draws[f"delta_y_minus_mt_{label}"]
        value_col = point_column[label]
        mt_point = float(adjustment.loc[adjustment.marker == "mtDNA", value_col].mean())
        y_point = float(adjustment.loc[adjustment.marker == "Y", value_col].mean())
        theta = y_point - mt_point
        q_low = float(delta.quantile(0.025))
        q_high = float(delta.quantile(0.975))
        bootstrap_median = float(delta.median())
        summary_rows.append(
            {
                "y_encoding": encoding,
                "adjustment": label,
                "mt_point_mean_tv": mt_point,
                "y_point_mean_tv": y_point,
                "delta_point_y_minus_mt": theta,
                "bootstrap_delta_median_diagnostic": bootstrap_median,
                "bootstrap_delta_mean_minus_point": float(delta.mean() - theta),
                "percentile_delta_ci_low": q_low,
                "percentile_delta_ci_high": q_high,
                "basic_delta_ci_low": 2 * theta - q_high,
                "basic_delta_ci_high": 2 * theta - q_low,
                "median_centered_delta_ci_low": theta + q_low - bootstrap_median,
                "median_centered_delta_ci_high": theta + q_high - bootstrap_median,
                "cluster_bootstrap_replicates": N_CLUSTER_BOOT,
                "rejected_empty_period_draws": rejected,
                "total_attempted_cluster_draws": N_CLUSTER_BOOT + rejected,
                "rejected_empty_period_fraction": rejected / (N_CLUSTER_BOOT + rejected),
                "interval_caveat": (
                    "Observed-data values are the point estimates. Percentile, basic and "
                    "median-centered intervals are all reported because sparse-bin ratio "
                    "statistics make the bootstrap distribution off-center. Empty-period "
                    "draws are rejected; the reported fraction quantifies this conditioning. "
                    "Correction terms are fixed, so their model uncertainty is not propagated."
                ),
            }
        )
    return draws, pd.DataFrame(summary_rows)


def site_balanced_composition(
    data: pd.DataFrame,
    marker_col: str,
    categories: list[str],
    group_col: str,
    levels: list[str],
    min_calls: int,
) -> pd.DataFrame:
    profiles = site_profile_table(data, marker_col, categories, min_calls, time_col=group_col)
    rows: list[dict[str, object]] = []
    for level in levels:
        subset = profiles[profiles["time_group"] == level]
        if subset.empty:
            p = np.full(len(categories), np.nan)
        else:
            p = subset[[f"p_{cat}" for cat in categories]].mean().to_numpy(float)
        row: dict[str, object] = {
            "time_group": level,
            "n_profiles": len(subset),
            "n_sites": subset["site"].nunique(),
            "n_calls": int(subset["n_calls"].sum()) if len(subset) else 0,
        }
        row.update({cat: val for cat, val in zip(categories, p)})
        rows.append(row)
    return pd.DataFrame(rows)


def time_grid_sensitivity(
    data: pd.DataFrame, marker_spec: dict[str, tuple[str, list[str]]]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    edges = np.arange(-3500, 1501, 500)
    labels = [f"E{i+1} {edges[i]} to {edges[i+1]-1} CE-axis" for i in range(len(edges) - 1)]
    data = data.copy()
    data["equal_500y_bin"] = pd.cut(
        pd.to_numeric(data["year_ce"]),
        bins=edges,
        labels=labels,
        right=False,
        include_lowest=True,
    ).astype(str).replace("nan", "")
    compositions: list[pd.DataFrame] = []
    turnover_rows: list[dict[str, object]] = []
    original_scaled: list[dict[str, object]] = []
    for marker, (column, categories) in marker_spec.items():
        for minimum in [1, 2]:
            comp = site_balanced_composition(
                data, column, categories, "equal_500y_bin", labels, minimum
            )
            comp.insert(0, "marker", marker)
            comp.insert(1, "minimum_calls", minimum)
            compositions.append(comp)
            for i in range(len(comp) - 1):
                a = comp.loc[i, categories].to_numpy(float)
                b = comp.loc[i + 1, categories].to_numpy(float)
                turnover_rows.append(
                    {
                        "marker": marker,
                        "minimum_calls": minimum,
                        "transition": f"{labels[i]} -> {labels[i+1]}",
                        "interval_years": 500,
                        "tv": total_variation(a, b),
                        "n_profiles_first": comp.loc[i, "n_profiles"],
                        "n_profiles_second": comp.loc[i + 1, "n_profiles"],
                        "interpretation": (
                            "Equal-width-bin dissimilarity; not a lineage-change velocity."
                        ),
                    }
                )
            original = site_balanced_composition(
                data, column, categories, "analysis_bin", BIN_LABELS, minimum
            )
            for i, (first, second) in enumerate(zip(BIN_LABELS[:-1], BIN_LABELS[1:])):
                tv = total_variation(
                    original.loc[i, categories].to_numpy(float),
                    original.loc[i + 1, categories].to_numpy(float),
                )
                gap = BIN_MIDPOINT[second] - BIN_MIDPOINT[first]
                original_scaled.append(
                    {
                        "marker": marker,
                        "minimum_calls": minimum,
                        "transition": f"{first} -> {second}",
                        "midpoint_gap_years": gap,
                        "tv": tv,
                        "tv_per_100_years_linear_scaling": tv / gap * 100,
                        "warning": (
                            "Descriptive linear scaling only; assumes TV accumulates linearly "
                            "and must not be called a biological rate."
                        ),
                    }
                )
    return (
        pd.concat(compositions, ignore_index=True),
        pd.DataFrame(turnover_rows),
        pd.DataFrame(original_scaled),
    )


def continuous_time_tests(
    profiles: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for marker in profiles["marker"].unique():
        for minimum in [1, 2, 3]:
            table = profiles[
                (profiles["marker"] == marker)
                & (profiles["minimum_calls"] == minimum)
            ].copy()
            x = table[["mean_year_ce"]].to_numpy(float)
            x = StandardScaler().fit_transform(x)
            linear = x
            spline = SplineTransformer(
                degree=3, n_knots=4, knots="quantile", include_bias=False
            ).fit_transform(x)
            for label, effect in [("linear year", linear), ("cubic B-spline year", spline)]:
                result = cluster_wild_test(
                    table,
                    effect,
                    N_WILD,
                    named_rng(f"continuous:{marker}:{minimum}:{label}"),
                    label,
                )
                result.update(
                    {
                        "marker": marker,
                        "minimum_calls": minimum,
                        "estimand": (
                            "country-adjusted association between mean site-profile year and "
                            "Hellinger lineage composition; association, not turnover velocity"
                        ),
                    }
                )
                rows.append(result)
    output = pd.DataFrame(rows)
    output["holm_within_minimum_and_time_basis_p"] = np.nan
    for (_, minimum), idx in output.groupby(["effect", "minimum_calls"]).groups.items():
        output.loc[idx, "holm_within_minimum_and_time_basis_p"] = holm(
            output.loc[idx, "raw_cluster_wild_p"]
        )
    return output


def write_readme(
    stable: pd.DataFrame,
    identifiability: pd.DataFrame,
    predictive: pd.DataFrame,
    adjusted_summary: pd.DataFrame,
    equal_tv: pd.DataFrame,
    continuous: pd.DataFrame,
    output_path: Path,
) -> None:
    stable_view = stable[
        [
            "marker",
            "minimum_calls",
            "analysis_role",
            "n_analyzed_profiles",
            "n_calls_in_analyzed_profiles",
            "n_sites",
            "partial_r2",
            "raw_cluster_wild_p",
            "holm_within_threshold_p",
        ]
    ]
    equal_summary = (
        equal_tv.groupby(["marker", "minimum_calls"])
        .agg(mean_tv=("tv", "mean"), max_tv=("tv", "max"))
        .reset_index()
    )
    text = f"""# Statistical extensions v4

This directory was generated from the schema-locked, project-coded AADR-
derived analytical input included with the release. The input SHA-256 is
`{sha256(INPUT)}`.

This run used seed {SEED}, {N_WILD:,} cluster-wild resamples,
{N_PARAMETRIC:,} parametric draws per transition, and
{N_CLUSTER_BOOT:,} paired-cluster bootstrap replicates. The frozen defaults
are seed {DEFAULT_SEED}, {DEFAULT_N_WILD:,}, {DEFAULT_N_PARAMETRIC:,}, and
{DEFAULT_N_CLUSTER_BOOT:,}, respectively. Runs with fewer replicates are for
execution checks and should not replace the frozen analysis for inference.

The input is project-coded rather than anonymous because combinations of
public AADR-derived attributes may remain linkable to the source resource. No
source person identifier, locality string, skeletal field, exact coordinate or
terminal haplogroup call is required by this extension.

## Decision rules

* **Revised primary profile estimand:** the country-adjusted difference among
  equally weighted site-period lineage profiles containing at least two calls.
  This excludes indicator-vector singleton profiles while retaining a usable
  number of localities. It is an archive/site-profile estimand, not a Central
  Asian population-frequency estimand.
* **Stringent sensitivity:** profiles containing at least three calls.
* The all-profile analysis remains an archival-scope context analysis so the
  revision does not erase the originally analysed record.
* Cluster-wild tests use a null-imposed HC2 Rademacher sign shared across all
  time profiles of a country-locality. Holm adjustment is applied across the
  two markers within each threshold. Because the thresholds are nested, the
  >=3 result is interpreted as sensitivity rather than a second independent
  confirmatory test.

## Stable-profile results

{markdown_table(stable_view, digits=6)}

The all-profile Y row is an extension rerun with a different fixed seed. The
frozen 25 July 2026 analysis reported raw/Holm *P*=0.0012; the extension gave
0.0009 under the frozen defaults. The result of the current run is shown above.
The discussion below describes the frozen-default analysis; changing the seed
or resampling counts requires reviewing the current tables before using it.

The mtDNA period association in the all-profile archive is not reproduced
after singleton exclusion. The Y result persists at >=2 calls, while the >=3
analysis has few residual degrees of freedom and should not be overinterpreted.

## Why no claimed hierarchical Dirichlet-multinomial result

The requested site/publication multilevel specification is structurally weak:
most sites occur in only one period, so a full site intercept is almost nested
with period. The identifiability audit is:

{markdown_table(identifiability, digits=4)}

Fitting a nominal mixed model would not manufacture independent longitudinal
information. Instead, `multinomial_predictive_sensitivity.csv` reports a
nested site-group cross-validation analysis using L2-shrunk country,
publication and site intercept terms, with and without period. A ridge penalty
is a Gaussian-prior/MAP analogue, but this remains a predictive sensitivity,
not a Dirichlet-multinomial hypothesis test.

{markdown_table(predictive, digits=5)}

Negative full-minus-reduced log loss means period improved prediction for
held-out sites. The site-bootstrap interval is descriptive and is not a P
value.

## Finite-sample TV sensitivity

For each adjacent transition and marker, two corrections are supplied:

1. `null_expected_tv_noise_floor`: expected positive TV when the two periods
   share one Jeffreys-smoothed composition but retain their observed per-site
   sample sizes.
2. `working_model_estimated_bias`: parametric plug-in bias when the periods
   have their own Jeffreys-smoothed compositions.

Neither correction is assumption-free. In particular, multinomial simulation
does not include extra site heterogeneity, and subtracting a null noise floor
does not generally produce an unbiased distance. Therefore corrected values
are sensitivities, not replacements for the raw estimator.

{markdown_table(adjusted_summary, digits=5)}

The previously displayed surrogate ordinary-bootstrap probabilities are
intentionally not reproduced because they are not null-hypothesis P values.
Interpretation should be based on effect estimates, intervals and sensitivity
across Y encodings.

## Time-grid sensitivity

The manuscript bins span unequal durations. `equal_500y_*` recomputes
site-balanced composition and adjacent TV in ten equal 500-y bins. Summary:

{markdown_table(equal_summary, digits=5)}

`original_bin_linear_scaled_tv.csv` also exposes TV divided by midpoint gap,
but labels it explicitly as a linear-scaling diagnostic rather than a
biological rate. `continuous_time_cluster_tests.csv` tests a country-adjusted
linear and flexible spline association between profile year and composition:

{markdown_table(continuous[['marker','minimum_calls','effect','n_analyzed_profiles','partial_r2','raw_cluster_wild_p','holm_within_minimum_and_time_basis_p','valid']], digits=5)}

These checks support only statements about archive-indexed compositional
dissimilarity/association. They do not identify a constant lineage turnover
velocity, migration rate, or sex-biased demographic mechanism.

## Files

* `stable_profile_cluster_tests.csv` and `stable_profile_vectors.csv`
* `hierarchical_model_identifiability_audit.csv`
* `multinomial_predictive_sensitivity.csv`
* `paired_male_tv_finite_sample_adjustments.csv`
* `paired_male_tv_adjusted_cluster_bootstrap_summary.csv`
* `paired_male_tv_adjusted_cluster_bootstrap_draws_*.csv`
* `equal_500y_composition.csv`, `equal_500y_turnover.csv`
* `original_bin_linear_scaled_tv.csv`
* `continuous_time_cluster_tests.csv`
* `run_manifest.json`
"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    global INPUT, OUT, SEED, N_WILD, N_PARAMETRIC, N_CLUSTER_BOOT
    args = parse_args(argv)
    INPUT, OUT, SEED = args.input, args.outdir, args.seed
    N_WILD = args.wild_resamples
    N_PARAMETRIC = args.parametric_draws
    N_CLUSTER_BOOT = args.paired_bootstrap
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT, keep_default_na=False)
    expected_public = {
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
    if set(data.columns) != expected_public:
        raise ValueError("Project-coded analytical input schema mismatch")
    data = data.rename(
        columns={
            "site_key": "locality",
            "study_key": "publication",
            "mt_category": "mt_l1_pooled",
            "y_category": "y_l1_pooled",
            "y_prefix_category": "y_isogg_prefix_family_pooled",
        }
    )
    data["year_ce"] = 1950 - pd.to_numeric(data["date_bp"], errors="raise")
    for column in ["mt_l1_pooled", "y_l1_pooled"]:
        data[column] = data[column].astype(str)
    y_family_categories = (
        data.loc[
            data["y_isogg_prefix_family_pooled"] != "",
            "y_isogg_prefix_family_pooled",
        ]
        .value_counts()
        .index.tolist()
    )
    marker_spec = {
        "mtDNA": ("mt_l1_pooled", data.loc[data.mt_l1_pooled != "", "mt_l1_pooled"].value_counts().index.tolist()),
        "Y": ("y_l1_pooled", data.loc[data.y_l1_pooled != "", "y_l1_pooled"].value_counts().index.tolist()),
    }

    stable, profiles = stable_profile_tests(data, marker_spec)
    stable.to_csv(OUT / "stable_profile_cluster_tests.csv", index=False)
    profiles.to_csv(OUT / "stable_profile_vectors.csv", index=False)

    identifiability = pd.DataFrame(
        [
            identifiability_audit(profiles, marker, minimum)
            for marker in marker_spec
            for minimum in [1, 2, 3]
        ]
    )
    identifiability.to_csv(OUT / "hierarchical_model_identifiability_audit.csv", index=False)

    predictive = pd.DataFrame(
        [
            weighted_multinomial_predictive_sensitivity(data, marker, column, minimum)
            for marker, (column, _) in marker_spec.items()
            for minimum in [2, 3]
        ]
    )
    predictive.to_csv(OUT / "multinomial_predictive_sensitivity.csv", index=False)

    paired = data[
        data["molecular_sex"].eq("M")
        & (data["mt_l1_pooled"] != "")
        & (data["y_l1_pooled"] != "")
        & (data["y_isogg_prefix_family_pooled"] != "")
    ].copy()
    mt_categories = marker_spec["mtDNA"][1]
    adjustment_frames: list[pd.DataFrame] = []
    bootstrap_summaries: list[pd.DataFrame] = []
    for encoding, y_col, y_categories in [
        ("broad_L1", "y_l1_pooled", marker_spec["Y"][1]),
        ("AADR_ISOGG_prefix_family", "y_isogg_prefix_family_pooled", y_family_categories),
    ]:
        mt_adjustment = finite_sample_tv_adjustments(
            paired, "mtDNA", "mt_l1_pooled", mt_categories, encoding
        )
        y_adjustment = finite_sample_tv_adjustments(
            paired, "Y", y_col, y_categories, encoding
        )
        adjustment = pd.concat([mt_adjustment, y_adjustment], ignore_index=True)
        adjustment_frames.append(adjustment)
        draws, summary = paired_cluster_bootstrap_adjusted_delta(
            paired,
            y_col,
            mt_categories,
            y_categories,
            adjustment,
            encoding,
        )
        safe = "broad" if encoding == "broad_L1" else "isogg_prefix"
        draws.to_csv(
            OUT / f"paired_male_tv_adjusted_cluster_bootstrap_draws_{safe}.csv",
            index=False,
        )
        bootstrap_summaries.append(summary)
    adjustments = pd.concat(adjustment_frames, ignore_index=True)
    adjusted_summary = pd.concat(bootstrap_summaries, ignore_index=True)
    adjustments.to_csv(OUT / "paired_male_tv_finite_sample_adjustments.csv", index=False)
    adjusted_summary.to_csv(
        OUT / "paired_male_tv_adjusted_cluster_bootstrap_summary.csv", index=False
    )

    equal_comp, equal_tv, original_scaled = time_grid_sensitivity(data, marker_spec)
    equal_comp.to_csv(OUT / "equal_500y_composition.csv", index=False)
    equal_tv.to_csv(OUT / "equal_500y_turnover.csv", index=False)
    original_scaled.to_csv(OUT / "original_bin_linear_scaled_tv.csv", index=False)
    continuous = continuous_time_tests(profiles)
    continuous.to_csv(OUT / "continuous_time_cluster_tests.csv", index=False)

    readme_path = args.readme_output
    write_readme(stable, identifiability, predictive, adjusted_summary, equal_tv,
                 continuous, readme_path)
    manifest = {
        "input": manifest_path(INPUT),
        "input_sha256": sha256(INPUT),
        "analysis_script": manifest_path(Path(__file__)),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "readme": manifest_path(readme_path),
        "readme_sha256": sha256(readme_path),
        "test_script": manifest_path(Path(__file__).resolve().parent / "test_results.py"),
        "test_script_sha256": sha256(Path(__file__).resolve().parent / "test_results.py"),
        "results_dir": manifest_path(OUT),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            **{package: importlib.metadata.version(package)
               for package in ("numpy", "pandas", "scipy", "scikit-learn")},
        },
        "seed": SEED,
        "wild_resamples": N_WILD,
        "parametric_tv_draws_per_transition": N_PARAMETRIC,
        "paired_cluster_bootstrap_replicates": N_CLUSTER_BOOT,
        "predictive_site_bootstrap_replicates": 20_000,
        "n_catalogue_rows": len(data),
        "n_paired_males": len(paired),
        "outputs": {},
    }
    for name in sorted(OUTPUT_NAMES):
        manifest["outputs"][name] = sha256(OUT / name)
    (OUT / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
