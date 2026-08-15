#!/usr/bin/env python3
"""Reproducible uniparental-lineage synthesis for ancient Central Asia.

The script deliberately estimates the composition of published archaeological
samples and sites. It does not treat those samples as a probability sample of
past populations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency, t as student_t


COUNTRIES = [
    "Kazakhstan",
    "Kyrgyzstan",
    "Tajikistan",
    "Turkmenistan",
    "Uzbekistan",
]

BIN_EDGES = np.array([-3500, -2500, -1800, -900, -200, 301, 651, 1001, 1501])
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

CONTEXT_EDGES = np.array(
    [-100000, -4500, -3500, -2500, -1800, -900, -200, 301, 651, 1001, 1501, 100000]
)
CONTEXT_LABELS = [
    "Earlier than 4500 BCE",
    "B0 4500-3501 BCE",
    *BIN_LABELS,
    "After 1500 CE",
]

MISSING = {"", ".", "..", "na", "n/a", "nan", "none", "unknown", "?"}

COUNTRY_COLORS = {
    "Kazakhstan": "#2C7FB8",
    "Kyrgyzstan": "#41AB5D",
    "Tajikistan": "#F03B20",
    "Turkmenistan": "#756BB1",
    "Uzbekistan": "#FDAE6B",
}

EXPECTED_INPUT_SHA256 = {
    "aadr": "98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c",
    "amtdb": "531e8ee8fae181124f5a9b77b6fe8d677e64e35b815be2a3965020244fe31057",
    "aychr": "e297110a18cba73d4044e8a95c0fae98d7f48633ad6ccfae1cd364e460eb1b3c",
}

# Manually verified aliases that cannot be resolved by literal identifier
# matching. Each link was checked against site, date, publication and marker
# call; these rows remain explicitly labelled "manual_verified".
MANUAL_CROSSWALK = {
    ("AmtDB v1.009", "DA223"): "I26001",
    ("AmtDB v1.009", "DA224"): "I10141",
    ("aYChr-DB v5", "TU45/ BOT14"): "BOT14",
    ("aYChr-DB v5", "I0563/Be11"): "I0563",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--aadr", required=True, type=Path)
    p.add_argument("--amtdb", required=True, type=Path)
    p.add_argument("--aychr", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--paired-bootstrap", type=int, default=50000)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--callability-resamples", type=int, default=99999)
    p.add_argument("--date-draws", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260725)
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def named_rng(seed: int, name: str) -> np.random.Generator:
    """Create a stable RNG stream that is independent of call order."""
    digest = hashlib.sha256(name.encode("utf-8")).digest()
    entropy = [seed]
    entropy.extend(
        int.from_bytes(digest[index : index + 4], "little")
        for index in range(0, 16, 4)
    )
    return np.random.default_rng(np.random.SeedSequence(entropy))


def validate_frozen_inputs(args: argparse.Namespace) -> dict[str, str]:
    observed = {}
    for name, expected in EXPECTED_INPUT_SHA256.items():
        path = getattr(args, name)
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen {name} input: {path}")
        digest = sha256(path)
        observed[name] = digest
        if digest != expected:
            raise ValueError(
                f"{name} SHA-256 mismatch: expected {expected}, observed {digest}"
            )
    return observed


def mask_coordinates(data: pd.DataFrame, digits: int = 1) -> pd.DataFrame:
    """Round public archaeological coordinates to documented coarse precision."""
    result = data.copy()
    for column in ("latitude", "longitude"):
        if column in result:
            result[column] = pd.to_numeric(
                result[column], errors="coerce"
            ).round(digits)
    return result


def valid_call(value: object) -> bool:
    # Marker calls are scalar catalogue fields.  Check scalar-ness before
    # ``pd.isna`` so an accidentally supplied Series/array cannot trigger the
    # ambiguous-truth-value error produced by ``if pd.isna(array)``.
    if not pd.api.types.is_scalar(value):
        return False
    if pd.isna(value):
        return False
    s = str(value).strip()
    low = s.lower()
    return low not in MISSING and not low.startswith("n/a")


def normalize_id(value: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def find_col(df: pd.DataFrame, prefix: str) -> str:
    cols = [c for c in df.columns if c.startswith(prefix)]
    if len(cols) != 1:
        raise ValueError(f"Expected one column beginning {prefix!r}; found {cols}")
    return cols[0]


def major_haplogroup(call: object, marker: str) -> str:
    """Create a deliberately broad, non-geographic lineage label.

    The mapping is deliberately explicit for compound roots. In particular,
    mitochondrial HV is not a subclade of H, and basal compound Y calls such
    as CF or IJK must not be assigned to their first letter.
    """
    if not valid_call(call):
        return ""
    s = str(call).strip().upper()
    s = s.replace("HAPLOGROUP", "").strip()
    if marker == "mt":
        if s.startswith("HV"):
            return "HV"
    if marker == "Y":
        basal = s.rstrip("*")
        if basal in {
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
    m = re.search(r"[A-Z]", s)
    if not m:
        return "Other/unresolved"
    return m.group(0)


def y_isogg_family_prefix(call: object) -> str:
    """Encode an AADR ISOGG-style Y call at a transparent prefix depth.

    The encoding retains the first letter, first integer and the immediately
    following branch letter (for example R1a1a1 -> R1a, J2a1a4b -> J2a).
    Single-letter or letter-number calls are retained, while basal compound
    roots use the same ``Basal/unresolved`` label as the primary L1 encoding.
    This is a nomenclature-prefix sensitivity, not a phylogenetic re-call; its
    depth is not guaranteed to be biologically uniform across haplogroups.
    """
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


def add_y_resolution_sensitivity_encoding(
    df: pd.DataFrame,
    minimum_category_count: int = 5,
) -> tuple[pd.DataFrame, list[str], str]:
    """Add a pooled ISOGG-prefix encoding without altering the primary Y call."""
    result = df.copy()
    source_column = (
        "y_isogg_call" if "y_isogg_call" in result.columns else "y_isogg"
    )
    encoded = result[source_column].map(y_isogg_family_prefix)
    # A family label is only eligible where the marker-specific Y call exists.
    encoded = encoded.where(result["y_call"].map(valid_call), "")
    result["y_isogg_prefix_family"] = encoded
    pooled, categories = pool_rare(encoded, minimum_category_count)
    result["y_isogg_prefix_family_pooled"] = pooled
    return result, categories, source_column


class UnionFind:
    def __init__(self, items: list[str]):
        self.parent = {x: x for x in items}

    def find(self, x: str) -> str:
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def derive_kin_components(df: pd.DataFrame, max_degree: float = 2.0) -> pd.Series:
    ids = df["individual_id"].astype(str).tolist()
    known = set(ids)
    uf = UnionFind(ids)
    for rel in df["family_relations"].fillna("").astype(str):
        if not valid_call(rel):
            continue
        for item in rel.split(","):
            parts = [p.strip() for p in item.split(":")]
            if len(parts) < 3:
                continue
            degree_match = re.match(r"([0-9.]+)d", parts[0])
            if not degree_match or float(degree_match.group(1)) > max_degree:
                continue
            a, b = parts[-2], parts[-1]
            if a in known and b in known:
                uf.union(a, b)
    roots = [uf.find(x) for x in ids]
    members: dict[str, list[str]] = defaultdict(list)
    for iid, root in zip(ids, roots):
        members[root].append(iid)
    label = {
        iid: ("KIN_" + "__".join(sorted(group)) if len(group) > 1 else iid)
        for group in members.values()
        for iid in group
    }
    return df["individual_id"].map(label)


def choose_canonical_aadr(aadr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    c = {
        "genetic_id": aadr.columns[0],
        "persistent_genetic_id": "Persistent Genetic ID",
        "individual_id": "Individual ID",
        "skeletal_code": "Skeletal code",
        "skeletal_element": "Skeletal element",
        "first_publication": find_col(aadr, "First publication:"),
        "publication": "Publication abbreviation",
        "doi": "doi for publication of this representation of the data",
        "repository": "Link to the most permanent repository hosting these data",
        "date_method": find_col(aadr, "Method for Determining Date"),
        "date_bp": find_col(aadr, "Date mean in BP"),
        "date_sd_bp": find_col(aadr, "Date standard deviation in BP"),
        "full_date": find_col(aadr, "Full Date"),
        "age_sex_osteo": "Age at death, Morphological sex from physical anthropology",
        "group_id": "Group ID",
        "locality": "Locality",
        "country": "Political Entity",
        "latitude": "Latitude",
        "longitude": "Longitude",
        "pulldown": "Pulldown Strategy",
        "data_suffix": "Suffices (indicating data types used for sources which can be a subset of that in bam)",
        "data_type": "Data type",
        "coverage": find_col(aadr, "Mean coverage on 1.15M"),
        "snp_2m": find_col(aadr, "SNPs hit on autosomal targets (Computed using easystats on enhance"),
        "molecular_sex": "Molecular Sex",
        "family_relations": "Family relations",
        "y_terminal": find_col(aadr, "Y haplogroup in terminal"),
        "y_isogg": find_col(aadr, "Y haplogroup  in ISOGG"),
        "y_manual": find_col(aadr, "Y haplogroup manually"),
        "mt_coverage": "mtDNA coverage (merged data)",
        "mt_haplogroup": "mtDNA haplogroup if >2x or published",
        "library_type": find_col(aadr, "Library type"),
        "assessment": "ASSESSMENT",
        "assessment_warnings": find_col(aadr, "ASSESSMENT WARNINGS"),
    }
    raw = aadr.rename(columns={v: k for k, v in c.items()})[list(c)].copy()
    raw = raw[raw["country"].isin(COUNTRIES)].copy()
    raw["date_bp"] = pd.to_numeric(raw["date_bp"], errors="coerce")
    raw["date_sd_bp"] = pd.to_numeric(raw["date_sd_bp"], errors="coerce")
    raw["snp_2m"] = pd.to_numeric(raw["snp_2m"], errors="coerce")
    raw["coverage"] = pd.to_numeric(raw["coverage"], errors="coerce")
    raw["mt_coverage"] = pd.to_numeric(raw["mt_coverage"], errors="coerce")
    raw["latitude"] = pd.to_numeric(raw["latitude"], errors="coerce")
    raw["longitude"] = pd.to_numeric(raw["longitude"], errors="coerce")
    raw = raw[raw["date_bp"] > 0].copy()
    raw["y_call_row"] = np.where(
        raw["y_manual"].map(valid_call),
        raw["y_manual"],
        np.where(raw["y_terminal"].map(valid_call), raw["y_terminal"], ""),
    )
    raw["mt_call_row"] = raw["mt_haplogroup"].where(
        raw["mt_haplogroup"].map(valid_call), ""
    )
    raw["representation_count"] = raw.groupby("individual_id")["individual_id"].transform(
        "size"
    )

    # AADR explicitly recommends using the representation with the greatest
    # number of hit autosomal targets when no representation is otherwise
    # privileged.
    metadata = (
        raw.assign(_snp=raw["snp_2m"].fillna(-1))
        .sort_values(["individual_id", "_snp"], ascending=[True, False])
        .drop_duplicates("individual_id")
        .drop(columns="_snp")
    )

    def best_marker(group: pd.DataFrame, marker: str) -> pd.Series:
        if marker == "mt":
            called = group[group["mt_call_row"].map(valid_call)].copy()
            if called.empty:
                return pd.Series(
                    {
                        "mt_call": "",
                        "mt_call_genetic_id": "",
                        "mt_call_coverage": np.nan,
                    }
                )
            called["_quality"] = called["mt_coverage"].fillna(-1)
            row = called.sort_values(
                ["_quality", "snp_2m"], ascending=[False, False]
            ).iloc[0]
            return pd.Series(
                {
                    "mt_call": row["mt_call_row"],
                    "mt_call_genetic_id": row["genetic_id"],
                    "mt_call_coverage": row["mt_coverage"],
                }
            )
        called = group[
            group["y_call_row"].map(valid_call)
            & group["molecular_sex"].astype(str).str.startswith("M")
        ].copy()
        if called.empty:
            return pd.Series({"y_call": "", "y_call_genetic_id": ""})
        called["_quality"] = called["snp_2m"].fillna(-1)
        row = called.sort_values("_quality", ascending=False).iloc[0]
        return pd.Series(
            {"y_call": row["y_call_row"], "y_call_genetic_id": row["genetic_id"]}
        )

    mt = raw.groupby("individual_id", sort=False).apply(
        lambda g: best_marker(g, "mt"), include_groups=False
    )
    y = raw.groupby("individual_id", sort=False).apply(
        lambda g: best_marker(g, "Y"), include_groups=False
    )
    canonical = metadata.drop(
        columns=["mt_call_row", "y_call_row"], errors="ignore"
    ).merge(mt, on="individual_id", how="left").merge(
        y, on="individual_id", how="left"
    )
    canonical["mt_call"] = canonical["mt_call"].fillna("")
    canonical["y_call"] = canonical["y_call"].fillna("")
    canonical["year_ce"] = 1950 - canonical["date_bp"]
    canonical["context_period"] = pd.cut(
        canonical["year_ce"],
        bins=CONTEXT_EDGES,
        labels=CONTEXT_LABELS,
        right=False,
    ).astype(str)
    canonical["analysis_bin"] = pd.cut(
        canonical["year_ce"],
        bins=BIN_EDGES,
        labels=BIN_LABELS,
        right=False,
    )
    canonical["analysis_included"] = canonical["analysis_bin"].notna()
    canonical["mt_l1"] = canonical["mt_call"].map(lambda x: major_haplogroup(x, "mt"))
    canonical["y_l1"] = canonical["y_call"].map(lambda x: major_haplogroup(x, "Y"))
    canonical["direct_date"] = canonical["date_method"].astype(str).str.startswith(
        "Direct", na=False
    )
    canonical["strict_qc"] = canonical["assessment"].isin(
        ["Pass", "PROVISIONAL_PASS"]
    )
    canonical["population_outlier"] = canonical["group_id"].astype(str).str.contains(
        r"(?:^|[-_])o(?:$|[-_])|outlier", case=False, regex=True
    )
    canonical["kin_component_2d"] = derive_kin_components(canonical, max_degree=2.0)
    canonical["kin_component_2_5d"] = derive_kin_components(
        canonical, max_degree=2.5
    )

    # Conflicts are retained, not silently overwritten.
    audit_rows = []
    for iid, g in raw.groupby("individual_id", sort=False):
        mt_calls = sorted(set(g.loc[g["mt_call_row"].map(valid_call), "mt_call_row"]))
        y_calls = sorted(set(g.loc[g["y_call_row"].map(valid_call), "y_call_row"]))
        audit_rows.append(
            {
                "individual_id": iid,
                "representation_count": len(g),
                "genetic_ids": ";".join(g["genetic_id"].astype(str)),
                "mt_calls_all": ";".join(mt_calls),
                "mt_call_discordant": len(mt_calls) > 1,
                "y_calls_all": ";".join(y_calls),
                "y_call_discordant": len(y_calls) > 1,
            }
        )
    audit = pd.DataFrame(audit_rows)
    return canonical, audit


def pool_rare(series: pd.Series, minimum: int = 5) -> tuple[pd.Series, list[str]]:
    counts = series[series != ""].value_counts()
    keep = counts[counts >= minimum].index.tolist()
    pooled = series.where(series.isin(keep), np.where(series.eq(""), "", "Other"))
    ordered = (
        pooled[pooled != ""].value_counts().sort_values(ascending=False).index.tolist()
    )
    return pooled, ordered


def profile(
    df: pd.DataFrame,
    marker_col: str,
    categories: list[str],
    site_balanced: bool,
) -> pd.DataFrame:
    use = df[df[marker_col] != ""].copy()
    rows = []
    for period in BIN_LABELS:
        g = use[use["analysis_bin"].astype(str) == period]
        if g.empty:
            p = np.full(len(categories), np.nan)
            n_sites = 0
        elif site_balanced:
            site_vectors = []
            for _, sg in g.groupby(["country", "locality"], dropna=False):
                counts = sg[marker_col].value_counts()
                vec = np.array([counts.get(k, 0) for k in categories], dtype=float)
                site_vectors.append(vec / vec.sum())
            p = np.mean(site_vectors, axis=0)
            n_sites = len(site_vectors)
        else:
            counts = g[marker_col].value_counts()
            p = np.array([counts.get(k, 0) for k in categories], dtype=float)
            p /= p.sum()
            n_sites = g["locality"].nunique()
        row = {
            "analysis_bin": period,
            "n_calls": len(g),
            "n_sites": n_sites,
            **{k: v for k, v in zip(categories, p)},
        }
        rows.append(row)
    return pd.DataFrame(rows)


def hill_numbers(p: np.ndarray) -> tuple[float, float]:
    p = p[np.isfinite(p) & (p > 0)]
    if p.size == 0:
        return np.nan, np.nan
    q1 = float(np.exp(-np.sum(p * np.log(p))))
    q2 = float(1.0 / np.sum(p**2))
    return q1, q2


def total_variation(a: np.ndarray, b: np.ndarray) -> float:
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return np.nan
    return float(0.5 * np.abs(a - b).sum())


def bootstrap_site_profiles(
    df: pd.DataFrame,
    marker_col: str,
    categories: list[str],
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Bootstrap country-locality clusters across their complete trajectories.

    One multinomial multiplicity is drawn for every country-locality cluster
    and reused in every period contributed by that cluster. This preserves
    dependence for localities represented in more than one chronological bin.
    Draws that omit every contributing locality from any observed bin are
    rejected so that every accepted replicate retains the eight-bin,
    seven-transition estimand.
    """
    called = df[df[marker_col] != ""].copy()
    called["_cluster"] = list(zip(called["country"], called["locality"]))
    clusters = list(dict.fromkeys(called["_cluster"]))
    cluster_index = {cluster: index for index, cluster in enumerate(clusters)}
    profiles: dict[str, dict[tuple[str, str], np.ndarray]] = {
        period: {} for period in BIN_LABELS
    }
    for (cluster, period), group in called.groupby(
        ["_cluster", "analysis_bin"], sort=False, observed=True
    ):
        counts = group[marker_col].value_counts()
        vector = np.array(
            [counts.get(category, 0) for category in categories], dtype=float
        )
        profiles[str(period)][cluster] = vector / vector.sum()

    div_draws: list[dict[str, object]] = []
    tv_draws: list[dict[str, object]] = []
    rejected_empty_period_draws = 0
    probabilities = np.full(len(clusters), 1.0 / len(clusters))
    for replicate in range(n_boot):
        while True:
            weights = rng.multinomial(len(clusters), probabilities)
            p_by_period: dict[str, np.ndarray] = {}
            complete = True
            for period in BIN_LABELS:
                period_clusters = list(profiles[period])
                period_weights = np.array(
                    [weights[cluster_index[cluster]] for cluster in period_clusters]
                )
                if period_weights.sum() == 0:
                    complete = False
                    break
                vectors = np.stack(
                    [profiles[period][cluster] for cluster in period_clusters]
                )
                p_by_period[period] = np.average(
                    vectors, axis=0, weights=period_weights
                )
            if complete:
                break
            rejected_empty_period_draws += 1

        for period in BIN_LABELS:
            q1, q2 = hill_numbers(p_by_period[period])
            div_draws.append(
                {
                    "replicate": replicate,
                    "analysis_bin": period,
                    "q1": q1,
                    "q2": q2,
                }
            )
        for first, second in zip(BIN_LABELS[:-1], BIN_LABELS[1:]):
            tv_draws.append(
                {
                    "replicate": replicate,
                    "transition": f"{first} -> {second}",
                    "tv": total_variation(
                        p_by_period[first], p_by_period[second]
                    ),
                }
            )
    diagnostics = {
        "n_clusters": len(clusters),
        "accepted_replicates": n_boot,
        "rejected_empty_period_draws": rejected_empty_period_draws,
    }
    return pd.DataFrame(div_draws), pd.DataFrame(tv_draws), diagnostics


def paired_marker_turnover_bootstrap(
    df: pd.DataFrame,
    mt_categories: list[str],
    y_categories: list[str],
    n_boot: int,
    rng: np.random.Generator,
    y_marker_col: str = "y_l1_pooled",
    y_encoding_label: str = "broad_L1",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Compare markers with one shared cluster draw across periods and markers."""
    paired = df[
        df["molecular_sex"].astype(str).str.startswith("M")
        & (df["mt_l1_pooled"] != "")
        & (df[y_marker_col] != "")
    ].copy()
    paired["_cluster"] = list(zip(paired["country"], paired["locality"]))
    clusters = list(dict.fromkeys(paired["_cluster"]))
    cluster_index = {cluster: index for index, cluster in enumerate(clusters)}
    profiles: dict[
        str, dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]
    ] = {period: {} for period in BIN_LABELS}
    for (cluster, period), group in paired.groupby(
        ["_cluster", "analysis_bin"], sort=False, observed=True
    ):
        mt_counts = group["mt_l1_pooled"].value_counts()
        y_counts = group[y_marker_col].value_counts()
        mt = np.array(
            [mt_counts.get(category, 0) for category in mt_categories],
            dtype=float,
        )
        y = np.array(
            [y_counts.get(category, 0) for category in y_categories],
            dtype=float,
        )
        profiles[str(period)][cluster] = (mt / mt.sum(), y / y.sum())

    draws: list[dict[str, float | int]] = []
    rejected_empty_period_draws = 0
    probabilities = np.full(len(clusters), 1.0 / len(clusters))
    for replicate in range(n_boot):
        while True:
            weights = rng.multinomial(len(clusters), probabilities)
            mt_by_period: dict[str, np.ndarray] = {}
            y_by_period: dict[str, np.ndarray] = {}
            complete = True
            for period in BIN_LABELS:
                period_clusters = list(profiles[period])
                period_weights = np.array(
                    [weights[cluster_index[cluster]] for cluster in period_clusters]
                )
                if period_weights.sum() == 0:
                    complete = False
                    break
                mt_vectors = np.stack(
                    [profiles[period][cluster][0] for cluster in period_clusters]
                )
                y_vectors = np.stack(
                    [profiles[period][cluster][1] for cluster in period_clusters]
                )
                mt_by_period[period] = np.average(
                    mt_vectors, axis=0, weights=period_weights
                )
                y_by_period[period] = np.average(
                    y_vectors, axis=0, weights=period_weights
                )
            if complete:
                break
            rejected_empty_period_draws += 1

        mt_tv = [
            total_variation(mt_by_period[first], mt_by_period[second])
            for first, second in zip(BIN_LABELS[:-1], BIN_LABELS[1:])
        ]
        y_tv = [
            total_variation(y_by_period[first], y_by_period[second])
            for first, second in zip(BIN_LABELS[:-1], BIN_LABELS[1:])
        ]
        mt_mean = float(np.mean(mt_tv))
        y_mean = float(np.mean(y_tv))
        draws.append(
            {
                "replicate": replicate,
                "mt_mean_adjacent_tv": mt_mean,
                "y_mean_adjacent_tv": y_mean,
                "delta_y_minus_mt": y_mean - mt_mean,
            }
        )
    draws_df = pd.DataFrame(draws)
    mt_point = mean_adjacent_tv(paired, "mt_l1_pooled", mt_categories)
    y_point = mean_adjacent_tv(paired, y_marker_col, y_categories)
    summary = pd.DataFrame(
        [
            {
                "subset": (
                    "Molecular males with both calls; shared country-locality "
                    "cluster bootstrap across periods and markers"
                ),
                "y_encoding": y_encoding_label,
                "y_marker_column": y_marker_col,
                "n_individuals": len(paired),
                "n_sites": len(clusters),
                "bootstrap_replicates": n_boot,
                "mt_mean_adjacent_tv": mt_point,
                "y_mean_adjacent_tv": y_point,
                "delta_y_minus_mt": y_point - mt_point,
                "bootstrap_median_delta": draws_df[
                    "delta_y_minus_mt"
                ].median(),
                "delta_ci_low": draws_df["delta_y_minus_mt"].quantile(0.025),
                "delta_ci_high": draws_df["delta_y_minus_mt"].quantile(0.975),
                "bootstrap_two_sided_sign_tail_probability": min(
                    1.0,
                    2
                    * min(
                        float((draws_df["delta_y_minus_mt"] <= 0).mean()),
                        float((draws_df["delta_y_minus_mt"] >= 0).mean()),
                    ),
                ),
                "bootstrap_sign_tail_interpretation": (
                    "Bootstrap sign-stability diagnostic; not a null-hypothesis "
                    "p-value"
                ),
            }
        ]
    )
    diagnostics = {
        "n_clusters": len(clusters),
        "accepted_replicates": n_boot,
        "rejected_empty_period_draws": rejected_empty_period_draws,
    }
    return draws_df, summary, diagnostics


def observed_profile_statistics(
    site_profile: pd.DataFrame, categories: list[str]
) -> tuple[dict[str, float], dict[str, float]]:
    """Return the statistics evaluated on the original site-balanced profiles."""
    indexed = site_profile.set_index("analysis_bin").reindex(BIN_LABELS)
    diversity = {
        period: hill_numbers(indexed.loc[period, categories].to_numpy(float))[0]
        for period in BIN_LABELS
    }
    turnover = {
        f"{first} -> {second}": total_variation(
            indexed.loc[first, categories].to_numpy(float),
            indexed.loc[second, categories].to_numpy(float),
        )
        for first, second in zip(BIN_LABELS[:-1], BIN_LABELS[1:])
    }
    return diversity, turnover


def summarize_bootstrap(
    draws: pd.DataFrame,
    value: str,
    group: str,
    observed: dict[str, float],
) -> pd.DataFrame:
    """Attach percentile intervals to an observed point estimate.

    The median of a nonparametric cluster-bootstrap distribution need not equal
    the statistic evaluated on the original sample. It is retained as a
    diagnostic rather than being mislabeled as the point estimate.
    """
    summary = (
        draws.groupby(group)[value]
        .agg(
            bootstrap_median="median",
            ci_low=lambda x: x.quantile(0.025),
            ci_high=lambda x: x.quantile(0.975),
        )
        .reset_index()
    )
    missing = set(summary[group].astype(str)) - set(observed)
    if missing:
        raise ValueError(f"Missing observed estimates for {group}: {sorted(missing)}")
    summary.insert(
        1,
        "estimate",
        summary[group].astype(str).map(observed).astype(float),
    )
    return summary


def site_profile_table(
    df: pd.DataFrame, marker_col: str, categories: list[str], min_calls: int = 1
) -> pd.DataFrame:
    rows = []
    use = df[(df["analysis_bin"].notna()) & (df[marker_col] != "")]
    for (country, site, period), g in use.groupby(
        ["country", "locality", "analysis_bin"], observed=True, dropna=False
    ):
        if len(g) < min_calls:
            continue
        counts = g[marker_col].value_counts()
        p = np.array([counts.get(k, 0) for k in categories], dtype=float)
        p /= p.sum()
        h = np.sqrt(p)
        rows.append(
            {
                "country": country,
                "site": site,
                "analysis_bin": str(period),
                "n_calls": len(g),
                **{f"h_{k}": v for k, v in zip(categories, h)},
            }
        )
    return pd.DataFrame(rows)


def design_matrix(
    table: pd.DataFrame, include_period: bool
) -> tuple[np.ndarray, list[str]]:
    country = pd.get_dummies(table["country"], prefix="country", drop_first=True)
    frames = [pd.Series(1.0, index=table.index, name="intercept"), country]
    if include_period:
        period = pd.get_dummies(
            pd.Categorical(table["analysis_bin"], categories=BIN_LABELS),
            prefix="period",
            drop_first=True,
        )
        period.index = table.index
        frames.append(period)
    x = pd.concat(frames, axis=1).astype(float)
    return x.to_numpy(), x.columns.tolist()


def _sse_from_residual_maker(
    residual_maker: np.ndarray, response: np.ndarray
) -> float:
    residual = residual_maker @ response
    return float(np.sum(residual**2))


def _residual_maker(design: np.ndarray) -> np.ndarray:
    return np.eye(len(design)) - design @ np.linalg.pinv(design)


def site_cluster_wild_period_test(
    site_table: pd.DataFrame,
    n_resamples: int,
    rng: np.random.Generator,
    leverage_adjustment: str = "HC2",
) -> dict[str, float]:
    """Test the period term with a null-imposed site-cluster wild bootstrap.

    Hellinger-transformed lineage proportions are regressed on country in the
    reduced model and on country plus period in the full model. Rademacher
    weights are drawn independently for country-locality clusters and shared
    by every site-period row from the same cluster. This preserves dependence
    when one locality contributes to more than one period.
    """
    input_profiles = len(site_table)
    country_counts = site_table["country"].value_counts()
    singleton_countries = country_counts[country_counts.lt(2)].index.tolist()
    site_table = site_table[
        ~site_table["country"].isin(singleton_countries)
    ].reset_index(drop=True)
    y_cols = [c for c in site_table.columns if c.startswith("h_")]
    y = site_table[y_cols].to_numpy(float)
    x0, _ = design_matrix(site_table, include_period=False)
    x1, _ = design_matrix(site_table, include_period=True)
    m0 = _residual_maker(x0)
    m1 = _residual_maker(x1)
    s0 = _sse_from_residual_maker(m0, y)
    s1 = _sse_from_residual_maker(m1, y)
    r0, r1 = np.linalg.matrix_rank(x0), np.linalg.matrix_rank(x1)
    df_effect = r1 - r0
    df_resid = len(site_table) - r1
    raw_leverage = 1.0 - np.diag(m0)
    max_leverage = float(raw_leverage.max())
    invalid_reasons = []
    if df_effect <= 0:
        invalid_reasons.append("non-positive effect degrees of freedom")
    if df_resid <= 0:
        invalid_reasons.append("non-positive residual degrees of freedom")
    elif df_resid < 5:
        invalid_reasons.append("fewer than five residual degrees of freedom")
    if s1 <= np.finfo(float).eps:
        invalid_reasons.append("zero or near-zero full-model residual SSE")
    if max_leverage >= 1.0 - 1e-8:
        invalid_reasons.append("reduced-model leverage is one")
    observed = (
        ((s0 - s1) / df_effect) / (s1 / df_resid)
        if not any(
            reason
            in {
                "non-positive effect degrees of freedom",
                "non-positive residual degrees of freedom",
                "zero or near-zero full-model residual SSE",
            }
            for reason in invalid_reasons
        )
        else np.nan
    )
    cluster = (
        site_table["country"].astype(str)
        + "|||"
        + site_table["site"].astype(str)
    ).to_numpy()
    unique_clusters, cluster_index = np.unique(cluster, return_inverse=True)
    cluster_sizes = pd.Series(cluster).value_counts()
    base_result = {
        "pseudo_f": observed,
        "partial_r2": (s0 - s1) / s0 if s0 > 0 else np.nan,
        "df_effect": df_effect,
        "df_resid": df_resid,
        "cluster_wild_p": np.nan,
        "n_input_site_period_profiles": input_profiles,
        "n_site_period_profiles": len(site_table),
        "n_dropped_singleton_country_profiles": (
            input_profiles - len(site_table)
        ),
        "dropped_singleton_countries": "|".join(singleton_countries),
        "n_site_clusters": len(unique_clusters),
        "n_multiperiod_site_clusters": int((cluster_sizes > 1).sum()),
        "cluster_wild_resamples": n_resamples,
        "leverage_adjustment": leverage_adjustment,
        "max_reduced_model_leverage": max_leverage,
        "inference_valid": not invalid_reasons,
        "invalid_reason": "; ".join(invalid_reasons),
    }
    if invalid_reasons:
        return base_result

    fitted_null = y - m0 @ y
    raw_residual = m0 @ y
    leverage = np.clip(raw_leverage, 0.0, 1.0 - 1e-10)
    if leverage_adjustment == "HC2":
        exponent = 0.5
    elif leverage_adjustment == "HC3":
        exponent = 1.0
    else:
        raise ValueError(
            "leverage_adjustment must be 'HC2' or 'HC3', "
            f"not {leverage_adjustment!r}"
        )
    residual_null = raw_residual / (1.0 - leverage)[:, None] ** exponent
    ge = 0
    for _ in range(n_resamples):
        weights = rng.choice((-1.0, 1.0), size=len(unique_clusters))
        y_boot = fitted_null + residual_null * weights[cluster_index, None]
        s0_boot = _sse_from_residual_maker(m0, y_boot)
        s1_boot = _sse_from_residual_maker(m1, y_boot)
        fp = ((s0_boot - s1_boot) / df_effect) / (s1_boot / df_resid)
        ge += fp >= observed
    base_result["cluster_wild_p"] = (ge + 1) / (n_resamples + 1)
    return base_result


def site_cluster_effect_jackknife(
    site_table: pd.DataFrame,
) -> dict[str, float]:
    """Leave-one-locality-cluster jackknife interval for partial R²."""

    def effect(table: pd.DataFrame) -> float:
        y_cols = [c for c in table.columns if c.startswith("h_")]
        y = table[y_cols].to_numpy(float)
        x0, _ = design_matrix(table, include_period=False)
        x1, _ = design_matrix(table, include_period=True)
        if np.linalg.matrix_rank(x1) <= np.linalg.matrix_rank(x0):
            return np.nan
        s0 = _sse_from_residual_maker(_residual_maker(x0), y)
        s1 = _sse_from_residual_maker(_residual_maker(x1), y)
        if s0 <= 0:
            return np.nan
        return float((s0 - s1) / s0)

    source = site_table.copy()
    source["_cluster"] = (
        source["country"].astype(str) + "|||" + source["site"].astype(str)
    )
    observed = effect(source)
    all_clusters = source["_cluster"].drop_duplicates().to_numpy()
    jackknife = np.asarray(
        [
            effect(source[source["_cluster"] != cluster].copy())
            for cluster in all_clusters
        ],
        dtype=float,
    )
    jackknife = jackknife[np.isfinite(jackknife)]
    n_clusters = len(jackknife)
    jackknife_mean = float(jackknife.mean())
    standard_error = float(
        np.sqrt(
            (n_clusters - 1)
            / n_clusters
            * np.sum((jackknife - jackknife_mean) ** 2)
        )
    )
    critical = float(student_t.ppf(0.975, max(1, n_clusters - 1)))
    return {
        "partial_r2_jackknife_se": standard_error,
        "partial_r2_jackknife_ci_low": float(
            observed - critical * standard_error
        ),
        "partial_r2_jackknife_ci_high": float(
            observed + critical * standard_error
        ),
        "partial_r2_jackknife_clusters": n_clusters,
    }


def dispersion_distance_table(site_table: pd.DataFrame) -> pd.DataFrame:
    """Distances to country-period centroids for a model-aligned diagnostic."""
    y_cols = [c for c in site_table.columns if c.startswith("h_")]
    rows = []
    for (_, _), group in site_table.groupby(
        ["country", "analysis_bin"], observed=True
    ):
        if len(group) < 2:
            continue
        y = group[y_cols].to_numpy(float)
        centroid = y.mean(axis=0)
        correction = np.sqrt(len(group) / (len(group) - 1))
        distances = np.linalg.norm(y - centroid, axis=1) * correction
        for (_, record), distance in zip(group.iterrows(), distances):
            rows.append(
                {
                    "country": record["country"],
                    "site": record["site"],
                    "analysis_bin": record["analysis_bin"],
                    "n_calls": record["n_calls"],
                    "h_distance": float(distance),
                }
            )
    return pd.DataFrame(rows)


def model_residual_diagnostics(
    site_table: pd.DataFrame, marker: str
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Profile-level leverage and residual-size diagnostics."""
    y_cols = [c for c in site_table.columns if c.startswith("h_")]
    y = site_table[y_cols].to_numpy(float)
    x0, _ = design_matrix(site_table, include_period=False)
    x1, _ = design_matrix(site_table, include_period=True)
    m0 = _residual_maker(x0)
    m1 = _residual_maker(x1)
    table = site_table[
        ["country", "site", "analysis_bin", "n_calls"]
    ].copy()
    table.insert(0, "marker", marker)
    table["log1p_n_calls"] = np.log1p(table["n_calls"].to_numpy(float))
    table["reduced_residual_norm"] = np.linalg.norm(m0 @ y, axis=1)
    table["full_residual_norm"] = np.linalg.norm(m1 @ y, axis=1)
    table["reduced_leverage"] = 1.0 - np.diag(m0)
    table["full_leverage"] = 1.0 - np.diag(m1)
    summary = {
        "marker": marker,
        "n_profiles": len(table),
        "max_reduced_leverage": float(table["reduced_leverage"].max()),
        "max_full_leverage": float(table["full_leverage"].max()),
        "pearson_full_residual_norm_vs_log1p_calls": float(
            table["full_residual_norm"].corr(
                table["log1p_n_calls"], method="pearson"
            )
        ),
        "spearman_full_residual_norm_vs_log1p_calls": float(
            table["full_residual_norm"].corr(
                table["log1p_n_calls"], method="spearman"
            )
        ),
    }
    return table, summary


def repeated_site_period_test(
    site_table: pd.DataFrame, n_permutations: int, rng: np.random.Generator
) -> dict[str, float]:
    """Conservative period test using only localities observed in >1 bin.

    Site fixed effects absorb all time-invariant locality differences. Period
    labels are then permuted only among rows from the same locality. The test
    asks a narrower within-locality question and has less information than the
    all-site cluster-wild analysis.
    """
    table = site_table.copy()
    table["_cluster"] = (
        table["country"].astype(str) + "|||" + table["site"].astype(str)
    )
    counts = table["_cluster"].value_counts()
    table = table[table["_cluster"].map(counts).gt(1)].reset_index(drop=True)
    empty_result = {
        "repeated_site_pseudo_f": np.nan,
        "repeated_site_df_effect": 0,
        "repeated_site_df_resid": 0,
        "repeated_site_permutation_p": np.nan,
        "repeated_site_profiles": len(table),
        "repeated_sites": table["_cluster"].nunique(),
        "repeated_site_permutations": n_permutations,
    }
    if table.empty:
        return empty_result
    y_cols = [c for c in table.columns if c.startswith("h_")]
    y = table[y_cols].to_numpy(float)

    site_design = pd.get_dummies(
        table["_cluster"], prefix="site", drop_first=False
    ).astype(float)

    def fixed_site_design(current: pd.DataFrame, include_period: bool) -> np.ndarray:
        frames = [site_design]
        if include_period:
            period = pd.get_dummies(
                pd.Categorical(
                    current["analysis_bin"], categories=BIN_LABELS
                ),
                prefix="period",
                drop_first=True,
            ).astype(float)
            period.index = current.index
            frames.append(period)
        return pd.concat(frames, axis=1).to_numpy()

    x0 = fixed_site_design(table, include_period=False)
    x1 = fixed_site_design(table, include_period=True)
    m0 = _residual_maker(x0)
    m1 = _residual_maker(x1)
    s0 = _sse_from_residual_maker(m0, y)
    s1 = _sse_from_residual_maker(m1, y)
    r0, r1 = np.linalg.matrix_rank(x0), np.linalg.matrix_rank(x1)
    df_effect = r1 - r0
    df_resid = len(table) - r1
    if df_effect <= 0 or df_resid <= 0 or s1 <= 0:
        return empty_result
    observed = ((s0 - s1) / df_effect) / (s1 / df_resid)

    periods = table["analysis_bin"].to_numpy(copy=True)
    clusters = table["_cluster"].to_numpy()
    ge = 0
    for _ in range(n_permutations):
        shuffled = periods.copy()
        for cluster in np.unique(clusters):
            idx = np.flatnonzero(clusters == cluster)
            shuffled[idx] = rng.permutation(shuffled[idx])
        permuted = table.copy()
        permuted["analysis_bin"] = shuffled
        xp = fixed_site_design(permuted, include_period=True)
        mp = _residual_maker(xp)
        sp = _sse_from_residual_maker(mp, y)
        fp = ((s0 - sp) / df_effect) / (sp / df_resid)
        ge += fp >= observed
    return {
        "repeated_site_pseudo_f": observed,
        "repeated_site_df_effect": df_effect,
        "repeated_site_df_resid": df_resid,
        "repeated_site_permutation_p": (ge + 1) / (n_permutations + 1),
        "repeated_site_profiles": len(table),
        "repeated_sites": table["_cluster"].nunique(),
        "repeated_site_permutations": n_permutations,
    }


def holm_adjust(values: pd.Series) -> pd.Series:
    order = np.argsort(values.to_numpy(float))
    ordered = values.iloc[order].reset_index(drop=True)
    adjusted = []
    running = 0.0
    m = len(ordered)
    for i, value in enumerate(ordered):
        running = max(running, (m - i) * float(value))
        adjusted.append(min(1.0, running))
    result = pd.Series(index=range(len(values)), dtype=float)
    for rank, original_index in enumerate(order):
        result.iloc[original_index] = adjusted[rank]
    return result


def unrelated_subset(df: pd.DataFrame, marker_call_col: str) -> pd.DataFrame:
    quality = "mt_call_coverage" if marker_call_col.startswith("mt") else "snp_2m"
    use = df[df[marker_call_col] != ""].copy()
    use["_quality"] = pd.to_numeric(use[quality], errors="coerce").fillna(-1)
    return (
        use.sort_values(
            ["kin_component_2d", "_quality"], ascending=[True, False]
        )
        .drop_duplicates("kin_component_2d")
        .drop(columns="_quality")
    )


def mean_adjacent_tv(
    df: pd.DataFrame, marker_col: str, categories: list[str]
) -> float:
    p = profile(df, marker_col, categories, site_balanced=True)
    vals = []
    for i in range(len(p) - 1):
        vals.append(
            total_variation(
                p.loc[i, categories].to_numpy(float),
                p.loc[i + 1, categories].to_numpy(float),
            )
        )
    vals = [v for v in vals if np.isfinite(v)]
    return float(np.mean(vals)) if vals else np.nan


def date_uncertainty(
    df: pd.DataFrame,
    marker_col: str,
    categories: list[str],
    draws: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Run an assumption-based chronological-bin assignment sensitivity.

    Direct dates are perturbed with a normal distribution in uncalibrated BP
    using the catalogue mean and standard deviation. Indirect dates are
    perturbed uniformly over a moment-matched interval. These draws are not
    samples from calibrated radiocarbon probability distributions and must not
    be interpreted as posterior chronological uncertainty.

    Calling this function for different markers with identically initialized
    RNG streams gives the markers the same joint individual-level date
    scenarios, allowing like-for-like descriptive comparison.
    """
    rows = []
    base = df.copy()
    mean = base["date_bp"].to_numpy(float)
    sd = base["date_sd_bp"].fillna(0).to_numpy(float)
    direct = base["direct_date"].to_numpy(bool)
    for i in range(draws):
        sampled = np.empty(len(base), dtype=float)
        sampled[direct] = rng.normal(mean[direct], sd[direct])
        half_range = np.sqrt(3.0) * sd[~direct]
        sampled[~direct] = rng.uniform(
            mean[~direct] - half_range, mean[~direct] + half_range
        )
        sampled = np.clip(sampled, 0, None)
        tmp = base.copy()
        tmp["analysis_bin"] = pd.cut(
            1950 - sampled,
            bins=BIN_EDGES,
            labels=BIN_LABELS,
            right=False,
        )
        rows.append(
            {
                "scenario_draw": i,
                "analysis_type": (
                    "chronological_bin_assignment_scenario_sensitivity"
                ),
                "date_scenario_model": (
                    "direct_normal_BP; indirect_uniform_BP_moment_matched"
                ),
                "is_calibrated_date_posterior": False,
                "mean_adjacent_tv": mean_adjacent_tv(
                    tmp, marker_col, categories
                ),
            }
        )
    return pd.DataFrame(rows)


def exact_crosswalk(
    canonical: pd.DataFrame, amtdb: pd.DataFrame, aychr: pd.DataFrame
) -> pd.DataFrame:
    alias_to_ids: dict[str, set[str]] = defaultdict(set)
    for _, r in canonical.iterrows():
        aliases = [
            r["individual_id"],
            r["genetic_id"].split(".")[0],
            r["skeletal_code"],
        ]
        for alias in aliases:
            key = normalize_id(alias)
            if key:
                alias_to_ids[key].add(r["individual_id"])

    rows = []
    for source, table, id_col, alt_col, country_col in [
        ("AmtDB v1.009", amtdb, "identifier", "alternative_identifiers", "country"),
        ("aYChr-DB v5", aychr, "Published ID", None, "Country"),
    ]:
        subset = table[table[country_col].isin(COUNTRIES)].copy()
        for _, r in subset.iterrows():
            aliases = [r[id_col]]
            if alt_col:
                aliases.extend(re.split(r"[;,|/]", str(r[alt_col])))
            candidates: set[str] = set()
            for alias in aliases:
                candidates.update(alias_to_ids.get(normalize_id(alias), set()))
            manual = MANUAL_CROSSWALK.get((source, str(r[id_col]).strip()))
            if manual:
                if manual not in set(canonical["individual_id"]):
                    raise ValueError(
                        f"Manual crosswalk target {manual!r} is not present in AADR"
                    )
                candidates = {manual}
                status = "manual_verified"
            else:
                status = (
                    "exact_unique"
                    if len(candidates) == 1
                    else "exact_ambiguous"
                    if len(candidates) > 1
                    else "not_exactly_matched"
                )
            rows.append(
                {
                    "source": source,
                    "source_id": r[id_col],
                    "country": r[country_col],
                    "match_status": status,
                    "aadr_individual_ids": ";".join(sorted(candidates)),
                }
            )
    return pd.DataFrame(rows)


def extended_legacy_mtdna(
    amtdb: pd.DataFrame, crosswalk: pd.DataFrame
) -> pd.DataFrame:
    """Return Central Asian AmtDB records not represented in AADR.

    These rows are kept as an extended evidence stratum. They are not pooled
    into AADR-based inferential tests because laboratory methods, QC fields and
    curation histories differ.
    """
    amt_cross = crosswalk[crosswalk["source"] == "AmtDB v1.009"][
        ["source_id", "match_status", "aadr_individual_ids"]
    ].copy()
    subset = amtdb[amtdb["country"].isin(COUNTRIES)].copy()
    subset = subset.merge(
        amt_cross, left_on="identifier", right_on="source_id", how="left"
    )
    legacy = subset[subset["match_status"] == "not_exactly_matched"].copy()
    legacy["year_from"] = pd.to_numeric(legacy["year_from"], errors="coerce")
    legacy["year_to"] = pd.to_numeric(legacy["year_to"], errors="coerce")
    legacy["year_mid_ce"] = (legacy["year_from"] + legacy["year_to"]) / 2
    legacy["analysis_bin"] = pd.cut(
        legacy["year_mid_ce"],
        bins=BIN_EDGES,
        labels=BIN_LABELS,
        right=False,
    ).astype(str)
    legacy["mt_l1"] = legacy["mt_hg"].map(lambda x: major_haplogroup(x, "mt"))
    legacy["evidence_stratum"] = "legacy_mtdna_not_in_aadr"
    legacy["pooling_rule"] = (
        "Catalogue and extended sensitivity only; excluded from AADR primary tests"
    )
    keep = [
        "evidence_stratum",
        "identifier",
        "alternative_identifiers",
        "country",
        "site",
        "site_detail",
        "latitude",
        "longitude",
        "culture",
        "epoch",
        "year_from",
        "year_to",
        "year_mid_ce",
        "analysis_bin",
        "sex",
        "mt_hg",
        "mt_l1",
        "sequence_source",
        "avg_coverage",
        "reference_name",
        "reference_link",
        "data_link",
        "match_status",
        "pooling_rule",
    ]
    return legacy[keep].sort_values(
        ["country", "year_mid_ce", "site", "identifier"]
    )


def _fixed_margin_chi_square_monte_carlo(
    yes: np.ndarray,
    total: np.ndarray,
    n_resamples: int,
    rng: np.random.Generator,
) -> dict[str, float | int | bool | str]:
    """Pearson statistic with a reproducible fixed-margin Monte Carlo p-value.

    This is a conditional simulation under independence, not an exact test.
    Empty rows are removed before evaluating the table. The +1 correction
    prevents a zero Monte Carlo p-value.
    """
    yes = np.asarray(yes, dtype=int)
    total = np.asarray(total, dtype=int)
    if np.any(yes < 0) or np.any(total < yes):
        raise ValueError("Callability counts must satisfy 0 <= yes <= total")
    keep = total > 0
    yes = yes[keep]
    total = total[keep]
    total_yes = int(yes.sum())
    grand_total = int(total.sum())
    if (
        len(total) < 2
        or grand_total == 0
        or total_yes == 0
        or total_yes == grand_total
    ):
        return {
            "chi2": np.nan,
            "df": max(0, len(total) - 1),
            "asymptotic_p": np.nan,
            "monte_carlo_p": np.nan,
            "monte_carlo_resamples": n_resamples,
            "min_expected_count": np.nan,
            "n_expected_lt_5": 0,
            "fraction_expected_lt_5": np.nan,
            "sparse_expected_cells": False,
            "p_value": np.nan,
            "p_value_method": "not_estimable",
        }

    contingency = np.column_stack([yes, total - yes])
    chi2, asymptotic_p, dof, expected = chi2_contingency(
        contingency, correction=False
    )
    expected_yes = expected[:, 0]
    expected_no = expected[:, 1]
    simulated_yes = rng.multivariate_hypergeometric(
        total, total_yes, size=n_resamples
    )
    simulated_no = total[np.newaxis, :] - simulated_yes
    simulated_statistics = (
        ((simulated_yes - expected_yes) ** 2 / expected_yes).sum(axis=1)
        + ((simulated_no - expected_no) ** 2 / expected_no).sum(axis=1)
    )
    monte_carlo_p = float(
        (1 + np.count_nonzero(simulated_statistics >= chi2 - 1e-12))
        / (n_resamples + 1)
    )
    sparse = bool(np.any(expected < 5))
    return {
        "chi2": float(chi2),
        "df": int(dof),
        "asymptotic_p": float(asymptotic_p),
        "monte_carlo_p": monte_carlo_p,
        "monte_carlo_resamples": int(n_resamples),
        "min_expected_count": float(expected.min()),
        "n_expected_lt_5": int((expected < 5).sum()),
        "fraction_expected_lt_5": float((expected < 5).mean()),
        "sparse_expected_cells": sparse,
        "p_value": monte_carlo_p if sparse else float(asymptotic_p),
        "p_value_method": (
            "fixed_margin_monte_carlo" if sparse else "asymptotic_pearson"
        ),
    }


def callability_table(
    df: pd.DataFrame,
    monte_carlo_resamples: int = 99999,
    seed: int = 20260725,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for period in BIN_LABELS:
        g = df[df["analysis_bin"].astype(str) == period]
        males = g[g["molecular_sex"].astype(str).str.startswith("M")]
        rows.append(
            {
                "analysis_bin": period,
                "n_individuals": len(g),
                "mt_calls": int((g["mt_call"] != "").sum()),
                "mt_call_rate": float((g["mt_call"] != "").mean()),
                "n_molecular_males": len(males),
                "y_calls": int((males["y_call"] != "").sum()),
                "y_call_rate": (
                    float((males["y_call"] != "").mean()) if len(males) else np.nan
                ),
            }
        )
    table = pd.DataFrame(rows)
    tests = []
    for marker, yes, total in [
        ("mtDNA", table["mt_calls"], table["n_individuals"]),
        ("Y", table["y_calls"], table["n_molecular_males"]),
    ]:
        result = _fixed_margin_chi_square_monte_carlo(
            yes.to_numpy(int),
            total.to_numpy(int),
            monte_carlo_resamples,
            named_rng(seed, f"{marker}:callability-fixed-margin-monte-carlo"),
        )
        result.update(
            {
                "marker": marker,
                "note": (
                    "Descriptive Pearson test; fixed-margin Monte Carlo used "
                    "for sparse expected cells. Monte Carlo, not exact; does "
                    "not account for site/study clustering."
                ),
            }
        )
        tests.append(result)
    return table, pd.DataFrame(tests)


def count_matrix(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for country in COUNTRIES:
        for period in BIN_LABELS:
            g = df[
                (df["country"] == country)
                & (df["analysis_bin"].astype(str) == period)
            ]
            rows.append(
                {
                    "country": country,
                    "analysis_bin": period,
                    "individuals": len(g),
                    "sites": g["locality"].nunique(),
                    "studies": g["publication"].nunique(),
                    "mt_calls": int((g["mt_call"] != "").sum()),
                    "y_calls": int((g["y_call"] != "").sum()),
                    "molecular_males": int(
                        g["molecular_sex"].astype(str).str.startswith("M").sum()
                    ),
                    "direct_dates": int(g["direct_date"].sum()),
                }
            )
    return pd.DataFrame(rows)


def analysis_cell_adequacy(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for country in COUNTRIES:
        for period in BIN_LABELS:
            g = df[
                (df["country"] == country)
                & (df["analysis_bin"].astype(str) == period)
            ]
            for marker, col in [
                ("mtDNA", "mt_l1_pooled"),
                ("Y", "y_l1_pooled"),
            ]:
                called = g[g[col] != ""]
                n_calls = len(called)
                n_sites = called["locality"].nunique()
                rows.append(
                    {
                        "country": country,
                        "analysis_bin": period,
                        "marker": marker,
                        "n_calls": n_calls,
                        "n_sites": n_sites,
                        "meets_descriptive_frequency_threshold": (
                            n_calls >= 10 and n_sites >= 3
                        ),
                        "threshold_rule": "At least 10 calls from at least 3 sites",
                    }
                )
    return pd.DataFrame(rows)


def site_dominance(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for period in BIN_LABELS:
        for marker, col in [("mtDNA", "mt_l1_pooled"), ("Y", "y_l1_pooled")]:
            g = df[
                (df["analysis_bin"].astype(str) == period) & (df[col] != "")
            ]
            if g.empty:
                continue
            by_site = (
                g.groupby(["country", "locality"], dropna=False)
                .size()
                .sort_values(ascending=False)
            )
            top_country, top_site = by_site.index[0]
            by_country = g.groupby("country").size().sort_values(ascending=False)
            rows.append(
                {
                    "analysis_bin": period,
                    "marker": marker,
                    "n_calls": len(g),
                    "n_sites": len(by_site),
                    "largest_site_country": top_country,
                    "largest_site": top_site,
                    "largest_site_calls": int(by_site.iloc[0]),
                    "largest_site_share": float(by_site.iloc[0] / len(g)),
                    "largest_country": by_country.index[0],
                    "largest_country_calls": int(by_country.iloc[0]),
                    "largest_country_share": float(by_country.iloc[0] / len(g)),
                }
            )
    return pd.DataFrame(rows)


def figure_sampling(df: pd.DataFrame, counts: pd.DataFrame, path: Path) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    fig = plt.figure(figsize=(13.5, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35])
    ax = fig.add_subplot(gs[0, 0])
    sites = (
        df.groupby(["country", "locality"], dropna=False)
        .agg(longitude=("longitude", "median"), latitude=("latitude", "median"), n=("individual_id", "size"))
        .reset_index()
    )
    for country in COUNTRIES:
        g = sites[sites["country"] == country]
        ax.scatter(
            g["longitude"],
            g["latitude"],
            s=18 + 7 * np.sqrt(g["n"]),
            color=COUNTRY_COLORS[country],
            alpha=0.78,
            edgecolor="white",
            linewidth=0.5,
            label=country,
        )
    ax.set(
        xlabel="Longitude",
        ylabel="Latitude",
        title="Published archaeological sites in the AADR core catalogue",
        xlim=(45, 90),
        ylim=(35, 56),
    )
    ax.legend(frameon=False, fontsize=8, loc="lower left")

    ax2 = fig.add_subplot(gs[0, 1])
    pivot = counts.pivot(
        index="country", columns="analysis_bin", values="individuals"
    ).reindex(index=COUNTRIES, columns=BIN_LABELS)
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".0f",
        cmap="YlGnBu",
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Unique individuals"},
        ax=ax2,
    )
    ax2.set(
        xlabel="Pre-defined chronological bin",
        ylabel="",
        title="Sampling is strongly uneven across country and time",
    )
    ax2.tick_params(axis="x", rotation=48)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_composition(
    site_profile: pd.DataFrame,
    categories: list[str],
    marker: str,
    path: Path,
) -> None:
    sns.set_theme(style="white", font_scale=0.95)
    palette = sns.color_palette("tab20", n_colors=len(categories))
    fig, ax = plt.subplots(figsize=(12.5, 5.8), constrained_layout=True)
    bottom = np.zeros(len(site_profile))
    x = np.arange(len(site_profile))
    for cat, color in zip(categories, palette):
        values = site_profile[cat].fillna(0).to_numpy(float)
        ax.bar(x, values, bottom=bottom, label=cat, color=color, width=0.82)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(BIN_LABELS, rotation=42, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Mean site-balanced composition")
    ax.set_title(
        f"{marker} L1 composition among published sites (equal weight per site)"
    )
    for i, r in site_profile.iterrows():
        ax.text(
            i,
            1.025,
            f"n={int(r.n_calls)}; s={int(r.n_sites)}",
            ha="center",
            va="bottom",
            fontsize=7,
            rotation=90,
        )
    ax.legend(
        title="Broad lineage label",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
        ncol=1,
    )
    sns.despine(ax=ax)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_diversity_turnover(
    div_summary: pd.DataFrame, tv_summary: pd.DataFrame, path: Path
) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), constrained_layout=True)
    colors = {"mtDNA": "#D95F02", "Y": "#1B9E77"}
    for marker, g in div_summary.groupby("marker"):
        g = g.set_index("analysis_bin").reindex(BIN_LABELS).reset_index()
        x = np.arange(len(g))
        axes[0].plot(
            x, g["estimate"], marker="o", lw=2, color=colors[marker], label=marker
        )
        axes[0].fill_between(
            x,
            g["ci_low"],
            g["ci_high"],
            color=colors[marker],
            alpha=0.18,
        )
    axes[0].set_xticks(np.arange(len(BIN_LABELS)))
    axes[0].set_xticklabels(BIN_LABELS, rotation=43, ha="right")
    axes[0].set_ylabel("Effective number of L1 lineages (Hill q=1)")
    axes[0].set_title("Observed site-balanced diversity; cluster-bootstrap interval")
    axes[0].legend(frameon=False)

    transitions = [f"{a} -> {b}" for a, b in zip(BIN_LABELS[:-1], BIN_LABELS[1:])]
    for marker, g in tv_summary.groupby("marker"):
        g = g.set_index("transition").reindex(transitions).reset_index()
        x = np.arange(len(g))
        axes[1].vlines(
            x,
            g["ci_low"],
            g["ci_high"],
            lw=1.8,
            color=colors[marker],
        )
        axes[1].scatter(
            x,
            g["estimate"],
            marker="o",
            color=colors[marker],
            label=marker,
            zorder=3,
        )
    axes[1].set_xticks(np.arange(len(transitions)))
    axes[1].set_xticklabels(transitions, rotation=43, ha="right")
    axes[1].set_ylabel("Total-variation turnover (0-1)")
    axes[1].set_title("Observed adjacent-bin turnover; cluster-bootstrap interval")
    axes[1].legend(frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_sensitivity(table: pd.DataFrame, path: Path) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    order = list(dict.fromkeys(table["analysis"].tolist()))
    markers = ["mtDNA", "Y"]
    offsets = {"mtDNA": -0.16, "Y": 0.16}
    colors = {"mtDNA": "#D95F02", "Y": "#1B9E77"}
    for marker in markers:
        g = table[table["marker"] == marker].set_index("analysis").reindex(order)
        y = np.arange(len(order)) + offsets[marker]
        ax.scatter(
            g["mean_adjacent_tv"],
            y,
            color=colors[marker],
            s=45,
            label=marker,
            zorder=3,
        )
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Mean adjacent-bin total-variation turnover")
    ax.set_title("Sensitivity of the descriptive turnover estimate")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_database_coverage(table: pd.DataFrame, path: Path) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.95)
    fig, ax = plt.subplots(figsize=(10.5, 5.7), constrained_layout=True)
    pivot = table.pivot(
        index="country", columns="source", values="records"
    ).reindex(COUNTRIES)
    pivot.plot(kind="bar", ax=ax, color=["#2C7FB8", "#F03B20", "#756BB1"])
    ax.set_ylabel("Database records before cross-database deduplication")
    ax.set_xlabel("")
    ax.set_title("Specialized legacy databases cover only a subset of the region")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(title="")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_hashes = validate_frozen_inputs(args)
    out = args.outdir
    tables_dir = out / "tables"
    figures_dir = out / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    aadr = pd.read_csv(args.aadr, sep="\t", dtype=str, keep_default_na=False)
    amtdb = pd.read_csv(args.amtdb, dtype=str, keep_default_na=False)
    aychr = pd.read_excel(args.aychr, dtype=str).fillna("")
    canonical, dedup_audit = choose_canonical_aadr(aadr)
    analysis = canonical[canonical["analysis_included"]].copy()

    analysis["mt_l1_pooled"], mt_categories = pool_rare(analysis["mt_l1"], 5)
    analysis["y_l1_pooled"], y_categories = pool_rare(analysis["y_l1"], 5)
    resolution_analysis, y_family_categories, y_family_source_column = (
        add_y_resolution_sensitivity_encoding(analysis, 5)
    )
    # Carry pooled labels back to the full catalogue where possible.
    mt_keep = set(mt_categories) - {"Other"}
    y_keep = set(y_categories) - {"Other"}
    canonical["mt_l1_pooled"] = canonical["mt_l1"].where(
        canonical["mt_l1"].isin(mt_keep),
        np.where(canonical["mt_l1"].eq(""), "", "Other"),
    )
    canonical["y_l1_pooled"] = canonical["y_l1"].where(
        canonical["y_l1"].isin(y_keep),
        np.where(canonical["y_l1"].eq(""), "", "Other"),
    )

    counts = count_matrix(analysis)
    adequacy = analysis_cell_adequacy(analysis)
    dominance = site_dominance(analysis)
    callability, callability_tests = callability_table(
        analysis,
        monte_carlo_resamples=args.callability_resamples,
        seed=args.seed,
    )
    crosswalk = exact_crosswalk(canonical, amtdb, aychr)
    legacy_mt = extended_legacy_mtdna(amtdb, crosswalk)

    db_rows = []
    for country in COUNTRIES:
        db_rows.extend(
            [
                {
                    "country": country,
                    "source": "AADR v66.p1",
                    "records": int((canonical["country"] == country).sum()),
                },
                {
                    "country": country,
                    "source": "AmtDB v1.009",
                    "records": int((amtdb["country"] == country).sum()),
                },
                {
                    "country": country,
                    "source": "aYChr-DB v5",
                    "records": int((aychr["Country"] == country).sum()),
                },
            ]
        )
    db_coverage = pd.DataFrame(db_rows)

    profiles = {}
    div_summaries = []
    tv_summaries = []
    global_tests = []
    dispersion_tests = []
    dispersion_summaries = []
    residual_diagnostic_tables = []
    residual_diagnostic_summaries = []
    bootstrap_diagnostics: dict[str, dict[str, int]] = {}
    for marker, col, cats in [
        ("mtDNA", "mt_l1_pooled", mt_categories),
        ("Y", "y_l1_pooled", y_categories),
    ]:
        ind = profile(analysis, col, cats, site_balanced=False)
        site = profile(analysis, col, cats, site_balanced=True)
        profiles[(marker, "individual")] = ind
        profiles[(marker, "site")] = site
        ind.to_csv(
            tables_dir / f"composition_{marker.lower()}_individual_weighted.csv",
            index=False,
        )
        site.to_csv(
            tables_dir / f"composition_{marker.lower()}_site_balanced.csv",
            index=False,
        )
        div_draws, tv_draws, marker_bootstrap_diagnostics = (
            bootstrap_site_profiles(
                analysis,
                col,
                cats,
                args.bootstrap,
                named_rng(args.seed, f"{marker}:site-cluster-bootstrap"),
            )
        )
        bootstrap_diagnostics[marker] = marker_bootstrap_diagnostics
        observed_diversity, observed_turnover = observed_profile_statistics(
            site, cats
        )
        div = summarize_bootstrap(
            div_draws,
            "q1",
            "analysis_bin",
            observed_diversity,
        )
        div["marker"] = marker
        tv = summarize_bootstrap(
            tv_draws,
            "tv",
            "transition",
            observed_turnover,
        )
        tv["marker"] = marker
        div_summaries.append(div)
        tv_summaries.append(tv)
        div_draws.to_csv(
            tables_dir / f"bootstrap_diversity_{marker.lower()}_draws.csv",
            index=False,
        )
        tv_draws.to_csv(
            tables_dir / f"bootstrap_turnover_{marker.lower()}_draws.csv",
            index=False,
        )
        site_table = site_profile_table(analysis, col, cats, min_calls=1)
        site_table.to_csv(
            tables_dir / f"site_profiles_{marker.lower()}.csv", index=False
        )
        result = site_cluster_wild_period_test(
            site_table,
            args.permutations,
            named_rng(args.seed, f"{marker}:primary-cluster-wild"),
        )
        result.update(site_cluster_effect_jackknife(site_table))
        result.update(
            repeated_site_period_test(
                site_table,
                args.permutations,
                named_rng(args.seed, f"{marker}:repeated-site-permutation"),
            )
        )
        result["marker"] = marker
        global_tests.append(result)

        residual_table, residual_summary = model_residual_diagnostics(
            site_table, marker
        )
        residual_diagnostic_tables.append(residual_table)
        residual_diagnostic_summaries.append(residual_summary)

        dispersion = dispersion_distance_table(site_table)
        dispersion.to_csv(
            tables_dir / f"dispersion_profiles_{marker.lower()}.csv",
            index=False,
        )
        dispersion_test = site_cluster_wild_period_test(
            dispersion,
            args.permutations,
            named_rng(args.seed, f"{marker}:dispersion-cluster-wild"),
        )
        dispersion_test["marker"] = marker
        dispersion_tests.append(dispersion_test)
        dispersion_summary = (
            dispersion.groupby("analysis_bin", observed=True)["h_distance"]
            .agg(
                n_profiles="size",
                median="median",
                q1=lambda x: x.quantile(0.25),
                q3=lambda x: x.quantile(0.75),
            )
            .reset_index()
        )
        dispersion_summary["marker"] = marker
        dispersion_summaries.append(dispersion_summary)

    diversity = pd.concat(div_summaries, ignore_index=True)
    turnover = pd.concat(tv_summaries, ignore_index=True)
    tests = pd.DataFrame(global_tests)
    tests["holm_cluster_wild_p"] = holm_adjust(tests["cluster_wild_p"])
    tests["holm_repeated_site_p"] = holm_adjust(
        tests["repeated_site_permutation_p"]
    )
    tests = tests.sort_values("marker").reset_index(drop=True)
    dispersion_test_table = pd.DataFrame(dispersion_tests)
    dispersion_test_table["holm_cluster_wild_p"] = holm_adjust(
        dispersion_test_table["cluster_wild_p"]
    )
    dispersion_test_table = dispersion_test_table.sort_values(
        "marker"
    ).reset_index(drop=True)
    dispersion_summary_table = pd.concat(
        dispersion_summaries, ignore_index=True
    )
    residual_diagnostic_table = pd.concat(
        residual_diagnostic_tables, ignore_index=True
    )
    residual_diagnostic_summary_table = pd.DataFrame(
        residual_diagnostic_summaries
    )

    sensitivity_rows = []
    filters = {
        "All marker-qualified calls": analysis,
        "AADR assessment-positive subset": analysis[analysis["strict_qc"]],
        "Exclude population-outlier labels": analysis[
            ~analysis["population_outlier"]
        ],
        "Direct dates only": analysis[analysis["direct_date"]],
    }
    for marker, call_col, pooled_col, cats in [
        ("mtDNA", "mt_call", "mt_l1_pooled", mt_categories),
        ("Y", "y_call", "y_l1_pooled", y_categories),
    ]:
        marker_filters = dict(filters)
        unrelated = unrelated_subset(analysis, call_col)
        marker_filters["One representative per <=2d kin component"] = unrelated
        paired = analysis[
            analysis["molecular_sex"].astype(str).str.startswith("M")
            & (analysis["mt_call"] != "")
            & (analysis["y_call"] != "")
        ]
        marker_filters["Male-paired subset"] = paired
        for name, subset in marker_filters.items():
            sensitivity_rows.append(
                {
                    "analysis": name,
                    "marker": marker,
                    "n_calls": int((subset[pooled_col] != "").sum()),
                    "n_sites": subset.loc[
                        subset[pooled_col] != "", "locality"
                    ].nunique(),
                    "mean_adjacent_tv": mean_adjacent_tv(
                        subset, pooled_col, cats
                    ),
                }
            )
        for country in COUNTRIES:
            subset = analysis[analysis["country"] != country]
            sensitivity_rows.append(
                {
                    "analysis": f"Leave out country: {country}",
                    "marker": marker,
                    "n_calls": int((subset[pooled_col] != "").sum()),
                    "n_sites": subset.loc[
                        subset[pooled_col] != "", "locality"
                    ].nunique(),
                    "mean_adjacent_tv": mean_adjacent_tv(
                        subset, pooled_col, cats
                    ),
                }
            )
    sensitivity = pd.DataFrame(sensitivity_rows)
    paired_draws, paired_summary, paired_bootstrap_diagnostics = (
        paired_marker_turnover_bootstrap(
            analysis,
            mt_categories,
            y_categories,
            args.paired_bootstrap,
            named_rng(args.seed, "paired-marker:site-cluster-bootstrap"),
            y_marker_col="y_l1_pooled",
            y_encoding_label="broad_L1",
        )
    )
    family_draws, family_summary, family_bootstrap_diagnostics = (
        paired_marker_turnover_bootstrap(
            resolution_analysis,
            mt_categories,
            y_family_categories,
            args.paired_bootstrap,
            # The identical stream isolates the effect of the encoding: the
            # same cluster multiplicities are used for both resolution rows.
            named_rng(args.seed, "paired-marker:site-cluster-bootstrap"),
            y_marker_col="y_isogg_prefix_family_pooled",
            y_encoding_label="AADR_ISOGG_prefix_family",
        )
    )
    paired_resolution_sensitivity = pd.concat(
        [paired_summary, family_summary], ignore_index=True
    )
    paired_resolution_sensitivity["category_count"] = [
        len(y_categories),
        len(y_family_categories),
    ]
    paired_resolution_sensitivity["categories"] = [
        ";".join(y_categories),
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

    date_draw_summaries = []
    for marker, col, cats in [
        ("mtDNA", "mt_l1_pooled", mt_categories),
        ("Y", "y_l1_pooled", y_categories),
    ]:
        draws = date_uncertainty(
            analysis,
            col,
            cats,
            args.date_draws,
            named_rng(args.seed, "shared:date-assignment-scenarios"),
        )
        draws["marker"] = marker
        draws.to_csv(
            tables_dir / f"date_uncertainty_{marker.lower()}_draws.csv",
            index=False,
        )
        date_draw_summaries.append(
            {
                "marker": marker,
                "analysis_type": (
                    "chronological_bin_assignment_scenario_sensitivity"
                ),
                "scenario_draws": args.date_draws,
                "observed_mean_adjacent_tv": mean_adjacent_tv(
                    analysis, col, cats
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
    date_uncertainty_summary = pd.DataFrame(date_draw_summaries)

    public_canonical = mask_coordinates(canonical)
    public_analysis = mask_coordinates(analysis)
    public_legacy_mt = mask_coordinates(legacy_mt)
    public_canonical.sort_values(
        ["country", "date_bp"], ascending=[True, False]
    ).to_csv(
        tables_dir / "aadr_central_asia_unique_individual_catalogue.csv", index=False
    )
    public_analysis.sort_values(
        ["country", "date_bp"], ascending=[True, False]
    ).to_csv(
        tables_dir / "aadr_primary_analysis_catalogue.csv", index=False
    )
    dedup_audit.to_csv(tables_dir / "aadr_deduplication_audit.csv", index=False)
    crosswalk.to_csv(tables_dir / "cross_database_exact_id_audit.csv", index=False)
    public_legacy_mt.to_csv(
        tables_dir / "extended_legacy_mtdna_not_in_aadr.csv", index=False
    )
    db_coverage.to_csv(tables_dir / "database_coverage_by_country.csv", index=False)
    counts.to_csv(tables_dir / "counts_country_by_bin.csv", index=False)
    adequacy.to_csv(tables_dir / "country_period_marker_adequacy.csv", index=False)
    dominance.to_csv(tables_dir / "site_and_country_dominance.csv", index=False)
    callability.to_csv(tables_dir / "marker_callability_by_bin.csv", index=False)
    callability_tests.to_csv(
        tables_dir / "marker_callability_tests.csv", index=False
    )
    diversity.to_csv(tables_dir / "diversity_site_bootstrap_summary.csv", index=False)
    turnover.to_csv(tables_dir / "turnover_site_bootstrap_summary.csv", index=False)
    stale_global = tables_dir / "global_composition_permutation_tests.csv"
    if stale_global.exists():
        stale_global.unlink()
    tests.to_csv(tables_dir / "global_composition_cluster_tests.csv", index=False)
    dispersion_test_table.to_csv(
        tables_dir / "composition_dispersion_cluster_tests.csv", index=False
    )
    dispersion_summary_table.to_csv(
        tables_dir / "composition_dispersion_by_period.csv", index=False
    )
    residual_diagnostic_table.to_csv(
        tables_dir / "cluster_model_residual_diagnostics.csv", index=False
    )
    residual_diagnostic_summary_table.to_csv(
        tables_dir / "cluster_model_diagnostic_summary.csv", index=False
    )
    sensitivity.to_csv(tables_dir / "turnover_sensitivity.csv", index=False)
    paired_draws.to_csv(
        tables_dir / "paired_male_turnover_bootstrap_draws.csv", index=False
    )
    family_draws.to_csv(
        tables_dir
        / "paired_male_y_resolution_family_bootstrap_draws.csv",
        index=False,
    )
    paired_summary.to_csv(
        tables_dir / "paired_male_turnover_bootstrap_summary.csv", index=False
    )
    paired_resolution_sensitivity.to_csv(
        tables_dir / "paired_male_y_resolution_sensitivity.csv", index=False
    )
    date_uncertainty_summary.to_csv(
        tables_dir / "date_uncertainty_summary.csv", index=False
    )

    figure_sampling(
        public_analysis, counts, figures_dir / "figure_1_sampling.png"
    )
    figure_composition(
        profiles[("mtDNA", "site")],
        mt_categories,
        "mtDNA",
        figures_dir / "figure_2_mtdna_composition.png",
    )
    figure_composition(
        profiles[("Y", "site")],
        y_categories,
        "Y chromosome",
        figures_dir / "figure_3_y_composition.png",
    )
    figure_diversity_turnover(
        diversity, turnover, figures_dir / "figure_4_diversity_turnover.png"
    )
    figure_sensitivity(
        sensitivity[
            ~sensitivity["analysis"].str.startswith("Leave out country:")
        ],
        figures_dir / "figure_5_sensitivity.png",
    )
    figure_database_coverage(
        db_coverage, figures_dir / "figure_6_database_coverage.png"
    )

    mt_calls = int((analysis["mt_call"] != "").sum())
    y_calls = int((analysis["y_call"] != "").sum())
    results = {
        "aadr_catalogue_unique_archaeological_individuals": int(len(canonical)),
        "primary_analysis_unique_individuals_3500BCE_to_1500CE": int(
            len(analysis)
        ),
        "primary_sites": int(analysis["locality"].nunique()),
        "primary_publication_labels": int(analysis["publication"].nunique()),
        "primary_mt_calls": mt_calls,
        "primary_y_calls_in_molecular_males": y_calls,
        "primary_molecular_males": int(
            analysis["molecular_sex"].astype(str).str.startswith("M").sum()
        ),
        "strict_qc_individuals": int(analysis["strict_qc"].sum()),
        "directly_dated_individuals": int(analysis["direct_date"].sum()),
        "reported_kin_components_le2d": int(
            analysis.groupby("kin_component_2d").size().gt(1).sum()
        ),
        "mt_l1_categories": mt_categories,
        "y_l1_categories": y_categories,
        "database_exact_crosswalk": [
            {
                "source": source,
                "match_status": status,
                "records": int(n),
            }
            for (source, status), n in crosswalk.groupby(
                ["source", "match_status"]
            ).size().items()
        ],
        "extended_legacy_mtdna_records_not_in_aadr": int(len(legacy_mt)),
        "extended_legacy_mtdna_valid_calls_not_in_aadr": int(
            legacy_mt["mt_hg"].map(valid_call).sum()
        ),
        "global_composition_tests": tests.to_dict("records"),
        "composition_dispersion_diagnostics": dispersion_test_table.to_dict(
            "records"
        ),
        "cluster_model_diagnostics": residual_diagnostic_summary_table.to_dict(
            "records"
        ),
        "paired_male_turnover_comparison": paired_summary.to_dict("records"),
        "paired_male_y_resolution_sensitivity": (
            paired_resolution_sensitivity.to_dict("records")
        ),
        "date_assignment_scenario_sensitivity": (
            date_uncertainty_summary.to_dict("records")
        ),
        "country_period_cells_meeting_threshold": [
            {
                "marker": marker,
                "cells": int(
                    g["meets_descriptive_frequency_threshold"].sum()
                ),
                "total_cells": int(len(g)),
            }
            for marker, g in adequacy.groupby("marker")
        ],
        "callability_tests": callability_tests.to_dict("records"),
        "important_estimand_warning": (
            "Estimates describe published individuals, equal-weighted sites "
            "within descriptive bins, or equal-weighted published "
            "site-period profiles in the regression; they are not population "
            "frequencies in the past."
        ),
    }
    (out / "results_summary.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    manifest = {
        "analysis_date": "2026-07-25",
        "code_revision_date": "2026-08-15",
        "random_seed": args.seed,
        "rng_streams": (
            "Stable SHA-256-named NumPy SeedSequence streams; each analysis "
            "is invariant to changes in unrelated random procedures"
        ),
        "bootstrap_replicates": args.bootstrap,
        "paired_bootstrap_replicates": args.paired_bootstrap,
        "site_cluster_bootstrap": {
            "cluster": "country + locality",
            "draw": (
                "one multinomial multiplicity per cluster, shared across all "
                "periods; paired analysis also shares it across both markers"
            ),
            "empty_period_handling": "reject complete draw and redraw",
            "marker_diagnostics": bootstrap_diagnostics,
            "paired_diagnostics": paired_bootstrap_diagnostics,
            "paired_y_resolution_diagnostics": {
                "broad_L1": paired_bootstrap_diagnostics,
                "AADR_ISOGG_prefix_family": family_bootstrap_diagnostics,
                "shared_cluster_draws_across_encodings": True,
                "warning": (
                    "Resolution comparison is an encoding sensitivity, not a "
                    "new demographic or phylogenetic test"
                ),
            },
        },
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
        "primary_inference": {
            "response": "Hellinger-transformed broad-lineage proportions",
            "unit": "country-locality-period profile",
            "reduced_model": "country",
            "full_model": "country + period",
            "resampling": (
                "Null-model residuals divided rowwise by sqrt(1-h_ii), "
                "with one shared Rademacher sign per country-locality cluster"
            ),
            "multiplicity": "Holm adjustment across mtDNA and Y",
        },
        "exploratory_within_locality_test": {
            "subset": "localities represented in more than one period",
            "model": "site fixed effects + period",
            "resampling": "period-label permutation within locality",
        },
        "haplogroup_harmonization": {
            "mtDNA_HV": "retained as HV, separate from H",
            "basal_Y_compounds": "Basal/unresolved, including CF, CT, F and IJK",
            "rare_pooling": "fewer than five primary calls pooled as Other",
            "Y_resolution_sensitivity": (
                "AADR ISOGG call prefix: first letter + first integer + "
                "immediately following branch letter; nomenclature sensitivity "
                "only, not a phylogenetic re-call"
            ),
        },
        "inputs": {
            str(args.aadr): input_hashes["aadr"],
            str(args.amtdb): input_hashes["amtdb"],
            str(args.aychr): input_hashes["aychr"],
        },
        "public_coordinate_policy": {
            "columns": ["latitude", "longitude"],
            "rounding_decimal_degrees": 1,
            "rationale": (
                "reduce archaeological site-location precision; coordinates "
                "are not used in statistical models"
            ),
        },
        "temporal_bins": {
            label: [int(start), int(end)]
            for label, start, end in zip(
                BIN_LABELS, BIN_EDGES[:-1], BIN_EDGES[1:]
            )
        },
        "source_code_sha256": sha256(Path(__file__)),
    }
    (out / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    print(json.dumps(results, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
