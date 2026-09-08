#!/usr/bin/env python3
"""Build the project-coded AADR-derived input used by the public workflow.

The output deliberately omits source identifiers, archaeological locality
strings, skeletal fields, exact coordinates and terminal haplogroup calls.  It
is project-coded rather than anonymous: combinations of public AADR-derived
attributes may still be linkable to the source catalogue.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from run_analysis import (
    add_y_resolution_sensitivity_encoding,
    apply_site_locality_aliases,
    unrelated_subset,
)


OUTPUT_COLUMNS = [
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
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def coded_mapping(values: list[object], prefix: str, width: int) -> dict[object, str]:
    return {
        value: f"{prefix}{index:0{width}d}"
        for index, value in enumerate(sorted(set(values)), start=1)
    }


def main() -> None:
    args = parse_args()
    source = pd.read_csv(args.input, keep_default_na=False)
    source = apply_site_locality_aliases(source)
    source, _, _ = add_y_resolution_sensitivity_encoding(source, 5)

    required = {
        "individual_id",
        "country",
        "locality",
        "analysis_bin",
        "publication",
        "molecular_sex",
        "mt_call",
        "y_call",
        "mt_l1_pooled",
        "y_l1_pooled",
        "y_isogg_prefix_family_pooled",
        "strict_qc",
        "population_outlier",
        "direct_date",
        "date_bp",
        "date_sd_bp",
        "kin_component_2d",
        "latitude",
        "longitude",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError(f"Source catalogue lacks required columns: {missing}")

    # Preserve the frozen catalogue order so named random streams allocate
    # draws identically in downstream procedures. Project codes, rather than
    # source identifiers, make that order stable in the public input.
    source = source.copy()
    source["project_record_id"] = [
        f"CAU{index:06d}" for index in range(1, len(source) + 1)
    ]

    site_values = list(zip(source["country"], source["locality"]))
    site_map = coded_mapping(site_values, "SITE", 4)
    study_map = coded_mapping(source["publication"].tolist(), "STUDY", 3)
    mt_representative_indices = set(unrelated_subset(source, "mt_call").index)
    y_representative_indices = set(unrelated_subset(source, "y_call").index)

    output = pd.DataFrame(
        {
            "project_record_id": source["project_record_id"],
            "country": source["country"],
            "site_key": [site_map[value] for value in site_values],
            "analysis_bin": source["analysis_bin"],
            "study_key": source["publication"].map(study_map),
            "molecular_sex": source["molecular_sex"],
            "mt_category": source["mt_l1_pooled"],
            "y_category": source["y_l1_pooled"],
            "y_prefix_category": source["y_isogg_prefix_family_pooled"],
            "mt_called": source["mt_l1_pooled"].ne(""),
            "y_called": source["y_l1_pooled"].ne(""),
            "strict_qc": source["strict_qc"],
            "population_outlier": source["population_outlier"],
            "direct_date": source["direct_date"],
            "date_bp": pd.to_numeric(source["date_bp"], errors="raise"),
            "date_sd_bp": pd.to_numeric(
                source["date_sd_bp"], errors="coerce"
            ).fillna(0),
            "kin_representative_mt": source.index.isin(
                mt_representative_indices
            ),
            "kin_representative_y": source.index.isin(
                y_representative_indices
            ),
            "latitude_0_1deg": pd.to_numeric(
                source["latitude"], errors="coerce"
            ).round(1),
            "longitude_0_1deg": pd.to_numeric(
                source["longitude"], errors="coerce"
            ).round(1),
        }
    )[OUTPUT_COLUMNS]

    if len(output) != 489 or not output["project_record_id"].is_unique:
        raise AssertionError("Public input must contain 489 uniquely coded rows")
    if output[["country", "site_key"]].drop_duplicates().shape[0] != 136:
        raise AssertionError("Public input must contain 136 normalized sites")
    if int(output["mt_called"].sum()) != 438:
        raise AssertionError("Public input must contain 438 mtDNA calls")
    if int(output["y_called"].sum()) != 229:
        raise AssertionError("Public input must contain 229 Y calls")
    if set(output["molecular_sex"]) - {"M", "F", "U"}:
        raise AssertionError("Unexpected molecular-sex code")
    if not output.loc[output["y_called"], "molecular_sex"].eq("M").all():
        raise AssertionError("Every released Y call must belong to an M-coded row")
    payload = output.to_csv(index=False)
    if re.search(r"CKZ00[1-4]", payload):
        raise AssertionError("Restricted post-freeze identifier in public input")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")
    print(
        f"Wrote {len(output)} rows, "
        f"{output[['country', 'site_key']].drop_duplicates().shape[0]} sites, "
        f"{int(output['mt_called'].sum())}/{int(output['y_called'].sum())} calls"
    )


if __name__ == "__main__":
    main()
