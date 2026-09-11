#!/usr/bin/env python3
"""Apply the documented public-coordinate policy to derived release files."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from run_analysis import figure_sampling, mask_coordinates


FILES = (
    "aadr_central_asia_unique_individual_catalogue.csv",
    "aadr_primary_analysis_catalogue.csv",
    "extended_legacy_mtdna_not_in_aadr.csv",
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
    parser.add_argument("--digits", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tables = args.analysis_output / "tables"
    audit = {}
    for name in FILES:
        path = tables / name
        before = sha256(path)
        data = pd.read_csv(path, keep_default_na=False)
        sanitized = mask_coordinates(data, args.digits)
        sanitized.to_csv(path, index=False)
        audit[name] = {
            "private_input_sha256": before,
            "public_output_sha256": sha256(path),
            "rows": len(sanitized),
        }

    primary = pd.read_csv(
        tables / "aadr_primary_analysis_catalogue.csv",
        keep_default_na=False,
    )
    # Preserve blank public coordinates in CSV, but convert the two columns to
    # numbers for plotting. Coordinates are not used by the statistical models.
    for column in ("latitude", "longitude"):
        primary[column] = pd.to_numeric(primary[column], errors="coerce")
    counts = pd.read_csv(tables / "counts_country_by_bin.csv")
    figure_sampling(
        primary,
        counts,
        args.analysis_output / "figures" / "figure_1_sampling.png",
    )

    manifest_path = args.analysis_output / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["public_coordinate_policy"] = {
        "columns": ["latitude", "longitude"],
        "rounding_decimal_degrees": args.digits,
        "rationale": (
            "reduce archaeological site-location precision; coordinates are "
            "not used in statistical models"
        ),
        "files": audit,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
