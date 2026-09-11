#!/usr/bin/env python3
"""Rebuild and verify the primary study input directly from frozen AADR."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from run_analysis import (
    EXPECTED_INPUT_SHA256,
    mask_coordinates,
    prepare_aadr_catalogues,
    sha256,
)


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_OUTPUT_SHA256 = "692a69cf38cc736f96ea5aa6b3b15024a9a49d50c03ae3e1470a8dc475504cc3"


def extract_aadr_input(aadr_path: Path, workdir: Path, output: Path) -> dict:
    """Verify the source first; expose a candidate only after exact hash match."""
    aadr_path, workdir, output = (
        aadr_path.resolve(), workdir.resolve(), output.resolve()
    )
    if workdir.is_relative_to(ROOT):
        raise ValueError("Source-labelled working catalogues require a workdir outside the repository")
    if output.is_relative_to(ROOT):
        raise ValueError("Output must be outside the repository")
    if output == aadr_path or output.exists():
        raise ValueError("Output must be a new file and must not replace the AADR source")
    if workdir.exists() and (not workdir.is_dir() or any(workdir.iterdir())):
        raise ValueError("Source workdir must be absent or empty")
    reserved = {"aadr_primary_analysis_catalogue.csv",
                "central_asia_analysis_input_v1.candidate.csv",
                "aadr_extraction_manifest.json"}
    if output.parent == workdir and output.name in reserved:
        raise ValueError("Output conflicts with a reserved extraction file")
    observed_source = sha256(aadr_path)
    if observed_source != EXPECTED_INPUT_SHA256["aadr"]:
        raise ValueError(
            "AADR SHA-256 mismatch: "
            f"expected {EXPECTED_INPUT_SHA256['aadr']}, observed {observed_source}"
        )
    aadr = pd.read_csv(aadr_path, sep="\t", dtype=str, keep_default_na=False)
    canonical, primary, _, mt_categories, y_categories = prepare_aadr_catalogues(aadr)
    workdir.mkdir(parents=True, exist_ok=True)
    catalogue = workdir / "aadr_primary_analysis_catalogue.csv"
    mask_coordinates(primary).to_csv(catalogue, index=False)
    candidate = workdir / "central_asia_analysis_input_v1.candidate.csv"
    subprocess.run(
        [sys.executable, str(ROOT / "analysis/build_public_analysis_input.py"),
         "--input", str(catalogue), "--output", str(candidate)],
        check=True,
        capture_output=True,
        text=True,
    )
    observed_output = sha256(candidate)
    if observed_output != EXPECTED_OUTPUT_SHA256:
        raise ValueError(
            "Rebuilt analytical input differs from frozen input: "
            f"expected {EXPECTED_OUTPUT_SHA256}, observed {observed_output}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation also protects against a destination appearing while
    # the source is being parsed. Never replace a source or release file.
    with output.open("xb") as destination, candidate.open("rb") as source:
        shutil.copyfileobj(source, destination)
    report = {
        "route": "verified_aadr_to_frozen_primary_input",
        "source_sha256": observed_source,
        "output_sha256": observed_output,
        "matches_frozen_input_byte_for_byte": True,
        "aadr_canonical_records": len(canonical),
        "primary_records": len(primary),
        "primary_sites": len(primary[["country", "locality"]].drop_duplicates()),
        "mt_calls": int(primary["mt_call"].ne("").sum()),
        "y_calls": int(primary["y_call"].ne("").sum()),
        "mt_l1_categories": mt_categories,
        "y_l1_categories": y_categories,
        "source_code_sha256": sha256(ROOT / "analysis/run_analysis.py"),
        "extractor_code_sha256": sha256(Path(__file__).resolve()),
        "builder_code_sha256": sha256(ROOT / "analysis/build_public_analysis_input.py"),
        "cross_database_audits_recomputed": False,
        "scope": "Rebuilds AADR primary study input; AmtDB/aYChr legacy cross-database audits require the separate full route.",
    }
    (workdir / "aadr_extraction_manifest.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aadr", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = extract_aadr_input(args.aadr, args.workdir, args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
