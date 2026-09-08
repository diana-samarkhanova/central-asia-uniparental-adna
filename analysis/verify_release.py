#!/usr/bin/env python3
"""Fail closed unless the repository satisfies the public-release contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path, PurePosixPath

# The verifier must not create a cache file that its own release policy rejects.
sys.dont_write_bytecode = True

from generate_checksums import (
    CACHE_DIRECTORY_NAMES,
    CACHE_FILE_NAMES,
    CACHE_SUFFIXES,
    DEFAULT_CHECKSUM_NAME,
    release_files,
    sha256,
)


DEFAULT_RESULTS = Path("results/aadr-v66p1_2026-07-25")
FORBIDDEN_GRANULAR_TABLES = (
    "aadr_central_asia_unique_individual_catalogue.csv",
    "aadr_deduplication_audit.csv",
    "aadr_primary_analysis_catalogue.csv",
    "cross_database_exact_id_audit.csv",
    "extended_legacy_mtdna_not_in_aadr.csv",
    "site_profiles_mtdna.csv",
    "site_profiles_y.csv",
    "cluster_model_residual_diagnostics.csv",
    "dispersion_profiles_mtdna.csv",
    "dispersion_profiles_y.csv",
)
INDIVIDUAL_LEVEL_HEADER_FIELDS = frozenset(
    {
        "individual_id",
        "genetic_id",
        "persistent_genetic_id",
        "source_id",
        "aadr_individual_ids",
        "mt_call",
        "y_call",
        "mt_hg",
        "mt_haplogroup",
        "y_haplogroup",
    }
)
EXPECTED_COUNTS = {
    "catalogue": 501,
    "primary": 489,
    "mt_calls": 438,
    "y_calls": 229,
    "paired_bootstrap": 50000,
}
PUBLIC_INPUT = Path("data/derived/central_asia_analysis_input_v1.csv")
EXTENSION_PROFILE_VECTORS = Path(
    "analysis_extensions/statistical_extensions_v4/results/"
    "stable_profile_vectors.csv"
)
PUBLIC_INPUT_COLUMNS = (
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
)
RESTRICTED_IDENTIFIERS = tuple(f"CKZ{index:03d}" for index in range(1, 5))

FORBIDDEN_DOCUMENT_SUFFIXES = {
    ".doc",
    ".docx",
    ".odt",
    ".pdf",
    ".ppt",
    ".pptx",
    ".rtf",
    ".xls",
    ".xlsx",
}
FORBIDDEN_ARCHIVE_SUFFIXES = {
    ".7z",
    ".bz2",
    ".gz",
    ".rar",
    ".tar",
    ".tgz",
    ".xz",
    ".zip",
}
FORBIDDEN_RAW_SUFFIXES = {
    ".anno",
    ".bai",
    ".bam",
    ".bcf",
    ".cram",
    ".crai",
    ".fa",
    ".fasta",
    ".fastq",
    ".fq",
    ".geno",
    ".sam",
    ".vcf",
}
FORBIDDEN_TEMP_SUFFIXES = {".bak", ".orig", ".swp", ".tmp"}
FORBIDDEN_RAW_FILENAMES = {
    "v66.p1_2M.aadr.PUB.anno",
    "amtdb_v1.009_metadata.csv",
    "a-YChr-DB_V5.xlsx",
}
CHECKSUM_LINE = re.compile(r"^([0-9a-fA-F]{64}) ([ *])(.+)$")
SESSION_PATH_PATTERNS = (
    re.compile(r"/(?:workspace)/scratch/"),
    re.compile(r"/(?:root)/[.]codex/(?:sessions)/"),
    re.compile(r"sandbox:/(?:workspace)/"),
    re.compile(r"/tmp/[0-9a-f]{12,}(?:[-_/])", re.IGNORECASE),
    re.compile(r"[A-Za-z]:[\\/](?:Users|workspace)[\\/].*?[\\/]scratch[\\/]", re.IGNORECASE),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify checksums, provenance, counts and release sanitization."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root (default: parent of the analysis directory).",
    )
    parser.add_argument(
        "--checksums",
        type=Path,
        default=Path(DEFAULT_CHECKSUM_NAME),
        help="Checksum manifest, relative to --root unless absolute.",
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_RESULTS,
        help="Released result directory, relative to --root unless absolute.",
    )
    return parser.parse_args()


def resolve_inside(root: Path, value: Path, label: str) -> Path:
    path = value if value.is_absolute() else root / value
    path = path.resolve()
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must be inside the repository root: {path}") from error
    return path


def verify_checksums(root: Path, checksum_path: Path, failures: list[str]) -> int:
    if not checksum_path.is_file():
        failures.append(f"missing checksum manifest: {checksum_path.relative_to(root)}")
        return 0

    declared: dict[str, str] = {}
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        match = CHECKSUM_LINE.fullmatch(line)
        if not match:
            failures.append(
                f"invalid SHA256SUMS line {line_number}; expected '<hash>  <path>'"
            )
            continue
        digest, _, name = match.groups()
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or name.startswith("./"):
            failures.append(f"unsafe SHA256SUMS path on line {line_number}: {name!r}")
            continue
        canonical_name = relative.as_posix()
        if canonical_name in declared:
            failures.append(f"duplicate SHA256SUMS entry: {canonical_name}")
            continue
        declared[canonical_name] = digest.lower()

    try:
        current_files = release_files(root, checksum_path)
    except ValueError as error:
        failures.append(str(error))
        return len(declared)
    current = {
        path.relative_to(root).as_posix(): path for path in current_files
    }
    for name in sorted(set(current) - set(declared)):
        failures.append(f"file missing from SHA256SUMS: {name}")
    for name in sorted(set(declared) - set(current)):
        failures.append(f"SHA256SUMS lists a missing or excluded file: {name}")
    for name in sorted(set(current) & set(declared)):
        observed = sha256(current[name])
        if observed != declared[name]:
            failures.append(
                f"checksum mismatch: {name} (expected {declared[name]}, observed {observed})"
            )
    return len(declared)


def verify_forbidden_files(root: Path, failures: list[str]) -> None:
    for path in root.rglob("*"):
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue
        if ".git" in relative.parts:
            continue
        relative_name = relative.as_posix()
        if path.is_symlink():
            failures.append(f"symlink is not allowed in the release: {relative_name}")
            continue
        if any(part in CACHE_DIRECTORY_NAMES for part in relative.parts):
            if path.is_dir():
                failures.append(f"cache directory is present: {relative_name}")
            continue
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if path.name in FORBIDDEN_GRANULAR_TABLES:
            failures.append(f"granular table is forbidden in GitHub release: {relative_name}")
        if path.name in CACHE_FILE_NAMES or suffix in CACHE_SUFFIXES:
            failures.append(f"cache file is present: {relative_name}")
        if "raw" in {part.lower() for part in relative.parts[:-1]}:
            failures.append(f"raw-data directory contains a file: {relative_name}")
        if path.name in FORBIDDEN_RAW_FILENAMES or suffix in FORBIDDEN_RAW_SUFFIXES:
            failures.append(f"raw scientific data file is present: {relative_name}")
        if suffix in FORBIDDEN_DOCUMENT_SUFFIXES:
            failures.append(f"rendered/office document is present: {relative_name}")
        if suffix in FORBIDDEN_ARCHIVE_SUFFIXES:
            failures.append(f"archive is present: {relative_name}")
        if suffix in FORBIDDEN_TEMP_SUFFIXES or path.name.endswith("~"):
            failures.append(f"temporary/backup file is present: {relative_name}")


def verify_sensitive_content(
    root: Path, checksum_path: Path, failures: list[str]
) -> None:
    try:
        files = release_files(root, checksum_path)
    except ValueError:
        return
    restricted = {
        identifier: identifier.encode("ascii")
        for identifier in RESTRICTED_IDENTIFIERS
    }
    for path in files:
        payload = path.read_bytes()
        relative = path.relative_to(root).as_posix()
        for identifier, encoded in restricted.items():
            if encoded in payload:
                failures.append(
                    f"restricted person identifier {identifier} is present in {relative}"
                )
        text = payload.decode("utf-8", errors="ignore")
        for pattern in SESSION_PATH_PATTERNS:
            match = pattern.search(text)
            if match:
                failures.append(
                    f"scratch/session absolute path is present in {relative}: "
                    f"{match.group(0)!r}"
                )


def verify_aggregate_csv_headers(root: Path, failures: list[str]) -> None:
    """Reject CSV files that expose individual-level identifiers or marker calls."""
    for path in root.rglob("*.csv"):
        if ".git" in path.parts:
            continue
        try:
            with path.open("r", encoding="utf-8-sig", newline="") as stream:
                header = next(csv.reader(stream), [])
        except (OSError, csv.Error, UnicodeError) as error:
            failures.append(f"cannot inspect CSV header {path}: {error}")
            continue
        prohibited = sorted(INDIVIDUAL_LEVEL_HEADER_FIELDS.intersection(header))
        if prohibited:
            failures.append(
                "individual-level CSV header field(s) "
                f"{prohibited} in {path.relative_to(root).as_posix()}"
            )


def verify_public_analysis_input(root: Path, failures: list[str]) -> None:
    """Validate the sole schema-locked, project-coded individual-level input."""
    path = root / PUBLIC_INPUT
    rows = read_catalogue_without_coordinates(path, failures)
    if not rows:
        return
    observed_columns = tuple(rows[0])
    if observed_columns != PUBLIC_INPUT_COLUMNS:
        failures.append(
            "public analytical-input schema mismatch: "
            f"{observed_columns!r}"
        )
        return
    if len(rows) != EXPECTED_COUNTS["primary"]:
        failures.append(
            f"public analytical input has {len(rows)} rows, expected 489"
        )
    record_ids = [row["project_record_id"] for row in rows]
    if len(set(record_ids)) != len(record_ids) or any(
        re.fullmatch(r"CAU[0-9]{6}", value) is None for value in record_ids
    ):
        failures.append("project_record_id values are not unique CAU###### codes")
    if any(
        re.fullmatch(r"SITE[0-9]{4}", row["site_key"]) is None for row in rows
    ):
        failures.append("site_key values are not project SITE#### codes")
    if any(
        re.fullmatch(r"STUDY[0-9]{3}", row["study_key"]) is None for row in rows
    ):
        failures.append("study_key values are not project STUDY### codes")
    sites = {(row["country"], row["site_key"]) for row in rows}
    if len(sites) != 136:
        failures.append(f"public analytical input has {len(sites)} sites, expected 136")
    if {row["molecular_sex"] for row in rows} - {"M", "F", "U"}:
        failures.append("public analytical input has an unexpected molecular-sex code")
    boolean_fields = (
        "mt_called",
        "y_called",
        "strict_qc",
        "population_outlier",
        "direct_date",
        "kin_representative_mt",
        "kin_representative_y",
    )
    if any(
        row[field].lower() not in {"true", "false"}
        for row in rows
        for field in boolean_fields
    ):
        failures.append("public analytical input contains a non-Boolean flag")
    mt_calls = sum(row["mt_called"].lower() == "true" for row in rows)
    y_calls = sum(row["y_called"].lower() == "true" for row in rows)
    if mt_calls != EXPECTED_COUNTS["mt_calls"] or y_calls != EXPECTED_COUNTS["y_calls"]:
        failures.append(
            f"public analytical input has {mt_calls}/{y_calls} mtDNA/Y calls, "
            "expected 438/229"
        )
    if any(
        row["y_called"].lower() == "true" and row["molecular_sex"] != "M"
        for row in rows
    ):
        failures.append("a released Y call is assigned to a non-M row")
    if any(
        (row["latitude_0_1deg"] == "")
        != (row["longitude_0_1deg"] == "")
        for row in rows
    ):
        failures.append(
            "public analytical input has an unpaired missing coordinate"
        )
    for coordinate in ("latitude_0_1deg", "longitude_0_1deg"):
        for row in rows:
            if row[coordinate] == "":
                continue
            try:
                value = float(row[coordinate])
            except ValueError:
                failures.append(f"non-numeric {coordinate} in public analytical input")
                break
            if not math.isclose(value, round(value, 1), abs_tol=1e-12):
                failures.append(
                    f"{coordinate} is not rounded to 0.1 degree in public analytical input"
                )
                break


def verify_project_coded_extension_profiles(
    root: Path, failures: list[str]
) -> None:
    rows = read_catalogue_without_coordinates(
        root / EXTENSION_PROFILE_VECTORS, failures
    )
    if not rows:
        return
    if any(
        re.fullmatch(r"SITE[0-9]{4}", row.get("site", "")) is None
        for row in rows
    ):
        failures.append("extension profile vectors contain a non-project site key")
    publications = {
        token
        for row in rows
        for token in row.get("publication", "").split("|")
        if token
    }
    if any(re.fullmatch(r"STUDY[0-9]{3}", value) is None for value in publications):
        failures.append(
            "extension profile vectors contain a non-project study key"
        )


def load_json(path: Path, failures: list[str]) -> dict[str, object]:
    if not path.is_file():
        failures.append(f"missing JSON file: {path}")
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        failures.append(f"invalid JSON file {path}: {error}")
        return {}
    if not isinstance(value, dict):
        failures.append(f"JSON root must be an object: {path}")
        return {}
    return value


def integer_value(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric) or not numeric.is_integer():
        return None
    return int(numeric)


def verify_manifest_hash(
    manifest: dict[str, object], key: str, script: Path, label: str, failures: list[str]
) -> None:
    expected = manifest.get(key)
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
        failures.append(f"{label} has no valid {key}")
        return
    observed = sha256(script)
    if expected != observed:
        failures.append(
            f"{label} {key} mismatch for {script.name}: expected {expected}, observed {observed}"
        )


def verify_manifests_and_counts(
    root: Path, results: Path, failures: list[str]
) -> None:
    summary = load_json(results / "results_summary.json", failures)
    analysis_manifest = load_json(results / "analysis_manifest.json", failures)
    sensitivity_manifest = load_json(
        results / "global_sensitivity_manifest.json", failures
    )

    summary_fields = {
        "aadr_catalogue_unique_archaeological_individuals": "catalogue",
        "primary_analysis_unique_individuals_3500BCE_to_1500CE": "primary",
        "primary_mt_calls": "mt_calls",
        "primary_y_calls_in_molecular_males": "y_calls",
    }
    for field, expected_name in summary_fields.items():
        observed = integer_value(summary.get(field))
        expected = EXPECTED_COUNTS[expected_name]
        if observed != expected:
            failures.append(
                f"results_summary.json {field} is {observed!r}, expected {expected}"
            )
    if integer_value(summary.get("primary_sites")) != 136:
        failures.append("results_summary.json primary_sites must equal 136")

    run_analysis = root / "analysis" / "run_analysis.py"
    run_sensitivity = root / "analysis" / "run_global_sensitivities.py"
    recompute = root / "analysis" / "recompute_from_catalogue.py"
    verify_manifest_hash(
        analysis_manifest,
        "source_code_sha256",
        run_analysis,
        "analysis_manifest.json",
        failures,
    )
    verify_manifest_hash(
        sensitivity_manifest,
        "source_code_sha256",
        run_sensitivity,
        "global_sensitivity_manifest.json",
        failures,
    )
    verify_manifest_hash(
        sensitivity_manifest,
        "primary_analysis_code_sha256",
        run_analysis,
        "global_sensitivity_manifest.json",
        failures,
    )
    sensitivity_input = sensitivity_manifest.get("analytical_input")
    public_input_path = root / PUBLIC_INPUT
    if not isinstance(sensitivity_input, dict):
        failures.append(
            "global_sensitivity_manifest.json must record its analytical input"
        )
    else:
        if sensitivity_input.get("file") != PUBLIC_INPUT.as_posix():
            failures.append(
                "global_sensitivity_manifest.json analytical_input.file must "
                f"be {PUBLIC_INPUT.as_posix()}"
            )
        if sensitivity_input.get("public_project_coded_input") is not True:
            failures.append(
                "global_sensitivity_manifest.json must record use of the "
                "project-coded public input"
            )
        if (
            public_input_path.is_file()
            and sensitivity_input.get("sha256") != sha256(public_input_path)
        ):
            failures.append(
                "global_sensitivity_manifest.json analytical-input SHA-256 "
                "does not match the released file"
            )

    recomputed = analysis_manifest.get("recomputed_from")
    if isinstance(recomputed, dict):
        verify_manifest_hash(
            recomputed,
            "script_sha256",
            recompute,
            "analysis_manifest.json recomputed_from",
            failures,
        )
        if recomputed.get("script") != "analysis/recompute_from_catalogue.py":
            failures.append(
                "analysis_manifest.json recomputed_from.script must be the "
                "repository-relative path analysis/recompute_from_catalogue.py"
            )
        if recomputed.get("file") != PUBLIC_INPUT.as_posix():
            failures.append(
                "analysis_manifest.json recomputed_from.file must be "
                f"{PUBLIC_INPUT.as_posix()}"
            )
        expected_input_hash = recomputed.get("sha256")
        if not isinstance(expected_input_hash, str) or (
            public_input_path.is_file()
            and expected_input_hash != sha256(public_input_path)
        ):
            failures.append(
                "analysis_manifest.json public analytical-input SHA-256 "
                "does not match the released file"
            )

    paired_expected = EXPECTED_COUNTS["paired_bootstrap"]
    manifest_paired = integer_value(
        analysis_manifest.get("paired_bootstrap_replicates")
    )
    if manifest_paired != paired_expected:
        failures.append(
            "analysis_manifest.json paired_bootstrap_replicates is "
            f"{manifest_paired!r}, expected {paired_expected}"
        )
    public_section = analysis_manifest.get("public_project_coded_input")
    if not isinstance(public_section, dict) or not public_section.get(
        "used_for_this_run"
    ):
        failures.append(
            "analysis_manifest.json must record recomputation from the "
            "project-coded public analytical input"
        )
    surrogate_policy = analysis_manifest.get("paired_surrogate_probability_policy")
    if not isinstance(surrogate_policy, dict) or not surrogate_policy.get(
        "removed_from_user_facing_outputs"
    ):
        failures.append(
            "analysis_manifest.json must record removal of the surrogate "
            "bootstrap probability from user-facing outputs"
        )
    bootstrap_section = analysis_manifest.get("site_cluster_bootstrap")
    paired_diagnostics = (
        bootstrap_section.get("paired_diagnostics")
        if isinstance(bootstrap_section, dict)
        else None
    )
    accepted = (
        integer_value(paired_diagnostics.get("accepted_replicates"))
        if isinstance(paired_diagnostics, dict)
        else None
    )
    if accepted != paired_expected:
        failures.append(
            "analysis_manifest.json paired accepted_replicates is "
            f"{accepted!r}, expected {paired_expected}"
        )

    comparisons = summary.get("paired_male_turnover_comparison")
    summary_paired = None
    if (
        isinstance(comparisons, list)
        and len(comparisons) == 1
        and isinstance(comparisons[0], dict)
    ):
        summary_paired = integer_value(comparisons[0].get("bootstrap_replicates"))
    if summary_paired != paired_expected:
        failures.append(
            "results_summary.json paired bootstrap count is "
            f"{summary_paired!r}, expected {paired_expected}"
        )

    paired_path = results / "tables" / "paired_male_turnover_bootstrap_summary.csv"
    paired_rows = read_catalogue_without_coordinates(paired_path, failures)
    csv_paired = (
        integer_value(paired_rows[0].get("bootstrap_replicates"))
        if len(paired_rows) == 1
        else None
    )
    if csv_paired != paired_expected:
        failures.append(
            f"paired bootstrap summary count is {csv_paired!r}, expected {paired_expected}"
        )


def read_catalogue_without_coordinates(
    path: Path, failures: list[str]
) -> list[dict[str, str]]:
    if not path.is_file():
        failures.append(f"missing CSV file: {path}")
        return []
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            return list(csv.DictReader(stream))
    except (OSError, csv.Error, UnicodeError) as error:
        failures.append(f"cannot read CSV file {path}: {error}")
        return []


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"Repository root does not exist: {root}")
    try:
        checksum_path = resolve_inside(root, args.checksums, "checksum manifest")
        results = resolve_inside(root, args.results, "results directory")
    except ValueError as error:
        raise SystemExit(str(error)) from error

    failures: list[str] = []
    checksum_entries = verify_checksums(root, checksum_path, failures)
    verify_forbidden_files(root, failures)
    verify_sensitive_content(root, checksum_path, failures)
    verify_aggregate_csv_headers(root, failures)
    verify_public_analysis_input(root, failures)
    verify_project_coded_extension_profiles(root, failures)
    verify_manifests_and_counts(root, results, failures)

    if failures:
        unique_failures = list(dict.fromkeys(failures))
        print(
            f"Release verification failed with {len(unique_failures)} issue(s):",
            file=sys.stderr,
        )
        for failure in unique_failures:
            print(f"- {failure}", file=sys.stderr)
        raise SystemExit(1)
    print(
        "Project-coded-input public release verification passed: "
        f"{checksum_entries} checksums; 501/489 individuals; "
        "438 mtDNA calls; 229 Y calls; 50,000 paired bootstrap replicates."
    )


if __name__ == "__main__":
    main()
