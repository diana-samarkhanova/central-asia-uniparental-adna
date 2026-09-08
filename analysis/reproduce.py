#!/usr/bin/env python3
"""Run the frozen analysis into an empty, external directory with an audit trail.

The derived route inherits only explicitly listed source-audit metadata and
category definitions. Statistical tables are regenerated, never seeded by
copying an old results directory. The raw route requires all three exact
inputs and also rebuilds the project-coded input used by the v4 extensions.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results/aadr-v66p1_2026-07-25"
EXTENSION = ROOT / "analysis_extensions/statistical_extensions_v4"
PUBLIC_INPUT = ROOT / "data/derived/central_asia_analysis_input_v1.csv"
FROZEN_INPUT_SHA256 = "692a69cf38cc736f96ea5aa6b3b15024a9a49d50c03ae3e1470a8dc475504cc3"
FROZEN_OPTIONS = ["--bootstrap", "2000", "--paired-bootstrap", "50000",
                  "--permutations", "9999", "--callability-resamples", "99999",
                  "--date-draws", "5000", "--seed", "20260725"]
INHERITED_TABLES = {"database_coverage_by_country.csv"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def environment_check() -> dict:
    packages = {}
    errors = []
    if sys.version_info[:2] != (3, 12):
        errors.append("Python 3.12 is required")
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, expected = line.split("==")
        try:
            installed = version(name)
        except PackageNotFoundError:
            installed = None
        packages[name] = {"expected": expected, "installed": installed}
        if installed != expected:
            errors.append(f"{name}: expected {expected}, found {installed}")
    if errors:
        raise ValueError("Frozen environment mismatch: " + "; ".join(errors)
                         + ". Install requirements.txt before reproducing.")
    return {"python": platform.python_version(), "platform": platform.platform(),
            "packages": packages}


def verify_analytical_input(path: Path) -> str:
    observed = sha256(path)
    if observed != FROZEN_INPUT_SHA256:
        raise ValueError("Analytical input SHA-256 differs from the frozen v4 input")
    return observed


def reserve_output(output: Path, repo_root: Path = ROOT) -> None:
    if output.is_relative_to(repo_root.resolve()):
        raise ValueError("Use an output directory outside the repository; generated "
                         "source rows and results must not replace release files.")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output must be absent or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)


def compare_csv(expected: Path, observed: Path) -> dict:
    """Compare every row/column; tolerate only floating-point roundoff."""
    import numpy as np
    import pandas as pd
    if not observed.is_file():
        return {"file": expected.name, "status": "missing"}
    left = pd.read_csv(expected, dtype=str, keep_default_na=False)
    right = pd.read_csv(observed, dtype=str, keep_default_na=False)
    if left.shape != right.shape or list(left) != list(right):
        return {"file": expected.name, "status": "schema_mismatch",
                "expected_shape": list(left.shape), "observed_shape": list(right.shape)}
    differences = []
    max_error = 0.0
    for name in left:
        a, b = left[name], right[name]
        try:
            av = pd.to_numeric(a.mask(a.eq(""), float("nan")), errors="raise").to_numpy(float)
            bv = pd.to_numeric(b.mask(b.eq(""), float("nan")), errors="raise").to_numpy(float)
        except (ValueError, TypeError):
            if not a.equals(b):
                differences.append(name)
            continue
        finite = np.isfinite(av) & np.isfinite(bv)
        if finite.any():
            max_error = max(max_error, float(np.max(np.abs(av[finite] - bv[finite]))))
        if not np.isclose(av, bv, rtol=1e-9, atol=1e-12, equal_nan=True).all():
            differences.append(name)
    return {"file": expected.name, "status": "match" if not differences else "different",
            "rows": len(left), "different_columns": differences,
            "maximum_absolute_numeric_difference": max_error,
            "byte_identical": sha256(expected) == sha256(observed)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["derived", "aadr", "raw"], required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-root", type=Path, default=ROOT,
                        help="Root containing the data/raw paths in SOURCES.tsv")
    parser.add_argument("--fetch", action="store_true",
                        help="AADR/raw mode: retrieve exact inputs before analysis")
    parser.add_argument("--amtdb", type=Path,
                        help="Raw mode with --fetch: import original hash-matching CSV")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output.expanduser().resolve()
    try:
        if args.mode == "derived" and (args.fetch or args.amtdb):
            raise ValueError("--fetch applies to --mode aadr/raw; --amtdb only to raw")
        if args.mode != "raw" and args.amtdb:
            raise ValueError("--amtdb applies only to --mode raw")
        if args.amtdb and not args.fetch:
            raise ValueError("Use --fetch with --amtdb to verify and import the CSV")
        reserve_output(output)
    except ValueError as error:
        print(error, file=sys.stderr)
        return 2
    report = {"schema_version": "1.0", "mode": args.mode, "status": "running",
              "started_at": datetime.now(timezone.utc).isoformat(), "stages": [],
              "inherited_fixtures": [], "comparison_tolerance": {"rtol": 1e-9, "atol": 1e-12},
              "raw_sources_reproduced": False, "primary_aadr_input_rebuilt": False}
    report_path = output / "reproduction_report.json"
    logs = output / "logs"
    logs.mkdir()
    started = time.monotonic()

    def save() -> None:
        temporary = report_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
        temporary.replace(report_path)

    def stage(name: str, arguments: list[str]) -> None:
        entry = {"name": name, "arguments": arguments, "status": "running",
                 "log": f"logs/{name}.log"}
        report["stages"].append(entry)
        save()
        print(f"START {name}", flush=True)
        env = os.environ.copy()
        env.update({"MPLCONFIGDIR": str(output / "matplotlib_cache"),
                    "MPLBACKEND": "Agg", "PYTHONDONTWRITEBYTECODE": "1",
                    "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1"})
        begin = time.monotonic()
        with (logs / f"{name}.log").open("w") as log:
            result = subprocess.run([sys.executable, *arguments], cwd=ROOT, env=env,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, check=False)
        entry.update({"returncode": result.returncode,
                      "seconds": round(time.monotonic() - begin, 3),
                      "status": "passed" if result.returncode == 0 else "failed"})
        save()
        if result.returncode:
            raise RuntimeError(f"{name} failed; see {entry['log']}")
        print(f"PASS {name} ({entry['seconds']} s)", flush=True)

    try:
        report["environment"] = environment_check()
        report["frozen_analytical_input_sha256"] = verify_analytical_input(PUBLIC_INPUT)
        result_dir = output / "analysis"
        result_dir.mkdir()
        tables = result_dir / "tables"
        tables.mkdir()
        extension_input = PUBLIC_INPUT
        analytical_input = PUBLIC_INPUT
        if args.mode == "raw":
            source_root = args.source_root.expanduser().resolve()
            if args.fetch:
                arguments = ["analysis/fetch_sources.py", "--output-root", str(source_root),
                             "--report", str(output / "source_acquisition.json")]
                if args.amtdb:
                    arguments.extend(["--local", "amtdb=" + str(args.amtdb.resolve())])
                stage("fetch_sources", arguments)
            import csv
            with (ROOT / "data/SOURCES.tsv").open(newline="") as stream:
                registry = list(csv.DictReader(stream, delimiter="\t"))
            names = {"AADR": "aadr", "AmtDB": "amtdb", "aYChr-DB": "aychr"}
            source_options = []
            for row in registry:
                source_options.extend(["--" + names[row["resource"]],
                                       str(source_root / row["expected_path"])])
            stage("verify_sources", ["analysis/verify_inputs.py", *source_options])
            stage("raw_analysis", ["analysis/run_analysis.py", *source_options,
                                    "--outdir", str(result_dir), *FROZEN_OPTIONS])
            stage("global_sensitivities", ["analysis/run_global_sensitivities.py",
                  "--analysis-output", str(result_dir), "--permutations", "1999",
                  "--seed", "20260726"])
            stage("full_output_checks", ["analysis/test_full_outputs.py", str(result_dir)])
            extension_input = output / "rebuilt_analysis_input.csv"
            stage("rebuild_analysis_input", ["analysis/build_public_analysis_input.py",
                  "--input", str(tables / "aadr_primary_analysis_catalogue.csv"),
                  "--output", str(extension_input)])
            if sha256(extension_input) != sha256(PUBLIC_INPUT):
                raise RuntimeError("Rebuilt analytical input differs from frozen input")
            report["raw_sources_reproduced"] = True
            report["primary_aadr_input_rebuilt"] = True
            report["not_reproduced"] = ["literature search and source-scope audits",
                                        "separate wild-bootstrap calibration simulation"]
        else:
            if args.mode == "aadr":
                source_root = args.source_root.expanduser().resolve()
                if args.fetch:
                    stage("fetch_aadr", ["analysis/fetch_sources.py", "--resource", "aadr",
                          "--output-root", str(source_root), "--report",
                          str(output / "source_acquisition.json")])
                analytical_input = output / "rebuilt_analysis_input.csv"
                stage("extract_aadr_input", ["analysis/extract_aadr_input.py", "--aadr",
                      str(source_root / "data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno"),
                      "--workdir", str(output / "source_extraction"),
                      "--output", str(analytical_input)])
                if sha256(analytical_input) != sha256(PUBLIC_INPUT):
                    raise RuntimeError("Rebuilt AADR input differs from frozen input")
                extension_input = analytical_input
                report["primary_aadr_input_rebuilt"] = True
            # Only initialization metadata and the non-AADR database audit are
            # inherited. No statistical result CSV or figure is copied.
            for relative in ["results_summary.json", "analysis_manifest.json",
                             "tables/database_coverage_by_country.csv"]:
                source = REFERENCE / relative
                destination = result_dir / relative
                shutil.copyfile(source, destination)
                report["inherited_fixtures"].append({"file": relative, "sha256": sha256(source)})
            report["not_reproduced"] = [
                "raw database downloads, parsing and person-level deduplication",
                "cross-database coverage audit, including Figure 6",
                "catalogue total of 501 (the derived input contains only 489 primary rows)",
                "literature/source-scope audits and separate wild-bootstrap calibration simulation"]
            if args.mode == "aadr":
                report["not_reproduced"] = [
                    "AmtDB/aYChr cross-database coverage audit, including Figure 6",
                    "literature/source-scope audits and separate wild-bootstrap calibration simulation"]
            stage("derived_analysis", ["analysis/recompute_from_catalogue.py",
                  "--catalogue", str(analytical_input), "--analysis-output", str(result_dir),
                  *FROZEN_OPTIONS, "--aggregate-only"])
            stage("global_sensitivities", ["analysis/run_global_sensitivities.py",
                  "--catalogue", str(analytical_input), "--analysis-output", str(result_dir),
                  "--permutations", "1999", "--seed", "20260726"])
            stage("aggregate_output_checks", ["analysis/test_aggregate_release.py", str(result_dir),
                  "--input", str(analytical_input)])
        extension_results = output / "extensions"
        extension_readme = output / "extension_results.md"
        stage("statistical_extensions", [str(EXTENSION / "run_stat_extensions.py"),
              "--input", str(extension_input), "--outdir", str(extension_results),
              "--readme-output", str(extension_readme)])
        stage("extension_output_checks", [str(EXTENSION / "test_results.py"),
              "--results", str(extension_results), "--input", str(extension_input),
              "--readme", str(extension_readme)])
        stage("plot_v4_sensitivity", ["analysis/plot_v4_sensitivity.py",
              "--analysis-output", str(result_dir),
              "--extension-results", str(extension_results),
              "--output", str(result_dir / "figures/figure_5_sensitivity_v4.png")])
        comparisons = []
        for source in sorted((REFERENCE / "tables").glob("*.csv")):
            if args.mode != "raw" and source.name in INHERITED_TABLES:
                continue
            comparisons.append({"group": "main", **compare_csv(source, tables / source.name)})
        for source in sorted((EXTENSION / "results").glob("*.csv")):
            comparisons.append({"group": "extensions",
                                **compare_csv(source, extension_results / source.name)})
        report["comparisons"] = comparisons
        if any(item["status"] != "match" for item in comparisons):
            raise RuntimeError("Recomputed tables differ from the frozen reference; inspect comparisons")
        report["status"] = "passed"
        returncode = 0
    except (Exception, KeyboardInterrupt) as error:
        report["status"] = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        report["error"] = str(error)
        print(f"STOP: {error}", file=sys.stderr, flush=True)
        returncode = 1
    report["finished_at"] = datetime.now(timezone.utc).isoformat()
    report["elapsed_seconds"] = round(time.monotonic() - started, 3)
    report["generated_files"] = {str(path.relative_to(output)): sha256(path)
                                  for folder in [output / "analysis", output / "extensions"]
                                  for path in sorted(folder.rglob("*")) if path.is_file()}
    save()
    print(f"Report: {report_path}", flush=True)
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
