#!/usr/bin/env python3
"""Acquire frozen metadata with strict hashes; never replace a conflicting file.

The registry identifies bytes, not merely release names. AmtDB v1.009 has no
verified download locator: supply its original CSV with --local amtdb=PATH.
Files are staged beside their destination, verified, then installed atomically
without clobbering an existing path. Failed/interrupted staging files are removed.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import http.client
import json
import os
from pathlib import Path
import socket
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from verify_inputs import EXPECTED, sha256

ROOT = Path(__file__).resolve().parents[1]
RESOURCE_KEYS = {"AADR": "aadr", "AmtDB": "amtdb", "aYChr-DB": "aychr"}
TRANSIENT_HTTP = {408, 429, 500, 502, 503, 504}
MAX_BYTES = 128 * 1024 * 1024
BLOCK = 1024 * 1024


def load_sources(registry: Path) -> dict[str, dict[str, str]]:
    with registry.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    sources = {}
    for row in rows:
        name = RESOURCE_KEYS.get(row.get("resource", ""))
        if name is None or name in sources:
            raise ValueError("Unknown or duplicate resource in source registry")
        if row.get("sha256") != EXPECTED[name]:
            raise ValueError(f"{name}: registry hash differs from frozen verify_inputs.py")
        path = Path(row.get("expected_path", ""))
        if path.is_absolute() or ".." in path.parts or path.parts[:2] != ("data", "raw"):
            raise ValueError(f"{name}: expected_path must be a relative path below data/raw")
        url = row.get("download_url", "").strip()
        parsed = urlsplit(url)
        if url and (parsed.scheme != "https" or not parsed.netloc or parsed.username or parsed.password):
            raise ValueError(f"{name}: download_url must be a credential-free HTTPS URL")
        sources[name] = {**row, "key": name, "download_url": url}
    if set(sources) != set(EXPECTED):
        raise ValueError("Registry must contain all three frozen sources")
    return sources


def is_transient(error: BaseException) -> bool:
    if isinstance(error, HTTPError):
        return error.code in TRANSIENT_HTTP
    if isinstance(error, URLError):
        return isinstance(error.reason, BaseException) and is_transient(error.reason)
    if isinstance(error, socket.gaierror):
        return error.errno == socket.EAI_AGAIN
    return isinstance(error, (TimeoutError, ConnectionError, http.client.IncompleteRead))


def _install(stream, destination: Path, expected: str, result: dict, timeout: float) -> None:
    """Copy, hash and install with an atomic, no-overwrite hard link."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    deadline = time.monotonic() + timeout
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent, prefix=destination.name + ".", suffix=".part", delete=False
        ) as output:
            temporary = Path(output.name)
            digest, length = hashlib.sha256(), 0
            while True:
                chunk = stream.read(BLOCK)
                if time.monotonic() > deadline:
                    raise TimeoutError("Transfer exceeded the per-attempt time limit")
                if not chunk:
                    break
                length += len(chunk)
                if length > MAX_BYTES:
                    raise ValueError("Metadata file exceeds the 128 MiB download limit")
                output.write(chunk)
                digest.update(chunk)
            output.flush()
            os.fsync(output.fileno())
        result.update(observed_sha256=digest.hexdigest(), bytes=length)
        if result["observed_sha256"] != expected:
            result.update(status="HASH_MISMATCH", message="Acquired bytes do not match the frozen SHA-256; destination was not written.")
            return
        try:
            os.link(temporary, destination)
        except FileExistsError:
            # A concurrent invocation may have completed during this transfer.
            if destination.is_file() and not destination.is_symlink() and sha256(destination) == expected:
                result.update(status="READY", message="Another process installed the same verified bytes.")
            else:
                result.update(status="CONFLICT", message="Destination appeared during acquisition; left unchanged.")
            return
        result.update(status="READY", message="Frozen SHA-256 verified; file installed atomically.")
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def acquire(source: dict, output_root: Path, *, local: Path | None = None,
            timeout: float = 30, retries: int = 2, offline: bool = False) -> dict:
    destination = output_root / source["expected_path"]
    result = {
        "resource": source["key"], "release": source["release"],
        "path": str(destination), "expected_sha256": source["sha256"],
        "observed_sha256": None, "bytes": None, "attempts": 0,
        "source": str(local) if local else source["download_url"] or None,
        "mode": "local" if local else "download", "status": "NOT_ATTEMPTED",
    }
    try:
        if not destination.resolve().is_relative_to(output_root.resolve()):
            raise ValueError("Destination resolves outside output root")
        if destination.exists() or destination.is_symlink():
            if destination.is_symlink() or not destination.is_file():
                result.update(status="CONFLICT", message="Existing destination is not a regular file; left unchanged.")
            else:
                result.update(observed_sha256=sha256(destination), bytes=destination.stat().st_size, mode="existing")
                matching = result["observed_sha256"] == source["sha256"]
                result.update(status="READY" if matching else "CONFLICT",
                              message="Existing file matches frozen SHA-256." if matching else
                              "Existing file has a different SHA-256; move it aside explicitly before retrying.")
            return result
        if local is not None:
            result["attempts"] = 1
            with local.open("rb") as stream:
                _install(stream, destination, source["sha256"], result, timeout)
            return result
        if not source["download_url"]:
            result.update(status="BLOCKED", message=(
                "No verified immutable AmtDB v1.009 metadata CSV URL is available. Obtain the original "
                "v1.009 (2024-02-28) export from a retained snapshot or the database authors and use "
                "--local amtdb=/path/to/export.csv. Its SHA-256 must match expected_sha256. "
                "The live v1.010 export cannot substitute for this frozen input."
            ))
            return result
        if offline:
            result.update(status="MISSING", message="File is absent; offline mode prevents downloading it.")
            return result
        for attempt in range(retries + 1):
            result["attempts"] = attempt + 1
            try:
                request = Request(source["download_url"], headers={"User-Agent": "central-asia-uniparental-adna/frozen-sources"})
                with urlopen(request, timeout=timeout) as stream:
                    if urlsplit(stream.geturl()).scheme != "https":
                        raise ValueError("Refusing a redirect to non-HTTPS transport")
                    _install(stream, destination, source["sha256"], result, timeout)
                return result
            except (OSError, http.client.HTTPException) as error:
                if not is_transient(error) or attempt == retries:
                    raise
                time.sleep(min(2 ** attempt, 4))
    except KeyboardInterrupt:
        result.update(status="INTERRUPTED", message="Acquisition interrupted; incomplete staging file removed.")
    except (OSError, ValueError, http.client.HTTPException) as error:
        result.update(status="ERROR", message=f"{type(error).__name__}: {error}")
    return result


def _bounded_timeout(value: str) -> float:
    number = float(value)
    if not 1 <= number <= 120:
        raise argparse.ArgumentTypeError("timeout must be between 1 and 120 seconds")
    return number


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=ROOT / "data/SOURCES.tsv")
    parser.add_argument("--output-root", type=Path, default=ROOT, help="Root under which data/raw/... is created")
    parser.add_argument("--resource", action="append", choices=sorted(EXPECTED), help="Repeat to select resources; default: all")
    parser.add_argument("--local", action="append", default=[], metavar="RESOURCE=PATH", help="Import and verify retained bytes, e.g. amtdb=/path/export.csv")
    parser.add_argument("--offline", action="store_true", help="Only check existing files or import local files")
    parser.add_argument("--timeout", type=_bounded_timeout, default=30, help="Socket timeout and streaming time limit per attempt (1-120 seconds)")
    parser.add_argument("--retries", type=int, choices=range(6), default=2, help="Additional attempts for transient network failures (0-5)")
    parser.add_argument("--report", type=Path, help="Also save the JSON report here; stdout always contains the report")
    args = parser.parse_args(argv)
    local = {}
    for item in args.local:
        name, separator, path = item.partition("=")
        if not separator or name not in EXPECTED or not path or name in local:
            parser.error("--local requires a unique known RESOURCE=PATH")
        local[name] = Path(path).expanduser()
    selected = list(dict.fromkeys(args.resource or EXPECTED))
    if set(local) - set(selected):
        parser.error("A --local resource was not selected by --resource")
    try:
        sources = load_sources(args.registry)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if args.report:
        protected = [args.registry, *local.values(), *(args.output_root / s["expected_path"] for s in sources.values())]
        if args.report.resolve() in {p.resolve() for p in protected}:
            parser.error("Report path must differ from the registry, source and destination files")
    results, interrupted = [], False
    for name in selected:
        if interrupted:
            results.append({"resource": name, "status": "NOT_ATTEMPTED", "message": "Stopped after interruption."})
            continue
        result = acquire(sources[name], args.output_root, local=local.get(name),
                         timeout=args.timeout, retries=args.retries, offline=args.offline)
        results.append(result)
        interrupted = result["status"] == "INTERRUPTED"
    ready = all(result["status"] == "READY" for result in results)
    report = {"schema_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
              "status": "READY" if ready else "INCOMPLETE", "selected_resources": selected,
              "all_three_sources_ready": ready and set(selected) == set(EXPECTED), "sources": results}
    serialized = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    print(serialized, end="")
    if args.report:
        try:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(serialized, encoding="utf-8")
        except OSError as error:
            parser.exit(2, f"Could not save JSON report: {error}\n")
    return 130 if interrupted else 0 if ready else 1


if __name__ == "__main__":
    raise SystemExit(main())
