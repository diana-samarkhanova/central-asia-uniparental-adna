#!/usr/bin/env python3
"""Fail fast unless the three inputs match the frozen analysis release."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


EXPECTED = {
    "aadr": "98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c",
    "amtdb": "531e8ee8fae181124f5a9b77b6fe8d677e64e35b815be2a3965020244fe31057",
    "aychr": "e297110a18cba73d4044e8a95c0fae98d7f48633ad6ccfae1cd364e460eb1b3c",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aadr", required=True, type=Path)
    parser.add_argument("--amtdb", required=True, type=Path)
    parser.add_argument("--aychr", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    failures = []
    for name, expected in EXPECTED.items():
        path = getattr(args, name)
        if not path.is_file():
            failures.append(f"{name}: missing file: {path}")
            continue
        observed = sha256(path)
        if observed != expected:
            failures.append(
                f"{name}: SHA-256 mismatch; expected {expected}, observed {observed}"
            )
        else:
            print(f"{name}: OK")
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
