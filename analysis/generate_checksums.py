#!/usr/bin/env python3
"""Generate a deterministic, GNU-compatible SHA256SUMS.txt release manifest."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path


DEFAULT_CHECKSUM_NAME = "SHA256SUMS.txt"
CACHE_DIRECTORY_NAMES = {
    "__pycache__",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".ipynb_checkpoints",
    ".mplconfig",
}
CACHE_FILE_NAMES = {".DS_Store", "Thumbs.db"}
CACHE_SUFFIXES = {".pyc", ".pyo"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def excluded_from_checksums(path: Path, root: Path, checksum_path: Path) -> bool:
    """Return whether a path is intentionally outside the checksum manifest."""
    try:
        relative = path.relative_to(root)
    except ValueError:
        return True
    if path == checksum_path:
        return True
    if ".git" in relative.parts:
        return True
    if any(part in CACHE_DIRECTORY_NAMES for part in relative.parts):
        return True
    if path.name in CACHE_FILE_NAMES or path.suffix.lower() in CACHE_SUFFIXES:
        return True
    return False


def release_files(root: Path, checksum_path: Path) -> list[Path]:
    """List checksummed files in stable repository-relative order."""
    files: list[Path] = []
    for path in root.rglob("*"):
        if excluded_from_checksums(path, root, checksum_path):
            continue
        if path.is_symlink():
            raise ValueError(f"Release symlinks are not supported: {path}")
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if any(character in relative for character in ("\n", "\r", "\\")):
            raise ValueError(
                "Release paths may not contain newlines or backslashes: "
                f"{relative!r}"
            )
        files.append(path)
    return sorted(files, key=lambda item: item.relative_to(root).as_posix())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Write a standard two-column SHA-256 manifest for the public "
            "release tree."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root (default: parent of the analysis directory).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_CHECKSUM_NAME),
        help="Manifest path, relative to --root unless absolute.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"Repository root does not exist: {root}")
    output = args.output if args.output.is_absolute() else root / args.output
    output = output.resolve()
    try:
        output.relative_to(root)
    except ValueError as error:
        raise SystemExit("Checksum manifest must be inside the repository root") from error

    files = release_files(root, output)
    lines = [
        f"{sha256(path)}  {path.relative_to(root).as_posix()}"
        for path in files
    ]
    payload = "\n".join(lines) + ("\n" if lines else "")

    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(payload, encoding="utf-8", newline="\n")
    os.replace(temporary, output)
    print(f"Wrote {output} with {len(files)} files")


if __name__ == "__main__":
    main()
