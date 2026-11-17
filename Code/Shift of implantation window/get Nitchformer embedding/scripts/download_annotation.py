#!/usr/bin/env python3
"""Download and checksum the exact Ensembl annotation used for the RIF run."""

from __future__ import annotations

import argparse
import os
import tempfile
import urllib.request
from pathlib import Path

from common import atomic_json, setup_logging, sha256_file


URL = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
SHA256 = "8c87436bab973c871447743c2f5cbb2b557be0523fe6efe3ff33307a4c836002"
FILENAME = "Homo_sapiens.GRCh38.112.gtf.gz"


def download(output_root: Path, local_gtf: Path | None = None) -> Path:
    """Verify an existing annotation or stream an atomic, verified download."""
    logger = setup_logging(output_root / "logs/download_annotation.log")
    destination = output_root / "model_assets/annotations" / FILENAME
    destination.parent.mkdir(parents=True, exist_ok=True)
    if local_gtf is not None:
        destination = local_gtf.resolve(strict=True)
    if not destination.is_file():
        fd, temporary = tempfile.mkstemp(dir=destination.parent, prefix=".annotation.")
        try:
            with os.fdopen(fd, "wb") as target:
                with urllib.request.urlopen(URL, timeout=120) as source:
                    while block := source.read(8 * 1024 * 1024):
                        target.write(block)
            if sha256_file(Path(temporary)) != SHA256:
                raise ValueError("Downloaded annotation SHA256 differs from the reference")
            os.replace(temporary, destination)
        finally:
            Path(temporary).unlink(missing_ok=True)
    actual = sha256_file(destination)
    if actual != SHA256:
        raise ValueError(f"Annotation SHA256 mismatch: {actual}; expected {SHA256}")
    atomic_json(output_root / "model_assets/annotations/manifest.json", {
        "source_url": URL, "sha256": actual, "release": "Ensembl 112 / GRCh38",
        "path": str(destination.resolve()), "supplied_locally": local_gtf is not None,
    })
    logger.info("Verified annotation %s", destination)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--local-gtf", type=Path)
    args = parser.parse_args()
    download(args.output_root, args.local_gtf)
