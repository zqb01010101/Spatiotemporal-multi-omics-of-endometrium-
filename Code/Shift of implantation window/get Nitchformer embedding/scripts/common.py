"""Shared utilities for the Nicheformer preprocessing pipeline."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import tempfile
from pathlib import Path
from typing import Any, Iterator

import numpy as np
from scipy import sparse


def setup_logging(log_path: Path, verbose: bool = False) -> logging.Logger:
    """Configure file and console logging."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(log_path.stem)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for handler in (logging.FileHandler(log_path), logging.StreamHandler()):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    """Return a streaming SHA256 digest without loading the file into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write a JSON document."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, default=json_default)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def json_default(value: Any) -> Any:
    """Convert common NumPy and path values for JSON serialization."""
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)!r}")


def csr_chunk(matrix: Any, start: int, end: int) -> sparse.csr_matrix:
    """Read a bounded row slice and ensure CSR representation."""
    chunk = matrix[start:end]
    if not sparse.issparse(chunk):
        raise TypeError(f"Expected a sparse row chunk, got {type(chunk)!r}")
    return chunk.tocsr(copy=False)


def iter_csr_chunks(matrix: Any, n_rows: int, chunk_size: int) -> Iterator[tuple[int, int, sparse.csr_matrix]]:
    """Yield bounded CSR row chunks."""
    for start in range(0, n_rows, chunk_size):
        end = min(start + chunk_size, n_rows)
        yield start, end, csr_chunk(matrix, start, end)


def describe(values: np.ndarray) -> dict[str, float | int]:
    """Return compact finite numeric distribution metrics."""
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"count": int(values.size), "finite": 0}
    return {
        "count": int(values.size),
        "finite": int(finite.size),
        "min": float(np.min(finite)),
        "q01": float(np.quantile(finite, 0.01)),
        "q25": float(np.quantile(finite, 0.25)),
        "median": float(np.median(finite)),
        "mean": float(np.mean(finite)),
        "q75": float(np.quantile(finite, 0.75)),
        "q99": float(np.quantile(finite, 0.99)),
        "max": float(np.max(finite)),
    }


def runtime_info() -> dict[str, str]:
    """Return core runtime provenance."""
    import sys

    return {
        "python": sys.version.replace("\n", " "),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "uname": " ".join(platform.uname()),
    }
