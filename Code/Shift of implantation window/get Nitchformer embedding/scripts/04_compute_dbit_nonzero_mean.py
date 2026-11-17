#!/usr/bin/env python3
"""Compute DBiT-specific per-model-gene nonzero means from raw counts."""

from __future__ import annotations

import argparse
import importlib.metadata
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

from common import atomic_json, runtime_info, setup_logging, sha256_file
from gene_utils import load_mapping, mapping_vector, normalize_nonempty_rows, remap_csr


def read_whitelist(value: str | None) -> set[str] | None:
    """Parse a comma-separated whitelist or a one-value-per-line file."""
    if not value:
        return None
    path = Path(value)
    entries = path.read_text(encoding="utf-8").splitlines() if path.is_file() else value.split(",")
    result = {entry.strip() for entry in entries if entry.strip()}
    if not result:
        raise ValueError("Whitelist was supplied but contained no values")
    return result


def package_versions() -> dict[str, str]:
    names = ["anndata", "numpy", "scipy", "pandas", "h5py"]
    return {name: importlib.metadata.version(name) for name in names}


def compute(
    input_h5ad: Path,
    mapping_path: Path,
    model_manifest_path: Path,
    output_root: Path,
    chunk_size: int,
    sample_whitelist: set[str] | None,
    donor_whitelist: set[str] | None,
) -> dict[str, Any]:
    """Compute sparse normalized sums and nonzero denominators."""
    logger = setup_logging(output_root / "logs/04_compute_dbit_nonzero_mean.log")
    input_digest = sha256_file(input_h5ad)
    mapping_digest = sha256_file(mapping_path)
    manifest = __import__("json").loads(model_manifest_path.read_text(encoding="utf-8"))
    n_model_genes = int(manifest["model_gene_count"])
    mapping = load_mapping(mapping_path)
    vector = mapping_vector(mapping)
    sums = np.zeros(n_model_genes, dtype=np.float64)
    nonzero = np.zeros(n_model_genes, dtype=np.int64)

    adata = ad.read_h5ad(input_h5ad, backed="r")
    try:
        selected = np.ones(adata.n_obs, dtype=bool)
        if sample_whitelist is not None:
            selected &= adata.obs["sample_id"].astype(str).isin(sample_whitelist).to_numpy()
        if donor_whitelist is not None:
            selected &= adata.obs["donor_id"].astype(str).isin(donor_whitelist).to_numpy()
        row_indices = np.flatnonzero(selected)
        if not len(row_indices):
            raise ValueError("Whitelists selected zero spots")
        counts = adata.layers["counts"]
        zero_mapped_libraries = 0
        for offset in range(0, len(row_indices), chunk_size):
            rows = row_indices[offset : offset + chunk_size]
            raw = counts[rows].tocsr(copy=False)
            mapped = remap_csr(raw, vector, n_model_genes)
            normalized, totals = normalize_nonempty_rows(mapped, 10_000.0)
            zero_mapped_libraries += int(np.count_nonzero(totals == 0))
            sums += np.asarray(normalized.sum(axis=0)).ravel()
            nonzero += np.bincount(normalized.indices, minlength=n_model_genes)
            logger.info("Processed selected rows [%d, %d) of %d", offset, min(offset + chunk_size, len(row_indices)), len(row_indices))
        sample_count = int(adata.obs.iloc[row_indices]["sample_id"].nunique(dropna=True))
        donor_count = int(adata.obs.iloc[row_indices]["donor_id"].nunique(dropna=True))
    finally:
        adata.file.close()

    means = np.ones(n_model_genes, dtype=np.float32)
    expressed = nonzero > 0
    means[expressed] = (sums[expressed] / nonzero[expressed]).astype(np.float32)
    mean_path = output_root / "technology_means/dbit_global_nonzero_mean.npy"
    count_path = output_root / "technology_means/dbit_global_nonzero_count.npy"
    mean_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(mean_path, means, allow_pickle=False)
    np.save(count_path, nonzero, allow_pickle=False)
    metadata: dict[str, Any] = {
        "input_h5ad": str(input_h5ad.resolve()),
        "input_sha256": input_digest,
        "selected_spots": int(len(row_indices)),
        "sample_count": sample_count,
        "donor_count": donor_count,
        "sample_whitelist": sorted(sample_whitelist) if sample_whitelist else None,
        "donor_whitelist": sorted(donor_whitelist) if donor_whitelist else None,
        "model_revision": manifest["resolved_revision"],
        "model_gene_count": n_model_genes,
        "mapping_path": str(mapping_path.resolve()),
        "mapping_sha256": mapping_digest,
        "normalization_target": 10_000.0,
        "normalization_library": "sum of raw counts retained in mapped Nicheformer gene space",
        "zero_values_excluded_from_denominator": True,
        "never_expressed_mean_value": 1.0,
        "expressed_model_genes": int(expressed.sum()),
        "zero_mapped_library_spots": zero_mapped_libraries,
        "chunk_size": chunk_size,
        "full_matrix_densified": False,
        "computed_at_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": runtime_info(),
        "packages": package_versions(),
        "mean_sha256": sha256_file(mean_path),
        "nonzero_count_sha256": sha256_file(count_path),
    }
    atomic_json(output_root / "technology_means/dbit_global_nonzero_mean_metadata.json", metadata)
    logger.info("Wrote means for %d genes; %d were expressed", n_model_genes, int(expressed.sum()))
    return metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--sample-whitelist")
    parser.add_argument("--donor-whitelist")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    compute(
        args.input_h5ad,
        args.mapping,
        args.model_manifest,
        args.output_root,
        args.chunk_size,
        read_whitelist(args.sample_whitelist),
        read_whitelist(args.donor_whitelist),
    )
