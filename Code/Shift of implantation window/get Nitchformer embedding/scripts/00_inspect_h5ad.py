#!/usr/bin/env python3
"""Inspect the source AnnData using sparse, read-only row chunks."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from common import atomic_json, describe, iter_csr_chunks, runtime_info, setup_logging, sha256_file


def categorical_counts(frame: pd.DataFrame, column: str) -> dict[str, int]:
    """Return value counts including missing values."""
    if column not in frame:
        return {}
    return {str(key): int(value) for key, value in frame[column].value_counts(dropna=False).items()}


def inspect(input_h5ad: Path, output_root: Path, chunk_size: int) -> dict[str, Any]:
    """Profile source metadata and raw counts without densifying full matrices."""
    logger = setup_logging(output_root / "logs/00_inspect_h5ad.log")
    digest = sha256_file(input_h5ad)
    logger.info("Input SHA256: %s", digest)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    try:
        if "counts" not in adata.layers:
            raise KeyError("Required layer adata.layers['counts'] is absent")
        counts = adata.layers["counts"]
        if not sparse.issparse(counts) and "SparseDataset" not in type(counts).__name__:
            raise TypeError(f"counts must be sparse, got {type(counts)!r}")
        n_obs, n_vars = adata.shape
        per_spot_total = np.zeros(n_obs, dtype=np.float64)
        per_spot_genes = np.zeros(n_obs, dtype=np.int32)
        per_gene_nnz = np.zeros(n_vars, dtype=np.int64)
        negative = 0
        noninteger = 0
        counts_nnz = 0
        counts_min = np.inf
        counts_max = -np.inf
        for start, end, chunk in iter_csr_chunks(counts, n_obs, chunk_size):
            data = chunk.data
            counts_nnz += int(chunk.nnz)
            negative += int(np.count_nonzero(data < 0))
            noninteger += int(np.count_nonzero(~np.isclose(data, np.rint(data), atol=1e-6)))
            if data.size:
                counts_min = min(counts_min, float(data.min()))
                counts_max = max(counts_max, float(data.max()))
            per_spot_total[start:end] = np.asarray(chunk.sum(axis=1)).ravel()
            per_spot_genes[start:end] = np.diff(chunk.indptr)
            per_gene_nnz += np.bincount(chunk.indices, minlength=n_vars)
            logger.info("Inspected rows [%d, %d)", start, end)

        obs = adata.obs.copy()
        var = adata.var.copy()
        obs.index = obs.index.astype(str)
        var.index = var.index.astype(str)
        spatial_exists = "spatial" in adata.obsm
        spatial = np.asarray(adata.obsm["spatial"]) if spatial_exists else None
        spot_id_equal = bool("spot_id" in obs and np.array_equal(obs["spot_id"].astype(str).to_numpy(), obs.index.to_numpy()))
        xy_equal = False
        if spatial is not None and {"x", "y"}.issubset(obs.columns) and spatial.shape == (n_obs, 2):
            xy_equal = bool(np.allclose(obs[["x", "y"]].to_numpy(dtype=float), spatial, equal_nan=True))

        x_obj = adata.X
        metrics: dict[str, Any] = {
            "input_path": str(input_h5ad.resolve()),
            "input_sha256": digest,
            "shape": [int(n_obs), int(n_vars)],
            "x_present": x_obj is not None,
            "x_type": type(x_obj).__name__,
            "x_sparse": bool(sparse.issparse(x_obj) or type(x_obj).__name__ in {"CSRDataset", "CSCDataset"}),
            "x_dtype": str(getattr(x_obj, "dtype", "unknown")),
            "x_nnz": int(getattr(x_obj, "nnz", x_obj.group["data"].shape[0] if hasattr(x_obj, "group") else -1)),
            "counts_present": True,
            "counts_type": type(counts).__name__,
            "counts_sparse": True,
            "counts_dtype": str(counts.dtype),
            "counts_nnz": counts_nnz,
            "counts_data_min": None if not np.isfinite(counts_min) else counts_min,
            "counts_data_max": None if not np.isfinite(counts_max) else counts_max,
            "counts_negative_values": negative,
            "counts_noninteger_values": noninteger,
            "obs_names_unique": bool(adata.obs_names.is_unique),
            "spot_id_present": "spot_id" in obs,
            "spot_id_unique": bool(obs["spot_id"].is_unique) if "spot_id" in obs else False,
            "spot_id_matches_obs_names": spot_id_equal,
            "var_names_unique": bool(adata.var_names.is_unique),
            "spatial_present": spatial_exists,
            "spatial_shape": list(spatial.shape) if spatial is not None else None,
            "xy_matches_spatial": xy_equal,
            "sample_count": int(obs["sample_id"].nunique(dropna=True)) if "sample_id" in obs else None,
            "donor_count": int(obs["donor_id"].nunique(dropna=True)) if "donor_id" in obs else None,
            "stage_distribution": categorical_counts(obs, "stage"),
            "cell_type_major_distribution": categorical_counts(obs, "cell_type_major"),
            "condition_distribution": categorical_counts(obs, "condition"),
            "batch_distribution": categorical_counts(obs, "batch"),
            "zero_library_spots": int(np.count_nonzero(per_spot_total == 0)),
            "empty_spots": int(np.count_nonzero(per_spot_genes == 0)),
            "per_spot_total_counts": describe(per_spot_total),
            "per_spot_detected_genes": describe(per_spot_genes),
            "per_gene_nonzero_spots": describe(per_gene_nnz),
            "chunk_size": chunk_size,
            "full_matrix_densified": False,
            "runtime": runtime_info(),
        }
        (output_root / "metadata").mkdir(parents=True, exist_ok=True)
        obs.to_parquet(output_root / "metadata/input_obs.parquet")
        var.to_parquet(output_root / "metadata/input_var.parquet")
        atomic_json(output_root / "reports/input_h5ad_metrics.json", metrics)
        report = f"""# Input H5AD Quality Report

- Input: `{input_h5ad.resolve()}`
- SHA256: `{digest}`
- Shape: **{n_obs:,} spots × {n_vars:,} genes**
- X: `{metrics['x_type']}`, `{metrics['x_dtype']}`, sparse={metrics['x_sparse']}, nnz={metrics['x_nnz']:,}
- counts: `{metrics['counts_type']}`, `{metrics['counts_dtype']}`, nnz={counts_nnz:,}
- Negative / non-integer nonzero counts: **{negative:,} / {noninteger:,}**
- Empty or zero-library spots: **{metrics['empty_spots']:,} / {metrics['zero_library_spots']:,}**
- Unique obs names / spot IDs / exact match: **{metrics['obs_names_unique']} / {metrics['spot_id_unique']} / {spot_id_equal}**
- Unique var names: **{metrics['var_names_unique']}**
- Spatial shape / x-y agreement: **{metrics['spatial_shape']} / {xy_equal}**
- Samples / donors: **{metrics['sample_count']} / {metrics['donor_count']}**

## Stage Distribution

{pd.Series(metrics['stage_distribution']).to_markdown()}

## Major Cell-Type Distribution

{pd.Series(metrics['cell_type_major_distribution']).to_markdown()}

## Numeric Summaries

- Library size: `{metrics['per_spot_total_counts']}`
- Detected genes per spot: `{metrics['per_spot_detected_genes']}`
- Nonzero spots per gene: `{metrics['per_gene_nonzero_spots']}`

The full expression matrices were never converted to dense form; inspection used {chunk_size}-spot CSR row chunks.
"""
        (output_root / "reports/input_h5ad_report.md").write_text(report, encoding="utf-8")
        return metrics
    finally:
        adata.file.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=256)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    inspect(args.input_h5ad, args.output_root, args.chunk_size)
