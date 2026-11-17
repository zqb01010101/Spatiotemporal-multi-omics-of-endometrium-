"""Sparse gene-space mapping helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse


def load_mapping(path: Path, n_input_genes: int | None = None) -> pd.DataFrame:
    """Load and validate the input-to-model mapping table."""
    frame = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = {"input_index", "model_index", "gene_token_id", "mapped"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Mapping table lacks columns: {sorted(missing)}")
    if n_input_genes is not None and len(frame) != n_input_genes:
        raise ValueError(f"Mapping has {len(frame)} rows but input has {n_input_genes} genes")
    if not np.array_equal(frame["input_index"].to_numpy(), np.arange(len(frame))):
        raise ValueError("Mapping input_index must be complete and ordered")
    return frame


def mapping_vector(frame: pd.DataFrame) -> np.ndarray:
    """Return input-column to model-column indices with -1 for unmapped."""
    vector = np.full(len(frame), -1, dtype=np.int64)
    mapped = frame["mapped"].astype(str).str.lower().isin(["true", "1"])
    vector[mapped.to_numpy()] = frame.loc[mapped, "model_index"].astype(int).to_numpy()
    return vector


def gene_token_vector(frame: pd.DataFrame, n_model_genes: int) -> np.ndarray:
    """Return a unique token ID for every model gene represented by the mapping."""
    result = np.full(n_model_genes, -1, dtype=np.int64)
    mapped = frame["mapped"].astype(str).str.lower().isin(["true", "1"])
    pairs = frame.loc[mapped, ["model_index", "gene_token_id"]].drop_duplicates()
    for model_index, group in pairs.groupby("model_index"):
        values = group["gene_token_id"].astype(int).unique()
        if len(values) != 1:
            raise ValueError(f"Model gene {model_index} has conflicting token IDs: {values}")
        result[int(model_index)] = values[0]
    return result


def remap_csr(chunk: sparse.csr_matrix, input_to_model: np.ndarray, n_model_genes: int) -> sparse.csr_matrix:
    """Map columns and sum duplicate mappings using sparse COO canonicalization."""
    chunk = chunk.tocsr(copy=False)
    rows = np.repeat(np.arange(chunk.shape[0], dtype=np.int64), np.diff(chunk.indptr))
    cols = input_to_model[chunk.indices]
    keep = cols >= 0
    mapped = sparse.coo_matrix(
        (chunk.data[keep], (rows[keep], cols[keep])),
        shape=(chunk.shape[0], n_model_genes),
        dtype=np.result_type(chunk.dtype, np.float32),
    ).tocsr()
    mapped.sum_duplicates()
    mapped.eliminate_zeros()
    return mapped


def normalize_nonempty_rows(matrix: sparse.csr_matrix, target_sum: float = 10_000.0) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Library-size normalize CSR rows without densifying the matrix."""
    result = matrix.astype(np.float64, copy=True)
    totals = np.asarray(result.sum(axis=1)).ravel()
    scales = np.divide(target_sum, totals, out=np.zeros_like(totals), where=totals > 0)
    result.data *= np.repeat(scales, np.diff(result.indptr))
    return result, totals
