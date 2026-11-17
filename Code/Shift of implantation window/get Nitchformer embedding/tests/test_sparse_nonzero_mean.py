"""Sparse nonzero mean correctness tests."""

import numpy as np
from scipy import sparse

from gene_utils import normalize_nonempty_rows


def test_sparse_nonzero_mean_matches_dense_reference():
    counts = sparse.csr_matrix([[1, 0, 3], [0, 2, 2], [0, 0, 0]], dtype=np.float64)
    normalized, _ = normalize_nonempty_rows(counts, 10_000)
    sums = np.asarray(normalized.sum(axis=0)).ravel()
    nonzero = np.bincount(normalized.indices, minlength=3)
    observed = np.ones(3)
    observed[nonzero > 0] = sums[nonzero > 0] / nonzero[nonzero > 0]
    dense = counts.toarray()
    totals = dense.sum(axis=1)
    dense = dense * np.divide(10_000, totals, out=np.zeros_like(totals), where=totals > 0)[:, None]
    expected = np.array([dense[dense[:, j] > 0, j].mean() if np.any(dense[:, j] > 0) else 1 for j in range(3)])
    np.testing.assert_allclose(observed, expected)
