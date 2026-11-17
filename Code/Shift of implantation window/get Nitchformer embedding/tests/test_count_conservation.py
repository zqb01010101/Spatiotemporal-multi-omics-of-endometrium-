"""Duplicate gene aggregation and count conservation tests."""

import numpy as np
from scipy import sparse

from gene_utils import remap_csr


def test_duplicate_mapping_sums_and_conserves_retained_counts():
    counts = sparse.csr_matrix([[1, 2, 7], [3, 4, 0]], dtype=np.int32)
    vector = np.array([0, 0, -1])
    mapped = remap_csr(counts, vector, 2)
    np.testing.assert_array_equal(mapped.toarray(), [[3, 0], [7, 0]])
    assert mapped.sum() == counts[:, :2].sum()
