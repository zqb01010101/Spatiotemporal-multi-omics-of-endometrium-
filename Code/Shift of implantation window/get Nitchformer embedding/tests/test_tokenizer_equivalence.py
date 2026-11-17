"""Sparse tokenizer equivalence against a safe dense reference."""

import numpy as np
from scipy import sparse


def test_sparse_ranking_matches_dense_reference(load_script):
    module = load_script("05_tokenize_nicheformer.py")
    matrix = sparse.csr_matrix([[0, 2, 5, 0], [4, 0, 1, 3]], dtype=np.float64)
    means = np.array([1, 2, 5, 1], dtype=np.float64)
    tokens = np.array([30, 31, 32, 33], dtype=np.int64)
    ids, mask, lengths = module.tokenize_mapped_chunk(matrix, means, tokens, 7, [5, 11, 4], [True] * 3)
    dense = matrix.toarray()
    dense *= (10_000 / dense.sum(axis=1))[:, None]
    dense = np.log1p(dense / means)
    expected = np.zeros((2, 7), dtype=np.int32)
    expected[:, :3] = [5, 11, 4]
    for row in range(2):
        nz = np.nonzero(dense[row])[0]
        order = nz[np.argsort(-dense[row, nz])]
        expected[row, 3 : 3 + len(order)] = tokens[order]
    np.testing.assert_array_equal(ids, expected)
    np.testing.assert_array_equal(mask, expected != 0)
    np.testing.assert_array_equal(lengths, [2, 3])
