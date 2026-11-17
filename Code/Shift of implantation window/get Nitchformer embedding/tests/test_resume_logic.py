"""Chunk checksum behavior used by resume verification."""

import json

import numpy as np


def test_completed_chunk_is_verified_before_skip(tmp_path, load_script):
    module = load_script("05_tokenize_nicheformer.py")
    ids = np.array([[5, 0, 4, 30, 0]], dtype=np.int32)
    mask = ids != 0
    checksum = module.chunk_checksum(ids, mask)
    marker = tmp_path / "chunk.json"
    marker.write_text(json.dumps({"checksum": checksum, "complete": True}))
    before = ids.copy()
    assert json.loads(marker.read_text())["checksum"] == module.chunk_checksum(ids, mask)
    np.testing.assert_array_equal(ids, before)
    ids[0, 3] = 31
    assert json.loads(marker.read_text())["checksum"] != module.chunk_checksum(ids, mask)
