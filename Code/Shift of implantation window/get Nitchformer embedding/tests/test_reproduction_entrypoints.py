"""Test release checkpoints, fail-fast gating, and embedding inference resume."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import torch


def test_checkpoint_rejects_modified_product(tmp_path, load_script):
    module = load_script("run_pipeline.py")
    path = tmp_path / "product.txt"
    path.write_text("original")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    module.verify_products(tmp_path, {path.name: digest})
    path.write_text("changed")
    with pytest.raises(RuntimeError, match="Checkpoint product"):
        module.verify_products(tmp_path, {path.name: digest})


def test_runner_stops_before_full_if_validation_fails(tmp_path, load_script, monkeypatch):
    module = load_script("run_pipeline.py")
    inp = tmp_path / "input.h5ad"
    inp.write_text("synthetic input; subprocesses are mocked")
    output = tmp_path / "output"
    args = argparse.Namespace(input_h5ad=inp, output_root=output, annotation_gtf=None,
                              seed=20260806, chunk_size=256, resume=False,
                              dry_run=False, stop_after=None)
    calls = []
    monkeypatch.setattr(module.sys, "version_info", (3, 10, 20))
    monkeypatch.setattr(module, "runtime_info", lambda: {"python": "synthetic test"})
    monkeypatch.setattr(module, "stage_commands", lambda _: [
        ("validate", ["mock", "validate"], ["reports/dev_validation_metrics.json"]),
        ("full_no_assay", ["mock", "full"], []),
    ])

    def fake_run(command, **kwargs):
        calls.append(command)
        (output / "reports/dev_validation_metrics.json").write_text(json.dumps({"all_passed": False}))

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="validation failed"):
        module.run(args)
    assert calls == [["mock", "validate"]]


def test_embedding_resume_verifies_arrays_and_preserves_order(tmp_path, load_script, monkeypatch):
    module = load_script("09_infer_embeddings.py")
    tokens = tmp_path / "tokens"
    model_dir = tmp_path / "model"
    tokens.mkdir()
    model_dir.mkdir()
    weights = model_dir / "model.safetensors"
    weights.write_bytes(b"synthetic weights; model loader mocked")
    digest = hashlib.sha256(weights.read_bytes()).hexdigest()
    manifest = tmp_path / "model_manifest.json"
    manifest.write_text(json.dumps({"resolved_revision": module.REVISION,
        "source_audit_completed_before_remote_code": True,
        "hidden_dimension": 2, "files": {"model.safetensors": digest}}))
    ids = np.array([[5, 0, 4, 30, 0], [5, 0, 4, 31, 32], [5, 0, 4, 33, 0]], dtype=np.int32)
    np.save(tokens / "input_ids.npy", ids)
    np.save(tokens / "attention_mask.npy", ids != 0)
    np.save(tokens / "obs_names.npy", np.array(["z", "a", "b"]))
    (tokens / "manifest.json").write_text(json.dumps({"complete": True, "shape": list(ids.shape), "model_revision": module.REVISION}))

    class Model:
        calls = 0

        def to(self, device):
            return self

        def eval(self):
            return self

        def get_embeddings(self, input_ids, **kwargs):
            self.calls += 1
            return input_ids[:, :2].float()

    model = Model()
    monkeypatch.setattr(module.AutoModel, "from_pretrained", lambda *a, **kw: model)
    args = argparse.Namespace(output_dir=tmp_path / "embeddings", token_dir=tokens,
        model_dir=model_dir, model_manifest=manifest, device="cpu", batch_size=2,
        threads=1, seed=20260806, max_spots=None, resume=False)
    module.infer(args)
    assert model.calls == 2
    np.testing.assert_array_equal(np.load(args.output_dir / "embeddings.npy"), ids[:, :2])
    args.resume = True
    module.infer(args)
    assert model.calls == 2
    result = np.load(args.output_dir / "embeddings.npy", mmap_mode="r+")
    result[0, 0] = -123
    result.flush()
    with pytest.raises(ValueError, match="checksum mismatch"):
        module.infer(args)
