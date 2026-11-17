#!/usr/bin/env python3
"""Extract per-spot embeddings with the pinned official pooling implementation."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModel

from common import atomic_json, setup_logging, sha256_file


REVISION = "0ba4ba162f9ad3d6c8b713761954c31be2fa69ed"


def infer(args: argparse.Namespace) -> None:
    """Stream float32 embeddings to disk, checking completed batches on resume."""
    if args.batch_size < 1 or args.threads < 1:
        raise ValueError("Batch size and thread count must be positive")
    root = args.output_dir.resolve()
    if root.exists() and any(root.iterdir()) and not args.resume:
        raise FileExistsError("Embedding output exists; choose a new directory or --resume")
    model_manifest = json.loads(args.model_manifest.read_text())
    token_manifest = json.loads((args.token_dir / "manifest.json").read_text())
    if not token_manifest.get("complete"):
        raise ValueError("Tokenization must be complete before embedding inference")
    if model_manifest["resolved_revision"] != REVISION or token_manifest["model_revision"] != REVISION:
        raise ValueError("This implementation is audited for the fixed paper revision only")
    if not model_manifest.get("source_audit_completed_before_remote_code"):
        raise ValueError("Run 02_inspect_model_assets.py before executing model code")
    for file, digest in model_manifest["files"].items():
        if file.endswith(".py") or file in {"config.json", "model.safetensors"}:
            if sha256_file(args.model_dir / file) != digest:
                raise ValueError(f"Audited model asset changed: {file}")
    ids = np.load(args.token_dir / "input_ids.npy", mmap_mode="r")
    mask = np.load(args.token_dir / "attention_mask.npy", mmap_mode="r")
    names = np.load(args.token_dir / "obs_names.npy", allow_pickle=False)
    if ids.shape != mask.shape or list(ids.shape) != token_manifest["shape"] or len(names) != len(ids):
        raise ValueError("Token shape or observation count differs from manifest")
    count = min(args.max_spots, len(ids)) if args.max_spots is not None else len(ids)
    if count < 1:
        raise ValueError("At least one spot must be selected")
    root.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(root / "inference.log")
    device = "cuda" if args.device == "auto" and torch.cuda.is_available() else args.device
    device = "cpu" if device == "auto" else device
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    torch.set_num_threads(args.threads)
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
    model = AutoModel.from_pretrained(str(args.model_dir), trust_remote_code=True, local_files_only=True)
    model = model.to(device).eval()
    dimension = int(model_manifest["hidden_dimension"])
    identity = {
        "revision": REVISION, "shape": [count, dimension], "batch_size": args.batch_size,
        "seed": args.seed, "device": device, "torch_version": torch.__version__,
        "with_context": False, "layer": -1, "precision": "float32",
        "pooling": "official get_embeddings; mean over all positions after first three, including PAD query positions",
        "token_files": {file: sha256_file(args.token_dir / file) for file in
                        ["input_ids.npy", "attention_mask.npy", "obs_names.npy", "manifest.json"]},
        "model_weights_sha256": model_manifest["files"]["model.safetensors"],
    }
    state_path = root / "manifest.json"
    array_path = root / "embeddings.npy"
    if args.resume:
        state = json.loads(state_path.read_text())
        if state["identity"] != identity:
            raise ValueError("Inference resume input, device, runtime, or parameters changed")
        result = np.lib.format.open_memmap(array_path, mode="r+")
        if list(result.shape) != identity["shape"] or result.dtype != np.float32:
            raise ValueError("Embedding memmap differs from expected shape/dtype")
        if not np.array_equal(np.load(root / "obs_names.npy", allow_pickle=False), names[:count]):
            raise ValueError("Embedding observation order changed")
    else:
        result = np.lib.format.open_memmap(array_path, mode="w+", dtype=np.float32, shape=(count, dimension))
        np.save(root / "obs_names.npy", names[:count], allow_pickle=False)
        state = {"identity": identity, "completed_batches": {}, "complete": False}
        atomic_json(state_path, state)
    for start in range(0, count, args.batch_size):
        end = min(start + args.batch_size, count)
        key = str(start)
        import hashlib
        if key in state["completed_batches"]:
            digest = hashlib.sha256(np.ascontiguousarray(result[start:end]).tobytes()).hexdigest()
            if digest != state["completed_batches"][key]:
                raise ValueError(f"Completed embedding batch checksum mismatch: {start}:{end}")
            continue
        with torch.inference_mode():
            embedding = model.get_embeddings(
                input_ids=torch.tensor(np.array(ids[start:end]), dtype=torch.long, device=device),
                attention_mask=torch.tensor(np.array(mask[start:end]), dtype=torch.bool, device=device),
                layer=-1, with_context=False,
            ).cpu().float().numpy()
        if embedding.shape != (end - start, dimension) or not np.isfinite(embedding).all():
            raise ValueError(f"Invalid embedding batch {start}:{end}")
        result[start:end] = embedding
        result.flush()
        state["completed_batches"][key] = hashlib.sha256(embedding.tobytes()).hexdigest()
        state["complete"] = len(state["completed_batches"]) == (count + args.batch_size - 1) // args.batch_size
        atomic_json(state_path, state)
        logger.info("Embedded rows [%s, %s) / %s", start, end, count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--token-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260806)
    parser.add_argument("--max-spots", type=int, help="Bounded smoke test only; omit for full inference")
    parser.add_argument("--resume", action="store_true")
    infer(parser.parse_args())
