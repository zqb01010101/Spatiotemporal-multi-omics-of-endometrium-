#!/usr/bin/env python3
"""Sparse, chunked and resumable Nicheformer tokenization from raw counts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import time
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import psutil

from common import atomic_json, describe, setup_logging, sha256_file
from gene_utils import load_mapping, mapping_vector, normalize_nonempty_rows, remap_csr


CONTEXT = {
    "no_assay": ([5, 0, 4], [True, False, True]),
    "10x3v3_proxy": ([5, 11, 4], [True, True, True]),
}


def load_vocab(path: Path) -> dict[str, int]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or not all(isinstance(value, int) for value in raw.values()):
        raise ValueError("Expected vocab.json to contain token->integer mappings")
    return {str(key): int(value) for key, value in raw.items()}


def model_gene_tokens(model_dir: Path) -> np.ndarray:
    """Resolve gene tokens exactly from model.h5ad var order and vocab.json."""
    reference = ad.read_h5ad(model_dir / "model.h5ad", backed="r")
    try:
        genes = reference.var_names.astype(str).tolist()
    finally:
        reference.file.close()
    vocab = load_vocab(model_dir / "vocab.json")
    missing = [gene for gene in genes if gene not in vocab]
    if missing:
        raise ValueError(f"Vocabulary lacks {len(missing)} exact model genes; refusing to guess IDs")
    result = np.array([vocab[gene] for gene in genes], dtype=np.int64)
    if len(np.unique(result)) != len(result):
        raise ValueError("Model gene token IDs are not unique")
    return result


def tokenize_mapped_chunk(
    mapped: Any,
    means: np.ndarray,
    gene_tokens: np.ndarray,
    max_length: int,
    context_ids: list[int],
    context_mask: list[bool],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tokenize one CSR chunk with official ranking/tie semantics."""
    normalized, _ = normalize_nonempty_rows(mapped, 10_000.0)
    normalized.data /= means[normalized.indices]
    np.log1p(normalized.data, out=normalized.data)
    n_rows = normalized.shape[0]
    input_ids = np.zeros((n_rows, max_length), dtype=np.int32)
    attention = np.zeros((n_rows, max_length), dtype=bool)
    input_ids[:, :3] = np.asarray(context_ids, dtype=np.int32)
    attention[:, :3] = np.asarray(context_mask, dtype=bool)
    capacity = max_length - 3
    lengths = np.zeros(n_rows, dtype=np.int32)
    for row in range(n_rows):
        left, right = normalized.indptr[row : row + 2]
        indices = normalized.indices[left:right]
        values = normalized.data[left:right]
        order = np.argsort(-values)[:capacity]
        ranked = indices[order]
        count = len(ranked)
        if count:
            input_ids[row, 3 : 3 + count] = gene_tokens[ranked].astype(np.int32, copy=False)
            attention[row, 3 : 3 + count] = True
        lengths[row] = count
    return input_ids, attention, lengths


def chunk_checksum(ids: np.ndarray, mask: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(ids).view(np.uint8))
    digest.update(np.ascontiguousarray(mask).view(np.uint8))
    return digest.hexdigest()


def tokenize(
    input_h5ad: Path,
    mapping_path: Path,
    means_path: Path,
    model_dir: Path,
    model_manifest_path: Path,
    output_dir: Path,
    strategy: str,
    chunk_size: int,
    resume: bool,
    start_chunk: int | None,
    end_chunk: int | None,
    overwrite: bool,
    export_metadata: bool,
    output_root: Path,
) -> dict[str, Any]:
    """Write token arrays incrementally with atomic per-chunk completion records."""
    logger = setup_logging(output_root / f"logs/05_tokenize_{output_dir.name}.log")
    if strategy not in CONTEXT:
        raise ValueError(f"Unknown strategy {strategy!r}")
    model_manifest = json.loads(model_manifest_path.read_text(encoding="utf-8"))
    max_length = int(model_manifest["max_sequence_length"])
    n_model_genes = int(model_manifest["model_gene_count"])
    if max_length != 1500:
        logger.warning("Model max length is %d, not the historical 1500", max_length)
    means = np.load(means_path, mmap_mode="r")
    if means.shape != (n_model_genes,) or not np.all(np.isfinite(means)) or np.any(means <= 0):
        raise ValueError(f"Invalid technology means shape/values: {means.shape}")
    mapping = load_mapping(mapping_path)
    input_to_model = mapping_vector(mapping)
    gene_tokens = model_gene_tokens(model_dir)
    if gene_tokens.shape != (n_model_genes,):
        raise ValueError("Model gene token vector shape mismatch")

    identity = {
        "input_h5ad": str(input_h5ad.resolve()),
        "input_sha256": sha256_file(input_h5ad),
        "mapping_sha256": sha256_file(mapping_path),
        "technology_mean_sha256": sha256_file(means_path),
        "model_revision": model_manifest["resolved_revision"],
        "model_safetensors_sha256": model_manifest["files"]["model.safetensors"],
        "strategy": strategy,
        "context_order": ["species", "technology", "modality"],
        "context_ids": CONTEXT[strategy][0],
        "context_attention_mask": CONTEXT[strategy][1],
        "normalization_source": "layers[counts]",
        "normalization_target": 10_000.0,
        "max_length": max_length,
        "gene_capacity": max_length - 3,
        "chunk_size": chunk_size,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    marker_dir = output_dir / ".chunks"
    manifest_path = output_dir / "manifest.json"
    array_paths = {name: output_dir / name for name in ("input_ids.npy", "attention_mask.npy", "obs_names.npy")}
    if overwrite:
        for path in [*array_paths.values(), manifest_path]:
            path.unlink(missing_ok=True)
        if marker_dir.exists():
            shutil.rmtree(marker_dir)
    marker_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_h5ad, backed="r")
    try:
        if len(mapping) != adata.n_vars:
            raise ValueError("Mapping gene count differs from input H5AD")
        n_obs = adata.n_obs
        shape = (n_obs, max_length)
        n_chunks = math.ceil(n_obs / chunk_size)
        if resume:
            if not manifest_path.is_file():
                raise FileNotFoundError("--resume requested but manifest.json is absent")
            old = json.loads(manifest_path.read_text(encoding="utf-8"))
            for key, value in identity.items():
                if old.get(key) != value:
                    raise RuntimeError(f"Resume identity mismatch for {key}: {old.get(key)!r} != {value!r}")
            ids_mm = np.lib.format.open_memmap(array_paths["input_ids.npy"], mode="r+")
            mask_mm = np.lib.format.open_memmap(array_paths["attention_mask.npy"], mode="r+")
            if ids_mm.shape != shape or mask_mm.shape != shape:
                raise RuntimeError("Resume array shape mismatch")
            manifest = old
        else:
            if any(path.exists() for path in array_paths.values()) or manifest_path.exists():
                raise FileExistsError("Output exists; use --resume or explicit --overwrite")
            ids_mm = np.lib.format.open_memmap(array_paths["input_ids.npy"], mode="w+", dtype=np.int32, shape=shape)
            mask_mm = np.lib.format.open_memmap(array_paths["attention_mask.npy"], mode="w+", dtype=bool, shape=shape)
            names_text = adata.obs_names.astype(str).tolist()
            names = np.asarray(names_text, dtype=f"<U{max(map(len, names_text), default=1)}")
            np.save(array_paths["obs_names.npy"], names, allow_pickle=False)
            manifest = {**identity, "shape": list(shape), "input_ids_dtype": "int32", "attention_mask_dtype": "bool", "n_chunks": n_chunks, "completed_chunks": [], "complete": False}
            atomic_json(manifest_path, manifest)

        first = 0 if start_chunk is None else start_chunk
        last = n_chunks if end_chunk is None else min(end_chunk, n_chunks)
        if not 0 <= first <= last <= n_chunks:
            raise ValueError(f"Invalid chunk bounds [{first}, {last}) for {n_chunks} chunks")
        process = psutil.Process(os.getpid())
        for chunk_index in range(first, last):
            marker_path = marker_dir / f"chunk_{chunk_index:06d}.json"
            start = chunk_index * chunk_size
            end = min(start + chunk_size, n_obs)
            if marker_path.is_file():
                if not resume:
                    raise RuntimeError(f"Unexpected existing chunk marker {marker_path}")
                marker = json.loads(marker_path.read_text(encoding="utf-8"))
                observed = chunk_checksum(ids_mm[start:end], mask_mm[start:end])
                if marker.get("checksum") != observed:
                    raise RuntimeError(f"Checksum mismatch for completed chunk {chunk_index}")
                logger.info("Skipping verified chunk %d", chunk_index)
                continue
            began = time.perf_counter()
            rss_before = process.memory_info().rss
            raw = adata.layers["counts"][start:end].tocsr(copy=False)
            mapped = remap_csr(raw, input_to_model, n_model_genes)
            ids, mask, lengths = tokenize_mapped_chunk(mapped, means, gene_tokens, max_length, *CONTEXT[strategy])
            ids_mm[start:end] = ids
            mask_mm[start:end] = mask
            ids_mm.flush()
            mask_mm.flush()
            checksum = chunk_checksum(ids_mm[start:end], mask_mm[start:end])
            marker = {
                "chunk_index": chunk_index,
                "row_start": start,
                "row_end": end,
                "first_spot_id": str(adata.obs_names[start]),
                "last_spot_id": str(adata.obs_names[end - 1]),
                "elapsed_seconds": time.perf_counter() - began,
                "rss_before_bytes": rss_before,
                "rss_peak_observed_bytes": max(rss_before, process.memory_info().rss),
                "local_dense_created": False,
                "chunk_size": end - start,
                "gene_token_count": describe(lengths),
                "checksum": checksum,
                "complete": True,
            }
            atomic_json(marker_path, marker)
            completed = sorted({*manifest.get("completed_chunks", []), chunk_index})
            manifest["completed_chunks"] = completed
            manifest["complete"] = len(completed) == n_chunks
            atomic_json(manifest_path, manifest)
            logger.info("Completed chunk %d rows [%d, %d) in %.2fs", chunk_index, start, end, marker["elapsed_seconds"])

        if export_metadata:
            wanted = ["spot_id", "sample_id", "donor_id", "x", "y", "stage", "stage_original", "stage_numeric", "condition", "cell_type_major", "cell_type_sub", "total_counts", "n_genes", "pct_mt", "pct_ribo"]
            missing = [column for column in wanted if column not in adata.obs]
            if missing:
                raise KeyError(f"Cannot export required obs columns: {missing}")
            metadata_obs = adata.obs[wanted].copy()
            stage_rif = metadata_obs["stage"].astype(str).str.startswith("RIF_")
            sample_rif = metadata_obs["sample_id"].astype(str).str.startswith("RIF_")
            if not np.array_equal(stage_rif.to_numpy(), sample_rif.to_numpy()):
                raise ValueError("RIF markers disagree between stage and sample_id during metadata export")
            metadata_obs["analysis_condition"] = np.where(stage_rif, "RIF", "normal")
            metadata_obs.to_parquet(output_root / "metadata/obs.parquet")
            np.save(output_root / "metadata/spatial.npy", np.asarray(adata.obsm["spatial"]), allow_pickle=False)
    finally:
        adata.file.close()
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--technology-mean", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--strategy", choices=sorted(CONTEXT), required=True)
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--start-chunk", type=int)
    parser.add_argument("--end-chunk", type=int)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--export-metadata", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    tokenize(
        args.input_h5ad, args.mapping, args.technology_mean, args.model_dir,
        args.model_manifest, args.output_dir, args.strategy, args.chunk_size,
        args.resume, args.start_chunk, args.end_chunk, args.overwrite,
        args.export_metadata, args.output_root,
    )
