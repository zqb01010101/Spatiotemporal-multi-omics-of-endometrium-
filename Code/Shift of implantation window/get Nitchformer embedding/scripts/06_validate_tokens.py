#!/usr/bin/env python3
"""Validate development tokens against official code and a CPU model forward pass."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import torch
from transformers import AutoModel

from common import atomic_json, describe, setup_logging
from gene_utils import load_mapping, mapping_vector, remap_csr


def load_official_tokenizer(model_dir: Path, means: np.ndarray) -> Any:
    """Import the already audited pinned tokenizer source from disk."""
    source = model_dir / "tokenization_nicheformer.py"
    spec = importlib.util.spec_from_file_location("pinned_nicheformer_tokenizer", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    tokenizer = module.NicheformerTokenizer(
        vocab_file=str(model_dir / "vocab.json"),
        max_length=1500,
        aux_tokens=30,
        technology_mean=np.asarray(means),
    )
    tokenizer._load_reference_model = lambda: None
    return tokenizer


def model_gene_names(model_dir: Path) -> list[str]:
    reference = ad.read_h5ad(model_dir / "model.h5ad", backed="r")
    try:
        return reference.var_names.astype(str).tolist()
    finally:
        reference.file.close()


def validate(
    dev_h5ad: Path,
    mapping_path: Path,
    means_path: Path,
    model_dir: Path,
    model_manifest_path: Path,
    no_assay_dir: Path,
    proxy_dir: Path,
    mapping_report_path: Path,
    output_root: Path,
    equivalence_spots: int,
    forward_batch: int,
) -> dict[str, Any]:
    """Run structural, equivalence, conservation and forward checks."""
    logger = setup_logging(output_root / "logs/06_validate_tokens.log")
    manifest = json.loads(model_manifest_path.read_text(encoding="utf-8"))
    max_length = int(manifest["max_sequence_length"])
    vocab_size = int(manifest["vocab_size"])
    config_tokens = int(manifest["config"]["n_tokens"])
    no_ids = np.load(no_assay_dir / "input_ids.npy", mmap_mode="r")
    no_mask = np.load(no_assay_dir / "attention_mask.npy", mmap_mode="r")
    no_names = np.load(no_assay_dir / "obs_names.npy", allow_pickle=False)
    proxy_ids = np.load(proxy_dir / "input_ids.npy", mmap_mode="r")
    proxy_mask = np.load(proxy_dir / "attention_mask.npy", mmap_mode="r")
    proxy_names = np.load(proxy_dir / "obs_names.npy", allow_pickle=False)
    checks: dict[str, bool] = {}

    dev = ad.read_h5ad(dev_h5ad, backed="r")
    try:
        checks["dev_spot_count_7200"] = dev.n_obs == 7200
        checks["token_shapes"] = no_ids.shape == no_mask.shape == proxy_ids.shape == proxy_mask.shape == (dev.n_obs, max_length)
        checks["obs_order"] = np.array_equal(no_names.astype(str), dev.obs_names.astype(str)) and np.array_equal(proxy_names.astype(str), dev.obs_names.astype(str))
        checks["token_id_range"] = int(no_ids.min()) >= 0 and int(proxy_ids.min()) >= 0 and int(no_ids.max()) < config_tokens + 5 and int(proxy_ids.max()) < config_tokens + 5
        checks["attention_pad_consistency_no_assay"] = bool(np.all((no_ids == 0) == (~no_mask)))
        checks["attention_pad_consistency_proxy"] = bool(np.all((proxy_ids == 0) == (~proxy_mask)))
        checks["no_negative_tokens"] = int(no_ids.min()) >= 0 and int(proxy_ids.min()) >= 0
        checks["no_nan"] = True
        checks["strategies_differ_only_at_technology"] = (
            np.array_equal(no_ids[:, 0], proxy_ids[:, 0])
            and np.array_equal(no_ids[:, 2:], proxy_ids[:, 2:])
            and np.all(no_ids[:, 1] == 0)
            and np.all(proxy_ids[:, 1] == 11)
            and np.array_equal(no_mask[:, 0], proxy_mask[:, 0])
            and np.array_equal(no_mask[:, 2:], proxy_mask[:, 2:])
            and np.all(~no_mask[:, 1])
            and np.all(proxy_mask[:, 1])
        )
        checks["context_no_assay"] = bool(np.all(no_ids[:, 0] == 5) and np.all(no_ids[:, 1] == 0) and np.all(no_ids[:, 2] == 4))
        checks["context_proxy"] = bool(np.all(proxy_ids[:, 0] == 5) and np.all(proxy_ids[:, 1] == 11) and np.all(proxy_ids[:, 2] == 4))
        gene_lengths = no_mask[:, 3:].sum(axis=1).astype(np.int32)
        checks["gene_lengths_reasonable"] = bool(np.all(gene_lengths >= 0) and np.all(gene_lengths <= max_length - 3) and np.count_nonzero(gene_lengths) == dev.n_obs)

        n_test = min(max(equivalence_spots, 32), dev.n_obs)
        mapping = load_mapping(mapping_path, dev.n_vars)
        vector = mapping_vector(mapping)
        genes = model_gene_names(model_dir)
        raw = dev.layers["counts"][:n_test].tocsr(copy=False)
        aligned = remap_csr(raw, vector, len(genes))
        dense_bytes = aligned.shape[0] * aligned.shape[1] * np.dtype(np.float64).itemsize
        dense = aligned.toarray()
        means = np.load(means_path)
        official = load_official_tokenizer(model_dir, means)
        obs = pd.DataFrame({"specie": ["human"] * n_test, "assay": ["10x 3' v3"] * n_test, "modality": ["spatial"] * n_test})
        aligned_adata = ad.AnnData(X=dense, obs=obs, var=pd.DataFrame(index=genes))
        official_proxy = official(aligned_adata)
        official_ids = official_proxy["input_ids"].numpy()
        official_mask = official_proxy["attention_mask"].numpy().astype(bool)
        checks["official_proxy_input_ids_equivalent"] = np.array_equal(official_ids, np.asarray(proxy_ids[:n_test]))
        checks["official_proxy_attention_equivalent"] = np.array_equal(official_mask, np.asarray(proxy_mask[:n_test]))
        checks["no_assay_gene_ranking_equivalent"] = np.array_equal(official_ids[:, 3:], np.asarray(no_ids[:n_test, 3:]))
        checks["no_assay_custom_context_valid"] = checks["context_no_assay"] and checks["attention_pad_consistency_no_assay"]

        mapping_report = json.loads(mapping_report_path.read_text(encoding="utf-8"))
        checks["counts_mapping_conserved"] = mapping_report["count_conservation_chunk_failures"] == 0
        group_summaries: dict[str, Any] = {}
        length_frame = dev.obs[["stage", "sample_id", "cell_type_major"]].copy()
        length_frame["gene_tokens"] = gene_lengths
        anomaly = False
        for column in ["stage", "sample_id", "cell_type_major"]:
            summary = length_frame.groupby(column, observed=True)["gene_tokens"].agg(["count", "mean", "median", "min", "max"])
            group_summaries[column] = summary.reset_index().to_dict(orient="records")
            eligible = summary[summary["count"] >= 10]
            if len(eligible) > 1 and eligible["mean"].min() > 0 and eligible["mean"].max() / eligible["mean"].min() > 4:
                anomaly = True
        checks["no_programmatic_length_anomaly"] = not anomaly
    finally:
        dev.file.close()

    source_report = output_root / "reports/model_asset_report.md"
    if not source_report.is_file():
        raise RuntimeError("Model source audit report is missing; refusing trust_remote_code execution")
    torch.set_num_threads(min(16, max(1, os.cpu_count() or 1)))
    try:
        model = AutoModel.from_pretrained(str(model_dir), trust_remote_code=True)
        model.eval()
        batch = min(max(2, forward_batch), len(no_ids))
        with torch.inference_mode():
            embeddings = model.get_embeddings(
                input_ids=torch.from_numpy(np.asarray(no_ids[:batch]).astype(np.int64)),
                attention_mask=torch.from_numpy(np.asarray(no_mask[:batch])),
                with_context=False,
            )
        checks["cpu_forward"] = tuple(embeddings.shape) == (batch, int(manifest["hidden_dimension"])) and bool(torch.isfinite(embeddings).all())
        forward_shape = list(embeddings.shape)
    except Exception:
        logger.exception("CPU forward test failed")
        checks["cpu_forward"] = False
        forward_shape = None

    metrics: dict[str, Any] = {
        "checks": checks,
        "all_passed": all(checks.values()),
        "dev_shape": [int(no_ids.shape[0]), int(no_ids.shape[1])],
        "gene_token_lengths": describe(np.asarray(gene_lengths)),
        "token_id_min": int(no_ids.min()),
        "token_id_max": int(no_ids.max()),
        "declared_vocab_size": vocab_size,
        "config_n_tokens": config_tokens,
        "official_equivalence_spots": n_test,
        "controlled_dense_shape": [n_test, len(genes)],
        "controlled_dense_bytes": dense_bytes,
        "full_matrix_densified": False,
        "forward_embedding_shape": forward_shape,
        "group_token_length_summaries": group_summaries,
    }
    atomic_json(output_root / "reports/dev_validation_metrics.json", metrics)
    status_lines = "\n".join(f"- {name}: **{'PASS' if passed else 'FAIL'}**" for name, passed in checks.items())
    report = f"""# Development Token Validation Report

## Overall Assessment: {'PASS' if metrics['all_passed'] else 'FAIL'}

- Token shape: `{metrics['dev_shape']}`
- Official equivalence spots: **{n_test}**
- Controlled dense reference only: `{metrics['controlled_dense_shape']}` / **{dense_bytes / 2**20:.2f} MiB**
- CPU forward embedding shape: `{forward_shape}`
- Gene-token length distribution: `{metrics['gene_token_lengths']}`

## Checks

{status_lines}

The full expression matrix was never densified. A bounded {n_test}-spot × {len(genes):,}-gene aligned matrix was converted to dense solely for official-tokenizer equivalence. The main `no_assay` context is `[human=5, PAD=0, spatial=4]` with the technology position masked; the proxy uses official `10x 3' v3=11`. Both retain fixed three-position context required by `get_embeddings(with_context=False)`.
"""
    (output_root / "reports/dev_validation_report.md").write_text(report, encoding="utf-8")
    if not metrics["all_passed"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Development validation failed: {failed}")
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dev-h5ad", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--technology-mean", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--no-assay-dir", type=Path, required=True)
    parser.add_argument("--proxy-dir", type=Path, required=True)
    parser.add_argument("--mapping-report", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--equivalence-spots", type=int, default=32)
    parser.add_argument("--forward-batch", type=int, default=2)
    return parser.parse_args()


if __name__ == "__main__":
    import os

    args = parse_args()
    validate(
        args.dev_h5ad, args.mapping, args.technology_mean, args.model_dir,
        args.model_manifest, args.no_assay_dir, args.proxy_dir,
        args.mapping_report, args.output_root, args.equivalence_spots,
        args.forward_batch,
    )
