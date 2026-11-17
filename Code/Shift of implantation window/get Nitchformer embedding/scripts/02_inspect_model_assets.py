#!/usr/bin/env python3
"""Download a pinned Nicheformer snapshot and audit code/config without executing remote code."""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from huggingface_hub import HfApi, snapshot_download
from safetensors import safe_open

from common import atomic_json, setup_logging, sha256_file


DEFAULT_REVISION = "0ba4ba162f9ad3d6c8b713761954c31be2fa69ed"


def literal_assignments(source: str) -> dict[str, Any]:
    """Collect safe literal module/class assignments from Python source."""
    tree = ast.parse(source)
    values: dict[str, Any] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            value_node = node.value
            for target in targets:
                if isinstance(target, ast.Name):
                    try:
                        values[target.id] = ast.literal_eval(value_node)
                    except (ValueError, TypeError):
                        pass
    return values


def source_excerpt(source: str, pattern: str, radius: int = 12) -> str:
    """Return numbered lines around the first source pattern."""
    lines = source.splitlines()
    match = next((i for i, line in enumerate(lines) if pattern in line), None)
    if match is None:
        return f"Pattern `{pattern}` not found."
    start, end = max(0, match - radius), min(len(lines), match + radius + 1)
    return "\n".join(f"{i + 1:4d}: {lines[i]}" for i in range(start, end))


def infer_identifier_type(values: list[str]) -> dict[str, Any]:
    """Classify model gene identifiers using explicit regex rates."""
    ens = re.compile(r"^ENSG\d+(?:\.\d+)?$")
    ensembl_rate = sum(bool(ens.fullmatch(x)) for x in values) / max(len(values), 1)
    symbol_rate = sum(bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", x)) for x in values) / max(len(values), 1)
    label = "ensembl_gene_id" if ensembl_rate >= 0.8 else "gene_symbol_or_name"
    return {"inferred_type": label, "ensembl_rate": ensembl_rate, "symbol_like_rate": symbol_rate}


def inspect_assets(repo_id: str, revision: str, output_root: Path) -> dict[str, Any]:
    """Download exact revision, hash assets, and statically inspect model rules."""
    logger = setup_logging(output_root / "logs/02_inspect_model_assets.log")
    destination = output_root / "model_assets/Nicheformer"
    destination.mkdir(parents=True, exist_ok=True)
    info = HfApi().model_info(repo_id, revision=revision, files_metadata=True)
    resolved = info.sha
    if resolved != revision:
        raise RuntimeError(f"Requested revision resolved to {resolved}, expected {revision}")
    snapshot_download(repo_id=repo_id, revision=resolved, local_dir=destination)
    required = ["config.json", "model.safetensors", "vocab.json", "model.h5ad", "tokenization_nicheformer.py", "modeling_nicheformer.py"]
    missing = [name for name in required if not (destination / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Downloaded snapshot is missing required files: {missing}")

    file_hashes = {str(path.relative_to(destination)): sha256_file(path) for path in sorted(destination.rglob("*")) if path.is_file()}
    config = json.loads((destination / "config.json").read_text(encoding="utf-8"))
    vocab = json.loads((destination / "vocab.json").read_text(encoding="utf-8"))
    tokenizer_source = (destination / "tokenization_nicheformer.py").read_text(encoding="utf-8")
    modeling_source = (destination / "modeling_nicheformer.py").read_text(encoding="utf-8")
    literals = literal_assignments(tokenizer_source)

    with safe_open(destination / "model.safetensors", framework="np") as handle:
        tensor_shapes = {name: list(handle.get_slice(name).get_shape()) for name in handle.keys()}
    parameter_count = int(sum(np.prod(shape, dtype=np.int64) for shape in tensor_shapes.values()))

    model_data = ad.read_h5ad(destination / "model.h5ad", backed="r")
    try:
        model_genes = model_data.var_names.astype(str).tolist()
        model_h5ad_shape = list(model_data.shape)
        model_var_columns = model_data.var.columns.astype(str).tolist()
    finally:
        model_data.file.close()
    id_info = infer_identifier_type(model_genes)
    gene_vocab = {str(key): int(value) for key, value in vocab.items() if str(key) in set(model_genes)}
    token_ids = list(vocab.values()) if isinstance(vocab, dict) else []
    pad_id = literals.get("PAD_TOKEN", vocab.get("[PAD]", vocab.get("", 0)))
    mask_id = literals.get("MASK_TOKEN", vocab.get("[MASK]", 1))
    modality = literals.get("MODALITY_DICT", literals.get("modality_dict", {}))
    species = literals.get("SPECIES_DICT", literals.get("SPECIE_DICT", literals.get("specie_dict", {})))
    technology = literals.get("TECHNOLOGY_DICT", literals.get("technology_dict", {}))
    max_length = config.get("context_length", config.get("max_seq_len", config.get("max_position_embeddings", config.get("max_length"))))
    if max_length is None:
        raise ValueError("No maximum sequence/context length found in config.json")
    hidden_dim = config.get("dim_model", config.get("hidden_size", config.get("d_model")))

    manifest: dict[str, Any] = {
        "repo_id": repo_id,
        "requested_revision": revision,
        "resolved_revision": resolved,
        "download_path": str(destination.resolve()),
        "files": file_hashes,
        "config": config,
        "parameter_count": parameter_count,
        "tensor_count": len(tensor_shapes),
        "vocab_size": len(vocab),
        "vocab_token_id_min": int(min(token_ids)) if token_ids else None,
        "vocab_token_id_max": int(max(token_ids)) if token_ids else None,
        "model_h5ad_shape": model_h5ad_shape,
        "model_gene_count": len(model_genes),
        "model_var_columns": model_var_columns,
        "gene_identifiers": id_info,
        "model_genes_found_in_vocab": len(gene_vocab),
        "pad_token_id": pad_id,
        "mask_token_id": mask_id,
        "modality_tokens": modality,
        "species_tokens": species,
        "technology_tokens": technology,
        "max_sequence_length": max_length,
        "hidden_dimension": hidden_dim,
        "trust_remote_code_executed": False,
        "source_audit_completed_before_remote_code": True,
    }
    atomic_json(output_root / "model_assets/model_download_manifest.json", manifest)

    report = f"""# Nicheformer Model Asset Report

## Pinned Snapshot

- Repository: `{repo_id}`
- Immutable revision: `{resolved}`
- Parameter count from safetensors shapes: **{parameter_count:,}**
- Vocabulary size / token-ID range: **{len(vocab):,} / {manifest['vocab_token_id_min']}–{manifest['vocab_token_id_max']}**
- Model gene count: **{len(model_genes):,}**
- Identifier inference: `{id_info}`
- Model H5AD shape: `{model_h5ad_shape}`
- Maximum sequence length: **{max_length}**
- Hidden dimension: **{hidden_dim}**

## Context and Special Tokens

- PAD: `{pad_id}`; MASK: `{mask_id}`
- Modality: `{modality}`
- Species: `{species}`
- Technology/assay: `{technology}`

No downloaded Python source was executed and `trust_remote_code=True` was not used during this audit. The following excerpts are preserved from the pinned files so context construction and embedding behavior are decided from source rather than memory.

## Tokenizer Context Excerpt

```python
{source_excerpt(tokenizer_source, 'technology', 28)}
```

## Tokenization Excerpt

```python
{source_excerpt(tokenizer_source, 'tokenize', 32)}
```

## `get_embeddings()` Excerpt

```python
{source_excerpt(modeling_source, 'get_embeddings', 45)}
```

## `with_context=False` Interpretation

The exact implementation is shown above. Validation code must confirm its context slicing and pooling semantics with a CPU forward pass before full tokenization.

## Required SHA256 Values

- `model.safetensors`: `{file_hashes['model.safetensors']}`
- `vocab.json`: `{file_hashes['vocab.json']}`
- `model.h5ad`: `{file_hashes['model.h5ad']}`
- `tokenization_nicheformer.py`: `{file_hashes['tokenization_nicheformer.py']}`
- `modeling_nicheformer.py`: `{file_hashes['modeling_nicheformer.py']}`
"""
    (output_root / "reports/model_asset_report.md").write_text(report, encoding="utf-8")
    logger.info("Downloaded and audited revision %s", resolved)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-id", default="theislab/Nicheformer")
    parser.add_argument("--revision", default=DEFAULT_REVISION)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    inspect_assets(args.repo_id, args.revision, args.output_root)
