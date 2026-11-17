#!/usr/bin/env python3
"""Prepare a checksummed Windows inference export after successful validation."""

from __future__ import annotations

import argparse
import json
import shutil
import tarfile
from pathlib import Path
from typing import Any

from common import atomic_json, setup_logging, sha256_file


TOKEN_DIRS = ["full_no_assay", "dev_no_assay", "dev_10x3v3_proxy"]
TOKEN_FILES = ["input_ids.npy", "attention_mask.npy", "obs_names.npy", "manifest.json"]


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def prepare(output_root: Path, project_root: Path) -> dict[str, Any]:
    """Copy validated artifacts, hash them, and create an uncompressed tar archive."""
    logger = setup_logging(output_root / "logs/07_prepare_windows_export.log")
    validation = json.loads((output_root / "reports/dev_validation_metrics.json").read_text(encoding="utf-8"))
    if not validation.get("all_passed"):
        raise RuntimeError("Development validation did not pass; export is forbidden")
    full_manifest = json.loads((output_root / "tokenized/full_no_assay/manifest.json").read_text(encoding="utf-8"))
    if not full_manifest.get("complete"):
        raise RuntimeError("Full tokenization is incomplete; export is forbidden")

    export = output_root / "windows_export"
    export.mkdir(parents=True, exist_ok=True)
    for token_dir in TOKEN_DIRS:
        for filename in TOKEN_FILES:
            source = output_root / "tokenized" / token_dir / filename
            if not source.is_file():
                raise FileNotFoundError(source)
            copy_file(source, export / token_dir / filename)
    extras = {
        output_root / "metadata/obs.parquet": export / "obs.parquet",
        output_root / "metadata/spatial.npy": export / "spatial.npy",
        output_root / "technology_means/dbit_global_nonzero_mean.npy": export / "dbit_global_nonzero_mean.npy",
        output_root / "reports/gene_mapping_report.md": export / "gene_mapping_report.md",
        output_root / "model_assets/model_download_manifest.json": export / "model_download_manifest.json",
    }
    for source, destination in extras.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        copy_file(source, destination)

    model_revision = full_manifest["model_revision"]
    readme = f"""# Windows Nicheformer Inference

This package contains pre-tokenized DBiT-seq data; it intentionally excludes the original H5AD and model weights.

## Exact model

- Repository: `theislab/Nicheformer`
- Revision: `{model_revision}`
- Main context: `[human=5, PAD/no-assay=0, spatial=4]`; technology attention is `False`.
- Sensitivity context: `[human=5, 10x 3' v3=11, spatial=4]` (development set only).

Download the identical model revision on Windows PowerShell:

```powershell
py -m pip install huggingface_hub torch transformers safetensors numpy pandas pyarrow
py -c "from huggingface_hub import snapshot_download; snapshot_download('theislab/Nicheformer', revision='{model_revision}', local_dir='model_assets/Nicheformer')"
Get-FileHash -Algorithm SHA256 .\\full_no_assay\\input_ids.npy
```

Load tokens and run a GPU batch:

```python
import numpy as np
import torch
from transformers import AutoModel

ids = np.load("full_no_assay/input_ids.npy", mmap_mode="r")
mask = np.load("full_no_assay/attention_mask.npy", mmap_mode="r")
model = AutoModel.from_pretrained("model_assets/Nicheformer", trust_remote_code=True).cuda().eval()
with torch.inference_mode():
    emb = model.get_embeddings(
        torch.as_tensor(ids[:32], dtype=torch.long, device="cuda"),
        torch.as_tensor(mask[:32], dtype=torch.bool, device="cuda"),
        with_context=False,
    )
print(emb.shape)
```

Inspect the pinned Python source before enabling `trust_remote_code`. Note that the official `get_embeddings()` mean-pools padded positions without applying an attention-weighted pooling denominator; preserve the validated behavior for reproducibility unless intentionally analyzing an alternative pooling method.
"""
    (export / "WINDOWS_INFERENCE_README.md").write_text(readme, encoding="utf-8")

    files = sorted(path for path in export.rglob("*") if path.is_file() and path.name != "SHA256SUMS.txt")
    checksum_lines = [f"{sha256_file(path)}  {path.relative_to(export).as_posix()}" for path in files]
    (export / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    archive = output_root / "nicheformer_windows_export.tar"
    with tarfile.open(archive, mode="w") as handle:
        handle.add(export, arcname="windows_export", recursive=True)
    manifest = {
        "export_directory": str(export.resolve()),
        "archive": str(archive.resolve()),
        "archive_sha256": sha256_file(archive),
        "archive_bytes": archive.stat().st_size,
        "file_count": len(files) + 1,
        "model_revision": model_revision,
        "source_h5ad_included": False,
        "transfer_command": f"rsync -avh --progress {export.resolve()}/ USER@WINDOWS_OR_TRANSFER_TARGET:/target/path/",
    }
    atomic_json(output_root / "reports/windows_export_manifest.json", manifest)
    logger.info("Prepared %s (%d bytes)", archive, archive.stat().st_size)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    prepare(args.output_root, args.project_root)
