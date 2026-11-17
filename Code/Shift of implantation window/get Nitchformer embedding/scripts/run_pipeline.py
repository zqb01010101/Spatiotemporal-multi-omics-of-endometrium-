#!/usr/bin/env python3
"""Run the paper's RIF preprocessing workflow with verified stage checkpoints."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

from common import atomic_json, runtime_info, setup_logging, sha256_file
from download_annotation import FILENAME


SCRIPTS = Path(__file__).resolve().parent
PACKAGE = SCRIPTS.parent
REVISION = "0ba4ba162f9ad3d6c8b713761954c31be2fa69ed"


def stage_commands(args: argparse.Namespace) -> list[tuple[str, list[str], list[str]]]:
    """Return ordered commands and checkpoint products, using configurable paths."""
    out = args.output_root.resolve()
    inp = args.input_h5ad.resolve()
    model = out / "model_assets/Nicheformer"
    manifest = out / "model_assets/model_download_manifest.json"
    mapping = out / "mapping/input_to_nicheformer_gene.tsv"
    means = out / "technology_means/dbit_global_nonzero_mean.npy"
    dev = out / "dev_subset/dev_spots.h5ad"
    annotation = args.annotation_gtf.resolve() if args.annotation_gtf else out / "model_assets/annotations" / FILENAME
    base = ["--output-root", str(out)]
    source = ["--input-h5ad", str(inp), *base]
    chunk = ["--chunk-size", str(args.chunk_size)]

    def command(script: str, options: list[str]) -> list[str]:
        return [sys.executable, str(SCRIPTS / script), *options]

    def tokens(name: str, h5ad: Path, strategy: str, full: bool = False) -> tuple[str, list[str], list[str]]:
        directory = out / "tokenized" / name
        options = ["--input-h5ad", str(h5ad), "--mapping", str(mapping),
                   "--technology-mean", str(means), "--model-dir", str(model),
                   "--model-manifest", str(manifest), "--output-dir", str(directory),
                   "--strategy", strategy, *base, *chunk]
        if args.resume and (directory / "manifest.json").is_file():
            options.append("--resume")
        products = [f"tokenized/{name}/{file}" for file in
                    ["manifest.json", "input_ids.npy", "attention_mask.npy", "obs_names.npy"]]
        if full:
            options.append("--export-metadata")
            products.extend(["metadata/obs.parquet", "metadata/spatial.npy"])
        return name, command("05_tokenize_nicheformer.py", options), products

    return [
        ("inspect", command("00_inspect_h5ad.py", [*source, *chunk]),
         ["reports/input_h5ad_metrics.json", "reports/input_h5ad_report.md", "metadata/input_obs.parquet", "metadata/input_var.parquet"]),
        ("subset", command("01_create_dev_subset.py", [*source, "--seed", str(args.seed), "--target", "7200"]),
         ["dev_subset/dev_spots.h5ad", "dev_subset/dev_selection_manifest.json", "dev_subset/dev_spot_ids.txt", "reports/dev_subset_report.md"]),
        ("model", command("02_inspect_model_assets.py", [*base, "--revision", REVISION]),
         ["model_assets/model_download_manifest.json", "reports/model_asset_report.md", *[
             f"model_assets/Nicheformer/{name}" for name in ["config.json", "model.safetensors", "model.h5ad", "vocab.json", "__init__.py", "configuration_nicheformer.py", "tokenization_nicheformer.py", "modeling_nicheformer.py", "masking.py"]]]),
        ("annotation", command("download_annotation.py", [*base, *(["--local-gtf", str(annotation)] if args.annotation_gtf else [])]),
         ["model_assets/annotations/manifest.json"]),
        ("mapping", command("03_build_gene_mapping.py", [*source, *chunk, "--model-dir", str(model), "--annotation-gtf", str(annotation)]),
         ["reports/gene_mapping_report.json", "reports/gene_mapping_report.md", "mapping/input_to_nicheformer_gene.tsv", "mapping/duplicate_mapping_groups.tsv", "mapping/unmapped_input_genes.tsv", "mapping/model_genes_not_present.tsv", "metadata/var_mapping.parquet"]),
        ("mean", command("04_compute_dbit_nonzero_mean.py", [*source, *chunk, "--mapping", str(mapping), "--model-manifest", str(manifest)]),
         ["technology_means/dbit_global_nonzero_mean.npy", "technology_means/dbit_global_nonzero_count.npy", "technology_means/dbit_global_nonzero_mean_metadata.json"]),
        tokens("dev_no_assay", dev, "no_assay"),
        tokens("dev_10x3v3_proxy", dev, "10x3v3_proxy"),
        ("tests", [sys.executable, "-m", "pytest", "-q", str(PACKAGE / "tests")], []),
        ("validate", command("06_validate_tokens.py", ["--dev-h5ad", str(dev), "--mapping", str(mapping), "--technology-mean", str(means), "--model-dir", str(model), "--model-manifest", str(manifest), "--no-assay-dir", str(out / "tokenized/dev_no_assay"), "--proxy-dir", str(out / "tokenized/dev_10x3v3_proxy"), "--mapping-report", str(out / "reports/gene_mapping_report.json"), *base]),
         ["reports/dev_validation_metrics.json", "reports/dev_validation_report.md"]),
        tokens("full_no_assay", inp, "no_assay", full=True),
        ("export", command("07_prepare_windows_export.py", [*base, "--project-root", str(inp.parent)]),
         ["reports/windows_export_manifest.json", "windows_export/SHA256SUMS.txt", "nicheformer_windows_export.tar"]),
        ("summary", command("08_final_summary.py", [*base, "--project-root", str(inp.parent), "--input-h5ad", str(inp)]),
         ["reports/final_report.md"]),
    ]


def verify_products(root: Path, products: dict[str, str]) -> None:
    """Reject missing or changed checkpoint outputs before skipping a stage."""
    for relative, expected in products.items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Checkpoint product changed or is missing: {path}")


def run(args: argparse.Namespace) -> None:
    """Execute all stages; a failed validation prevents full tokenization."""
    stages = stage_commands(args)
    if args.stop_after and args.stop_after not in {name for name, _, _ in stages}:
        raise ValueError(f"Unknown --stop-after {args.stop_after!r}")
    if args.dry_run:
        for name, command, _ in stages:
            print(f"[{name}] {shlex.join(command)}")
            if name == args.stop_after:
                break
        return
    if sys.version_info[:2] not in {(3, 9), (3, 10)}:
        raise RuntimeError("Use an isolated Python 3.10 (preferred) or 3.9 environment")
    if not args.input_h5ad.is_file():
        raise FileNotFoundError(args.input_h5ad)
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive")
    root = args.output_root.resolve()
    state_path = root / "pipeline_state.json"
    if root.exists() and any(root.iterdir()) and not (args.resume and state_path.is_file()):
        raise FileExistsError("Choose an empty output directory, or --resume a run created by this runner")
    for directory in ["logs", "reports", "metadata", "dev_subset", "mapping", "technology_means", "tokenized", "environment"]:
        (root / directory).mkdir(parents=True, exist_ok=True)
    identity: dict[str, Any] = {
        "input_h5ad": str(args.input_h5ad.resolve()), "input_sha256": sha256_file(args.input_h5ad),
        "revision": REVISION, "chunk_size": args.chunk_size, "seed": args.seed,
        "annotation_sha256": sha256_file(args.annotation_gtf) if args.annotation_gtf else None,
        "code_sha256": {str(p.relative_to(PACKAGE)): sha256_file(p) for p in
                        sorted([*SCRIPTS.glob("*.py"), *(PACKAGE / "tests").glob("*.py")])},
        "packages": {name: importlib.metadata.version(name) for name in
                     ["numpy", "scipy", "anndata", "pandas", "torch", "transformers", "huggingface_hub"]},
    }
    state = json.loads(state_path.read_text()) if args.resume and state_path.is_file() else {"identity": identity, "stages": {}}
    if state["identity"] != identity:
        raise RuntimeError("Input, code, parameters, or runtime changed; use a new output directory")
    atomic_json(state_path, state)
    atomic_json(root / "environment/runtime.json", {"runtime": runtime_info(), **identity})
    logger = setup_logging(root / "logs/pipeline.log")
    for name, command, products in stages:
        if name in state["stages"]:
            verify_products(root, state["stages"][name])
            logger.info("Verified checkpoint %s", name)
        else:
            if name == "full_no_assay":
                validation = json.loads((root / "reports/dev_validation_metrics.json").read_text())
                if not validation.get("all_passed"):
                    raise RuntimeError("Development validation failed; full tokenization forbidden")
            logger.info("Running %s: %s", name, shlex.join(command))
            with (root / f"logs/stage_{name}.log").open("a", encoding="utf-8") as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            state["stages"][name] = {p: sha256_file(root / p) for p in products}
            atomic_json(state_path, state)
        if name == args.stop_after:
            break


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--annotation-gtf", type=Path, help="Optional local Ensembl 112 GTF; checksum is enforced")
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260806)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--stop-after", help="Stop at a named stage, e.g. validate; resume later")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without reading data or creating outputs")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
