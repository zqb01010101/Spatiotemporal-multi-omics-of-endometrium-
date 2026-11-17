#!/usr/bin/env python3
"""Reproduce the embedding kNN graph and UMAP parameters in Untitled.ipynb."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from common import atomic_json, setup_logging, sha256_file


def analyze(args: argparse.Namespace) -> None:
    """Attach aligned embeddings to metadata, compute kNN/UMAP and save results."""
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Choose a new empty graph output directory")
    if args.n_neighbors < 2:
        raise ValueError("--n-neighbors must be at least 2")
    embeddings = np.load(args.embedding_dir / "embeddings.npy", mmap_mode="r")
    names = np.load(args.embedding_dir / "obs_names.npy", allow_pickle=False).astype(str)
    manifest = json.loads((args.embedding_dir / "manifest.json").read_text())
    if not manifest.get("complete") or list(embeddings.shape) != manifest["identity"]["shape"]:
        raise ValueError("Embedding inference is not complete or shape differs from manifest")
    obs = pd.read_parquet(args.obs_parquet)
    obs.index = obs.index.astype(str)
    if not obs.index.is_unique or len(np.unique(names)) != len(names):
        raise ValueError("Observation names must be unique")
    if len(embeddings) != len(obs) or not np.array_equal(names, obs.index.to_numpy()):
        raise ValueError("Embedding and metadata row order must match exactly")
    if not np.isfinite(embeddings).all():
        raise ValueError("Embeddings contain nonfinite values")
    if len(obs) <= args.n_neighbors:
        raise ValueError("Number of spots must exceed --n-neighbors")
    args.output_dir.mkdir(parents=True)
    logger = setup_logging(args.output_dir / "neighbors_umap.log")
    # Expression is unnecessary: only the small spot x 512 embedding is materialized.
    data = ad.AnnData(X=sparse.csr_matrix((len(obs), 0), dtype=np.float32), obs=obs)
    data.obsm["X_nicheformer"] = np.array(embeddings, dtype=np.float32)
    sc.pp.neighbors(data, use_rep="X_nicheformer", n_neighbors=args.n_neighbors,
                    metric="euclidean", random_state=args.seed)
    sc.tl.umap(data, random_state=args.seed)
    if args.spatial_npy:
        spatial = np.load(args.spatial_npy, allow_pickle=False)
        if spatial.shape != (len(obs), 2):
            raise ValueError("Spatial coordinates must have shape (n_spots, 2)")
        data.obsm["spatial"] = spatial
    data.write_h5ad(args.output_dir / "nicheformer_graph.h5ad", compression="gzip")
    frame = obs.copy()
    frame[["UMAP1", "UMAP2"]] = data.obsm["X_umap"]
    frame.to_parquet(args.output_dir / "umap_metadata.parquet")
    atomic_json(args.output_dir / "manifest.json", {
        "embedding_sha256": sha256_file(args.embedding_dir / "embeddings.npy"),
        "obs_sha256": sha256_file(args.obs_parquet), "n_spots": len(obs),
        "n_neighbors": args.n_neighbors, "metric": "euclidean", "seed": args.seed,
        "scanpy_version": sc.__version__, "graph_type": "embedding-space kNN, not gene regulatory or spatial-adjacency network",
    })
    if args.interactive_html:
        import plotly.express as px
        if args.color not in frame.columns:
            raise KeyError(f"Missing plot color field {args.color}")
        frame["spot_id_plot"] = frame.index
        figure = px.scatter(frame, x="UMAP1", y="UMAP2", color=args.color,
                            hover_name="spot_id_plot", render_mode="webgl", opacity=0.7)
        figure.update_traces(marker={"size": 3})
        figure.update_layout(width=1000, height=800, title="Nicheformer UMAP")
        figure.write_html(args.output_dir / "umap.html", include_plotlyjs=True)
    logger.info("Saved kNN graph and UMAP for %d spots", len(obs))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding-dir", type=Path, required=True)
    parser.add_argument("--obs-parquet", type=Path, required=True)
    parser.add_argument("--spatial-npy", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--interactive-html", action="store_true")
    parser.add_argument("--color", default="sample_id")
    analyze(parser.parse_args())
