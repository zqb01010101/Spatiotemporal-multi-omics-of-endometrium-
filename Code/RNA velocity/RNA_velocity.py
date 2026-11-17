import os
from pathlib import Path
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scvelo as scv
import scanpy.external as sce
import numpy as np
import harmonypy as hm



scv.settings.verbosity = 3
scv.set_figure_params("scvelo")


adata = sc.read_h5ad(
    '/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/RNA_velocity/adata_scvelo_stochastic.h5ad'
)

adata_v = adata.copy()
adata_v = adata_v[adata_v.obs["super_cell_type"].astype(str) == "stroma"].copy()
adata_v = adata_v[adata_v.obs['stage'].astype(str).isin(['sec_PO11'])].copy()

seurat_umap = pd.read_csv(
    "sec_late_stroma_harmony_umap_with_cell_id.csv",
    index_col=0
)

seurat_harmony = pd.read_csv(
    "sec_late_stroma_harmony_embedding_with_cell_id.csv",
    index_col=0
)

adata_v.obs["cell_id"] = adata_v.obs["cell_id"].astype(str)
seurat_umap["cell_id"] = seurat_umap["cell_id"].astype(str)
seurat_harmony["cell_id"] = seurat_harmony["cell_id"].astype(str)

seurat_umap = seurat_umap.set_index("cell_id")
seurat_harmony = seurat_harmony.set_index("cell_id")

common_ids = (
    pd.Index(adata_v.obs["cell_id"])
    .intersection(seurat_umap.index)
    .intersection(seurat_harmony.index)
)

keep_mask = adata_v.obs["cell_id"].isin(common_ids)

adata_v = adata_v[keep_mask].copy()
cell_id_order = adata_v.obs["cell_id"].astype(str).values
seurat_umap_aligned = seurat_umap.loc[cell_id_order]
seurat_harmony_aligned = seurat_harmony.loc[cell_id_order]

umap_cols = seurat_umap_aligned.select_dtypes(include=[np.number]).columns[:2]

print("UMAP columns:", umap_cols.tolist())

adata_v.obsm["X_umap"] = seurat_umap_aligned[umap_cols].to_numpy()

harmony_cols = seurat_harmony_aligned.select_dtypes(include=[np.number]).columns

print("Harmony dimensions:", len(harmony_cols))
print("Harmony columns head:", harmony_cols[:5].tolist())

adata_v.obsm["X_pca_harmony"] = seurat_harmony_aligned[harmony_cols].to_numpy()

print("adata_v:", adata_v.shape)
print("X_umap:", adata_v.obsm["X_umap"].shape)
print("X_pca_harmony:", adata_v.obsm["X_pca_harmony"].shape)

adata_v.obsm["X_pca_harmony_30"] = adata_v.obsm["X_pca_harmony"][:, :30]

sc.pp.neighbors(
    adata_v,
    n_neighbors=30,
    use_rep="X_pca_harmony_30"
)
scv.pp.moments(
    adata_v,
    n_pcs=None,
    n_neighbors=None
)
scv.tl.velocity(
    adata_v,
    mode="stochastic"
)

scv.tl.velocity_graph(adata_v)

scv.tl.terminal_states(adata_v)
scv.pl.scatter(
    adata_v,
    basis="umap",
    color=["exclusive_cell_type_4", "root_cells", "end_points"],
    cmap="viridis"
)


scv.tl.terminal_states(adata_v)
scv.pl.scatter(
    adata_v,
    basis="umap",
    color=["exclusive_cell_type_4", "root_cells", "end_points"],
    cmap="viridis"
)