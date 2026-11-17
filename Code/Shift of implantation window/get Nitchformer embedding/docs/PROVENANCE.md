# Code provenance and organization

Packaged on 2026-10-21 from the completed `pro_sec_with_RIF.h5ad` analysis.

- `scripts/00_inspect_h5ad.py` through `08_final_summary.py`, `common.py` and
  `gene_utils.py` are byte-identical copies of the RIF run's source files.
- The seven original `test_*.py` modules and `conftest.py` are copied unchanged.
- `provenance/source_inventory.json` records original relative locations, SHA256,
  sizes, and the source notebook hash. Its inventory can verify the copies.
- `environment/requirements-lock.txt` retains the successful run's dependency
  versions, replacing local packaging/PySocks build paths with version pins.
- `environment/original-linux-environment.yml` is the historical conda export with
  its server-specific environment prefix removed.
- `download_annotation.py` adds a fixed URL/hash download, eliminating the previous
  manual annotation-copy step.
- `run_pipeline.py` supplies directory creation, execution order, validation gate,
  logs and SHA256-based stage recovery. It does not change the scientific
  preprocessing implementation.
- `09_infer_embeddings.py` is a new executable implementation of the pinned model's
  audited last-layer embedding interface. Historical full GPU inference code was
  not available among the source scripts; identical historical GPU results cannot
  be claimed solely from this packaging step.
- `10_neighbors_umap.py` adapts code cells 4 and 7 of `Untitled.ipynb`: kNN=15,
  Euclidean metric, UMAP seed=0 and optional WebGL scatter colored by sample.
  The notebook's 34 MB rendered outputs were not copied.
- No original H5AD, per-spot metadata, model binary, token arrays or generated
  analysis result is copied into this source directory.

The directory title requested by the author contains `Nitchformer`; model names
and Python documentation consistently use official `Nicheformer` spelling.
No GitHub upload or author commit was performed by the organization task.

See `VALIDATION.md` for checks actually run on the organized package.
