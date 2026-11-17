# RIF implantation: Nicheformer embeddings and neighborhood analysis

Reproducible code for the `pro_sec_with_RIF.h5ad` normal + recurrent implantation
failure (RIF) dataset. The requested directory name
`Script_10_RIF_implantation_Nitchformer_network` is retained; the model's official
name is **Nicheformer**.

## Scope and organization

| Path | Purpose |
| --- | --- |
| `scripts/00_*.py`–`08_*.py` | Original RIF preprocessing, model audit, mapping, means, tokens, validation, export and summary |
| `scripts/common.py`, `gene_utils.py` | Shared logging, checksums and sparse matrix helpers |
| `scripts/download_annotation.py` | Verified Ensembl release 112 annotation download |
| `scripts/run_pipeline.py` | Portable ordered preprocessing runner with stage checkpoints |
| `scripts/09_infer_embeddings.py` | Added CPU/CUDA embedding inference and batch resume |
| `scripts/10_neighbors_umap.py` | kNN graph and UMAP from the original exploratory notebook |
| `tests/` | Synthetic sparse references, sampling, tokenization and resume tests |
| `environment/` | Portable environment, reference dependency lock and optional analysis dependencies |
| `docs/` | Methods, provenance, expected results and GitHub release notes |
| `provenance/` | Source-code checksums and fixed resource hashes |

The original `00`–`08` scripts and helpers are copied unchanged from
`nicheformer_embedding_with_RIF/scripts/`. New entry points organize their
execution. The historical full-cohort GPU embedding command was not present in
that project; script `09` makes the audited official inference interface runnable.
Script `10` preserves the notebook's k=15, Euclidean metric and UMAP seed 0.
Here, “network” means a **spot embedding kNN graph**. This code does not claim to
reconstruct a gene regulatory network or spatial-adjacency graph.

## Input data

Obtain `pro_sec_with_RIF.h5ad` separately and keep it outside the code directory.
The exact reference input has SHA256
`2960e65257c74312b196e172e6226c9a632cc613c16ce0b62a1f06254ca2544b`,
146,414 spots × 38,299 genes, 24 samples and 21 donors. It includes 108,643 normal
and 37,771 RIF spots. See [expected results](docs/EXPECTED_RESULTS.md).

Required AnnData contents include sparse raw integer `layers["counts"]`, `X`,
`obsm["spatial"]`, the original gene annotation columns, and the metadata schema
documented in [METHODS.md](docs/METHODS.md). All expression tokenization reads
`layers["counts"]`. The full expression matrix is never made dense; the official
tokenizer comparison densifies only 32 aligned spots.

The source `condition` is uniformly `unspecified`. Derived `analysis_condition`
labels require concordant `RIF_` prefixes in `stage` and `sample_id`; both fields
must agree or the pipeline stops. Neither source H5AD nor its labels are rewritten.

The input H5AD, donor-level metadata, model weights, downloaded annotations,
generated tokens and generated results are excluded from this code release.
An input data accession/download link must be supplied with the paper to make
the entire workflow reproducible. This package starts from the prepared H5AD;
upstream construction of that H5AD is outside its scope.

## Install an isolated CPU environment

Run commands from this directory. Use mamba or micromamba equivalently to conda.
Python 3.10 is the reference version; do not install into an existing base env.

```bash
conda env create -f environment/environment.yml
conda activate rif_nicheformer
python --version
python -m pip install torch==2.4.1+cpu --index-url https://download.pytorch.org/whl/cpu
python -m pip install -r environment/requirements-lock.txt
python -m pytest -q tests
```

`requirements-lock.txt` records the original resolved versions. Two unportable
`file://` references were converted to their exact version pins. The original
Linux conda dependency export is preserved in
`environment/original-linux-environment.yml` (without the server-specific prefix).
It is an audit record; the small portable `environment.yml` is the install entry.
The complete environment was reused for validation on this server; installation
on a brand-new machine has not been tested during this packaging step.

## Reproduce preprocessing

```bash
python scripts/run_pipeline.py \
  --input-h5ad ../pro_sec_with_RIF.h5ad \
  --output-root outputs/rif_reference \
  --chunk-size 256 --seed 20260806
```

The runner creates output directories, verifies/downloads fixed resources,
performs all preprocessing, runs tests, validates the development set against
the official tokenizer and a CPU forward pass, then tokenizes the full cohort.
It creates the Windows export and final report. Commands and stdout/stderr are
saved under `outputs/rif_reference/logs/`; stage identity, product hashes and
completion records are saved in `pipeline_state.json`.

Download access is required for the official Hugging Face model and Ensembl GTF.
An existing annotation can be supplied with `--annotation-gtf /path/to/file.gtf.gz`;
the same fixed SHA256 is always enforced. The model revision is
`theislab/Nicheformer@0ba4ba162f9ad3d6c8b713761954c31be2fa69ed`.

```bash
# Print every command without creating outputs.
python scripts/run_pipeline.py --input-h5ad ../pro_sec_with_RIF.h5ad \
  --output-root outputs/rif_reference --dry-run

# Run only through development validation; resume later with the second command.
python scripts/run_pipeline.py --input-h5ad ../pro_sec_with_RIF.h5ad \
  --output-root outputs/rif_reference --stop-after validate
python scripts/run_pipeline.py --input-h5ad ../pro_sec_with_RIF.h5ad \
  --output-root outputs/rif_reference --resume
```

Use a **new empty output directory** for the runner. Its resume checkpoints apply
to runs created by this runner, not to the historical output folder. Resume
rejects changed input, code, key runtime versions, seed, chunk size or recorded
products. Individual tokenization scripts retain `--resume`, `--start-chunk`,
`--end-chunk` and explicit `--overwrite`; use their `--help` for manual recovery.
Do not run concurrent writers against the same output directory.

## Extract per-spot embeddings

On Linux CPU (slow for the full cohort), test two spots first:

```bash
python scripts/09_infer_embeddings.py \
  --token-dir outputs/rif_reference/tokenized/dev_no_assay \
  --model-dir outputs/rif_reference/model_assets/Nicheformer \
  --model-manifest outputs/rif_reference/model_assets/model_download_manifest.json \
  --output-dir outputs/cpu_smoke --device cpu --batch-size 2 --max-spots 2
```

For full CUDA inference, use a separate GPU environment with Python 3.10 and an
appropriate CUDA-enabled **torch 2.4.1** build, followed by
`python -m pip install -r environment/requirements-gpu.txt`.
`requirements-gpu.txt` pins direct dependencies, rather than claiming a complete
platform-independent CUDA lock. Review the downloaded Python code before execution;
the inference script checks its SHA256 against the audited model manifest.

```bash
python scripts/09_infer_embeddings.py \
  --token-dir outputs/rif_reference/tokenized/full_no_assay \
  --model-dir outputs/rif_reference/model_assets/Nicheformer \
  --model-manifest outputs/rif_reference/model_assets/model_download_manifest.json \
  --output-dir outputs/full_embeddings --device cuda --batch-size 8
```

Add `--resume` with the same parameters after interruption. Output is a float32
`embeddings.npy` of shape `(146414, 512)`, `obs_names.npy` and a batch-checksummed
`manifest.json`. On Windows, copy this code directory and `windows_export/`;
replace `--token-dir` with `windows_export/full_no_assay`. Download the exact model
revision shown in `windows_export/WINDOWS_INFERENCE_README.md` and use its model
manifest. Full GPU inference remains to be run on the target GPU machine.

## Reproduce the embedding neighborhood graph and UMAP

```bash
# Optional, only needed for interactive HTML.
python -m pip install plotly==5.24.1
python scripts/10_neighbors_umap.py \
  --embedding-dir outputs/full_embeddings \
  --obs-parquet outputs/rif_reference/metadata/obs.parquet \
  --spatial-npy outputs/rif_reference/metadata/spatial.npy \
  --output-dir outputs/neighbors_umap --n-neighbors 15 --seed 0 \
  --interactive-html --color sample_id
```

The resulting `nicheformer_graph.h5ad` contains `X_nicheformer`, `X_umap`, sparse
neighbor graphs and metadata, with a zero-column expression matrix. It does not
duplicate the source RNA matrix. Coordinates and metadata are matched by exact
spot order. No sample, donor, condition, stage, cell type or coordinate is encoded
into Nicheformer input tokens. Data-dependent batch effects can still be present
in expression and embeddings; no batch correction is performed here.
