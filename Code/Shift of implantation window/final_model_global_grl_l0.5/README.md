# Window-of-implantation model

This workflow fits separate stromal and epithelial models on Nicheformer embeddings. It includes `pro_stage` with a temporal target of `-1`, uses a global gradient reversal layer for the method branch (maximum adversarial weight 0.5), and weights the temporal loss by 0.3. The analysis settings are in [config.yaml](config.yaml).

## Inputs and environment

Run from this directory. Set `input_h5ad` to the annotated H5AD containing the `X_nicheformer` embedding and `method_metadata` to the sample-level metadata table. Use an environment with PyTorch, AnnData, NumPy, pandas, scikit-learn, joblib, and PyYAML. The scripts also import shared modules from the sibling directory `../stage_prediction/src/`, which must be supplied separately. The commands below use the original Conda environment name.

## Commands

```bash
conda run -n nicheformer_gpu python scripts/01_train.py --cell-type stroma
conda run -n nicheformer_gpu python scripts/01_train.py --cell-type epithelial
conda run -n nicheformer_gpu python scripts/03_inference.py --cell-type stroma
conda run -n nicheformer_gpu python scripts/03_inference.py --cell-type epithelial
conda run -n nicheformer_gpu python scripts/02_evaluate.py --cell-type stroma
conda run -n nicheformer_gpu python scripts/02_evaluate.py --cell-type epithelial
conda run -n nicheformer_gpu python scripts/04_combine_celltypes.py
```

Training creates model artifacts under `models/`. Inference and evaluation create tables under `results/`. The combined analysis uses the stromal and epithelial prediction tables from the preceding steps. The repository includes the scripts and configuration for generating these outputs.
