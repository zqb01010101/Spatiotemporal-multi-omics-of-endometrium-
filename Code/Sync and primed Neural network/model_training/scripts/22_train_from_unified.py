#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""22_train_from_unified.py — Train the v5 network on the unified pair table

Read peak_gene_label_app/data_cache/pair_features_all.csv:
  - Keep rows with all four scores.
  - Use only the 12 peak_* and gene_* features as inputs.
  - Omit metadata such as is_organoid, context, and batch.

Save the trained model to outputs/model/ and print validation metrics.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODEL_DIR = Path(__file__).resolve().parents[1]
APP_DIR = MODEL_DIR.parent / "peak_gene_label_app"
DATA_FILE = APP_DIR / "data_cache" / "pair_features_all.csv"
OUT_DIR = MODEL_DIR / "outputs"
MODEL_OUT = OUT_DIR / "model"
FIG_OUT = OUT_DIR / "figures"
MODEL_OUT.mkdir(parents=True, exist_ok=True)
FIG_OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(MODEL_DIR / "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("train_mod", MODEL_DIR / "scripts" / "02_train.py")
train_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train_mod)

STAGES = ["PO3", "PO5", "PO7", "PO9", "PO11", "PO13"]
FEATURES = [f"peak_{s}" for s in STAGES] + [f"gene_{s}" for s in STAGES]
SCORE_COLS = ["score_pre_open", "score_sync_open", "score_sync_repress", "score_unrelated"]


def main():
    df = pd.read_csv(DATA_FILE)
    scored = df.dropna(subset=SCORE_COLS).reset_index(drop=True)
    print("统一文件行数：", len(df), "；有完整打分的行数：", len(scored))

    X = scored[FEATURES].to_numpy(dtype=float)
    Y = scored[SCORE_COLS].to_numpy(dtype=float) / 10.0
    hard = scored.get("hard_label")

    from sklearn.model_selection import train_test_split
    strat = None
    if hard is not None:
        strat = hard.astype(str)
    X_tr, X_val, Y_tr, Y_val = train_test_split(
        X, Y, test_size=0.2, random_state=42, stratify=strat)
    print("train/val:", len(Y_tr), "/", len(Y_val))

    model = train_mod.ScoreMLP(
        X.shape[1], train_mod.HIDDEN, Y.shape[1], seed=42,
        keep_prob=train_mod.DROPOUT_KEEP,
        loss_weights=train_mod.OUTPUT_LOSS_WEIGHTS)
    train_mod.train(model, X_tr, Y_tr, X_val, Y_val)

    pred = model.forward(X_val) * 10.0
    y = Y_val * 10.0
    print("\n验证集（0-1 尺度）：")
    for i, c in enumerate(train_mod.CLASSES):
        rmse = float(np.sqrt(np.mean((y[:, i] - pred[:, i]) ** 2)))
        corr = float(np.corrcoef(y[:, i], pred[:, i])[0, 1])
        print(f"  {c:12s} RMSE={rmse:.4f} corr={corr:.4f}")
    np.savez(MODEL_OUT / "model.npz", **model.params())
    print("\n模型已保存到 outputs/model/")


if __name__ == "__main__":
    main()
