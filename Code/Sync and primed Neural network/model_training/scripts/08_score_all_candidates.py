#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""08_score_all_candidates.py — Score all candidates with the v3 regression model

Outputs:
  outputs/all_candidate_scores.csv
      Columns: peak, gene, score_pre_open, score_sync_open, score_sync_repress, score_unrelated
      Scores range from 0 to 10. Pairs present in both contexts appear on two rows.
  outputs/all_candidate_scores_with_context.csv
      Full export with candidate_id and context
"""

import numpy as np
import pandas as pd
from pathlib import Path

MODEL_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = MODEL_DIR / "data"
MODEL_OUT = MODEL_DIR / "outputs" / "model"
OUT_DIR = MODEL_DIR / "outputs"

STAGES = ["PO3", "PO5", "PO7", "PO9", "PO11", "PO13"]
FEATURES = [f"peak_{s}" for s in STAGES] + [f"gene_{s}" for s in STAGES]
CLASSES = ["pre_open", "sync_open", "sync_repress", "unrelated"]
SCORE_COLS = [f"score_{c}" for c in CLASSES]


def relu(x):
    return np.maximum(0.0, x)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def main():
    feat = pd.read_csv(DATA_DIR / "all_candidates_features.csv")
    print("候选总数：", len(feat))
    X = feat[FEATURES].to_numpy(dtype=float)

    Xs = X    # Use raw fractions as in training

    w = np.load(MODEL_OUT / "model.npz")
    layer_keys = sorted([k for k in w.files if k.startswith("W")],
                        key=lambda k: int(k[1:]))
    z = Xs
    for i, wk in enumerate(layer_keys):
        z = z @ w[wk] + w["b" + wk[1:]]
        z = relu(z) if i < len(layer_keys) - 1 else sigmoid(z)
    out = z * 10.0     # Scores on the 0–10 scale

    for i, c in enumerate(SCORE_COLS):
        feat[c] = out[:, i]

    table = feat[["peak", "gene"] + SCORE_COLS].copy()
    table = table.sort_values(["peak", "gene"]).reset_index(drop=True)
    table.to_csv(OUT_DIR / "all_candidate_scores.csv", index=False)
    print("已写出：outputs/all_candidate_scores.csv（", len(table), " 行）")
    print(table.head(5).to_string(index=False))

    full = feat[["candidate_id", "context", "peak", "gene"] + SCORE_COLS].copy()
    full = full.sort_values(["context", "peak", "gene"]).reset_index(drop=True)
    full.to_csv(OUT_DIR / "all_candidate_scores_with_context.csv", index=False)
    print("已写出：outputs/all_candidate_scores_with_context.csv（", len(full), " 行）")

    print("\n预测类别分布（四个分数取最大）：")
    pred = np.array(CLASSES)[out.argmax(axis=1)]
    print(pd.Series(pred).value_counts().to_string())


if __name__ == "__main__":
    main()
