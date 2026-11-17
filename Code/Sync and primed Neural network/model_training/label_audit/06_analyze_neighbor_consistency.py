#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""06_analyze_neighbor_consistency.py — Check nearest-neighbor agreement independently of the models

Find k nearest neighbors within each context using Euclidean distance on the 12 unscaled features.
Compare the four human scores with neighborhood medians:
  neighbor_score_mae / neighbor_top_class / neighbor_top_class_fraction
Output: outputs/neighbor_consistency.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("06_analyze_neighbor_consistency")


def main():
    audit = common.load_audit_dataset()
    k = common.cfg("neighbor_k", 15)
    rows = []
    for ctx in ["stroma", "epithelial"]:
        sub = audit[audit["context"] == ctx].reset_index(drop=True)
        X = sub[common.FEATURES].to_numpy(dtype=float)
        D = cdist(X, X)
        np.fill_diagonal(D, np.inf)
        for i in range(len(sub)):
            nbr = np.argpartition(D[i], min(k, len(sub) - 1))[:k]
            nbr_scores = sub.iloc[nbr][common.HUMAN_SCORE_COLS].to_numpy(dtype=float)
            med = np.median(nbr_scores, axis=0)
            target = sub.iloc[i][common.HUMAN_SCORE_COLS].to_numpy(dtype=float)
            target_class = common.CLASSES[int(target.argmax())]
            nbr_class = np.array(common.CLASSES)[nbr_scores.argmax(axis=1)]
            rows.append({
                "sample_id": sub.iloc[i]["sample_id"],
                "neighbor_k": len(nbr),
                "neighbor_score_mae": float(np.abs(target - med).mean()),
                "neighbor_top_class": common.CLASSES[int(med.argmax())],
                "neighbor_top_class_fraction": float(np.mean(nbr_class == target_class)),
            })
    out = pd.DataFrame(rows)
    out.to_csv(common.OUT_DIR / "neighbor_consistency.csv", index=False)
    logger.info("已写出 neighbor_consistency.csv（%d 行，k=%d）", len(out), k)


if __name__ == "__main__":
    main()
