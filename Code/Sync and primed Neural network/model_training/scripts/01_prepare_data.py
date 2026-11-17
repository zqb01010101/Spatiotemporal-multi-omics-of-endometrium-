#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""01_prepare_data.py — Build regression training, validation, and prediction sets

Fit four independent 0–10 scores, scaled to [0, 1] by dividing by 10.
Use hard_label only for stratified splitting and evaluation.

Inputs:
  peak_gene_label_app/labels/labeled_training_data.csv   (labeled samples)
  data/queue_features.csv                                (labeling queue features)
Outputs:
  data/annotated_samples.csv   Labeled sample manifest
  data/train_test.npz          X_train/X_val (12 features), Y_train/Y_val (four normalized scores),
                               hard_train/hard_val (hard_label indices with -1 for uncertain), and X_queue
  data/class_meta.json         Label and split statistics
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

MODEL_DIR = Path(__file__).resolve().parents[1]
APP_DIR = MODEL_DIR.parent / "peak_gene_label_app"
DATA_DIR = MODEL_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

STAGES = ["PO3", "PO5", "PO7", "PO9", "PO11", "PO13"]
FEATURES = [f"peak_{s}" for s in STAGES] + [f"gene_{s}" for s in STAGES]
CLASSES = ["pre_open", "sync_open", "sync_repress", "unrelated"]
SCORE_COLS = ["score_pre_open", "score_sync_open", "score_sync_repress", "score_unrelated"]

RANDOM_STATE = 42


def main():
    labeled = pd.read_csv(APP_DIR / "labels" / "labeled_training_data.csv")
    print("已标注样本：", len(labeled))

    X = labeled[FEATURES].to_numpy(dtype=float)
    scores = labeled[SCORE_COLS].to_numpy(dtype=float)
    Y = scores / 10.0                     # Scale scores from 0–10 to [0, 1]

    # ---------- Export the labeled sample manifest ----------
    argmax_idx = scores.argmax(axis=1)
    second_max = np.partition(scores, -2, axis=1)[:, -2]
    ties = second_max == scores.max(axis=1)

    out = labeled[["candidate_id", "context", "peak", "gene", "distance",
                   "mandatory", "source_tri_pre", "source_sync_open",
                   "source_sync_repress_heuristic", "S_pre_shape", "S_final",
                   "overall_confidence", "notes", "updated_at"]].copy()
    for c, name in zip(SCORE_COLS, ["pre_open", "sync_open", "sync_repress", "unrelated"]):
        out[f"score_{name}"] = labeled[c].to_numpy()
    out["argmax_label"] = [CLASSES[i] for i in argmax_idx]
    out["score_tie"] = ties
    out["hard_label_rule"] = labeled["hard_label"].to_numpy()
    out = out.sort_values("candidate_id").reset_index(drop=True)
    out.to_csv(DATA_DIR / "annotated_samples.csv", index=False)
    print("已写出 data/annotated_samples.csv")

    # ---------- Encode hard_label for splitting and evaluation (-1 = uncertain) ----------
    label_str = labeled["hard_label"].to_numpy()
    hard_idx = np.array([CLASSES.index(c) if c in CLASSES else -1 for c in label_str])

    # ---------- Split by hard_label, with uncertain as its own stratum ----------
    strat = np.where(hard_idx >= 0, hard_idx, len(CLASSES))
    X_train, X_val, Y_train, Y_val, h_tr, h_va = train_test_split(
        X, Y, hard_idx, test_size=0.2, random_state=RANDOM_STATE, stratify=strat
    )
    print(f"train: {len(Y_train)}  val: {len(Y_val)}")
    print("val hard_label 分布:", {c: int((h_va == i).sum()) for i, c in enumerate(CLASSES)},
          "uncertain:", int((h_va < 0).sum()))

    # ---------- Prediction set: full labeling queue ----------
    queue_feat = pd.read_csv(DATA_DIR / "queue_features.csv")
    X_queue = queue_feat[FEATURES].to_numpy(dtype=float)
    queue_feat.to_csv(DATA_DIR / "queue_features_clean.csv", index=False)
    print("预测集（标注队列）候选数：", len(queue_feat))

    np.savez(
        DATA_DIR / "train_test.npz",
        X_train=X_train, Y_train=Y_train, hard_train=h_tr,
        X_val=X_val, Y_val=Y_val, hard_val=h_va,
        X_queue=X_queue,
    )
    meta = {
        "classes": CLASSES,
        "features": FEATURES,
        "stages": STAGES,
        "random_state": RANDOM_STATE,
        "n_labeled": int(len(labeled)),
        "n_train": int(len(Y_train)),
        "n_val": int(len(Y_val)),
        "n_queue": int(len(queue_feat)),
        "score_targets": ["score_pre_open", "score_sync_open", "score_sync_repress", "score_unrelated"],
        "val_hardlabel_counts": {c: int((h_va == i).sum()) for i, c in enumerate(CLASSES)},
        "val_uncertain": int((h_va < 0).sum()),
        "score_stats": {
            c: {"mean": float(scores[:, i].mean()), "sd": float(scores[:, i].std())}
            for i, c in enumerate(["pre_open", "sync_open", "sync_repress", "unrelated"])
        },
    }
    with open(DATA_DIR / "class_meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)
    print("完成：data/train_test.npz 与 data/class_meta.json")


if __name__ == "__main__":
    main()
