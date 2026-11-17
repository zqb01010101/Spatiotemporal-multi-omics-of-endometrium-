#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""10_compare_before_after.py — Compare model performance before and after label review on shared grouped CV splits

Reuse the audit cv_splits.csv (repeat 1 by default, see config.compare_repeats).
Train the Original MLP separately on original and revised labels. Compare:
  Four-score regression MSE, RMSE, MAE, and correlation, plus hard-label accuracy, balanced accuracy, and macro-F1
Outputs: outputs/before_after_comparison.csv and figures/before_after.png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             f1_score, mean_squared_error, mean_absolute_error)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("10_compare_before_after")
train_mod = common.import_train_module()


def evaluate_split(audit, cv, repeat, label_suffix=""):
    X_all = audit[common.FEATURES].to_numpy(dtype=float)
    score_cols = common.HUMAN_SCORE_COLS
    Y_all = audit[score_cols].to_numpy(dtype=float) / 10.0
    hard_all = audit["hard_label"].astype(str).to_numpy()
    groups = audit["group_id"].to_numpy()

    preds = np.zeros_like(Y_all)
    n_folds = common.cfg("n_folds", 5)
    seed = common.cfg("random_seed", 42)
    for f in range(n_folds):
        te_sid = set(cv.loc[(cv["repeat"] == repeat) & (cv["fold"] == f), "sample_id"])
        tr_idx = np.array([i for i, sid in enumerate(audit["sample_id"]) if sid not in te_sid])
        te_idx = np.array([i for i, sid in enumerate(audit["sample_id"]) if sid in te_sid])
        tr_g = groups[tr_idx]
        tr_s = hard_all[tr_idx]
        iv_tr, iv_te = common.inner_val_split(tr_g, tr_s, 0.1, seed + f)
        m = train_mod.ScoreMLP(12, [128, 64, 32], 4, seed=seed + f,
                               keep_prob=common.cfg("mlp_keep_prob", 0.8),
                               loss_weights=train_mod.OUTPUT_LOSS_WEIGHTS)
        train_mod.train(m, X_all[tr_idx][iv_tr], Y_all[tr_idx][iv_tr],
                        X_all[tr_idx][iv_te], Y_all[tr_idx][iv_te])
        preds[te_idx] = m.forward(X_all[te_idx]) * 10.0

    y = Y_all * 10.0
    rows = {}
    for i, c in enumerate(common.CLASSES):
        rows[f"{c}_mse"] = float(mean_squared_error(y[:, i], preds[:, i]))
        rows[f"{c}_rmse"] = float(np.sqrt(rows[f"{c}_mse"]))
        rows[f"{c}_mae"] = float(mean_absolute_error(y[:, i], preds[:, i]))
        rows[f"{c}_pearson"] = float(np.corrcoef(y[:, i], preds[:, i])[0, 1])
        rows[f"{c}_spearman"] = float(pd.Series(y[:, i]).corr(pd.Series(preds[:, i]), method="spearman"))
    pred_label, _, _, _ = train_mod.derive_hard_label(preds)
    pred_label_str = np.array([common.CLASSES[i] if i >= 0 else "uncertain" for i in pred_label])
    known = hard_all != "uncertain"
    rows["hard_accuracy"] = float(accuracy_score(hard_all[known], pred_label_str[known]))
    rows["hard_balanced"] = float(balanced_accuracy_score(hard_all[known], pred_label_str[known]))
    rows["hard_macro_f1"] = float(f1_score(
        hard_all[known], pred_label_str[known], average="macro",
        labels=list(common.CLASSES), zero_division=0))
    return rows, pred_label_str, hard_all


def main():
    clean_path = common.OUT_DIR / "labels_clean_v1.csv"
    if not clean_path.exists():
        logger.info("尚未生成 labels_clean_v1.csv（需先完成人工复核并运行 09），跳过对比。")
        return

    audit_orig = common.load_audit_dataset()
    audit_clean = pd.read_csv(clean_path)
    cv = pd.read_csv(common.DATA_DIR / "cv_splits.csv")
    repeat = common.cfg("compare_repeats", 1)

    logger.info("开始 before/after 对比（repeat=%d，5-fold grouped CV）", repeat)
    before, _, _ = evaluate_split(audit_orig, cv, repeat)
    after, _, _ = evaluate_split(audit_clean, cv, repeat)

    comp = pd.DataFrame({"metric": list(before.keys()),
                         "before": list(before.values()),
                         "after": list(after.values())})
    comp["delta"] = comp["after"] - comp["before"]
    comp.to_csv(common.OUT_DIR / "before_after_comparison.csv", index=False)
    logger.info("已写出 before_after_comparison.csv")
    print(comp.to_string(index=False))

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(comp))
    ax.bar(x - 0.2, comp["before"], 0.4, label="before")
    ax.bar(x + 0.2, comp["after"], 0.4, label="after")
    ax.set_xticks(x, comp["metric"], rotation=60, ha="right")
    ax.legend()
    ax.set_title("Before vs After label cleaning")
    fig.tight_layout()
    fig.savefig(common.FIG_DIR / "before_after.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
