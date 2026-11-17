#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""03_run_oof_models.py — Generate repeated grouped out-of-fold predictions for three models

Model A: Original MLP (v4.1 with 12 features, no z-scoring, and dropout)
Model B: Context-aware MLP (same architecture with 12 features plus binary context)
Model C: ExtraTreesRegressor (12 features, multi-output regression)

Each OOF prediction comes from a model trained without that sample's group.
Outputs:
  outputs/oof_predictions_long.csv      sample_id, repeat, fold, model, pred_*
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("03_run_oof_models")
train_mod = common.import_train_module()


def train_mlp_oof(X_tr, Y_tr, groups_tr, strat_tr, seed, n_in):
    tr_idx, iv_idx = common.inner_val_split(
        groups_tr, strat_tr, common.cfg("inner_val_fraction", 0.1), seed)
    model = train_mod.ScoreMLP(
        n_in, [128, 64, 32], 4, seed=seed,
        keep_prob=common.cfg("mlp_keep_prob", 0.8),
        loss_weights=train_mod.OUTPUT_LOSS_WEIGHTS)
    train_mod.train(model, X_tr[tr_idx], Y_tr[tr_idx], X_tr[iv_idx], Y_tr[iv_idx])
    return model


def main():
    audit = common.load_audit_dataset()
    cv = pd.read_csv(common.DATA_DIR / "cv_splits.csv")
    groups = audit.set_index("sample_id")["group_id"].to_dict()
    strat = audit.set_index("sample_id")["hard_label"].astype(str).to_dict()

    X_all = audit[common.FEATURES].to_numpy(dtype=float)
    Y_all = audit[common.HUMAN_SCORE_COLS].to_numpy(dtype=float) / 10.0
    ctx_bin = (audit["context"] == "epithelial").astype(float).to_numpy()
    X_ctx_all = np.column_stack([X_all, ctx_bin])

    n_folds = common.cfg("n_folds", 5)
    n_repeats = common.cfg("n_repeats", 5)
    seed = common.cfg("random_seed", 42)

    rows = []
    t0 = time.time()
    for r in range(1, n_repeats + 1):
        for f in range(n_folds):
            te_sid = set(cv.loc[(cv["repeat"] == r) & (cv["fold"] == f), "sample_id"])
            tr_idx = np.array([i for i, sid in enumerate(audit["sample_id"]) if sid not in te_sid])
            te_idx = np.array([i for i, sid in enumerate(audit["sample_id"]) if sid in te_sid])
            groups_tr = np.array([groups[s] for s in audit["sample_id"].iloc[tr_idx]])
            strat_tr = np.array([strat[s] for s in audit["sample_id"].iloc[tr_idx]])
            seed_rf = seed + 1000 * r + 10 * f

            # Model A: Original MLP
            mA = train_mlp_oof(X_all[tr_idx], Y_all[tr_idx], groups_tr, strat_tr, seed_rf, 12)
            predA = mA.forward(X_all[te_idx]) * 10.0
            # Model B: Context-aware MLP
            mB = train_mlp_oof(X_ctx_all[tr_idx], Y_all[tr_idx], groups_tr, strat_tr, seed_rf + 1, 13)
            predB = mB.forward(X_ctx_all[te_idx]) * 10.0
            # Model C: ExtraTrees
            from sklearn.ensemble import ExtraTreesRegressor
            mC = ExtraTreesRegressor(
                n_estimators=common.cfg("tree_n_estimators", 300),
                min_samples_leaf=2, n_jobs=-1,
                random_state=seed_rf + 2)
            mC.fit(X_all[tr_idx], Y_all[tr_idx])
            predC = mC.predict(X_all[te_idx]) * 10.0

            for i, sid in enumerate(audit["sample_id"].iloc[te_idx]):
                rows.append({
                    "sample_id": sid, "repeat": r, "fold": f, "model": "original_mlp",
                    "pred_pre_open": predA[i, 0], "pred_sync_open": predA[i, 1],
                    "pred_sync_repress": predA[i, 2], "pred_unrelated": predA[i, 3]})
                rows.append({
                    "sample_id": sid, "repeat": r, "fold": f, "model": "context_mlp",
                    "pred_pre_open": predB[i, 0], "pred_sync_open": predB[i, 1],
                    "pred_sync_repress": predB[i, 2], "pred_unrelated": predB[i, 3]})
                rows.append({
                    "sample_id": sid, "repeat": r, "fold": f, "model": "extratrees",
                    "pred_pre_open": predC[i, 0], "pred_sync_open": predC[i, 1],
                    "pred_sync_repress": predC[i, 2], "pred_unrelated": predC[i, 3]})
            logger.info("[r=%d/%d f=%d/%d] done, elapsed=%.0fs",
                        r, n_repeats, f + 1, n_folds, time.time() - t0)

    long_df = pd.DataFrame(rows)
    long_df.to_csv(common.OUT_DIR / "oof_predictions_long.csv", index=False)
    logger.info("已写出 oof_predictions_long.csv（%d 行，耗时 %.0fs）", len(long_df), time.time() - t0)


if __name__ == "__main__":
    main()
