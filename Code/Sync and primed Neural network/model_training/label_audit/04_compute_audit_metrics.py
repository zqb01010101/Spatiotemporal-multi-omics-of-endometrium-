#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""04_compute_audit_metrics.py — Summarize OOF predictions and compute audit metrics

Outputs:
  outputs/oof_predictions_summary.csv  OOF mean and SD for each sample and model
  outputs/audit_metrics.csv            Per-sample metrics based on the Original MLP:
    vector_mae / vector_rmse / max_abs_error / residual for each class
    predicted_strong_label / hard_label_conflict / contradiction_gap
    prediction_sd_mean / prediction_sd_max
    model_consensus_n / model_consensus_fraction / model_prediction_dispersion
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("04_compute_audit_metrics")
train_mod = common.import_train_module()


def main():
    audit = common.load_audit_dataset()
    long_df = pd.read_csv(common.OUT_DIR / "oof_predictions_long.csv")

    pred_cols = [f"pred_{c}" for c in common.CLASSES]
    human = audit.set_index("sample_id")[common.HUMAN_SCORE_COLS]

    # ---------- OOF mean and SD per sample and model ----------
    summ = long_df.groupby(["sample_id", "model"])[pred_cols].agg(
        ["mean", "std"]).reset_index()
    summ.columns = ["sample_id", "model"] + [
        f"pred_{c}_{stat}" for c in common.CLASSES for stat in ("mean", "std")]
    summ.to_csv(common.OUT_DIR / "oof_predictions_summary.csv", index=False)
    logger.info("oof_predictions_summary.csv：%d 行", len(summ))

    # ---------- Derive strong labels by model ----------
    model_rows = []
    for model in ["original_mlp", "context_mlp", "extratrees"]:
        sub = summ[summ["model"] == model].copy()
        Xm = sub[[f"pred_{c}_mean" for c in common.CLASSES]].to_numpy(dtype=float)
        lab, _, _, _ = train_mod.derive_hard_label(Xm)
        sub["predicted_strong_label"] = [
            common.CLASSES[i] if i >= 0 else "uncertain" for i in lab]
        model_rows.append(sub)
    per_model = pd.concat(model_rows, ignore_index=True)
    per_model.to_csv(common.OUT_DIR / "oof_summary_by_model.csv", index=False)

    # ---------- Compute primary metrics using the Original MLP ----------
    orig = per_model[per_model["model"] == "original_mlp"].copy()
    ctx_mlp = per_model[per_model["model"] == "context_mlp"].copy()
    tree = per_model[per_model["model"] == "extratrees"].copy()

    out = audit[["sample_id", "context", "peak", "gene", "batch_id", "source",
                 "is_repeat", "repeat_of", "hard_label", "human_confidence"]].copy()
    out = out.merge(human.reset_index(), on="sample_id", how="left")
    for c in common.CLASSES:
        out[f"human_{c}"] = out[c]
    out = out.drop(columns=common.HUMAN_SCORE_COLS)

    po = orig.set_index("sample_id")
    for c in common.CLASSES:
        out[f"oof_{c}"] = po[f"pred_{c}_mean"].values
        out[f"oof_{c}_sd"] = po[f"pred_{c}_std"].values
    out["prediction_sd_mean"] = po[[f"pred_{c}_std" for c in common.CLASSES]].mean(axis=1).values
    out["prediction_sd_max"] = po[[f"pred_{c}_std" for c in common.CLASSES]].max(axis=1).values
    out["predicted_strong_label"] = po["predicted_strong_label"].values

    y = out[[f"human_{c}" for c in common.CLASSES]].to_numpy(dtype=float)
    yhat = out[[f"oof_{c}" for c in common.CLASSES]].to_numpy(dtype=float)
    err = y - yhat
    out["vector_mae"] = np.abs(err).mean(axis=1)
    out["vector_rmse"] = np.sqrt((err ** 2).mean(axis=1))
    out["max_abs_error"] = np.abs(err).max(axis=1)
    for i, c in enumerate(common.CLASSES):
        out[f"{c}_residual"] = err[:, i]

    out["hard_label_conflict"] = (
        out["predicted_strong_label"] != "uncertain") & (
        out["predicted_strong_label"] != out["hard_label"])
    # contradiction_gap: predicted top-class score minus the human score for that class
    top_idx = yhat.argmax(axis=1)
    gap = yhat[np.arange(len(yhat)), top_idx] - y[np.arange(len(y)), top_idx]
    out["contradiction_gap"] = gap

    # ---------- Compare agreement across models ----------
    def class_map(s):
        return common.CLASSES.index(s) if s != "uncertain" else -1
    for name, dfm in [("context_mlp", ctx_mlp), ("extratrees", tree)]:
        tmp = dfm.set_index("sample_id")["predicted_strong_label"]
        out[f"{name}_strong_label"] = out["sample_id"].map(tmp).values
    top3 = np.stack([
        out["predicted_strong_label"].map(class_map).to_numpy(),
        out["context_mlp_strong_label"].map(class_map).to_numpy(),
        out["extratrees_strong_label"].map(class_map).to_numpy(),
    ], axis=1)
    valid = top3 >= 0
    out["model_consensus_n"] = np.where(
        valid.sum(axis=1) > 0,
        np.array([np.bincount(top3[i][valid[i]], minlength=4).max()
                  for i in range(len(top3))]),
        0)
    out["model_consensus_fraction"] = (
        out["model_consensus_n"] / np.maximum(valid.sum(axis=1), 1)).round(4)
    # Prediction dispersion: mean across four classes of the score SD across three models
    stack = np.stack([
        po[[f"pred_{c}_mean" for c in common.CLASSES]].to_numpy(),
        ctx_mlp.set_index("sample_id")[[f"pred_{c}_mean" for c in common.CLASSES]].to_numpy(),
        tree.set_index("sample_id")[[f"pred_{c}_mean" for c in common.CLASSES]].to_numpy(),
    ])
    out["model_prediction_dispersion"] = stack.std(axis=0).mean(axis=1)

    out.to_csv(common.OUT_DIR / "audit_metrics.csv", index=False)
    logger.info("已写出 audit_metrics.csv（%d 行）", len(out))
    logger.info("hard_label_conflict 数: %d", int(out["hard_label_conflict"].sum()))


if __name__ == "__main__":
    main()
