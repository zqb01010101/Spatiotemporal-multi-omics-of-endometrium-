#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""07_rank_label_issues.py — Combine evidence, compute audit_score, and assign tiers A–D

audit_score ranks samples within a tier. See the README and metadata for the formula.
  Average the 0–100 percentile ranks of the component signals:
  vector_mae_rank, max_abs_error_rank, contradiction_gap_rank,
  neighbor_score_mae_rank, (1 - model_prediction_dispersion)_rank,
  (1 - prediction_sd_mean)_rank, model_consensus_fraction_rank,
  repeat_score_mae_rank (contributes zero if there is no repeat)

Tier rules (see README):
  Tier A —— strong review candidate
  Tier B —— probable review candidate
  Tier C — ambiguous or difficult labels, without a confirmed labeling error
  Tier D —— currently consistent

Outputs:
  outputs/label_audit_ranked.csv
  outputs/batch_audit_report.csv
  outputs/audit_meta.json
  figures/ (human vs. OOF scatter plots, residuals, residuals vs. uncertainty,
           batch residuals, tier proportions by batch, model agreement, and audit scores)
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("07_rank_label_issues")


def pct_rank(s):
    return s.rank(pct=True) * 100.0


def make_tier(df):
    ma = df["vector_mae_pct"]
    conflict = df["hard_label_conflict"]
    consensus = df["model_consensus_fraction"]
    sd_low = df["prediction_sd_mean_pct"] <= 50
    disp_high = df["model_prediction_dispersion_pct"] >= 60
    gap_big = df["contradiction_gap"] >= 5.0
    rep_conflict = df.get("repeat_hard_label_conflict", pd.Series(0, index=df.index)).fillna(0).astype(int) == 1

    tier = pd.Series("D", index=df.index)
    reason = [""] * len(df)
    for i in df.index:
        parts = []
        if conflict.loc[i] and ma.loc[i] >= 70:
            parts.append("hard-label conflict with OOF")
        if gap_big.loc[i]:
            parts.append(f"large contradiction gap ({df.loc[i, 'contradiction_gap']:.1f})")
        if consensus.loc[i] >= 2 / 3:
            parts.append("3/3 or 2/3 model consensus")
        if sd_low.loc[i]:
            parts.append("stable repeated predictions")
        if ma.loc[i] >= 80:
            parts.append("large OOF residual")
        if rep_conflict.loc[i]:
            parts.append("repeat annotation conflict")
        if disp_high.loc[i]:
            parts.append("models strongly disagree")
        if ma.loc[i] >= 80 and sd_low.loc[i]:
            parts.append("large residual but low variance")
        reason[i] = "; ".join(parts)

        strong = (conflict.loc[i] and consensus.loc[i] >= 2 / 3 and
                  (sd_low.loc[i] or gap_big.loc[i]))
        if strong or rep_conflict.loc[i]:
            tier.loc[i] = "A"
        elif conflict.loc[i] and ma.loc[i] >= 70:
            tier.loc[i] = "B"
        elif ma.loc[i] >= 80 and (not sd_low.loc[i] or disp_high.loc[i]):
            tier.loc[i] = "C"
    return tier, pd.Series(reason, index=df.index)


def main():
    audit = common.load_audit_dataset()
    metrics = pd.read_csv(common.OUT_DIR / "audit_metrics.csv")
    neigh = pd.read_csv(common.OUT_DIR / "neighbor_consistency.csv")
    rep = pd.read_csv(common.OUT_DIR / "repeat_consistency.csv")

    df = metrics.merge(
        neigh[["sample_id", "neighbor_score_mae", "neighbor_top_class",
               "neighbor_top_class_fraction"]], on="sample_id", how="left")
    rep_map = None
    if len(rep):
        rep_map = rep.groupby("original_sample_id").agg(
            repeat_score_mae=("repeat_score_mae", "mean"),
            repeat_hard_label_conflict=("repeat_hard_label_conflict", "max")).reset_index()
        df = df.merge(rep_map, left_on="sample_id", right_on="original_sample_id", how="left")
    for c in ["neighbor_score_mae", "repeat_score_mae", "repeat_hard_label_conflict"]:
        if c not in df.columns:
            df[c] = 0
        else:
            df[c] = df[c].fillna(0)

    # ---------- Convert scores to percentile ranks ----------
    df["vector_mae_pct"] = pct_rank(df["vector_mae"])
    df["max_abs_error_pct"] = pct_rank(df["max_abs_error"])
    df["contradiction_gap_pct"] = pct_rank(df["contradiction_gap"].clip(lower=0))
    df["neighbor_score_mae_pct"] = pct_rank(df["neighbor_score_mae"])
    df["model_prediction_dispersion_pct"] = pct_rank(df["model_prediction_dispersion"])
    df["prediction_sd_mean_pct"] = pct_rank(df["prediction_sd_mean"])
    df["repeat_score_mae_pct"] = pct_rank(df["repeat_score_mae"])
    df["stability_pct"] = 100 - df["prediction_sd_mean_pct"]
    df["consensus_pct"] = pct_rank(df["model_consensus_fraction"])

    score_parts = [
        df["vector_mae_pct"], df["max_abs_error_pct"],
        df["contradiction_gap_pct"], df["neighbor_score_mae_pct"],
        100 - df["model_prediction_dispersion_pct"],
        df["stability_pct"], df["consensus_pct"], df["repeat_score_mae_pct"],
    ]
    df["audit_score"] = np.mean(np.column_stack(score_parts), axis=1)

    tier, reason = make_tier(df)
    df["audit_tier"] = tier
    df["audit_reason"] = reason
    df = df.sort_values(["audit_tier", "audit_score"], ascending=[True, False])
    df.to_csv(common.OUT_DIR / "label_audit_ranked.csv", index=False)

    logger.info("Tier 分布:\n%s", df["audit_tier"].value_counts().to_string())
    logger.info("Tier A 主要可疑原因（top5）：\n%s",
                df[df["audit_tier"] == "A"]["audit_reason"].value_counts().head(5).to_string())

    # ---------- Summarize by batch ----------
    batch = df.groupby("batch_id").agg(
        n=("sample_id", "count"),
        median_vector_mae=("vector_mae", "median"),
        strong_contradiction_rate=("hard_label_conflict", "mean"),
        tier_A=("audit_tier", lambda s: (s == "A").sum()),
        tier_B=("audit_tier", lambda s: (s == "B").sum()),
        tier_C=("audit_tier", lambda s: (s == "C").sum()),
        tier_D=("audit_tier", lambda s: (s == "D").sum()),
    ).reset_index()
    batch.to_csv(common.OUT_DIR / "batch_audit_report.csv", index=False)

    # ---------- figures ----------
    for i, c in enumerate(common.CLASSES):
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(df[f"human_{c}"], df[f"oof_{c}"], s=8, alpha=0.4)
        lim = [0, 10]
        ax.plot(lim, lim, "k--", lw=0.8)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("human"); ax.set_ylabel("OOF (original MLP)")
        ax.set_title(f"Human vs OOF: {c}")
        fig.tight_layout()
        fig.savefig(common.FIG_DIR / f"human_vs_oof_{c}.png", dpi=130)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(df["vector_mae"], bins=40)
    ax.set_xlabel("vector MAE (0-10)"); ax.set_ylabel("count")
    ax.set_title("OOF residual distribution")
    fig.tight_layout(); fig.savefig(common.FIG_DIR / "residual_distribution.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(df["prediction_sd_mean"], df["vector_mae"], s=8, alpha=0.4)
    ax.set_xlabel("OOF prediction sd (mean over 4)"); ax.set_ylabel("vector MAE")
    ax.set_title("Residual vs prediction uncertainty")
    fig.tight_layout(); fig.savefig(common.FIG_DIR / "residual_vs_uncertainty.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4))
    order = batch.sort_values("median_vector_mae")["batch_id"]
    ax.barh(order, batch.set_index("batch_id").loc[order, "median_vector_mae"])
    ax.set_xlabel("median vector MAE"); ax.set_title("Per-batch OOF residual")
    fig.tight_layout(); fig.savefig(common.FIG_DIR / "batch_residual.png", dpi=130)
    plt.close(fig)

    tier_plot = df.replace("", "none").groupby(["batch_id", "audit_tier"]).size().unstack(fill_value=0)
    tier_plot = tier_plot.div(tier_plot.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(9, 4))
    tier_plot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("proportion"); ax.set_title("Per-batch Tier proportion")
    fig.tight_layout(); fig.savefig(common.FIG_DIR / "batch_tier.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(df["model_consensus_fraction"], bins=20)
    ax.set_xlabel("model consensus fraction"); ax.set_title("Model-model agreement")
    fig.tight_layout(); fig.savefig(common.FIG_DIR / "model_agreement.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(df["audit_score"], bins=40)
    ax.set_xlabel("audit score"); ax.set_title("Audit score distribution")
    fig.tight_layout(); fig.savefig(common.FIG_DIR / "audit_score_distribution.png", dpi=130)
    plt.close(fig)

    meta = {
        "audit_score_formula": (
            "mean of percentile ranks: vector_mae, max_abs_error, contradiction_gap, "
            "neighbor_score_mae, 100-dispersion, 100-pred_sd_mean, consensus_fraction, "
            "repeat_score_mae (0 if none)"),
        "tier_rules": (
            "A: (conflict & consensus>=2/3 & (stable|gap>=5)) or repeat conflict; "
            "B: conflict & vector_mae_pct>=70; "
            "C: vector_mae_pct>=80 & (unstable|models disagree); D: rest"),
        "n_tier_a": int((df["audit_tier"] == "A").sum()),
        "n_tier_b": int((df["audit_tier"] == "B").sum()),
        "n_tier_c": int((df["audit_tier"] == "C").sum()),
        "n_tier_d": int((df["audit_tier"] == "D").sum()),
    }
    with open(common.OUT_DIR / "audit_meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)
    logger.info("已写出 label_audit_ranked.csv / batch_audit_report.csv / audit_meta.json 与 figures")


if __name__ == "__main__":
    main()
