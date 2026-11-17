#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""16_organoid_score.py — Score organoid features with the existing v5 model

Output directory: model_training/organoid_run/
  all_candidate_scores.csv / all_candidate_scores_with_context.csv
  queue_predictions.csv
  review_pdfs_organoid/ (four PDFs with the _organoid suffix)
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

MODEL_DIR = Path(__file__).resolve().parents[1]
ORUN = MODEL_DIR / "organoid_run"
MODEL_NPZ = MODEL_DIR / "outputs" / "model" / "model.npz"
PDF_DIR = ORUN / "review_pdfs_organoid"
PDF_DIR.mkdir(parents=True, exist_ok=True)

STAGES = ["PO3", "PO5", "PO7", "PO9", "PO11", "PO13"]
FEATURES = [f"peak_{s}" for s in STAGES] + [f"gene_{s}" for s in STAGES]
CLASSES = ["pre_open", "sync_open", "sync_repress", "unrelated"]
SCORE_COLS = [f"score_{c}" for c in CLASSES]
TOP_N, LOWEST_N, PER_PAGE = 20, 50, 6


def relu(x):
    return np.maximum(0.0, x)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def predict(X):
    w = np.load(MODEL_NPZ)
    keys = sorted([k for k in w.files if k.startswith("W")], key=lambda k: int(k[1:]))
    z = X
    for i, wk in enumerate(keys):
        z = z @ w[wk] + w["b" + wk[1:]]
        z = relu(z) if i < len(keys) - 1 else sigmoid(z)
    return z * 10.0


def render_candidate(ax, row):
    x = np.arange(len(STAGES))
    ax.plot(x, row[FEATURES[:6]].to_numpy(dtype=float), "o-", color="#c0392b", lw=1.3, ms=3, label="peak")
    ax.plot(x, row[FEATURES[6:]].to_numpy(dtype=float), "s-", color="#2471a3", lw=1.3, ms=3, label="gene")
    ax.set_xticks(x, STAGES, fontsize=6)
    ax.set_ylim(0, 1)
    t = f"{row['peak']} | {row['gene']}\n"
    t += ("pre={:.1f} open={:.1f} rep={:.1f} unrel={:.1f} conf={:.2f}"
          .format(row["score_pre_open"], row["score_sync_open"],
                  row["score_sync_repress"], row["score_unrelated"], row["confidence"]))
    if pd.notna(row.get("ann_scores")):
        t += f"\nann={row['ann_scores']}"
    ax.set_title(t, fontsize=6.2, loc="left")
    if ax.get_subplotspec().is_first_col() and ax.get_subplotspec().is_last_row():
        ax.legend(fontsize=6, loc="upper left")


def make_pdf(rows, name, title):
    if not len(rows):
        return
    pages = int(np.ceil(len(rows) / PER_PAGE))
    with PdfPages(PDF_DIR / name) as pdf:
        for pg in range(pages):
            fig = plt.figure(figsize=(11.69, 8.27))
            fig.suptitle(f"{title} (page {pg + 1}/{pages})", fontsize=11)
            for j in range(PER_PAGE):
                idx = pg * PER_PAGE + j
                if idx >= len(rows):
                    break
                ax = fig.add_subplot(2, 3, j + 1)
                render_candidate(ax, rows.iloc[idx])
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)
    print(f"已生成：{PDF_DIR / name}（{len(rows)} 条，{pages} 页）")


def main():
    feat = pd.read_csv(ORUN / "all_candidates_features.csv")
    X = feat[FEATURES].to_numpy(dtype=float)
    out = predict(X)
    for i, c in enumerate(SCORE_COLS):
        feat[c] = out[:, i]

    table = feat[["peak", "gene"] + SCORE_COLS].sort_values(["peak", "gene"])
    table.to_csv(ORUN / "all_candidate_scores.csv", index=False)
    full = feat[["candidate_id", "context", "peak", "gene"] + SCORE_COLS].sort_values(
        ["context", "peak", "gene"])
    full.to_csv(ORUN / "all_candidate_scores_with_context.csv", index=False)
    print("已写出 all_candidate_scores.csv / with_context（", len(feat), " 行）")

    queue = pd.read_csv(ORUN / "queue_features.csv")
    qpred = full[full["candidate_id"].isin(queue["candidate_id"])].copy()
    top = np.array(CLASSES)[out[feat["candidate_id"].isin(queue["candidate_id"]).to_numpy()].argmax(axis=1)]
    qscores = qpred[SCORE_COLS].to_numpy(dtype=float)
    order = np.argsort(-qscores, axis=1)
    qpred["pred_class"] = [CLASSES[i] for i in order[:, 0]]
    qpred["confidence"] = qscores[np.arange(len(qscores)), order[:, 0]] / 10.0
    qpred["margin"] = (qscores[np.arange(len(qscores)), order[:, 0]] -
                       qscores[np.arange(len(qscores)), order[:, 1]]) / 10.0
    qpred = qpred.sort_values("confidence")
    qpred.to_csv(ORUN / "queue_predictions.csv", index=False)
    print("已写出 queue_predictions.csv（", len(qpred), " 行）")

    ann = pd.read_csv(MODEL_DIR / "data" / "annotated_samples.csv")
    ann_map = {
        r.candidate_id: "{:.0f}/{:.0f}/{:.0f}/{:.0f}".format(
            r.score_pre_open, r.score_sync_open, r.score_sync_repress, r.score_unrelated)
        for r in ann.itertuples()
    }
    qpred["ann_scores"] = qpred["candidate_id"].map(ann_map)
    qpred = qpred.merge(queue[["candidate_id"] + FEATURES], on="candidate_id", how="left")

    for cls in ["pre_open", "sync_open", "sync_repress"]:
        sel = qpred.sort_values(f"score_{cls}", ascending=False).head(TOP_N)
        make_pdf(sel, f"top20_{cls}_organoid.pdf", f"Top-{TOP_N} {cls} (organoid)")
    low = qpred.head(LOWEST_N)
    make_pdf(low, "lowest50_confidence_organoid.pdf", f"Lowest-{LOWEST_N} confidence (organoid)")


if __name__ == "__main__":
    main()
