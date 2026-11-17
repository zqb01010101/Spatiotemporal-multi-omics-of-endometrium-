#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""03_generate_review_pdfs.py — Create review PDFs for the v3 regression model

Outputs in outputs/review_pdfs/:
  top20_pre_open.pdf      Top 20 predicted pre-opening scores
  top20_sync_open.pdf     Top 20 predicted synchronous opening scores
  top20_sync_repress.pdf  Top 20 predicted synchronous repression scores
  lowest50_confidence.pdf 50 lowest-confidence predictions
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

MODEL_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = MODEL_DIR / "data"
REVIEW_DIR = MODEL_DIR / "outputs" / "review_pdfs"
REVIEW_DIR.mkdir(parents=True, exist_ok=True)

STAGES = ["PO3", "PO5", "PO7", "PO9", "PO11", "PO13"]
CLASSES = ["pre_open", "sync_open", "sync_repress", "unrelated"]
SCORE_COLS = [f"score_{c}" for c in CLASSES]
PEAK_COLS = [f"peak_{s}" for s in STAGES]
GENE_COLS = [f"gene_{s}" for s in STAGES]
TOP_N = 20
LOWEST_N = 50
PER_PAGE = 6


def render_candidate(ax, row):
    x = np.arange(len(STAGES))
    peak = row[PEAK_COLS].to_numpy(dtype=float)
    gene = row[GENE_COLS].to_numpy(dtype=float)
    ax.plot(x, peak, "o-", color="#c0392b", lw=1.4, ms=3, label="peak")
    ax.plot(x, gene, "s-", color="#2471a3", lw=1.4, ms=3, label="gene")
    ax.set_xticks(x, STAGES, fontsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax.set_ylim(0, 1)
    title = f"{row['peak']} | {row['gene']}\n"
    title += ("pre={:.1f} open={:.1f} rep={:.1f} unrel={:.1f}  conf={:.2f}"
              .format(row["score_pre_open"], row["score_sync_open"],
                      row["score_sync_repress"], row["score_unrelated"],
                      row["confidence"]))
    if pd.notna(row.get("ann_scores")):
        title += f"\nann={row['ann_scores']}"
    ax.set_title(title, fontsize=6.2, loc="left")
    if ax.get_subplotspec().is_first_col() and ax.get_subplotspec().is_last_row():
        ax.legend(fontsize=6, loc="upper left")


def make_pdf(rows, out_name, title):
    n = len(rows)
    if n == 0:
        print(f"跳过 {out_name}：无候选")
        return
    pages = int(np.ceil(n / PER_PAGE))
    with PdfPages(REVIEW_DIR / out_name) as pdf:
        for pg in range(pages):
            fig = plt.figure(figsize=(11.69, 8.27))
            fig.suptitle(f"{title}  (page {pg + 1}/{pages})", fontsize=11)
            for i in range(PER_PAGE):
                idx = pg * PER_PAGE + i
                if idx >= n:
                    break
                ax = fig.add_subplot(2, 3, i + 1)
                render_candidate(ax, rows.iloc[idx])
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)
    print(f"已生成：{REVIEW_DIR / out_name}（{n} 条，{pages} 页）")


def main():
    pred = pd.read_csv(DATA_DIR / "queue_predictions.csv")
    feat = pd.read_csv(DATA_DIR / "queue_features_clean.csv")
    df = pred.merge(feat, on=["candidate_id", "context", "peak", "gene", "distance"],
                    how="left", validate="one_to_one")

    ann = pd.read_csv(DATA_DIR / "annotated_samples.csv")
    ann_map = {
        r.candidate_id: "{:.0f}/{:.0f}/{:.0f}/{:.0f}".format(
            r.score_pre_open, r.score_sync_open, r.score_sync_repress, r.score_unrelated)
        for r in ann.itertuples()
    }
    df["ann_scores"] = df["candidate_id"].map(ann_map)

    selections = {}
    for cls in ["pre_open", "sync_open", "sync_repress"]:
        sel = df.sort_values(f"score_{cls}", ascending=False).head(TOP_N).reset_index(drop=True)
        selections[f"top20_{cls}"] = sel
        make_pdf(sel, f"top20_{cls}.pdf", f"Top-{TOP_N} predicted score '{cls}'")
        sel.to_csv(DATA_DIR / f"review_top20_{cls}.csv", index=False)

    low = df.sort_values("confidence", ascending=True).head(LOWEST_N).reset_index(drop=True)
    selections["lowest50_confidence"] = low
    make_pdf(low, "lowest50_confidence.pdf", f"Lowest-{LOWEST_N} model confidence")
    low.to_csv(DATA_DIR / "review_lowest50_confidence.csv", index=False)

    for name, sel in selections.items():
        n_ann = int(sel["ann_scores"].notna().sum())
        print(f"  {name}: {len(sel)} 条，其中已标注 {n_ann} 条")


if __name__ == "__main__":
    main()
