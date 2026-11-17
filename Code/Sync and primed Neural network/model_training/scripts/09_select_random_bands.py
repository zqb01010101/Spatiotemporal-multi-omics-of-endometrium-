#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""09_select_random_bands.py — Sample review candidates within predicted score bands

By default, randomly sample 50 unlabeled candidates from each of three bands:
  1) Pre-opening: score_pre_open in [7.0, 8.5] on the 0–10 scale
  2) Synchronous opening: score_sync_open in [7.5, 8.5]
  3) Synchronous repression: score_sync_repress in [8.5, 9.5]

Override score bands with command-line options:
  python 09_select_random_bands.py \
      --pre 7.0 8.5 --open 7.5 8.5 --repress 8.5 9.5 \
      --prefix review_v3 --out review_batch_v3.csv --seed 42

Output: data/review_batch_*.csv for import with 06_import_review_batch.R
"""

import argparse
import pandas as pd
from pathlib import Path

MODEL_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = MODEL_DIR / "data"
OUT_DIR = MODEL_DIR / "outputs"
APP_LABELS = MODEL_DIR.parent / "peak_gene_label_app" / "labels"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre", nargs=2, type=float, default=[7.0, 8.5])
    ap.add_argument("--open", nargs=2, type=float, default=[7.5, 8.5])
    ap.add_argument("--repress", nargs=2, type=float, default=[8.5, 9.5])
    ap.add_argument("--prefix", default="review_v3")
    ap.add_argument("--out", default="review_batch_v3.csv")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n", type=int, default=50)
    args = ap.parse_args()
    n_per_band = args.n

    bands = [
        (f"{args.prefix}_pre_open_{args.pre[0]:.2f}_{args.pre[1]:.2f}", "score_pre_open", *args.pre),
        (f"{args.prefix}_sync_open_{args.open[0]:.2f}_{args.open[1]:.2f}", "score_sync_open", *args.open),
        (f"{args.prefix}_sync_repress_{args.repress[0]:.2f}_{args.repress[1]:.2f}", "score_sync_repress", *args.repress),
    ]

    scores = pd.read_csv(OUT_DIR / "all_candidate_scores_with_context.csv")
    ann = pd.read_csv(APP_LABELS / "annotations_export.csv")
    annotated = set(ann["candidate_id"].astype(str))
    unannotated = ~scores["candidate_id"].astype(str).isin(annotated)
    print("全部候选：", len(scores), "；已标注：", len(annotated))

    parts = []
    for bid, col, lo, hi in bands:
        sub = scores[unannotated & (scores[col] >= lo) & (scores[col] <= hi)]
        print(f"{bid}: 可用 {len(sub)} 条（{col} ∈ [{lo}, {hi}]）")
        if len(sub) < n_per_band:
            print(f"  !! 不足 {n_per_band} 条，按可用数量抽取")
        sel = sub.sample(n=min(n_per_band, len(sub)), random_state=args.seed)
        sel = sel[["candidate_id", "context", "peak", "gene"]].copy()
        sel["batch_id"] = bid
        sel["batch_score"] = sel["candidate_id"].map(
            sub.set_index("candidate_id")[col])
        sel["note"] = f"{col} in [{lo}, {hi}]"
        parts.append(sel)

    batch = pd.concat(parts, ignore_index=True)
    batch = batch[["candidate_id", "context", "peak", "gene",
                   "batch_id", "batch_score", "note"]]
    print("重复 candidate_id 数：", int(batch["candidate_id"].duplicated().sum()))
    out = DATA_DIR / args.out
    batch.to_csv(out, index=False)
    print("已写出：", out, "（", len(batch), " 条）")
    print(batch["batch_id"].value_counts().to_string())


if __name__ == "__main__":
    main()
