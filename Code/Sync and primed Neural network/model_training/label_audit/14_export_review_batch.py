#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""14_export_review_batch.py — Export selected candidates for app import and back up prior scores

Usage:
  python label_audit/14_export_review_batch.py \
      --selection data/review_lowest50_confidence.csv \
      --batch-id review_v5_low_confidence \
      --score-col confidence \
      --out review_batch_v5_lowconf.csv

Outputs:
  model_training/data/<out>                          Input to 06_import_review_batch.R
  label_audit/outputs/<out prefix>_old_scores_backup.csv   Snapshot of prior scores
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

MODEL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True)
    ap.add_argument("--batch-id", required=True)
    ap.add_argument("--score-col", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sel = pd.read_csv(MODEL_DIR / args.selection)
    assert args.score_col in sel.columns, f"缺少分数列: {args.score_col}"

    batch = pd.DataFrame({
        "candidate_id": sel["candidate_id"],
        "context": sel["context"],
        "peak": sel["peak"],
        "gene": sel["gene"],
        "batch_id": args.batch_id,
        "batch_score": sel[args.score_col].round(4),
        "note": f"{args.batch_id} (score={args.score_col})",
    })
    assert not batch["candidate_id"].duplicated().any()
    dest = MODEL_DIR / "data" / args.out
    batch.to_csv(dest, index=False)
    print("已写出导入批次：", dest, "（", len(batch), " 条）")

    # Back up prior scores from the latest export while preserving the source data
    labeled = pd.read_csv(common.APP_LABELS / "labeled_training_data.csv").set_index("candidate_id")
    backup = batch[["candidate_id"]].rename(columns={"candidate_id": "sample_id"})
    for c in common.CLASSES:
        backup[f"old_{c}"] = backup["sample_id"].map(labeled[f"score_{c}"])
    backup["old_confidence"] = backup["sample_id"].map(labeled["overall_confidence"])
    backup["old_hard_label"] = backup["sample_id"].map(labeled["hard_label"])
    bpath = common.OUT_DIR / (Path(args.out).stem + "_old_scores_backup.csv")
    backup.to_csv(bpath, index=False)
    print("已写出旧分备份：", bpath)


if __name__ == "__main__":
    main()
