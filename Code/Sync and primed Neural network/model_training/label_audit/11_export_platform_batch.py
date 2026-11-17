#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""11_export_platform_batch.py — Export the audit review queue for import into the labeling app

Input: outputs/review_queue.csv (50 Tier A, 30 Tier B, and 20 random controls)
Output: model_training/data/review_batch_audit_v1.csv
      Columns: candidate_id, context, peak, gene, batch_id, batch_score, and note
Then run:
  Rscript scripts/06_import_review_batch.R review_batch_audit_v1.csv
This prioritizes the batch in the labeling app queue and displays batch_id in the status bar.
"""

import sys
from pathlib import Path

import pandas as pd

MODEL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="audit_v1")
    args = ap.parse_args()
    prefix = args.prefix
    group_to_batch = {
        "tier_a": f"{prefix}_tierA",
        "tier_b": f"{prefix}_tierB",
        "random_control": f"{prefix}_control",
    }

    queue = pd.read_csv(common.OUT_DIR / "review_queue.csv")
    audit = common.load_audit_dataset().set_index("sample_id")
    out = pd.DataFrame({
        "candidate_id": queue["sample_id"],
        "context": queue["context"],
        "peak": queue["peak"],
        "gene": queue["gene"],
        "batch_id": queue["review_group"].map(group_to_batch),
        "batch_score": queue["audit_score"].round(3),
        "note": ("audit tier=" + queue["audit_tier"].astype(str) +
                 " | " + queue["audit_reason"].fillna("")),
    })
    assert out["batch_id"].notna().all()
    assert not out["candidate_id"].duplicated().any()
    dest = MODEL_DIR / "data" / f"review_batch_{prefix}.csv"
    out.to_csv(dest, index=False)
    print("已写出：", dest)
    print(out["batch_id"].value_counts().to_string())

    # ---------- Save prior scores for comparison and preserve the original data ----------
    backup = queue[["sample_id", "context", "peak", "gene", "review_group",
                    "audit_tier", "audit_score", "audit_reason"]].copy()
    for c in common.HUMAN_SCORE_COLS:
        backup[f"old_{c}"] = backup["sample_id"].map(audit[c])
    backup["old_confidence"] = backup["sample_id"].map(audit["human_confidence"])
    backup["old_hard_label"] = backup["sample_id"].map(audit["hard_label"])
    backup_path = common.OUT_DIR / f"review_queue_{prefix}_old_scores_backup.csv"
    backup.to_csv(backup_path, index=False)
    print("已写出旧分备份：", backup_path)


if __name__ == "__main__":
    main()
