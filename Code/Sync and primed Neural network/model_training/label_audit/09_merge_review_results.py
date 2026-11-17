#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""09_merge_review_results.py — Merge manual review decisions while preserving original labels

Read outputs/review_results.csv with manually entered new_* and decision:
  decision ∈ {keep_original, replace_with_new, mark_uncertain, exclude, pending}
Outputs:
  outputs/labels_clean_v1.csv   Revised training set with original_label, reviewed_label,
                                review_status, and revision_timestamp. Original data are preserved.
  outputs/review_results_merged.csv
"""

import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("09_merge_review_results")


def main():
    results = pd.read_csv(common.OUT_DIR / "review_results.csv")
    audit = common.load_audit_dataset().set_index("sample_id")

    clean = audit.copy()
    clean["original_label"] = clean["hard_label"]
    clean["reviewed_label"] = clean["hard_label"]
    clean["review_status"] = "not_reviewed"
    clean["revision_timestamp"] = ""

    now = datetime.now().isoformat(timespec="seconds")
    n_apply = 0
    for _, r in results.iterrows():
        sid = r["sample_id"]
        if sid not in clean.index:
            logger.warning("review_results 中的 sample_id 不在 audit dataset: %s", sid)
            continue
        decision = str(r.get("decision", "pending")).strip().lower()
        clean.loc[sid, "review_status"] = decision
        clean.loc[sid, "revision_timestamp"] = now
        if decision == "keep_original":
            clean.loc[sid, "reviewed_label"] = r.get("old_hard_label", clean.loc[sid, "hard_label"])
        elif decision == "replace_with_new":
            for c in common.CLASSES:
                v = r.get(f"new_{c}")
                if pd.notna(v):
                    clean.loc[sid, c] = float(v)
            conf = r.get("new_confidence")
            if pd.notna(conf):
                clean.loc[sid, "human_confidence"] = float(conf)
            scores = clean.loc[sid, common.HUMAN_SCORE_COLS].to_numpy(dtype=float)
            clean.loc[sid, "hard_label"] = common.human_hard_label(
                scores, clean.loc[sid, "human_confidence"])
            clean.loc[sid, "reviewed_label"] = clean.loc[sid, "hard_label"]
            n_apply += 1
        elif decision == "mark_uncertain":
            clean.loc[sid, "hard_label"] = "uncertain"
            clean.loc[sid, "reviewed_label"] = "uncertain"
        elif decision == "exclude":
            clean.loc[sid, "review_status"] = "exclude"
            clean.loc[sid, "reviewed_label"] = "excluded"
        elif decision == "pending":
            pass
        else:
            logger.warning("未知 decision: %s（%s）", decision, sid)

    clean = clean[clean["review_status"] != "exclude"]
    clean.reset_index().to_csv(common.OUT_DIR / "labels_clean_v1.csv", index=False)
    results.to_csv(common.OUT_DIR / "review_results_merged.csv", index=False)
    logger.info("已写出 labels_clean_v1.csv（%d 行，应用 replace_with_new %d 条）",
                len(clean), n_apply)
    logger.info("review_status 分布:\n%s", clean["review_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
