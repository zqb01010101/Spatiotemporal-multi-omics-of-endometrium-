#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""13_autofill_review_results.py — Fill review_results.csv from updated app scores

Inputs:
  outputs/review_queue_old_scores_backup.csv   Snapshot of prior scores
  peak_gene_label_app/labels/labeled_training_data.csv   Latest export with revised scores
Outputs:
  outputs/review_results.csv   old_*, new_*, decision, and review_note for each candidate
    Set decision to keep_original if all scores and confidence are unchanged. Otherwise use replace_with_new.
    Manually verify or edit the decisions, then rerun step 09.
  outputs/review_change_summary.csv   Summary of changes
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODEL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common


def main():
    backup = pd.read_csv(common.OUT_DIR / "review_queue_old_scores_backup.csv")
    labeled = pd.read_csv(common.APP_LABELS / "labeled_training_data.csv")
    new = labeled.set_index("candidate_id")

    rows = []
    for _, r in backup.iterrows():
        sid = r["sample_id"]
        if sid not in new.index:
            print(f"!! {sid} 不在最新导出中，跳过")
            continue
        nr = new.loc[sid]
        old = [r[f"old_{c}"] for c in common.CLASSES]
        nw = [nr[f"score_{c}"] for c in common.CLASSES]
        old_conf = r["old_confidence"]
        new_conf = nr["overall_confidence"]
        changed = (not np.allclose(old, nw, atol=1e-9)) or (not np.allclose(old_conf, new_conf, atol=1e-9))
        decision = "replace_with_new" if changed else "keep_original"
        new_hard = common.human_hard_label(np.array(nw), float(new_conf))
        rows.append({
            "sample_id": sid,
            "context": r["context"],
            "peak": r["peak"],
            "gene": r["gene"],
            **{f"old_{c}": old[i] for i, c in enumerate(common.CLASSES)},
            "old_confidence": old_conf,
            "old_hard_label": r["old_hard_label"],
            **{f"new_{c}": nw[i] for i, c in enumerate(common.CLASSES)},
            "new_confidence": new_conf,
            "new_hard_label": new_hard,
            "decision": decision,
            "review_note": "platform re-score" if changed else "unchanged",
        })

    results = pd.DataFrame(rows)
    results.to_csv(common.OUT_DIR / "review_results.csv", index=False)
    print("review_results.csv 已生成：", len(results), "条")
    print(results["decision"].value_counts().to_string())

    summary = results.groupby("review_group" if "review_group" in results else "decision").size()
    changed = results[results["decision"] == "replace_with_new"]
    if len(changed):
        hl_flip = (changed["old_hard_label"] != changed["new_hard_label"])
        print("hard_label 变化的条数：", int(hl_flip.sum()))
        print(changed["old_hard_label"].value_counts().to_string())
        print("→ new：\n", changed["new_hard_label"].value_counts().to_string())
    changed.to_csv(common.OUT_DIR / "review_change_summary.csv", index=False)


if __name__ == "__main__":
    main()
