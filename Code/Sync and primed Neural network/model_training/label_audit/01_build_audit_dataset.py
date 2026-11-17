#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""01_build_audit_dataset.py — Build a unified audit dataset

Read-only inputs:
  peak_gene_label_app/labels/labeled_training_data.csv  (12 features, four scores, confidence, and hard_label)
  peak_gene_label_app/labels/annotations_export.csv     (batch_id, is_repeat, repeat_of, and source)
  model_training/data/annotated_samples.csv             (optional reference values such as argmax)

Outputs:
  data/audit_dataset.csv            One row per labeled sample
  data/audit_dataset_summary.json   Counts, distributions, anomalies, and input hashes
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("01_build_audit_dataset")


def main():
    labeled = pd.read_csv(common.APP_LABELS / "labeled_training_data.csv")
    ann_export = pd.read_csv(common.APP_LABELS / "annotations_export.csv")
    logger.info("labeled_training_data.csv 行数: %d", len(labeled))
    logger.info("annotations_export.csv 行数: %d", len(ann_export))

    # ---------- Join batch, repeat, and source metadata ----------
    meta_cols = ["candidate_id", "batch_id", "is_repeat", "repeat_of",
                 "mandatory", "source_tri_pre", "source_sync_open",
                 "source_sync_repress_heuristic", "queue_priority"]
    meta = ann_export[meta_cols].drop_duplicates("candidate_id")
    df = labeled.merge(meta, on="candidate_id", how="left")

    # ---------- Standardize column names ----------
    score_map = {f"score_{c}": c for c in common.CLASSES}
    df = df.rename(columns={**score_map, "overall_confidence": "human_confidence"})
    df = df.rename(columns={"hard_label": "hard_label"})

    def derive_source(r):
        if pd.notna(r.get("batch_id")):
            return str(r["batch_id"])
        if int(r.get("mandatory") or 0) == 1:
            return "mandatory"
        if int(r.get("source_sync_repress_heuristic") or 0) == 1:
            return "sync_repress_heuristic"
        if pd.notna(r.get("queue_priority")) and int(r["queue_priority"]) == 3:
            return "background"
        return "other"

    df["source"] = df.apply(derive_source, axis=1)
    df["is_repeat"] = df["is_repeat"].fillna(0).astype(int)
    df["repeat_of"] = df["repeat_of"].where(df["repeat_of"].notna(), "")
    df["batch_id"] = df["batch_id"].where(df["batch_id"].notna(), "")

    keep_cols = (["candidate_id", "context", "peak", "gene"] + common.FEATURES +
                 common.HUMAN_SCORE_COLS + ["human_confidence", "hard_label",
                                            "batch_id", "source", "is_repeat", "repeat_of"])
    audit = df[keep_cols].copy()
    audit = audit.rename(columns={"candidate_id": "sample_id"})
    audit = audit.sort_values("sample_id").reset_index(drop=True)

    # ---------- Log anomalous scores without changing the data ----------
    scores = audit[common.HUMAN_SCORE_COLS].to_numpy(dtype=float)
    anomalies = {}
    anomalies["missing_score"] = int(pd.isna(scores).any(axis=1).sum())
    anomalies["score_out_of_range"] = int(((scores < 0) | (scores > 10)).any(axis=1).sum())
    anomalies["all_zero_scores"] = int((scores.sum(axis=1) == 0).sum())
    recomputed = [
        common.human_hard_label(scores[i], audit["human_confidence"].iloc[i])
        for i in range(len(audit))
    ]
    anomalies["hard_label_mismatch"] = int(
        (np.array(recomputed) != audit["hard_label"].to_numpy()).sum())
    audit["recomputed_hard_label"] = recomputed
    audit["hard_label_mismatch"] = (
        audit["recomputed_hard_label"] != audit["hard_label"]).astype(int)
    audit["dup_sample_id"] = audit["sample_id"].duplicated().astype(int)

    logger.info("异常评分统计: %s", json.dumps(anomalies, ensure_ascii=False))

    # ---------- Write outputs ----------
    audit.to_csv(common.DATA_DIR / "audit_dataset.csv", index=False)
    summary = {
        "generated_at": pd.Timestamp.now().isoformat(),
        "n_samples": int(len(audit)),
        "n_groups": int(audit["group_id"].nunique()) if "group_id" in audit else None,
        "context_counts": audit["context"].value_counts().to_dict(),
        "hard_label_counts": audit["hard_label"].value_counts().to_dict(),
        "batch_counts": audit["batch_id"].replace("", "none").value_counts().to_dict(),
        "n_repeats": int((audit["is_repeat"] == 1).sum()),
        "anomalies": anomalies,
        "n_dup_sample_id": int(audit["dup_sample_id"].sum()),
        "score_stats": {
            c: {"mean": float(audit[c].mean()), "sd": float(audit[c].std()),
                "min": float(audit[c].min()), "max": float(audit[c].max())}
            for c in common.HUMAN_SCORE_COLS
        },
        "input_hashes": {
            "labeled_training_data.csv": common.sha256_file(
                common.APP_LABELS / "labeled_training_data.csv"),
            "annotations_export.csv": common.sha256_file(
                common.APP_LABELS / "annotations_export.csv"),
        },
    }
    with open(common.DATA_DIR / "audit_dataset_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    logger.info("已写出 data/audit_dataset.csv（%d 行）与 data/audit_dataset_summary.json", len(audit))


if __name__ == "__main__":
    main()
