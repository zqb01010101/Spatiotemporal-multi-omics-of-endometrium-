#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""02_build_cv_splits.py — Group related samples and create repeated grouped CV splits

Outputs:
  data/group_mapping.csv       sample_id -> group_id
  data/cv_splits.csv           sample_id, group_id, repeat, fold
  outputs/cv_fold_summary.csv  Hard-label, context, and batch distributions by repeat and fold
  outputs/cv_meta.json         Split strategy, group counts, and seeds
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("02_build_cv_splits")


def main():
    audit = common.load_audit_dataset()
    audit["group_id"] = common.build_group_mapping(audit)
    audit[["sample_id", "group_id"]].to_csv(common.DATA_DIR / "group_mapping.csv", index=False)
    audit.to_csv(common.DATA_DIR / "audit_dataset.csv", index=False)   # Save group assignments for downstream scripts
    logger.info("样本数=%d，组数=%d", len(audit), audit["group_id"].nunique())

    n_folds = common.cfg("n_folds", 5)
    n_repeats = common.cfg("n_repeats", 5)
    seed = common.cfg("random_seed", 42)
    strat = audit["hard_label"].astype(str).to_numpy()
    groups = audit["group_id"].to_numpy()

    rows = []
    methods = {}
    for r in range(1, n_repeats + 1):
        splits, method = common.group_stratified_splitter(
            groups, strat, n_folds, seed + 1000 * r)
        methods[r] = method
        for f, (_, te_idx) in enumerate(splits):
            rows.append(pd.DataFrame({
                "sample_id": audit["sample_id"].iloc[te_idx].to_numpy(),
                "group_id": groups[te_idx],
                "repeat": r,
                "fold": f,
            }))
    cv = pd.concat(rows, ignore_index=True)
    cv.to_csv(common.DATA_DIR / "cv_splits.csv", index=False)

    # ---------- Summarize fold distributions ----------
    merged = cv.merge(audit[["sample_id", "context", "hard_label", "batch_id"]],
                      on="sample_id", how="left")
    merged["batch_id"] = merged["batch_id"].replace("", "none")
    summ = merged.groupby(["repeat", "fold"]).agg(
        n=("sample_id", "count"),
        hard_label_counts=("hard_label", lambda s: json.dumps(s.value_counts().to_dict(), ensure_ascii=False)),
        context_counts=("context", lambda s: json.dumps(s.value_counts().to_dict(), ensure_ascii=False)),
        batch_counts=("batch_id", lambda s: json.dumps(s.value_counts().to_dict(), ensure_ascii=False)),
    ).reset_index()
    summ.to_csv(common.OUT_DIR / "cv_fold_summary.csv", index=False)

    meta = {
        "n_samples": len(audit),
        "n_groups": int(audit["group_id"].nunique()),
        "n_folds": n_folds,
        "n_repeats": n_repeats,
        "seed": seed,
        "split_method_per_repeat": methods,
        "min_group_size": int(audit["group_id"].value_counts().min()),
        "max_group_size": int(audit["group_id"].value_counts().max()),
        "group_size_distribution": audit["group_id"].value_counts().value_counts().to_dict(),
    }
    with open(common.OUT_DIR / "cv_meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)
    logger.info("划分方式: %s", set(methods.values()))
    logger.info("已写出 cv_splits.csv / cv_fold_summary.csv / cv_meta.json")


if __name__ == "__main__":
    main()
