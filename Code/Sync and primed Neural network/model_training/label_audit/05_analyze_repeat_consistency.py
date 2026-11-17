#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""05_analyze_repeat_consistency.py — Check repeat annotation agreement

Compare each is_repeat=1 sample with its repeat_of original:
  repeat_score_mae / repeat_max_abs_difference / repeat_hard_label_conflict /
  repeat_confidence_diff / repeat_margin_diff
Output: outputs/repeat_consistency.csv (empty table if no repeats are available)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("05_analyze_repeat_consistency")


def main():
    audit = common.load_audit_dataset()
    reps = audit[audit["is_repeat"] == 1].copy()
    if not len(reps):
        empty = pd.DataFrame(columns=[
            "repeat_sample_id", "original_sample_id", "repeat_score_mae",
            "repeat_max_abs_difference", "repeat_hard_label_conflict",
            "repeat_confidence_diff", "repeat_margin_diff"])
        empty.to_csv(common.OUT_DIR / "repeat_consistency.csv", index=False)
        logger.info("无重复标注样本（is_repeat=1），已生成空 repeat_consistency.csv")
        return

    orig = audit.set_index("sample_id")
    rows = []
    for _, r in reps.iterrows():
        oid = r["repeat_of"]
        if oid not in orig.index:
            continue
        o = orig.loc[oid]
        diff = (r[common.HUMAN_SCORE_COLS].astype(float) -
                o[common.HUMAN_SCORE_COLS].astype(float)).to_numpy()
        rows.append({
            "repeat_sample_id": r["sample_id"],
            "original_sample_id": oid,
            "repeat_score_mae": float(np.abs(diff).mean()),
            "repeat_max_abs_difference": float(np.abs(diff).max()),
            "repeat_hard_label_conflict": int(r["hard_label"] != o["hard_label"]),
            "repeat_confidence_diff": float(r["human_confidence"] - o["human_confidence"]),
        })
    out = pd.DataFrame(rows)
    out.to_csv(common.OUT_DIR / "repeat_consistency.csv", index=False)
    logger.info("重复样本对 %d 组，mean repeat_score_mae=%.3f",
                len(out), out["repeat_score_mae"].mean() if len(out) else float("nan"))


if __name__ == "__main__":
    main()
