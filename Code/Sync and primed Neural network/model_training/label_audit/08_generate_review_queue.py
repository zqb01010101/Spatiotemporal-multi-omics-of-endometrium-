#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""08_generate_review_queue.py — Build full and blinded manual review queues

Defaults from config.json: top 50 Tier A, top 30 Tier B, and 20 random Tier D controls
Outputs:
  outputs/review_queue.csv            Full queue with prior scores, OOF predictions, and audit evidence
  outputs/review_queue_blinded.csv    Blinded queue with peak, gene, context, and 12 features
                                      Prior scores, predictions, and audit evidence are hidden.
  outputs/review_results.csv          Review template with prior scores filled in and new_* and decision blank
  figures/review_queue_blinded.pdf    Blinded plots (six candidates per page)
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

logger = common.get_logger("08_generate_review_queue")


def make_blinded_pdf(df, path):
    per_page = 6
    n = len(df)
    pages = int(np.ceil(n / per_page))
    with PdfPages(path) as pdf:
        for pg in range(pages):
            fig = plt.figure(figsize=(11.69, 8.27))
            fig.suptitle(f"Blinded review (page {pg + 1}/{pages})", fontsize=11)
            for i in range(per_page):
                idx = pg * per_page + i
                if idx >= n:
                    break
                row = df.iloc[idx]
                ax = fig.add_subplot(2, 3, i + 1)
                x = np.arange(len(common.STAGES))
                peak = row[common.FEATURES[:6]].to_numpy(dtype=float)
                gene = row[common.FEATURES[6:]].to_numpy(dtype=float)
                ax.plot(x, peak, "o-", color="#c0392b", lw=1.4, ms=3, label="peak")
                ax.plot(x, gene, "s-", color="#2471a3", lw=1.4, ms=3, label="gene")
                ax.set_xticks(x, common.STAGES, fontsize=6)
                ax.set_ylim(0, 1)
                ax.set_title(f"{row['peak']} | {row['gene']} | {row['context']}",
                             fontsize=7, loc="left")
                if i == 0:
                    ax.legend(fontsize=6, loc="upper left")
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)


def main():
    ranked = pd.read_csv(common.OUT_DIR / "label_audit_ranked.csv")
    audit = common.load_audit_dataset().set_index("sample_id")
    seed = common.cfg("random_seed", 42)
    n_a = common.cfg("tier_a_n", 50)
    n_b = common.cfg("tier_b_n", 30)
    n_ctrl = common.cfg("random_control_n", 20)

    tier_a = ranked[ranked["audit_tier"] == "A"].head(n_a)
    tier_b = ranked[ranked["audit_tier"] == "B"].head(n_b)
    pool_d = ranked[ranked["audit_tier"] == "D"]
    control = pool_d.sample(n=min(n_ctrl, len(pool_d)), random_state=seed)
    queue = pd.concat([tier_a, tier_b, control], ignore_index=True)
    queue["review_group"] = (["tier_a"] * len(tier_a) + ["tier_b"] * len(tier_b) +
                             ["random_control"] * len(control))
    logger.info("review queue: Tier A %d, Tier B %d, random control %d",
                len(tier_a), len(tier_b), len(control))

    feat = audit.loc[queue["sample_id"], common.FEATURES].reset_index()
    queue = queue.merge(feat, on="sample_id", how="left")
    queue.to_csv(common.OUT_DIR / "review_queue.csv", index=False)

    blinded_cols = ["sample_id", "context", "peak", "gene"] + common.FEATURES
    blinded = queue[blinded_cols].copy()
    blinded.to_csv(common.OUT_DIR / "review_queue_blinded.csv", index=False)
    make_blinded_pdf(blinded, common.FIG_DIR / "review_queue_blinded.pdf")

    # Review template: prefill prior scores and leave new_* and decision for manual entry
    tmpl_cols = (["sample_id", "context", "peak", "gene"] +
                 [f"old_{c}" for c in common.CLASSES] + ["old_confidence", "old_hard_label"] +
                 [f"new_{c}" for c in common.CLASSES] + ["new_confidence"] +
                 ["decision", "review_note"])
    tmpl = pd.DataFrame(index=range(len(queue)), columns=tmpl_cols)
    for i, sid in enumerate(queue["sample_id"]):
        r = audit.loc[sid]
        tmpl.loc[i, "sample_id"] = sid
        tmpl.loc[i, "context"] = r["context"]
        tmpl.loc[i, "peak"] = r["peak"]
        tmpl.loc[i, "gene"] = r["gene"]
        for c in common.CLASSES:
            tmpl.loc[i, f"old_{c}"] = r[c]
        tmpl.loc[i, "old_confidence"] = r["human_confidence"]
        tmpl.loc[i, "old_hard_label"] = r["hard_label"]
        tmpl.loc[i, "decision"] = "pending"
    tmpl.to_csv(common.OUT_DIR / "review_results.csv", index=False)
    logger.info("已写出 review_queue.csv / review_queue_blinded.csv / review_results.csv / blinded PDF")


if __name__ == "__main__":
    main()
