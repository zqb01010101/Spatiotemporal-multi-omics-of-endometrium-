"""Evaluation summary for the final Global GRL model."""

from __future__ import annotations

import argparse

import pandas as pd
from sklearn.metrics import accuracy_score, f1_score

import common


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-type", required=True, choices=["stroma", "epithelial"])
    args = p.parse_args()
    out = common.HERE / "results" / args.cell_type
    df = pd.read_csv(out / "all_spot_predictions.csv.gz")
    normal = df[df["analysis_condition"] == "normal"]
    po7 = normal[normal["true_stage"] == "sec_PO7"]
    sm = pd.read_csv(out / "all_sample_predictions.csv")
    normal_sm = sm[sm["analysis_condition"] == "normal"]
    held = common.utils.load_json(out / "metrics_pooled.json")
    eval_json = {
        "cell_type": args.cell_type,
        "deployment_view": {
            "n_spots": int(len(normal)), "n_samples": int(len(normal_sm)),
            "spot_accuracy": float(accuracy_score(normal["true_stage"],
                                                   normal["predicted_stage"])),
            "spot_macro_f1": float(f1_score(normal["true_stage"],
                                            normal["predicted_stage"],
                                            average="macro", zero_division=0)),
            "sample_accuracy": float(accuracy_score(normal_sm["true_stage"],
                                                     normal_sm["majority_vote_stage"])),
            "PO7_to_PO9": float((po7["predicted_stage"] == "sec_PO9").mean()),
        },
        "held_out_pooled": held,
    }
    common.utils.save_json(eval_json, out / "evaluation.json")
    lines = [
        f"# Final Global GRL lambda=0.5 results ({args.cell_type})",
        "",
        "## Deployment view",
        f"- spot accuracy: {eval_json['deployment_view']['spot_accuracy']:.4f}",
        f"- spot macro-F1: {eval_json['deployment_view']['spot_macro_f1']:.4f}",
        f"- sample accuracy: {eval_json['deployment_view']['sample_accuracy']:.4f}",
        f"- PO7->PO9: {eval_json['deployment_view']['PO7_to_PO9']:.4f}",
        "",
        "## Held-out pooled view",
        f"- spot accuracy: {held['accuracy']:.4f}",
        f"- macro-F1: {held['macro_F1']:.4f}",
        f"- sample accuracy: {held.get('sample_accuracy_majority_vote', float('nan')):.4f}",
        f"- PO7->PO9: {held.get('PO7_to_PO9_error_rate', float('nan')):.4f}",
    ]
    (out / "RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
