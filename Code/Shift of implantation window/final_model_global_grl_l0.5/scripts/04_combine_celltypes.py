"""Combine stroma + epithelial results (spot-count weighted) into one view."""

from __future__ import annotations

import pandas as pd
from sklearn.metrics import accuracy_score, f1_score

import common


def main():
    root = common.HERE / "results"
    out = root / "combined"
    out.mkdir(parents=True, exist_ok=True)

    epi = pd.read_csv(root / "epithelial" / "all_spot_predictions.csv.gz")
    strm = pd.read_csv(root / "stroma" / "all_spot_predictions.csv.gz")
    all_spots = pd.concat([epi, strm], ignore_index=True)
    all_spots.to_csv(out / "all_spot_predictions.csv.gz",
                     index=False, compression="gzip")
    common.sample_summary(all_spots, "classification").to_csv(
        out / "all_sample_predictions.csv", index=True)
    common.sample_summary(all_spots, "temporal_nearest").to_csv(
        out / "all_sample_predictions_temporal_nearest.csv", index=True)

    normal = all_spots[all_spots["analysis_condition"] == "normal"]
    po7 = normal[normal["true_stage"] == "sec_PO7"]
    sm = pd.read_csv(out / "all_sample_predictions.csv")
    normal_sm = sm[sm["analysis_condition"] == "normal"]

    # held-out pooled: combine the two cell-type pooled test sets
    pe = pd.read_csv(root / "epithelial" / "spot_predictions_pooled.csv.gz")
    ps = pd.read_csv(root / "stroma" / "spot_predictions_pooled.csv.gz")
    pooled = pd.concat([pe, ps], ignore_index=True)
    pooled_po7 = pooled[pooled["true_stage"] == "sec_PO7"]

    eval_json = {
        "scope": "combined epithelial + stroma (spot-count weighted)",
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
        "held_out_pooled_view": {
            "n_spots": int(len(pooled)),
            "spot_accuracy": float(accuracy_score(pooled["true_stage"],
                                                  pooled["predicted_stage"])),
            "macro_f1": float(f1_score(pooled["true_stage"],
                                       pooled["predicted_stage"],
                                       average="macro", zero_division=0)),
            "sample_accuracy": float(accuracy_score(
                pooled.groupby("sample_id")["true_stage"].first(),
                pooled.groupby("sample_id")["predicted_stage"].agg(
                    lambda s: s.mode().iloc[0] if len(s.mode()) else None))),
            "PO7_to_PO9": float((pooled_po7["predicted_stage"] == "sec_PO9").mean())
            if len(pooled_po7) else None,
        },
    }
    common.utils.save_json(eval_json, out / "combined_evaluation.json")
    lines = [
        "# Combined epithelial + stroma results (Global GRL lambda=0.5)",
        "",
        "Spot-count weighted by concatenating the two cell-type spot sets.",
        "",
        "## Deployment view",
        f"- n spots: {eval_json['deployment_view']['n_spots']}",
        f"- spot accuracy: {eval_json['deployment_view']['spot_accuracy']:.4f}",
        f"- spot macro-F1: {eval_json['deployment_view']['spot_macro_f1']:.4f}",
        f"- sample accuracy: {eval_json['deployment_view']['sample_accuracy']:.4f}",
        f"- PO7->PO9: {eval_json['deployment_view']['PO7_to_PO9']:.4f}",
        "",
        "## Held-out pooled view",
        f"- n spots: {eval_json['held_out_pooled_view']['n_spots']}",
        f"- spot accuracy: {eval_json['held_out_pooled_view']['spot_accuracy']:.4f}",
        f"- macro-F1: {eval_json['held_out_pooled_view']['macro_f1']:.4f}",
        f"- sample accuracy: {eval_json['held_out_pooled_view']['sample_accuracy']:.4f}",
        f"- PO7->PO9: {eval_json['held_out_pooled_view']['PO7_to_PO9']:.4f}",
    ]
    (out / "combined_RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
