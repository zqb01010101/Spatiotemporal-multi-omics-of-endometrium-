#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""14_compare_architectures.py — Train and compare two deeper networks

Architecture A: 256 -> 128 -> 64 -> 32
Architecture B: 256 -> 64 -> 32

Both use the same raw inputs, loss weights (5:5:3.5:2), dropout keep probability (0.8), seed, and early stopping.
Outputs:
  outputs/model/arch_256_128_64_32/   model.npz + metrics.json
  outputs/model/arch_256_64_32/       model.npz + metrics.json
  outputs/architecture_comparison.csv
  figures/architecture_comparison.png
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, mean_squared_error

MODEL_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = MODEL_DIR / "data"
OUT_DIR = MODEL_DIR / "outputs"
FIG_DIR = OUT_DIR / "figures"
MODEL_OUT = OUT_DIR / "model"

sys.path.insert(0, str(MODEL_DIR / "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("train_mod", MODEL_DIR / "scripts" / "02_train.py")
train_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train_mod)

ARCHS = {
    "256_128_64_32": [256, 128, 64, 32],
    "256_64_32": [256, 64, 32],
}
SEED = 42


def param_count(hidden, n_in=12, n_out=4):
    sizes = [n_in] + hidden + [n_out]
    return sum(sizes[i] * sizes[i + 1] + sizes[i + 1] for i in range(len(sizes) - 1))


def evaluate_arch(name, hidden):
    data = np.load(DATA_DIR / "train_test.npz")
    X_tr, Y_tr, hard_tr = data["X_train"], data["Y_train"], data["hard_train"]
    X_val, Y_val, hard_val = data["X_val"], data["Y_val"], data["hard_val"]

    model = train_mod.ScoreMLP(
        X_tr.shape[1], hidden, Y_tr.shape[1],
        seed=SEED, keep_prob=train_mod.DROPOUT_KEEP,
        loss_weights=train_mod.OUTPUT_LOSS_WEIGHTS)
    hist, best_val = train_mod.train(model, X_tr, Y_tr, X_val, Y_val)

    pred = model.forward(X_val) * 10.0
    y = Y_val * 10.0
    rmse = {f"rmse_{c}": float(np.sqrt(mean_squared_error(y[:, i], pred[:, i])))
            for i, c in enumerate(train_mod.CLASSES)}
    corr = {f"corr_{c}": float(np.corrcoef(y[:, i], pred[:, i])[0, 1])
            for i, c in enumerate(train_mod.CLASSES)}

    pred_label, _, _, _ = train_mod.derive_hard_label(pred)
    pred_label_str = np.array([train_mod.CLASSES[i] if i >= 0 else "uncertain"
                               for i in pred_label])
    known = hard_val >= 0
    acc = float(accuracy_score(hard_val[known], pred_label[known]))
    bal = float(balanced_accuracy_score(hard_val[known], pred_label[known]))

    out_dir = MODEL_OUT / f"arch_{name}"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / "model.npz", **model.params())
    metrics = {
        "name": name,
        "hidden": hidden,
        "n_params": param_count(hidden),
        "best_val_weighted_mse": float(best_val),
        "epochs_used": len(hist["train"]),
        "hard_accuracy": acc,
        "hard_balanced": bal,
        **rmse, **corr,
    }
    with open(out_dir / "metrics.json", "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, ensure_ascii=False, indent=2)
    print(f"\n=== {name} === params={metrics['n_params']} epochs={metrics['epochs_used']}")
    print(f"  val weighted MSE = {best_val:.6f}")
    for c in train_mod.CLASSES:
        print(f"  {c:12s} RMSE={metrics[f'rmse_{c}']:.4f}  corr={metrics[f'corr_{c}']:.4f}")
    print(f"  hard_label acc={acc:.4f} balanced={bal:.4f}")
    return metrics


def main():
    rows = []
    for name, hidden in ARCHS.items():
        rows.append(evaluate_arch(name, hidden))
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT_DIR / "architecture_comparison.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for i, (key, label) in enumerate([("rmse_", "RMSE (lower better)"),
                                      ("corr_", "Pearson corr (higher better)")]):
        cols = [c for c in comp.columns if c.startswith(key)]
        x = np.arange(len(cols))
        w = 0.38
        axes[i].bar(x - w / 2, comp.iloc[0][cols], w, label=comp.iloc[0]["name"])
        axes[i].bar(x + w / 2, comp.iloc[1][cols], w, label=comp.iloc[1]["name"])
        axes[i].set_xticks(x, [c.split("_", 1)[1] for c in cols], rotation=20)
        axes[i].set_title(label)
        axes[i].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "architecture_comparison.png", dpi=130)
    plt.close(fig)
    print("\n已写出 outputs/architecture_comparison.csv 与 figures/architecture_comparison.png")


if __name__ == "__main__":
    main()
