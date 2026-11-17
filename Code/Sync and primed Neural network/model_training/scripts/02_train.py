#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""02_train.py — Train an MLP to regress four independent scores

Network architecture (v3):
  12 inputs (six peak and six gene features, as raw fractions without z-scoring)
    -> Dense 128 (ReLU)
    -> Dense 64  (ReLU)
    -> Dense 32  (ReLU)
    -> Four independent sigmoid scores in [0, 1], multiplied by 10 for reporting
  Confidence is the maximum predicted score. Also report the gap between the top two scores.
  v4: dropout after each hidden-layer ReLU (keep probability 0.8)
  v4.1: use unscaled fractions so flat profiles across stages retain their ambiguity
  v4.2: weight the four losses 5:5:3.5:2 to downweight unrelated and sync_repress
  v5: hidden layer widths of 256, 128, 64, and 32

Use weighted MSE across the four outputs, Adam, and early stopping.

Outputs:
  outputs/model/model.npz, scaler.npz, metrics.json, class_map.json
  outputs/figures/（loss_curve, confusion_matrix, score_calibration）
  data/queue_predictions.csv
"""

import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             confusion_matrix, classification_report)

MODEL_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = MODEL_DIR / "data"
OUT_DIR = MODEL_DIR / "outputs"
MODEL_OUT = OUT_DIR / "model"
FIG_OUT = OUT_DIR / "figures"
for d in (MODEL_OUT, FIG_OUT):
    d.mkdir(parents=True, exist_ok=True)

HIDDEN = [256, 128, 64, 32]
DROPOUT_KEEP = 0.8          # Hidden-layer keep probability (dropout = 0.2)
LR = 1e-3
EPOCHS = 500
PATIENCE = 40
BATCH = 32
RANDOM_STATE = 42
CLASSES = ["pre_open", "sync_open", "sync_repress", "unrelated"]
SCORE_OUT = ["score_pre_open", "score_sync_open", "score_sync_repress", "score_unrelated"]
OUTPUT_LOSS_WEIGHTS = np.array([5.0, 5.0, 3.5, 2.0])   # pre_open:sync_open:sync_repress:unrelated


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def relu(x):
    return np.maximum(0.0, x)


class ScoreMLP:
    """Three hidden layers with dropout and four independent sigmoid outputs"""

    def __init__(self, n_in, hidden, n_out, seed=0, keep_prob=1.0, loss_weights=None):
        self.rng = np.random.default_rng(seed)
        self.keep_prob = keep_prob
        if loss_weights is None:
            self.loss_weights = np.ones(n_out)
        else:
            self.loss_weights = np.asarray(loss_weights, dtype=float)
            self.loss_weights = self.loss_weights / self.loss_weights.mean()
        sizes = [n_in] + list(hidden) + [n_out]
        self.W = []
        self.b = []
        for a, c in zip(sizes[:-1], sizes[1:]):
            self.W.append(self.rng.standard_normal((a, c)) * np.sqrt(2.0 / a))
            self.b.append(np.zeros(c))
        self.act = []   # Save each layer's activations and the inputs for backpropagation
        self.masks = [] # Save hidden-layer dropout masks for the backward pass

    def forward(self, X, training=False):
        self.act = [X]
        self.masks = []
        h = X
        for W, b in zip(self.W[:-1], self.b[:-1]):
            h = relu(h @ W + b)
            if training and self.keep_prob < 1.0:
                m = (self.rng.random(h.shape) < self.keep_prob) / self.keep_prob
                h = h * m
                self.masks.append(m)
            else:
                self.masks.append(None)
            self.act.append(h)
        z = self.act[-1] @ self.W[-1] + self.b[-1]
        self.out = sigmoid(z)
        return self.out

    def backward(self, X, Y):
        self.forward(X, training=True)
        n = X.shape[0]
        dout = 2.0 * (self.out - Y) * self.loss_weights / n
        dz = dout * self.out * (1.0 - self.out)
        dW = [None] * len(self.W)
        db = [None] * len(self.b)
        dW[-1] = self.act[-1].T @ dz
        db[-1] = dz.sum(axis=0)
        for k in range(len(self.W) - 2, -1, -1):
            da = dz @ self.W[k + 1].T
            da[self.act[k + 1] <= 0] = 0.0
            if self.masks[k] is not None:
                da *= self.masks[k]
            dz = da
            dW[k] = self.act[k].T @ dz
            db[k] = dz.sum(axis=0)
        return dW, db

    def params(self):
        return {"W" + str(i): w for i, w in enumerate(self.W)} | \
               {"b" + str(i): b for i, b in enumerate(self.b)}


class Adam:
    def __init__(self, lr=1e-3, beta1=0.9, beta2=0.999, eps=1e-8):
        self.lr, self.b1, self.b2, self.eps = lr, beta1, beta2, eps
        self.m, self.v, self.t = {}, {}, 0

    def step(self, params, grads):
        self.t += 1
        for k, p in params.items():
            g = grads[k]
            m = self.m.get(k, np.zeros_like(p))
            v = self.v.get(k, np.zeros_like(p))
            m = self.b1 * m + (1 - self.b1) * g
            v = self.b2 * v + (1 - self.b2) * g * g
            mh = m / (1 - self.b1 ** self.t)
            vh = v / (1 - self.b2 ** self.t)
            params[k] -= self.lr * mh / (np.sqrt(vh) + self.eps)
            self.m[k], self.v[k] = m, v


def mse(y_true, y_pred):
    return float(np.mean((y_true - y_pred) ** 2))


def weighted_mse(y_true, y_pred, weights):
    return float(np.mean(((y_true - y_pred) ** 2) * weights))


def train(model, X_tr, Y_tr, X_val, Y_val):
    adam = Adam(lr=LR)
    best = None
    best_val = np.inf
    patience = 0
    hist = {"train": [], "val": []}
    n = len(Y_tr)
    rng = np.random.default_rng(RANDOM_STATE)

    for epoch in range(EPOCHS):
        idx = rng.permutation(n)
        for start in range(0, n, BATCH):
            sel = idx[start:start + BATCH]
            dW, db = model.backward(X_tr[sel], Y_tr[sel])
            grads = {("W" + str(i)): w for i, w in enumerate(dW)}
            grads.update({("b" + str(i)): b for i, b in enumerate(db)})
            adam.step(model.params(), grads)

        l_tr = weighted_mse(Y_tr, model.forward(X_tr), model.loss_weights)
        l_val = weighted_mse(Y_val, model.forward(X_val), model.loss_weights)
        hist["train"].append(l_tr)
        hist["val"].append(l_val)

        if l_val < best_val - 1e-5:
            best_val = l_val
            best = {k: v.copy() for k, v in model.params().items()}
            patience = 0
        else:
            patience += 1
            if patience >= PATIENCE:
                print(f"early stop @ epoch {epoch + 1}, best val MSE = {best_val:.6f}")
                break
    for k, v in best.items():
        if k.startswith("W"):
            model.W[int(k[1:])] = v
        else:
            model.b[int(k[1:])] = v
    return hist, best_val


def derive_hard_label(scores10):
    """For scores on the 0–10 scale, assign the top class if its score is ≥7 and leads by ≥2. Otherwise assign uncertain (-1)."""
    order = np.argsort(-scores10, axis=1)
    top = scores10[np.arange(len(scores10)), order[:, 0]]
    second = scores10[np.arange(len(scores10)), order[:, 1]]
    margin = top - second
    label = np.where((top >= 7.0) & (margin >= 2.0), order[:, 0], -1)
    return label, top, second, margin


def main():
    data = np.load(DATA_DIR / "train_test.npz")
    X_train, Y_train = data["X_train"], data["Y_train"]
    X_val, Y_val = data["X_val"], data["Y_val"]
    hard_val = data["hard_val"]
    X_queue = data["X_queue"]
    with open(DATA_DIR / "class_meta.json", encoding="utf-8") as fh:
        meta = json.load(fh)

    # Feed raw fractions directly to the model.
    # Stage-wise z-scoring can distort flat profiles
    # and obscure their ambiguity.
    X_train_s = X_train
    X_val_s = X_val
    X_queue_s = X_queue

    model = ScoreMLP(X_train.shape[1], HIDDEN, Y_train.shape[1],
                     seed=RANDOM_STATE, keep_prob=DROPOUT_KEEP,
                     loss_weights=OUTPUT_LOSS_WEIGHTS)
    hist, best_val = train(model, X_train_s, Y_train, X_val_s, Y_val)

    # ---------- Evaluate regression ----------
    pred = model.forward(X_val_s)
    rmse = np.sqrt(np.mean((Y_val - pred) ** 2, axis=0))
    mae = np.mean(np.abs(Y_val - pred), axis=0)
    corr = [float(np.corrcoef(Y_val[:, k], pred[:, k])[0, 1]) for k in range(4)]
    print("\n==== 验证集回归评估（分数 0-1 尺度） ====")
    for k, c in enumerate(["pre_open", "sync_open", "sync_repress", "unrelated"]):
        print(f"  {c:12s} RMSE={rmse[k]:.4f}  MAE={mae[k]:.4f}  corr={corr[k]:.4f}")

    # ---------- Evaluate derived hard labels ----------
    pred_label, ptop, psec, pmargin = derive_hard_label(pred * 10.0)
    true_label = hard_val
    known = true_label >= 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        acc = accuracy_score(true_label[known], pred_label[known]) if known.sum() else float("nan")
        bal = balanced_accuracy_score(true_label[known], pred_label[known]) if known.sum() else float("nan")
    print(f"hard_label（预测分 vs 标注分）：accuracy={acc:.4f}  balanced={bal:.4f}  "
          f"(评估样本 {int(known.sum())}，其中预测 uncertain {int((pred_label[known] < 0).sum())})")
    pred_label_vis = np.where(pred_label >= 0, pred_label, 4)   # 4 denotes predicted uncertain
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cm = confusion_matrix(true_label[known], pred_label_vis[known], labels=list(range(5)))
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(cm, cmap="Blues")
    tick_labels = CLASSES + ["uncertain"]
    ax.set_xticks(range(5), tick_labels, rotation=30)
    ax.set_yticks(range(5), tick_labels)
    for i in range(5):
        for j in range(5):
            ax.text(j, i, cm[i, j], ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set_xlabel("Predicted hard label")
    ax.set_ylabel("True hard label")
    ax.set_title(f"Hard-label confusion (val, acc={acc:.3f})")
    fig.colorbar(im)
    fig.tight_layout()
    fig.savefig(FIG_OUT / "confusion_matrix.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, len(hist["train"]) + 1), hist["train"], label="train")
    ax.plot(range(1, len(hist["val"]) + 1), hist["val"], label="val")
    ax.set_xlabel("epoch")
    ax.set_ylabel("MSE")
    ax.legend()
    ax.set_title("Loss curves")
    fig.tight_layout()
    fig.savefig(FIG_OUT / "loss_curve.png", dpi=150)
    plt.close(fig)

    # Calibration plot per output: mean observed vs. mean predicted score
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    for k, ax in enumerate(axes.ravel()):
        bins = np.linspace(0, 1, 11)
        centers, truth, counts = [], [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = (pred[:, k] >= lo) & (pred[:, k] < hi)
            if m.sum() > 0:
                centers.append((lo + hi) / 2)
                truth.append(Y_val[m, k].mean())
                counts.append(int(m.sum()))
        ax.plot(centers, truth, "o-")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
        ax.set_title(CLASSES[k])
        ax.set_xlabel("predicted")
        ax.set_ylabel("true")
    fig.suptitle("Score calibration (val)")
    fig.tight_layout()
    fig.savefig(FIG_OUT / "score_calibration.png", dpi=150)
    plt.close(fig)

    # ---------- Cross-check with scikit-learn ----------
    from sklearn.neural_network import MLPRegressor
    sk = MLPRegressor(hidden_layer_sizes=tuple(HIDDEN), activation="relu",
                      max_iter=500, early_stopping=True, validation_fraction=0.15,
                      random_state=RANDOM_STATE)
    sk.fit(X_train_s, Y_train)
    sk_rmse = float(np.sqrt(np.mean((Y_val - sk.predict(X_val_s)) ** 2)))
    print(f"[核对] sklearn MLPRegressor 验证集 RMSE = {sk_rmse:.4f}")

    # ---------- Save outputs ----------
    params = model.params()
    np.savez(MODEL_OUT / "model.npz", **params)
    with open(MODEL_OUT / "class_map.json", "w", encoding="utf-8") as fh:
        json.dump({"classes": CLASSES, "score_outputs": SCORE_OUT}, fh, ensure_ascii=False)
    metrics = {
        "rmse_per_output": dict(zip(SCORE_OUT, np.round(rmse, 5).tolist())),
        "mae_per_output": dict(zip(SCORE_OUT, np.round(mae, 5).tolist())),
        "corr_per_output": dict(zip(SCORE_OUT, np.round(corr, 5).tolist())),
        "hardlabel_accuracy": acc,
        "hardlabel_balanced_accuracy": bal,
        "val_mse": best_val,
        "sklearn_mlp_rmse": sk_rmse,
        "n_train": int(len(Y_train)),
        "n_val": int(len(Y_val)),
        "architecture": {"input": X_train.shape[1], "hidden": HIDDEN, "output": 4,
                         "output_activation": "sigmoid (no softmax)",
                         "dropout_keep_prob": DROPOUT_KEEP,
                         "normalization": "none (raw 0-1 inputs)",
                         "loss_weights_ratio": {"pre_open": 5, "sync_open": 5,
                                                "sync_repress": 3.5, "unrelated": 2}},
        "epochs_used": len(hist["train"]),
    }
    with open(MODEL_OUT / "metrics.json", "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, ensure_ascii=False, indent=2)
    print("模型已保存：", MODEL_OUT)

    # ---------- Score the full labeling queue ----------
    p_queue = model.forward(X_queue_s) * 10.0
    queue = pd.read_csv(DATA_DIR / "queue_features_clean.csv")
    pred_df = queue[["candidate_id", "context", "peak", "gene", "distance"]].copy()
    for i, c in enumerate(SCORE_OUT):
        pred_df[c] = p_queue[:, i]
    qlabel, qtop, qsec, qmargin = derive_hard_label(p_queue)
    pred_df["pred_class"] = [CLASSES[i] if i >= 0 else "uncertain" for i in qlabel]
    pred_df["confidence"] = qtop / 10.0
    pred_df["margin"] = qmargin / 10.0
    pred_df = pred_df.sort_values("confidence").reset_index(drop=True)
    pred_df.to_csv(DATA_DIR / "queue_predictions.csv", index=False)
    print(f"队列预测已保存：data/queue_predictions.csv（{len(pred_df)} 行）")
    print("预测 hard_label 分布：")
    print(pred_df["pred_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
