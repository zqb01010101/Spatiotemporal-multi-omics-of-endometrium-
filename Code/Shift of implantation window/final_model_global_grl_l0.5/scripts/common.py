"""Shared helpers for the final Global GRL (lambda=0.5) template."""

from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE.parent / "stage_prediction"))

from src import config as scfg, utils  # noqa: E402


CFG = scfg.load_config(HERE / "config.yaml")
CFG["_project_root"] = str(HERE)
CFG["random_seed"] = int(CFG["random_seed"])

CLASSES = ["pro_stage", "sec_PO3", "sec_PO5", "sec_PO7",
           "sec_PO9", "sec_PO11", "sec_PO13"]
TIMES = np.array([-1, 3, 5, 7, 9, 11, 13], dtype=float)


def load_final_model(cell_type: str):
    model_dir = HERE / "models" / cell_type
    meta = utils.load_json(model_dir / "metadata.json")
    from src.model_method_invariant import MethodInvariantStageMLP

    model = MethodInvariantStageMLP(
        input_dim=int(meta["input_dim"]), n_classes=len(meta["classes"]),
        hidden=CFG["hidden_dim"], latent=CFG["latent_dim"],
        dropout=CFG["dropout"], conditional=bool(meta.get("conditional", False)))
    model.load_state_dict(torch.load(model_dir / "model.pt", map_location="cpu"))
    model.to(utils.device()).eval()
    scaler = joblib.load(model_dir / "scaler.joblib")
    return model, scaler, meta


def predict_ct(cfg, model, scaler, meta, cell_type: str, batch: int = 4096):
    import anndata as ad

    adata = ad.read_h5ad(cfg["input_h5ad"])
    adata = adata[adata.obs["cell_type_major"].astype(str) == cell_type]
    X = np.asarray(adata.obsm[cfg["embedding_key"]], dtype=np.float32)
    X_s = scaler.transform(X)
    classes = list(meta["classes"])
    tmean, tstd = float(meta["temporal_mean"]), float(meta["temporal_std"])
    dev = utils.device()
    logits_all, tpred_all = [], []
    with torch.inference_mode():
        for i in range(0, len(X_s), batch):
            xb = torch.from_numpy(X_s[i:i + batch]).to(dev)
            with torch.autocast("cuda", dtype=torch.float16, enabled=dev.type == "cuda"):
                logits, tpred, _, _ = model(xb, with_method=False)
            logits_all.append(logits.float().cpu())
            tpred_all.append(tpred.float().cpu())
    logits = torch.cat(logits_all)
    tpred = torch.cat(tpred_all).numpy() * tstd + tmean
    probs = torch.softmax(logits, dim=1).numpy()
    pred = probs.argmax(1)
    df = pd.DataFrame({"spot_index": np.arange(len(X)), "spot_id": adata.obs.index})
    df["true_stage"] = adata.obs["stage"].astype(str).values
    df["predicted_stage"] = [classes[p] for p in pred]
    df["classification_confidence"] = probs.max(1)
    df["predicted_temporal_value"] = tpred
    for c, name in enumerate(classes):
        df[f"prob_{name}"] = probs[:, c]
    for col in ("sample_id", "donor_id", "cell_type_major", "cell_type_sub",
                "analysis_condition", "Method", "x", "y"):
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values
    return df


def sample_summary(df: pd.DataFrame, rule: str = "classification") -> pd.DataFrame:
    if rule == "classification":
        label = df["predicted_stage"]
    else:
        idx = np.argmin(np.abs(df["predicted_temporal_value"].to_numpy()[:, None]
                               - TIMES[None, :]), axis=1)
        label = pd.Series(np.asarray(CLASSES)[idx], index=df.index)
    df = df.assign(_label=label)
    g = df.groupby("sample_id", observed=True)
    out = pd.DataFrame({"true_stage": g["true_stage"].first(),
                        "analysis_condition": g["analysis_condition"].first()})
    for c in CLASSES:
        out[f"pct_pred_{c}"] = (
            df.groupby("sample_id", observed=True)["_label"]
            .apply(lambda s: (s == c).mean()))
    out["majority_vote_stage"] = (
        df.groupby("sample_id", observed=True)["_label"]
        .apply(lambda s: s.mode().iloc[0] if len(s.mode()) else None))
    out["mean_predicted_temporal_value"] = g["predicted_temporal_value"].mean()
    out["median_predicted_temporal_value"] = g["predicted_temporal_value"].median()
    return out.reset_index().sort_values("sample_id").set_index("sample_id")
