"""Reproduce the final Global GRL (lambda_adv=0.5, lambda_time=0.3) model."""

from __future__ import annotations

import argparse

import joblib
import numpy as np
import pandas as pd
import torch

import common
from src import config, utils
from src.evaluate import sample_level, spot_metrics
from src.load_data import load_adata, prepare_normal_data
from src.method_data import add_method_column, audit_method
from src.split_samples import make_splits
from src.stage_mapping import build_targets
from src.train_method_invariant import fit_method_model


def all_data_split(obs_norm, cfg):
    sample_ids = obs_norm["sample_id"].to_numpy()
    all_idx = np.arange(len(obs_norm))
    rng = np.random.default_rng(int(cfg["random_seed"]) + 999)
    samp = np.unique(sample_ids)
    val_samp = rng.choice(samp, size=max(1, int(round(len(samp) * 0.15))), replace=False)
    val_idx = np.flatnonzero(np.isin(sample_ids, val_samp))
    return type("S", (), {"train_idx": all_idx, "val_idx": val_idx,
                         "test_idx": np.array([], dtype=int)})()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-type", required=True, choices=["stroma", "epithelial"])
    p.add_argument("--smoke", action="store_true")
    args = p.parse_args()

    cfg = config.load_config(common.HERE / "config.yaml")
    cfg["_project_root"] = str(common.HERE)
    cfg["cell_type_filter"] = args.cell_type
    cfg["random_seed"] = int(cfg["random_seed"])
    utils.set_seed(cfg["random_seed"])
    if args.smoke:
        cfg["max_epochs"] = 1

    adata = load_adata(cfg)
    add_method_column(adata, cfg)
    audit = audit_method(adata, cfg)
    X, X_norm, obs_norm, _ = prepare_normal_data(adata, cfg)
    y_class, y_time, classes, temporal_source = build_targets(obs_norm, cfg, "v1")
    splits = make_splits(obs_norm, cfg)
    dual = [s for s in audit["dual_method_stages"] if s in classes]
    n_classes = len(classes)
    out_root = common.HERE / ("checkpoints" if args.smoke else "results") / args.cell_type
    logger = utils.setup_logger(f"train_{args.cell_type}", out_root / "train.log")
    logger.info(f"classes={classes}; dual={dual}; repeats={len(splits['splits'])}")

    pooled = []
    for split in splits["splits"]:
        res, *_ = fit_method_model(
            X_norm, y_class, y_time, obs_norm, split, cfg, "v1", out_root, logger,
            f"repeat{split.repeat}", n_classes, classes, temporal_source,
            conditional=False, lambda_max=0.5, dual_stages=dual)
        pooled.append(res["test_predictions"])
    pooled_df = pd.concat(pooled, ignore_index=True).drop_duplicates(
        "spot_index", keep="first")
    pooled_df.to_csv(out_root / "spot_predictions_pooled.csv.gz",
                     index=False, compression="gzip")
    m = spot_metrics(pooled_df, classes, "v1", cfg, out_root, "mlp")
    m.update(sample_level(pooled_df, classes, out_root, "mlp"))
    utils.save_json(m, out_root / "metrics.json")

    final_split = all_data_split(obs_norm, cfg)
    f_res, f_model, f_scaler, f_tstats, f_hist = fit_method_model(
        X_norm, y_class, y_time, obs_norm, final_split, cfg, "v1", out_root, logger,
        "final", n_classes, classes, temporal_source,
        conditional=False, lambda_max=0.5, dual_stages=dual)
    model_dir = common.HERE / "models" / args.cell_type
    model_dir.mkdir(parents=True, exist_ok=True)
    torch.save(f_model.state_dict(), model_dir / "model.pt")
    joblib.dump(f_scaler, model_dir / "scaler.joblib")
    f_hist.to_csv(model_dir / "training_history.csv", index=False)
    utils.save_json({
        "classes": classes, "input_dim": int(X_norm.shape[1]),
        "lambda_adv": 0.5, "conditional": False, "grl_schedule": "dann",
        "dual_method_stages": dual, "method_mapping_raw": cfg.get("method_mapping_raw"),
        "temporal_mean": f_tstats[0], "temporal_std": f_tstats[1],
        "temporal_source": temporal_source, "random_seed": cfg["random_seed"],
        "cell_type": args.cell_type, "pro_temporal_value": -1,
        "lambda_time": cfg["lambda_time"],
    }, model_dir / "metadata.json")
    logger.info("final Global GRL model saved")


if __name__ == "__main__":
    main()
