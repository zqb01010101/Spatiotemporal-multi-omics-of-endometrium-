# -*- coding: utf-8 -*-
"""label_audit/common.py — Shared paths, configuration, logging, grouping, CV, and models

Peak–gene label audit rules:
  * Identify potential labeling errors from out-of-fold predictions only.
  * Put both contexts of each peak–gene pair and its repeat chain in the same CV group.
  * Use batch_id and human confidence only as audit metadata.
  * Reuse the v4.1 models in scripts/02_train.py.
"""

import hashlib
import importlib.util
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

MODULE_DIR = Path(__file__).resolve().parent
MODEL_DIR = MODULE_DIR.parent
APP_DIR = MODEL_DIR.parent
APP_LABELS = APP_DIR / "peak_gene_label_app" / "labels"

CONFIG = json.loads((MODULE_DIR / "config.json").read_text(encoding="utf-8"))

DATA_DIR = MODULE_DIR / CONFIG.get("data_dir", "data")
OUT_DIR = MODULE_DIR / CONFIG.get("outputs_dir", "outputs")
FIG_DIR = MODULE_DIR / CONFIG.get("figures_dir", "figures")
LOG_DIR = OUT_DIR / CONFIG.get("log_dir", "logs")
for d in (DATA_DIR, OUT_DIR, FIG_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

STAGES = ["PO3", "PO5", "PO7", "PO9", "PO11", "PO13"]
FEATURES = [f"peak_{s}" for s in STAGES] + [f"gene_{s}" for s in STAGES]
CLASSES = ["pre_open", "sync_open", "sync_repress", "unrelated"]
HUMAN_SCORE_COLS = ["pre_open", "sync_open", "sync_repress", "unrelated"]
DB_SCORE_COLS = ["score_pre_open", "score_sync_open", "score_sync_repress", "score_unrelated"]


def cfg(key, default=None):
    return CONFIG.get(key, default)


def get_logger(name):
    log_file = LOG_DIR / f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger(f"label_audit.{name}")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(fh)
        logger.addHandler(sh)
    return logger


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def import_train_module():
    """Reuse ScoreMLP, Adam, train, and derive_hard_label from scripts/02_train.py."""
    spec = importlib.util.spec_from_file_location(
        "label_audit_train_mod", MODEL_DIR / "scripts" / "02_train.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_audit_dataset():
    return pd.read_csv(DATA_DIR / "audit_dataset.csv")


def human_hard_label(scores10, conf, min_score=None, min_margin=None, min_conf=None):
    min_score = min_score if min_score is not None else cfg("hard_label_min_score", 7.0)
    min_margin = min_margin if min_margin is not None else cfg("hard_label_min_margin", 2.0)
    min_conf = min_conf if min_conf is not None else cfg("hard_label_min_confidence", 6.0)
    s = np.asarray(scores10, dtype=float)
    order = np.argsort(-s)
    top, second = s[order[0]], s[order[1]]
    if top >= min_score and (top - second) >= min_margin and conf >= min_conf:
        return CLASSES[int(order[0])]
    return "uncertain"


def build_group_mapping(df):
    """group_id joins both contexts of each peak–gene pair and repeat chains using union-find."""
    id_col = "sample_id" if "sample_id" in df.columns else "candidate_id"
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for sid in df[id_col].astype(str):
        parent[sid] = sid
    pair_rep = {}
    for sid, pk, gn in zip(df[id_col].astype(str),
                           df["peak"].astype(str), df["gene"].astype(str)):
        key = pk + "|" + gn
        if key in pair_rep:
            union(sid, pair_rep[key])
        else:
            pair_rep[key] = sid
    for sid, ro in zip(df[id_col].astype(str),
                       df["repeat_of"].fillna("").astype(str)):
        if ro and ro in parent:
            union(sid, ro)

    root2id, groups = {}, []
    for sid in df[id_col].astype(str):
        r = find(sid)
        if r not in root2id:
            root2id[r] = f"g{len(root2id) + 1:06d}"
        groups.append(root2id[r])
    return pd.Series(groups, index=df.index, name="group_id")


def group_stratified_splitter(groups, strat, n_splits, seed):
    """Use StratifiedGroupKFold when available. Otherwise, use StratifiedKFold on groups labeled by majority vote.
    Return (splits, used_method), with splits as a list of (train_idx, test_idx) sample indices."""
    groups = np.asarray(groups)
    strat = np.asarray(strat)
    try:
        from sklearn.model_selection import StratifiedGroupKFold
        sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        splits = list(sgkf.split(np.zeros(len(groups)), strat, groups))
        return splits, "StratifiedGroupKFold"
    except Exception:
        pass
    from sklearn.model_selection import StratifiedKFold
    gdf = pd.DataFrame({"group": groups, "strat": strat})
    gm = gdf.groupby("group")["strat"].agg(
        lambda s: s.value_counts().index[0]).reset_index()
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    splits = []
    for tr_g, te_g in skf.split(gm["group"], gm["strat"]):
        tr_set = set(gm["group"].iloc[tr_g])
        te_set = set(gm["group"].iloc[te_g])
        tr_idx = np.where(np.isin(groups, list(tr_set)))[0]
        te_idx = np.where(np.isin(groups, list(te_set)))[0]
        splits.append((tr_idx, te_idx))
    return splits, "fallback: group-majority StratifiedKFold"


def inner_val_split(groups, strat, frac, seed):
    """Hold out a small stratified set of training groups for MLP early stopping."""
    uniq = np.unique(groups)
    if len(uniq) < 10:
        n = max(2, int(round(1.0 / max(frac, 0.05))))
        splits, _ = group_stratified_splitter(groups, strat, n, seed)
        tr_idx, iv_idx = splits[-1]
        return tr_idx, iv_idx
    n = int(round(1.0 / frac))
    splits, _ = group_stratified_splitter(groups, strat, n, seed)
    tr_idx, iv_idx = splits[-1]
    return tr_idx, iv_idx
