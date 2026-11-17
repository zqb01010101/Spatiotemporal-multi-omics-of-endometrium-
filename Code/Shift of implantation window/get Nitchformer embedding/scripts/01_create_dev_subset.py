#!/usr/bin/env python3
"""Create a reproducible, joint-stratified development subset."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from common import atomic_json, setup_logging, sha256_file


STRATA = ["analysis_condition", "stage", "sample_id", "cell_type_major"]


def add_analysis_condition(obs: pd.DataFrame) -> pd.DataFrame:
    """Derive auditable normal/RIF labels from concordant stage and sample IDs."""
    if not {"stage", "sample_id"}.issubset(obs.columns):
        raise KeyError("stage and sample_id are required to identify RIF observations")
    stage_rif = obs["stage"].astype(str).str.startswith("RIF_")
    sample_rif = obs["sample_id"].astype(str).str.startswith("RIF_")
    if not np.array_equal(stage_rif.to_numpy(), sample_rif.to_numpy()):
        mismatch = int(np.count_nonzero(stage_rif.to_numpy() != sample_rif.to_numpy()))
        raise ValueError(f"RIF markers disagree between stage and sample_id for {mismatch} rows")
    result = obs.copy()
    result["analysis_condition"] = np.where(stage_rif, "RIF", "normal")
    return result


def allocate_quotas(sizes: np.ndarray, target: int) -> np.ndarray:
    """Allocate exact near-proportional quotas while retaining each stratum."""
    sizes = np.asarray(sizes, dtype=np.int64)
    if target > int(sizes.sum()):
        raise ValueError("Target exceeds available observations")
    if target < len(sizes):
        raise ValueError("Target is too small to retain every joint stratum")
    ideal = sizes * (target / sizes.sum())
    quota = np.minimum(sizes, np.maximum(1, np.floor(ideal).astype(np.int64)))
    while int(quota.sum()) < target:
        eligible = quota < sizes
        priority = np.where(eligible, ideal - quota, -np.inf)
        idx = int(np.argmax(priority))
        quota[idx] += 1
    while int(quota.sum()) > target:
        eligible = quota > 1
        priority = np.where(eligible, quota - ideal, -np.inf)
        idx = int(np.argmax(priority))
        quota[idx] -= 1
    return quota


def distribution(frame: pd.DataFrame, column: str) -> dict[str, int]:
    return {str(k): int(v) for k, v in frame[column].value_counts(dropna=False).items()}


def select_indices(obs: pd.DataFrame, target: int, seed: int) -> tuple[np.ndarray, pd.DataFrame]:
    """Select rows using seeded sampling within joint strata."""
    if "analysis_condition" not in obs.columns:
        obs = add_analysis_condition(obs)
    missing = sorted(set(STRATA) - set(obs.columns))
    if missing:
        raise KeyError(f"Missing stratification columns: {missing}")
    work = obs[STRATA].astype(object).where(obs[STRATA].notna(), "<NA>").astype(str)
    work["_row"] = np.arange(len(work), dtype=np.int64)
    groups = list(work.groupby(STRATA, sort=True, dropna=False))
    sizes = np.array([len(group) for _, group in groups], dtype=np.int64)
    quotas = allocate_quotas(sizes, target)
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    manifest_rows: list[dict[str, Any]] = []
    for ((analysis_condition, stage, sample_id, cell_type), group), quota in zip(groups, quotas):
        rows = group["_row"].to_numpy()
        chosen = rows if quota == len(rows) else rng.choice(rows, size=int(quota), replace=False)
        selected.append(chosen)
        manifest_rows.append({
            "analysis_condition": analysis_condition,
            "stage": stage,
            "sample_id": sample_id,
            "cell_type_major": cell_type,
            "available": int(len(rows)),
            "selected": int(quota),
        })
    indices = np.sort(np.concatenate(selected)).astype(np.int64)
    if len(indices) != target or len(np.unique(indices)) != target:
        raise RuntimeError("Stratified selection did not produce the exact unique target")
    return indices, pd.DataFrame(manifest_rows)


def create_subset(input_h5ad: Path, output_root: Path, target: int, seed: int) -> dict[str, Any]:
    """Write and re-read the development AnnData subset."""
    logger = setup_logging(output_root / "logs/01_create_dev_subset.log")
    digest = sha256_file(input_h5ad)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    try:
        sampling_obs = add_analysis_condition(adata.obs)
        indices, strata = select_indices(sampling_obs, target, seed)
        subset = adata[indices, :].to_memory(copy=True)
    finally:
        adata.file.close()
    subset.obs["full_data_row_index"] = indices
    subset.obs["analysis_condition"] = sampling_obs.iloc[indices]["analysis_condition"].to_numpy()
    subset_path = output_root / "dev_subset/dev_spots.h5ad"
    subset_path.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(subset_path, compression="gzip")
    np.savetxt(output_root / "dev_subset/dev_spot_ids.txt", subset.obs_names.to_numpy(), fmt="%s")

    check = ad.read_h5ad(subset_path, backed="r")
    try:
        if check.shape != subset.shape or "counts" not in check.layers or "spatial" not in check.obsm:
            raise RuntimeError("Re-read validation of dev_spots.h5ad failed")
        if not np.array_equal(check.obs["full_data_row_index"].to_numpy(), indices):
            raise RuntimeError("full_data_row_index changed after serialization")
    finally:
        check.file.close()

    full_obs = ad.read_h5ad(input_h5ad, backed="r")
    try:
        source_obs = add_analysis_condition(full_obs.obs)
    finally:
        full_obs.file.close()
    selected_obs = subset.obs
    columns = ["analysis_condition", "condition", "stage", "sample_id", "donor_id", "cell_type_major"]
    distributions = {
        column: {"full": distribution(source_obs, column), "dev": distribution(selected_obs, column)}
        for column in columns
    }
    manifest: dict[str, Any] = {
        "input_h5ad": str(input_h5ad.resolve()),
        "input_sha256": digest,
        "target_spots": target,
        "actual_spots": int(subset.n_obs),
        "genes": int(subset.n_vars),
        "seed": seed,
        "strata": STRATA,
        "algorithm": "one-per-joint-stratum floor allocation followed by largest-deficit proportional redistribution; seeded sampling within strata; final indices sorted",
        "relative_order_preserved": bool(np.all(np.diff(indices) > 0)),
        "stage_coverage": int(selected_obs["stage"].nunique(dropna=True)),
        "sample_coverage": int(selected_obs["sample_id"].nunique(dropna=True)),
        "donor_coverage": int(selected_obs["donor_id"].nunique(dropna=True)),
        "cell_type_major_coverage": int(selected_obs["cell_type_major"].nunique(dropna=True)),
        "analysis_condition_coverage": int(selected_obs["analysis_condition"].nunique(dropna=True)),
        "analysis_condition_rule": "RIF iff stage and sample_id both start with RIF_; source condition is retained unchanged",
        "distributions": distributions,
        "joint_strata": strata.to_dict(orient="records"),
        "subset_sha256": sha256_file(subset_path),
    }
    atomic_json(output_root / "dev_subset/dev_selection_manifest.json", manifest)

    sections = []
    for column in columns:
        full_counts = pd.Series(distributions[column]["full"], name="full_count")
        dev_counts = pd.Series(distributions[column]["dev"], name="dev_count")
        table = pd.concat([full_counts, dev_counts], axis=1).fillna(0).astype(int)
        table["full_pct"] = table["full_count"] / len(source_obs) * 100
        table["dev_pct"] = table["dev_count"] / len(selected_obs) * 100
        sections.append(f"## {column}\n\n{table.to_markdown(floatfmt='.3f')}")
    distribution_sections = "\n\n".join(sections)
    report = f"""# Development Subset Report

- Target / actual: **{target:,} / {subset.n_obs:,} spots**
- Shape: **{subset.n_obs:,} × {subset.n_vars:,}**
- Random seed: **{seed}**
- Input SHA256: `{digest}`
- Stage / sample / donor / major-cell-type coverage: **{manifest['stage_coverage']} / {manifest['sample_coverage']} / {manifest['donor_coverage']} / {manifest['cell_type_major_coverage']}**
- Relative source order preserved: **{manifest['relative_order_preserved']}**

Selection jointly stratified by `analysis_condition`, `stage`, `sample_id`, and `cell_type_major`. Because source `condition` is uniformly `unspecified`, `analysis_condition` is derived as RIF only where both `stage` and `sample_id` start with `RIF_`; disagreement is a hard error. The source `condition` is retained unchanged. Every joint stratum receives one row, then remaining capacity is assigned approximately in proportion to source abundance; underfilled strata are capped and remaining capacity is redistributed. Rows are sampled without replacement and sorted by original row index. This preserves rare groups, including low-quality spots, without inflating them beyond the stratified allocation.

{distribution_sections}
"""
    (output_root / "reports/dev_subset_report.md").write_text(report, encoding="utf-8")
    logger.info("Wrote %s with shape %s", subset_path, subset.shape)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--target", type=int, default=7200)
    parser.add_argument("--seed", type=int, default=20260806)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    create_subset(args.input_h5ad, args.output_root, args.target, args.seed)
