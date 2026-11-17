"""Inference with the final Global GRL model + reference-format sample summaries."""

from __future__ import annotations

import argparse

import common


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-type", required=True, choices=["stroma", "epithelial"])
    args = p.parse_args()
    cfg = dict(common.CFG)
    cfg["cell_type_filter"] = args.cell_type
    model, scaler, meta = common.load_final_model(args.cell_type)
    df = common.predict_ct(cfg, model, scaler, meta, args.cell_type)
    out = common.HERE / "results" / args.cell_type
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "all_spot_predictions.csv.gz", index=False, compression="gzip")
    common.sample_summary(df, "classification").to_csv(
        out / "all_sample_predictions.csv", index=True)
    common.sample_summary(df, "temporal_nearest").to_csv(
        out / "all_sample_predictions_temporal_nearest.csv", index=True)
    print(f"[{args.cell_type}] saved {len(df)} spots and sample summaries under {out}")


if __name__ == "__main__":
    main()
