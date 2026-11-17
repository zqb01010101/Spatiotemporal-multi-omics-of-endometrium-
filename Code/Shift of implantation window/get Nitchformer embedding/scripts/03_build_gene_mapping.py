#!/usr/bin/env python3
"""Build an explicit, auditable input-to-Nicheformer gene mapping."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from common import atomic_json, describe, iter_csr_chunks, setup_logging, sha256_file
from gene_utils import mapping_vector, remap_csr


INPUT_PRIORITY = ["gene_id_without_version", "gene_id", "gene_symbol", "gene_name_original"]
MODEL_FIELDS = ["gene_id_without_version", "gene_id", "gene_symbol", "gene_name_original", "gene_name", "gene_id"]


def clean_value(value: Any, strip_version: bool) -> str | None:
    """Normalize missingness and optionally strip an Ensembl version suffix."""
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return None
    return re.sub(r"\.\d+$", "", text) if strip_version else text


def vocab_mapping(raw: Any) -> dict[str, int]:
    """Normalize supported vocabulary JSON orientations to token->ID."""
    if not isinstance(raw, dict):
        raise TypeError("vocab.json must be a JSON object")
    if all(isinstance(value, int) for value in raw.values()):
        return {str(key): int(value) for key, value in raw.items()}
    if all(str(key).isdigit() and isinstance(value, str) for key, value in raw.items()):
        return {value: int(key) for key, value in raw.items()}
    raise ValueError("Unsupported vocab.json orientation")


def build_model_lookup(model_var: pd.DataFrame, model_names: list[str]) -> dict[str, list[int]]:
    """Build exact model identifier lookup, retaining ambiguity."""
    values: dict[str, set[int]] = {}
    fields = [field for field in MODEL_FIELDS if field in model_var.columns]
    for idx, name in enumerate(model_names):
        candidates = [clean_value(name, name.startswith("ENSG"))]
        for field in fields:
            candidates.append(clean_value(model_var.iloc[idx][field], field in {"gene_id", "gene_id_without_version"}))
        for candidate in candidates:
            if candidate:
                values.setdefault(candidate, set()).add(idx)
    return {key: sorted(indices) for key, indices in values.items()}



def augment_lookup_from_ensembl_gtf(lookup: dict[str, list[int]], annotation_gtf: Path, model_names: list[str]) -> dict[str, list[int]]:
    """Add exact Ensembl gene_name-to-model-ID links from a fixed GTF."""
    values = {key: set(indices) for key, indices in lookup.items()}
    model_index = {name.split(".")[0]: idx for idx, name in enumerate(model_names)}
    with gzip.open(annotation_gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
            gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not gene_id_match or not gene_name_match:
                continue
            gene_id = gene_id_match.group(1).split(".")[0]
            if gene_id in model_index:
                values.setdefault(gene_name_match.group(1), set()).add(model_index[gene_id])
    return {key: sorted(indices) for key, indices in values.items()}


def token_ids_for_model_genes(model_names: list[str], vocab: dict[str, int]) -> np.ndarray:
    """Resolve exact model gene token IDs and reject guessed offsets."""
    token_ids = np.array([vocab.get(name, -1) for name in model_names], dtype=np.int64)
    if np.any(token_ids < 0):
        missing = np.flatnonzero(token_ids < 0)
        raise ValueError(
            f"vocab.json lacks exact entries for {len(missing)} model genes; first indices {missing[:10].tolist()}. "
            "Refusing to guess a token offset."
        )
    if len(np.unique(token_ids)) != len(token_ids):
        raise ValueError("Model gene token IDs are not unique")
    return token_ids


def build_mapping(
    input_h5ad: Path, model_dir: Path, annotation_gtf: Path | None,
    output_root: Path, chunk_size: int, min_retention: float,
) -> dict[str, Any]:
    """Map genes, profile retained counts, and enforce a retention gate."""
    logger = setup_logging(output_root / "logs/03_build_gene_mapping.log")
    input_digest = sha256_file(input_h5ad)
    model_h5ad = ad.read_h5ad(model_dir / "model.h5ad", backed="r")
    try:
        model_names = model_h5ad.var_names.astype(str).tolist()
        model_var = model_h5ad.var.copy()
    finally:
        model_h5ad.file.close()
    raw_vocab = json.loads((model_dir / "vocab.json").read_text(encoding="utf-8"))
    vocab = vocab_mapping(raw_vocab)
    model_token_ids = token_ids_for_model_genes(model_names, vocab)
    lookup = build_model_lookup(model_var, model_names)
    annotation_sha256 = None
    if annotation_gtf is not None:
        if not annotation_gtf.is_file():
            raise FileNotFoundError(annotation_gtf)
        lookup = augment_lookup_from_ensembl_gtf(lookup, annotation_gtf, model_names)
        annotation_sha256 = sha256_file(annotation_gtf)

    adata = ad.read_h5ad(input_h5ad, backed="r")
    try:
        input_var = adata.var.copy()
        input_names = adata.var_names.astype(str).tolist()
        records: list[dict[str, Any]] = []
        for idx, input_name in enumerate(input_names):
            matches: list[int] = []
            method = ""
            matched_value = ""
            ambiguous = False
            for field in INPUT_PRIORITY:
                if field not in input_var.columns:
                    continue
                value = clean_value(input_var.iloc[idx][field], field in {"gene_id", "gene_id_without_version"})
                if not value:
                    continue
                candidates = lookup.get(value, [])
                if len(candidates) == 1:
                    matches = candidates
                    method = field
                    matched_value = value
                    break
                if len(candidates) > 1:
                    ambiguous = True
            mapped = len(matches) == 1
            model_index = matches[0] if mapped else -1
            records.append({
                "input_index": idx,
                "input_var_name": input_name,
                "gene_id_without_version": clean_value(input_var.iloc[idx].get("gene_id_without_version"), True) or "",
                "gene_id": clean_value(input_var.iloc[idx].get("gene_id"), True) or "",
                "gene_symbol": clean_value(input_var.iloc[idx].get("gene_symbol"), False) or "",
                "gene_name_original": clean_value(input_var.iloc[idx].get("gene_name_original"), False) or "",
                "mapped": mapped,
                "mapping_method": method,
                "matched_value": matched_value,
                "model_index": model_index,
                "model_gene": model_names[model_index] if mapped else "",
                "gene_token_id": int(model_token_ids[model_index]) if mapped else -1,
                "unmapped_reason": "ambiguous_model_identifier" if ambiguous and not mapped else ("no_exact_match" if not mapped else ""),
                "is_mt": bool(input_var.iloc[idx].get("is_mt", False)),
                "is_ribo": bool(input_var.iloc[idx].get("is_ribo", False)),
            })
        mapping = pd.DataFrame(records)
        vector = mapping_vector(mapping)
        model_present = np.zeros(len(model_names), dtype=bool)
        model_present[vector[vector >= 0]] = True
        expressed = np.zeros(adata.n_vars, dtype=bool)
        input_gene_counts = np.zeros(adata.n_vars, dtype=np.float64)
        mapped_nnz = np.zeros(adata.n_obs, dtype=np.int32)
        total_counts = 0.0
        retained_counts = 0.0
        conservation_failures = 0
        for start, end, chunk in iter_csr_chunks(adata.layers["counts"], adata.n_obs, chunk_size):
            expressed[np.unique(chunk.indices)] = True
            input_gene_counts += np.asarray(chunk.sum(axis=0)).ravel()
            mapped_chunk = remap_csr(chunk, vector, len(model_names))
            original_sum = float(chunk.sum(dtype=np.float64))
            retained_sum = float(mapped_chunk.sum(dtype=np.float64))
            expected = float(chunk[:, vector >= 0].sum(dtype=np.float64))
            conservation_failures += int(not np.isclose(retained_sum, expected, rtol=0, atol=1e-6))
            total_counts += original_sum
            retained_counts += retained_sum
            mapped_nnz[start:end] = np.diff(mapped_chunk.indptr)
            logger.info("Profiled mapping rows [%d, %d)", start, end)
    finally:
        adata.file.close()

    mapping_dir = output_root / "mapping"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    mapping_path = mapping_dir / "input_to_nicheformer_gene.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)
    duplicates = mapping[mapping["mapped"]].groupby(["model_index", "model_gene"], as_index=False).filter(lambda x: len(x) > 1)
    duplicates.to_csv(mapping_dir / "duplicate_mapping_groups.tsv", sep="\t", index=False)
    unmapped = mapping[~mapping["mapped"]].copy()
    unmapped["total_input_counts"] = input_gene_counts[unmapped["input_index"].to_numpy()]
    unmapped.sort_values("total_input_counts", ascending=False).to_csv(mapping_dir / "unmapped_input_genes.tsv", sep="\t", index=False)
    pd.DataFrame({"model_index": np.flatnonzero(~model_present), "model_gene": np.asarray(model_names)[~model_present], "gene_token_id": model_token_ids[~model_present]}).to_csv(
        mapping_dir / "model_genes_not_present.tsv", sep="\t", index=False
    )
    var_mapping = mapping.copy()
    var_mapping.to_parquet(output_root / "metadata/var_mapping.parquet", index=False)
    mapping_digest = sha256_file(mapping_path)
    mapped_count = int(mapping["mapped"].sum())
    expressed_count = int(expressed.sum())
    mapped_expressed = int((mapping["mapped"].to_numpy() & expressed).sum())
    retention = retained_counts / total_counts if total_counts else 0.0
    metrics: dict[str, Any] = {
        "input_h5ad": str(input_h5ad.resolve()),
        "input_sha256": input_digest,
        "mapping_sha256": mapping_digest,
        "annotation_gtf": str(annotation_gtf.resolve()) if annotation_gtf else None,
        "annotation_release": "Ensembl 112 / GRCh38" if annotation_gtf else None,
        "annotation_sha256": annotation_sha256,
        "input_gene_count": len(mapping),
        "mapped_input_gene_count": mapped_count,
        "mapped_input_gene_rate": mapped_count / len(mapping),
        "expressed_input_gene_count": expressed_count,
        "mapped_expressed_input_gene_count": mapped_expressed,
        "mapped_expressed_input_gene_rate": mapped_expressed / max(expressed_count, 1),
        "model_gene_count": len(model_names),
        "covered_model_gene_count": int(model_present.sum()),
        "model_gene_coverage_rate": float(model_present.mean()),
        "total_input_counts": total_counts,
        "retained_counts": retained_counts,
        "counts_retention_rate": retention,
        "mapped_nonzero_genes_per_spot": describe(mapped_nnz),
        "many_to_one_model_gene_groups": int(duplicates["model_index"].nunique()) if len(duplicates) else 0,
        "count_conservation_chunk_failures": conservation_failures,
        "mitochondrial": {
            "total": int(mapping["is_mt"].sum()),
            "mapped": int((mapping["is_mt"] & mapping["mapped"]).sum()),
        },
        "ribosomal": {
            "total": int(mapping["is_ribo"].sum()),
            "mapped": int((mapping["is_ribo"] & mapping["mapped"]).sum()),
        },
        "minimum_required_retention": min_retention,
        "retention_gate_passed": retention >= min_retention,
    }
    atomic_json(output_root / "reports/gene_mapping_report.json", metrics)
    top_unmapped = unmapped.head(20)[["input_var_name", "gene_symbol", "total_input_counts"]]
    report = f"""# Gene Mapping Report

- Input genes: **{len(mapping):,}**
- Mapped input genes: **{mapped_count:,} ({metrics['mapped_input_gene_rate']:.2%})**
- Expressed input genes mapped: **{mapped_expressed:,}/{expressed_count:,} ({metrics['mapped_expressed_input_gene_rate']:.2%})**
- Model genes covered: **{metrics['covered_model_gene_count']:,}/{len(model_names):,} ({metrics['model_gene_coverage_rate']:.2%})**
- Raw-count retention: **{retention:.4%}**
- Many-to-one model genes: **{metrics['many_to_one_model_gene_groups']:,}**
- Count-conservation chunk failures: **{conservation_failures}**
- Mitochondrial mapped/total: **{metrics['mitochondrial']['mapped']}/{metrics['mitochondrial']['total']}**
- Ribosomal mapped/total: **{metrics['ribosomal']['mapped']}/{metrics['ribosomal']['total']}**
- Mapping SHA256: `{mapping_digest}`
- Annotation: `Ensembl release 112 / GRCh38`, SHA256 `{annotation_sha256}`
- Retention gate (≥ {min_retention:.1%}): **{metrics['retention_gate_passed']}**

Mapping used exact identifiers only, in the required priority order. Ensembl version suffixes were removed only from Ensembl ID fields. No case folding, fuzzy matching, or unique-name suffix interpretation was used. Multiple input genes mapping to one model gene are summed by sparse COO canonicalization and checked against the exact retained-column sum for every chunk.

## Highest-Count Unmapped Genes

{top_unmapped.to_markdown(index=False, floatfmt='.0f')}

## Mapped Nonzero Genes per Spot

`{metrics['mapped_nonzero_genes_per_spot']}`
"""
    (output_root / "reports/gene_mapping_report.md").write_text(report, encoding="utf-8")
    if conservation_failures:
        raise RuntimeError(f"Counts conservation failed in {conservation_failures} chunks")
    if retention < min_retention:
        raise RuntimeError(f"Counts retention {retention:.2%} is below required gate {min_retention:.2%}")
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--annotation-gtf", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--min-count-retention", type=float, default=0.5)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_mapping(args.input_h5ad, args.model_dir, args.annotation_gtf, args.output_root, args.chunk_size, args.min_count_retention)
