# Computational methods and limits

## Data and quality checks

The reference input is `pro_sec_with_RIF.h5ad`, 146,414 spots × 38,299 genes.
Its `X` is log-normalized expression, while `layers["counts"]` contains raw RNA
counts. Only counts enter mapping, mean estimation and tokenization. Spatial
coordinates are copied for downstream use. Checks include sparse types, integer
and nonnegative counts, zero libraries, unique spot IDs, exact index agreement,
spatial agreement and population distributions.

Required observation columns for the full workflow/export are `spot_id`,
`sample_id`, `donor_id`, `x`, `y`, `stage`, `stage_original`, `stage_numeric`,
`condition`, `cell_type_major`, `cell_type_sub`, `total_counts`, `n_genes`, `pct_mt`
and `pct_ribo`. Mapping prioritizes `gene_id_without_version`, `gene_id`,
`gene_symbol` and `gene_name_original`; `is_mt`/`is_ribo` support mapping QC.

## Sampling

The development subset contains 7,200 spots with NumPy seed 20260806. Joint strata
are `analysis_condition`, `stage`, `sample_id` and `cell_type_major`. Quotas are
approximately proportional, keep at least one spot per joint stratum, cap at the
available size and redistribute the remainder. Selected indices are sorted to
preserve full-data order, saved as `full_data_row_index`, and the H5AD is re-read.
Derived RIF labels require matching `RIF_` markers in both stage and sample ID.

## Fixed gene space and mapping

The pinned official model defines 20,310 genes, mostly Ensembl gene IDs. Mapping
uses exact case-sensitive identifiers and Ensembl release 112 / GRCh38 gene-name
links; there is no fuzzy or case-insensitive gene matching. Only Ensembl ID fields
have version suffixes stripped. Ambiguity and unmapped genes are recorded. If
several input genes match one model gene, counts are summed sparsely and compared
with the counts in the retained input columns. Unmapped counts are explicitly
reported; conservation does not imply retention of all input genes.

## DBiT nonzero mean and ranking

The released code retains the original successful implementation: first map
counts into model gene space; then normalize each spot by **the sum of its retained
mapped counts** to 10,000. For each model gene, average only positive normalized
expression; never-expressed genes receive a mean of 1. The full normal + RIF cohort
contributes to the global exploratory reference.

Tokenization applies the same mapping/library normalization, divides each
positive gene value by this gene-specific mean, applies the official monotonic
`log1p` transform, then sorts descending with NumPy's original `argsort` tie
behavior. At most 1,497 positive gene tokens are retained.

This precise normalization denominator matters: it differs from using the
pre-mapping library total. Both mean calculation and tokenization use the same
denominator, and the aligned 32-spot official tokenizer comparison passes.

Global means use both study groups. For held-out donor/sample evaluations,
estimate the reference from training data alone using
`04_compute_dbit_nonzero_mean.py --sample-whitelist ...` or `--donor-whitelist ...`
in a separate output directory, and use those means for all folds' relevant
tokenization. The default runner reproduces the exploratory full-cohort analysis,
not a leakage-free supervised cross-validation experiment.

## Tokens and model inference

Fixed context order is `[species, technology, modality]`. Main `no_assay` tokens
are `[5, 0, 4]` with mask `[True, False, True]`; development sensitivity tokens are
`[5, 11, 4]` with all context positions visible. Sequence length is 1,500 with
PAD ID 0. These constants correspond to the pinned audited revision only.
All gene IDs are read explicitly from the downloaded vocabulary/model gene order.

Input does not encode sample ID, donor ID, stage, disease label, cell type, QC
values, spatial coordinates, or a fabricated batch label. The shared technology
mean is a gene-expression reference, not per-donor batch correction.

Embedding inference calls the official last-layer
`model.get_embeddings(layer=-1, with_context=False)`. Context tokens participate
in attention; the first three output positions are removed before pooling.
The official code then averages **all remaining positions, including PAD query
positions**. The attention mask controls transformer key padding, not this pooling
denominator. Script `09` intentionally reproduces this behavior without substituting
attention-weighted pooling. Output dimension is 512 and inference uses float32.

## Neighborhood analysis

The notebook used `scanpy.pp.neighbors(use_rep="X_nicheformer", n_neighbors=15,
metric="euclidean")` followed by `scanpy.tl.umap(random_state=0)`. Script `10`
exports these operations with explicit seeds. This is a spot embedding-space kNN
graph. It is not a physical neighborhood graph, an inferred gene regulatory
network, a batch correction, or evidence of causal interactions.

## Computational limits

Full expression is always sparse. The AnnData version used here loads the sparse
counts layer into RAM even with a backed read; “backed” therefore does not mean
fully disk-backed layers. The matrix is still sparse, and computations use
256-spot slices. Sufficient RAM is required. Token arrays use disk-backed NumPy
memmaps. Only the official equivalence reference densifies 32 × 20,310 values.

Floating-point GPU inference and approximate neighbor/UMAP routines can differ
slightly across hardware and runtime versions. Use the supplied original version
lock, exact input/resource hashes and saved inference manifests for comparison.
Synthetic tests and a bounded real-model CPU smoke test validate the packaged
entry points; they do not establish that every historical embedding or figure
was regenerated on the original Windows GPU hardware.
