# Code-package validation

Validated on 2026-10-21 with the existing isolated Python 3.10.20 environment,
NumPy 1.26.4, torch 2.4.1+cpu, transformers 4.49.0 and scanpy 1.10.4.
All generated smoke-test data/results were written to a temporary directory
outside this source package. No historical analysis result was overwritten.

| Check | Result |
| --- | --- |
| Original scripts/helpers and original tests | 19 copied files verified byte-for-byte against the RIF source |
| Python syntax | All 24 Python files parsed successfully |
| Synthetic test suite | **12 passed** |
| Runner dry-run | All 13 stage commands generated; no output directory created |
| Runner real-data smoke test | `inspect` completed against the reference H5AD |
| Runner checkpoint recovery | `--resume --stop-after inspect` verified products and skipped recomputation |
| Failed-development gate | Test confirms failed validation blocks full tokenization |
| Changed checkpoint product | Test confirms modified products are rejected |
| Fixed local annotation | Release 112 GTF verified against expected SHA256 |
| Real model CPU inference | Two real development spots produced finite float32 `(2, 512)` embeddings |
| Real inference recovery | Same command with `--resume` completed without recomputing the batch |
| Synthetic inference recovery | Test verifies row order, skipped batches and corrupted-batch rejection |
| Synthetic graph/UMAP | 60 spots × 12 artificial embedding dimensions; k=15, Euclidean, seed=0 |
| Graph output round trip | H5AD re-read successfully; shape `(60, 0)`, embeddings `(60, 12)`, UMAP `(60, 2)`, connectivities `(60, 60)` |

The test suite command, from this package directory:

```bash
python -m pytest -q -p no:cacheprovider tests
```

The real inference smoke command uses `09_infer_embeddings.py` with the original
RIF run's `tokenized/dev_no_assay`, audited local model and model manifest, plus
`--device cpu --batch-size 2 --max-spots 2 --threads 8`; its resume uses the same
arguments and `--resume`.

This validation did not reinstall the complete environment on a clean host,
redownload the model/annotation, rerun full-cohort preprocessing, perform full GPU
inference, or render optional Plotly HTML. The original preprocessing run's
reported validation is documented separately in `EXPECTED_RESULTS.md`.
UMAP smoke tests validate execution and file structure, not replication of a
particular historical figure. Hardware-dependent numerical differences remain
possible; see `METHODS.md`.
