# Reference preprocessing results

These values are transcribed from the successful 2026-08-08 RIF preprocessing run.
They are comparison targets, not the result of re-running the complete cohort
during code organization.

| Check | Reference result |
| --- | --- |
| Input SHA256 | `2960e65257c74312b196e172e6226c9a632cc613c16ce0b62a1f06254ca2544b` |
| Input shape | 146,414 × 38,299 |
| Sparse counts dtype / nnz | int32 / 304,082,467 |
| Negative / noninteger counts / empty spots | 0 / 0 / 0 |
| Normal / RIF spots | 108,643 / 37,771 |
| Stages / samples / donors | 8 / 24 / 21 |
| Development shape | 7,200 × 38,299 |
| Development normal / RIF | 5,343 / 1,857 |
| Development stage / sample / donor coverage | 8 / 24 / 21 |
| Mapped input genes | 19,000 / 38,299 (49.61%) |
| Model genes covered | 19,000 / 20,310 (93.55%) |
| Counts retention | 94.5135561622% |
| Mapping conservation failures | 0 |
| Expressed model genes in global mean | 18,527 |
| Development arrays (each strategy) | 7,200 × 1,500 |
| Full no-assay arrays | 146,414 × 1,500 |
| Full token chunks at size 256 | 572 |
| Embedding dimension | 512 |
| Official tokenizer equivalence | 32 spots; exact token/mask equality |
| CPU forward validation | 2 × 512 finite embeddings |

Historical full token files, checksums and complete run reports remain in the
original analysis directory; this public code package contains no per-spot data.
Only code can be shared until the study's data release location is supplied.
