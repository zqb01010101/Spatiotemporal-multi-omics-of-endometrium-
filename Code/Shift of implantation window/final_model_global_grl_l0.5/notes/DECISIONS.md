# Decision log - Global GRL (lambda_adv=0.5), pro-inclusive

## Architecture

512D Nicheformer -> StandardScaler(train) -> LayerNorm(512)
-> Linear(512->256) -> GELU -> Dropout(0.2)
-> ResidualBlock(256) -> Linear(256->128) -> GELU -> Dropout(0.2)
-> z_stage (128D)
   ├── Stage head (128 -> 7)
   ├── Temporal head (128 -> 1), pro_stage target = -1
   └── GRL(lambda=0.5, DANN) -> Method discriminator (128 -> 64 -> 2)

L = CE(stage) + 0.3 * Huber(time) + CE(method | dual stages, via GRL)

## Key choices

- Method branch: Global GRL (no stage one-hot conditioning).
- lambda_adv = 0.5.
- lambda_time = 0.3 (default).
- pro_stage included, temporal target = -1.
- Per cell type: stroma / epithelial.
- Method loss only on dual-method stages: pro_stage, sec_PO3, sec_PO7, sec_PO11.

## Reproducibility

- random seed 42, sample-level grouped holdout (4 repeats, full coverage),
  scaler/class-weights/temporal stats fitted on training splits only.
