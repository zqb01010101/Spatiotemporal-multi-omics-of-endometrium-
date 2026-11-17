# RIF 着床窗口 Nicheformer 复现代码

本目录整理了 `pro_sec_with_RIF.h5ad` 对应的完整预处理代码、测试和依赖，
并补充了可执行的 embedding 推理入口及原 notebook 的近邻图/UMAP 步骤。
完整参数、Windows 使用方法和各输出的含义见 [英文 README](README.md)。

## 目录

- `scripts/`：原来的 `00`–`08`、公共函数、新总入口、注释下载、embedding 推理、近邻图/UMAP。
- `tests/`：人工稀疏矩阵、映射、抽样、官方排序语义及断点恢复测试。
- `environment/`：Python 3.10 环境、原依赖锁和 GPU/可视化依赖。
- `docs/`：方法说明、参考结果、来源和本次验证记录。
- `provenance/`：原始源码 SHA256、固定模型和注释版本。

## 快速运行

在本目录执行，并将输入 H5AD 放在上一级目录或自行指定路径：

```bash
conda env create -f environment/environment.yml
conda activate rif_nicheformer
python -m pip install torch==2.4.1+cpu --index-url https://download.pytorch.org/whl/cpu
python -m pip install -r environment/requirements-lock.txt
python -m pytest -q tests
python scripts/run_pipeline.py \
  --input-h5ad ../pro_sec_with_RIF.h5ad \
  --output-root outputs/rif_reference
```

新运行应使用空输出目录；中断后用相同参数加 `--resume`。
使用 `--stop-after validate` 可停在开发验证通过之后；
使用 `--dry-run` 可先查看全部命令。

CPU 预处理完成后，可在 GPU 上执行 `09_infer_embeddings.py`。
`10_neighbors_umap.py` 复现 k=15、Euclidean、UMAP seed=0 的 spot 近邻图分析。
这里的 network 指 embedding 空间的 spot 近邻图。

## 论文方法中的关键细节

1. 表达输入严格来自 `layers["counts"]`。
2. 先做精确基因映射，再按保留基因的 counts 总和归一化到 10,000。
3. 当前全局 DBiT nonzero mean 使用全部 normal + RIF spots；严格交叉验证时需单独用训练组计算。
4. sample、donor、stage、疾病标签、细胞类型和空间坐标不编码进 token。
5. 主版本保留 `[human=5, PAD=0, spatial=4]` 上下文，技术位置 mask=False。
6. 官方 `with_context=False` 会在池化前移除三个上下文位置，但仍对剩余 PAD query 位置取平均。

本目录只整理公开代码及汇总说明，原始数据、逐 spot 元数据、权重、token 和旧结果均未复制。
原来的分析工程保留在原处。正式公开时，请补充论文的数据获取地址、引用信息和代码许可证。

新增入口的验证结果见 [VALIDATION.md](docs/VALIDATION.md)。
