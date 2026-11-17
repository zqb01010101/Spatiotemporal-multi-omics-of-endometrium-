#!/usr/bin/env python3
"""Print and persist the required final pipeline summary."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def summarize(project_root: Path, input_h5ad: Path, output_root: Path) -> str:
    inspect = load(output_root / "reports/input_h5ad_metrics.json")
    dev = load(output_root / "dev_subset/dev_selection_manifest.json")
    model = load(output_root / "model_assets/model_download_manifest.json")
    mapping = load(output_root / "reports/gene_mapping_report.json")
    validation = load(output_root / "reports/dev_validation_metrics.json")
    full = load(output_root / "tokenized/full_no_assay/manifest.json")
    export = load(output_root / "reports/windows_export_manifest.json")
    archive_gib = export["archive_bytes"] / 2**30
    lines = [
        "# Final Nicheformer Preprocessing Summary",
        "",
        f"- 项目根目录: `{project_root.resolve()}`",
        f"- 输入H5AD绝对路径: `{input_h5ad.resolve()}`",
        f"- 输入H5AD SHA256: `{inspect['input_sha256']}`",
        f"- Python解释器路径和版本: `{sys.executable}` / `{sys.version.split()[0]}`",
        f"- 模型revision: `{model['resolved_revision']}`",
        f"- 模型最大序列长度: **{model['max_sequence_length']}**",
        f"- 模型embedding维度: **{model['hidden_dimension']}**",
        f"- 开发子集shape: **{dev['actual_spots']:,} × {dev['genes']:,}**",
        f"- 开发子集stage/sample/donor覆盖: **{dev['stage_coverage']} / {dev['sample_coverage']} / {dev['donor_coverage']}**",
        f"- 输入基因数: **{mapping['input_gene_count']:,}**",
        f"- 成功映射基因数和比例: **{mapping['mapped_input_gene_count']:,} ({mapping['mapped_input_gene_rate']:.2%})**",
        f"- counts保留率: **{mapping['counts_retention_rate']:.4%}**",
        f"- DBiT nonzero mean输出路径: `{(output_root / 'technology_means/dbit_global_nonzero_mean.npy').resolve()}`",
        f"- 开发集token shape: **{validation['dev_shape']}**",
        f"- 开发集验证是否通过: **{validation['all_passed']}**",
        f"- 完整数据token shape: **{full['shape']}**",
        f"- 完整tokenization是否完成: **{full['complete']}**",
        f"- Windows导出目录: `{(output_root / 'windows_export').resolve()}`",
        f"- Windows导出包大小: **{archive_gib:.3f} GiB**",
        f"- 下一步Windows端命令: `py -c \"from huggingface_hub import snapshot_download; snapshot_download('theislab/Nicheformer', revision='{model['resolved_revision']}', local_dir='model_assets/Nicheformer')\"`",
        f"- 最终状态: **{'SUCCESS' if validation['all_passed'] and full['complete'] else 'FAILURE'}**",
        "",
        "完整复现命令、模型源码审计、映射证据、分块日志和 SHA256 清单均保存在本输出目录。",
    ]
    result = "\n".join(lines) + "\n"
    (output_root / "reports/final_report.md").write_text(result, encoding="utf-8")
    print(result)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    summarize(args.project_root, args.input_h5ad, args.output_root)
