#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import subprocess
from pathlib import Path

def find_trimal_launcher(base_dir: Path) -> str:
    """
    优先使用同级目录 trimal_win/trimal.exe；
    否则退回到 PATH 中的 'trimal'（适用于已配置好可执行的环境）。
    """
    win_exe = base_dir / "trimal_win" / "trimal.exe"
    if win_exe.exists():
        return str(win_exe)
    return "trimal"

def run_trimal_on_folder(input_dir: str, output_dir: str, mode: str = "automated1"):
    in_dir = Path(input_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 解析 trimal 可执行文件位置
    trimal_bin = find_trimal_launcher(base_dir=Path(__file__).parent)

    # 遍历 *_aligned.fasta
    for fasta_file in sorted(in_dir.glob("*_aligned.fasta")):
        # 生成 *_trimmed.fasta
        out_name = fasta_file.name.replace("_aligned.fasta", "_trimmed.fasta")
        output_file = out_dir / out_name

        cmd = [
            trimal_bin,
            "-in", str(fasta_file),
            "-out", str(output_file),
            f"-{mode}"
        ]
        print(f"Trimming alignment: {fasta_file.name} → {output_file.name}")
        # 与原版一致：失败时抛异常
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run trimAl on aligned FASTA files.")
    parser.add_argument("-i", "--input_dir", required=True, help="输入文件夹，包含 *_aligned.fasta")
    parser.add_argument("-o", "--output_dir", required=True, help="输出文件夹，将生成 *_trimmed.fasta")
    parser.add_argument("-m", "--mode", default="automated1", help="trimAl 模式（默认为 automated1）")
    args = parser.parse_args()

    run_trimal_on_folder(args.input_dir, args.output_dir, args.mode)
