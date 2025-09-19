# run_trimal.py
import os
import argparse
import subprocess
from pathlib import Path

def run_trimal_on_folder(input_dir: str, output_dir: str, mode: str = "automated1"):
    os.makedirs(output_dir, exist_ok=True)

    for fasta_file in Path(input_dir).glob("*_aligned.fasta"):
        output_file = Path(output_dir) / fasta_file.name.replace("_aligned.fasta", "_trimmed.fasta")
        cmd = [
            "trimal",
            "-in", str(fasta_file),
            "-out", str(output_file),
            f"-{mode}"
        ]
        print(f"Trimming alignment: {fasta_file.name} → {output_file.name}")
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run trimAl on aligned FASTA files.")
    parser.add_argument("-i", "--input_dir", required=True, help="输入文件夹，包含 *_aligned.fasta")
    parser.add_argument("-o", "--output_dir", required=True, help="输出文件夹，将生成 *_trimmed.fasta")
    parser.add_argument("-m", "--mode", default="automated1", help="trimAl 模式（默认为 automated1）")

    args = parser.parse_args()
    run_trimal_on_folder(args.input_dir, args.output_dir, args.mode)
