#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
per_gene_filter.py —— 已修改版，输出逻辑简化为 filtered/{gene}.fq
"""

import os
import sys
import argparse
import subprocess
import shutil
import tempfile

def run_cmd(cmd, cwd=None):
    print(f"\n>>> Running: {' '.join(cmd)}")
    try:
        subprocess.check_call(cmd, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Error: command {' '.join(cmd)} failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(1)

def per_gene_filter(args):
    fastq_ext = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    files = [f for f in os.listdir(args.sample_dir) if f.endswith(fastq_ext)]
    if not files:
        sys.exit(f"Error: {args.sample_dir} 中未找到 FASTQ 文件")

    # 先尝试找 *_1.* 作为 R1；若不存在，则把第 1 个文件当作单端
    r1_name = next((f for f in files if "_1." in f), files[0])
    paired  = "_1." in r1_name          # 只有带 _1.* 才可能是双端

    r1_path = os.path.join(args.sample_dir, r1_name)
    r2_path = None
    if paired:
        r2_candidate = r1_name.replace("_1.", "_2.")
        if r2_candidate in files:
            r2_path = os.path.join(args.sample_dir, r2_candidate)
        else:
            print(f"Warning: 未找到 {r2_candidate}，回退为单端模式。")
            paired = False

    gene_files = [f for f in os.listdir(args.gene_fasta_dir) if f.lower().endswith(".fasta")]
    if not gene_files:
        print(f"Error: 在 {args.gene_fasta_dir} 下未找到任何 .fasta 文件。", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)
    for gene_file in gene_files:
        gene_name = os.path.splitext(gene_file)[0]
        gene_fasta = os.path.join(args.gene_fasta_dir, gene_file)
        print(f"\n=== 现在处理基因 {gene_name} ===")

        tmpdir = tempfile.mkdtemp(prefix=f"tmp_filter_{gene_name}_")
        try:
            ext1 = None
            if r1_path.endswith(".fastq.gz"):
                ext1 = ".fastq.gz"
            elif r1_path.endswith(".fq.gz"):
                ext1 = ".fq.gz"
            elif r1_path.endswith(".fastq"):
                ext1 = ".fastq"
            elif r1_path.endswith(".fq"):
                ext1 = ".fq"
            else:
                print(f"Error: 无法识别 R1 文件后缀：{r1_path}", file=sys.stderr)
                sys.exit(1)

            # -------- 复制 FASTQ 开始 --------
            ext1 = os.path.splitext(r1_path)[1]          # .fastq / .fq / .gz
            tmp_r1 = os.path.join(tmpdir, f"{gene_name}_1{ext1}")
            shutil.copy(r1_path, tmp_r1)

            if paired:
                tmp_r2 = os.path.join(tmpdir, f"{gene_name}_2{ext1}")
                shutil.copy(r2_path, tmp_r2)
            # -------- 复制 FASTQ 结束 --------
            print(f"  已在 {tmpdir} 创建临时输入文件")

            # 调用 main_refilter_new.py
            gene_outdir = os.path.join(tmpdir, "refilter_out")
            os.makedirs(gene_outdir, exist_ok=True)

            cmd = [
                sys.executable, "main_refilter_new.py",
                ("-qd" if paired else "-qs"), tmpdir,
                "-r", args.gene_fasta_dir,
                "-o", gene_outdir,
                "-p", str(args.threads),
                "--min-depth", str(args.min_depth),
                "--max-depth", str(args.max_depth),
                "--max-size", str(args.max_size),
                "-kf", str(args.kmer_size),
            ]

            if args.keep_temporaries:
                cmd.append("--keep-temporaries")


            run_cmd(cmd)

            # 从 gene_outdir 里找出输出结果，移动到目标 out_dir
            output_found = False
            for fname in os.listdir(gene_outdir):
                if fname.lower().endswith((".fastq", ".fq", ".fasta")):
                    full_path = os.path.join(gene_outdir, fname)
                    if os.path.getsize(full_path) > 0:
                        dest_path = os.path.join(args.out_dir, f"{gene_name}.fq")
                        shutil.move(full_path, dest_path)
                        print(f"  【OK】基因 {gene_name} 过滤有结果，输出至 {dest_path}")
                        output_found = True
                        break
            if not output_found:
                print(f"  【Info】基因 {gene_name} 没有通过过滤的 reads。")
        finally:
            shutil.rmtree(tmpdir)

    print("\n=== per-gene 过滤完成。所有结果集中保存在目录：", args.out_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Per-gene 调用 main_refilter_new.py，简化输出结构")
    parser.add_argument("--sample-dir", required=True, help="样本目录")
    parser.add_argument("--gene-fasta-dir", required=True, help="基因FASTA参考目录")
    parser.add_argument("--out-dir", required=True, help="最终输出目录（filtered/）")
    parser.add_argument("-p","--threads", type=int, default=1, help="线程数")
    parser.add_argument("-d","--min-depth", type=int, default=50, help="最低覆盖度")
    parser.add_argument("-D","--max-depth", type=int, default=768, help="最高覆盖度")
    parser.add_argument("--max-size", type=int, default=6, help="最大允许文件大小（MB）")
    parser.add_argument("-kf","--kmer-size", type=int, default=31, help="kmer 大小")
    parser.add_argument("--keep-temporaries", action="store_true", help="是否保留中间文件")

    args = parser.parse_args()
    per_gene_filter(args)

