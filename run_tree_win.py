#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import argparse
import subprocess
from pathlib import Path
import pandas as pd
import shutil


def find_iqtree_launcher(base_dir: Path) -> str:
    """
    优先 iqtree_win/iqtree3.exe；否则依次尝试 PATH 里的 iqtree3、iqtree2、iqtree。
    """
    cand = base_dir / "iqtree_win" / "iqtree3.exe"
    if cand.exists():
        return str(cand)
    for name in ("iqtree3", "iqtree2", "iqtree"):
        p = shutil.which(name)
        if p:
            return p
    # 兜底（系统再去解析）
    return "iqtree2"



def run_iqtree(input_fasta: Path, output_dir: Path, threads: int, launcher: str):
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        launcher,
        "-s", str(input_fasta),
        "-m", "MFP",
        "-bb", "1000",
        "-nt", str(threads),
        "-pre", str(output_dir / "iqtree_result")
    ]
    print("Running IQ-TREE:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_iqtree_batch(report_csv: Path, alignment_dir: Path, output_dir: Path,
                     threads: int, launcher: str):
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(report_csv)
    gene_list = df["Gene"].tolist()

    for gene in gene_list:
        # 与原脚本保持一致：输入是 {gene}.fasta
        fasta_path = alignment_dir / f"{gene}.fasta"
        if not fasta_path.exists():
            print(f"[WARNING] Alignment not found: {fasta_path.name}, skipped.")
            continue

        out_prefix = output_dir / gene.replace("_trimmed", "")
        cmd = [
            launcher,
            "-s", str(fasta_path),
            "-m", "MFP",
            "-bb", "1000",
            "-nt", str(threads),
            "-pre", str(out_prefix)
        ]
        print("Running IQ-TREE:", " ".join(map(str, cmd)))
        subprocess.run(cmd, check=True)


def run_astral(tree_dir, output_dir):
    from collections import defaultdict
    import re
    from pathlib import Path

    tree_dir   = Path(tree_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_tree_file = output_dir / "all_gene_trees.tre"
    mapping_file  = output_dir / "sample_gene_mapping.tsv"
    sample_to_genes = defaultdict(set)

    # —— 1) 列出基因名（从 *.treefile 的前缀推断）——
    gene_names = []
    for tf in tree_dir.glob("*.treefile"):
        base = tf.name[:-9] if tf.name.endswith(".treefile") else tf.stem
        if base.endswith(".iqtree_result"):
            base = base[:-15]
        gene_names.append(base)
    if not gene_names:
        raise RuntimeError(f"No .treefile found in {tree_dir}")

    # IQ-TREE 标签净化规则：
    # - 非 [A-Za-z0-9_\-] 的字符 → '_' （注意：破折号 '-' 会保留）
    # - 末尾可能会多一个 '_'，也可能不会
    def iqtree_sanitize(name: str) -> str:
        s = re.sub(r"[^A-Za-z0-9_\-]+", "_", name)  # 保留 '-'
        s = re.sub(r"_+", "_", s)                   # 多余下划线合并
        return s

    variants = []
    for g in set(gene_names):
        raw = re.escape(g)
        san_keepdash = re.escape(iqtree_sanitize(g))          # 破折号保留
        san_nodash   = re.escape(iqtree_sanitize(g).replace("-", "_"))  # 破折号也转为下划线

        # 三种形态 + 末尾下划线可选
        variants.extend([
            raw + "_?",          # 原始
            san_keepdash + "_?", # 括号转 '_'，保留 '-'
            san_nodash + "_?",   # 括号转 '_'，'-' 也转 '_'
        ])

    variants = sorted(set(variants), key=len, reverse=True)
    gene_pat = "|".join(variants)

    label_regex = re.compile(
        r'([\(:,])'                      # 1: 前导 (, : 等
        r'(?P<sample>[^:,\(\)]+?)'       # 样本 ID
        r'_(?P<gene>(' + gene_pat + r'))'  # 基因名变体
        r'(?=[:\),])',                   # 后面必须是 :,) 或 ,
        flags=re.IGNORECASE
    )

    # —— 2) 合并并清洗基因树 —— 
    with open(all_tree_file, "w", encoding="utf-8", newline="\n") as out_f:
        for tf in sorted(tree_dir.glob("*.treefile")):
            with open(tf, "r", encoding="utf-8") as fh:
                tree_str = fh.read().strip()

            # 统计 sample->gene（便于排查哪些没被识别）
            for m in label_regex.finditer(tree_str):
                sample_to_genes[m.group("sample")].add(m.group("gene"))

            # 去掉 “_基因名” 只保留样本 ID
            cleaned = label_regex.sub(lambda m: m.group(1) + m.group("sample"), tree_str)
            out_f.write(cleaned + "\n")

    # —— 3) 输出映射表，便于核查 —— 
    with open(mapping_file, "w", encoding="utf-8", newline="\n") as mf:
        mf.write("SampleID\tGeneNames\n")
        for sample in sorted(sample_to_genes):
            mf.write(f"{sample}\t{','.join(sorted(sample_to_genes[sample]))}\n")

    # —— 4) 运行 ASTRAL —— 
    jar_path = Path(__file__).parent / "ASTRAL-master" / "astral.5.7.8.jar"
    cmd = ["java", "-jar", str(jar_path), "-i", str(all_tree_file), "-o", str(output_dir / "species_tree.tre")]
    print("Running ASTRAL:", " ".join(cmd))
    subprocess.run(cmd, check=True)


# ---------------- main ----------------
def main():
    parser = argparse.ArgumentParser(description="Phylogenetic tree builder (Windows-compatible).")
    parser.add_argument("-m", "--mode", choices=["concat", "coalescent"], required=True,
                        help="Tree construction mode")
    parser.add_argument("-r", "--report", required=True, help="Filtered CSV report file")
    parser.add_argument("-a", "--alignment_dir", required=True, help="Directory with gene alignments")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_dir = Path(__file__).parent
    iqtree_launcher = find_iqtree_launcher(base_dir)

    if args.mode == "concat":
        # 按你现有流程：先用 combin_genes.py 串联，再跑 IQ-TREE
        concat_fasta = output_dir / "concatenated_alignment.fasta"
        combin_path = base_dir / "combin_genes.py"
        cmd = [
            sys.executable, str(combin_path),
            "-r", args.report,
            "-m", args.alignment_dir,
            "-o", str(concat_fasta)
        ]
        print("Concatenating alignments:", " ".join(cmd))
        subprocess.run(cmd, check=True)

        run_iqtree(concat_fasta, output_dir, args.threads, iqtree_launcher)

    elif args.mode == "coalescent":
        single_tree_dir = output_dir / "gene_trees"
        run_iqtree_batch(Path(args.report), Path(args.alignment_dir), single_tree_dir,
                         args.threads, iqtree_launcher)
        run_astral(single_tree_dir, output_dir)


if __name__ == "__main__":
    main()
