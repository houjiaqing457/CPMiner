import os
import subprocess
import argparse
from pathlib import Path
import pandas as pd
import sys
import re  

def run_iqtree(input_fasta, output_dir, threads=4):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "iqtree2",
        "-s", str(input_fasta),
        "-m", "MFP",
        "-bb", "1000",
        "-nt", str(threads),
        "-pre", str(output_dir / "iqtree_result")
    ]
    subprocess.run(cmd, check=True)


def get_shared_prefix(labels):
    """
    获取标签列表的最长公共前缀（以 _ 分割）
    如 ['503_atpI', '503_ccsA'] -> '503'
    """
    split_labels = [label.split("_") for label in labels]
    shared = []
    for parts in zip(*split_labels):
        if all(part == parts[0] for part in parts):
            shared.append(parts[0])
        else:
            break
    return "_".join(shared) if shared else None


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
        s = re.sub(r"_+", "_", s)                   # 合并多余下划线
        return s

    variants = []
    for g in set(gene_names):
        raw = re.escape(g)
        san_keepdash = re.escape(iqtree_sanitize(g))          # 保留 '-'
        san_nodash   = re.escape(iqtree_sanitize(g).replace("-", "_"))  # '-' 转 '_'

        # 三种形态 + 末尾下划线可选
        variants.extend([
            raw + "_?",          
            san_keepdash + "_?", 
            san_nodash + "_?",   
        ])

    variants = sorted(set(variants), key=len, reverse=True)
    gene_pat = "|".join(variants)

    label_regex = re.compile(
        r'([\(:,])'                      # 前导 (, : 等
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

            # 统计 sample->gene（便于排查）
            for m in label_regex.finditer(tree_str):
                sample_to_genes[m.group("sample")].add(m.group("gene"))

            # 去掉 “_基因名”，只保留样本 ID
            cleaned = label_regex.sub(lambda m: m.group(1) + m.group("sample"), tree_str)
            out_f.write(cleaned + "\n")

    # —— 3) 输出映射表 —— 
    with open(mapping_file, "w", encoding="utf-8", newline="\n") as mf:
        mf.write("SampleID\tGeneNames\n")
        for sample in sorted(sample_to_genes):
            mf.write(f"{sample}\t{','.join(sorted(sample_to_genes[sample]))}\n")

    # —— 4) 运行 ASTRAL —— 
    jar_path = Path(__file__).parent / "ASTRAL-master" / "astral.5.7.8.jar"
    cmd = ["java", "-jar", str(jar_path), "-i", str(all_tree_file), "-o", str(output_dir / "species_tree.tre")]
    print("Running ASTRAL:", " ".join(cmd))
    subprocess.run(cmd, check=True)



def run_iqtree_batch(report_csv, alignment_dir, output_dir, threads=4):
    alignment_dir = Path(alignment_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(report_csv)
    gene_list = df['Gene'].tolist()

    for gene in gene_list:
        fasta_path = alignment_dir / f"{gene}.fasta"
        if not fasta_path.exists():
            print(f"[WARNING] Alignment not found: {fasta_path.name}, skipped.")
            continue

        out_prefix = output_dir / gene.replace("_trimmed", "")
        cmd = [
            "iqtree2",
            "-s", str(fasta_path),
            "-m", "MFP",
            "-bb", "1000",
            "-nt", str(threads),
            "-pre", str(out_prefix)
        ]
        subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Phylogenetic tree builder")
    parser.add_argument("-m", "--mode", choices=["concat", "coalescent"], required=True, help="Tree construction mode")
    parser.add_argument("-r", "--report", required=True, help="Filtered CSV report file")
    parser.add_argument("-a", "--alignment_dir", required=True, help="Directory with gene alignments")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    if args.mode == "concat":
        concat_fasta = Path(args.output) / "concatenated_alignment.fasta"

        # 调用 combin_genes.py 串联
        # 调用 combin_genes.py 串联
        combin_path = Path(__file__).parent / "combin_genes.py"
        cmd = [
            sys.executable,
            str(combin_path),
            "-r", args.report,
            "-m", args.alignment_dir,
            "-o", str(concat_fasta)
        ]
        subprocess.run(cmd, check=True)


        # 构建串联树
        run_iqtree(concat_fasta, args.output, args.threads)

    elif args.mode == "coalescent":
        single_tree_dir = Path(args.output) / "gene_trees"
        run_iqtree_batch(args.report, args.alignment_dir, single_tree_dir, args.threads)
        run_astral(single_tree_dir, args.output)


if __name__ == "__main__":
    main()
