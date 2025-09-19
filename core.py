#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import gzip
import shutil
from Bio import SeqIO
import csv
from pathlib import Path
from refgb_database import load_species_database, find_best_reference
from report_filtered import filter_markers


base_dir = os.path.dirname(os.path.realpath(__file__))

def run_command(cmd):
    print(f"\n>>> Running in {base_dir}:\n  {' '.join(cmd)}\n")
    try:
        subprocess.check_call(cmd, cwd=base_dir)
    except subprocess.CalledProcessError as e:
        print(f"Error: command {' '.join(cmd)} failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(1)

def prepare_reference(args, ref_gb, workdir):
    split_dir = os.path.join(workdir, "split_genes")
    org_seq_dir = os.path.join(split_dir, "org_seq")
    gene_db_dir = os.path.join(workdir, "gene_db")
    os.makedirs(split_dir, exist_ok=True)
    os.makedirs(gene_db_dir, exist_ok=True)

    cmd_split = [
        sys.executable,
        "split_genes.py",
        "-input", ref_gb,
        "-out_dir", split_dir,
        "-soft_boundary", f"{args.soft_boundary[0]},{args.soft_boundary[1]}",
        "-max_seq_length", str(args.max_seq_length),
        "-min_seq_length", str(args.min_seq_length),
        *(["--intron_only"] if args.intron_only else []),
        "-t", str(args.threads),
    ]
    run_command(cmd_split)

    # —— 可选：拆分间隔区 —— #
    if args.enable_intergenic:
        inter_dir = os.path.join(split_dir, "intergenic_seq")
        os.makedirs(inter_dir, exist_ok=True)
        cmd_inter = [
            sys.executable,
            "split_intergenic_2.py",
            "-input", ref_gb,
            "-out_dir", inter_dir,
            "--soft_boundary", str(args.soft_boundary[0]), str(args.soft_boundary[1])
        ]
        if args.alias_mode:
            cmd_inter.append("--enable_alias")
        run_command(cmd_inter)

        # 把所有物种子目录下的 .fasta（间隔区）复制到 org_seq/
        for root, dirs, files in os.walk(inter_dir):
            for fn in files:
                if fn.lower().endswith(".fasta"):
                    src = os.path.join(root, fn)
                    dst = os.path.join(org_seq_dir, fn)
                    shutil.copy2(src, dst)
    # —— 拆分完成 —— #

    try:
        with open(ref_gb, "r") as handle:
            rec = next(SeqIO.parse(handle, "genbank"))
            organism = rec.annotations.get("organism", "").strip()
            if not organism:
                raise ValueError("无法从 GenBank 注释中提取 organism。")
            species_str = organism.replace(" ", "_")
    except Exception as e:
        print(f"Error: 读取 GenBank 文件失败：{e}", file=sys.stderr)
        sys.exit(1)

    illegal_chars = set('#|:/\\')
    for fname in os.listdir(org_seq_dir):
        path = os.path.join(org_seq_dir, fname)
        if not os.path.isfile(path):
            continue
        root, ext = os.path.splitext(path)
        if ext.lower() not in (".fa", ".fas", ".fasta"):
            continue
        if ext.lower() != ".fasta":
            new_path = root + ".fasta"
            os.replace(path, new_path)
        else:
            new_path = path

        new_lines = []
        with open(new_path, "r") as fin:
            for line in fin:
                if line.startswith(">"):
                    header_rest = line[1:].rstrip()
                    sanitized = "".join("_" if c in illegal_chars else c for c in header_rest)
                    new_header = f">Species_{species_str}_Repository_RefSeq {sanitized}\n"
                    new_lines.append(new_header)
                else:
                    new_lines.append(line)
        with open(new_path, "w") as fout:
            fout.writelines(new_lines)

    cmd_db = [
        sys.executable,
        "build_database.py",
        "-i", org_seq_dir,
        "-o", gene_db_dir
    ]
    run_command(cmd_db)

    return org_seq_dir, gene_db_dir

def decompress_fastq(original_dir, sample_name, sample_out_dir):
    gz_suffixes = ['.fq.gz', '.fastq.gz']
    for suf in gz_suffixes:
        p1 = os.path.join(original_dir, f"{sample_name}_1{suf}")
        p2 = os.path.join(original_dir, f"{sample_name}_2{suf}")
        if os.path.isfile(p1) and os.path.isfile(p2):
            temp_dir = os.path.join(sample_out_dir, "raw_reads")
            if os.path.isdir(temp_dir):
                shutil.rmtree(temp_dir)
            os.makedirs(temp_dir, exist_ok=True)

            for gz_path in (p1, p2):
                base = os.path.basename(gz_path)
                out_name = base[:-3]
                out_path = os.path.join(temp_dir, out_name)
                print(f"解压 {gz_path} → {out_path}")
                with gzip.open(gz_path, 'rt') as src, open(out_path, 'w') as dst:
                    shutil.copyfileobj(src, dst)
            return temp_dir
    return original_dir


def process_sample(args, sample_name, read_paths, ref_dir, workdir):
    sample_out = os.path.join(workdir, sample_name)
    raw_dir = os.path.join(sample_out, "raw_reads")
    
    # 增加删除后的存在检查
    if os.path.exists(raw_dir):
        print(f"Removing existing: {raw_dir}")
        if os.path.isfile(raw_dir):
            os.remove(raw_dir)
        else:
            shutil.rmtree(raw_dir, ignore_errors=True)  # 强制删除
        
        # 等待确保删除完成
        import time
        retry = 0
        while os.path.exists(raw_dir) and retry < 5:
            time.sleep(0.2)
            retry += 1
        if os.path.exists(raw_dir):
            raise RuntimeError(f"Failed to remove {raw_dir}")

    os.makedirs(raw_dir, exist_ok=True)  # 安全创建目录


    # 把原始文件复制/解压到 raw_reads，并统一命名
    for idx, src in enumerate(read_paths, 1):
        dst = os.path.join(raw_dir, f"{sample_name}_{idx}.fastq")
        if src.endswith('.gz'):       # 解压
            with gzip.open(src, 'rt') as fin, open(dst, 'w') as fout:
                shutil.copyfileobj(fin, fout)
        else:                         # 直接复制
            shutil.copy2(src, dst)

    # 交给 per_gene_filter
    filtered_dir = os.path.join(sample_out, "filtered")
    cmd = [
        sys.executable, "per_gene_filter.py",
        "--sample-dir", raw_dir,
        "--gene-fasta-dir", ref_dir,
        "--out-dir", filtered_dir,
        "--threads", str(args.threads),
        *(["--log-file", args.log_file] if args.log_file else []),
        "--min-depth", str(args.min_depth),
        "--max-depth", str(args.max_depth),
        "--max-size", str(args.max_size),
        * (["--keep-temporaries"] if args.keep_temporaries else []),
        "-kf", str(args.kmer_size),
        "-p", str(args.processes_assembler),   # 或者专门改个名称，避免与 assembler 冲突
    ]
    run_command(cmd)


def run_assembly(args, ref_dir, workdir):
    cmd_assemble = [
        sys.executable,
        "main_assembler.py",
        "-r", ref_dir,
        "-o", workdir,
        "-ka", str(args.ka),
        "-k_max", str(args.k_max),
        "-k_min", str(args.k_min),
        "-limit_count", str(args.limit_count),
        "-iteration", str(args.iteration),
        "-sb", str(args.soft_boundary[0]), str(args.soft_boundary[1]),
        "-p", str(args.processes_assembler)
    ]
    run_command(cmd_assemble)

def run_cat_mafft(result_dir):
    """Step 4: 对所有样本的基因片段做合并 + MAFFT 比对"""
    cmd = [
        sys.executable,
        "cat_mafft.py",
        "--result_dir", result_dir
    ]
    run_command(cmd)

def run_trimal(input_dir, output_dir, mode="automated1"):
    """Step 5: 对每个 *_aligned.fasta 执行 trimAl 切齐"""
    cmd = [
        sys.executable,
        "run_trimal.py",
        "-i", input_dir,
        "-o", output_dir,
        "-m", mode
    ]
    run_command(cmd)

def run_barcoding_report(args, alignment_dir, report_dir):
    """Step 5: 对每个 *_aligned.fasta 生成变异位点报告"""
    # barcoding_report.py 会自行创建/覆盖 output 目录
    cmd = [
        sys.executable,
        "barcoding_report.py",
        "-input", alignment_dir,
        "-pattern", args.pattern,
        "-output", report_dir
    ]
    run_command(cmd)

def run_filter_summary_table(args, workdir):

    report_dir = Path(workdir) / "variant_reports"
    raw_summary = report_dir / "summary_table.csv"
    filtered_summary = report_dir / "summary_table_filtered.csv"

    if not raw_summary.exists():
        print(f"[ERROR] 缺失 summary_table.csv：{raw_summary}")
        return

    filter_markers(
        str(raw_summary), str(filtered_summary), 
        sample_size=args.sample_size,
        seq_pct=args.seq_pct,
        aln_len_min=args.aln_len_min,
        gap_pct_max=args.gap_pct_max,
        var_sites=(args.var_sites_min, args.var_sites_max),
        var_pct=(args.var_pct_min, args.var_pct_max),
        pi_pct_min=args.pi_pct_min
    )
    print(f"[INFO] 已生成过滤结果：{filtered_summary}")


def run_tree_building(args, workdir):

    # ① 从 args 拿到 mode 和 threads
    mode = args.tree_mode
    threads = args.threads

    tree_script   = os.path.join(os.path.dirname(__file__), "run_tree.py")
    alignment_dir = Path(workdir) / "trimmed"
    report_csv    = Path(args.report) if args.report \
                    else Path(workdir) / "variant_reports" / "summary_table_filtered.csv"
    output_dir    = Path(workdir) / "tree_results"

    if not report_csv.exists():
        print(f"[ERROR] 缺失过滤结果 summary_table_filtered.csv：{report_csv}")
        return

    # ② 在 cmd 中使用上面定义的 mode 和 threads
    cmd = [
        sys.executable,
        tree_script,
        "-m", mode,
        "-r", str(report_csv),
        "-a", str(alignment_dir),
        "-o", str(output_dir),
        "-t", str(threads)
    ]
    print(f"[INFO] 正在运行 run_tree.py 构建系统发育树...\n{' '.join(cmd)}")
    subprocess.run(cmd, check=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", help="输入 GenBank 文件路径；若未指定，将通过 --species-name 从数据库中检索")
    parser.add_argument("--file-list", required=True, help="逗号分隔样品清单 CSV 路径，每行格式为：SampleName,Read1[,Read2]")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--threads", type=int, default=4, help="线程数")
    parser.add_argument("--species-name", help="待匹配参考的物种名称（用于数据库自动检索）")


    #split_genes.py
    parser.add_argument("-soft_boundary", nargs=2, type=int, default=[200,200], metavar=("LEFT","RIGHT"), help="基因区域软边界（格式：左,右）")
    parser.add_argument("-max_seq_length", type=int, default=5000, help="最大基因序列长度")
    parser.add_argument("-min_seq_length", type=int, default=200, help="最小基因序列长度")
    parser.add_argument("--intron_only", action="store_true", help="是否提取 intron 区段")

    #split_intergenic_2.py
    parser.add_argument("--enable_intergenic", action="store_true", help="是否提取经典间隔区区域")
    parser.add_argument("--alias_mode", action="store_true", help="是否启用 split_intergenic_2.py 的模糊匹配模式")

    #per_gene_filter.py
    parser.add_argument('--log-file', default=None, help='Log file')
    parser.add_argument('--min-depth', default=50, help='Min allowed coverage', type=int)
    parser.add_argument('--max-depth', default=768, help='Max allowed coverage', type=int)
    parser.add_argument('--max-size', default=6, help='Max allowed size in million bases', type=int)
    parser.add_argument('--keep-temporaries', action='store_true', help='Keep temporary files')
    parser.add_argument('-kf', '--kmer-size', default=31, help='K-mer size', type=int)
    parser.add_argument('-p', '--processes', default=1, help='Number of parallel processes', type=int)


    #main_assembler.py
    parser.add_argument('-ka', metavar='<int>', type=int, help='''kmer of assemble''',  default=39)
    parser.add_argument('-k_max', metavar='<int>', type=int, help='''max kmer of assemble''',  default=39)
    parser.add_argument('-k_min', metavar='<int>', type=int, help='''max kmer of assemble''',  default=21)
    parser.add_argument('-limit_count', metavar='<int>', type=int, help='''limit of kmer count''', required=False, default=2)
    parser.add_argument('-iteration', metavar='<int>', type=int, help='''iteration''', required=False, default=8192)
    parser.add_argument("--processes-assembler", type=int, default=1, help="用于组装步骤的线程数")

    parser.add_argument("--trim-mode", type=str, default="automated1", help="trimAl 模式参数（默认为 automated1）")


    #barcoding_report.py
    parser.add_argument("--pattern", type=str, default="10", help="统计中基因出现在样本中的最小样本数")
    parser.add_argument("--skip-report", action="store_true", help="跳过 barcoding_report 统计")

    #report_filtered.py
    parser.add_argument("--filtered", action="store_true", help="Run report_filtered before tree building")
    parser.add_argument("--sample_size", type=int, default=None, help="Total sample count (optional)")
    parser.add_argument("--seq_pct", type=float, default=0.7, help="Minimum fraction of samples per marker [0–1]")
    parser.add_argument("--aln_len_min", type=int, default=500, help="Minimum alignment length")
    parser.add_argument("--gap_pct_max", type=float, default=60.0, help="Maximum gap percentage")
    parser.add_argument("--var_sites_min", type=int, default=0, help="Minimum variable sites")
    parser.add_argument("--var_sites_max", type=int, default=200, help="Maximum variable sites")
    parser.add_argument("--var_pct_min", type=float, default=0.0, help="Minimum variable sites percentage")
    parser.add_argument("--var_pct_max", type=float, default=10.0, help="Maximum variable sites percentage")
    parser.add_argument("--pi_pct_min", type=float, default=1.0, help="Minimum parsimony informative site %")

    #run_tree.py
    parser.add_argument("--tree", action="store_true", help="Run phylogenetic tree building after filtering")
    parser.add_argument("--tree-mode", choices=["concat", "coalescent"], default="concat", help="系统发育树构建模式")
    parser.add_argument("--report", default=None, help="仅在单独运行建树时指定过滤后的报告文件路径；" "未提供则使用默认 variant_reports/summary_table_filtered.csv")

    args = parser.parse_args()

    def parse_csv(csv_path):
        sample_dict = {}
        with open(csv_path, newline='') as csvfile:
            reader = csv.reader(csvfile)
            for ln, row in enumerate(reader, 1):
                if not row or row[0].startswith('#'):
                    continue
                if len(row) < 2:
                    raise ValueError(f"CSV 第 {ln} 行少于 2 列")
                sample, *reads = row
                sample_dict[sample] = reads
        return sample_dict

    # 若未指定参考但提供了物种名，则从数据库选择
    if not args.ref:
        if not args.species_name:
            raise ValueError("未提供 --ref，也未指定 --species-name，无法确定参考文件。")
        db_path = os.path.join(base_dir, "info_list_cp.tsv")  # 数据库描述文件
        db = load_species_database(db_path)
        ref_gb = find_best_reference(args.species_name, db)
        print(f"自动选取的参考 GenBank 文件：{ref_gb}")
    else:
        ref_gb = os.path.abspath(args.ref)

    workdir = os.path.abspath(args.outdir)
    os.makedirs(workdir, exist_ok=True)

    print("Step 1: Splitting reference ...")
    split_ref_dir, _ = prepare_reference(args, ref_gb, workdir)

    sample_table = parse_csv(os.path.abspath(args.file_list))

    print("Step 2: Filtering samples ...")
    for name, reads in sample_table.items():
        process_sample(args, name, reads, split_ref_dir, workdir)

    print("Step 3: Assembling ...")
    for name in sample_table:
        sample_out_dir = os.path.join(workdir, name)
        run_assembly(args, split_ref_dir, sample_out_dir)

    
    print("Step 4: Merging per-gene sequences & running MAFFT ...")
    run_cat_mafft(workdir)                      # workdir == args.outdir

    print("Step 5: Trimming alignments with trimAl ...")
    run_trimal(os.path.join(workdir, "alignments"),
               os.path.join(workdir, "trimmed"), args.trim_mode)

    print("Step 6: Calculating variable-site statistics ...")
    run_barcoding_report(args,
                         os.path.join(workdir, "trimmed"),
                         os.path.join(workdir, "variant_reports"))

    if args.filtered:
        print("Step 7: Filtering summary table ...")
        run_filter_summary_table(args, workdir)
    else:
        print("[跳过] Step 7: 不执行 summary_table_filtered.csv 过滤")

    if args.tree:
        print("Step 8: Building phylogenetic tree ...")
        run_tree_building(args, workdir)
    else:
        print("[跳过] Step 8: 不执行系统发育树构建")


    print(f"All finished.")


if __name__ == "__main__":
    main()
