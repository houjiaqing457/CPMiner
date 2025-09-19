#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from refgb_database import load_species_database, find_best_reference
from report_filtered import filter_markers
from split_genes import Extract_reference
import shutil
from typing import Tuple

# ======== Win/Linux 外部程序适配工具 ========

def _which(cmd: str) -> str:
    """shutil.which 的轻封装：若找到返回绝对路径，否则原样返回（交给系统 PATH）"""
    p = shutil.which(cmd)
    return p if p else cmd

def find_tool(name: str,
              win_subdir: str | None = None,
              env_var: str | None = None,
              prefer_exts: Tuple[str, ...] = (".exe", ".bat")) -> str:
    """
    跨平台查找外部程序。
    - Windows: 优先 <script_dir>/<win_subdir>/<name>.exe|.bat；再看 env_var 指定目录；再走系统 PATH (name.exe / name.bat)。
    - Linux/Mac: 直接用 name（走 PATH）。
    """
    base = Path(__file__).parent
    if os.name == "nt":
        candidates = []
        if win_subdir:
            for ext in prefer_exts:
                candidates.append(base / win_subdir / f"{name}{ext}")
        if env_var and os.environ.get(env_var):
            for ext in prefer_exts:
                candidates.append(Path(os.environ[env_var]) / f"{name}{ext}")
        # PATH fallback
        for ext in prefer_exts:
            candidates.append(f"{name}{ext}")
        for c in candidates:
            c = str(c)
            if Path(c).exists():
                return c
            w = shutil.which(c)
            if w:
                return w
        return name  # 最后兜底
    else:
        return _which(name)

def run_mafft_xplat(input_fa: Path, output_fa: Path, threads: int):
    """
    与原始 Linux 行为一致：mafft --auto --thread N input > output
    - Win：若同级目录有 mafft-win/mafft.bat，用 cmd 重定向；否则走 PATH 里的 mafft(.exe)。
    - *不改参数*。
    """
    output_fa.parent.mkdir(parents=True, exist_ok=True)
    script_dir = Path(__file__).parent
    if os.name == "nt":
        bat = script_dir / "mafft-win" / "mafft.bat"
        if bat.exists():
            # .bat 必须 shell=True 才能使用 ">" 重定向
            cmd = f'"{bat}" --auto --thread {threads} "{input_fa}" > "{output_fa}"'
            subprocess.run(cmd, shell=True, check=True)
            return
        # 没有 bat 就尝试 PATH 里的 mafft(.exe)
        exe = find_tool("mafft")
        with open(output_fa, "w", encoding="utf-8") as fout:
            subprocess.run([exe, "--auto", "--thread", str(threads), str(input_fa)],
                           stdout=fout, stderr=subprocess.DEVNULL, check=True)
    else:
        # Linux/macOS：直接 list 传参 + stdout 重定向句柄亦可
        with open(output_fa, "w", encoding="utf-8") as fout:
            subprocess.run(["mafft", "--auto", "--thread", str(threads), str(input_fa)],
                           stdout=fout, stderr=subprocess.DEVNULL, check=True)

def run_trimal_xplat(aligned_fa: Path, output_fa: Path, mode: str = "automated1"):
    """
    与原始 Linux 行为一致：trimal -in x -out y -automated1
    - Win：优先 trimal_win/trimal.exe；否则走 PATH。
    """
    output_fa.parent.mkdir(parents=True, exist_ok=True)
    trimal = find_tool("trimal", win_subdir="trimal_win")
    subprocess.run([trimal, "-in", str(aligned_fa), "-out", str(output_fa), f"-{mode}"], check=True)

def blast_exe(name: str) -> str:
    """
    给 BLAST+ 可执行找路径（可选但推荐）：
    - Windows：支持 <script_dir>/blast_win/ 、环境变量 BLAST_BIN，或 PATH。
    - *不改变参数*，只是找得到正确的可执行。
    """
    return find_tool(name, win_subdir="blast_win", env_var="BLAST_BIN")
# ============================================


# -----------------------------------------
def safe_mkdir(path: Path):
    if path.exists() and not path.is_dir():
        print(f"[警告] {path} 已存在但不是目录，类型为：{'file' if path.is_file() else 'unknown'}")
        bak = path.with_name(path.name + ".bak")
        idx = 1
        while bak.exists():
            bak = path.with_name(f"{path.name}.bak{idx}")
            idx += 1
        path.rename(bak)
    path.mkdir(parents=True, exist_ok=True)


def match_blast_test(fasta_path, blast_tsv, output_tsv, missing_list):
    # 直接从fasta读取基因长度
    id_to_length = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        id_to_length[record.id] = len(record.seq)
    ref_df = pd.DataFrame({"id": list(id_to_length.keys()), "length_bp": list(id_to_length.values())})

    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_df = pd.read_csv(blast_tsv, sep="\t", names=blast_cols)

    blast_df["qlen"] = blast_df["qseqid"].map(id_to_length)
    blast_df["coverage"] = blast_df["length"] / blast_df["qlen"]
    valid_hits = blast_df[(blast_df["pident"] >= 90) & (blast_df["coverage"] >= 0.8)]
    hit_ids = set(valid_hits["qseqid"])
    ref_df["status"] = ref_df["id"].apply(lambda x: "present" if x in hit_ids else "missing")
    # 生成完整带BLAST过滤信息的输出
    merged_df = ref_df.merge(blast_df[["qseqid", "pident", "coverage"]].drop_duplicates(subset="qseqid"),
                            left_on="id", right_on="qseqid", how="left")

    merged_df = merged_df.drop(columns=["qseqid"])
    merged_df.to_csv(output_tsv, sep="\t", index=False)

    missing_ids = ref_df[ref_df["status"] == "missing"]["id"].tolist()

    with open(missing_list, "w") as f:
        for gene_id in missing_ids:
            f.write(f"{gene_id}\n")

def rename_fasta_header(fasta_path, gene_name):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    for rec in records:
        rec.id = gene_name
        rec.name = gene_name
        rec.description = ""
    SeqIO.write(records, fasta_path, "fasta")
    return fasta_path


def target_gene(reference_fasta, blast_tsv, scaffold_fasta, output_fasta, gene_name):
    """修改后的 target_gene 函数，统一FASTA头格式"""
    id_to_length = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        id_to_length[record.id] = len(record.seq)
    ref_df = pd.DataFrame({"id": list(id_to_length.keys()), "length_bp": list(id_to_length.values())})
    # 更稳健的版本，明确用非正则方式进行子串匹配，忽略大小写
    target_gene_info = ref_df[ref_df['id'].str.contains(gene_name, case=False, regex=False)]
    if target_gene_info.empty:
        print(f"警告: 未找到基因 {gene_name} 的信息（序列 ID 列表：{list(ref_df['id'])}）")
        return


    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_df = pd.read_csv(blast_tsv, sep="\t", names=blast_cols)

    target_blast_hits = blast_df[blast_df['qseqid'].isin(target_gene_info['id'])]
    if target_blast_hits.empty:
        print(f"警告: 未找到基因 {gene_name} 的 BLAST 匹配结果")
        return
        
    top_hits = target_blast_hits.sort_values(by=["qseqid", "bitscore"], ascending=[True, False]).groupby("qseqid").head(1)

    def extract_prefix(contig_id):
        return contig_id.split("_length")[0] if "_length" in contig_id else contig_id

    fasta_records = {record.id: record for record in SeqIO.parse(scaffold_fasta, "fasta")}
    fasta_prefix_map = {}
    for record_id in fasta_records:
        prefix = extract_prefix(record_id)
        fasta_prefix_map[prefix] = record_id

    with open(output_fasta, "w") as out_file:
        for index, row in top_hits.iterrows():
            contig_id = row['sseqid']
            contig_prefix = extract_prefix(contig_id)
            
            # 获取样本ID (output_fasta的父目录的父目录名称)
            sample_id = Path(output_fasta).parent.name
            
            # 尝试直接匹配完整ID
            if contig_id in fasta_records:
                fasta_contig_id = contig_id
            # 尝试匹配前缀
            elif contig_prefix in fasta_prefix_map:
                fasta_contig_id = fasta_prefix_map[contig_prefix]
            else:
                print(f"警告: 无法找到 contig ID {contig_id} 或前缀 {contig_prefix}")
                continue
                
            record = fasta_records[fasta_contig_id]
            start = row['sstart']
            end = row['send']
            seq_start = min(start, end)
            seq_end = max(start, end)
            
            if seq_start < 1:
                seq_start = 1
            if seq_end > len(record.seq):
                seq_end = len(record.seq)
                
            gene_seq = record.seq[seq_start - 1: seq_end]
            if start > end:
                gene_seq = gene_seq.reverse_complement()
                
            # 修改FASTA头格式为: >sampleID_geneName_geneID|contig_info
            gene_id = gene_name  # 获取基因基本名称 (如matK)
            out_file.write(f">{sample_id}_{gene_name}|{fasta_contig_id}_{seq_start}_{seq_end}\n{gene_seq}\n")

def assert_writable_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    testfile = path / "_write_test.tmp"
    try:
        with open(testfile, "w") as f:
            f.write("ok")
        testfile.unlink(missing_ok=True)
        return path
    except Exception as e:
        raise PermissionError(f"[BLAST] 目录不可写：{path}（{e}）")

# -----------------------------------------
# 核心流程

class Pipeline:
    def __init__(self, args):
        self.args = args
        self.output_root = Path(args.output)
        # ① 仅当用户提供 --genes 时才读取
        if args.genes and Path(args.genes).exists():
            self.genes = self.load_gene_list()
        else:
            self.genes = []      # 先空置；稍后 collect_all_genes 填充

        self.enable_intergenic = args.enable_intergenic

        self.species_table = None
        if args.use_ref_db:
            self.species_table = pd.read_csv(
                args.species_csv, sep=None, engine="python", index_col=0
            )

    # ---------- 新增：真正把 collect_all_genes 放到类里 ----------
    def collect_all_genes(self):
        """
        扫描 complete/ 与 incomplete/ 两端结果，自动汇总出现过的基因/间隔区
        """
        genes = set()
        for path in (self.output_root / "complete").glob("*/*/*.fasta"):
            genes.add(path.stem.split("_", 1)[-1])
        for path in (self.output_root / "incomplete").glob("*/*/*.fasta"):
            genes.add(path.stem.split("_", 1)[-1])
        self.genes = sorted(genes)


    def load_gene_list(self):
        """
        若 --genes 未提供/文件不存在，则返回 []
        """
        if not self.args.genes or not Path(self.args.genes).exists():
            return []
        with open(self.args.genes) as f:
            return [line.strip() for line in f if line.strip()]

    def create_dirs(self):
        safe_mkdir(self.output_root / "complete" / "shared")
        safe_mkdir(self.output_root / "incomplete" / "shared")
        for gene in self.genes:
            safe_mkdir(self.output_root / "complete" / gene)
            safe_mkdir(self.output_root / "incomplete" / gene)
        safe_mkdir(self.output_root / "merged_fasta")
        safe_mkdir(self.output_root / "alignments")
        safe_mkdir(self.output_root / "trimmed")
        safe_mkdir(self.output_root / "blast_dbs")


    def process_complete_samples(self):
        print("\nProcessing complete samples...")
        complete_dir = Path(self.args.complete)
        
        gb_files = list(complete_dir.glob("*.gb"))
        if not gb_files:
            print("未在 complete 目录中发现 .gb 文件")
            return

        for gb_file in gb_files:
            sample_id = gb_file.stem  # e.g., 401 → 作为 sample ID
            shared_dir = self.output_root / "complete" / "shared" / sample_id
            safe_mkdir(shared_dir)

            config = {
                "out": shared_dir.as_posix(),
                "input": str(gb_file),
                "soft_boundary": ",".join(map(str, self.args.soft_boundary)),
                "max_seq_length": self.args.max_seq_length,
                "min_seq_length": self.args.min_seq_length,
                "thread_number": self.args.threads,
                "intron_only": self.args.intron_only
            }

            extractor = Extract_reference(config)
            extractor.extract_reference_from_gb_parallel()
            org_seq_dir = Path(config["out"]) / "org_seq"

            for gene_file in org_seq_dir.glob("*.fasta"):
                gene_name = gene_file.stem
                # ——全部复制——
                output_dir = self.output_root / "complete" / gene_name / sample_id
                safe_mkdir(output_dir)
                shutil.copy(gene_file, output_dir / f"{sample_id}_{gene_name}.fasta")

            
            # -------- 如果启用提取间隔区 --------
            if self.args.enable_intergenic:
                intergenic_tmp = shared_dir / "intergenic_seq"
                safe_mkdir(intergenic_tmp)
                cmd = [
                    sys.executable,
                    str(Path(__file__).with_name("split_intergenic_2.py")),
                    "-input", str(gb_file),                 
                    "-out_dir", str(intergenic_tmp),        
                    "--soft_boundary", *map(str, self.args.soft_boundary),  
                ]
                if self.args.alias_mode:
                    cmd.append("--enable_alias")        

                subprocess.run(cmd, check=True)

                for fasta in intergenic_tmp.rglob("*.fasta"):
                    region_name = fasta.stem
                    out_dir = self.output_root / "complete" / region_name / sample_id
                    safe_mkdir(out_dir)
                    shutil.copy(fasta, out_dir / f"{sample_id}_{region_name}.fasta")
            # ----------------------------------



    def process_incomplete_samples(self):
        print("\nProcessing incomplete samples...")
        incomplete_dir = Path(self.args.incomplete)
        sample_dirs = [d for d in incomplete_dir.iterdir() if d.is_dir()]

        # 若启用数据库模式，则提前加载参考数据库与物种表
        if self.args.use_ref_db and self.species_table is None:
            raise ValueError("启用 --use-ref-db 时，无法读取 --species-csv")

        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            scaffold_fasta = sample_dir / "finalscaffold.fasta"
            reference_gb  = sample_dir / "reference.gb"
            shared_dir    = self.output_root / "incomplete" / "shared" / sample_id
            safe_mkdir(shared_dir)

            # 如果 reference.gb 不存在，则尝试从数据库中检索
            if not reference_gb.exists():
                if self.args.use_ref_db:
                    if sample_id not in self.species_table.index:
                        print(f"[WARNING] 物种表中缺少样本 {sample_id} 的记录，跳过")
                        continue
                    try:
                        species_name = self.species_table.loc[sample_id, "species_name"]
                        best_ref = find_best_reference(species_name, self.ref_db)
                        shutil.copy(best_ref, reference_gb)
                        print(f"[INFO] 样本 {sample_id} 自动引入参考：{os.path.basename(best_ref)}")
                    except Exception as e:
                        print(f"[WARNING] 样本 {sample_id} 引入数据库参考失败：{e}")
                        continue
                else:
                    print(f"[ERROR] {reference_gb} 不存在，且未启用数据库参考，跳过样本 {sample_id}")
                    continue

            # 1) 拆分参考
            config = {
                "out": shared_dir.as_posix(),
                "input": str(reference_gb),
                "soft_boundary": ",".join(map(str, self.args.soft_boundary)),
                "max_seq_length": self.args.max_seq_length,
                "min_seq_length": self.args.min_seq_length,
                "thread_number": self.args.threads,
                "intron_only": self.args.intron_only
            }
            extractor = Extract_reference(config)
            extractor.extract_reference_from_gb_parallel()
            org_seq_dir = Path(config["out"]) / "org_seq"

            # 2) 拆分间隔区（若启用）
            if self.args.enable_intergenic:
                intergenic_tmp = shared_dir / "intergenic_seq"
                safe_mkdir(intergenic_tmp)
                cmd = [
                    sys.executable,
                    str(Path(__file__).with_name("split_intergenic_2.py")),
                    "-input", str(reference_gb),
                    "-out_dir", str(intergenic_tmp),
                    "--soft_boundary", *map(str, self.args.soft_boundary),
                ]
                if self.args.alias_mode:
                    cmd.append("--enable_alias")
                subprocess.run(cmd, check=True)

            # 3) 建库（scaffold）
            blast_db_dir = self.output_root / "blast_dbs"
            assert_writable_dir(blast_db_dir)  # 确认可写；不可写直接抛错更清楚

            db_prefix = f"{sample_id}_db"      # 只保留前缀名
            makeblastdb = blast_exe("makeblastdb")

            # 关键点：把工作目录切到 blast_db_dir，-out 只给前缀名
            subprocess.run(
                [
                    makeblastdb,
                    "-in", str(scaffold_fasta),
                    "-dbtype", "nucl",
                    "-out", db_prefix
                ],
                check=True,
                cwd=str(blast_db_dir)
            )

            # 统一得到库的绝对前缀路径，供后续 blastn 使用
            blast_db_prefix = blast_db_dir / db_prefix


            # 4) BLAST + match + extract target
            all_fragments = list(org_seq_dir.glob("*.fasta"))
            if self.args.enable_intergenic:
                all_fragments += list((shared_dir / "intergenic_seq").rglob("*.fasta"))

            for gene_file in all_fragments:
                gene_name = gene_file.stem
                renamed_path = rename_fasta_header(gene_file, gene_name)
                output_dir = self.output_root / "incomplete" / gene_name / sample_id
                safe_mkdir(output_dir)

                blast_output = output_dir / "blast_results.tsv"
                blastn = blast_exe("blastn")
                subprocess.run([
                    blastn,
                    "-query", str(renamed_path),
                    "-db", str(blast_db_prefix),   # 注意这里
                    "-evalue", "1e-5",
                    "-outfmt", "6",
                    "-out", str(blast_output),
                    "-max_target_seqs", "5",
                    "-num_threads", str(self.args.threads)
                ], check=True)


                match_tsv   = output_dir / "match_result.tsv"
                missing_txt = output_dir / "missing_list.txt"
                match_blast_test(
                    str(renamed_path),
                    str(blast_output),
                    str(match_tsv),
                    str(missing_txt)
                )

                target_fasta = output_dir / f"{sample_id}_{gene_name}.fasta"
                target_gene(
                    reference_fasta=str(renamed_path),
                    blast_tsv=str(blast_output),
                    scaffold_fasta=str(scaffold_fasta),
                    output_fasta=str(target_fasta),
                    gene_name=gene_name
                )




    def merge_fasta(self):
        """
        合并所有样本的同名基因，并对缺失样本自动用 gap 补齐，
        输出目录：{output_root}/merged_fasta/all_samples_{gene}.fasta
        """
        print("\nMerging FASTA files (with gap padding for missing genes)...")

        #  关键：再保险，确保 merged_fasta 存在
        merged_root = self.output_root / "merged_fasta"
        merged_root.mkdir(parents=True, exist_ok=True)

        # ----------- 1. 收集全局样本 ID 列表 -----------
        all_sample_ids = set()

        # 完整样本
        shared_complete = self.output_root / "complete" / "shared"
        if shared_complete.exists():
            for d in shared_complete.iterdir():
                if d.is_dir():
                    all_sample_ids.add(d.name)

        # 不完整样本
        shared_incomplete = self.output_root / "incomplete" / "shared"
        if shared_incomplete.exists():
            for d in shared_incomplete.iterdir():
                if d.is_dir():
                    all_sample_ids.add(d.name)

        # ----------- 2. 按基因逐一处理 -----------
        for gene in self.genes:
            # 用字典保存每个样本的 SeqRecord（缺失先置 None）
            sample_rec_map = {sid: None for sid in all_sample_ids}
            max_len = 0      # 记录该基因最长序列长度

            # 2-1 读取「完整样本」已存在序列
            for fasta_path in (self.output_root / "complete" / gene).rglob(f"*_{gene}.fasta"):
                sample_id = fasta_path.parent.name          # 上一级目录就是 sample_id
                rec = next(SeqIO.parse(fasta_path, "fasta"))
                rec.id = f"{sample_id}_{gene}"
                rec.description = ""
                sample_rec_map[sample_id] = rec
                max_len = max(max_len, len(rec.seq))

            # 2-2 读取「不完整样本」已存在序列
            for fasta_path in (self.output_root / "incomplete" / gene).rglob(f"*_{gene}.fasta"):
                sample_id = fasta_path.parent.name
                rec = next(SeqIO.parse(fasta_path, "fasta"))
                rec.id = f"{sample_id}_{gene}"
                rec.description = ""
                sample_rec_map[sample_id] = rec
                max_len = max(max_len, len(rec.seq))

            # 2-3 对缺失样本生成 gap 序列
            for sid, rec in sample_rec_map.items():
                if rec is None:
                    gap_seq = "-" * max_len if max_len > 0 else "-"
                    sample_rec_map[sid] = SeqRecord(
                        Seq(gap_seq),
                        id=f"{sid}_{gene}",
                        description=""
                    )

            # ----------- 3. 写出合并文件 -----------
            merged_file = merged_root / f"all_samples_{gene}.fasta"
            # 关键：写之前确保父目录存在（防止被外部清空/误删）
            merged_file.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(sample_rec_map.values(), merged_file, "fasta")
            print(f"完成合并（含 gap 补齐）: {merged_file}")


    def generate_summary(self):
        all_records = []
        for gene in self.genes:
            incomplete_dirs = (self.output_root / "incomplete" / gene).glob("*")
            for sample_dir in incomplete_dirs:
                sample_id = sample_dir.name
                match_file = sample_dir / "match_result.tsv"
                if match_file.exists():
                    df = pd.read_csv(match_file, sep="\t")
                    for _, row in df.iterrows():
                        all_records.append({
                            "Sample": sample_id,
                            "Gene": row['id'],
                            "Status": row['status'],
                            "Identity": row.get('pident', ''),
                            "Coverage": row.get('coverage', '')
                        })
        pd.DataFrame(all_records).to_csv(self.output_root / "summary_report.tsv", sep="\t", index=False)


    def run_alignments(self):
        print("\nRunning MAFFT alignments...")
        merged_dir = self.output_root / "merged_fasta"
        alignment_dir = self.output_root / "alignments"
        for fasta_file in merged_dir.glob("*.fasta"):
            gene_name = fasta_file.stem.replace("all_samples_", "")
            output_file = alignment_dir / f"{gene_name}_aligned.fasta"
            run_mafft_xplat(fasta_file, output_file, self.args.threads)


    def trim_alignments(self):
        print("\nTrimming alignment ends using trimAl (automated1 mode)...")
        alignment_dir = self.output_root / "alignments"
        trimmed_dir = self.output_root / "trimmed"
        safe_mkdir(trimmed_dir)

        for aligned_file in alignment_dir.glob("*_aligned.fasta"):
            gene_name = aligned_file.stem.replace("_aligned", "")
            output_file = trimmed_dir / f"{gene_name}_trimmed.fasta"
            try:
                run_trimal_xplat(aligned_file, output_file, mode="automated1")
                print(f"切齐完成: {output_file}")
            except subprocess.CalledProcessError as e:
                print(f"trimAl 处理失败: {aligned_file}（错误信息: {e}）")



    def run_barcoding_report(self):
        """
        Step-X: 对所有 *_aligned.fasta 运行 barcoding_report.py，
        生成各类变异位点报告与 summary_table.csv
        """
        # 优先用切齐后的结果；若不存在则退回 alignments/
        align_dir = self.output_root / "trimmed"
        if not align_dir.exists() or not any(align_dir.glob("*.fasta")):
            align_dir = self.output_root / "alignments"

        out_dir = self.output_root / "variant_reports"
        safe_mkdir(out_dir)

        script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "barcoding_report.py")  # 使用绝对路径防止找不到文件
        cmd = [
            sys.executable,
            script_path,
            "-input", str(align_dir),
            "-pattern", str(self.args.pattern),
            "-output", str(out_dir)
        ]
        print(f"\nRunning barcoding_report.py on {align_dir} → {out_dir}")
        subprocess.run(cmd, check=True)

    def run_filter_summary(self):
        """
        Step-X: 调用 report_filtered.py，过滤 summary_table.csv → summary_table_filtered.csv
        """
        report_dir = self.output_root / "variant_reports"
        raw_summary = report_dir / "summary_table.csv"
        filtered_summary = report_dir / "summary_table_filtered.csv"

        if not raw_summary.exists():
            print(f"[ERROR] summary_table.csv not found in {report_dir}")
            sys.exit(1)

        print(f"\nFiltering summary_table.csv → summary_table_filtered.csv")
        filter_markers(
            str(raw_summary), str(filtered_summary),
            sample_size=self.args.sample_size,
            seq_pct=self.args.seq_pct,
            aln_len_min=self.args.aln_len_min,
            gap_pct_max=self.args.gap_pct_max,
            var_sites=(self.args.var_sites_min, self.args.var_sites_max),
            var_pct=(self.args.var_pct_min, self.args.var_pct_max),
            pi_pct_min=self.args.pi_pct_min
        )

    def run_phylogenetic_tree(self, mode="coalescent"):
        """
        Step-X: 调用 run_tree.py 构建系统发育树，输入为 trimmed + summary_table_filtered.csv
        """

        trimmed_dir = self.output_root / "trimmed"
        report_dir = self.output_root / "variant_reports"
        if self.args.report:
            filtered_summary = Path(self.args.report)
        else:
            filtered_summary = report_dir / "summary_table_filtered.csv"

        tree_output_dir = self.output_root / "tree_results"
        safe_mkdir(tree_output_dir)

        base_dir = os.path.dirname(os.path.realpath(__file__))
        run_tree_script = os.path.join(base_dir, "run_tree_win.py" if os.name == "nt" else "run_tree.py")

        cmd = [
            sys.executable,
            run_tree_script,
            "-m", mode,
            "-r", str(filtered_summary),
            "-a", str(trimmed_dir),
            "-o", str(tree_output_dir),
            "-t", str(self.args.threads)
        ]
        print(f"\nRunning run_tree.py to construct species tree → {tree_output_dir}")
        subprocess.run(cmd, check=True)


# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="完整融合版 CPTOOLS")
    parser.add_argument("--complete", required=True, help="完整样本目录")
    parser.add_argument("--incomplete", required=True, help="不完整样本目录")
    parser.add_argument("--output", required=True, help="输出目录")
    parser.add_argument("--genes", required=False, default=None, help="基因列表文件")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip-extract", action="store_true")
    parser.add_argument("--skip-alignment", action="store_true")
    parser.add_argument("--skip-trimming", action="store_true")
    parser.add_argument("--use-ref-db", action="store_true", help="若样本目录中未提供参考 gb 文件，则自动从数据库中选择参考")
    parser.add_argument("--species-csv", type=str, default=None, help="提供样本物种名称表格（CSV），包含 sample_id 和 species_name 两列")
    

    #split_genes.py
    parser.add_argument("--soft_boundary", nargs=2, type=int, default=[200,200], metavar=("LEFT","RIGHT"), help="基因区域软边界（格式：左,右）")
    parser.add_argument("--max_seq_length", type=int, default=5000, help="最大基因序列长度")
    parser.add_argument("--min_seq_length", type=int, default=200, help="最小基因序列长度")
    parser.add_argument("--intron_only", action="store_true", help="是否提取 intron 区段")

    #split_intergenic_2.py
    parser.add_argument("--enable_intergenic", action="store_true", help="是否提取经典间隔区区域")
    parser.add_argument("--alias_mode", action="store_true", help="是否启用 split_intergenic_2.py 的模糊匹配模式")

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

    safe_mkdir(Path(args.output))
    pipeline = Pipeline(args)

    if args.use_ref_db:
        if not args.species_csv:
            raise ValueError("启用 --use-ref-db 时，必须提供 --species-csv 文件。")
        if not os.path.exists(args.species_csv):
            raise FileNotFoundError(f"指定的物种信息文件不存在：{args.species_csv}")
        ref_db = load_species_database("info_list_cp.tsv")  # 数据库路径可视情况修改
        pipeline.ref_db = ref_db

        pipeline.create_dirs()

    if not args.skip_extract:
        pipeline.process_complete_samples()
        pipeline.process_incomplete_samples()
        pipeline.collect_all_genes()
        pipeline.merge_fasta()
        pipeline.generate_summary()
    if not args.skip_alignment:
        pipeline.run_alignments()
    if not args.skip_trimming:
        pipeline.trim_alignments()
    if not args.skip_report:
        pipeline.run_barcoding_report()

    if args.filtered:
        pipeline.run_filter_summary()
    if args.tree:
        pipeline.run_phylogenetic_tree(args.tree_mode)

    print("\n 全流程完成")

if __name__ == "__main__":
    main()
