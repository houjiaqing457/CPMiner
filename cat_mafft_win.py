#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_gene_names_from_org_seq(org_seq_dir: Path):
    return sorted([Path(f).stem for f in glob.glob(str(org_seq_dir / "*.fasta"))])

def parse_species_name(sample_dir: Path):
    # 与原逻辑一致：用样本目录名（空格转下划线）
    return sample_dir.name.replace(" ", "_")

def read_gene_sequence(file_path: Path, new_id: str):
    records = list(SeqIO.parse(str(file_path), "fasta"))
    if not records:
        return None
    record = records[0]
    record.id = new_id
    record.name = ""
    record.description = ""
    return record

def create_gap_record(species_id: str, length: int):
    return SeqRecord(Seq("-" * length), id=species_id, name="", description="")

def find_mafft_launcher(base_dir: Path):
    """
    优先返回同级目录 mafft-win/mafft.bat 的调用方式；
    若不存在，则回退 'mafft'（要求用户已配置 PATH）。
    返回 (invoke_type, launcher)
      - invoke_type == "bat"  -> 需要 shell=True 并用重定向
      - invoke_type == "exe"  -> 直接用可执行（或 PATH 中的 'mafft'），可用 stdout 重定向句柄
    """
    bat = base_dir / "mafft-win" / "mafft.bat"
    if bat.exists():
        return "bat", str(bat)
    # 允许用户已将 mafft.exe/mafft（WSL/其他环境）加入 PATH
    return "exe", "mafft"

def run_mafft(invoke_type: str, launcher: str, input_fasta: Path, output_fasta: Path):
    """
    保持与原 Linux 逻辑一致：使用 '--auto'。不新增线程或额外参数。
    Windows 下 .bat 必须 shell=True，且用 '>' 做输出重定向。
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    if invoke_type == "bat":
        # 用 cmd.exe /C 并重定向输出
        cmd = f'"{launcher}" --auto "{input_fasta}" > "{output_fasta}"'
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"[Warning] MAFFT failed on {input_fasta}")
    else:
        # 走可执行/PATH 中的 mafft，直接把 stdout 写文件
        try:
            with open(output_fasta, "w", encoding="utf-8") as fout:
                subprocess.run([launcher, "--auto", str(input_fasta)],
                               stdout=fout,
                               stderr=subprocess.DEVNULL,
                               check=True)
        except subprocess.CalledProcessError:
            print(f"[Warning] MAFFT failed on {input_fasta}")

def main(result_dir: str):
    root_dir = Path(result_dir).resolve()

    # 输入：参考基因目录
    org_seq_dir = root_dir / "split_genes" / "org_seq"
    if not org_seq_dir.is_dir():
        raise SystemExit(f"[Error] Not found: {org_seq_dir}")

    # 输出：与 Linux 版一致的两个目录
    merged_dir = root_dir / "merged_fasta"
    aligned_dir = root_dir / "alignments"
    merged_dir.mkdir(exist_ok=True)
    aligned_dir.mkdir(exist_ok=True)

    # 发现基因
    gene_names = get_gene_names_from_org_seq(org_seq_dir)
    if not gene_names:
        raise SystemExit(f"[Error] No *.fasta genes found under {org_seq_dir}")

    # 发现样本（规则：root 下含有 results 子目录）
    sample_dirs = [p for p in root_dir.iterdir() if (p / "results").is_dir()]
    if not sample_dirs:
        raise SystemExit(f"[Error] No sample 'results' dirs found under {root_dir}")

    # 准备 MAFFT 启动方式
    invoke_type, launcher = find_mafft_launcher(base_dir=Path(__file__).parent)

    print("Merging FASTA files (with gap padding for missing genes)...")
    for gene in gene_names:
        # 收集各样本的该基因序列
        records = []
        max_len = 0

        for sample_dir in sample_dirs:
            species = parse_species_name(sample_dir)
            # 原脚本路径结构：<sample>/results/<gene>.fasta
            f = sample_dir / "results" / f"{gene}.fasta"
            if f.exists():
                rec = read_gene_sequence(f, species)
                if rec:
                    records.append(rec)
                    max_len = max(max_len, len(rec.seq))
                else:
                    # 空文件按缺失处理
                    records.append(None)
            else:
                records.append(None)

        # 按最大长度补齐/填充缺失
        filled_records = []
        for i, rec in enumerate(records):
            species = parse_species_name(sample_dirs[i])
            if rec is None:
                filled_records.append(create_gap_record(species, max_len))
            else:
                if len(rec.seq) < max_len:
                    rec.seq = rec.seq + "-" * (max_len - len(rec.seq))
                filled_records.append(rec)

        # 生成与 Linux 版一致的“合并后（未比对）”文件
        merged_fa = merged_dir / f"all_samples_{gene}.fasta"
        SeqIO.write(filled_records, merged_fa, "fasta")
        print(f"已合并（gap 填充）: {merged_fa}")

    print("\nRunning MAFFT alignments...")
    for gene in gene_names:
        merged_fa = merged_dir / f"all_samples_{gene}.fasta"
        aligned_out = aligned_dir / f"{gene}_aligned.fasta"

        # 调用 MAFFT 做比对
        run_mafft(invoke_type, launcher, merged_fa, aligned_out)

        # 仅当比对文件成功生成时，删除中间合并文件
        if aligned_out.exists():
            try:
                merged_fa.unlink()  # 删除 all_samples_<gene>.fasta
            except Exception as e:
                print(f"[Warn] 删除中间文件失败：{merged_fa} ({e})")
        else:
            print(f"[Warn] 未生成对齐文件，保留中间合并文件用于排查：{merged_fa}")

        print(f"完成比对: {aligned_out}")
    # --- 可选清理：如果 merged_fasta 目录已经空了，就删除该目录 ---
    try:
        if not any(merged_dir.iterdir()):
            merged_dir.rmdir()
            print(f"[Clean] 已删除空目录：{merged_dir}")
    except Exception as e:
        print(f"[Warn] 无法删除目录 {merged_dir}：{e}")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge per-gene FASTAs across samples and run MAFFT alignment (Windows-compatible).")
    # 与原版保持一致：只保留 --result_dir
    parser.add_argument("--result_dir", required=True, help="整个结果目录，例如 D:/Atest/GM5/result3")
    args = parser.parse_args()
    main(args.result_dir)
