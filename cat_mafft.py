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

def get_gene_names_from_org_seq(org_seq_dir):
    return sorted([Path(f).stem for f in glob.glob(f"{org_seq_dir}/*.fasta")])

def parse_species_name(result_dir):
    return Path(result_dir).name.replace(" ", "_")

def read_gene_sequence(file_path, new_id):
    records = list(SeqIO.parse(file_path, "fasta"))
    if not records:
        return None
    record = records[0]
    record.id = new_id
    record.name = ""
    record.description = ""
    return record

def create_gap_record(gene_name, species_id, length):
    return SeqRecord(Seq("-" * length), id=species_id, name="", description="")

def run_mafft(input_fasta, output_fasta):
    try:
        subprocess.run(["mafft", "--auto", input_fasta], stdout=open(output_fasta, "w"), stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        print(f"[Warning] MAFFT failed on {input_fasta}")

def main(root_dir):
    root_dir = Path(root_dir).resolve()
    org_seq_dir = root_dir / "split_genes" / "org_seq"
    alignments_dir = root_dir / "alignments"
    alignments_dir.mkdir(exist_ok=True)

    gene_names = get_gene_names_from_org_seq(org_seq_dir)
    sample_dirs = [p for p in root_dir.iterdir() if (p / "results").is_dir()]

    for gene in gene_names:
        merged_records = []
        max_len = 0

        for sample_dir in sample_dirs:
            species = parse_species_name(sample_dir)
            result_file = sample_dir / "results" / f"{gene}.fasta"
            if result_file.exists():
                rec = read_gene_sequence(result_file, f"{species}_{gene}")
                if rec:
                    merged_records.append(rec)
                    max_len = max(max_len, len(rec.seq))
                else:
                    print(f"[Warning] {result_file} is empty")
            else:
                merged_records.append(None)

        # 填充缺失序列
        filled_records = []
        for i, rec in enumerate(merged_records):
            if rec is None:
                species = parse_species_name(sample_dirs[i])
                filled = create_gap_record(gene, f"{species}_{gene}", max_len)
                filled_records.append(filled)
            else:
                # 对于较短序列补 gap（防止不一致）
                padded_seq = rec.seq + "-" * (max_len - len(rec.seq))
                rec.seq = padded_seq
                filled_records.append(rec)

        # 合并临时文件
        temp_merge = alignments_dir / f"{gene}_temp.fasta"
        SeqIO.write(filled_records, temp_merge, "fasta")

        # MAFFT 比对
        aligned_output = alignments_dir / f"{gene}_aligned.fasta"
        run_mafft(str(temp_merge), str(aligned_output))

        # 删除临时合并文件
        temp_merge.unlink()

        print(f"[Done] {gene} alignment saved to {aligned_output}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge and align gene sequences across samples")
    parser.add_argument("--result_dir", required=True, help="路径为整个结果目录，例如 D:/Atest/GM5/result3")
    args = parser.parse_args()

    main(args.result_dir)
