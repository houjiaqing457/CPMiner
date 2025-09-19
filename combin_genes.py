#!/usr/bin/env python3

import os
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Concatenate MSA files based on sample name.")
    parser.add_argument("-r", "--report", required=True, help="Filtered CSV file with Gene column")
    parser.add_argument("-m", "--msa_dir", required=True, help="Directory containing *_trimmed.fasta files")
    parser.add_argument("-o", "--output", required=True, help="Output concatenated alignment fasta file")
    return parser.parse_args()

def extract_sample_name(fasta_id: str, gene_name: str, valid_sample_set: set) -> str:
    parts = fasta_id.split("_")
    # 尝试倒数第1个下划线之前
    if "_".join(parts[:-1]) in valid_sample_set:
        return "_".join(parts[:-1])
    # 尝试倒数第2个下划线之前
    elif len(parts) > 2 and "_".join(parts[:-2]) in valid_sample_set:
        return "_".join(parts[:-2])
    else:
        return "_".join(parts[:-1])  # 默认策略（即使匹配不到）


def concatenate_msa(report_csv, msa_dir, output_file):
    df = pd.read_csv(report_csv)
    genes = df['Gene'].dropna().tolist()

    all_sample_ids = []
    sample_to_concat_seq = defaultdict(str)
    first = True

    for gene in genes:
        msa_path = Path(msa_dir) / f"{gene}.fasta"
        if not msa_path.exists():
            print(f"[WARNING] Missing file: {msa_path}")
            continue

        records = list(SeqIO.parse(msa_path, "fasta"))
        seq_dict = {}
        valid_sample_set = set(all_sample_ids)
        for record in records:
            sample_name = extract_sample_name(record.id, gene, valid_sample_set)
            seq_dict[sample_name] = str(record.seq)


        # 锁定样本顺序
        if first:
            all_sample_ids = list(seq_dict.keys())
            first = False

        aln_len = len(next(iter(seq_dict.values())))
        for sid in all_sample_ids:
            if sid in seq_dict:
                sample_to_concat_seq[sid] += seq_dict[sid]
            else:
                print(f"[WARNING] Sample {sid} missing in {gene}, filling with gaps.")
                sample_to_concat_seq[sid] += "-" * aln_len

    # 输出拼接文件
    with open(output_file, "w") as f:
        for sid in all_sample_ids:
            f.write(f">{sid}\n{sample_to_concat_seq[sid]}\n")


    print(f"[INFO] Concatenated alignment written to: {output_file}")

def main():
    args = parse_args()
    concatenate_msa(args.report, args.msa_dir, args.output)

if __name__ == "__main__":
    main()
