# find_reference.py

import os
import pandas as pd

def load_species_database(tsv_file):
    return pd.read_csv(tsv_file, sep="\t")

def find_best_reference(species_name, db, db_root="cp_gb"):
    # 精确匹配种名
    matches = db[db["ORGANISM_GENUS_SPECIES"].str.lower() == species_name.lower()]
    
    # 若找不到，则尝试属名匹配
    if matches.empty:
        genus = species_name.split()[0]
        matches = db[db["ORGANISM_GENUS_SPECIES"].str.lower().str.startswith(genus.lower())]
        if matches.empty:
            raise ValueError(f"未找到匹配物种或属：{species_name}")
    
    # 优先选择 NC_ 开头的 accession
    if any(matches["ACCESSION"].str.startswith("NC_")):
        matches = matches[matches["ACCESSION"].str.startswith("NC_")]

    # 按序列长度降序排序
    matches = matches.sort_values(by="SEQ_LEN", ascending=False)
    acc = matches.iloc[0]["ACCESSION"]

    # 构造路径
    prefix = acc[:2]
    subdir = acc[:6]
    gb_path = os.path.join(db_root, prefix, subdir, f"{acc}.gb")

    if not os.path.exists(gb_path):
        raise FileNotFoundError(f"文件未找到：{gb_path}")
    return gb_path
