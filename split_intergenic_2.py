#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
import re

def normalize_trn_name(name):
    name = name.replace("(", "").replace(")", "").replace("_", "-")
    match = re.match(r"^(trn\w)-?(\w{3})$", name, re.IGNORECASE)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    return name.strip()


def get_feature_name(feature):
    for key in ['gene', 'locus_tag', 'product']:
        if key in feature.qualifiers:
            return feature.qualifiers[key][0].strip()
    return None

def sanitize_name(name):
    return name.replace(" ", "_").replace("|", "_").replace("/", "_")

def extract_region(seq, start, end, strand):
    if strand == -1:
        return seq[start:end].reverse_complement()
    else:
        return seq[start:end]

def soft_boundary(start, end, total_length, left, right):
    new_start = (start - left) % total_length
    new_end = (end + right) % total_length
    if new_start < new_end:
        return new_start, new_end
    else:
        return new_start, new_end + total_length  # wrap around origin

def calculate_circular_distance(pos1, pos2, total_length):
    """计算环形基因组中两个位置的最短距离"""
    linear_dist = abs(pos1 - pos2)
    circular_dist = total_length - linear_dist
    return min(linear_dist, circular_dist)

def find_fallback_gene(base, gene_map, anchor_pos, total_length, max_distance=3000):
    candidates = []
    base_norm = normalize_trn_name(base).lower()
    print(f"[Debug] Looking for fallback of '{base}' (normalized: '{base_norm}') around position {anchor_pos}")
    
    # 先收集所有可能匹配的基因，不考虑距离
    potential_matches = []
    for name, info in gene_map.items():
        name_norm = normalize_trn_name(name).lower()
        raw_name = info.get("raw_name", "").lower()
        
        match_info = []
        if base_norm in name_norm:
            match_info.append("base in name")
        if name_norm in base_norm:
            match_info.append("name in base")
        if base_norm in raw_name:
            match_info.append("base in raw")
        if raw_name in base_norm:
            match_info.append("raw in base")
        
        if match_info:
            center = (info["start"] + info["end"]) // 2
            dist = calculate_circular_distance(center, anchor_pos, total_length)
            potential_matches.append((name, center, dist, match_info))
            print(f"  [Potential] {name} (start={info['start']}, end={info['end']}, center={center}, circular_dist={dist}) matched by: {', '.join(match_info)}")

    # 从潜在匹配中筛选距离合适的候选
    for name, center, dist, match_info in potential_matches:
        if dist <= max_distance:
            print(f"  [Candidate] {name} (center={center}, circular_dist={dist}) matched by: {', '.join(match_info)}")
            candidates.append((name, center, dist))

    if not candidates:
        print(f"[Debug] No fallback candidates found within {max_distance}bp. Total potential matches: {len(potential_matches)}")
        if potential_matches:
            closest = min(potential_matches, key=lambda x: x[2])
            print(f"[Debug] Closest match was: {closest[0]} at distance {closest[2]}")
        return None

    best = min(candidates, key=lambda x: x[2])
    print(f"[Debug] → Selected fallback: {best[0]} (circular_dist={best[2]})")
    return best[0]


def find_unique_gene_match(base, gene_map):
    matches = []
    base_norm = normalize_trn_name(base).lower()
    for name, info in gene_map.items():
        name_norm = normalize_trn_name(name).lower()
        raw_name = info.get("raw_name", "").lower()

        if base_norm in name_norm or base_norm in raw_name:
            matches.append(name)

    if matches:
        print(f"[Unique Match] base='{base}' → found {len(matches)} candidates: {matches}")
    else:
        print(f"[Unique Match] base='{base}' → no match found")
    return matches[0] if len(matches) == 1 else None


INTERGENIC_REGIONS = {
    "trnH-psbA": ("trnH-GUG", "psbA"),
    "atpF-atpH": ("atpF", "atpH"),
    "psbK-psbI": ("psbK", "psbI"),
    "rpoB-trnC(GCA)": ("rpoB", "trnC-GCA"),
    "rpl32-trnL(UAG)": ("rpl32", "trnL-UAG"),
    "ndhF-rpl32": ("ndhF", "rpl32"),
    "petA-psbJ": ("petA", "psbJ"),
    "psaA-ycf3": ("psaA", "ycf3"),
    "ycf4-cemA": ("ycf4", "cemA"),
    "psaC-ndhE": ("psaC", "ndhE"),
    "rps16-trnQ(UUG)": ("rps16", "trnQ-UUG"),
    "trnS(UGA)-trnG(UCC)": ("trnS-UGA", "trnG-UCC"),
    "trnL(UAA)-trnF(GAA)": ("trnL-UAA", "trnF-GAA"),
    "ndhC-trnV(UAC)": ("ndhC", "trnV-UAC"),
    "rpl16-rps3": ("rpl16", "rps3"),
}

def strict_match(target, gene_map):
    norm_target = normalize_trn_name(target)
    for key in gene_map:
        if normalize_trn_name(key) == norm_target:
            print(f"[Strict Match] '{target}' matched strictly as '{key}'")
            return key
    print(f"[Strict Match] '{target}' NOT found in gene_map (normalized as '{norm_target}')")
    return None


def process_genbank_file(file_path, out_dir, soft_left, soft_right, enable_alias):
    records = list(SeqIO.parse(file_path, "genbank"))
    if not records:
        print(f"[Error] No GenBank records found in {file_path}.", file=sys.stderr)
        return

    record = records[0]
    seq = record.seq
    total_len = len(seq)
    species = sanitize_name(record.annotations.get("organism", "Unknown"))

    missing_regions = []

    gene_map = {}
    for feature in record.features:
        if feature.type == "gene":
            name = get_feature_name(feature)
            if name:
                raw_name = name
                norm_name = normalize_trn_name(name)
                gene_map[norm_name] = {
                    "start": int(feature.location.start),
                    "end": int(feature.location.end),
                    "strand": feature.location.strand,
                    "raw_name": raw_name
                }

    for region_name, (g1, g2) in INTERGENIC_REGIONS.items():
        g1_key = strict_match(g1, gene_map)
        g2_key = strict_match(g2, gene_map)

        fallback_used = False

        # 如果启用别名模式，进行模糊匹配
        if enable_alias:
            # 提取基因名的基础部分（去掉反密码子部分）
            g1_base = g1.split("-")[0] if "-" in g1 else g1.replace("(", "").replace(")", "")
            g2_base = g2.split("-")[0] if "-" in g2 else g2.replace("(", "").replace(")", "")

            if g1_key and not g2_key:
                # g1严格匹配到，尝试模糊匹配g2
                anchor_pos = (gene_map[g1_key]["start"] + gene_map[g1_key]["end"]) // 2
                g2_key_fallback = find_fallback_gene(g2_base, gene_map, anchor_pos, total_len, max_distance=5000)
                if g2_key_fallback:
                    g2_key = g2_key_fallback
                    fallback_used = True
                    print(f"[Debug] Found fallback for g2: {g2_base} → {g2_key}")

            elif g2_key and not g1_key:
                # g2严格匹配到，尝试模糊匹配g1
                anchor_pos = (gene_map[g2_key]["start"] + gene_map[g2_key]["end"]) // 2
                g1_key_fallback = find_fallback_gene(g1_base, gene_map, anchor_pos, total_len, max_distance=5000)
                if g1_key_fallback:
                    g1_key = g1_key_fallback
                    fallback_used = True
                    print(f"[Debug] Found fallback for g1: {g1_base} → {g1_key}")

            elif not g1_key and not g2_key:
                # 两端都没有严格匹配到，尝试先确定唯一命中的一端
                g1_candidate = find_unique_gene_match(g1_base, gene_map)
                g2_candidate = find_unique_gene_match(g2_base, gene_map)
                
                if g1_candidate and not g2_candidate:
                    # 找到g1的唯一候选，以它为锚点查找g2
                    anchor_pos = (gene_map[g1_candidate]["start"] + gene_map[g1_candidate]["end"]) // 2
                    g2_key_fallback = find_fallback_gene(g2_base, gene_map, anchor_pos, total_len, max_distance=5000)
                    if g2_key_fallback:
                        g1_key = g1_candidate
                        g2_key = g2_key_fallback
                        fallback_used = True
                        print(f"[Debug] Using unique g1 candidate '{g1_candidate}' as anchor for g2 fallback")
                
                elif g2_candidate and not g1_candidate:
                    # 找到g2的唯一候选，以它为锚点查找g1
                    anchor_pos = (gene_map[g2_candidate]["start"] + gene_map[g2_candidate]["end"]) // 2
                    g1_key_fallback = find_fallback_gene(g1_base, gene_map, anchor_pos, total_len, max_distance=5000)
                    if g1_key_fallback:
                        g2_key = g2_candidate
                        g1_key = g1_key_fallback
                        fallback_used = True
                        print(f"[Debug] Using unique g2 candidate '{g2_candidate}' as anchor for g1 fallback")
                
                elif g1_candidate and g2_candidate:
                    # 两个都有唯一候选，直接使用
                    g1_key = g1_candidate
                    g2_key = g2_candidate
                    fallback_used = True
                    print(f"[Debug] Using both unique candidates: g1='{g1_candidate}', g2='{g2_candidate}'")

        # 输出两端基因中心位置和距离
        if g1_key and g2_key:
            g1_center = (gene_map[g1_key]["start"] + gene_map[g1_key]["end"]) // 2
            g2_center = (gene_map[g2_key]["start"] + gene_map[g2_key]["end"]) // 2
            distance = calculate_circular_distance(g1_center, g2_center, total_len)
        else:
            g1_center = g2_center = distance = "N/A"

        print(f"[Fallback Summary] Region: {region_name}, g1_key: {g1_key}, g2_key: {g2_key}, used fallback: {fallback_used}, centers: {g1_center}~{g2_center}, distance: {distance}")

        if not g1_key or not g2_key:
            # 提取基因名的基础部分用于相似性检查
            g1_base = g1.split("-")[0] if "-" in g1 else g1.replace("(", "").replace(")", "")
            g2_base = g2.split("-")[0] if "-" in g2 else g2.replace("(", "").replace(")", "")
            similar = [k for k in gene_map if g1_base.lower() in normalize_trn_name(k).lower() or g2_base.lower() in normalize_trn_name(k).lower()]
            msg = f"[Warning] Skipping {region_name}: missing gene {g1} or {g2}; similar in file: {similar}"
            if not similar:
                msg += " (no similar genes found)"
            print(msg, file=sys.stderr)
            missing_regions.append(msg)
            continue
        elif fallback_used:
            print(f"[Info] Fallback used for {region_name}: matched {g1_key}, {g2_key} instead of {g1}, {g2} [alias mode enabled]")

        start = gene_map[g1_key]["end"]
        end = gene_map[g2_key]["start"]
        if start > end:
            start, end = end, start

        inter_start, inter_end = soft_boundary(start, end, total_len, soft_left, soft_right)

        if inter_end <= total_len:
            region_seq = extract_region(seq, inter_start, inter_end, 1)
        else:
            part1 = extract_region(seq, inter_start, total_len, 1)
            part2 = extract_region(seq, 0, inter_end % total_len, 1)
            region_seq = part1 + part2

        species_dir = os.path.join(out_dir, species)
        os.makedirs(species_dir, exist_ok=True)
        out_path = os.path.join(species_dir, f"{region_name}.fasta")
        header = f">{species}|{region_name}|{g1_key}~{g2_key}|{inter_start}_{inter_end}"
        with open(out_path, "w") as f:
            f.write(f"{header}\n{region_seq}\n")

        print(f"[Info] Extracted: {region_name} → {out_path}")

    if missing_regions:
        species_dir = os.path.join(out_dir, species)
        log_path = os.path.join(species_dir, "missing_regions.log")
        with open(log_path, "w") as logf:
            logf.write("\n".join(missing_regions) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract chloroplast intergenic regions from GenBank files or folders")
    parser.add_argument("-input", required=True, nargs='+', help="One or more GenBank files or folders containing them")
    parser.add_argument("-out_dir", required=True, help="Output directory")
    parser.add_argument("--soft_boundary", nargs=2, type=int, default=[200, 200], help="Soft boundary left and right")
    parser.add_argument("--enable_alias", action="store_true", help="Enable fallback to nearby similarly named genes when exact match is missing")

    global args
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    gb_files = []
    for path in args.input:
        if os.path.isfile(path) and path.endswith((".gb", ".gbk")):
            gb_files.append(path)
        elif os.path.isdir(path):
            for fname in os.listdir(path):
                full_path = os.path.join(path, fname)
                if os.path.isfile(full_path) and fname.endswith((".gb", ".gbk")):
                    gb_files.append(full_path)
        else:
            print(f"[Warning] Ignoring invalid input: {path}", file=sys.stderr)

    if not gb_files:
        sys.exit("[Error] No valid GenBank files found.")

    print(f"[Info] Total input files: {len(gb_files)}")
    if args.enable_alias:
        print("[Info] Fallback alias mode is ENABLED (±3kbp, partial name match)")
    else:
        print("[Info] Fallback alias mode is DISABLED (only exact gene names used)")

    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = {
            executor.submit(
                process_genbank_file,
                gb_file,
                args.out_dir,
                args.soft_boundary[0],
                args.soft_boundary[1],
                args.enable_alias
            ): gb_file
            for gb_file in gb_files
        }
        for future in as_completed(futures):
            gb_file = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"[Error] Failed processing {gb_file}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()