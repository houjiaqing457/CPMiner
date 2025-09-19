#!/usr/bin/env python3
# cli.py


import argparse
import sys
import os
from core import main as core_main
from CPTOOLS8 import main as cptools_main

def main():
    parser = argparse.ArgumentParser(
        description="CPMiner: A toolkit for phylogenetic reconstruction from raw reads or assembled sequences."
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # 模块一：从测序数据出发
    parser_reads = subparsers.add_parser(
        "reads2tree", help="Full pipeline from raw reads to phylogenetic tree"
    )
    parser_reads.add_argument("--ref", help="Input GenBank file (or auto-select with --species-name)")
    parser_reads.add_argument("-f", "--file-list", required=True, help="Sample list in TSV format")
    parser_reads.add_argument("--outdir", required=True, help="Output directory")
    parser_reads.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser_reads.add_argument("--enable-intergenic", action="store_true", help="Enable intergenic region extraction")
    parser_reads.add_argument("--soft-boundary", default="200,200", help="Soft boundary for intergenic extraction")
    parser_reads.add_argument("--alias-mode", action="store_true", help="Enable fuzzy matching for intergenic")
    parser_reads.add_argument("--species-name", help="Species name for auto-selecting reference")
    parser_reads.add_argument("--filtered", action="store_true", help="Run report_filtered.py")
    parser_reads.add_argument("--tree", action="store_true", help="Build phylogenetic tree")
    parser_reads.add_argument("--tree-mode", choices=["concat", "parallel"], default="concat", help="Tree mode")

    # 模块二：从组装数据出发
    parser_assembled = subparsers.add_parser(
        "assembled2tree", help="Build tree from assembled and aligned sequences"
    )
    parser_assembled.add_argument("--summary", required=True, help="Filtered summary_table CSV file")
    parser_assembled.add_argument("--alignment", required=True, help="Alignment folder containing *_trimmed.fasta")
    parser_assembled.add_argument("--outdir", required=True, help="Output directory for tree results")
    parser_assembled.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser_assembled.add_argument("--mode", choices=["concat", "parallel"], default="concat", help="Tree construction mode")

    args = parser.parse_args()

    if args.command == "reads2tree":
        sys.argv = ["core.py"]
        for key, value in vars(args).items():
            if value is True:
                sys.argv.append(f"--{key.replace('_', '-')}")
            elif value not in [False, None]:
                sys.argv.extend([f"--{key.replace('_', '-')}", str(value)])
        core_main()

    elif args.command == "assembled2tree":
        sys.argv = ["CPTOOLS8.py"]
        for key, value in vars(args).items():
            sys.argv.extend([f"--{key.replace('_', '-')}", str(value)])
        cptools_main()

if __name__ == "__main__":
    main()
