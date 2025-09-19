#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import shutil
import platform
from pathlib import Path
import re
from collections import defaultdict, Counter
import pandas as pd


def get_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='分析FASTA文件中的变异位点')
    parser.add_argument('-input', default='input', help='输入目录 (默认: input)')
    parser.add_argument('-pattern', default='10', help='分析模式 (1-10, 默认: 10)')
    parser.add_argument('-output', default='output', help='输出目录 (默认: output)')
    return parser.parse_args()

def default_value(default_val, key):
    """返回默认值函数，模拟原始Perl的default函数"""
    args = get_arguments()
    return getattr(args, key, default_val) or default_val

def setup_output_directory(output_directory):
    """设置输出目录"""
    osname = platform.system().lower()
    
    if os.path.exists(output_directory):
        if osname == "windows":
            os.system(f'del /f /s /q "{output_directory}"')
        elif osname in ["linux", "darwin"]:
            shutil.rmtree(output_directory)
        elif osname == "cygwin":
            os.system(f'rm -rf "{output_directory}"')
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

def find_fasta_files(input_directory):
    """查找FASTA文件"""
    fasta_extensions = ['.fasta', '.fas', '.fa', '.fsa']
    filenames = []
    
    for root, dirs, files in os.walk(input_directory):
        for file in files:
            if any(file.endswith(ext) for ext in fasta_extensions):
                filenames.append(os.path.join(root, file))
    
    return filenames

def read_fasta_sequences(filename):
    """读取FASTA文件并返回序列矩阵"""
    sequences = []
    ids = []
    current_seq = []
    current_id = None
    
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append(list(''.join(current_seq)))
                    ids.append(current_id)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        # 处理最后一个序列
        if current_id is not None:
            sequences.append(list(''.join(current_seq)))
            ids.append(current_id)
    
    return sequences, ids

def transpose_matrix(matrix):
    """转置矩阵"""
    if not matrix or not matrix[0]:
        return []
    return list(map(list, zip(*matrix)))

def count_unique_bases(column, include_gaps=False):
    """统计一列中唯一碱基的数量"""
    if include_gaps:
        unique_bases = set(column)
    else:
        unique_bases = set(base for base in column if base not in ['-', 'N', '?'])
    return len(unique_bases)

def is_parsimony_informative(column):
    """判断是否为简约信息位点"""
    # 统计每个碱基的频次（不包括gap）
    base_counts = Counter(base for base in column if base not in ['-', 'N', '?'])
    
    # 至少有2种碱基，每种至少出现2次
    informative_bases = [base for base, count in base_counts.items() if count >= 2]
    return len(informative_bases) >= 2

def has_gaps(column):
    """判断是否包含gap字符"""
    return any(base in ['-', 'N', '?'] for base in column)

def analyze_sites_by_pattern(sequences, pattern):
    """根据模式分析位点"""
    if not sequences:
        return {}
    
    # 转置矩阵，每列代表一个位点
    transposed = transpose_matrix(sequences)
    
    results = {
        'variable_array': [],                    # pattern 1: 所有变异位点
        'two_bases_array': [],                   # pattern 2: 2种碱基的位点
        'three_bases_array': [],                 # pattern 3: 3种碱基的位点
        'four_bases_array': [],                  # pattern 4: 4种碱基的位点
        'parsimony_informative_sites_array': [], # pattern 5: 简约信息位点
        'non_parsimony_informative_sites_array': [], # pattern 6: 非信息性变异位点
        'invariable_array': [],                  # pattern 7: 保守位点
        'gap_array': []                          # pattern 8: 缺失位点
    }
    
    for col_idx, column in enumerate(transposed):
        unique_count = count_unique_bases(column, include_gaps=False)
        unique_count_with_gaps = count_unique_bases(column, include_gaps=True)
        has_gap = has_gaps(column)
        is_pis = is_parsimony_informative(column)
        
        # Pattern 1: 所有变异位点（至少2个物种有不同）
        if pattern in ['1', '9', '10'] and unique_count >= 2:
            results['variable_array'].append((col_idx, column))
        
        # Pattern 2: 仅有2种碱基的变异位点
        if pattern in ['2', '10'] and unique_count == 2:
            results['two_bases_array'].append((col_idx, column))
        
        # Pattern 3: 有3种碱基的变异位点
        if pattern in ['3', '10'] and unique_count == 3:
            results['three_bases_array'].append((col_idx, column))
        
        # Pattern 4: 有4种碱基的变异位点
        if pattern in ['4', '10'] and unique_count == 4:
            results['four_bases_array'].append((col_idx, column))
        
        # Pattern 5: 简约信息位点
        if pattern in ['5', '9', '10'] and is_pis:
            results['parsimony_informative_sites_array'].append((col_idx, column))
        
        # Pattern 6: 非信息性变异位点
        if pattern in ['6', '10'] and unique_count >= 2 and not is_pis:
            results['non_parsimony_informative_sites_array'].append((col_idx, column))
        
        # Pattern 7: 保守位点
        if pattern in ['7', '10'] and unique_count <= 1:
            results['invariable_array'].append((col_idx, column))
        
        # Pattern 8: 缺失位点
        if pattern in ['8', '10'] and has_gap:
            results['gap_array'].append((col_idx, column))
    
    return results

def write_output_files(results, ids, filename, output_directory, pattern, seq_length):
    """写入输出文件"""
    gene_output_dir = os.path.join(output_directory, filename)
    os.makedirs(gene_output_dir, exist_ok=True)

    def write_array_to_files(array_data, array_name, pattern_num, seq_length):
        if array_data:
            transposed_array = transpose_matrix([col for _, col in array_data])
            alignment_length = seq_length  # 来自 analyze_sites_by_pattern 里返回的总长
            num_sites = len(array_data)
            percent_sites = num_sites / alignment_length if alignment_length > 0 else 0
            positions = [str(col_idx + 1) for col_idx, _ in array_data]

            # 文本输出
            txt_filename = f"{gene_output_dir}/{filename}_p{pattern_num}_{array_name}.txt"
            with open(txt_filename, 'w', encoding='utf-8') as f:
                f.write(f"Total sequence length of alignment matrix {filename}:\n")
                f.write(f"{alignment_length}\n\n")

                f.write(f"Number of {array_name.replace('_', ' ')} in alignment matrix {filename}:\n")
                f.write(f"{num_sites}\n\n")

                f.write(f"Percentage of {array_name.replace('_', ' ')} in alignment matrix {filename}:\n")
                f.write(f"{percent_sites:.4f}\n\n")

                f.write(f"Position of {array_name.replace('_', ' ')} in alignment matrix {filename}:\n")
                f.write(' '.join(positions) + '\n\n')

                f.write(f"{array_name.replace('_', ' ').capitalize()} in alignment matrix {filename}:\n")
                for i, seq_id in enumerate(ids):
                    f.write(f">{seq_id}\n")
                    f.write(''.join(transposed_array[i]) + "\n")

            
            # 写入FASTA文件
            fasta_filename = f"{gene_output_dir}/{filename}_p{pattern_num}_{array_name}.fasta"
            with open(fasta_filename, 'w', encoding='utf-8') as f:
                for i, seq_id in enumerate(ids):
                    f.write(f">{seq_id}\n")
                    f.write(''.join(transposed_array[i]) + "\n")
    
    # 根据模式写入相应的文件
    if pattern in ['1', '9', '10']:
        write_array_to_files(results['variable_array'], 'variable_sites_allnucleotides', '1', seq_length)
    
    if pattern in ['2', '10']:
        write_array_to_files(results['two_bases_array'], 'two_bases_sites', '2', seq_length)
    
    if pattern in ['3', '10']:
        write_array_to_files(results['three_bases_array'], 'three_bases_sites', '3', seq_length)
    
    if pattern in ['4', '10']:
        write_array_to_files(results['four_bases_array'], 'four_bases_sites', '4', seq_length)
    
    if pattern in ['5', '9', '10']:
        write_array_to_files(results['parsimony_informative_sites_array'], 'parsimony_informative_sites', '5', seq_length)
    
    if pattern in ['6', '10']:
        write_array_to_files(results['non_parsimony_informative_sites_array'], 'non-parsimony_informative_sites', '6', seq_length)
    
    if pattern in ['7', '10']:
        write_array_to_files(results['invariable_array'], 'invariable_sites', '7', seq_length)
    
    if pattern in ['8', '10']:
        write_array_to_files(results['gap_array'], 'gap_sites', '8', seq_length)

def preprocess_fasta(input_filename, output_filename):
    """预处理FASTA文件，模拟原始Perl的格式化"""
    with open(input_filename, 'r', encoding='utf-8') as input_file, \
         open(output_filename, 'w', encoding='utf-8') as output_file:
        
        first_line = input_file.readline()
        output_file.write(first_line)
        
        for line in input_file:
            line = line.replace('\r', '').replace('\n', '')
            if line.startswith('>'):
                output_file.write('\n' + line + '\n')
            else:
                output_file.write(line)
        output_file.write('\n')

def calculate_statistics(results, total_length):
    """计算统计信息"""
    stats = {}
    
    for key, array_data in results.items():
        count = len(array_data)
        percentage = (count / total_length * 100) if total_length > 0 else 0
        stats[key] = {
            'count': count,
            'percentage': percentage
        }
    
    return stats

def main():
    """主函数"""
    # 解析参数
    args = get_arguments()
    
    input_directory = args.input or "input"
    pattern = args.pattern or "10"
    output_directory = args.output or "output"
    
    # 验证pattern参数
    if pattern not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']:
        print(f"错误: 不支持的分析模式 {pattern}")
        print("支持的模式:")
        print("1 - 所有变异位点")
        print("2 - 仅有2种碱基的变异位点")
        print("3 - 有3种碱基的变异位点")
        print("4 - 有4种碱基的变异位点")
        print("5 - 简约信息位点")
        print("6 - 非信息性变异位点")
        print("7 - 保守位点")
        print("8 - 缺失位点")
        print("9 - 模式1+5")
        print("10 - 全部模式")
        return
    
    # 设置输出目录
    setup_output_directory(output_directory)
    
    # 查找FASTA文件
    filenames = find_fasta_files(input_directory)
    
    if not filenames:
        print(f"在目录 {input_directory} 中未找到FASTA文件")
        return
    
    # 统计变量
    percentage_variable = {}
    percentage_pis = {}
    number_variable = {}
    number_pis = {}
    
    print(f"找到 {len(filenames)} 个FASTA文件")
    print(f"使用分析模式: {pattern}")
    print()
    
    # 处理每个FASTA文件
    for filename_fasta in filenames:
        print(f"处理文件: {filename_fasta}")
        
        # 获取文件名（不含扩展名）
        target_name = os.path.splitext(filename_fasta)[0]
        target_name = target_name.replace(' ', '_')
        filename = os.path.basename(target_name)
        
        # 创建临时文件
        output_filename_temp = f"{input_directory}/{filename}_temp"
        
        try:
            # 预处理FASTA文件
            preprocess_fasta(filename_fasta, output_filename_temp)
            
            # 读取序列
            sequences, ids = read_fasta_sequences(output_filename_temp)
            
            if not sequences:
                print(f"警告: 文件 {filename_fasta} 中没有有效序列")
                continue
            
            # 获取序列长度
            total_length = len(sequences[0]) if sequences else 0
            
            # 分析位点
            results = analyze_sites_by_pattern(sequences, pattern)
            
            # 写入输出文件
            write_output_files(results, ids, filename, output_directory, pattern, total_length)
            
            # 计算统计信息
            stats = calculate_statistics(results, total_length)
            
            # 收集统计信息（用于兼容原始脚本）
            if 'variable_array' in stats:
                percentage_variable[filename] = stats['variable_array']['percentage']
                number_variable[filename] = stats['variable_array']['count']
            
            if 'parsimony_informative_sites_array' in stats:
                percentage_pis[filename] = stats['parsimony_informative_sites_array']['percentage']
                number_pis[filename] = stats['parsimony_informative_sites_array']['count']
            
            # 输出统计信息
            print(f"  序列数量: {len(sequences)}")
            print(f"  序列长度: {total_length}")
            
            if pattern in ['1', '9', '10']:
                print(f"  变异位点: {stats['variable_array']['count']} ({stats['variable_array']['percentage']:.2f}%)")
            if pattern in ['2', '10']:
                print(f"  2种碱基位点: {stats['two_bases_array']['count']} ({stats['two_bases_array']['percentage']:.2f}%)")
            if pattern in ['3', '10']:
                print(f"  3种碱基位点: {stats['three_bases_array']['count']} ({stats['three_bases_array']['percentage']:.2f}%)")
            if pattern in ['4', '10']:
                print(f"  4种碱基位点: {stats['four_bases_array']['count']} ({stats['four_bases_array']['percentage']:.2f}%)")
            if pattern in ['5', '9', '10']:
                print(f"  简约信息位点: {stats['parsimony_informative_sites_array']['count']} ({stats['parsimony_informative_sites_array']['percentage']:.2f}%)")
            if pattern in ['6', '10']:
                print(f"  非信息性变异位点: {stats['non_parsimony_informative_sites_array']['count']} ({stats['non_parsimony_informative_sites_array']['percentage']:.2f}%)")
            if pattern in ['7', '10']:
                print(f"  保守位点: {stats['invariable_array']['count']} ({stats['invariable_array']['percentage']:.2f}%)")
            if pattern in ['8', '10']:
                print(f"  缺失位点: {stats['gap_array']['count']} ({stats['gap_array']['percentage']:.2f}%)")
            
            print()
            
        except Exception as e:
            print(f"处理文件 {filename_fasta} 时出错: {str(e)}")
            continue
        
        finally:
            # 清理临时文件
            if os.path.exists(output_filename_temp):
                os.remove(output_filename_temp)
    
    # 输出总体统计
    print("=" * 20)
    print("处理完成!")
    print(f"共处理了 {len(filenames)} 个文件")
    print(f"分析模式: {pattern}")
    print(f"结果保存在目录: {output_directory}")

    # 根据所选 pattern 决定是否输出平均比例
    if pattern in ['1', '9', '10'] and percentage_variable:
        print(f"\n平均变异位点比例: {sum(percentage_variable.values())/len(percentage_variable):.2f}%")

    if pattern in ['5', '9', '10'] and percentage_pis:
        print(f"平均简约信息位点比例: {sum(percentage_pis.values())/len(percentage_pis):.2f}%")

    # 生成汇总统计表（只包含当前分析模式涉及的统计项）
    summary_data = []
    headers = ['Gene', 'Sequence Count', 'Alignment Length']

    # 映射模式编号到对应的统计键及人类可读名称
    pattern_stat_keys = {
        '1': [('variable_array', 'Variable Sites')],
        '2': [('two_bases_array', 'Two Base Sites')],
        '3': [('three_bases_array', 'Three Base Sites')],
        '4': [('four_bases_array', 'Four Base Sites')],
        '5': [('parsimony_informative_sites_array', 'Parsimony Informative Sites')],
        '6': [('non_parsimony_informative_sites_array', 'Non-Parsimony Sites')],
        '7': [('invariable_array', 'Invariable Sites')],
        '8': [('gap_array', 'Gap Sites')],
        '9': [
            ('variable_array', 'Variable Sites'),
            ('parsimony_informative_sites_array', 'Parsimony Informative Sites')
        ],
        '10': [
            ('variable_array', 'Variable Sites'),
            ('two_bases_array', 'Two Base Sites'),
            ('three_bases_array', 'Three Base Sites'),
            ('four_bases_array', 'Four Base Sites'),
            ('parsimony_informative_sites_array', 'Parsimony Informative Sites'),
            ('non_parsimony_informative_sites_array', 'Non-Parsimony Sites'),
            ('invariable_array', 'Invariable Sites'),
            ('gap_array', 'Gap Sites')
        ]
    }

    selected_keys = pattern_stat_keys.get(pattern, [])
    for _, readable_name in selected_keys:
        headers.extend([f"{readable_name} Count", f"{readable_name} %"])

    for filename_fasta in filenames:
        target_name = os.path.splitext(filename_fasta)[0]
        filename = os.path.basename(target_name)

        output_filename_temp = f"{input_directory}/{filename}_temp"
        preprocess_fasta(filename_fasta, output_filename_temp)
        sequences, ids = read_fasta_sequences(output_filename_temp)
        total_length = len(sequences[0]) if sequences else 0
        results = analyze_sites_by_pattern(sequences, pattern)
        stats = calculate_statistics(results, total_length)
        row = [filename, len(sequences), total_length]

        for key, _ in selected_keys:
            stat = stats.get(key, {'count': 0, 'percentage': 0})
            row.append(stat['count'])
            row.append(f"{stat['percentage']:.2f}")
        summary_data.append(row)

        if os.path.exists(output_filename_temp):
            os.remove(output_filename_temp)

    # 写出CSV文件
    summary_csv_path = os.path.join(output_directory, 'summary_table.csv')
    df_summary = pd.DataFrame(summary_data, columns=headers)
    df_summary.to_csv(summary_csv_path, index=False, encoding='utf-8-sig')
    print(f"\n已输出汇总统计表: {summary_csv_path}")


if __name__ == "__main__":
    main()