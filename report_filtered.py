import pandas as pd

def filter_markers(
    input_file,
    output_file,
    sample_size=None,
    seq_pct=0.7,
    aln_len_min=500,
    gap_pct_max=60.0,
    var_sites=(0, 200),
    var_pct=(0.0, 10.0),
    pi_pct_min=1.0,
    verbose=False
):
    df = pd.read_csv(input_file)
    original_len = len(df)

    if sample_size is None:
        if "Sequence Count" in df.columns:
            sample_size = df["Sequence Count"].max()
        else:
            raise ValueError("Cannot infer sample size: 'Sequence Count' column missing.")

    filters = []

    if "Sequence Count" in df.columns:
        min_seq_count = int(sample_size * seq_pct)
        filters.append(df["Sequence Count"] >= min_seq_count)

    if "Alignment Length" in df.columns:
        filters.append(df["Alignment Length"] >= aln_len_min)

    if "Gap Sites %" in df.columns:
        filters.append(df["Gap Sites %"] <= gap_pct_max)

    if "Variable Sites Count" in df.columns:
        filters.append(df["Variable Sites Count"].between(*var_sites))

    if "Variable Sites %" in df.columns:
        filters.append(df["Variable Sites %"].between(*var_pct))

    if "Parsimony Informative Sites %" in df.columns:
        filters.append(df["Parsimony Informative Sites %"] >= pi_pct_min)

    if filters:
        combined_filter = filters[0]
        for f in filters[1:]:
            combined_filter &= f
        df_filtered = df[combined_filter]
    else:
        df_filtered = df.copy()

    df_filtered.to_csv(output_file, index=False)  # 用逗号写出 CSV

    if verbose:
        print(f"Input rows: {original_len}")
        print(f"Filtered rows: {len(df_filtered)}")
        print(f"Output written to: {output_file}")

    return df_filtered


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Filter phylogenetic markers based on summary report.")
    parser.add_argument("-i", "--input", required=True, help="Input summary_table.csv file")
    parser.add_argument("-o", "--output", required=True, help="Output filtered marker table")
    parser.add_argument("--sample_size", type=int, default=None, help="Total sample count (optional)")
    parser.add_argument("--seq_pct", type=float, default=0.7, help="Minimum fraction of samples per marker [0–1]")
    parser.add_argument("--aln_len_min", type=int, default=500, help="Minimum alignment length")
    parser.add_argument("--gap_pct_max", type=float, default=60.0, help="Maximum gap percentage")
    parser.add_argument("--var_sites_min", type=int, default=0, help="Minimum variable sites")
    parser.add_argument("--var_sites_max", type=int, default=200, help="Maximum variable sites")
    parser.add_argument("--var_pct_min", type=float, default=0.0, help="Minimum variable sites percentage")
    parser.add_argument("--var_pct_max", type=float, default=10.0, help="Maximum variable sites percentage")
    parser.add_argument("--pi_pct_min", type=float, default=1.0, help="Minimum parsimony informative site %")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print summary")
    args = parser.parse_args()

    filter_markers(
        input_file=args.input,
        output_file=args.output,
        sample_size=args.sample_size,
        seq_pct=args.seq_pct,
        aln_len_min=args.aln_len_min,
        gap_pct_max=args.gap_pct_max,
        var_sites=(args.var_sites_min, args.var_sites_max),
        var_pct=(args.var_pct_min, args.var_pct_max),
        pi_pct_min=args.pi_pct_min,
        verbose=args.verbose
    )
