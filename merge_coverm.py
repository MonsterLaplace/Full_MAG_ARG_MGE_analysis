#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
import sys
import glob
import re


def detect_genome_column(df):
    """
    自动识别Genome列
    """
    possible = ["Genome", "genome", "Bin", "bin", "MAG", "mag"]
    for c in df.columns:
        if c in possible:
            return c

    for c in df.columns:
        cl = c.lower()
        if "genome" in cl or "bin" in cl or "mag" in cl:
            return c

    return None


def detect_value_column(df, metric):
    """
    自动识别要提取的数值列
    """
    # 优先精确匹配
    for c in df.columns:
        if c == metric:
            return c

    # 再宽松匹配
    for c in df.columns:
        if c.lower() == metric.lower():
            return c

    return None


def extract_sample_name(filepath, sample_regex=None):
    """
    从文件名提取 sample_id
    默认去掉扩展名
    可用正则自定义
    """
    fname = os.path.basename(filepath)
    sample = re.sub(r'\.tsv$|\.txt$|\.csv$', '', fname)

    if sample_regex:
        m = re.search(sample_regex, sample)
        if m:
            if m.groups():
                return m.group(1)
            else:
                return m.group(0)

    return sample


def read_coverm_file(filepath, metric, sample_regex=None):
    """
    读取单个 coverm 文件并转成长表
    输出列:
      sample_id, MAG, value
    """
    try:
        df = pd.read_csv(filepath, sep="\t")
    except Exception as e:
        raise RuntimeError(f"Failed to read {filepath}: {e}")

    genome_col = detect_genome_column(df)
    value_col = detect_value_column(df, metric)

    if genome_col is None:
        raise ValueError(f"[ERROR] Cannot detect genome column in {filepath}. Columns: {list(df.columns)}")
    if value_col is None:
        raise ValueError(f"[ERROR] Cannot detect metric column '{metric}' in {filepath}. Columns: {list(df.columns)}")

    sample_id = extract_sample_name(filepath, sample_regex=sample_regex)

    sub = df[[genome_col, value_col]].copy()
    sub.columns = ["MAG", "value"]
    sub["sample_id"] = sample_id

    # 数值化
    sub["value"] = pd.to_numeric(sub["value"], errors="coerce").fillna(0)

    return sub[["sample_id", "MAG", "value"]]


def merge_coverm_files(file_list, metric="relative_abundance", sample_regex=None):
    """
    合并多个 CoverM 文件，输出：
      long_df: sample_id, MAG, value
      wide_df: sample x MAG
    """
    all_long = []

    for f in file_list:
        print(f"[INFO] Reading: {f}", file=sys.stderr)
        sub = read_coverm_file(f, metric=metric, sample_regex=sample_regex)
        all_long.append(sub)

    if not all_long:
        raise ValueError("[ERROR] No input files parsed successfully.")

    long_df = pd.concat(all_long, ignore_index=True)

    wide_df = long_df.pivot_table(
        index="sample_id",
        columns="MAG",
        values="value",
        aggfunc="sum",
        fill_value=0
    ).reset_index()

    return long_df, wide_df


def main():
    parser = argparse.ArgumentParser(
        description="Merge CoverM genome output files into sample x MAG matrix."
    )
    parser.add_argument("inputs", nargs="+", help="Input CoverM TSV files, or glob-expanded file list")
    parser.add_argument("-m", "--metric", default="relative_abundance",
                        help="Metric column to extract from CoverM output (default: relative_abundance)")
    parser.add_argument("--sample_regex", default=None,
                        help="Optional regex to extract sample_id from filename")
    parser.add_argument("--long_out", default=None,
                        help="Optional long-format output file")
    parser.add_argument("--sort_mags", action="store_true",
                        help="Sort MAG columns alphabetically")
    parser.add_argument("-o", "--output", default=None,
                        help="Wide-format output file. Default: stdout")
    args = parser.parse_args()

    # 支持 shell 未展开时自己glob
    file_list = []
    for x in args.inputs:
        if any(ch in x for ch in ["*", "?", "["]):
            file_list.extend(glob.glob(x))
        else:
            file_list.append(x)

    file_list = sorted(set(file_list))

    if len(file_list) == 0:
        print("[ERROR] No input files found.", file=sys.stderr)
        sys.exit(1)

    long_df, wide_df = merge_coverm_files(
        file_list=file_list,
        metric=args.metric,
        sample_regex=args.sample_regex
    )

    if args.sort_mags:
        mag_cols = sorted([c for c in wide_df.columns if c != "sample_id"])
        wide_df = wide_df[["sample_id"] + mag_cols]

    if args.long_out:
        long_df.to_csv(args.long_out, sep="\t", index=False)
        print(f"[INFO] Long-format output written to: {args.long_out}", file=sys.stderr)

    if args.output:
        wide_df.to_csv(args.output, sep="\t", index=False)
        print(f"[INFO] Wide-format output written to: {args.output}", file=sys.stderr)
    else:
        wide_df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
