#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
import sys
import numpy as np


def find_fixed_isescan_tsv(input_path):
    """
    在输入路径中寻找优先解析的 ISEScan tsv 文件
    优先顺序：
    1. *.fixed.tsv
    2. *.clean.fa.tsv
    3. 任意 .tsv
    """
    if os.path.isfile(input_path):
        return input_path

    candidates = []
    for root, dirs, files in os.walk(input_path):
        for f in files:
            if f.endswith(".tsv"):
                candidates.append(os.path.join(root, f))

    if not candidates:
        return None

    fixed_files = [x for x in candidates if x.endswith(".fixed.tsv")]
    if fixed_files:
        return fixed_files[0]

    clean_files = [x for x in candidates if x.endswith(".clean.fa.tsv")]
    if clean_files:
        return clean_files[0]

    return candidates[0]


def safe_numeric(series):
    return pd.to_numeric(series, errors="coerce")


def parse_isescan_tsv(tsv_file):
    """
    解析 ISEScan tsv（你当前 clean.fa.tsv / fixed.tsv 格式）
    """
    df = pd.read_csv(tsv_file, sep="\t", comment="#")

    required_cols = {"seqID", "family", "cluster", "score", "E-value", "type", "orfLen", "tir"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"[ERROR] Missing required columns in {tsv_file}: {missing}")

    # 数值列转换
    numeric_cols = [
        "score", "E-value", "orfLen", "isLen", "ncopy4is",
        "irId", "irLen", "nGaps", "ov"
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = safe_numeric(df[col])

    # 标准化明细表
    detail = df.copy()
    detail = detail.rename(columns={
        "seqID": "Contig",
        "family": "IS_family",
        "cluster": "IS_cluster",
        "type": "IS_type",
        "score": "IS_score",
        "E-value": "IS_evalue",
        "orfLen": "IS_orfLen",
        "tir": "IS_tir"
    })

    # 统一字符串列
    for col in ["IS_family", "IS_cluster", "IS_type", "IS_tir"]:
        if col in detail.columns:
            detail[col] = detail[col].fillna("").astype(str)

    # 是否存在TIR
    detail["tir_present"] = detail["IS_tir"].apply(lambda x: 0 if x.strip() == "" or x.strip().lower() == "nan" else 1)

    return detail


def summarize_isescan(detail_df):
    """
    从明细表汇总为 contig-level 表
    """
    grouped = detail_df.groupby("Contig")

    summary = grouped.agg(
        is_element_count=("Contig", "size"),
        is_family_list=("IS_family", lambda x: ";".join(sorted(set([i for i in map(str, x) if i])))),
        is_family_count=("IS_family", lambda x: len(set([i for i in map(str, x) if i]))),
        is_cluster_list=("IS_cluster", lambda x: ";".join(sorted(set([i for i in map(str, x) if i])))),
        is_cluster_count=("IS_cluster", lambda x: len(set([i for i in map(str, x) if i]))),
        isescan_type_list=("IS_type", lambda x: ";".join(sorted(set([i for i in map(str, x) if i])))),
        isescan_type_count=("IS_type", lambda x: len(set([i for i in map(str, x) if i]))),
        isescan_score_max=("IS_score", "max"),
        isescan_score_mean=("IS_score", "mean"),
        isescan_score_min=("IS_score", "min"),
        isescan_evalue_min=("IS_evalue", "min"),
        isescan_evalue_mean=("IS_evalue", "mean"),
        isescan_orf_len_mean=("IS_orfLen", "mean"),
        isescan_orf_len_max=("IS_orfLen", "max"),
        tir_present_count=("tir_present", "sum")
    ).reset_index()

    # 在你当前 ISEScan 结果中，每条 IS 预测记录可以近似作为一个 transposase-associated event
    summary["transposase_count"] = summary["is_element_count"]

    # MGE摘要标签
    def make_summary(row):
        tags = []
        if row["is_element_count"] > 0:
            tags.append("IS")
        if row["transposase_count"] > 0:
            tags.append("transposase")
        if row["tir_present_count"] > 0:
            tags.append("TIR")
        return ";".join(tags) if tags else "none"

    summary["mge_summary"] = summary.apply(make_summary, axis=1)

    return summary


def main():
    parser = argparse.ArgumentParser(description="Parse fixed ISEScan TSV and generate contig-level summary.")
    parser.add_argument("-i", "--input", required=True, help="ISEScan tsv file or directory")
    parser.add_argument("-o", "--output", required=True, help="output summary tsv")
    parser.add_argument("--detail", required=False, help="optional output detail tsv")
    args = parser.parse_args()

    target = find_fixed_isescan_tsv(args.input)
    if target is None:
        print("[ERROR] No ISEScan TSV file found.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Parsing ISEScan TSV: {target}")
    detail_df = parse_isescan_tsv(target)
    summary_df = summarize_isescan(detail_df)

    summary_df.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Summary written to: {args.output}")

    if args.detail:
        detail_df.to_csv(args.detail, sep="\t", index=False)
        print(f"[INFO] Detail written to: {args.detail}")

    print(f"[INFO] Parsed contigs with IS hits: {summary_df.shape[0]}")
    print(f"[INFO] Total IS records: {detail_df.shape[0]}")


if __name__ == "__main__":
    main()
