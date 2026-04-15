#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
import sys
import glob


def find_summary_files(input_dir):
    """
    在输入目录中寻找:
    - *_plasmid_summary.tsv
    - *_virus_summary.tsv
    """
    plasmid_files = glob.glob(os.path.join(input_dir, "*_plasmid_summary.tsv"))
    virus_files = glob.glob(os.path.join(input_dir, "*_virus_summary.tsv"))

    return plasmid_files, virus_files


def detect_contig_column(df):
    """
    自动识别 contig 列名
    常见候选:
    seq_name, contig, sequence_name, name
    """
    for c in df.columns:
        cl = c.lower()
        if cl in ["seq_name", "contig", "sequence_name", "contig_name", "name"]:
            return c
    return None


def detect_score_column(df, mode="plasmid"):
    """
    自动识别 score 列
    """
    for c in df.columns:
        cl = c.lower()
        if mode == "plasmid":
            if "plasmid" in cl and "score" in cl:
                return c
        elif mode == "virus":
            if "virus" in cl and "score" in cl:
                return c

    # 若没有精确匹配，再尝试更宽松匹配
    for c in df.columns:
        cl = c.lower()
        if "score" in cl:
            return c
    return None


def parse_plasmid_file(plasmid_file):
    df = pd.read_csv(plasmid_file, sep="\t", comment="#")
    contig_col = detect_contig_column(df)
    score_col = detect_score_column(df, mode="plasmid")

    if contig_col is None:
        raise ValueError(f"[ERROR] Cannot detect contig column in {plasmid_file}")
    if score_col is None:
        raise ValueError(f"[ERROR] Cannot detect plasmid score column in {plasmid_file}")

    out = df[[contig_col, score_col]].copy()
    out.columns = ["Contig", "plasmid_score"]
    return out


def parse_virus_file(virus_file):
    df = pd.read_csv(virus_file, sep="\t", comment="#")
    contig_col = detect_contig_column(df)
    score_col = detect_score_column(df, mode="virus")

    if contig_col is None:
        raise ValueError(f"[ERROR] Cannot detect contig column in {virus_file}")
    if score_col is None:
        raise ValueError(f"[ERROR] Cannot detect virus score column in {virus_file}")

    out = df[[contig_col, score_col]].copy()
    out.columns = ["Contig", "virus_score"]
    return out


def infer_label(row, threshold=0.5):
    p = row["plasmid_score"]
    v = row["virus_score"]

    if p >= threshold and p >= v:
        return "plasmid"
    elif v >= threshold and v > p:
        return "virus"
    else:
        return "chromosome"


def main():
    parser = argparse.ArgumentParser(
        description="Parse geNomad plasmid and virus summary files and merge into a contig-level MGE table."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="geNomad sample output directory, e.g. 11_mge/genomad/sample01/")
    parser.add_argument("-o", "--output", required=True,
                        help="output parsed tsv file")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="score threshold for assigning plasmid/virus label (default: 0.5)")
    args = parser.parse_args()

    plasmid_files, virus_files = find_summary_files(args.input)

    if len(plasmid_files) == 0:
        print("[ERROR] No *_plasmid_summary.tsv found", file=sys.stderr)
        sys.exit(1)
    if len(virus_files) == 0:
        print("[ERROR] No *_virus_summary.tsv found", file=sys.stderr)
        sys.exit(1)

    # 默认各取第一个
    plasmid_file = plasmid_files[0]
    virus_file = virus_files[0]

    print(f"[INFO] Using plasmid summary: {plasmid_file}")
    print(f"[INFO] Using virus summary:   {virus_file}")

    p_df = parse_plasmid_file(plasmid_file)
    v_df = parse_virus_file(virus_file)

    merged = pd.merge(p_df, v_df, on="Contig", how="outer")

    # 缺失值填0
    merged["plasmid_score"] = merged["plasmid_score"].fillna(0.0)
    merged["virus_score"] = merged["virus_score"].fillna(0.0)

    # 构建一个近似 chromosome_score
    merged["chromosome_score"] = 1 - merged[["plasmid_score", "virus_score"]].max(axis=1)
    merged["chromosome_score"] = merged["chromosome_score"].clip(lower=0)

    # 分类标签
    merged["genomad_label"] = merged.apply(lambda row: infer_label(row, threshold=args.threshold), axis=1)

    # source 列
    merged["source"] = os.path.basename(args.input.rstrip("/"))

    # 排序
    merged = merged.sort_values("Contig")

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Done. Output written to: {args.output}")
    print(f"[INFO] Total contigs parsed: {merged.shape[0]}")


if __name__ == "__main__":
    main()
