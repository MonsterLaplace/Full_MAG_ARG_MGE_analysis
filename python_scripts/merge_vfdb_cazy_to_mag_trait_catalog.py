#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import argparse
import pandas as pd


def extract_vfg_id(s):
    """
    从第二列中提取 VFG 编号，例如：
    VFG018662(gb|WP_012130507) -> VFG018662
    """
    if pd.isna(s):
        return None
    m = re.match(r'^(VFG\d+)', str(s))
    return m.group(1) if m else None


def parse_vfdb_file(fp):
    """
    解析单个 VFDB tsv 文件
    假定无表头，至少前2列：
      col1 = query gene
      col2 = VFDB/VFG ID
    """
    mag = os.path.basename(fp).rsplit(".tsv", 1)[0]

    try:
        df = pd.read_csv(fp, sep='\t', header=None, comment='#')
    except pd.errors.EmptyDataError:
        return {
            "MAG": mag,
            "vfdb_hit_count": 0,
            "vfdb_gene_count": 0,
            "vfdb_vfg_count": 0,
            "vfdb_vfg_list": "",
            "has_vfdb_hit": 0,
        }

    if df.empty or df.shape[1] < 2:
        return {
            "MAG": mag,
            "vfdb_hit_count": 0,
            "vfdb_gene_count": 0,
            "vfdb_vfg_count": 0,
            "vfdb_vfg_list": "",
            "has_vfdb_hit": 0,
        }

    query_genes = df.iloc[:, 0].astype(str).dropna()
    vfg_ids = df.iloc[:, 1].astype(str).map(extract_vfg_id).dropna()

    unique_genes = sorted(set(query_genes))
    unique_vfgs = sorted(set(vfg_ids))

    return {
        "MAG": mag,
        "vfdb_hit_count": int(len(df)),
        "vfdb_gene_count": int(len(unique_genes)),
        "vfdb_vfg_count": int(len(unique_vfgs)),
        "vfdb_vfg_list": "; ".join(unique_vfgs),
        "has_vfdb_hit": 1 if len(df) > 0 else 0,
    }


def load_vfdb_summary(vfdb_dir):
    """
    汇总目录下所有 VFDB 结果文件
    """
    vfdb_files = sorted(glob.glob(os.path.join(vfdb_dir, "*.tsv")))
    if not vfdb_files:
        print(f"[WARN] No TSV files found in: {vfdb_dir}")

    records = [parse_vfdb_file(fp) for fp in vfdb_files]
    vfdb_df = pd.DataFrame(records)

    if vfdb_df.empty:
        vfdb_df = pd.DataFrame(columns=[
            "MAG",
            "vfdb_hit_count",
            "vfdb_gene_count",
            "vfdb_vfg_count",
            "vfdb_vfg_list",
            "has_vfdb_hit"
        ])

    return vfdb_df, vfdb_files


def main():
    parser = argparse.ArgumentParser(
        description="Merge VFDB and CAZy results into mag_trait_catalog.tsv"
    )
    parser.add_argument(
        "-i", "--input_table",
        default="13_tables/mag_trait_catalog.tsv",
        help="Input MAG trait table"
    )
    parser.add_argument(
        "-v", "--vfdb_dir",
        default="10_annotation/vfdb",
        help="Directory containing VFDB result TSV files"
    )
    parser.add_argument(
        "-c", "--cazy_table",
        default="10_annotation/cazy_summary/mag_cazy_catalog.tsv",
        help="CAZy summary table"
    )
    parser.add_argument(
        "-o", "--output_table",
        default="13_tables/mag_trait_catalog.with_vfdb_cazy.tsv",
        help="Output merged table"
    )
    args = parser.parse_args()

    # 读取主表
    mag_df = pd.read_csv(args.input_table, sep="\t", dtype=str)
    if "MAG" not in mag_df.columns:
        raise ValueError("Input table must contain a 'MAG' column.")

    # 读取 VFDB 汇总
    vfdb_df, vfdb_files = load_vfdb_summary(args.vfdb_dir)

    # 合并 VFDB
    merged = mag_df.merge(vfdb_df, on="MAG", how="left")

    # 填充 VFDB 缺失值
    vfdb_zero_cols = [
        "vfdb_hit_count",
        "vfdb_gene_count",
        "vfdb_vfg_count",
        "has_vfdb_hit"
    ]
    for col in vfdb_zero_cols:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    if "vfdb_vfg_list" in merged.columns:
        merged["vfdb_vfg_list"] = merged["vfdb_vfg_list"].fillna("")

    # 读取 CAZy 表
    cazy_df = pd.read_csv(args.cazy_table, sep="\t", dtype=str)
    if "MAG" not in cazy_df.columns:
        raise ValueError("CAZy table must contain a 'MAG' column.")

    # 避免与已有列冲突：只合并 cazy 表中除 MAG 外的列
    cazy_cols = [col for col in cazy_df.columns if col != "MAG"]
    merged = merged.merge(cazy_df[["MAG"] + cazy_cols], on="MAG", how="left")

    # CAZy 数值列填 0
    cazy_zero_cols = [
        "CAZy_count",
        "CAZy_family_count",
        "GH_count",
        "GT_count",
        "CE_count",
        "CBM_count",
        "PL_count",
        "AA_count",
        "GH_family_count",
        "GT_family_count",
        "CE_family_count",
        "CBM_family_count",
        "PL_family_count",
        "AA_family_count",
    ]
    for col in cazy_zero_cols:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    # CAZy 列表字段填空字符串
    if "CAZy_family_list" in merged.columns:
        merged["CAZy_family_list"] = merged["CAZy_family_list"].fillna("")

    # 输出
    merged.to_csv(args.output_table, sep="\t", index=False)

    print(f"[OK] Output written to: {args.output_table}")
    print(f"[INFO] Input MAG table rows: {len(mag_df)}")
    print(f"[INFO] VFDB files parsed: {len(vfdb_files)}")
    print(f"[INFO] CAZy table rows: {len(cazy_df)}")


if __name__ == "__main__":
    main()
