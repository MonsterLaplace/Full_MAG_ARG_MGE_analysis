#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
import numpy as np


def safe_read_tsv(file_path, required=True):
    if not os.path.exists(file_path):
        if required:
            raise FileNotFoundError(f"[ERROR] File not found: {file_path}")
        return pd.DataFrame()
    return pd.read_csv(file_path, sep="\t")


def ensure_columns(df, required_cols, df_name):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"[ERROR] Missing columns in {df_name}: {missing}")


def fill_arg_defaults(df):
    """
    针对当前整合后的 contig_arg.tsv 补默认值
    """
    str_cols = [
        "protein", "ARG", "ARG_name", "Class", "Subclass",
        "Resistance_Mechanism", "AMR_Gene_Family", "Source",
        "ARG_std", "ARG_name_std", "ARG_group", "Standardization_source",
        "AMRFinder_ARG", "AMRFinder_ARG_name", "AMRFinder_Type",
        "AMRFinder_Subtype", "AMRFinder_Class", "AMRFinder_Subclass",
        "AMRFinder_Method", "RGI_ARO", "RGI_Drug_Class",
        "RGI_Resistance_Mechanism", "RGI_AMR_Gene_Family",
        "RGI_Cut_Off", "RGI_Model_type", "RGI_Antibiotic"
    ]
    for col in str_cols:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str)

    int_cols = [
        "Validated_by_AMRFinder",
        "Validated_by_RGI",
        "Validated_by_both",
        "High_confidence"
    ]
    for col in int_cols:
        if col not in df.columns:
            df[col] = 0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    float_cols = [
        "AMRFinder_Identity", "AMRFinder_Coverage",
        "RGI_Identity", "RGI_Coverage"
    ]
    for col in float_cols:
        if col not in df.columns:
            df[col] = np.nan
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def fill_mge_defaults(df):
    int_cols = [
        "plasmid_like", "virus_like", "chromosome_like",
        "integron_present", "integrase_present",
        "attc_count", "complete_integron", "in0_count",
        "is_element_count", "is_family_count", "is_cluster_count",
        "isescan_type_count", "transposase_count",
        "tir_present_count", "IS_present", "TIR_present"
    ]
    for col in int_cols:
        if col not in df.columns:
            df[col] = 0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    float_cols = [
        "plasmid_score", "virus_score", "chromosome_score",
        "isescan_score_max", "isescan_score_mean", "isescan_score_min",
        "isescan_evalue_min", "isescan_evalue_mean",
        "isescan_orf_len_mean", "isescan_orf_len_max"
    ]
    for col in float_cols:
        if col not in df.columns:
            df[col] = np.nan
        df[col] = pd.to_numeric(df[col], errors="coerce")

    str_cols = [
        "genomad_label", "mge_summary",
        "is_family_list", "is_cluster_list", "isescan_type_list"
    ]
    for col in str_cols:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str)

    return df


def harmonize_genomad_flags(df):
    """
    根据 genomad_label 纠正 plasmid_like / virus_like / chromosome_like
    支持:
      plasmid, plasmid_like
      virus, viral, virus_like
      chromosome, chromosome_like
    """
    if "genomad_label" not in df.columns:
        return df

    labels = df["genomad_label"].fillna("").astype(str).str.strip().str.lower()

    if "plasmid_like" not in df.columns:
        df["plasmid_like"] = 0
    if "virus_like" not in df.columns:
        df["virus_like"] = 0
    if "chromosome_like" not in df.columns:
        df["chromosome_like"] = 0

    plasmid_mask = labels.isin(["plasmid", "plasmid_like"])
    virus_mask = labels.isin(["virus", "viral", "virus_like"])
    chromosome_mask = labels.isin(["chromosome", "chromosome_like"])

    df.loc[plasmid_mask, "plasmid_like"] = 1
    df.loc[virus_mask, "virus_like"] = 1
    df.loc[chromosome_mask, "chromosome_like"] = 1

    df["plasmid_like"] = pd.to_numeric(df["plasmid_like"], errors="coerce").fillna(0).astype(int)
    df["virus_like"] = pd.to_numeric(df["virus_like"], errors="coerce").fillna(0).astype(int)
    df["chromosome_like"] = pd.to_numeric(df["chromosome_like"], errors="coerce").fillna(0).astype(int)

    return df


def build_loose_flag(row):
    cond = (
        row["plasmid_like"] == 1 or
        row["virus_like"] == 1 or
        row["integron_present"] == 1 or
        row["integrase_present"] == 1 or
        row["IS_present"] == 1 or
        row["transposase_count"] > 0 or
        row["TIR_present"] == 1
    )
    return 1 if cond else 0


def build_strict_flag(row):
    cond = (
        row["integron_present"] == 1 or
        row["integrase_present"] == 1 or
        (row["plasmid_like"] == 1 and (
            row["IS_present"] == 1 or row["transposase_count"] > 0 or row["TIR_present"] == 1
        )) or
        (row["virus_like"] == 1 and (
            row["IS_present"] == 1 or row["transposase_count"] > 0 or row["TIR_present"] == 1
        ))
    )
    return 1 if cond else 0


def build_loose_evidence(row):
    tags = []
    if row["plasmid_like"] == 1:
        tags.append("plasmid_like")
    if row["virus_like"] == 1:
        tags.append("virus_like")
    if row["integron_present"] == 1:
        tags.append("integron")
    if row["integrase_present"] == 1:
        tags.append("integrase")
    if row["IS_present"] == 1:
        tags.append("IS")
    if row["transposase_count"] > 0:
        tags.append("transposase")
    if row["TIR_present"] == 1:
        tags.append("TIR")
    return ";".join(tags) if tags else "none"


def build_strict_evidence(row):
    tags = []
    if row["integron_present"] == 1:
        tags.append("integron")
    if row["integrase_present"] == 1:
        tags.append("integrase")
    if row["plasmid_like"] == 1 and (
        row["IS_present"] == 1 or row["transposase_count"] > 0 or row["TIR_present"] == 1
    ):
        tags.append("plasmid+IS/transposase/TIR")
    if row["virus_like"] == 1 and (
        row["IS_present"] == 1 or row["transposase_count"] > 0 or row["TIR_present"] == 1
    ):
        tags.append("virus+IS/transposase/TIR")
    return ";".join(tags) if tags else "none"


def unique_join(values):
    vals = [str(x) for x in values if pd.notnull(x) and str(x) != "" and str(x) != "nan"]
    vals = sorted(set(vals))
    return ";".join(vals)


def make_contig_summary(df):
    rows = []

    for (sample, contig), subdf in df.groupby(["Sample", "Contig"]):
        row = {
            "Sample": sample,
            "Contig": contig,

            "protein_count": subdf["protein"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "protein_list": unique_join(subdf["protein"]),

            "ARG_count": subdf["ARG"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "ARG_list": unique_join(subdf["ARG"]),
            "ARG_name_list": unique_join(subdf["ARG_name"]),

            "ARG_std_count": subdf["ARG_std"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "ARG_std_list": unique_join(subdf["ARG_std"]),
            "ARG_name_std_list": unique_join(subdf["ARG_name_std"]),
            "ARG_group_count": subdf["ARG_group"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "ARG_group_list": unique_join(subdf["ARG_group"]),

            "Class_count": subdf["Class"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "Class_list": unique_join(subdf["Class"]),
            "Subclass_count": subdf["Subclass"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "Subclass_list": unique_join(subdf["Subclass"]),
            "Resistance_Mechanism_list": unique_join(subdf["Resistance_Mechanism"]),
            "AMR_Gene_Family_list": unique_join(subdf["AMR_Gene_Family"]),

            "high_conf_ARG_count": int((subdf["High_confidence"] == 1).sum()),
            "validated_by_amrfinder_any": int((subdf["Validated_by_AMRFinder"] == 1).any()),
            "validated_by_rgi_any": int((subdf["Validated_by_RGI"] == 1).any()),
            "validated_by_both_any": int((subdf["Validated_by_both"] == 1).any()),

            "putative_mobile_ARG_loose": int((subdf["putative_mobile_ARG_loose"] == 1).any()),
            "putative_mobile_ARG_strict": int((subdf["putative_mobile_ARG_strict"] == 1).any()),
        }

        # 汇总 mobile evidence
        loose_evs = sorted(set([x for x in subdf["mobile_evidence_loose"].astype(str) if x != "none"]))
        strict_evs = sorted(set([x for x in subdf["mobile_evidence_strict"].astype(str) if x != "none"]))

        row["mobile_evidence_loose"] = ";".join(loose_evs) if loose_evs else "none"
        row["mobile_evidence_strict"] = ";".join(strict_evs) if strict_evs else "none"

        # MGE字段（同一contig基本一致）
        mge_cols = [
            "plasmid_score", "virus_score", "chromosome_score", "genomad_label",
            "plasmid_like", "virus_like", "chromosome_like",
            "integron_present", "integrase_present", "attc_count",
            "complete_integron", "in0_count",
            "is_element_count", "is_family_list", "is_family_count",
            "is_cluster_list", "is_cluster_count",
            "isescan_type_list", "isescan_type_count",
            "isescan_score_max", "isescan_score_mean", "isescan_score_min",
            "isescan_evalue_min", "isescan_evalue_mean",
            "transposase_count", "isescan_orf_len_mean", "isescan_orf_len_max",
            "tir_present_count", "IS_present", "TIR_present",
            "mge_summary"
        ]

        for col in mge_cols:
            if col in subdf.columns:
                row[col] = subdf.iloc[0][col]

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Merge standardized contig ARG table with contig MGE table."
    )
    parser.add_argument("--arg", required=True, help="standardized contig_arg.tsv")
    parser.add_argument("--mge", required=True, help="all_samples.contig_mge.tsv")
    parser.add_argument("-o", "--output", required=True, help="output contig_arg_mge.tsv")
    parser.add_argument("--summary", required=False, help="optional output contig_arg_mge_summary.tsv")
    args = parser.parse_args()

    # 读取
    arg_df = safe_read_tsv(args.arg, required=True)
    mge_df = safe_read_tsv(args.mge, required=True)

    # 检查关键列
    ensure_columns(arg_df, ["Sample", "Contig", "protein", "ARG"], "contig_arg.tsv")
    ensure_columns(mge_df, ["Sample", "Contig"], "contig_mge.tsv")

    arg_df = fill_arg_defaults(arg_df)
    mge_df = fill_mge_defaults(mge_df)

    # 先对MGE表做一次genomad标签修正
    mge_df = harmonize_genomad_flags(mge_df)

    # 合并
    merged = pd.merge(arg_df, mge_df, on=["Sample", "Contig"], how="left")

    # 对没有匹配到MGE的ARG记录补默认值
    merged = fill_arg_defaults(merged)
    merged = fill_mge_defaults(merged)

    # 合并后再次根据genomad_label纠正 plasmid/virus/chromosome 标记
    merged = harmonize_genomad_flags(merged)

    # ARG存在标志
    merged["ARG_present"] = 1

    # mobile ARG 判定
    merged["putative_mobile_ARG_loose"] = merged.apply(build_loose_flag, axis=1)
    merged["putative_mobile_ARG_strict"] = merged.apply(build_strict_flag, axis=1)

    merged["mobile_evidence_loose"] = merged.apply(build_loose_evidence, axis=1)
    merged["mobile_evidence_strict"] = merged.apply(build_strict_evidence, axis=1)

    # 输出列顺序
    output_cols = [
        "Sample", "Contig", "protein",
        "ARG", "ARG_name", "Class", "Subclass",
        "Resistance_Mechanism", "AMR_Gene_Family",

        "ARG_std", "ARG_name_std", "ARG_group", "Standardization_source",

        "Validated_by_AMRFinder", "Validated_by_RGI",
        "Validated_by_both", "High_confidence", "Source",

        "AMRFinder_ARG", "AMRFinder_ARG_name", "AMRFinder_Type",
        "AMRFinder_Subtype", "AMRFinder_Class", "AMRFinder_Subclass",
        "AMRFinder_Method", "AMRFinder_Identity", "AMRFinder_Coverage",

        "RGI_ARO", "RGI_Drug_Class", "RGI_Resistance_Mechanism",
        "RGI_AMR_Gene_Family", "RGI_Identity", "RGI_Coverage",
        "RGI_Cut_Off", "RGI_Model_type", "RGI_Antibiotic",

        "ARG_present",

        "plasmid_score", "virus_score", "chromosome_score", "genomad_label",
        "plasmid_like", "virus_like", "chromosome_like",

        "integron_present", "integrase_present", "attc_count",
        "complete_integron", "in0_count",

        "is_element_count", "is_family_list", "is_family_count",
        "is_cluster_list", "is_cluster_count",
        "isescan_type_list", "isescan_type_count",
        "isescan_score_max", "isescan_score_mean", "isescan_score_min",
        "isescan_evalue_min", "isescan_evalue_mean",
        "transposase_count", "isescan_orf_len_mean", "isescan_orf_len_max",
        "tir_present_count", "IS_present", "TIR_present",

        "mge_summary",
        "putative_mobile_ARG_loose",
        "putative_mobile_ARG_strict",
        "mobile_evidence_loose",
        "mobile_evidence_strict"
    ]

    for col in output_cols:
        if col not in merged.columns:
            merged[col] = np.nan

    merged = merged[output_cols].sort_values(["Sample", "Contig", "protein", "ARG"])
    merged.to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Output written to: {args.output}")
    print(f"[INFO] Total ARG records: {merged.shape[0]}")

    if args.summary:
        summary_df = make_contig_summary(merged)
        summary_df = summary_df.sort_values(["Sample", "Contig"])
        summary_df.to_csv(args.summary, sep="\t", index=False)
        print(f"[INFO] Summary written to: {args.summary}")
        print(f"[INFO] Total ARG-carrying contigs: {summary_df.shape[0]}")


if __name__ == "__main__":
    main()
