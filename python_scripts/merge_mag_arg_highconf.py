#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import glob
import pandas as pd
import numpy as np


def read_amrfinder_file(file_path):
    mag = os.path.basename(file_path).replace(".tsv", "")
    df = pd.read_csv(file_path, sep="\t")

    if df.empty:
        return pd.DataFrame()

    required = [
        "Protein id", "Element symbol", "Element name",
        "Type", "Subtype", "Class", "Subclass", "Method"
    ]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"[ERROR] Missing column '{col}' in AMRFinder file: {file_path}")

    out = pd.DataFrame({
        "MAG": mag,
        "Protein_ID": df["Protein id"].astype(str),
        "AMRFinder_ARG": df["Element symbol"].astype(str),
        "AMRFinder_ARG_name": df["Element name"].astype(str),
        "AMRFinder_Type": df["Type"].astype(str),
        "AMRFinder_Subtype": df["Subtype"].astype(str),
        "AMRFinder_Class": df["Class"].astype(str),
        "AMRFinder_Subclass": df["Subclass"].astype(str),
        "AMRFinder_Method": df["Method"].astype(str),
        "AMRFinder_Identity": pd.to_numeric(df.get("% Identity to reference", np.nan), errors="coerce"),
        "AMRFinder_Coverage": pd.to_numeric(df.get("% Coverage of reference", np.nan), errors="coerce"),
        "Validated_by_AMRFinder": 1
    })

    out = out.drop_duplicates(subset=["MAG", "Protein_ID", "AMRFinder_ARG"])
    return out


def read_rgi_file(file_path, keep_loose=False, min_loose_identity=80.0):
    mag = os.path.basename(file_path).replace(".txt", "")
    df = pd.read_csv(file_path, sep="\t")

    if df.empty:
        return pd.DataFrame()

    required = [
        "ORF_ID", "ARO", "Drug Class", "Resistance Mechanism",
        "AMR Gene Family", "Best_Identities", "Cut_Off", "Model_type"
    ]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"[ERROR] Missing column '{col}' in RGI file: {file_path}")

    df["Best_Identities"] = pd.to_numeric(df["Best_Identities"], errors="coerce")
    df["Percentage Length of Reference Sequence"] = pd.to_numeric(
        df.get("Percentage Length of Reference Sequence", np.nan), errors="coerce"
    )

    if keep_loose:
        df = df[
            (df["Cut_Off"].isin(["Strict", "Perfect"])) |
            ((df["Cut_Off"] == "Loose") & (df["Best_Identities"] >= min_loose_identity))
        ]
    else:
        df = df[df["Cut_Off"].isin(["Strict", "Perfect"])]

    if df.empty:
        return pd.DataFrame()

    out = pd.DataFrame({
        "MAG": mag,
        "Protein_ID": df["ORF_ID"].astype(str).str.split().str[0],   # 修正：只保留ORF ID本体
        "RGI_ARO": df["ARO"].astype(str),
        "RGI_Drug_Class": df["Drug Class"].astype(str),
        "RGI_Resistance_Mechanism": df["Resistance Mechanism"].astype(str),
        "RGI_AMR_Gene_Family": df["AMR Gene Family"].astype(str),
        "RGI_Identity": pd.to_numeric(df["Best_Identities"], errors="coerce"),
        "RGI_Coverage": pd.to_numeric(df.get("Percentage Length of Reference Sequence", np.nan), errors="coerce"),
        "RGI_Cut_Off": df["Cut_Off"].astype(str),
        "RGI_Model_type": df["Model_type"].astype(str),
        "RGI_Antibiotic": df.get("Antibiotic", pd.Series([""] * len(df))).astype(str),
        "Validated_by_RGI": 1
    })

    out = out.sort_values(
        by=["MAG", "Protein_ID", "RGI_Identity"],
        ascending=[True, True, False]
    ).drop_duplicates(subset=["MAG", "Protein_ID"], keep="first")

    return out


def load_all_amrfinder(amr_dir):
    dfs = []
    files = sorted(glob.glob(os.path.join(amr_dir, "*.tsv")))
    print(f"[INFO] Found {len(files)} AMRFinder files")
    for f in files:
        try:
            df = read_amrfinder_file(f)
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            print(f"[WARNING] Failed parsing AMRFinder file: {f}\n  {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def load_all_rgi(rgi_dir, keep_loose=False, min_loose_identity=80.0):
    dfs = []
    files = sorted(glob.glob(os.path.join(rgi_dir, "*.txt")))
    print(f"[INFO] Found {len(files)} RGI files")
    for f in files:
        try:
            df = read_rgi_file(f, keep_loose=keep_loose, min_loose_identity=min_loose_identity)
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            print(f"[WARNING] Failed parsing RGI file: {f}\n  {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def merge_amr_rgi(amr_df, rgi_df):
    if amr_df.empty and rgi_df.empty:
        return pd.DataFrame()

    if amr_df.empty:
        merged = rgi_df.copy()
    elif rgi_df.empty:
        merged = amr_df.copy()
    else:
        merged = pd.merge(amr_df, rgi_df, on=["MAG", "Protein_ID"], how="outer")

    if "Validated_by_AMRFinder" not in merged.columns:
        merged["Validated_by_AMRFinder"] = 0
    merged["Validated_by_AMRFinder"] = merged["Validated_by_AMRFinder"].fillna(0).astype(int)

    if "Validated_by_RGI" not in merged.columns:
        merged["Validated_by_RGI"] = 0
    merged["Validated_by_RGI"] = merged["Validated_by_RGI"].fillna(0).astype(int)

    merged["Validated_by_both"] = (
        (merged["Validated_by_AMRFinder"] == 1) &
        (merged["Validated_by_RGI"] == 1)
    ).astype(int)

    # 修正后的主注释字段
    merged["ARG"] = merged["AMRFinder_ARG"]
    merged.loc[
        merged["ARG"].isna() | (merged["ARG"] == ""),
        "ARG"
    ] = merged["RGI_AMR_Gene_Family"]

    merged["ARG_name"] = merged["AMRFinder_ARG_name"]
    merged.loc[
        merged["ARG_name"].isna() | (merged["ARG_name"] == ""),
        "ARG_name"
    ] = merged["RGI_AMR_Gene_Family"]

    merged["Class"] = merged["AMRFinder_Class"]
    merged.loc[
        merged["Class"].isna() | (merged["Class"] == ""),
        "Class"
    ] = merged["RGI_Drug_Class"]

    merged["Subclass"] = merged["AMRFinder_Subclass"]
    merged.loc[
        merged["Subclass"].isna() | (merged["Subclass"] == ""),
        "Subclass"
    ] = merged["RGI_Drug_Class"]

    merged["Resistance_Mechanism"] = merged["AMRFinder_Subtype"]
    merged.loc[
        merged["Resistance_Mechanism"].isna() | (merged["Resistance_Mechanism"] == ""),
        "Resistance_Mechanism"
    ] = merged["RGI_Resistance_Mechanism"]

    merged["AMR_Gene_Family"] = merged["AMRFinder_ARG"]
    merged.loc[
        merged["AMR_Gene_Family"].isna() | (merged["AMR_Gene_Family"] == ""),
        "AMR_Gene_Family"
    ] = merged["RGI_AMR_Gene_Family"]

    # 如果 RGI_AMR_Gene_Family 也空，再回退到 RGI_ARO
    merged.loc[
        (merged["ARG"].isna() | (merged["ARG"] == "")) & (~merged["RGI_ARO"].isna()),
        "ARG"
    ] = merged["RGI_ARO"]

    merged.loc[
        (merged["ARG_name"].isna() | (merged["ARG_name"] == "")) & (~merged["RGI_ARO"].isna()),
        "ARG_name"
    ] = merged["RGI_ARO"]

    def build_source(row):
        if row["Validated_by_both"] == 1:
            return "AMRFinder;RGI"
        elif row["Validated_by_AMRFinder"] == 1:
            return "AMRFinder"
        elif row["Validated_by_RGI"] == 1:
            return "RGI"
        else:
            return "Unknown"

    merged["Source"] = merged.apply(build_source, axis=1)

    merged["High_confidence"] = (
        (merged["Validated_by_AMRFinder"] == 1) |
        (
            (merged["Validated_by_RGI"] == 1) &
            (merged.get("RGI_Cut_Off", pd.Series([""] * len(merged))).isin(["Strict", "Perfect"]))
        ) |
        (merged["Validated_by_both"] == 1)
    ).astype(int)

    cols = [
        "MAG", "Protein_ID",
        "ARG", "ARG_name", "Class", "Subclass",
        "Resistance_Mechanism", "AMR_Gene_Family",
        "Validated_by_AMRFinder", "Validated_by_RGI", "Validated_by_both",
        "High_confidence", "Source",
        "AMRFinder_ARG", "AMRFinder_ARG_name",
        "AMRFinder_Type", "AMRFinder_Subtype",
        "AMRFinder_Class", "AMRFinder_Subclass",
        "AMRFinder_Method", "AMRFinder_Identity", "AMRFinder_Coverage",
        "RGI_ARO", "RGI_Drug_Class", "RGI_Resistance_Mechanism",
        "RGI_AMR_Gene_Family", "RGI_Identity", "RGI_Coverage",
        "RGI_Cut_Off", "RGI_Model_type", "RGI_Antibiotic"
    ]

    for c in cols:
        if c not in merged.columns:
            merged[c] = np.nan

    merged = merged[cols].sort_values(["MAG", "Protein_ID"])
    return merged


def main():
    parser = argparse.ArgumentParser(description="Merge AMRFinder and RGI by MAG + Protein_ID with corrected ARG/ARG_name fields.")
    parser.add_argument("--amrfinder_dir", required=True)
    parser.add_argument("--rgi_dir", required=True)
    parser.add_argument("--out_merged", required=True)
    parser.add_argument("--out_highconf", required=True)
    parser.add_argument("--keep_rgi_loose", action="store_true")
    parser.add_argument("--min_rgi_loose_identity", type=float, default=80.0)
    args = parser.parse_args()

    amr_df = load_all_amrfinder(args.amrfinder_dir)
    rgi_df = load_all_rgi(args.rgi_dir, keep_loose=args.keep_rgi_loose,
                          min_loose_identity=args.min_rgi_loose_identity)

    merged = merge_amr_rgi(amr_df, rgi_df)
    merged.to_csv(args.out_merged, sep="\t", index=False)

    highconf = merged[merged["High_confidence"] == 1].copy()
    highconf.to_csv(args.out_highconf, sep="\t", index=False)

    print(f"[INFO] Merged table written to: {args.out_merged}")
    print(f"[INFO] High-confidence table written to: {args.out_highconf}")
    print(f"[INFO] Total merged records: {merged.shape[0]}")
    print(f"[INFO] Total high-confidence records: {highconf.shape[0]}")


if __name__ == "__main__":
    main()
