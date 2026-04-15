#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
import numpy as np


def safe_read_tsv(path, required=True):
    if not os.path.exists(path):
        if required:
            raise FileNotFoundError(f"[ERROR] File not found: {path}")
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def ensure_columns(df, required_cols, df_name):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"[ERROR] Missing columns in {df_name}: {missing}")


def unique_join(series):
    vals = [str(x) for x in series if pd.notnull(x) and str(x) != "" and str(x) != "nan"]
    vals = sorted(set(vals))
    return ";".join(vals)


def summarize_mge_types(subdf):
    tags = []
    if "plasmid_like" in subdf.columns and (subdf["plasmid_like"] == 1).any():
        tags.append("plasmid_like")
    if "virus_like" in subdf.columns and (subdf["virus_like"] == 1).any():
        tags.append("virus_like")
    if "integron_present" in subdf.columns and (subdf["integron_present"] == 1).any():
        tags.append("integron")
    if "integrase_present" in subdf.columns and (subdf["integrase_present"] == 1).any():
        tags.append("integrase")
    if "IS_present" in subdf.columns and (subdf["IS_present"] == 1).any():
        tags.append("IS")
    if "transposase_count" in subdf.columns and (pd.to_numeric(subdf["transposase_count"], errors="coerce").fillna(0) > 0).any():
        tags.append("transposase")
    if "TIR_present" in subdf.columns and (subdf["TIR_present"] == 1).any():
        tags.append("TIR")
    return ";".join(tags) if tags else "none"


def ensure_optional_columns(df, cols):
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    return df


def load_and_prepare(contig_to_rep, contig_arg_mge, mag_arg_highconf):
    c2r = safe_read_tsv(contig_to_rep, required=True)
    cargm = safe_read_tsv(contig_arg_mge, required=True)
    magarg = safe_read_tsv(mag_arg_highconf, required=True)

    ensure_columns(c2r, ["Sample", "Contig", "raw_bin"], "contig_to_repMAG.tsv")
    ensure_columns(cargm, ["Sample", "Contig", "protein", "ARG"], "contig_arg_mge.tsv")
    ensure_columns(magarg, ["MAG", "Protein_ID", "ARG"], "mag_arg_highconf.tsv")

    if "repMAG" not in c2r.columns:
        c2r["repMAG"] = np.nan
    if "cluster" not in c2r.columns:
        c2r["cluster"] = np.nan

    for df in [c2r, cargm]:
        df["Sample"] = df["Sample"].astype(str)
        df["Contig"] = df["Contig"].astype(str)

    c2r["raw_bin"] = c2r["raw_bin"].astype(str)
    c2r["repMAG"] = c2r["repMAG"].astype(str)
    c2r["cluster"] = c2r["cluster"].astype(str)

    for col in ["repMAG", "cluster"]:
        c2r[col] = c2r[col].replace(["", "nan", "None", "NA", "NaN"], np.nan)

    cargm["protein"] = cargm["protein"].astype(str)

    magarg["MAG"] = magarg["MAG"].astype(str)
    magarg["Protein_ID"] = magarg["Protein_ID"].astype(str)

    merged = pd.merge(cargm, c2r, on=["Sample", "Contig"], how="left")

    merged["mapped_to_repMAG"] = (~merged["repMAG"].isna()).astype(int)

    mapped = merged[merged["mapped_to_repMAG"] == 1].copy()
    unmapped = merged[merged["mapped_to_repMAG"] == 0].copy()

    return c2r, cargm, merged, mapped, unmapped, magarg


def build_mapping_summary(merged_df):
    summary = {
        "total_ARG_records": merged_df.shape[0],
        "mapped_ARG_records": int((merged_df["mapped_to_repMAG"] == 1).sum()),
        "unmapped_ARG_records": int((merged_df["mapped_to_repMAG"] == 0).sum()),
        "mapped_ARG_contigs": merged_df.loc[merged_df["mapped_to_repMAG"] == 1, ["Sample", "Contig"]].drop_duplicates().shape[0],
        "unmapped_ARG_contigs": merged_df.loc[merged_df["mapped_to_repMAG"] == 0, ["Sample", "Contig"]].drop_duplicates().shape[0],
        "mapped_raw_bins": merged_df.loc[merged_df["mapped_to_repMAG"] == 1, "raw_bin"].dropna().astype(str).replace("nan", np.nan).dropna().nunique(),
        "unmapped_raw_bins": merged_df.loc[merged_df["mapped_to_repMAG"] == 0, "raw_bin"].dropna().astype(str).replace("nan", np.nan).dropna().nunique(),
    }
    return pd.DataFrame([summary])


def build_mag_arg_catalog(mapped_df, mag_arg_highconf):
    cols_needed = [
        "repMAG", "cluster", "raw_bin", "Sample", "Contig", "protein",
        "ARG", "ARG_name", "Class", "Subclass",
        "Resistance_Mechanism", "AMR_Gene_Family",
        "ARG_std", "ARG_name_std", "ARG_group", "Standardization_source",
        "Validated_by_AMRFinder", "Validated_by_RGI", "Validated_by_both",
        "High_confidence", "Source",
        "putative_mobile_ARG_loose", "putative_mobile_ARG_strict"
    ]
    for c in cols_needed:
        if c not in mapped_df.columns:
            mapped_df[c] = np.nan

    mag_arg_catalog = mapped_df[cols_needed].copy()
    mag_arg_catalog = mag_arg_catalog.rename(columns={
        "repMAG": "MAG",
        "protein": "Protein_ID"
    })

    mag_arg_catalog = pd.merge(
        mag_arg_catalog,
        mag_arg_highconf[[
            "MAG", "Protein_ID",
            "ARG", "ARG_name", "Class", "Subclass",
            "Resistance_Mechanism", "AMR_Gene_Family",
            "Validated_by_AMRFinder", "Validated_by_RGI",
            "Validated_by_both", "High_confidence", "Source"
        ]].drop_duplicates(),
        on=["MAG", "Protein_ID"],
        how="left",
        suffixes=("", "_MAG")
    )

    for col in [
        "ARG", "ARG_name", "Class", "Subclass",
        "Resistance_Mechanism", "AMR_Gene_Family",
        "Validated_by_AMRFinder", "Validated_by_RGI",
        "Validated_by_both", "High_confidence", "Source"
    ]:
        alt = f"{col}_MAG"
        if alt in mag_arg_catalog.columns:
            mag_arg_catalog[col] = mag_arg_catalog[col].fillna(mag_arg_catalog[alt])

    mag_arg_catalog = mag_arg_catalog[[c for c in mag_arg_catalog.columns if not c.endswith("_MAG")]]

    mag_arg_catalog = mag_arg_catalog.drop_duplicates(
        subset=["MAG", "Protein_ID", "ARG", "ARG_std", "ARG_group"]
    )

    return mag_arg_catalog.sort_values(["MAG", "Protein_ID"])


def build_mag_mge_catalog(mapped_df):
    rows = []

    optional_cols = [
        "ARG_std", "ARG_group", "Class", "Subclass",
        "Resistance_Mechanism", "AMR_Gene_Family",
        "High_confidence", "Validated_by_AMRFinder",
        "Validated_by_RGI", "Validated_by_both",
        "putative_mobile_ARG_loose", "putative_mobile_ARG_strict",
        "plasmid_like", "virus_like", "chromosome_like",
        "integron_present", "integrase_present", "attc_count",
        "complete_integron", "in0_count",
        "IS_present", "is_element_count", "transposase_count", "TIR_present"
    ]
    mapped_df = ensure_optional_columns(mapped_df, optional_cols)

    for mag, subdf in mapped_df.groupby("repMAG"):
        row = {
            "MAG": mag,
            "cluster": unique_join(subdf["cluster"]),
            "raw_bin_count": subdf["raw_bin"].dropna().astype(str).replace("nan", np.nan).dropna().nunique(),
            "raw_bin_list": unique_join(subdf["raw_bin"]),

            "n_arg_records": subdf.shape[0],
            "n_contigs": subdf["Contig"].nunique(),
            "n_samples": subdf["Sample"].nunique(),

            "arg_contig_count": subdf.loc[subdf["ARG"].astype(str).replace("nan", "") != "", "Contig"].nunique(),
            "ARG_count": subdf["ARG"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "ARG_list": unique_join(subdf["ARG"]),

            "ARG_std_count": subdf["ARG_std"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "ARG_std_list": unique_join(subdf["ARG_std"]),
            "ARG_group_count": subdf["ARG_group"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "ARG_group_list": unique_join(subdf["ARG_group"]),

            "Class_count": subdf["Class"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "Class_list": unique_join(subdf["Class"]),
            "Subclass_count": subdf["Subclass"].astype(str).replace("nan", "").loc[lambda x: x != ""].nunique(),
            "Subclass_list": unique_join(subdf["Subclass"]),
            "Resistance_Mechanism_list": unique_join(subdf["Resistance_Mechanism"]),
            "AMR_Gene_Family_list": unique_join(subdf["AMR_Gene_Family"]),

            "high_conf_ARG_count": int((subdf["High_confidence"] == 1).sum()),
            "validated_by_amrfinder_count": int((subdf["Validated_by_AMRFinder"] == 1).sum()),
            "validated_by_rgi_count": int((subdf["Validated_by_RGI"] == 1).sum()),
            "validated_by_both_count": int((subdf["Validated_by_both"] == 1).sum()),

            "mobile_ARG_loose_count": int((subdf["putative_mobile_ARG_loose"] == 1).sum()),
            "mobile_ARG_strict_count": int((subdf["putative_mobile_ARG_strict"] == 1).sum()),
            "mobile_ARG_loose_contig_count": subdf.loc[subdf["putative_mobile_ARG_loose"] == 1, "Contig"].nunique(),
            "mobile_ARG_strict_contig_count": subdf.loc[subdf["putative_mobile_ARG_strict"] == 1, "Contig"].nunique(),

            "plasmid_like_contig_count": subdf.loc[subdf["plasmid_like"] == 1, "Contig"].nunique(),
            "virus_like_contig_count": subdf.loc[subdf["virus_like"] == 1, "Contig"].nunique(),
            "chromosome_like_contig_count": subdf.loc[subdf["chromosome_like"] == 1, "Contig"].nunique(),

            "integron_contig_count": subdf.loc[subdf["integron_present"] == 1, "Contig"].nunique(),
            "integrase_contig_count": subdf.loc[subdf["integrase_present"] == 1, "Contig"].nunique(),
            "attc_total": int(pd.to_numeric(subdf["attc_count"], errors="coerce").fillna(0).sum()),
            "complete_integron_count": int(pd.to_numeric(subdf["complete_integron"], errors="coerce").fillna(0).sum()),
            "in0_count": int(pd.to_numeric(subdf["in0_count"], errors="coerce").fillna(0).sum()),

            "IS_contig_count": subdf.loc[subdf["IS_present"] == 1, "Contig"].nunique(),
            "is_element_total": int(pd.to_numeric(subdf["is_element_count"], errors="coerce").fillna(0).sum()),
            "transposase_total": int(pd.to_numeric(subdf["transposase_count"], errors="coerce").fillna(0).sum()),
            "tir_contig_count": subdf.loc[subdf["TIR_present"] == 1, "Contig"].nunique(),

            "MGE_summary": summarize_mge_types(subdf)
        }

        rows.append(row)

    out = pd.DataFrame(rows)
    return out.sort_values("MAG")


def build_arg_mge_mag_edges(mapped_df):
    edge_cols = [
        "repMAG", "cluster", "raw_bin", "Sample", "Contig", "protein",
        "ARG", "ARG_std", "ARG_group",
        "Class", "Subclass", "Resistance_Mechanism",
        "High_confidence",
        "putative_mobile_ARG_loose", "putative_mobile_ARG_strict",
        "mobile_evidence_loose", "mobile_evidence_strict",
        "plasmid_like", "virus_like", "integron_present", "integrase_present",
        "IS_present", "transposase_count", "TIR_present",
        "mge_summary"
    ]
    for c in edge_cols:
        if c not in mapped_df.columns:
            mapped_df[c] = np.nan

    edges = mapped_df[edge_cols].copy()
    edges = edges.rename(columns={
        "repMAG": "MAG",
        "protein": "Protein_ID"
    })

    edges = edges.drop_duplicates()
    return edges.sort_values(["MAG", "Contig", "Protein_ID", "ARG"])


def build_mag_trait_catalog(mag_mge_catalog):
    trait = mag_mge_catalog.copy()

    trait["has_ARG"] = (trait["ARG_count"] > 0).astype(int)
    trait["has_mobile_ARG_loose"] = (trait["mobile_ARG_loose_count"] > 0).astype(int)
    trait["has_mobile_ARG_strict"] = (trait["mobile_ARG_strict_count"] > 0).astype(int)
    trait["has_integron"] = (trait["integron_contig_count"] > 0).astype(int)
    trait["has_plasmid_like_contig"] = (trait["plasmid_like_contig_count"] > 0).astype(int)
    trait["has_virus_like_contig"] = (trait["virus_like_contig_count"] > 0).astype(int)
    trait["has_IS"] = (trait["IS_contig_count"] > 0).astype(int)

    trait["mobile_ARG_loose_ratio"] = trait["mobile_ARG_loose_count"] / trait["n_arg_records"].replace(0, np.nan)
    trait["mobile_ARG_strict_ratio"] = trait["mobile_ARG_strict_count"] / trait["n_arg_records"].replace(0, np.nan)
    trait["arg_contig_ratio"] = trait["arg_contig_count"] / trait["n_contigs"].replace(0, np.nan)

    return trait.sort_values("MAG")


def build_all_mag_trait_catalog(c2r, cargm):
    all_contigs = c2r.copy()
    all_contigs = all_contigs[~all_contigs["repMAG"].isna()].copy()

    needed_cols = [
        "Sample", "Contig", "ARG", "ARG_std", "ARG_group",
        "Class", "Subclass", "Resistance_Mechanism", "AMR_Gene_Family",
        "High_confidence", "Validated_by_AMRFinder", "Validated_by_RGI", "Validated_by_both",
        "putative_mobile_ARG_loose", "putative_mobile_ARG_strict",
        "plasmid_like", "virus_like", "chromosome_like",
        "integron_present", "integrase_present", "attc_count",
        "complete_integron", "in0_count",
        "IS_present", "is_element_count", "transposase_count", "TIR_present"
    ]
    cargm2 = cargm.copy()
    cargm2 = ensure_optional_columns(cargm2, needed_cols)

    merged = pd.merge(
        all_contigs,
        cargm2[needed_cols].copy(),
        on=["Sample", "Contig"],
        how="left"
    )

    rows = []
    for mag, subdf in merged.groupby("repMAG", dropna=True):
        arg_mask = subdf["ARG"].notna() & (subdf["ARG"].astype(str).replace("nan", "") != "")
        arg_std_mask = subdf["ARG_std"].notna() & (subdf["ARG_std"].astype(str).replace("nan", "") != "")
        arg_group_mask = subdf["ARG_group"].notna() & (subdf["ARG_group"].astype(str).replace("nan", "") != "")
        class_mask = subdf["Class"].notna() & (subdf["Class"].astype(str).replace("nan", "") != "")
        subclass_mask = subdf["Subclass"].notna() & (subdf["Subclass"].astype(str).replace("nan", "") != "")

        row = {
            "MAG": mag,
            "cluster": unique_join(subdf["cluster"]),
            "raw_bin_count": subdf["raw_bin"].dropna().astype(str).replace("nan", np.nan).dropna().nunique(),
            "raw_bin_list": unique_join(subdf["raw_bin"]),

            "n_contigs": subdf["Contig"].nunique(),
            "n_samples": subdf["Sample"].nunique(),

            "arg_contig_count": subdf.loc[arg_mask, "Contig"].nunique(),
            "ARG_count": subdf.loc[arg_mask, "ARG"].astype(str).nunique(),
            "ARG_list": unique_join(subdf.loc[arg_mask, "ARG"]),

            "ARG_std_count": subdf.loc[arg_std_mask, "ARG_std"].astype(str).nunique(),
            "ARG_std_list": unique_join(subdf["ARG_std"]),
            "ARG_group_count": subdf.loc[arg_group_mask, "ARG_group"].astype(str).nunique(),
            "ARG_group_list": unique_join(subdf["ARG_group"]),

            "Class_count": subdf.loc[class_mask, "Class"].astype(str).nunique(),
            "Class_list": unique_join(subdf["Class"]),
            "Subclass_count": subdf.loc[subclass_mask, "Subclass"].astype(str).nunique(),
            "Subclass_list": unique_join(subdf["Subclass"]),
            "Resistance_Mechanism_list": unique_join(subdf["Resistance_Mechanism"]),
            "AMR_Gene_Family_list": unique_join(subdf["AMR_Gene_Family"]),

            "high_conf_ARG_count": int((subdf["High_confidence"] == 1).sum()),
            "validated_by_amrfinder_count": int((subdf["Validated_by_AMRFinder"] == 1).sum()),
            "validated_by_rgi_count": int((subdf["Validated_by_RGI"] == 1).sum()),
            "validated_by_both_count": int((subdf["Validated_by_both"] == 1).sum()),

            "mobile_ARG_loose_count": int((subdf["putative_mobile_ARG_loose"] == 1).sum()),
            "mobile_ARG_strict_count": int((subdf["putative_mobile_ARG_strict"] == 1).sum()),
            "mobile_ARG_loose_contig_count": subdf.loc[subdf["putative_mobile_ARG_loose"] == 1, "Contig"].nunique(),
            "mobile_ARG_strict_contig_count": subdf.loc[subdf["putative_mobile_ARG_strict"] == 1, "Contig"].nunique(),

            "plasmid_like_contig_count": subdf.loc[subdf["plasmid_like"] == 1, "Contig"].nunique(),
            "virus_like_contig_count": subdf.loc[subdf["virus_like"] == 1, "Contig"].nunique(),
            "chromosome_like_contig_count": subdf.loc[subdf["chromosome_like"] == 1, "Contig"].nunique(),

            "integron_contig_count": subdf.loc[subdf["integron_present"] == 1, "Contig"].nunique(),
            "integrase_contig_count": subdf.loc[subdf["integrase_present"] == 1, "Contig"].nunique(),
            "attc_total": int(pd.to_numeric(subdf["attc_count"], errors="coerce").fillna(0).sum()),
            "complete_integron_count": int(pd.to_numeric(subdf["complete_integron"], errors="coerce").fillna(0).sum()),
            "in0_count": int(pd.to_numeric(subdf["in0_count"], errors="coerce").fillna(0).sum()),

            "IS_contig_count": subdf.loc[subdf["IS_present"] == 1, "Contig"].nunique(),
            "is_element_total": int(pd.to_numeric(subdf["is_element_count"], errors="coerce").fillna(0).sum()),
            "transposase_total": int(pd.to_numeric(subdf["transposase_count"], errors="coerce").fillna(0).sum()),
            "tir_contig_count": subdf.loc[subdf["TIR_present"] == 1, "Contig"].nunique(),

            "MGE_summary": summarize_mge_types(subdf)
        }

        rows.append(row)

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    out["has_ARG"] = (out["ARG_count"] > 0).astype(int)
    out["has_mobile_ARG_loose"] = (out["mobile_ARG_loose_count"] > 0).astype(int)
    out["has_mobile_ARG_strict"] = (out["mobile_ARG_strict_count"] > 0).astype(int)
    out["has_integron"] = (out["integron_contig_count"] > 0).astype(int)
    out["has_plasmid_like_contig"] = (out["plasmid_like_contig_count"] > 0).astype(int)
    out["has_virus_like_contig"] = (out["virus_like_contig_count"] > 0).astype(int)
    out["has_IS"] = (out["IS_contig_count"] > 0).astype(int)

    out["mobile_ARG_loose_ratio"] = out["mobile_ARG_loose_contig_count"] / out["n_contigs"].replace(0, np.nan)
    out["mobile_ARG_strict_ratio"] = out["mobile_ARG_strict_contig_count"] / out["n_contigs"].replace(0, np.nan)
    out["arg_contig_ratio"] = out["arg_contig_count"] / out["n_contigs"].replace(0, np.nan)

    return out.sort_values("MAG")


def main():
    parser = argparse.ArgumentParser(
        description="Build MAG-level ARG/MGE catalogs with support for unmapped contigs and all-MAG trait table."
    )
    parser.add_argument("--contig_to_repMAG", required=True, help="13_tables/contig_to_repMAG.tsv")
    parser.add_argument("--contig_arg_mge", required=True, help="13_tables/contig_arg_mge.tsv")
    parser.add_argument("--mag_arg_highconf", required=True, help="13_tables/mag_arg_highconf.tsv")

    parser.add_argument("--out_mag_arg_catalog", required=True, help="output mag_arg_catalog.tsv")
    parser.add_argument("--out_mag_mge_catalog", required=True, help="output mag_mge_catalog.tsv")
    parser.add_argument("--out_arg_mge_mag_edges", required=True, help="output arg_mge_mag_edges.tsv")
    parser.add_argument("--out_mag_trait_catalog", required=True, help="output mag_trait_catalog.tsv")
    parser.add_argument("--out_all_mag_trait_catalog", required=True, help="output all.mag_trait_catalog.tsv")

    parser.add_argument("--out_unmapped_contig_arg_mge", required=True, help="output unmapped_contig_arg_mge.tsv")
    parser.add_argument("--out_mapping_summary", required=True, help="output mapping_summary.tsv")
    args = parser.parse_args()

    c2r, cargm, merged_all, mapped, unmapped, magarg = load_and_prepare(
        args.contig_to_repMAG,
        args.contig_arg_mge,
        args.mag_arg_highconf
    )

    unmapped.to_csv(args.out_unmapped_contig_arg_mge, sep="\t", index=False)

    mapping_summary = build_mapping_summary(merged_all)
    mapping_summary.to_csv(args.out_mapping_summary, sep="\t", index=False)

    mag_arg_catalog = build_mag_arg_catalog(mapped, magarg)
    mag_mge_catalog = build_mag_mge_catalog(mapped)
    arg_mge_mag_edges = build_arg_mge_mag_edges(mapped)
    mag_trait_catalog = build_mag_trait_catalog(mag_mge_catalog)
    all_mag_trait_catalog = build_all_mag_trait_catalog(c2r, cargm)

    mag_arg_catalog.to_csv(args.out_mag_arg_catalog, sep="\t", index=False)
    mag_mge_catalog.to_csv(args.out_mag_mge_catalog, sep="\t", index=False)
    arg_mge_mag_edges.to_csv(args.out_arg_mge_mag_edges, sep="\t", index=False)
    mag_trait_catalog.to_csv(args.out_mag_trait_catalog, sep="\t", index=False)
    all_mag_trait_catalog.to_csv(args.out_all_mag_trait_catalog, sep="\t", index=False)

    print(f"[INFO] mag_arg_catalog written to: {args.out_mag_arg_catalog}")
    print(f"[INFO] mag_mge_catalog written to: {args.out_mag_mge_catalog}")
    print(f"[INFO] arg_mge_mag_edges written to: {args.out_arg_mge_mag_edges}")
    print(f"[INFO] mag_trait_catalog written to: {args.out_mag_trait_catalog}")
    print(f"[INFO] all.mag_trait_catalog written to: {args.out_all_mag_trait_catalog}")
    print(f"[INFO] unmapped_contig_arg_mge written to: {args.out_unmapped_contig_arg_mge}")
    print(f"[INFO] mapping_summary written to: {args.out_mapping_summary}")

    print(f"[INFO] mapped ARG records: {mapped.shape[0]}")
    print(f"[INFO] unmapped ARG records: {unmapped.shape[0]}")
    print(f"[INFO] MAGs retained in mapped trait table: {mag_mge_catalog.shape[0]}")
    print(f"[INFO] MAGs retained in all trait table: {all_mag_trait_catalog.shape[0]}")


if __name__ == "__main__":
    main()
