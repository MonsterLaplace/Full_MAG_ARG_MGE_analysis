#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import pandas as pd
from collections import Counter, defaultdict

############################################
# 1. Config
############################################

# 修改成你的真实目录
INPUT_DIR = "/data/xb/MAG/10_annotation/cazy"
OUT_DIR = "/data/xb/MAG/10_annotation/cazy_summary"

os.makedirs(OUT_DIR, exist_ok=True)

OUT_GENE = os.path.join(OUT_DIR, "mag_cazy_gene_level.tsv")
OUT_LONG = os.path.join(OUT_DIR, "mag_cazy_long.tsv")
OUT_WIDE = os.path.join(OUT_DIR, "mag_cazy_catalog.tsv")
OUT_STATUS = os.path.join(OUT_DIR, "mag_cazy_file_status.tsv")

############################################
# 2. Regex
############################################

MAIN_FAMILY_PATTERN = re.compile(r'\b(GH\d+|GT\d+|CE\d+|CBM\d+|PL\d+|AA\d+)\b')
EXTENDED_PATTERN = re.compile(
    r'\b(GH\d+(?:_[A-Za-z0-9]+)?|GT\d+(?:_[A-Za-z0-9]+)?|CE\d+(?:_[A-Za-z0-9]+)?|CBM\d+(?:_[A-Za-z0-9]+)?|PL\d+(?:_[A-Za-z0-9]+)?|AA\d+(?:_[A-Za-z0-9]+)?)\b'
)

############################################
# 3. Helper functions
############################################

def classify_family(fam):
    if fam.startswith("GH"):
        return "GH"
    elif fam.startswith("GT"):
        return "GT"
    elif fam.startswith("CE"):
        return "CE"
    elif fam.startswith("CBM"):
        return "CBM"
    elif fam.startswith("PL"):
        return "PL"
    elif fam.startswith("AA"):
        return "AA"
    else:
        return "Other"


def normalize_family(token):
    m = MAIN_FAMILY_PATTERN.search(token)
    if m:
        return m.group(1)
    return None


def extract_families_from_text(text):
    if pd.isna(text):
        return []

    text = str(text).strip()
    if text in ["-", "", "NA", "nan", "None"]:
        return []

    # 去掉坐标区间
    text2 = re.sub(r'\([^)]*\)', '', text)

    # 以 + 和 | 分隔
    parts = re.split(r'[+|]', text2)

    fams = []
    for p in parts:
        p = p.strip()
        if not p:
            continue

        ext_hits = EXTENDED_PATTERN.findall(p)
        if len(ext_hits) > 0:
            for h in ext_hits:
                fam = normalize_family(h)
                if fam:
                    fams.append(fam)
        else:
            fam = normalize_family(p)
            if fam:
                fams.append(fam)

    return fams


def parse_overview_file(filepath, mag_name):
    """
    Parse a dbCAN overview file.
    """
    try:
        df = pd.read_csv(filepath, sep="\t", dtype=str)
    except Exception as e:
        print(f"[WARN] Cannot parse file: {filepath}; error={e}")
        return []

    if df.shape[0] == 0:
        return []

    # Gene ID 列
    gene_col = "Gene ID" if "Gene ID" in df.columns else df.columns[0]

    parse_cols = [c for c in ["dbCAN_hmm", "dbCAN_sub", "DIAMOND", "Recommend Results"] if c in df.columns]
    if len(parse_cols) == 0:
        print(f"[WARN] No CAZy columns found in: {filepath}")
        return []

    records = []

    for _, row in df.iterrows():
        gene_id = str(row[gene_col]).strip()
        fams = set()

        for c in parse_cols:
            fams.update(extract_families_from_text(row[c]))

        for fam in sorted(fams):
            records.append({
                "MAG": mag_name,
                "Gene_ID": gene_id,
                "CAZy_family": fam,
                "Class": classify_family(fam)
            })

    return records


def find_overview_files(input_dir):
    """
    Recursively scan input directory for overview.tsv or overview.txt.
    """
    file_map = []
    status_records = []

    # 找所有子目录中的 overview.tsv / overview.txt
    for root, dirs, files in os.walk(input_dir):
        for fname in files:
            if fname in ["overview.tsv", "overview.txt"]:
                fp = os.path.join(root, fname)
                mag_name = os.path.basename(root)
                file_map.append((mag_name, fp))
                status_records.append({
                    "MAG": mag_name,
                    "overview_file": fp,
                    "status": "found"
                })

    # 去重
    uniq = {}
    for mag, fp in file_map:
        uniq[(mag, fp)] = 1
    file_map = list(uniq.keys())

    return file_map, status_records


############################################
# 4. Main
############################################

def main():
    file_map, status_records = find_overview_files(INPUT_DIR)
    print(f"[INFO] Found {len(file_map)} overview files under {INPUT_DIR}")

    all_records = []

    for mag_name, fp in file_map:
        recs = parse_overview_file(fp, mag_name)
        all_records.extend(recs)

    status_df = pd.DataFrame(status_records).drop_duplicates()
    if status_df.shape[0] > 0:
        status_df.to_csv(OUT_STATUS, sep="\t", index=False)

    if len(all_records) == 0:
        print("[ERROR] No CAZy records extracted. Please check directory structure or overview file format.")
        return

    # gene-level
    gene_df = pd.DataFrame(all_records).drop_duplicates()
    gene_df.to_csv(OUT_GENE, sep="\t", index=False)

    # long table: per MAG per family
    long_df = (
        gene_df.groupby(["MAG", "CAZy_family", "Class"])
        .size()
        .reset_index(name="Count")
    )
    long_df.to_csv(OUT_LONG, sep="\t", index=False)

    # wide table: per MAG summary
    mag_summary = []
    for mag, sub in long_df.groupby("MAG"):
        total_count = int(sub["Count"].sum())
        fam_count = int(sub["CAZy_family"].nunique())

        class_counter = Counter()
        class_fam_counter = defaultdict(set)

        for _, row in sub.iterrows():
            class_counter[row["Class"]] += int(row["Count"])
            class_fam_counter[row["Class"]].add(row["CAZy_family"])

        fam_list = sorted(sub["CAZy_family"].unique().tolist())

        mag_summary.append({
            "MAG": mag,
            "CAZy_count": total_count,
            "CAZy_family_count": fam_count,
            "GH_count": int(class_counter.get("GH", 0)),
            "GT_count": int(class_counter.get("GT", 0)),
            "CE_count": int(class_counter.get("CE", 0)),
            "CBM_count": int(class_counter.get("CBM", 0)),
            "PL_count": int(class_counter.get("PL", 0)),
            "AA_count": int(class_counter.get("AA", 0)),
            "GH_family_count": len(class_fam_counter.get("GH", set())),
            "GT_family_count": len(class_fam_counter.get("GT", set())),
            "CE_family_count": len(class_fam_counter.get("CE", set())),
            "CBM_family_count": len(class_fam_counter.get("CBM", set())),
            "PL_family_count": len(class_fam_counter.get("PL", set())),
            "AA_family_count": len(class_fam_counter.get("AA", set())),
            "CAZy_family_list": ";".join(fam_list)
        })

    wide_df = pd.DataFrame(mag_summary)
    wide_df.to_csv(OUT_WIDE, sep="\t", index=False)

    print("[INFO] Finished successfully.")
    print(f"[INFO] Gene-level table : {OUT_GENE}")
    print(f"[INFO] Long table       : {OUT_LONG}")
    print(f"[INFO] Wide table       : {OUT_WIDE}")
    print(f"[INFO] Status table     : {OUT_STATUS}")


if __name__ == "__main__":
    main()
