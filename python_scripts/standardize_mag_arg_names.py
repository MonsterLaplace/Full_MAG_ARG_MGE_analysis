#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import re
import numpy as np


def clean_str(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def choose_best_text(row):
    """
    为标准化选择最可靠的源文本
    优先级：
    1. AMRFinder_ARG
    2. AMRFinder_ARG_name
    3. AMR_Gene_Family
    4. ARG_name
    5. ARG
    6. RGI_ARO
    """
    candidates = [
        ("AMRFinder_ARG", clean_str(row.get("AMRFinder_ARG", ""))),
        ("AMRFinder_ARG_name", clean_str(row.get("AMRFinder_ARG_name", ""))),
        ("AMR_Gene_Family", clean_str(row.get("AMR_Gene_Family", ""))),
        ("ARG_name", clean_str(row.get("ARG_name", ""))),
        ("ARG", clean_str(row.get("ARG", ""))),
        ("RGI_ARO", clean_str(row.get("RGI_ARO", "")))
    ]
    for source, text in candidates:
        if text != "" and text.lower() != "nan":
            return source, text
    return "fallback_original", ""


def extract_gene_symbol(text):
    """
    从文本中提取更标准的ARG symbol
    """
    t = text.strip()

    # 常见明确gene模式
    patterns = [
        r"(aac\([^)]+\)(?:-[A-Za-z0-9]+)?)",
        r"(aph\([^)]+\)(?:-[A-Za-z0-9]+)?)",
        r"(erm\([^)]+\)|erm[A-Za-z0-9_-]+)",
        r"(tet\([^)]+\)|tet[A-Za-z0-9_-]+)",
        r"(blaOXA|blaTEM|blaSHV|blaCTX-M|bla[A-Za-z0-9_-]*)",
        r"(van[A-Z])",
        r"(van[a-z])",
        r"(optrA)",
        r"(poxtA)",
        r"(cfr\(?[A-Za-z0-9-]*\)?)",
        r"(sul[0-9A-Za-z_-]+)",
        r"(dfr[A-Za-z0-9_-]+)",
        r"(mecA)",
        r"(macB)",
        r"(bcrA)",
        r"(novA)",
        r"(mupA|mupB)",
        r"(gyrA|gyrB|rpoB2|rpoB|fusA|fusE)",
        r"(IreK)",
        r"(cls)",
        r"(van ligase|Van ligase)",
        r"(vanH|vanR|vanS|vanT|vanW|vanY|vanU|vanL|vanX|vanZ)"
    ]

    for pat in patterns:
        m = re.search(pat, t, flags=re.IGNORECASE)
        if m:
            return m.group(1)

    # 处理 family 风格名称
    low = t.lower()

    if "aminoglycoside o-phosphotransferase" in low:
        return "aph-like"
    if "aminoglycoside 6'-n-acetyltransferase" in low or "aac(6')" in low:
        return "aac(6')-like"
    if "ribosomal protection protein" in low and "tetracycline" in low:
        return "tet-like"
    if "abc antibiotic efflux pump" in low:
        return "ABC_efflux"
    if "mfs antibiotic efflux pump" in low:
        return "MFS_efflux"
    if "rnd antibiotic efflux pump" in low:
        return "RND_efflux"
    if "smr antibiotic efflux pump" in low:
        return "SMR_efflux"
    if "mate transporter" in low:
        return "MATE_efflux"
    if "glycopeptide resistance gene cluster" in low:
        return "glycopeptide_cluster_gene"
    if "beta-lactamase" in low:
        return "beta-lactamase-like"
    if "23s ribosomal rna methyltransferase" in low:
        return "23S_rRNA_methyltransferase"

    return t


def standardize_arg_name(text):
    """
    把描述统一成更规范的 ARG_name_std
    """
    t = clean_str(text)
    low = t.lower()

    if re.search(r"van[a-z]", low):
        m = re.search(r"(van[a-z])", low)
        if m:
            gene = m.group(1)
            return f"glycopeptide resistance gene cluster; {gene}"

    if "van ligase" in low:
        return "glycopeptide resistance gene cluster; Van ligase"

    if "aph(3')" in low or "aminoglycoside o-phosphotransferase" in low:
        return "APH(3') family aminoglycoside O-phosphotransferase"

    if "aac(6')" in low or "aminoglycoside 6'-n-acetyltransferase" in low:
        return "aminoglycoside 6'-N-acetyltransferase"

    if "tetracycline-resistant ribosomal protection protein" in low or "ribosomal protection protein" in low:
        return "tetracycline-resistant ribosomal protection protein"

    if "abc antibiotic efflux pump" in low:
        return "ABC antibiotic efflux pump"

    if "mfs antibiotic efflux pump" in low:
        return "MFS antibiotic efflux pump"

    if "rnd antibiotic efflux pump" in low:
        return "RND antibiotic efflux pump"

    if "smr antibiotic efflux pump" in low:
        return "SMR antibiotic efflux pump"

    if "mate transporter" in low:
        return "MATE multidrug efflux transporter"

    if "beta-lactamase" in low:
        return t

    if "23s ribosomal rna methyltransferase" in low:
        return "23S ribosomal RNA methyltransferase"

    return t


def assign_arg_group(arg_std, name_std, row):
    """
    归为更高层 group
    """
    x = f"{clean_str(arg_std)} {clean_str(name_std)} {clean_str(row.get('Resistance_Mechanism', ''))}".lower()

    if any(k in x for k in ["vanh", "vanr", "vans", "vant", "vanw", "vany", "vanu", "van ligase", "glycopeptide resistance gene cluster"]):
        return "glycopeptide_cluster"

    if any(k in x for k in ["aac(", "aph(", "aminoglycoside o-phosphotransferase", "aminoglycoside 6'-n-acetyltransferase"]):
        return "aminoglycoside_modifying_enzyme"

    if any(k in x for k in ["tet(", "tet-like", "ribosomal protection protein", "optra", "poxta"]):
        return "ribosomal_protection"

    if any(k in x for k in ["erm", "23s ribosomal rna methyltransferase", "cfr", "llma", "emta", "tsnr"]):
        return "rRNA_methyltransferase"

    if any(k in x for k in ["beta-lactamase", "aim", "pln", "lra", "nmc", "cfxa", "bla"]):
        return "beta_lactamase"

    if any(k in x for k in ["ab c_efflux", "abc antibiotic efflux pump", "mfs antibiotic efflux pump",
                            "rnd antibiotic efflux pump", "smr antibiotic efflux pump",
                            "mate multidrug efflux transporter", "efflux"]):
        return "efflux_pump"

    if any(k in x for k in ["gyr", "rpob", "fusa", "fuse", "meca", "mura", "rpsa", "cls", "irek"]):
        return "target_alteration"

    if any(k in x for k in ["sul", "dfr", "nova", "mupa", "mupb"]):
        return "target_replacement"

    return "other"


def standardize_row(row):
    source, text = choose_best_text(row)

    arg_std = extract_gene_symbol(text)
    name_std = standardize_arg_name(text)
    group = assign_arg_group(arg_std, name_std, row)

    return pd.Series({
        "ARG_std": arg_std,
        "ARG_name_std": name_std,
        "ARG_group": group,
        "Standardization_source": source
    })


def main():
    parser = argparse.ArgumentParser(description="Standardize ARG names in mag_arg_highconf.tsv")
    parser.add_argument("-i", "--input", required=True, help="input mag_arg_highconf.tsv")
    parser.add_argument("-o", "--output", required=True, help="output standardized tsv")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    std_df = df.apply(standardize_row, axis=1)
    out = pd.concat([df, std_df], axis=1)

    out.to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Output written to: {args.output}")
    print("[INFO] ARG_group counts:")
    print(out["ARG_group"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
