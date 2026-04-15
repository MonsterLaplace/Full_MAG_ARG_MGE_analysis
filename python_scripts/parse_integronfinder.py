#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
import sys
from collections import defaultdict


def find_integrons_file(input_path):
    """
    支持输入：
    1) .integrons / .integrons.tsv 文件
    2) IntegronFinder结果目录
    自动寻找包含 '.integrons' 的文件
    """
    if os.path.isfile(input_path):
        return input_path

    candidates = []
    for root, dirs, files in os.walk(input_path):
        for f in files:
            if ".integrons" in f:
                candidates.append(os.path.join(root, f))

    if len(candidates) == 0:
        return None
    elif len(candidates) == 1:
        return candidates[0]
    else:
        print("[WARNING] Multiple integrons-like files found, using the first one:", file=sys.stderr)
        for c in candidates:
            print("  ", c, file=sys.stderr)
        return candidates[0]


def read_integrons_table(file_path):
    """
    读取 IntegronFinder 2.0.6 明细表
    自动跳过 # 注释行
    """
    try:
        df = pd.read_csv(file_path, sep="\t", comment="#", dtype=str)
    except Exception as e:
        raise RuntimeError(f"Failed to read file {file_path}: {e}")

    required_cols = [
        "ID_integron", "ID_replicon", "element", "type_elt",
        "annotation", "model", "type", "default", "considered_topology"
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    return df


def parse_integronfinder_detail(df):
    """
    针对 IntegronFinder 2.0.6 .integrons 明细表按 contig 汇总
    """

    summary = defaultdict(lambda: {
        "integron_present": 0,
        "integrase_present": 0,
        "attc_count": 0,
        "protein_count": 0,
        "complete_integron": 0,
        "in0_count": 0,
        "CALIN_count": 0,
        "integron_type_list": set(),
        "integron_id_list": set(),
        "considered_topology": set(),
        "default_yes_count": 0,
        "default_no_count": 0
    })

    for _, row in df.iterrows():
        contig = str(row["ID_replicon"]).strip()
        integron_id = str(row["ID_integron"]).strip()
        type_elt = str(row["type_elt"]).strip()
        annotation = str(row["annotation"]).strip()
        model = str(row["model"]).strip()
        int_type = str(row["type"]).strip()
        default = str(row["default"]).strip()
        topo = str(row["considered_topology"]).strip()

        summary[contig]["integron_present"] = 1
        summary[contig]["integron_id_list"].add(integron_id)
        summary[contig]["integron_type_list"].add(int_type)
        if topo and topo != "nan":
            summary[contig]["considered_topology"].add(topo)

        # 统计 default
        if default.lower() == "yes":
            summary[contig]["default_yes_count"] += 1
        elif default.lower() == "no":
            summary[contig]["default_no_count"] += 1

        # 统计 attC / protein
        if type_elt.lower() == "attc":
            summary[contig]["attc_count"] += 1
        elif type_elt.lower() == "protein":
            summary[contig]["protein_count"] += 1

        # integron 类型统计
        int_type_low = int_type.lower()
        if int_type_low == "complete":
            summary[contig]["complete_integron"] += 1
        elif int_type_low == "in0":
            summary[contig]["in0_count"] += 1
        elif int_type_low == "calin":
            summary[contig]["CALIN_count"] += 1

        # integrase 判断
        # 你的文件中没有单独 integrase 列，所以只能从 annotation / model / element 推断
        # 若 annotation/model/element 中出现 integrase 或 intl 则记为1
        text = " ".join([
            str(row.get("element", "")),
            annotation,
            model
        ]).lower()

        if ("integrase" in text) or ("intl" in text) or ("inti1" in text) or ("int1" in text):
            summary[contig]["integrase_present"] = 1

    rows = []
    for contig, info in summary.items():
        rows.append({
            "Contig": contig,
            "integron_present": info["integron_present"],
            "integrase_present": info["integrase_present"],
            "attc_count": info["attc_count"],
            "protein_count": info["protein_count"],
            "complete_integron": info["complete_integron"],
            "in0_count": info["in0_count"],
            "CALIN_count": info["CALIN_count"],
            "integron_type_list": ";".join(sorted(info["integron_type_list"])) if info["integron_type_list"] else "",
            "integron_id_list": ";".join(sorted(info["integron_id_list"])) if info["integron_id_list"] else "",
            "considered_topology": ";".join(sorted(info["considered_topology"])) if info["considered_topology"] else "",
            "default_yes_count": info["default_yes_count"],
            "default_no_count": info["default_no_count"]
        })

    out_df = pd.DataFrame(rows)

    # 保证类型正确
    int_cols = [
        "integron_present", "integrase_present", "attc_count", "protein_count",
        "complete_integron", "in0_count", "CALIN_count",
        "default_yes_count", "default_no_count"
    ]
    for c in int_cols:
        out_df[c] = out_df[c].astype(int)

    return out_df.sort_values("Contig")


def main():
    parser = argparse.ArgumentParser(
        description="Parse IntegronFinder 2.0.6 detailed .integrons table into contig-level summary."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="IntegronFinder result directory or .integrons/.integrons.tsv file"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output parsed TSV file"
    )
    args = parser.parse_args()

    integrons_file = find_integrons_file(args.input)
    if integrons_file is None:
        print("[ERROR] No .integrons file found.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Parsing file: {integrons_file}", file=sys.stderr)

    df = read_integrons_table(integrons_file)
    print(f"[INFO] Input shape: {df.shape}", file=sys.stderr)
    print(f"[INFO] Columns: {list(df.columns)}", file=sys.stderr)

    out_df = parse_integronfinder_detail(df)
    out_df.to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Done. Output written to: {args.output}", file=sys.stderr)
    print(f"[INFO] Total contigs with integron signal: {out_df.shape[0]}", file=sys.stderr)


if __name__ == "__main__":
    main()
