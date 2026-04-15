#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd


def load_id_map(id_map_file):
    """
    读取id映射表：
    第1列 = 修改后的id
    第2列 = 原始id
    """
    id_map = {}
    with open(id_map_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"[WARNING] Skip malformed line {line_num} in {id_map_file}: {line}", file=sys.stderr)
                continue
            new_id = parts[0]
            original_id = parts[1]
            id_map[new_id] = original_id
    return id_map


def replace_tsv_first_col(tsv_file, out_file, id_map):
    """
    替换TSV第一列（通常是seqID）
    保留表头
    """
    replaced = 0
    missing = set()

    with open(tsv_file, "r", encoding="utf-8") as fin, open(out_file, "w", encoding="utf-8") as fout:
        for i, line in enumerate(fin):
            line = line.rstrip("\n")
            if i == 0:
                # 表头原样输出
                fout.write(line + "\n")
                continue

            if not line:
                fout.write("\n")
                continue

            parts = line.split("\t")
            old_id = parts[0]

            if old_id in id_map:
                parts[0] = id_map[old_id]
                replaced += 1
            else:
                missing.add(old_id)

            fout.write("\t".join(parts) + "\n")

    print(f"[INFO] TSV replaced lines: {replaced}")
    if missing:
        print(f"[WARNING] TSV IDs not found in id_map: {len(missing)}", file=sys.stderr)
        for x in sorted(list(missing))[:20]:
            print(f"  {x}", file=sys.stderr)
        if len(missing) > 20:
            print("  ...", file=sys.stderr)


def replace_gff_first_col(gff_file, out_file, id_map):
    """
    替换GFF第一列seqid
    注释行(#开头)保持不变
    """
    replaced = 0
    missing = set()

    with open(gff_file, "r", encoding="utf-8") as fin, open(out_file, "w", encoding="utf-8") as fout:
        for line in fin:
            line = line.rstrip("\n")

            if not line:
                fout.write("\n")
                continue

            if line.startswith("#"):
                fout.write(line + "\n")
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                fout.write(line + "\n")
                continue

            old_id = parts[0]
            if old_id in id_map:
                parts[0] = id_map[old_id]
                replaced += 1
            else:
                missing.add(old_id)

            fout.write("\t".join(parts) + "\n")

    print(f"[INFO] GFF replaced lines: {replaced}")
    if missing:
        print(f"[WARNING] GFF IDs not found in id_map: {len(missing)}", file=sys.stderr)
        for x in sorted(list(missing))[:20]:
            print(f"  {x}", file=sys.stderr)
        if len(missing) > 20:
            print("  ...", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Replace ISEScan-modified contig IDs in TSV and GFF using id_map.tsv"
    )
    parser.add_argument("--id_map", required=True, help="sample.id_map.tsv")
    parser.add_argument("--tsv", required=True, help="sample.clean.fa.tsv")
    parser.add_argument("--gff", required=True, help="sample.clean.fa.gff")
    parser.add_argument("--out_tsv", required=True, help="output fixed TSV")
    parser.add_argument("--out_gff", required=True, help="output fixed GFF")

    args = parser.parse_args()

    print(f"[INFO] Loading id map: {args.id_map}")
    id_map = load_id_map(args.id_map)
    print(f"[INFO] Loaded {len(id_map)} ID mappings")

    print(f"[INFO] Processing TSV: {args.tsv}")
    replace_tsv_first_col(args.tsv, args.out_tsv, id_map)

    print(f"[INFO] Processing GFF: {args.gff}")
    replace_gff_first_col(args.gff, args.out_gff, id_map)

    print("[INFO] Done.")


if __name__ == "__main__":
    main()
