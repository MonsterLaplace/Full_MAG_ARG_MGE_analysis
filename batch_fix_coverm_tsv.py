#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path


def fix_one_file(input_file: Path, output_file: Path):
    with input_file.open("r", encoding="utf-8", newline="") as fin, \
         output_file.open("w", encoding="utf-8", newline="") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        # 读取并跳过原表头
        try:
            next(reader)
        except StopIteration:
            print(f"跳过空文件: {input_file}")
            return

        # 写入统一的新表头
        writer.writerow([
            "Genome",
            "relative_abundance",
            "trimmed_mean",
            "covered_fraction",
            "mean"
        ])

        # 写入数据，跳过 unmapped
        for row in reader:
            if not row:
                continue
            if row[0] == "unmapped":
                continue
            writer.writerow(row)


def main():
    coverm_dir = Path("09_mapping/coverm")

    if not coverm_dir.exists():
        print(f"目录不存在: {coverm_dir}")
        return

    tsv_files = sorted(coverm_dir.glob("*.tsv"))

    if not tsv_files:
        print(f"目录中没有找到 .tsv 文件: {coverm_dir}")
        return

    count = 0
    for input_file in tsv_files:
        # 跳过已经修复过的文件
        if input_file.name.endswith(".fixed.tsv"):
            continue

        output_file = input_file.with_name(input_file.stem + ".fixed.tsv")
        fix_one_file(input_file, output_file)
        print(f"已生成: {output_file}")
        count += 1

    print(f"\n处理完成，共生成 {count} 个 fixed 文件。")


if __name__ == "__main__":
    main()
