#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import glob
import os


def main():
    parser = argparse.ArgumentParser(description="Merge all sample contig_mge tables.")
    parser.add_argument("-i", "--input_dir", required=True, help="directory of sample *.contig_mge.tsv")
    parser.add_argument("-o", "--output", required=True, help="merged output file")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.input_dir, "*.contig_mge.tsv"))
    dfs = []

    for f in files:
        sample = os.path.basename(f).replace(".contig_mge.tsv", "")
        df = pd.read_csv(f, sep="\t")
        df["Sample"] = sample
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Done. Output: {args.output}")


if __name__ == "__main__":
    main()
