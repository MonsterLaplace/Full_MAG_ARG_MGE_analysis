#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Generate contig -> raw_bin -> representative MAG mapping"
    )
    parser.add_argument("--contig_to_rawbin", required=True, help="contig_to_rawbin.tsv")
    parser.add_argument("--bin_to_repMAG", required=True, help="bin_to_repMAG.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output contig_to_repMAG.tsv")
    args = parser.parse_args()

    c2b = pd.read_csv(args.contig_to_rawbin, sep="\t")
    b2r = pd.read_csv(args.bin_to_repMAG, sep="\t")

    required_c2b = {"Sample", "Contig", "raw_bin", "ContigKey"}
    required_b2r = {"raw_bin", "repMAG", "cluster"}

    if not required_c2b.issubset(set(c2b.columns)):
        raise ValueError(f"contig_to_rawbin.tsv must contain: {required_c2b}")

    if not required_b2r.issubset(set(b2r.columns)):
        raise ValueError(f"bin_to_repMAG.tsv must contain: {required_b2r}")

    merged = pd.merge(c2b, b2r, on="raw_bin", how="left")

    missing = merged["repMAG"].isna().sum()
    if missing > 0:
        print(f"[WARNING] {missing} contigs have no matched representative MAG", file=sys.stderr)

    merged = merged[["Sample", "Contig", "raw_bin", "repMAG", "cluster", "ContigKey"]]

    merged.to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Output written to: {args.output}")
    print(f"[INFO] Total contigs: {merged.shape[0]}")
    print(f"[INFO] Total raw bins: {merged['raw_bin'].nunique()}")
    print(f"[INFO] Total repMAGs: {merged['repMAG'].nunique()}")


if __name__ == "__main__":
    main()
