#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import re
import sys


def strip_fasta_suffix(x):
    x = str(x)
    x = re.sub(r'\.(fa|fasta|fna|fas)$', '', x)
    return x


def detect_cdb_columns(df):
    genome_col = None
    cluster_col = None

    for c in df.columns:
        cl = c.lower()
        if cl == "genome":
            genome_col = c
        elif "secondary_cluster" in cl:
            cluster_col = c

    if cluster_col is None:
        for c in df.columns:
            cl = c.lower()
            if cl == "cluster":
                cluster_col = c
                break

    return genome_col, cluster_col


def detect_wdb_columns(df):
    genome_col = None
    cluster_col = None

    for c in df.columns:
        cl = c.lower()
        if cl == "genome":
            genome_col = c
        elif "cluster" in cl:
            cluster_col = c

    return genome_col, cluster_col


def main():
    parser = argparse.ArgumentParser(
        description="Generate raw_bin -> representative MAG mapping from dRep Cdb.csv and Wdb.csv"
    )
    parser.add_argument("--cdb", required=True, help="Path to dRep Cdb.csv")
    parser.add_argument("--wdb", required=True, help="Path to dRep Wdb.csv")
    parser.add_argument("-o", "--output", required=True, help="Output bin_to_repMAG.tsv")
    args = parser.parse_args()

    cdb = pd.read_csv(args.cdb)
    wdb = pd.read_csv(args.wdb)

    c_genome, c_cluster = detect_cdb_columns(cdb)
    w_genome, w_cluster = detect_wdb_columns(wdb)

    if c_genome is None or c_cluster is None:
        print("[ERROR] Cannot detect genome/cluster columns in Cdb.csv", file=sys.stderr)
        print("[INFO] Cdb columns:", list(cdb.columns), file=sys.stderr)
        sys.exit(1)

    if w_genome is None or w_cluster is None:
        print("[ERROR] Cannot detect genome/cluster columns in Wdb.csv", file=sys.stderr)
        print("[INFO] Wdb columns:", list(wdb.columns), file=sys.stderr)
        sys.exit(1)

    cdb2 = cdb[[c_genome, c_cluster]].copy()
    cdb2.columns = ["raw_bin", "cluster"]
    cdb2["raw_bin"] = cdb2["raw_bin"].map(strip_fasta_suffix)

    wdb2 = wdb[[w_genome, w_cluster]].copy()
    wdb2.columns = ["repMAG", "cluster"]
    wdb2["repMAG"] = wdb2["repMAG"].map(strip_fasta_suffix)

    merged = pd.merge(cdb2, wdb2, on="cluster", how="left")

    missing = merged["repMAG"].isna().sum()
    if missing > 0:
        print(f"[WARNING] {missing} raw bins have no matched representative MAG", file=sys.stderr)

    merged = merged[["raw_bin", "repMAG", "cluster"]].drop_duplicates()

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Output written to: {args.output}")
    print(f"[INFO] Total raw bins: {merged['raw_bin'].nunique()}")
    print(f"[INFO] Total representative MAGs: {merged['repMAG'].nunique()}")


if __name__ == "__main__":
    main()
