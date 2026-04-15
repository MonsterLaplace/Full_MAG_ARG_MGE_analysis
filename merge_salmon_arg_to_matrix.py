#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Merge Salmon ARG quant results into sample x ARG matrix")
    parser.add_argument("--salmon_dir", required=True, help="Directory containing per-sample Salmon outputs")
    parser.add_argument("--arg_map", required=True, help="ARG mapping table from build_arg_reference_from_ffn.py")
    parser.add_argument("--output", required=True, help="Output sample_arg_abundance.tsv")
    parser.add_argument("--feature_level", default="ARG_std", choices=["ARG_std", "ARG_group", "ARG"], help="Feature level")
    parser.add_argument("--value_col", default="TPM", choices=["TPM", "NumReads"], help="Quant value")
    args = parser.parse_args()

    arg_map = pd.read_csv(args.arg_map, sep="\t", dtype=str).fillna("")
    required_cols = ["ref_id", args.feature_level]
    for c in required_cols:
        if c not in arg_map.columns:
            raise ValueError(f"Missing required column in arg_map: {c}")

    ref2feature = arg_map.set_index("ref_id")[args.feature_level].to_dict()

    sample_tables = []

    quant_files = sorted(glob.glob(os.path.join(args.salmon_dir, "*", "quant.sf")))
    if not quant_files:
        raise FileNotFoundError(f"No quant.sf files found under {args.salmon_dir}")

    for qf in quant_files:
        sample = os.path.basename(os.path.dirname(qf))
        df = pd.read_csv(qf, sep="\t")

        if "Name" not in df.columns or args.value_col not in df.columns:
            raise ValueError(f"quant.sf missing required columns in {qf}")

        df["feature"] = df["Name"].map(ref2feature)
        df = df[df["feature"].notna() & (df["feature"] != "")].copy()

        agg = df.groupby("feature")[args.value_col].sum()
        agg.name = sample
        sample_tables.append(agg)

    mat = pd.concat(sample_tables, axis=1).fillna(0).T
    mat.index.name = "sample_id"
    mat.reset_index().to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Output written to: {args.output}")
    print(f"[INFO] Matrix shape: {mat.shape}")

if __name__ == "__main__":
    main()
