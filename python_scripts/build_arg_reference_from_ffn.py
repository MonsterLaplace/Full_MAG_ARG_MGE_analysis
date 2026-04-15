#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Build ARG nucleotide reference fasta from Prodigal FFN files")
    parser.add_argument("--arg_table", required=True, help="contig_arg_mge.tsv")
    parser.add_argument("--ffn_dir", required=True, help="Directory containing sample .ffn files")
    parser.add_argument("--out_fasta", required=True, help="Output ARG reference fasta")
    parser.add_argument("--out_map", required=True, help="Output mapping table")
    args = parser.parse_args()

    arg_df = pd.read_csv(args.arg_table, sep="\t", dtype=str).fillna("")

    required_cols = ["Sample", "protein", "ARG_std", "ARG_group"]
    for c in required_cols:
        if c not in arg_df.columns:
            raise ValueError(f"Missing required column: {c}")

    arg_df = arg_df[arg_df["protein"] != ""].copy()

    # protein -> ARG info
    prot2info = {}
    for _, row in arg_df.iterrows():
        key = (row["Sample"], row["protein"])
        prot2info[key] = {
            "ARG_std": row["ARG_std"],
            "ARG_group": row["ARG_group"],
            "ARG": row.get("ARG", ""),
            "ARG_name": row.get("ARG_name", "")
        }

    records_out = []
    map_rows = []

    ffn_files = glob.glob(os.path.join(args.ffn_dir, "*.ffn"))
    if not ffn_files:
        raise FileNotFoundError(f"No .ffn files found in {args.ffn_dir}")

    found = 0

    for ffn in ffn_files:
        sample = os.path.basename(ffn).replace(".ffn", "")
        for rec in SeqIO.parse(ffn, "fasta"):
            gene_id = rec.id
            key = (sample, gene_id)
            if key in prot2info:
                info = prot2info[key]
                arg_std = info["ARG_std"] if info["ARG_std"] else "UnknownARG"
                arg_group = info["ARG_group"] if info["ARG_group"] else "UnknownGroup"

                new_id = f"{sample}|{gene_id}|{arg_std}"
                rec.id = new_id
                rec.description = ""

                records_out.append(rec)
                map_rows.append({
                    "ref_id": new_id,
                    "sample": sample,
                    "gene_id": gene_id,
                    "ARG_std": arg_std,
                    "ARG_group": arg_group,
                    "ARG": info["ARG"],
                    "ARG_name": info["ARG_name"]
                })
                found += 1

    if found == 0:
        raise RuntimeError("No ARG gene sequences found. Check Sample names and protein IDs.")

    SeqIO.write(records_out, args.out_fasta, "fasta")
    pd.DataFrame(map_rows).to_csv(args.out_map, sep="\t", index=False)

    print(f"[INFO] ARG reference fasta written to: {args.out_fasta}")
    print(f"[INFO] ARG mapping table written to: {args.out_map}")
    print(f"[INFO] Total ARG sequences: {found}")

if __name__ == "__main__":
    main()
