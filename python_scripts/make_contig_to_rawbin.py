#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import glob
import re


def strip_fasta_suffix(x):
    return re.sub(r'\.(fa|fasta|fna|fas)$', '', x)


def parse_fasta_headers(fasta_file):
    headers = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                contig = line[1:].strip().split()[0]
                headers.append(contig)
    return headers


def main():
    parser = argparse.ArgumentParser(
        description="Build contig -> raw_bin table from raw bin fasta files"
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Root directory of raw bins")
    parser.add_argument("-o", "--output", required=True, help="Output contig_to_rawbin.tsv")
    parser.add_argument("--sample_from_parent_dir", action="store_true",
                        help="Infer sample name from parent directory name")
    parser.add_argument("--sample_regex", default=None,
                        help="Optional regex to extract sample from raw_bin filename")
    args = parser.parse_args()

    fasta_files = glob.glob(os.path.join(args.input_dir, "**", "*.fa"), recursive=True)
    fasta_files += glob.glob(os.path.join(args.input_dir, "**", "*.fasta"), recursive=True)
    fasta_files += glob.glob(os.path.join(args.input_dir, "**", "*.fna"), recursive=True)
    fasta_files = sorted(set(fasta_files))

    with open(args.output, "w") as out:
        out.write("Sample\tContig\traw_bin\tContigKey\n")

        for fa in fasta_files:
            raw_bin = strip_fasta_suffix(os.path.basename(fa))

            sample = "unknown"

            if args.sample_from_parent_dir:
                sample = os.path.basename(os.path.dirname(fa))

            if args.sample_regex:
                m = re.search(args.sample_regex, raw_bin)
                if m:
                    sample = m.group(1) if m.groups() else m.group(0)

            headers = parse_fasta_headers(fa)

            for contig in headers:
                contig_key = f"{sample}|{contig}"
                out.write(f"{sample}\t{contig}\t{raw_bin}\t{contig_key}\n")

    print(f"[INFO] Output written to: {args.output}")
    print(f"[INFO] Parsed {len(fasta_files)} raw bin fasta files")


if __name__ == "__main__":
    main()
