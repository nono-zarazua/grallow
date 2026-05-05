#!/usr/bin/env python3

import os
import pandas as pd

print("Loading GSA-1 Manifest...")

script_dir = os.path.dirname(os.path.abspath(__file__))

csv_path = os.path.join(script_dir, "../../docs/GSA-24v1-0_C2.csv")


# skiprows=7 skips the Illumina metadata so pandas reads row 8 as the headers
manifest = pd.read_csv(csv_path, skiprows=7, low_memory=False)

print("Extracting coordinates and writing to gsa_targets.tsv...")
with open("../../docs/gsa_targets.tsv", "w") as out_file:
    for index, row in manifest.iterrows():
        chrom = str(row['Chr'])

        # Prevent "chrNaN" or empty rows at the end of the file
        if chrom == 'nan' or pd.isna(row['MapInfo']):
            continue

        # Add 'chr' prefix to match your QUILT 2 VCF format
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"

        pos = int(row['MapInfo'])

        # bcftools expects: CHROM [tab] POS
        out_file.write(f"{chrom}\t{pos}\n")

print("Done! Target list generated: gsa_targets.tsv")
