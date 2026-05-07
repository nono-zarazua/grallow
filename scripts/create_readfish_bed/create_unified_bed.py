#!/usr/bin/env python3
"""
This script has the objective of reading downloading epigenetic clocks data from pyaging package to use them to filter EPIC-v2
manifest for further BED formatting and merging with GSA-1 BED,
"""

#!/usr/bin/env python3
import pandas as pd
import pyaging as pya
import os
import sys

# 1. Configuration & Paths
epic_manifest_path = 'EPIC-8v2-0_A2.csv'
gsa_1_bed_path = 'GSA-24v3-0_A1.bed'
OUTPUT_PATH = 'epiclock_target_hg38.bed'

valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M'] 
valid_chroms_prefixed = ['chr' + c for c in valid_chroms]

# 2. Setup pyaging environment
os.makedirs('pyaging_data', exist_ok=True)
logger = pya.logger.Logger('extract_bed_logger')
device = 'cpu'
data_dir = 'pyaging_data'

list_clocks = ['grimage2', 'phenoage', 'horvath2013', 'dunedinpace']
all_features = set()

# 3. Extract pyaging data
for clock in list_clocks:
    print(f"📦 Loading clock: {clock}")
    clock_obj = pya.pred.load_clock(clock, device, data_dir, logger)

    if hasattr(clock_obj, 'features'):
        all_features.update(clock_obj.features)
        print(f"   Added sites: {len(clock_obj.features)}")

print(f"\n🎯 Total unique CpGs for filtering: {len(all_features)}")

# 4. Load and filter EPIC manifest
print("📖 Reading EPIC v2 manifest...")
cols_to_use = ['Name', 'CHR', 'MAPINFO']
df_epic = pd.read_csv(epic_manifest_path, skiprows=7, usecols=cols_to_use, low_memory=False)

# Filter by pyaging IDs and remove non-mappable/garbage data
df_clocks = df_epic[df_epic['Name'].isin(all_features)].copy()
df_clocks = df_clocks[
    (df_clocks['CHR'].astype(str).isin(valid_chroms)) & 
    (df_clocks['MAPINFO'] > 0)
].copy()

# Create BED columns
df_clocks['Start'] = df_clocks['MAPINFO'].astype(int) - 1
df_clocks['End'] = df_clocks['MAPINFO'].astype(int)
df_clocks['chrom'] = 'chr' + df_clocks['CHR'].astype(str)
df_clocks = df_clocks[['chrom', 'Start', 'End', 'Name']]

print(f"✅ Filtered epigenomic sites: {len(df_clocks)}")

# 5. Load GSA SNPs
print(f"📖 Loading GSA manifest: {gsa_1_bed_path}")
df_gsa = pd.read_csv(gsa_1_bed_path, sep='\t', skiprows=1, header=None, 
                     names=['chrom', 'Start', 'End', 'Name'])
print(f"✅ Loaded GSA sites: {len(df_gsa)}")

# 6. Unify and Categorical Sort
print("🔀 Unifying and sorting target BED...")
df_unified = pd.concat([df_gsa, df_clocks], axis=0)

df_unified['chrom'] = pd.Categorical(
    df_unified['chrom'], 
    categories=valid_chroms_prefixed, 
    ordered=True
)

# Final sort: Chromosome (natural order) then numerical Start position
df_final = df_unified.sort_values(by=['chrom', 'Start'])

# 7. Output to disk
df_final.to_csv(OUTPUT_PATH, sep='\t', header=False, index=False)

print(f"\n🚀 DONE! Target BED ready for ReadFish.")
print(f"📊 Final file: {OUTPUT_PATH} ({len(df_final)} positions)")