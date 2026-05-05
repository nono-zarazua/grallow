#!/usr/bin/env python3
import pysam
import pandas as pd
import argparse
import sys

# --- 1. Set up Command Line Arguments ---
parser = argparse.ArgumentParser(description="Translate a VCF into an Illumina GSA Final Report.")
parser.add_argument("-i", "--input", required=True, help="Path to the input VCF file (e.g., quilt_gsa_only.vcf.gz)")
parser.add_argument("-m", "--manifest", default="../grallow/docs/GSA-24v1-0_C2.csv", help="Path to the Illumina CSV Manifest")
parser.add_argument("-o", "--output", default="mock_GSA_final_report.txt", help="Path to save the output text file")
args = parser.parse_args()

print(f"1. Loading Manifest from: {args.manifest}")
try:
    manifest = pd.read_csv(args.manifest, skiprows=7, low_memory=False)
except FileNotFoundError:
    print(f"Error: Could not find manifest at {args.manifest}")
    sys.exit(1)

# --- 2. Build the Lookup Dictionary ---
probe_dict = {}
for index, row in manifest.iterrows():
    chrom = str(row['Chr'])
    if chrom == 'nan' or pd.isna(row['MapInfo']):
        continue
    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    pos = str(int(row['MapInfo']))
    
    probe_dict[f"{chrom}:{pos}"] = {
        'Name': row['Name'],
        'Strand': row['RefStrand']
    }

def flip_strand(allele):
    """Returns the reverse complement of a DNA base."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(allele, allele)

# --- 3. Translate the VCF ---
print(f"2. Translating VCF: {args.input}")
try:
    vcf_in = pysam.VariantFile(args.input)
except OSError:
    print(f"Error: Could not open VCF at {args.input}. Is it bgzipped?")
    sys.exit(1)

samples = list(vcf_in.header.samples)
hits = 0

print(f"3. Writing output to: {args.output}")
with open(args.output, "w") as out_file:
    # Write the standard Illumina header
    out_file.write("SNP_Name\tSample_ID\tAllele1_Forward\tAllele2_Forward\n")
    
    for record in vcf_in:
        key = f"{record.chrom}:{record.pos}"
        
        if key in probe_dict:
            hits += 1
            probe_info = probe_dict[key]
            
            for sample_name in samples:
                gt = record.samples[sample_name]['GT']
                
                if gt == (None, None) or gt is None:
                    a1, a2 = "-", "-"
                else:
                    raw_a1 = record.alleles[gt[0]]
                    raw_a2 = record.alleles[gt[1]]
                    
                    if probe_info['Strand'] == '-':
                        a1 = flip_strand(raw_a1)
                        a2 = flip_strand(raw_a2)
                    else:
                        a1 = raw_a1
                        a2 = raw_a2
                
                out_file.write(f"{probe_info['Name']}\t{sample_name}\t{a1}\t{a2}\n")

print(f"\nDone! Successfully translated {hits} loci across {len(samples)} samples.")
