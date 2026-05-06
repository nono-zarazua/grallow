#!/usr/bin/env python3
import pysam
import pandas as pd
import argparse
import sys
import os
import json

def flip_strand(allele):
    """Returns the reverse complement of a DNA base."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(allele, allele)

def main():
    # --- 1. Set up Command Line Arguments ---
    parser = argparse.ArgumentParser(description="Hybrid VCF to Illumina GSA Final Report with QC, Caching, and INDEL support.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input VCF file")
    parser.add_argument("-m", "--manifest", default="../docs/GSA-24v1-0_C2.csv", help="Path to the Illumina CSV Manifest")
    parser.add_argument("-o", "--output", default="mock_GSA_final_report.txt", help="Path to save the output text file")
    parser.add_argument("-c", "--cache", default="gsa_manifest_cache.json", help="Path to save/load the compiled dictionary")
    parser.add_argument("--min_gp", type=float, default=0.90, help="Minimum Genotype Probability (GP) to keep a call (default: 0.90)")
    parser.add_argument("--rsid_map", help="Optional path to a TSV mapping Illumina names to rsIDs")
    parser.add_argument("--use_rsid", action="store_true", help="Toggle to output rsIDs instead of Illumina names")
    args = parser.parse_args()

    # --- 2. Build or Load the Lookup Dictionary ---
    # The dictionary maps loci to the ID name from the rsid.txt and the RefStrand orientation from the manifes
    # dict(f"{chrom}:{pos}" => dict("Name" => final_name, "Strand" => str(ref_strand)))
    probe_dict = {}

    if os.path.exists(args.cache):
        print(f"1. Loading pre-compiled manifest from cache: {args.cache}")
        with open(args.cache, 'r') as f:
            probe_dict = json.load(f)
    else:
        print(f"1. Compiling Manifest from CSV: {args.manifest}")
        try:
            # First 7 lines are header
            manifest = pd.read_csv(args.manifest, skiprows=7, low_memory=False)
        except FileNotFoundError:
            print(f"Error: Could not find manifest at {args.manifest}")
            sys.exit(1)


        # Translate the gsa-1 loci to its rsid
        # If rsid file is not given, gsa-1 names will be used
        # Maps illumina generic name to rsid
        name_translator = {}
        if args.rsid_map and os.path.exists(args.rsid_map):
            print(f"Loading rsID map from {args.rsid_map}...")
            with open(args.rsid_map, 'r') as f:
                for line in f:
                    if line.startswith("Name") or line.startswith("#"): continue
                    parts = line.strip().split()
                    if len(parts) >= 2 and parts[1] != ".":
                        name_translator[parts[0]] = parts[1]

        print("Building dictionary...")
        for index, row in manifest.iterrows():
            chrom = str(row['Chr'])
            if chrom == 'nan' or pd.isna(row['MapInfo']): continue
            if not chrom.startswith('chr'): chrom = f"chr{chrom}"
            
            pos = str(int(row['MapInfo']))
            orig_name = row['Name']
            rsid_name = name_translator.get(orig_name, orig_name)

            # Store BOTH names permanently in the dictionary
            probe_dict[f"{chrom}:{pos}"] = {
                'Illumina_Name': orig_name,
                'rsID': rsid_name,
                'Strand': str(row['RefStrand'])
            }
        
        print(f"Saving compiled dictionary to {args.cache} for future runs...")
        with open(args.cache, 'w') as f:
            json.dump(probe_dict, f)

    # --- 3. Translate the VCF ---
    print(f"\n2. Translating VCF: {args.input} with GP threshold: {args.min_gp}")
    try:
        # Load vcf into memory using pysam package
        vcf_in = pysam.VariantFile(args.input)
    except OSError:
        print(f"Error: Could not open VCF at {args.input}. Is it bgzipped?")
        sys.exit(1)

    samples = list(vcf_in.header.samples)
    hits = 0
    passed_gp = 0
    failed_gp = 0

    print(f"3. Writing output to: {args.output}")
    with open(args.output, "w") as out_file:
        out_file.write("SNP_Name\tSample_ID\tAllele1_Forward\tAllele2_Forward\n")
        
        # Extract from the vcf the chromosome and position of the snp to match gsa format
        for record in vcf_in:
            key = f"{record.chrom}:{record.pos}"
            
            # Check if the SNP from the VCF is in the dictionary
            if key in probe_dict:
                hits += 1
                # Use the built key to access dictionary as a double safeguard
                probe_info = probe_dict[key]
                
                # Get the SNP for a specific loci for every sample
                for sample_name in samples:
                    sample_data = record.samples[sample_name]
                    # Get the genotype ((0,0),(0,1),(1,1))
                    gt = sample_data.get('GT', (None, None))
                    # Get the statistical likelihood for every possible biological state
                    # First vaue: prob that it is homozygous
                    # Second value: prob that it is heterozygous
                    # third value: prob that it is homozygous-alternate
                    gp = sample_data.get('GP', None)
                    confident = True
                    
                    # Likelihood filtering
                    if gp is not None and args.min_gp > 0:
                        if max(gp) < args.min_gp:
                            confident = False
                            failed_gp += 1
                        else:
                            passed_gp += 1
                    
                    # If null, missing, or low confidence -> Write "-"
                    if gt == (None, None) or None in gt or not confident:
                        a1, a2 = "-", "-"
                    else:
                        raw_a1 = record.alleles[gt[0]]
                        raw_a2 = record.alleles[gt[1]]
                        
                        # --- NEW INDEL LOGIC ---
                        is_indel = False
                        if record.alts:
                            for alt in record.alts:
                                if len(alt) != len(record.ref):
                                    is_indel = True
                                    break
                        
                        if is_indel:
                            max_len = max([len(record.ref)] + [len(a) for a in record.alts])
                            a1 = "I" if len(raw_a1) == max_len else "D"
                            a2 = "I" if len(raw_a2) == max_len else "D"
                        # --- END INDEL LOGIC ---
                        
                        # Standard SNPs
                        else:
                            if probe_info['Strand'] == '-':
                                a1 = flip_strand(raw_a1)
                                a2 = flip_strand(raw_a2)
                            else:
                                a1 = raw_a1
                                a2 = raw_a2
                    
                    # Decide which name to use based on your command line flag
                    final_output_name = probe_info['rsID'] if args.use_rsid else probe_info['Illumina_Name']
                    
                    out_file.write(f"{final_output_name}\t{sample_name}\t{a1}\t{a2}\n")

    print(f"\nDone! Translated {hits} loci across {len(samples)} samples.")
    if args.min_gp > 0:
        total_calls = passed_gp + failed_gp
        if total_calls > 0:
            fail_percent = (failed_gp / total_calls) * 100
            print(f"QC Filter: {passed_gp} passed, {failed_gp} failed ({fail_percent:.2f}% removed for < {args.min_gp} GP).")

if __name__ == "__main__":
    main()