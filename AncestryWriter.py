#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import gzip
import sys
import os
import csv
import datetime
from zoneinfo import ZoneInfo

class AncestryWriter:
    # AncestryWriter format parser for vcf data
    def __init__(self, vcf_path, txt_path, compiled_tsv_path, gsa_manifest_path=None, rsid_map_path=None):
        self.vcf_path = vcf_path
        self.txt_path = txt_path
        self.compiled_tsv_path = compiled_tsv_path
        self.gsa_manifest_path = gsa_manifest_path
        self.rsid_map_path = rsid_map_path

        self.name_translator = {}
        self.target_snps = {}
        self.raw_snp_data = []
        self.chrom_map = {
            "chrX": "23",
            "chrY": "24",
            "chrXY": "25",
            "chrM": "26",
            "MT": "26"
        }
        for i in range(1,23):
            self.chrom_map[f"chr{i}"] = str(i)
        
        # Check if data compilation in tsv already exists
        if os.path.exists(self.compiled_tsv_path):
            self._load_compiled_tsv()
        else:
            if not self.gsa_manifest_path or not self.rsid_map_path:
                print(f"Error: '{self.compiled_tsv_path}' not found, and raw Illumina files not provided to generate it.")
                sys.exit(1)
            self._load_rsid_translator()
            self._load_gsa_manifest()
            self._save_compiled_tsv()
    
    def _load_compiled_tsv(self):
        """Instantly loads the pre-compiled target SNPs from our permanent TSV."""
        print(f"Loading pre-compiled reference from {self.compiled_tsv_path}...")
        comp_object, comp_filename = self._open_file(self.compiled_tsv_path)
        with comp_object as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 4:
                    rsid, illumina_id, chrom, pos = parts
                    self.target_snps[f"{chrom}_{pos}"] = rsid
        print(f"Successfully loaded {len(self.target_snps)} target SNPs.")

    def _save_compiled_tsv(self):
        """Writes the unified translation data to a permanent TSV file."""
        print(f"Saving compiled reference to {self.compiled_tsv_path}...")
        with open(self.compiled_tsv_path, 'w') as f:
            f.write("rsid\tillumina_id\tchr\tposition\n")
            for line in self.raw_snp_data:
                f.write(line)
        print("Permanent TSV reference created successfully!")

    def _load_rsid_translator(self):
        """Loads the Name to RsID cross-reference file."""
        print(f"Loading RsID cross-reference from {self.rsid_map_path}...")
        ref_object, ref_filename = self._open_file(self.rsid_map_path)
        with ref_object as f:
            for line in f:
                if line.startswith("Name") or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    illumina_name = parts[0]
                    real_rsid = parts[1]
                    if real_rsid != ".":
                        self.name_translator[illumina_name] = real_rsid
        print(f"Loaded {len(self.name_translator)} RsID translations.")
    
    def _load_gsa_manifest(self):
        """Parses the CSV manifest and applies the translations."""
        print(f"Loading GSA array template from {self.gsa_manifest_path}...")
        in_assay_section = False

        manifest_object, manifest_filename = self._open_file(self.gsa_manifest_path)
        with manifest_object as f:
            reader = csv.reader(f)
            headers = []
            for row in reader:
                if not row: continue
                if row[0] == '[Assay]':
                    in_assay_section = True
                    continue
                
                if in_assay_section:
                    if not headers:
                        headers = row
                        genome_build_idx = headers.index("GenomeBuild")
                        name_idx = headers.index('Name')
                        chr_idx = headers.index('Chr')
                        pos_idx = headers.index('MapInfo')
                        continue

                    if len(row) > pos_idx:
                        name = row[name_idx]
                        chrom = row[chr_idx]
                        pos = row[pos_idx]
                        gbuild = row[genome_build_idx]

                        # Only enforce build 38 if the column actually has data
                        if gbuild and gbuild != "38":
                            print(f"Error: this parser works exclusively with genome build 38. Found: {gbuild}")
                            sys.exit(1)
                        
                        if chrom == 'X': chrom = '23'
                        elif chrom == 'Y': chrom = '24'
                        elif chrom == 'XY': chrom = '25'
                        elif chrom == 'MT': chrom = '26'

                        # Apply translation if it exists
                        final_id = self.name_translator.get(name, name)
                        lookup_key = f"{chrom}_{pos}"
                        self.target_snps[lookup_key] = final_id

                        self.raw_snp_data.append(f"{final_id}\t{name}\t{chrom}\t{pos}\n")

        print(f"Successfully loaded and translated {len(self.target_snps)} target array SNPs.")

    def _open_file(self, filename):
        """Opens a file and returns the file object along with its uncompressed filename."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File '{filename}' not found.")

        # If it's a compressed .gz file
        if filename.endswith(".gz"):
            uncompressed_filename = os.path.splitext(filename)[0]
            file_obj = gzip.open(filename, 'rt')
            return file_obj, uncompressed_filename
        else:
            # If it's a normal file
            file_obj = open(filename, 'r')
            return file_obj, filename

    def process_genotype(self, ref, alt_string, gt_string, gp_string=None, min_confidence=0.90):
        """Converts VCF GT (e.g., 0|1) to actual DNA letters."""

        if gp_string and min_confidence > 0:
            try:
                # Convert "0.05,0.90,0.05" -> [0.05, 0.90, 0.05]
                gp_probs = [float(p) for p in gp_string.split(',')]

                if max(gp_probs) < min_confidence:
                    return "0", "0"
            except ValueError:
                return "0", "0"

        alleles = [ref] + alt_string.split(',')
        gt_clean = gt_string.replace('|', '/').split('/')

        if gt_clean[0] == ".":
            return "0", "0"
        
        try:
            a1 = alleles[int(gt_clean[0])]
            a2 = alleles[int(gt_clean[1])] if len(gt_clean) > 1 else a1
            return a1, a2
        except (IndexError, ValueError):
            return "0", "0"

    def convert(self):
        print(f"Converting {self.vcf_path} to Ancestry format...\n")

        vcf_obj, vcf_path = self._open_file(self.vcf_path)
        timestamp = datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %Z')
        variants_written = 0

        with vcf_obj as vcf, open(self.txt_path, 'w') as out:
            # Header aligned strictly to the left edge to prevent tab/space corruption
            out.write(f"""#AncestryDNA raw data download format
#This file was generated by Zengen at: {timestamp}
#Data was collected using low-coverage WGS imputation using QUILT2.
#Data is formatted using the AncestryWriter class.
#
#THIS INFORMATION IS FOR YOUR PERSONAL USE AND IS INTENDED FOR GENEALOGICAL RESEARCH 
#ONLY. IT IS NOT INTENDED FOR MEDICAL, DIAGNOSTIC, OR HEALTH PURPOSES. YOU ASSUME 
#ALL RISK OF STORING, SECURING, AND PROTECTING YOUR DOWNLOADED DATA.
#
#Genetic data is provided below as five TAB delimited columns. Each line 
#corresponds to a SNP. Column one provides the SNP identifier (rsID where 
#possible). Columns two and three contain the chromosome and basepair position 
#of the SNP using human reference build 38 coordinates. Columns four and five 
#contain the two alleles observed at this SNP (genotype). The genotype is reported 
#on the forward (+) strand with respect to the GRCh38 / hg38 human reference.
rsid\tchromosome\tposition\tallele1\tallele2
""")

            for line in vcf:
                if line.startswith('#'):
                    continue

                parts = line.strip().split("\t")
                chrom = parts[0]
                pos = parts[1]
                rsid = parts[2]
                ref = parts[3]
                alt = parts[4]

                if len(ref) > 1 or any(len(a) > 1 for a in alt.split(',')):
                    continue

                pos_int = int(pos)
                is_par = False
                if chrom == "chrX":
                    if (10000 <= pos_int <= 2781479) or (155701383 <= pos_int <= 156030895):
                        is_par = True
                
                if is_par:
                    mapped_chrom = "25"
                else:
                    mapped_chrom = self.chrom_map.get(chrom)
                
                if not mapped_chrom:
                    continue

                lookup_key = f"{mapped_chrom}_{pos}"

                if lookup_key in self.target_snps:
                    final_rsid = self.target_snps[lookup_key]
                else:
                    continue

                format_keys = parts[8].split(":")
                sample_data = parts[9].split(":")
                if "GT" not in format_keys: continue
                gt_string = sample_data[format_keys.index("GT")]

                gp_string = None
                if "GP" in format_keys:
                    gp_string = sample_data[format_keys.index("GP")]

                allele1, allele2 = self.process_genotype(ref, alt, gt_string, gp_string, min_confidence=0.90)
                out.write(f"{final_rsid}\t{mapped_chrom}\t{pos}\t{allele1}\t{allele2}\n")
                variants_written += 1
        
        print(f"Success! Filtered and wrote {variants_written} array-matched SNPs to {self.txt_path}")


if __name__ == "__main__":
    # We now require at least 3 arguments (VCF, TXT, Compiled_TSV) 
    # and up to 5 if the compiled TSV doesn't exist yet!
    if len(sys.argv) < 4:
        print("Usage: python vcf_to_ancestry.py <input.vcf.gz> <output.txt> <compiled_map.tsv> [gsa_manifest.csv] [rsid_map.txt]")
        sys.exit(1)
        
    vcf_in = sys.argv[1]
    txt_out = sys.argv[2]
    compiled_tsv = sys.argv[3]
    
    # Safely handle the optional raw Illumina files
    gsa_manifest = sys.argv[4] if len(sys.argv) >= 5 else None
    rsid_map = sys.argv[5] if len(sys.argv) == 6 else None
        
    converter = AncestryWriter(vcf_in, txt_out, compiled_tsv, gsa_manifest, rsid_map)
    converter.convert()