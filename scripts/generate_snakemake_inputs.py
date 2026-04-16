#!/usr/bin/env python3
import csv
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Generate Snakemake TSV inputs. Samples present in the CSV will be processed.")
    parser.add_argument('--qc-csv', required=True, help="Path to the input qc_summary.csv")
    parser.add_argument('--run-name', required=True, help="Run name for building BAM paths")
    parser.add_argument('--out-samples', required=True, help="Path to output samples.tsv")
    parser.add_argument('--out-sex', required=True, help="Path to output samples-sex.tsv")
    args = parser.parse_args()

    if not os.path.exists(args.qc_csv):
        print(f"Error: QC summary file '{args.qc_csv}' not found.")
        sys.exit(1)

    samples_to_process = []
    
    # Read the CSV. Every row remaining in this file will be sent to Snakemake.
    with open(args.qc_csv, mode='r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row['sample_id'].strip()
            if not sample_id:
                continue

            # Default logic for sex remains based on sample suffix
            if sample_id.endswith('_1'):
                sex = 'M'
            elif sample_id.endswith('_2'):
                sex = 'F'
            else:
                sex = 'U'

            samples_to_process.append({
                'sample_id': sample_id,
                'depth': row.get('depth', 'NA'),
                'sex': sex
            })

    if not samples_to_process:
        print("Warning: No samples found in the CSV. Snakemake files will be empty.")

    # Write samples.tsv (sampleid, bam, depth)
    with open(args.out_samples, mode='w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sampleid', 'bam', 'depth'])
        for s in samples_to_process:
            # Construct path matching your Snakemake expectations [cite: 79]
            bam_path = f"data/bams/{args.run_name}/{s['sample_id']}.sorted.aligned.bam"
            writer.writerow([s['sample_id'], bam_path, s['depth']])

    # Write samples-sex.tsv (sample, sex)
    with open(args.out_sex, mode='w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample', 'sex'])
        for s in samples_to_process:
            writer.writerow([s['sample_id'], s['sex']])

    print(f"Done. Prepared {len(samples_to_process)} samples for Snakemake.")

if __name__ == '__main__':
    main()
