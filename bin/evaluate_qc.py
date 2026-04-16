#!/usr/bin/env python3
import sys
import os
import glob
import pandas as pd
import argparse

BOLD = '\033[1m'
CYAN = '\033[36m'
GREEN = '\033[32m'
RED = '\033[31m'
RESET = '\033[0m'

class QCParser:
    def __init__(self, multiqc_file, mosdepth_dir, min_depth=0.1, min_mapping=75.0, min_n50=1000, min_q=10.0):
        self.multiqc_file = multiqc_file
        self.mosdepth_dir = mosdepth_dir
        
        # Thresholds
        self.min_depth = min_depth
        self.min_mapping = min_mapping
        self.min_n50 = min_n50
        self.min_q = min_q

    def parse_mosdepth(self, sample_id):
        """Extracts total depth and chromosome ratios from a sample's mosdepth summary."""
        # Find the summary file in the provided directory
        pattern = os.path.join(self.mosdepth_dir, f"*{sample_id}*.mosdepth.summary.txt")
        files = glob.glob(pattern)
        
        if not files:
            return None, None, None
            
        df = pd.read_csv(files[0], sep='\t')
        
        try:
            total_depth = df[df['chrom'] == 'total']['mean'].values[0]
            
            # Fetch sex chromosomes
            chrx_depth = df[df['chrom'] == 'chrX']['mean'].values[0] if 'chrX' in df['chrom'].values else 0
            chry_depth = df[df['chrom'] == 'chrY']['mean'].values[0] if 'chrY' in df['chrom'].values else 0
            
            # Calculate ratios against total depth (proxy for autosomes)
            if total_depth == 0:
                return 0.0, 0.0, 0.0
                
            chrx_ratio = chrx_depth / total_depth
            chry_ratio = chry_depth / total_depth
            
            return total_depth, chrx_ratio, chry_ratio
            
        except Exception as e:
            print(f"Error parsing mosdepth for {sample_id}: {e}")
            return None, None, None

    def evaluate_sex(self, sample_id, chrx_ratio, chry_ratio):
        """Compares declared sex (_1 / _2) against biological depth ratios."""
        if chrx_ratio is None or chry_ratio is None:
            return "Unknown"
            
        declared_male = sample_id.endswith('_1')
        declared_female = sample_id.endswith('_2')
        
        if not declared_male and not declared_female:
            return "No Declared Sex"

        # Biological definitions based on ratios
        bio_male = chrx_ratio < 0.8 and chry_ratio > 0.1
        bio_female = chrx_ratio >= 0.8 and chry_ratio <= 0.1
        
        inferred_sex = "M" if bio_male else "F" if bio_female else "U"

        if declared_male:
            if bio_male:
                return "Match", inferred_sex
            elif bio_female:
                return "Mismatch (Expected Male)", inferred_sex
            else:
                return "Mismatch (Inconclusive)", inferred_sex

        elif declared_female:
            if bio_female:
                return "Match", inferred_sex
            elif bio_male:
                return "Mismatch (Expected Female)", inferred_sex
            else:
                return "Mismatch (Inconclusive)", inferred_sex 
        return "Unknown Declaration", inferred_sex

    def run(self, output_csv):
        # 1. Load MultiQC Data
        mqc_df = pd.read_csv(self.multiqc_file, sep='\t')
        
        # Map MultiQC columns to clean internal names
        col_map = {
            'Sample': 'sample_id',
            'samtools_flagstat-mapped_passed_pct': 'mapped_rate',
            'nanostat-Read_length_N50_aligned': 'n50',
            'nanostat-Mean_read_quality_aligned': 'mean_q'
        }
        
        cols_to_keep = [c for c in col_map.keys() if c in mqc_df.columns]
        df = mqc_df[cols_to_keep].rename(columns=col_map)
        
        # Strip nextflow/file extensions from sample ID
        df['sample_id'] = df['sample_id'].str.replace('.fixed', '', regex=False)
        
        final_data = []
        
        # Print header for terminal display
        header_format = f"\n{'Sample':<12} | {'Status':<6} | {'Depth':<6} | {'Map %':<6} | {'N50':<6} | {'Mean Q':<6} | {'Real Sex':<6} | {'Sex Match':<20}"
        print(f"\n{header_format}")
        table_width = len(header_format) + 2
        print("-" * table_width)

        # 2. Iterate through samples and apply logic
        for _, row in df.iterrows():
            s_id = row['sample_id']
            reasons = []
            
            # Fetch Mosdepth stats
            depth, x_ratio, y_ratio = self.parse_mosdepth(s_id)
            
            depth_str = f"{depth:.2f}" if depth is not None else "N/A"
            if depth is None or depth < self.min_depth:
                reasons.append(f"Low Depth ({depth_str}x)")
                depth_str = f"{RED}{depth_str}{RESET}"
                
            # Evaluate Mapping Rate
            mapped = row.get('mapped_rate')
            mapped_str = f"{mapped:.1f}" if pd.notna(mapped) else "N/A"
            if pd.notna(mapped) and mapped < self.min_mapping:
                reasons.append(f"Low Mapping ({mapped_str}%)")
                mapped_str = f"{RED}{mapped_str}{RESET}"

            # Evaluate N50
            n50 = row.get('n50')
            n50_str = f"{n50:.0f}" if pd.notna(n50) else "N/A"
            if pd.notna(n50) and n50 < self.min_n50:
                reasons.append(f"Low N50 ({n50})")
                n50_str = f"{RED}{n50_str}{RESET}"
                
            # Evaluate Mean Q-Score
            q = row.get('mean_q')
            q_str = f"{q:.1f}" if pd.notna(q) else "N/A"
            if pd.notna(q) and q < self.min_q:
                reasons.append(f"Low Q-Score ({q:.1f})")
                q_str = f"{RED}{q_str}{RESET}"
                
            # Evaluate Sex Check
            sex_status = self.evaluate_sex(s_id, x_ratio, y_ratio)
            sex_str = sex_status[0]
            sex_real = sex_status[1]
            if "Mismatch" in sex_str or "Unknown" in sex_str:
                reasons.append(f"Sex {sex_str}")
                sex_str = f"{RED}{sex_str}{RESET}"
                sex_real = f"{RED}{sex_real}{RESET}"

                
            # Compile results
            status = "FAIL" if reasons else "PASS"
            status_str = f"{RED}FAIL{RESET}" if status == "FAIL" else f"{GREEN}PASS{RESET}"

            # Print the color-coded row to the terminal
            print(f"{s_id:<12} | {status_str:<15} | {depth_str:<6} | {mapped_str:<6} | {n50_str:<6} | {q_str:<6} | {sex_real:<8} | {sex_str:<20}")

            final_data.append({
                'sample_id': s_id,
                'status': status,
                'depth': depth,
                'mapped_rate': mapped,
                'n50': n50,
                'mean_q': q,
                'inferred_sex': sex_real,
                'sex_match': sex_str,
                'fail_reasons': " | ".join(reasons)

            })
            
        # 3. Save Report
        out_df = pd.DataFrame(final_data)
        out_df.to_csv(output_csv, index=False)
        print("-" * table_width)
        print(f"\nClean CSV report saved to: {output_csv}\n")

# --- CLI Entry Point ---
if __name__ == "__main__":
    # 2. Setup the Parser
    description = f"{BOLD}{GREEN}QC Evaluator:{RESET} Parses MultiQC and Mosdepth data."
    usage = f"{BOLD}python3 %(prog)s{RESET} {CYAN}<multiqc_data.tsv>{RESET} {CYAN}<mosdepth_dir>{RESET} {CYAN}<output.csv>{RESET}"

    parser = argparse.ArgumentParser(
        description=description,
        usage=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # --- THE MISSING PART: CHECK FOR ARGUMENTS ---
    # sys.argv[0] is the script name, so if len is 1, no args were passed
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    # ---------------------------------------------

    # If your script uses positional arguments via sys.argv[1], sys.argv[2], etc.
    # instead of parser.add_argument, you need at least 4 args (script + 3 files)
    if len(sys.argv) < 4:
        print(f"{BOLD}Error:{RESET} Missing required positional arguments.")
        parser.print_help(sys.stderr)
        sys.exit(1)

    evaluator = QCParser(sys.argv[1], sys.argv[2])
    evaluator.run(sys.argv[3])
