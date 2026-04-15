#!/usr/bin/env python3
import sys
import os
import glob
import pandas as pd

RED = '\033[91m'
GREEN = '\033[92m'
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
        
        if declared_male and not bio_male:
            return "Mismatch (Expected Male)"
        if declared_female and not bio_female:
            return "Mismatch (Expected Female)"
            
        return "Match"

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
        print(f"\n{'Sample':<12} | {'Status':<6} | {'Depth':<6} | {'Map %':<6} | {'N50':<6} | {'Mean Q':<6} | {'Sex Match':<20}")
        print("-" * 80)

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
            sex_str = sex_status
            if "Mismatch" in sex_status:
                reasons.append(f"Sex {sex_status}")
                sex_str = f"{RED}{sex_status}{RESET}"
                
            # Compile results
            status = "FAIL" if reasons else "PASS"
            status_str = f"{RED}FAIL{RESET}" if status == "FAIL" else f"{GREEN}PASS{RESET}"

            # Print the color-coded row to the terminal
            print(f"{s_id:<12} | {status_str:<15} | {depth_str:<15} | {mapped_str:<15} | {n50_str:<15} | {q_str:<15} | {sex_str}")

            final_data.append({
                'sample_id': s_id,
                'status': status,
                'depth': depth,
                'mapped_rate': mapped,
                'n50': n50,
                'mean_q': q,
                'sex_match': sex_status,
                'fail_reasons': " | ".join(reasons)

            })
            
        # 3. Save Report
        out_df = pd.DataFrame(final_data)
        out_df.to_csv(output_csv, index=False)
        print("-" * 80)
        print(f"\nClean CSV report saved to: {output_csv}\n")

# --- CLI Entry Point ---
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: evaluate_qc.py <multiqc_general_stats.txt> <mosdepth_dir> <output_report.csv>")
        sys.exit(1)
        
    parser = QCParser(multiqc_file=sys.argv[1], mosdepth_dir=sys.argv[2])
    parser.run(output_csv=sys.argv[3])
