#!/usr/bin/env python3

import os
import sys
import re
import gzip

class QCParser:
    def __init__(self, qc_dir):
        self.qc_dir = qc_dir
        self.metrics = {}

    # ----- PRIVATE METHODS -----
    def _open_file(self, filename):
        """Opens a file and returns the file object."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File '{filename}' not found.")

        if filename.endswith(".gz"):
            return gzip.open(filename, 'rt')
        else:
            return open(filename, 'r')

    # ----- PUBLIC METHODS -----
    def parse_flagstat(self, sample_id): 
        """
        Parses samtools flagstat output to get Total Reads and Mapping %.
        """
        filename = os.path.join(self.qc_dir, "flagstat", f"{sample_id}.flagstat")

        stats = {
            'total_reads': 0,
            'mapped_reads': 0,
            'mapped_pct': 0.0
        }

        if not os.path.exists(filename):
            print(f"Warning: File not found {filename}")
            return stats

        try:
            file_obj = self._open_file(filename)
            
            with file_obj as f:
                content = f.read() 

                total_match = re.search(r'^(\d+) \+ \d+ in total', content)
                if total_match:
                    stats['total_reads'] = int(total_match.group(1))

                mapped_match = re.search(r'(\d+) \+ \d+ mapped \(([\d\.]+)%', content)
                if mapped_match:
                    stats['mapped_reads'] = int(mapped_match.group(1))
                    stats['mapped_pct'] = float(mapped_match.group(2))

        except Exception as e:
            print(f"Error parsing flagstat for {sample_id}: {e}")

        for key,value in stats.items():
            self.metrics[key] = value

# --- TEST BLOCK ---
if __name__ == "__main__":
    qc_path = "/home/ec2-user/workspace/trials/qc"
    test_sample = "B025-02_1"

    parser = QCParser(qc_path)
    result = parser.parse_flagstat(test_sample)
    
    for key,value in parser.metrics.items():
        print(f"{key}: {value}")
