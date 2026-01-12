#!/bin/bash

set -u

# initialize Conda for this script's shell
source ~/miniconda3/etc/profile.d/conda.sh

# Activate env
conda activate bcftools

INPUT_VCF=""

# --- 1. STRICT PARSER ---
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      # Check if file exists before processing
      if [[ ! -f "$2" ]]; then
          echo "Error: Input VCF '$2' not found."
          exit 1
      fi
      INPUT_VCF=$(realpath "$2")
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: $0 --input <input_vcf>"
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1"
      exit 1
      ;;
  esac
done


# --- VALIDATION ---
if [[ -z "$INPUT_VCF" ]]; then
    echo "Error: You must provide an input file using -i or --input."
    exit 1
fi

# Determine output directory
BASE_DIR=$(dirname "$INPUT_VCF")
OUT_DIR="${BASE_DIR}/split_files"

echo "Input:  $INPUT_VCF"
echo "Output: $OUT_DIR"

mkdir -p "$OUT_DIR"

bcftools +split "$INPUT_VCF" -O z -o "${OUT_DIR}"

echo "Done! Individual VCFs are in: $OUT_DIR"
echo "Indexing split VCFs..."

for vcf in "$OUT_DIR"/*.vcf.gz; do
	[ -e "$vcf" ] || continue

	echo "Indexing $(basename "$vcf")..."
	bcftools index -t "$vcf"
done

echo "Spliting and indexing completed!"
	
