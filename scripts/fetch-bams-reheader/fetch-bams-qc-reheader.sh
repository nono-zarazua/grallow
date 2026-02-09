#!/bin/bash

source /home/ec2-user/workdir/scripts/utils.sh

set -u

# Point to Conda envs' bins
CONDA_ROOT="/home/ec2-user/miniconda3/envs"
SAMTOOLS="$CONDA_ROOT/bstools/bin/samtools"
BCFTOOLS="$CONDA_ROOT/bstools/bin/bcftools"
MOSDEPTH="$CONDA_ROOT/mosdepth/bin/mosdepth"
SOMALIER="$CONDA_ROOT/somalier/bin/somalier"
NANOPLOT="$CONDA_ROOT/multi-nano/bin/NanoPlot"
MULTIQC="$CONDA_ROOT/multi-nano/bin/multiqc"

export PATH="$CONDA_ROOT/bstools/bin:$PATH"

# Configuration
S3_INPUT_ROOT="s3://omics-data-002313225286/clients/gtria/prod/inputs"
S3_OUTPUT_ROOT="s3://omics-data-002313225286/clients/gtria/prod/outputs"
EC_OUTPUT_ROOT="/home/ec2-user/workdir/project-quilt-workdir/data/bams"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOMALIER_SITES="${SCRIPT_DIR}/sites.hg38.vcf.gz"  # <--- FIXED: Absolute path
FASTA_REF="/home/ec2-user/workdir/project-quilt-workdir/data/ref/ref.fa"

# --- 1. STRICT PARSER ---
while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--csv)
      # Check if file exists before processing
      if [[ ! -f "$2" ]]; then
          echo "Error: Input CSV '$2' not found."
          exit 1
      fi
      # Convert to Absolute Path immediately
      # 1. cd into the file's directory to get absolute dir path
      # 2. Add the filename back
      INPUT_CSV="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: $0 --csv <input_csv>"
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1"
      exit 1
      ;;
  esac
done


# --- 2. VALIDATION (Fail if no dir provided) ---
if [[ -z "$INPUT_CSV" ]]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 --csv <file.csv>"
    exit 1
fi

# CHECK IF CSV EXISTS
if [[ ! -f "$INPUT_CSV" ]]; then
    echo "Error: Could not find $INPUT_CSV."
    exit 1
fi

echo "Input CSV: $INPUT_CSV"

# CREATE OUTPUT DIRECTORY
BASENAME_FILE=$(basename $INPUT_CSV)
OUTPUT_DIR_NAME=$(substring "$BASENAME_FILE" 0 18 )
OUTPUT_DIR="${EC_OUTPUT_ROOT}/${OUTPUT_DIR_NAME}"

mkdir -p $OUTPUT_DIR
cp $INPUT_CSV $OUTPUT_DIR
cd $OUTPUT_DIR

# DEFINE THE HEADERS WE ARE LOOKING FOR
TARGET_HEADERS=("alias" "barcode" "sample_id")
RESULT_STRING=$(csvcols "${TARGET_HEADERS[*]}" "$INPUT_CSV" "csv")
read -r -a RESULT_INDICES <<< "$RESULT_STRING"

for ((i=0; i<${#TARGET_HEADERS[@]}; i++)); do
	header_name="${TARGET_HEADERS[i]}"
	col_index="${RESULT_INDICES[i]}"

	if [[ "$col_index" -eq 0 ]]; then
		echo "Warning: Column '$header_name' not foind in CSV."
	else
		echo "Found '$header_name' at column $col_index"
	fi
done

# Create qc folders
mkdir -p "${OUTPUT_DIR}/qc/mosdepth" "${OUTPUT_DIR}/qc/flagstat" "${OUTPUT_DIR}/qc/nanoplot" "${OUTPUT_DIR}/qc/somalier"

# -- 3. LOOP THROUGH THE CSV FILE --
tail -n +2 "$INPUT_CSV" | csvextract "$RESULT_STRING" | while IFS=, read -r alias run_id sample_id; do
	# Cleanup invisible character (important for CSVs from Excel/Windows)
        alias=$(echo "$alias" | tr -d '\r')
        run_id=$(echo "$run_id" | tr -d '\r')
        sample_id=$(echo "$sample_id" | tr -d '\r')

	if [[ -z $alias ]]; then continue; fi
	
	# Construct path
	S3_BAM="${S3_INPUT_ROOT}/${sample_id}/${alias}.sorted.aligned.bam"
	NAME_BAM="${alias}.sorted.aligned.bam"
	LOCAL_BAM="${OUTPUT_DIR}/${NAME_BAM}"
	echo "Processing $alias (Run: $run_id)..."
	
	# Copy BAM and BAI
	aws s3 cp "$S3_BAM" "$LOCAL_BAM" --quiet
	
	if [[ $? -ne 0 ]]; then
		echo " X Error downloading $alias. Skipping."
		continue
	fi
	

	# Reheader
	echo "  -> Fixing header (adding 'chr')..."
	bamrehead "$LOCAL_BAM" index
	
	# QC
	echo "Running QC for $alias."

	# Mosdepth
	$MOSDEPTH -n -x "${OUTPUT_DIR}/qc/mosdepth/sample_coverage_${alias}" "$LOCAL_BAM" &
	
	# Flagstat
	$SAMTOOLS flagstat "$LOCAL_BAM" > "${OUTPUT_DIR}/qc/flagstat/${alias}.flagstat" &

	# Nanoplot	
	$NANOPLOT --bam "$LOCAL_BAM" \
		  --outdir "${OUTPUT_DIR}/qc/nanoplot/${alias}" \
		  --maxlength 4000 \
		  -p "${alias}_" &

	# Somalier Extract
	export SOMALIER_SAMPLE_NAME="$alias"
	$SOMALIER extract -d "${OUTPUT_DIR}/qc/somalier" \
	       --sites $SOMALIER_SITES \
	       -f $FASTA_REF \
	       "$LOCAL_BAM" &

	wait
done

# Somalier Relate
echo "Running Somalier Relate on cohort..."
if ls ${OUTPUT_DIR}/qc/somalier/*.somalier 1> /dev/null 2>&1; then
	$SOMALIER relate --infer \
		--output-prefix cohort_somalier \
		"${OUTPUT_DIR}/qc/somalier"/*.somalier
fi

# MultiQC
$MULTIQC .

echo "-------------------------------------------------------------"
echo "✓ Done! All tasks finished."
