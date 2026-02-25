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
SOMALIER_SITES="${SCRIPT_DIR}/sites.hg38.vcf.gz"  
FASTA_REF="/home/ec2-user/workdir/project-quilt-workdir/data/ref/GRCh38_chr.fa"

# Config files
CONFIG_DIR="/home/ec2-user/workdir/project-quilt-workdir/config"
SAMPLES_TSV="${CONFIG_DIR}/samples.tsv"
SAMPLES_SEX_TSV="${CONFIG_DIR}/samples-sex.tsv"
CONFIG_YAML="${CONFIG_DIR}/config.yaml"

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

# INITIALIZE TSV FILES WITH HEADERS
echo "Initializing configuration TSVs..."
echo -e "sampleid\tbam\tdepth" > "$SAMPLES_TSV"
echo -e "sample\tsex" > "$SAMPLES_SEX_TSV"

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
	
	# Copy BAM
	echo "Downloading $S3_BAM ..."
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
	$MOSDEPTH -n -x "${OUTPUT_DIR}/qc/mosdepth/${alias}" "$LOCAL_BAM" &
	
	# Flagstat
	$SAMTOOLS flagstat "$LOCAL_BAM" > "${OUTPUT_DIR}/qc/flagstat/${alias}.flagstat" &

	# Nanoplot	
	$NANOPLOT --bam "$LOCAL_BAM" \
		  --outdir "${OUTPUT_DIR}/qc/nanoplot/${alias}" \
		  --maxlength 4000 \
		  --no_static \
		  -p "${alias}_" &

	# Somalier Extract
	export SOMALIER_SAMPLE_NAME="$alias"
	$SOMALIER extract -d "${OUTPUT_DIR}/qc/somalier" \
	       --sites $SOMALIER_SITES \
	       -f $FASTA_REF \
	       "$LOCAL_BAM" &

	wait

	echo "  -> Generating Snakemake config entries for $alias..."

	MOSDEPTH_SUMMARY="${OUTPUT_DIR}/qc/mosdepth/${alias}.mosdepth.summary.txt"
	if [[ -f "$MOSDEPTH_SUMMARY" ]]; then
		MEAN_DEPTH=$(awk '$1=="total" {printf "%.2f", $4}' "$MOSDEPTH_SUMMARY")
	else
                MEAN_DEPTH="NA"
		echo "     Warning: mosdepth summary not found."
	fi

	# Determine sex
	if [[ "$alias" == *_1 ]]; then
		SEX="M"
	elif [[ "$alias" == *_2 ]]; then
		SEX="F"
	else
		echo "     Warning: Alias does not end in _1 or _2."
		SEX="U"
	fi

	# Append to TSV files
	REL_BAM_PATH="data/bams/${OUTPUT_DIR_NAME}/${NAME_BAM}"

	echo -e "${alias}\t${REL_BAM_PATH}\t${MEAN_DEPTH}" >> "$SAMPLES_TSV"
	echo -e "${alias}\t${SEX}" >> "$SAMPLES_SEX_TSV"
done

# Somalier Relate
echo "Running Somalier Relate on cohort..."
if ls ${OUTPUT_DIR}/qc/somalier/*.somalier 1> /dev/null 2>&1; then
	$SOMALIER relate --infer \
		--output-prefix cohort_somalier \
		"${OUTPUT_DIR}/qc/somalier"/*.somalier
fi

# UPDATE SNAKEMAKE CONFIG.YAML
echo "Updating config.yaml with run_name: ${OUTPUT_DIR_NAME}..."
if [[ -f "$CONFIG_YAML" ]]; then
    # Use sed to safely replace the line starting with run_name:
    sed -i 's/^run_name:.*/run_name: "'"${OUTPUT_DIR_NAME}"'"/' "$CONFIG_YAML"
else
    echo "Warning: $CONFIG_YAML not found. Cannot update run_name."
fi

# MultiQC
$MULTIQC .

echo "-------------------------------------------------------------"
echo "✓ Done! All tasks finished."
