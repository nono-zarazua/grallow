#!/bin/bash

set -u

# initialize Conda for this script's shell
source ~/miniconda3/etc/profile.d/conda.sh

# Activate env
conda activate bcftools

# Configuration
S3_ROOT="s3://omics-data-002313225286/clients/gtria/prod/inputs"
OUTPUT_DIR="" # Empty by default to force validation

# --- 1. STRICT PARSER ---
while [[ $# -gt 0 ]]; do
  case $1 in
    -d|--dir|--directory)
      mkdir -p "$2"
      OUTPUT_DIR=$(cd "$2" && pwd)
      shift # past argument
      shift # past value
      ;;
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
      echo "Usage: $0 --dir <output_directory> --csv <input_csv>"
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1"
      exit 1
      ;;
  esac
done


# --- 2. VALIDATION (Fail if no dir provided) ---
if [[ -z "$OUTPUT_DIR" || -z "$INPUT_CSV" ]]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 --dir <directory_name> --csv <file.csv>"
    exit 1
fi


# Check if CSV exists
if [[ ! -f "$INPUT_CSV" ]]; then
    echo "Error: Could not find $INPUT_CSV. Run the submission script first!"
    exit 1
fi

echo "Input CSV:  $INPUT_CSV"
echo "Output Dir: $OUTPUT_DIR"
echo "Downloading BAM file to: $OUTPUT_DIR"

# Define the headers we are looking for
ALIAS_COL_NAME="alias"
RUNID_COL_NAME="run_id"
SAMPLE_COL_NAME="sample_id"

# Find index for alias
ALIAS_IDX=$(awk -F, -v target="$ALIAS_COL_NAME" 'NR==1{for(i=1;i<=NF;i++){if($i==target){print i; exit}}}' "$INPUT_CSV")
# Find index for barcode
RUN_IDX=$(awk -F, -v target="$RUNID_COL_NAME" 'NR==1{for(i=1;i<=NF;i++){if($i==target){print i; exit}}}' "$INPUT_CSV")
# Find index for barcode
SAMPLE_IDX=$(awk -F, -v target="$SAMPLE_COL_NAME" 'NR==1{for(i=1;i<=NF;i++){if($i==target){print i; exit}}}' "$INPUT_CSV")

# Safety check
if [[ -z "$ALIAS_IDX" || -z "$RUN_IDX" || -z "$SAMPLE_IDX" ]]; then
	echo "Error: Could not find columns '$ALIAS_IDX', '$RUN_IDX', or '$SAMPLE_IDX' in $INPUT_CSV"
	exit 1
fi

echo "Found '$ALIAS_COL_NAME' at column $ALIAS_IDX"
echo "Found '$RUNID_COL_NAME' at column $RUN_IDX"
echo "Found '$SAMPLE_COL_NAME' at column $SAMPLE_IDX"

# --- 3. Processing loop
# tail skips header
# awk extracts ONLY the two specific columns we found, ensuring consistent orderfor the 'read' command
tail -n +2 "$INPUT_CSV" | awk -F, -v a="$ALIAS_IDX" -v b="$RUN_IDX" -v c="$SAMPLE_IDX" '{print $a "," $b "," $c}' | while IFS=, read -r alias run_id sample_id; do
        # Cleanup invisible character (important for CSVs from Excel/Windows)
	alias=$(echo "$alias" | tr -d '\r')
	run_id=$(echo "$run_id" | tr -d '\r')
	sample_id=$(echo "$sample_id" | tr -d '\r')

	if [[ -z "$alias" ]]; then continue; fi

	# Construct Path
	S3_PATH="${S3_ROOT}/${sample_id}/${run_id}/definition/output/${alias}.sorted.aligned.bam"
	LOCAL_BAM="${OUTPUT_DIR}/${alias}.sorted.aligned.bam"

	echo "Processing $alias (Run: $run_id)..."

	# Copy BAM and BAI
	aws s3 cp "$S3_PATH" "$LOCAL_BAM" --quiet
	
	if [[ $? -ne 0 ]]; then
		echo "X Error downloading $alias. Skipping."
		continue
	fi
	
	# Rheader (Add 'chr' prefix)
	echo "  -> Fixing header (adding 'chr')..."

        # Dump header
	samtools view -H "$LOCAL_BAM" > header.tmp

	# Modify header (1 -> chr1, MT -> chrM)
	sed -i 's/SN:\([0-9XY]\)/SN:chr\1/g' header.tmp
    	sed -i 's/SN:MT/SN:chrM/g' header.tmp

	# Apply new header to a temp file
	# (reheader is fast becuasue it does not decompress the whole file)
	samtools reheader header.tmp "$LOCAL_BAM" > "${LOCAL_BAM}.tmp"

	# Overwrite original
	mv "${LOCAL_BAM}.tmp" "$LOCAL_BAM"

	# Indexing
	echo "  -> Generating new index (.bai)..."
	samtools index "${LOCAL_BAM}"

	# Cleanup
	rm -f header.tmp

	echo "  ✓ Complete: ${LOCAL_BAM}"

done

echo "-------------------------------------------------------------"
echo "✓ Done! All tasks finished."
