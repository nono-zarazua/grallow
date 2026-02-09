#!/bin/bash

source ../utils.sh

# --- CONFIGURATION ---
WORKFLOW_ID="9388837"
ROLE_ARN="arn:aws:iam::002313225286:role/HOMICS"
S3_BUCKET="omics-data-002313225286"
BASE_URI="s3://${S3_BUCKET}/clients/gtria/prod/inputs"
REFERENCES="s3://${S3_BUCKET}/clients/gtria/dev/references/GRCh38/GRCh38_genome.fa"
OUTPUT_CSV_ROOT="/home/ec2-user/workdir/project-quilt-workdir/data/bams"

# --- ARGUMENT PARSING ---
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--iuri|--input-uri)
      INPUT_URI="${2}"
      shift # past argument
      shift # past value
      ;;
    -c|--csv)
      INPUT_CSV="$2"
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: $0 --dir <s3_input_uri> --csv <s3_input_csv_file>"
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT_URI" || -z "$INPUT_CSV" ]]; then
    echo "Usage: $0 -i <input_directory_uri> -c <csv_filename>"
    exit 1
fi

# Construct full S3 paths
INPUT_CSV_S3="${BASE_URI}/${INPUT_URI}/${INPUT_CSV}"
LOCAL_CSV_FILE=$(basename "$INPUT_CSV")

# --- STEP 1: FETCH THE CSV AND FIND COLUMN INDICES ---
aws s3 cp "${INPUT_CSV_S3}" "$LOCAL_CSV_FILE" || { echo "Failed to download CSV"; exit 1; }

# Define the headers we are looking for
ALIAS_COL_NAME="alias"
BARCODE_COL_NAME="barcode"
SAMPLE_COL_NAME="sample_id"

# Find index for alias
ALIAS_IDX=$(awk -F, -v target="$ALIAS_COL_NAME" 'NR==1{for(i=1;i<=NF;i++){if($i==target){print i; exit}}}' "$LOCAL_CSV_FILE")
# Find index for barcode
BARCODE_IDX=$(awk -F, -v target="$BARCODE_COL_NAME" 'NR==1{for(i=1;i<=NF;i++){if($i==target){print i; exit}}}' "$LOCAL_CSV_FILE")
# Find index for batch
SAMPLE_IDX=$(awk -F, -v target="$SAMPLE_COL_NAME" 'NR==1{for(i=1;i<=NF;i++){if($i==target){print i; exit}}}' "$LOCAL_CSV_FILE")


# Validate we found the three
if [[ -z "$ALIAS_IDX" || -z "$BARCODE_IDX" || -z "$SAMPLE_IDX" ]]; then
    echo echo "Error: Could not find one of the required columns: alias, barcode, or sample_id"
    echo "Available columns: $(head -n 1 "$INPUT_CSV")"
    exit 1
fi

# Get Batch ID (Extract from the first data row)
# We assume the whole file belongs to the same batch/experiment.
FIRST_DATA_ROW=$(sed -n '2p' "$LOCAL_CSV_FILE")
BATCH_ID=$(echo "$FIRST_DATA_ROW" | awk -F, -v idx="$SAMPLE_IDX" '{print $idx}' | tr -d '\r')

if [[ -z "$BATCH_ID" ]]; then
    echo "Error: Could not extract Batch ID from the first row."
    exit 1
fi

# Create S3 Folder
S3_OUTPUT_URI="${BASE_URI}/${BATCH_ID}/bams"
# S3 Key (used for creating the folder marker - strip s3://bucket/)
S3_FOLDER_KEY="clients/gtria/prod/inputs/${BATCH_ID}/bams"
echo "---------------------------------------------------"
echo "Batch Setup: $BATCH_ID"
echo "Output URI:  $S3_OUTPUT_URI"
echo "Creating S3 folder marker..."

# Create the folder explicitly (0-byte object ending in /)
aws s3api put-object --bucket "$S3_BUCKET" --key "${S3_FOLDER_KEY}/" > /dev/null 2>&1 || true

echo "Found '$ALIAS_COL_NAME' at column #$ALIAS"
echo "Found '$BARCODE_COL_NAME' at column #$BARCODE_IDX"
echo "Found '$SAMPLE_COL_NAME' at column #$SAMPLE_IDX"
echo "Processing batch for $INPUT_URI"

mkdir -p ${OUTPUT_CSV_ROOT}/${BATCH_ID}
OUTPUT_CSV="${OUTPUT_CSV_ROOT}/${BATCH_ID}/batch_runs.csv"
echo "alias,barcode,run_id,sample_id" > "${OUTPUT_CSV}"

# --- STEP 2: PROCESS THE CSV ---
# We print "$alias,$barcode,sample" separated by a comma so we can read them into two variables
tail -n +2 "$LOCAL_CSV_FILE" |awk -F, -v a="$ALIAS_IDX" -v b="$BARCODE_IDX" -v c="$SAMPLE_IDX" '{print $a "," $b "," $c}' | while IFS=, read -r alias barcode sample_id; do
    # Cleanup: Remove carriage returns and spaces
    alias=$(echo "$alias" | tr -d '\r' | tr -d ' ')
    barcode=$(echo "$barcode" | tr -d '\r' | tr -d ' ')
    sample_id=$(echo "$sample_id" | tr -d '\r' | tr -d ' ')

    
    if [[ -z "$alias" || -z "$barcode" || -z "$sample_id" ]]; then continue; fi

    echo "------------------------------------------------------"
    echo "Processing Sample: $alias"
    echo "  Sample:   $sample_id"
    echo "  Barcode: $barcode"

    # --- STEP 3: CONSTRUCT PATHS ---
    # INPUT: Uses the BARCODE (e.g., barcode01.fastq.gz)
    filename="PBI25824_pass_${barcode}_263f1dba_00000000_0.fastq"
    fastq_path="${BASE_URI}/${INPUT_URI}/${barcode}/${filename}"

    # OUTPUT / NAMING: Uses the ALIAS (e.g., corrida_20251220_1)
    echo "Submitting $alias (File: $filename)..."

    # Submit to AWS Omics
    # Note: We use alias for the "sample" parameter
    run_id=$(aws omics start-run \
        --workflow-id "$WORKFLOW_ID" \
        --role-arn "$ROLE_ARN" \
        --storage-type "DYNAMIC" \
        --output-uri "$S3_OUTPUT_URI" \
        --parameters "{
            \"references\": \"$REFERENCES\",
            \"fastq\": \"$fastq_path\",
            \"sample\": \"$alias\",
            \"minimap_preset\": \"dna\"
        }" \
        --region us-east-1 \
        --output text \
        --query 'id')

    # Check if a run was correctly processed and right to CSV file
    if [[ -z "$run_id" ]]; then
        echo "  X Error submitting $alias"
    else
        echo "$alias,$barcode,$run_id,$sample_id" >> "$OUTPUT_CSV"
        echo "  ✓ Run ID: $run_id"
    fi

    # Sleep to prevent API throttling
    sleep 15
done

rm $LOCAL_CSV_FILE

echo ""
echo "Batch complete. IDs saved to batch_runs.csv"
