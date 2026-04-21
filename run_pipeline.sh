#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e
cd "$(dirname "$(readlink -f "$0")")"

BATCH_NAME=${1:-"trial"}
INPUT_DIR=$2

echo "======================================================"
echo " PHASE 1: PREPROCESSING & QC (NEXTFLOW) "
echo "======================================================"
# Run the Nextflow pipeline
nextflow run main.nf -profile local -resume --batch "$BATCH_NAME"

echo ""
echo "======================================================"
echo " PHASE 2: HUMAN QC REVIEW "
echo "======================================================"
echo "Nextflow QC is complete. The color-coded summary has been printed above."
echo "Your full report is located at: results/evaluation/qc_summary.csv"
echo ""
echo "ACTION REQUIRED: Please review the failures and prepare your"
echo "Snakemake input file (e.g., config/samples.tsv) with the approved samples."
echo ""

while true; do
    read -p "Have you confirmed the GUI selection? Ready to sync to AWS? (y/n): " yn
    case $yn in
        [Yy]* ) 
            echo "QC Confirmed. Proceeding to AWS Upload..."
            
            # --- TRIGGER AWS SSO LOGIN ---
            echo "Checking AWS SSO authentication..."
            aws sso login
            
            # Format the text file into a string of --include flags
            if [ -f "aws_includes.txt" ]; then
                INCLUDES=$(awk '{printf "--include \"%s\" ", $0}' aws_includes.txt)
            else
                echo "Warning: aws_includes.txt not found. Nothing to include."
                INCLUDES=""
            fi

            echo "Syncing $INPUT_DIR to S3..."
            
            # Note the order: --exclude "*" MUST come before the includes
            eval aws s3 sync "$INPUT_DIR" s3://omics-data-002313225286/clients/gtria/prod/inputs/${BATCH_NAME} --exclude \"*\" $INCLUDES
            
            echo "AWS Upload Complete!"
            break
            ;;
        [Nn]* ) 
            echo "Upload paused. Take your time. Type 'y' when ready."
            ;;
        * ) 
            echo "Please answer y (yes) or n (no)."
            ;;
    esac
done

echo ""
echo "======================================================"
echo " PIPELINE COMPLETE "
echo "======================================================"

