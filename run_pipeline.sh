#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e
cd "$(dirname "$(readlink -f "$0")")"

CONFIG_DIR="config"
SAMPLES_TSV="${CONFIG_DIR}/samples.tsv"
SAMPLES_SEX_TSV="${CONFIG_DIR}/samples-sex.tsv"
CONFIG_YAML="${CONFIG_DIR}/config.yaml"

RUN_ID="trial"
RUN_NAME="corrida_${RUN_ID}"
OUTPUT_DIR="results/${RUN_NAME}"

echo "======================================================"
echo " PHASE 1: PREPROCESSING & QC (NEXTFLOW) "
echo "======================================================"
# Run the Nextflow pipeline
nextflow run main.nf -profile local -resume --outdir "$OUTPUT_DIR"

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

# Loop until the user explicitly types 'y' or 'yes'
while true; do
    read -p "Have you reviewed the QC and are you ready to launch Snakemake? (y/n): " yn
    case $yn in
        [Yy]* ) 
            echo "QC Confirmed. Proceeding to Imputation..."
            break
            ;;
        [Nn]* ) 
            echo "Pipeline paused. Take your time. Type 'y' when ready."
            ;;
        * ) 
            echo "Please answer y (yes) or n (no)."
            ;;
    esac
done

echo "Updating config.yaml with run_name: ${OUTPUT_DIR_NAME}..."
# This finds the line starting with run_name and replaces it
sed -i 's/^run_name:.*/run_name: "'"${RUN_NAME}"'"/' "$CONFIG_YAML"

echo "Generating samples.tsv and samples-sex.tsv..."
python3 /home/ec2-user/workdir/grallow/scripts/generate_snakemake_inputs.py \
    --qc-csv "${OUTPUT_DIR}/evaluation/qc_summary.csv" \
    --run-name "${OUTPUT_DIR_NAME}" \
    --out-samples "${SAMPLES_TSV}" \
    --out-sex "${SAMPLES_SEX_TSV}"

echo ""
echo "======================================================"
echo " PHASE 3: IMPUTATION (SNAKEMAKE) "
echo "======================================================"
# Trigger the Snakemake workflow 
# (Adjust cores and conda flags to match your exact environment)
#snakemake --cores all --sdm conda --keep-going --rerun-incomplete

echo ""
echo "======================================================"
echo " PIPELINE COMPLETE "
echo "======================================================"

