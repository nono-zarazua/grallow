#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e
cd "$(dirname "$(readlink -f "$0")")"

CONFIG_DIR="config"
SAMPLES_TSV="${CONFIG_DIR}/samples.tsv"
SAMPLES_SEX_TSV="${CONFIG_DIR}/samples-sex.tsv"
CONFIG_YAML="${CONFIG_DIR}/config.yaml"

echo "======================================================"
echo " PHASE 1: PREPROCESSING & QC (NEXTFLOW) "
echo "======================================================"
# Run the Nextflow pipeline
nextflow run main.nf -profile local -resume

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
echo "======================================================"
echo " QC COMPLETE "
echo "======================================================"
