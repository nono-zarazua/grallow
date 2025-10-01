#!/usr/bin/env bash
set -euo pipefail

# -------------------------- #
# ONT → BAM → Clair3 → VCF
# -------------------------- #
# Requires: minimap2, samtools, bcftools, singularity
# Pull Clair3 once:
#   singularity pull $HOME/containers/clair3.sif docker://hkubal/clair3:latest

usage() {
  cat <<EOF
Usage: $0 -r REF.fa -i READS.fastq[.gz] -o OUTDIR -m MODEL_NAME
          [-s SAMPLE]
          [--mm2-threads N] [--samtools-threads N]
          [--sort-mem 4G] [--tmpdir /scratch]
Notes:
  - Relative paths can be used for the input files.
  - Runs Clair3 from the active conda env (expects run_clair3.sh on PATH).
  - MODEL_NAME examples: r941_prom_hac_g360+g422, r1041_e82_400bps_sup_v500
EOF
  exit 1
}

# -------------------------- #
# Defaults
# -------------------------- #
SAMPLE=""
GPU=0
MM2_THREADS=24          # minimap2
ST_THREADS=8           # samtools/bgzip/tabix/bcftools -@ workers
CL_THREADS="$ST_THREADS"
ST_MEM=4G              # per-thread mem for samtools sort (~(1+ST_THREADS)*ST_MEM total)
MODEL_PATH="/mnt/silo/hts2025/zarazuanav/bin/Clair3/models/"
TMPDIR_DEFAULT="${TMPDIR:-/scratch}"


# -------------------------- #
# Parse args
# -------------------------- #
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref) REF="$2"; shift 2;;
    -i|--reads) READS="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -m|--model) MODEL="$2"; shift 2;;
    -s|--sample) SAMPLE="$2"; shift 2;;
    --gpu) GPU=1; shift;;
    --mm2-threads) MM2_THREADS="$2"; shift 2;;
    --samtools-threads) ST_THREADS="$2"; CLAIR3_THREADS="$2"; shift 2;;
    --sort-mem) ST_MEM="$2"; shift 2;;
    --tmpdir) TMPDIR_DEFAULT="$2"; shift 2;;
    --primary) PRIMARY_ONLY=1; shift;;
    -h|--help) usage;;
    *) ARGS+=("$1"); shift;;
  esac
done

# -------------------------- #
# Checks
# -------------------------- #
[[ -z "${REF:-}" || -z "${READS:-}" || -z "${OUTDIR:-}" || -z "${MODEL:-}" ]] && usage
for x in minimap2 samtools bcftools run_clair3.sh; do
  command -v "$x" >/dev/null || { echo "ERROR: $x not found"; exit 2; }
done
command -v tabix >/dev/null || { echo "ERROR: tabix not found"; exit 2; }

# make absolute
REF="$(readlink -f "$REF")"
READS="$(readlink -f "$READS")"
OUTDIR="$(readlink -f "$OUTDIR")"

mkdir -p "$OUTDIR"/{logs,align,clair3,vcf}
CLAIR3_OUT="$OUTDIR/clair3"
VCF_OUT_DIR="$OUTDIR/vcf"

# sample name
if [[ -z "$SAMPLE" ]]; then
  bn="$(basename "$READS")"
  SAMPLE="${bn%%.*}"
fi

# files
MMI="${REF}.mmi"
FAIDX="${REF}.fai"
BAM="$OUTDIR/align/${SAMPLE}.minimap2.sorted.bam"
FLAGSTAT="$OUTDIR/align/${SAMPLE}.flagstat.txt"

echo ">>> [Indexing] reference: $REF"
[[ -f "$FAIDX" ]] || samtools faidx "$REF"
[[ -f "$MMI"  ]] || minimap2 -d "$MMI" "$REF"

# temp dir on local disk
# use local scratch for big temp files
TMPDIR_DEFAULT="/scratch"
TMPDIR_RUN="${TMPDIR_DEFAULT%/}/ont_${SAMPLE}_$RANDOM"
mkdir -p "$TMPDIR_RUN"
trap 'rm -rf "$TMPDIR_RUN"' EXIT


echo ">>> [Align] $READS -> $BAM"
echo "    minimap2 threads: $MM2_THREADS ; samtools sort threads: $ST_THREADS ; sort mem: $ST_MEM ; tmp: $TMPDIR_RUN"
# map-ont preset for ONT reads
minimap2 -t "$MM2_THREADS" -ax map-ont "$MMI" "$READS" 2> "$OUTDIR/logs/minimap2.err"  \
  | samtools sort -@ "$ST_THREADS" -m "$ST_MEM" -T "$TMPDIR_RUN/sort_${SAMPLE}" -o "$BAM" - 2> "$OUTDIR/logs/sort.err"
samtools index -@ "$ST_THREADS" "$BAM"

echo ">>> [QC] flagstat -> $FLAGSTAT"
samtools flagstat -@ "$ST_THREADS" "$BAM" > "$FLAGSTAT"

# Optional: restrict to primary chromosomes (GRCh38p14 NC_0000xx, X, Y, MT)
CTG_ARG=()
if [[ "$PRIMARY_ONLY" -eq 1 ]]; then
  CTGS=$(awk -F'\t' '
    ($1 ~ /^NC_0000(0[1-9]|1[0-9]|2[0-2])\./) || $1=="NC_000023.11" || $1=="NC_000024.10" || $1=="NC_012920.1" {print $1}
  ' "${FAIDX}" | paste -sd, -)
  CTG_ARG=(--ctg_name "$CTGS")
else
  CTG_ARG=(--include_all_ctgs)
fi

echo ">>> [Clair3] calling SNP/indels (model: $MODEL)"
conda activate clair3
# IMPORTANT: use absolute paths for binds
run_clair3.sh \
  --bam_fn="$BAM" \
  --ref_fn="$REF" \
  --threads="$CL_THREADS" \
  --platform="ont" \
  --model_path="$MODEL_PATH/$MODEL" \
  --include_all_ctgs \
  --output="$CLAIR3_OUT"  
  2>&1 | tee "$OUTDIR/logs/clair3.log"

echo ">>> [Post-process] bgzip/index VCF with bcftools"
VCF_RAW=""
if   [[ -s "$CLAIR3_OUT/merge_output.vcf.gz" ]]; then
  VCF_RAW="$CLAIR3_OUT/merge_output.vcf.gz"
elif [[ -s "$CLAIR3_OUT/merge_output.vcf" ]]; then
  VCF_RAW="$CLAIR3_OUT/merge_output.vcf"
else
  echo "ERROR: merged VCF not found"; exit 3
fi
echo ">>> "$VCF_RAW" will be processed"

mkdir -p "$VCF_OUT_DIR"
VCF_GZ="$VCF_OUT_DIR/${SAMPLE}.clair3.vcf.gz"
bcftools sort -@ "$ST_THREADS" -T "$TMPDIR_RUN" -Oz -o "$VCF_GZ" "$VCF_RAW"
tabix -@ "$ST_THREADS" -p vcf "$VCF_GZ"

conda deactivate

echo ">>> Done"
echo "BAM : $BAM"
echo "QC  : $FLAGSTAT"
echo "VCF : $VCF_GZ (indexed)"

