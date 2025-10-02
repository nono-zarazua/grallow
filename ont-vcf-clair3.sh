#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -r REF.fa[.gz] -i READS.fastq[.gz] -o OUTDIR -m MODEL_NAME
          [-s SAMPLE] [--mm2-threads N] [--samtools-threads N] [--clair3-threads N]
          [--sort-mem 4G] [--primary]
          [--model-dir /path/to/models]
Notes:
  + Relative paths can be used for the input files and output directory.
  * MODEL_NAME examples: r941_prom_hac_g360+g422, r1041_e82_400bps_sup_v500
  * Default engine is conda (since Singularity was problematic on your host).

Requirements:
  * samtools
  * minimap2
  * clair3 - via conda
  * bcftools
  * tabix
EOF
  exit 1
}

# --- require the clair3 conda env to be active ---
REQUIRED_ENV="${REQUIRED_ENV:-clair3}"   # override by exporting REQUIRED_ENV if you like
ACTIVE_ENV="${CONDA_DEFAULT_ENV:-${MAMBA_DEFAULT_ENV:-}}"

if [[ "$ACTIVE_ENV" != "$REQUIRED_ENV" ]]; then
  echo "ERROR: This script requires the '$REQUIRED_ENV' conda env to be active."
  echo "       Current env: '${ACTIVE_ENV:-<none>}'"
  echo "       Run:  conda activate $REQUIRED_ENV"
  exit 2
fi

# Also ensure Clair3 runner is on PATH (and from the right env)
if ! command -v run_clair3.sh >/dev/null 2>&1; then
  echo "ERROR: 'run_clair3.sh' not found on PATH (is '$REQUIRED_ENV' activated?)"
  exit 2
fi


# ---------- defaults ----------
SAMPLE=""
MM2_THREADS=24
ST_THREADS=8              # samtools/bgzip/tabix worker threads
CLAIR3_THREADS="$ST_THREADS"
ST_MEM=4G
TMPDIR_DEFAULT="${TMPDIR:-/scratch}"
PRIMARY_ONLY=0
MODEL_DIR="/mnt/silo/hts2025/zarazuanav/bin/Clair3/models"

# ---------- parse args ----------
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref) REF="$2"; shift 2;;
    -i|--reads) READS="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -m|--model) MODEL="$2"; shift 2;;
    -s|--sample) SAMPLE="$2"; shift 2;;
    --mm2-threads) MM2_THREADS="$2"; shift 2;;
    --samtools-threads) ST_THREADS="$2"; shift 2;;
    --clair3-threads) CLAIR3_THREADS="$2"; shift 2;;
    --sort-mem) ST_MEM="$2"; shift 2;;
    --tmpdir) TMPDIR_DEFAULT="$2"; shift 2;;
    --primary) PRIMARY_ONLY=1; shift;;
    --model-dir) MODEL_DIR="$2"; shift 2;;
    -h|--help) usage;;
    *) ARGS+=("$1"); shift;;
  esac
done
[[ -z "${REF:-}" || -z "${READS:-}" || -z "${OUTDIR:-}" || -z "${MODEL:-}" ]] && usage

# ---------- tool checks ----------
for x in minimap2 samtools bcftools tabix conda bgzip; do
  command -v "$x" >/dev/null || { echo "ERROR: $x not found"; exit 2; }
done

# ---------- normalize paths & dirs ----------
REF="$(readlink -f "$REF")"
READS="$(readlink -f "$READS")"
MODEL_DIR="$(readlink -f "$MODEL_DIR")"
[[ -d "${MODEL_DIR}/${MODEL}" ]] || { echo "ERROR: model path not found: ${MODEL_DIR}/${MODEL}"; exit 2; }
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"/{logs,align,clair3,vcf,tmp}
CLAIR3_OUT="$OUTDIR/clair3"
VCF_OUT_DIR="$OUTDIR/vcf"

# ---------- sample name ----------
if [[ -z "$SAMPLE" ]]; then
  bn="$(basename "$READS")"
  SAMPLE="${bn%.fastq.gz}"; SAMPLE="${SAMPLE%.fq.gz}"
  SAMPLE="${SAMPLE%.fastq}";  SAMPLE="${SAMPLE%.fq}"
fi

# ---------- derived file paths ----------
MMI="${REF}.mmi"
FAIDX="${REF}.fai"
BAM="$OUTDIR/align/${SAMPLE}.minimap2.sorted.bam"
FLAGSTAT="$OUTDIR/align/${SAMPLE}.flagstat.txt"

# Silence perl locale warnings
export LC_ALL=C.UTF-8 LANG=C.UTF-8

echo ">>> [Indexing] reference"
[[ -f "$FAIDX" ]] || samtools faidx "$REF"
[[ -f "$MMI"  ]] || minimap2 -d "$MMI" "$REF"

# temp working area
TMPDIR_BASE="$OUTDIR"/tmp
TMPDIR_RUN="${TMPDIR_BASE%/}/ont_${SAMPLE}_$RANDOM"
mkdir -p "$TMPDIR_RUN"
export TMPDIR="$TMPDIR_RUN"
trap 'rm -rf "$TMPDIR_RUN"' EXIT

echo ">>> [Align] $READS -> $BAM"
echo "    minimap2 t=$MM2_THREADS ; sort t=$ST_THREADS mem=$ST_MEM tmp=$TMPDIR_RUN"
minimap2 -t "$MM2_THREADS" -ax map-ont "$MMI" "$READS" 2> "$OUTDIR/logs/minimap2.err" \
  | samtools sort -@ "$ST_THREADS" -m "$ST_MEM" -T "$TMPDIR_RUN/sort_${SAMPLE}" -o "$BAM" - 2> "$OUTDIR/logs/sort.err"
samtools index -@ "$ST_THREADS" "$BAM"

echo ">>> [QC] flagstat -> $FLAGSTAT"
samtools flagstat -@ "$ST_THREADS" "$BAM" > "$FLAGSTAT"

# ---------- contig selection ----------
CTG_ARGS=()
if [[ "$PRIMARY_ONLY" -eq 1 ]]; then
  CTG_LIST=$(awk -F'\t' '
    ($1 ~ /^NC_0000(0[1-9]|1[0-9]|2[0-2])\./) || $1=="NC_000023.11" || $1=="NC_000024.10" || $1=="NC_012920.1" {print $1}
  ' "$FAIDX" | paste -sd, -)
  CTG_ARGS=(--ctg_name "$CTG_LIST")
else
  CTG_ARGS=(--include_all_ctgs)
fi

echo ">>> [Clair3] model=$MODEL"
# assume the script is run inside the clair3 conda env
run_clair3.sh \
--bam_fn="$BAM" \
--ref_fn="$REF" \
--threads="$CLAIR3_THREADS" \
--platform=ont \
--model_path="${MODEL_DIR}/${MODEL}" \
"${CTG_ARGS[@]}" \
--sample_name="$SAMPLE" \
--output="$CLAIR3_OUT" \
2>&1 | tee "$OUTDIR/logs/clair3.log"

echo ">>> [Pick merged VCF]"
VCF_RAW=""
if   [[ -s "$CLAIR3_OUT/merge_output.vcf.gz" ]]; then VCF_RAW="$CLAIR3_OUT/merge_output.vcf.gz"
elif [[ -s "$CLAIR3_OUT/merge_output.vcf"    ]]; then VCF_RAW="$CLAIR3_OUT/merge_output.vcf"
else
  echo "ERROR: merged VCF not found in $CLAIR3_OUT"; exit 3
fi

echo ">>> [Sort/bgzip/index]"
mkdir -p "$VCF_OUT_DIR"
VCF_GZ="$VCF_OUT_DIR/${SAMPLE}.clair3.vcf.gz"
# bcftools sort has no -@ ; parallelize compression with bgzip
bcftools sort -T "$TMPDIR_RUN" -Ov "$VCF_RAW" \
  | bgzip -@ "$ST_THREADS" > "$VCF_GZ"
tabix -f -@ "$ST_THREADS" -p vcf "$VCF_GZ"

echo ">>> Done"
echo "BAM : $BAM"
echo "QC  : $FLAGSTAT"
echo "VCF : $VCF_GZ (indexed)"

