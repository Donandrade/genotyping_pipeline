#!/usr/bin/env bash
#SBATCH --job-name=geno_merge
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deandradesilvae@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --output=./logs/geno_merge_%A_%a.txt
#SBATCH --account=munoz
#SBATCH --qos=munoz
#SBATCH --array=1-12

set -euo pipefail
pwd; hostname; date

mkdir -p logs

module load bcftools/1.22
module load samtools/1.20

source ./genotyping.conf

mkdir -p "$PILEUP_DIR" "$SPLIT_DIR" "$MERGE_DIR" "$REPORT_DIR"

TIMING_MERGE_TSV="$REPORT_DIR/timing_merge.tsv"
[[ -s "$TIMING_MERGE_TSV" ]] || echo -e "scope\tid\tstep\tseconds\tstart_epoch\tend_epoch\ttask_id\thost" > "$TIMING_MERGE_TSV"

log_timing() {
  local scope="$1" id="$2" step="$3" start_sec="$4" end_sec="$5" outfile="$6"
  local dur=$(( end_sec - start_sec ))
  echo -e "${scope}\t${id}\t${step}\t${dur}\t${start_sec}\t${end_sec}\t${SLURM_ARRAY_TASK_ID:-NA}\t$(hostname)" >> "$outfile"
}

# Map task -> chromosome
CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
CHR="${CHR_LIST[$CHR_INDEX]:-}"

if [[ -z "${CHR:-}" ]]; then
  echo "[MERGE][ERROR] No CHR for task $SLURM_ARRAY_TASK_ID"
  exit 1
fi

SAFE_CHR=${CHR//[:\-]/_}
echo "[MERGE] CHR=$CHR SAFE_CHR=$SAFE_CHR"

merge_start=$(date +%s)

# 1) split each per-sample VCF into this CHR
mapfile -t NOW_PILEUPS < <(find "$PILEUP_DIR" -maxdepth 1 -type f -name "*_sorted_norm_split.vcf.gz" | sort)

CURR_CHR_PILEUPS=()
for f in "${NOW_PILEUPS[@]}"; do
  base=$(basename "$f" .vcf.gz)
  split_vcf="$SPLIT_DIR/${base}.${SAFE_CHR}.vcf.gz"

  if [[ ! -s "$split_vcf" ]]; then
    echo "[MERGE][SPLIT] $f -> $split_vcf"
    # Se CHR é um "contig name" literal, -i CHROM==... funciona.
    # Se CHR for uma região tipo chr:start-end, use -r "$CHR".
    bcftools view -i "CHROM==\"$CHR\"" -Oz -o "$split_vcf" "$f"
  fi
  [[ -s "${split_vcf}.tbi" ]] || tabix -f -p vcf "$split_vcf"
  CURR_CHR_PILEUPS+=("$split_vcf")
done

# 2) include previous pileups for this CHR, if configured
PREV_CHR_PILEUPS=()
if [[ "$USE_PREV_PILEUPS" == true && -f "$PILEUP_TSV" ]]; then
  while IFS=$'\t' read -r p; do
    [[ -z "$p" || "$p" =~ ^# ]] && continue
    if [[ "$p" == *"$SAFE_CHR"* || "$p" == *"$CHR"* ]]; then
      PREV_CHR_PILEUPS+=("$p")
    fi
  done < "$PILEUP_TSV"
fi

ALL_PILEUPS_CHR=("${CURR_CHR_PILEUPS[@]}" "${PREV_CHR_PILEUPS[@]}")
(( ${#ALL_PILEUPS_CHR[@]} > 0 )) || { echo "[MERGE][ERROR] No pileups for $CHR"; exit 20; }

for f in "${ALL_PILEUPS_CHR[@]}"; do
  [[ -s "${f}.tbi" ]] || tabix -f -p vcf "$f"
done

MERGE_LIST_CHR="$MERGE_DIR/merge.${SAFE_CHR}.list"
printf '%s\n' "${ALL_PILEUPS_CHR[@]}" > "$MERGE_LIST_CHR"

# 3) decide probes per CHR
USE_PROBES=false
OUT_VCF=""

if [[ -s "$PROBES" ]]; then
  CHR_PROBES_BED="$MERGE_DIR/probes.${SAFE_CHR}.bed"
  awk -v c="$CHR" 'BEGIN{OFS="\t"} $1==c {print}' "$PROBES" > "$CHR_PROBES_BED"

  if [[ -s "$CHR_PROBES_BED" ]]; then
    USE_PROBES=true
    OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.probes.vcf.gz"
  else
    OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.all.vcf.gz"
  fi
else
  OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.all.vcf.gz"
fi

# 4) merge
if [[ "$USE_PROBES" == true ]]; then
  echo "[MERGE] bcftools merge (probes) -> $OUT_VCF"
  bcftools merge -Oz --threads "$THREADS" -m none \
    --file-list "$MERGE_LIST_CHR" \
    -R "$CHR_PROBES_BED" \
    -o "$OUT_VCF"
else
  echo "[MERGE] bcftools merge (all) -> $OUT_VCF"
  bcftools merge -Oz --threads "$THREADS" -m none \
    --file-list "$MERGE_LIST_CHR" \
    -o "$OUT_VCF"
fi
tabix -f -p vcf "$OUT_VCF"

merge_end=$(date +%s)
log_timing "merge" "$CHR" "merge_chr" "$merge_start" "$merge_end" "$TIMING_MERGE_TSV"

# 5) call
CALLED_VCF="$MERGE_DIR/merged.${SAFE_CHR}.called.vcf.gz"
call_start=$(date +%s)
bcftools call -mv -Oz -o "$CALLED_VCF" "$OUT_VCF"
tabix -f -p vcf "$CALLED_VCF"
call_end=$(date +%s)
log_timing "merge" "$CHR" "bcftools_call_chr" "$call_start" "$call_end" "$TIMING_MERGE_TSV"

echo "[MERGE] Done -> $OUT_VCF"
echo "[CALL]  Done -> $CALLED_VCF"
date

