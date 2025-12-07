#!/usr/bin/env bash
#SBATCH --job-name=trim_align_rg_qc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deandradesilvae@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --output=/blue/munoz/deandradesilvae/vccinoum_genomes/pipelines/genotyping_pipeline/out_trim10_%A_%a.txt
#SBATCH --array=1-12
#SBATCH --account=munoz
#SBATCH --qos=munoz-b

set -euo pipefail
pwd; hostname; date

# ===== módulos =====
module load trimmomatic   # 0.39
module load bcftools/1.22
module load bwa/0.7.17
module load samtools/1.20
module load bamutil/1.0.15
module load picard/3.2.0
# module load htslib        # para tabix

# ===== CONFIG =====

# PROBES (BED) para restringir o merge
# PROBES="simulated_probes.bed"
PROBES=""

# samples.tsv com header: sample_id  r1  r2  [rgid rglb rgpl rgpu]
SAMPLES_TSV="samples.tsv"

# Same OUTDIR of trimmed
OUTDIR="./out"

# Temporary/final BAM  (same “base” of OUTDIR)
BAM_TMP_DIR="out/bam_tmp"
BAM_FINAL_DIR="out/bam"

# Reference genome
REF="reference/subgenome_blue.multi.fa"

# Configuration for the number of threads and the number of samples to be processed
THREADS="${SLURM_CPUS_PER_TASK:-4}"
PER_TASK=256   # 12×256 = 3072 (ajuste se necessário)

# Trimmomatic
ADAPTER_PE="${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa"
TRIM_OPTS_COMMON="SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50"

# ---- pileup/merge ----
PILEUP_DIR="out/pileup"
MERGE_DIR="out/merge"

USE_PREV_PILEUPS=true   # <- se true, inclui pileups listados em PILEUP_TSV

# 
PILEUP_TSV="old_pileup.list"

# Lista de "cromossomos" (headers do FASTA)
CHR_LIST=(
  "VaccDscaff1:42640288-42650287"
  "VaccDscaff2:28801683-28811682"
  "VaccDscaff4:21204352-21214351"
  "VaccDscaff6:11534481-11544480"
  "VaccDscaff7:3282650-3292649"
  "VaccDscaff11:27652190-27662189"
  "VaccDscaff12:6795405-6805404"
  "VaccDscaff13:8612145-8622144"
  "VaccDscaff17:2834352-2844351"
  "VaccDscaff20:13144615-13154614"
  "VaccDscaff21:7405786-7415785"
  "VaccDscaff22:28175518-28185517"
)

mkdir -p "$PILEUP_DIR" "$MERGE_DIR"
mkdir -p "$OUTDIR" "$BAM_TMP_DIR" "$BAM_FINAL_DIR"

# ===== detectar header e calcular fatia =====
if head -n1 "$SAMPLES_TSV" | grep -q $'\tr1\t'; then
  HAS_HEADER=1
else
  HAS_HEADER=0
fi

TOTAL_LINES=$(wc -l < "$SAMPLES_TSV")
if (( HAS_HEADER )); then
  NSAMPLES=$(( TOTAL_LINES - 1 ))
else
  NSAMPLES=$TOTAL_LINES
fi

echo "Total de amostras (linhas no TSV): $NSAMPLES (header=$HAS_HEADER)"

START_NUM=$(( (SLURM_ARRAY_TASK_ID - 1) * PER_TASK + 1 ))
END_NUM=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
(( END_NUM > NSAMPLES )) && END_NUM=$NSAMPLES

if (( START_NUM > NSAMPLES )); then
  echo "Nada a processar nesta task ($SLURM_ARRAY_TASK_ID): START_NUM=$START_NUM > NSAMPLES=$NSAMPLES"
  # ainda assim vamos deixar o merge rodar, caso você já tenha pileups prévios
else
  echo "Task $SLURM_ARRAY_TASK_ID processará amostras $START_NUM..$END_NUM do samples.tsv"
fi

##### ==============================
##### ========== FUNÇÕES ==========
##### ==============================

trim_pair () {
  local sample="$1" r1="$2" r2="$3"
  local r1p="$OUTDIR/${sample}_R1_paired.fq.gz"
  local r1u="$OUTDIR/${sample}_R1_unpaired.fq.gz"
  local r2p="$OUTDIR/${sample}_R2_paired.fq.gz"
  local r2u="$OUTDIR/${sample}_R2_unpaired.fq.gz"

  if [[ -s "$r1p" && -s "$r2p" ]]; then
    echo "SKIP trim: $r1p / $r2p já existem"
    return 0
  fi

  trimmomatic PE -threads "$THREADS" -phred33 \
    "$r1" "$r2" \
    "$r1p" "$r1u" \
    "$r2p" "$r2u" \
    ILLUMINACLIP:${ADAPTER_PE}:2:30:10 \
    $TRIM_OPTS_COMMON
}

align_bwa () {
  local sample="$1"
  local r1p="$OUTDIR/${sample}_R1_paired.fq.gz"
  local r2p="$OUTDIR/${sample}_R2_paired.fq.gz"
  local bam_sorted="$BAM_TMP_DIR/${sample}.sorted.group.bam"

  if [[ -s "$bam_sorted" ]]; then
    echo "SKIP align: $bam_sorted já existe"
    return 0
  fi

  bwa mem -t "$THREADS" -M \
      "$REF" \
      "$r1p" \
      "$r2p" \
      -R "@RG\tID:${sample}\tLB:lib1\tPL:ILLUMINA\tPU:10K\tSM:${sample}" \
          | samtools view -hbS - \
          | samtools sort -@ "$THREADS" -o "$bam_sorted" -

  samtools index "$bam_sorted"
}

qc_bam () {
  local sample="$1"
  local bam_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam"
  samtools flagstat "$bam_rg" > "$BAM_TMP_DIR/${sample}.flagstat.txt"
  samtools stats    "$bam_rg" > "$BAM_TMP_DIR/${sample}.stats.txt"
  samtools idxstats "$bam_rg" > "$BAM_TMP_DIR/${sample}.idxstats.txt"
  bam validate --in "$bam_rg" --so_coord --verbose > "$BAM_TMP_DIR/${sample}.bamvalidate.txt" || true
}

call_variants () {
  local sample="$1"
  local bam_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam"

  local raw_vcf_gz="$BAM_FINAL_DIR/${sample}.vcf.gz"
  local norm_vcf_gz="$BAM_FINAL_DIR/${sample}_sorted_norm_split.vcf.gz"

  # mpileup + index VCF
  bcftools mpileup "$bam_rg" --annotate FORMAT/AD,FORMAT/DP -Oz \
    --min-MQ 20 -f "$REF" -o "$raw_vcf_gz"
  tabix -f -p vcf "$raw_vcf_gz"

  bcftools sort "$raw_vcf_gz" \
    | bcftools norm -O u --atomize -f "$REF" \
    | bcftools norm --multiallelics -any -f "$REF" -O z -o "$norm_vcf_gz"
  tabix -f -p vcf "$norm_vcf_gz"

  # Count for log
  echo "Total of variants — ${sample}.vcf.gz"
  bcftools view -H "$raw_vcf_gz" | wc -l || true
  echo "Normalized and split variants — ${sample}_sorted_norm_split.vcf.gz"
  bcftools view -H "$norm_vcf_gz" | wc -l || true

  # mover finais para 'PILEUP_DIR'
  mv -f "$norm_vcf_gz" "$PILEUP_DIR"
}

finalize_sample () {
  local sample="$1"
  local bam_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam"
  local bai_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam.bai"
  if [[ -s "$bam_rg" ]]; then
    mv -f "$bam_rg" "$BAM_FINAL_DIR/"
    [[ -s "$bai_rg" ]] && mv -f "$bai_rg" "$BAM_FINAL_DIR/"
    echo "Moved $sample → $BAM_FINAL_DIR"
  fi
  rm -f "$BAM_TMP_DIR/${sample}.sorted.bam" "$BAM_TMP_DIR/${sample}.sorted.bam.bai" 2>/dev/null || true
}

# ===== leitura do TSV (pula header se existir) =====
SRC_STREAM=""
if (( HAS_HEADER )); then
  SRC_STREAM=$(awk -v s="$START_NUM" -v e="$END_NUM" 'NR>1 {i=NR-1} (i>=s && i<=e)' "$SAMPLES_TSV")
else
  SRC_STREAM=$(awk -v s="$START_NUM" -v e="$END_NUM" 'NR>=s && NR<=e' "$SAMPLES_TSV")
fi

# processa o bloco selecionado de amostras (por enquanto tudo comentado)
while IFS=$'\t' read -r sample r1 r2 rgid rglb rgpl rgpu; do
  [[ -z "${sample:-}" || -z "${r1:-}" || -z "${r2:-}" ]] && { echo "Linha malformada no TSV"; exit 1; }
  [[ -f "$r1" && -f "$r2" ]] || { echo "Arquivos ausentes: $r1 / $r2"; exit 2; }

  echo "[$SLURM_ARRAY_TASK_ID] sample=$sample"
  echo "  R1: $r1"
  echo "  R2: $r2"

  # 1) Trim
  # trim_pair "$sample" "$r1" "$r2"
  # 2) Alinhamento + sort
  # align_bwa "$sample"
  # 3) QC
  # qc_bam "$sample"
  # 4) Variants
  # call_variants "$sample"
  # 5) mover BAM final e limpar intermediários
  # finalize_sample "$sample"

done <<< "$SRC_STREAM"

echo "Task $SLURM_ARRAY_TASK_ID (amostras) completed."

########################################
# MERGE POR CROMOSSOMO USANDO PROBES  #
########################################

########################################
# MERGE POR CROMOSSOMO USANDO PROBES  #
########################################

echo "[MERGE] Task $SLURM_ARRAY_TASK_ID iniciando merge por cromossomo + probes"

# mapeia ID da array -> cromossomo
CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
CHR="${CHR_LIST[$CHR_INDEX]}"

if [[ -z "${CHR:-}" ]]; then
  echo "[MERGE][ERROR] Nenhum cromossomo definido para SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
  exit 1
fi

echo "[MERGE] Cromossomo desta task: $CHR"

# coleta TODOS os pileups disponíveis
echo "[MERGE] Coletando pileups em $PILEUP_DIR..."
mapfile -t NOW_PILEUPS < <(find "$PILEUP_DIR" -maxdepth 1 -type f -name "*_sorted_norm_split.vcf.gz" | sort)
echo "[MERGE] Pileups atuais:"
printf '  %s\n' "${NOW_PILEUPS[@]}"

EXTRA_PILEUPS=()
if [[ "$USE_PREV_PILEUPS" == true ]]; then
  if [[ -f "$PILEUP_TSV" ]]; then
    while IFS=$'\t' read -r p; do
      [[ -n "$p" ]] && EXTRA_PILEUPS+=("$p")
    done < "$PILEUP_TSV"
  else
    echo "[MERGE][WARN] USE_PREV_PILEUPS=true mas PILEUP_TSV não existe: $PILEUP_TSV"
  fi
fi

ALL_PILEUPS=("${NOW_PILEUPS[@]}" "${EXTRA_PILEUPS[@]}")

if (( ${#ALL_PILEUPS[@]} == 0 )); then
  echo "[MERGE][ERROR] Nenhum *_sorted_norm_split.vcf.gz encontrado."
  exit 20
fi

# garante índices .tbi
for f in "${ALL_PILEUPS[@]}"; do
  [[ -s "${f}.tbi" ]] || tabix -f -p vcf "$f"
done

printf '%s\n' "${ALL_PILEUPS[@]}" > "$MERGE_DIR/merge.list"
echo "[MERGE] N de arquivos na merge.list: ${#ALL_PILEUPS[@]}"

SAFE_CHR=${CHR//[:\-]/_}

########################################
# LÓGICA: USAR PROBES SE POSSÍVEL      #
# CASO CONTRÁRIO, MERGE COMPLETO      #
########################################

OUT_VCF=""
USE_PROBES=false

# 1) Verifica se o arquivo de probes existe e tem conteúdo
if [[ -s "$PROBES" ]]; then
  CHR_PROBES_BED="$MERGE_DIR/probes.${SAFE_CHR}.bed"
  echo "[MERGE] Gerando BED de probes para cromossomo $CHR -> $CHR_PROBES_BED"
  awk -v c="$CHR" 'BEGIN{OFS="\t"} $1 == c {print}' "$PROBES" > "$CHR_PROBES_BED"

  if [[ -s "$CHR_PROBES_BED" ]]; then
    echo "[MERGE] Probes encontradas para $CHR. Merge será restrito às probes."
    USE_PROBES=true
    OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.probes.vcf.gz"
  else
    echo "[MERGE][WARN] Nenhuma probe encontrada para $CHR em $PROBES."
    echo "[MERGE] Merge será feito com TODO o conteúdo do cromossomo (sem -R)."
    OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.all.vcf.gz"
  fi
else
  echo "[MERGE][WARN] Arquivo de probes inexistente ou vazio: $PROBES"
  echo "[MERGE] Merge será feito com TODO o conteúdo do cromossomo (sem -R)."
  OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.all.vcf.gz"
fi

########################################
# CHAMADA DO BCFTOOLS MERGE           #
########################################

if [[ "$USE_PROBES" == true ]]; then
  echo "[MERGE] Rodando bcftools merge para $CHR usando probes em $CHR_PROBES_BED"
  bcftools merge -Oz --threads "$THREADS" -m none \
    --file-list "$MERGE_DIR/merge.list" \
    -R "$CHR_PROBES_BED" \
    -o "$OUT_VCF"
else
  # Aqui faço merge completo, mas restringindo ao cromossomo (se os VCFs tiverem todos os cromossomos)
  echo "[MERGE] Rodando bcftools merge COMPLETO para cromossomo $CHR (sem -R)"
  bcftools merge -Oz --threads "$THREADS" -m none \
    --file-list "$MERGE_DIR/merge.list" \
    -r "$CHR" \
    -o "$OUT_VCF"
fi

tabix -f -p vcf "$OUT_VCF"

echo "[MERGE] OK -> $OUT_VCF"
date

