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

# ===== modules =====
module load trimmomatic   # 0.39
module load bcftools/1.22
module load bwa/0.7.17
module load samtools/1.20
module load bamutil/1.0.15
module load picard/3.2.0
# module load htslib        # for tabix

# ===== CONFIGURATION =====

# PROBES (BED) used to:
#   - restrict bcftools mpileup (-T) when set and non-empty
#   - restrict bcctools merge per chromosome in the merge step

PROBES="probes.bed"
# PROBES=""   # <- if unset/empty or non-existent, pipeline runs genome-wide

# samples.tsv with header: sample_id  r1  r2  [rgid rglb rgpl rgpu]
SAMPLES_TSV="samples.tsv"

# Same OUTDIR used for trimmed reads
OUTDIR="./out"

# Temporary/final BAM directories (same base as OUTDIR)
BAM_TMP_DIR="out/bam_tmp"
BAM_FINAL_DIR="out/bam"

# Reports directory (read counts, flagstat summaries, timing, etc.)
REPORT_DIR="reports"
mkdir -p "$REPORT_DIR"


# Read-count report
READCOUNT_TSV="$REPORT_DIR/read_counts.tsv"
if [[ ! -s "$READCOUNT_TSV" ]]; then
  echo -e "sample\tr1_raw\tr2_raw\tr1_trimmed_paired\tr2_trimmed_paired" > "$READCOUNT_TSV"
fi

# Flagstat summary report
FLAGSTAT_TSV="$REPORT_DIR/flagstat_summary.tsv"
if [[ ! -s "$FLAGSTAT_TSV" ]]; then
  echo -e "sample\ttotal_reads\tmapped_percent\tpaired_in_sequencing\tproperly_paired_percent" > "$FLAGSTAT_TSV"
fi

# Timing reports
TIMING_SAMPLE_TSV="$REPORT_DIR/timing_samples.tsv"
if [[ ! -s "$TIMING_SAMPLE_TSV" ]]; then
  echo -e "scope\tid\tstep\tseconds\tstart_epoch\tend_epoch\ttask_id\thost" > "$TIMING_SAMPLE_TSV"
fi

TIMING_MERGE_TSV="$REPORT_DIR/timing_merge.tsv"
if [[ ! -s "$TIMING_MERGE_TSV" ]]; then
  echo -e "scope\tid\tstep\tseconds\tstart_epoch\tend_epoch\ttask_id\thost" > "$TIMING_MERGE_TSV"
fi

# Reference genome
REF="reference/subgenome_blue.multi.fa"

# Configuration for the number of threads and samples processed per array task
THREADS="${SLURM_CPUS_PER_TASK:-4}"
PER_TASK=256   # 12×256 = 3072 (adjust if needed)

# Trimmomatic
ADAPTER_PE="${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa"
TRIM_OPTS_COMMON="SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50"

# ---- pileup / merge ----
PILEUP_DIR="out/pileup"
MERGE_DIR="out/merge"
# Directory for VCF pileups split by chromosome/region
SPLIT_DIR="${PILEUP_DIR}/split_chr"

# If true, include previous pileups listed in PILEUP_TSV during merge
USE_PREV_PILEUPS=true

# TSV/list of previous pileups to include
PILEUP_TSV="old_pileup.list"

# List of "chromosomes" (FASTA headers or regions)
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

mkdir -p "$PILEUP_DIR" "$MERGE_DIR" "$SPLIT_DIR"
mkdir -p "$OUTDIR" "$BAM_TMP_DIR" "$BAM_FINAL_DIR"


# -------------------------------------
# ----- Config to count variant -------

# Expected format: CHROM <tab> SIZE
SCAFFOLD_SIZE_TSV="scaffold_size.txt"

# SNP summary table (mpileup vs call)
SNP_TABLE_TSV="$REPORT_DIR/table_snps_count_last_by_scaffold.tsv"
if [[ ! -s "$SNP_TABLE_TSV" ]]; then
    echo -e "mode\tchrom\tlen_info\tlast_pos\tn_snps" > "$SNP_TABLE_TSV"
fi



# ===== detect header and compute slice for this array task =====
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

echo "Total number of samples (lines in TSV): $NSAMPLES (header_present=$HAS_HEADER)"

START_NUM=$(( (SLURM_ARRAY_TASK_ID - 1) * PER_TASK + 1 ))
END_NUM=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
(( END_NUM > NSAMPLES )) && END_NUM=$NSAMPLES

if (( START_NUM > NSAMPLES )); then
  echo "No samples to process in this task ($SLURM_ARRAY_TASK_ID): START_NUM=$START_NUM > NSAMPLES=$NSAMPLES"
  # we still allow the merge step to run, in case previous pileups already exist
else
  echo "Task $SLURM_ARRAY_TASK_ID will process samples $START_NUM..$END_NUM from samples.tsv"
fi

##### ==============================
##### ========== FUNCTIONS =========
##### ==============================

# Helper para registrar timings
log_timing() {
  local scope="$1"   # ex: sample / merge
  local id="$2"      # ex: sample_id ou CHR
  local step="$3"    # ex: trim_pair, align_bwa, merge_chr
  local start_sec="$4"
  local end_sec="$5"
  local outfile="$6"

  local dur=$(( end_sec - start_sec ))
  local host
  host=$(hostname)

  echo -e "${scope}\t${id}\t${step}\t${dur}\t${start_sec}\t${end_sec}\t${SLURM_ARRAY_TASK_ID:-NA}\t${host}" >> "$outfile"
}

# Count the number of reads in a FASTQ (gzipped or plain)
count_fastq_reads () {
  local fq="$1"
  if [[ ! -s "$fq" ]]; then
    echo 0
    return 0
  fi

  local lines
  if [[ "$fq" == *.gz ]]; then
    lines=$(zcat "$fq" | wc -l)
  else
    lines=$(wc -l < "$fq")
  fi
  # FASTQ has 4 lines per read
  echo $(( lines / 4 ))
}

# Trim paired-end reads with Trimmomatic and record read counts before/after
trim_pair () {
  local sample="$1" r1="$2" r2="$3"
  local r1p="$OUTDIR/${sample}_R1_paired.fq.gz"
  local r1u="$OUTDIR/${sample}_R1_unpaired.fq.gz"
  local r2p="$OUTDIR/${sample}_R2_paired.fq.gz"
  local r2u="$OUTDIR/${sample}_R2_unpaired.fq.gz"

  # Run trimming only if paired outputs do not exist
  if [[ -s "$r1p" && -s "$r2p" ]]; then
    echo "SKIP trim: $r1p / $r2p already exist"
  else
    trimmomatic PE -threads "$THREADS" -phred33 \
      "$r1" "$r2" \
      "$r1p" "$r1u" \
      "$r2p" "$r2u" \
      ILLUMINACLIP:${ADAPTER_PE}:2:30:10 \
      $TRIM_OPTS_COMMON
  fi

  # Compute read counts before and after trimming
  local raw_r1 raw_r2 trim_r1 trim_r2
  raw_r1=$(count_fastq_reads "$r1")
  raw_r2=$(count_fastq_reads "$r2")
  trim_r1=$(count_fastq_reads "$r1p")
  trim_r2=$(count_fastq_reads "$r2p")

  # Avoid duplicating entries for the same sample in the TSV
  if ! grep -q "^${sample}\b" "$READCOUNT_TSV"; then
    echo -e "${sample}\t${raw_r1}\t${raw_r2}\t${trim_r1}\t${trim_r2}" >> "$READCOUNT_TSV"
  else
    echo "[REPORT] Read counts for sample ${sample} already present in ${READCOUNT_TSV}, skipping append."
  fi
}

# Align trimmed reads with BWA-MEM, sort, and index
align_bwa () {
  local sample="$1"
  local r1p="$OUTDIR/${sample}_R1_paired.fq.gz"
  local r2p="$OUTDIR/${sample}_R2_paired.fq.gz"
  local bam_sorted="$BAM_TMP_DIR/${sample}.sorted.group.bam"

  if [[ -s "$bam_sorted" ]]; then
    echo "SKIP align: $bam_sorted already exists"
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

# Basic BAM QC metrics and flagstat summary
qc_bam () {
  local sample="$1"
  local bam_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam"
  local flagstat_txt="$BAM_TMP_DIR/${sample}.flagstat.txt"

  # Standard QC outputs
  samtools flagstat "$bam_rg" > "$flagstat_txt"
  samtools stats    "$bam_rg" > "$BAM_TMP_DIR/${sample}.stats.txt"
  samtools idxstats "$bam_rg" > "$BAM_TMP_DIR/${sample}.idxstats.txt"
  bam validate --in "$bam_rg" --so_coord --verbose > "$BAM_TMP_DIR/${sample}.bamvalidate.txt" || true

  # Parse flagstat to extract summary metrics for the consolidated TSV
  local total_reads mapped_pct paired_in_seq properly_paired_pct

  total_reads=$(
    grep ' in total' "$flagstat_txt" 2>/dev/null | head -n1 | awk '{print $1}'
  )

  mapped_pct=$(
    grep ' mapped (' "$flagstat_txt" 2>/dev/null | head -n1 \
      | sed -E 's/.*\(([^[:space:]]+)%:.*/\1/' || true
  )

  paired_in_seq=$(
    grep ' paired in sequencing' "$flagstat_txt" 2>/dev/null | head -n1 | awk '{print $1}'
  )

  properly_paired_pct=$(
    grep ' properly paired (' "$flagstat_txt" 2>/dev/null | head -n1 \
      | sed -E 's/.*\(([^[:space:]]+)%:.*/\1/' || true
  )

  # Fallbacks in case patterns are not found
  total_reads=${total_reads:-0}
  mapped_pct=${mapped_pct:-0}
  paired_in_seq=${paired_in_seq:-0}
  properly_paired_pct=${properly_paired_pct:-0}

  # Append to the global flagstat summary TSV, avoiding duplicates
  if ! grep -q "^${sample}\b" "$FLAGSTAT_TSV"; then
    echo -e "${sample}\t${total_reads}\t${mapped_pct}\t${paired_in_seq}\t${properly_paired_pct}" >> "$FLAGSTAT_TSV"
  else
    echo "[REPORT] Flagstat summary for sample ${sample} already present in ${FLAGSTAT_TSV}, skipping append."
  fi
}

# Variant calling (mpileup + normalization + splitting multiallelic sites)
# If PROBES is set and non-empty, use -T PROBES to restrict mpileup to BED regions
call_variants () {
  local sample="$1"
  local bam_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam"

  local raw_vcf_gz="$BAM_FINAL_DIR/${sample}.vcf.gz"
  local norm_vcf_gz="$BAM_FINAL_DIR/${sample}_sorted_norm_split.vcf.gz"

  # Decide whether to restrict mpileup to probes
  local mpileup_region_opts=()
  if [[ -n "${PROBES:-}" && -s "$PROBES" ]]; then
    echo "[MPILEUP] Using probes file $PROBES for sample $sample (bcftools mpileup -T)"
    mpileup_region_opts=(-T "$PROBES")
  else
    echo "[MPILEUP] No valid probes file set; running genome-wide mpileup for sample $sample"
  fi

  # mpileup + index VCF
  bcftools mpileup \
    -f "$REF" \
    --annotate FORMAT/AD,FORMAT/DP \
    --min-MQ 20 \
    "${mpileup_region_opts[@]}" \
    "$bam_rg" \
    -Oz -o "$raw_vcf_gz"

  tabix -f -p vcf "$raw_vcf_gz"

  bcftools sort "$raw_vcf_gz" \
    | bcftools norm -O u --atomize -f "$REF" \
    | bcftools norm --multiallelics -any -f "$REF" -O z -o "$norm_vcf_gz"
  tabix -f -p vcf "$norm_vcf_gz"

  # Counts for logging
  echo "Total number of variants — ${sample}.vcf.gz"
  bcftools view -H "$raw_vcf_gz" | wc -l || true
  echo "Normalized and split variants — ${sample}_sorted_norm_split.vcf.gz"
  bcftools view -H "$norm_vcf_gz" | wc -l || true

  # move final per-sample normalized VCF to PILEUP_DIR
  mv -f "$norm_vcf_gz" "$PILEUP_DIR"
}

# Move final BAMs and clean up intermediate BAMs
finalize_sample () {
  local sample="$1"
  local bam_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam"
  local bai_rg="$BAM_TMP_DIR/${sample}.sorted.group.bam.bai"
  if [[ -s "$bam_rg" ]]; then
    mv -f "$bam_rg" "$BAM_FINAL_DIR/"
    [[ -s "$bai_rg" ]] && mv -f "$bai_rg" "$BAM_FINAL_DIR/"
    echo "Moved $sample -> $BAM_FINAL_DIR"
  fi
  rm -f "$BAM_TMP_DIR/${sample}.sorted.bam" "$BAM_TMP_DIR/${sample}.sorted.bam.bai" 2>/dev/null || true
}



# -------------------------------------------------------------------------
# Run bcftools call on a per-chromosome merged pileup,
# count variants before/after, and log summary statistics.
# -------------------------------------------------------------------------
call_bcftools_per_chr () {
    local chr="$1"         # chromosome/region name, e.g., VaccDscaff1:...
    local safe_chr="$2"    # sanitized string for filenames
    local in_vcf="$3"      # merged pileup VCF input
    local out_vcf="$4"     # bcftools call output

    echo "[CALL] Calling genotypes (bcftools call) for $chr"
    bcftools call -mv -Oz -o "$out_vcf" "$in_vcf"
    tabix -f -p vcf "$out_vcf"

    # ----------------------------------------------------
    # Variant counting (mpileup-style vs called genotypes)
    # ----------------------------------------------------
    echo "[COUNT] Number of variants (merged mpileup) for $chr:"
    bcftools view -H "$in_vcf" | wc -l

    echo "[COUNT] Number of variants after bcftools call for $chr:"
    bcftools view -H "$out_vcf" | wc -l

    # ----------------------------------------------------
    # SNP summary table (if scaffold_size.txt is available)
    # ----------------------------------------------------
    if [[ -f "$SCAFFOLD_SIZE_TSV" ]]; then
        local tmp_snps="$REPORT_DIR/snp_list_${safe_chr}.txt"

        # ---------- mpileup summary ----------
        bcftools query -f '%CHROM\t%POS\n' "$in_vcf" > "$tmp_snps"

        local len_info_mp
        len_info_mp=$(grep -m1 "$chr" "$SCAFFOLD_SIZE_TSV" || echo -e "$chr\tNA")

        local last_pos_mp
        last_pos_mp=$(tail -n 1 "$tmp_snps" 2>/dev/null || echo -e "$chr\t0")

        local n_snps_mp
        n_snps_mp=$(wc -l < "$tmp_snps")

        paste \
            <(echo "mpileup") \
            <(echo "$len_info_mp") \
            <(echo "$last_pos_mp") \
            <(echo "$n_snps_mp") >> "$SNP_TABLE_TSV"

        # ---------- call summary ----------
        bcftools query -f '%CHROM\t%POS\n' "$out_vcf" > "$tmp_snps"

        local len_info_call
        len_info_call=$(grep -m1 "$chr" "$SCAFFOLD_SIZE_TSV" || echo -e "$chr\tNA")

        local last_pos_call
        last_pos_call=$(tail -n 1 "$tmp_snps" 2>/dev/null || echo -e "$chr\t0")

        local n_snps_call
        n_snps_call=$(wc -l < "$tmp_snps")

        paste \
            <(echo "call") \
            <(echo "$len_info_call") \
            <(echo "$last_pos_call") \
            <(echo "$n_snps_call") >> "$SNP_TABLE_TSV"

        rm -f "$tmp_snps"

    else
        echo "[WARN] No scaffold_size.txt found; skipping SNP summary table for $chr."
    fi
}






# ===== read TSV (skip header if present) =====
SRC_STREAM=""
if (( HAS_HEADER )); then
  SRC_STREAM=$(awk -v s="$START_NUM" -v e="$END_NUM" 'NR>1 {i=NR-1} (i>=s && i<=e)' "$SAMPLES_TSV")
else
  SRC_STREAM=$(awk -v s="$START_NUM" -v e="$END_NUM" 'NR>=s && NR<=e' "$SAMPLES_TSV")
fi

# Process the selected block of samples for this array task
while IFS=$'\t' read -r sample r1 r2 rgid rglb rgpl rgpu; do
  [[ -z "${sample:-}" || -z "${r1:-}" || -z "${r2:-}" ]] && { echo "Malformed line in TSV"; exit 1; }
  [[ -f "$r1" && -f "$r2" ]] || { echo "Missing files: $r1 / $r2"; exit 2; }

  echo "[$SLURM_ARRAY_TASK_ID] sample=$sample"
  echo "  R1: $r1"
  echo "  R2: $r2"

  # 1) Trim reads + record read counts
  t_start=$(date +%s)
  trim_pair "$sample" "$r1" "$r2"
  t_end=$(date +%s)
  log_timing "sample" "$sample" "trim_pair" "$t_start" "$t_end" "$TIMING_SAMPLE_TSV"

  # 2) Alignment + sorting
  t_start=$(date +%s)
  align_bwa "$sample"
  t_end=$(date +%s)
  log_timing "sample" "$sample" "align_bwa" "$t_start" "$t_end" "$TIMING_SAMPLE_TSV"

  # 3) QC (including flagstat summary)
  t_start=$(date +%s)
  qc_bam "$sample"
  t_end=$(date +%s)
  log_timing "sample" "$sample" "qc_bam" "$t_start" "$t_end" "$TIMING_SAMPLE_TSV"

  # 4) Variant calling (with optional probes restriction)
  t_start=$(date +%s)
  call_variants "$sample"
  t_end=$(date +%s)
  log_timing "sample" "$sample" "call_variants" "$t_start" "$t_end" "$TIMING_SAMPLE_TSV"

  # 5) Move final BAM and clean intermediates
  t_start=$(date +%s)
  finalize_sample "$sample"
  t_end=$(date +%s)
  log_timing "sample" "$sample" "finalize_sample" "$t_start" "$t_end" "$TIMING_SAMPLE_TSV"

done <<< "$SRC_STREAM"

echo "Task $SLURM_ARRAY_TASK_ID (sample processing) completed."

########################################
# MERGE BY CHROMOSOME/REGION USING PROBES
########################################

echo "[MERGE] Task $SLURM_ARRAY_TASK_ID starting per-chromosome merge + probe filtering"
merge_start=$(date +%s)

# Map array ID to chromosome/region
CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
CHR="${CHR_LIST[$CHR_INDEX]}"

if [[ -z "${CHR:-}" ]]; then
  echo "[MERGE][ERROR] No chromosome/region defined for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
  exit 1
fi

echo "[MERGE] Chromosome/region for this task: $CHR"

# Safe representation of CHR for filenames
SAFE_CHR=${CHR//[:\-]/_}

########################################
# 1) SPLIT CURRENT PILEUPS BY CHR      #
########################################

echo "[MERGE] Collecting CURRENT pileups in $PILEUP_DIR..."
mapfile -t NOW_PILEUPS < <(find "$PILEUP_DIR" -maxdepth 1 -type f -name "*_sorted_norm_split.vcf.gz" | sort)

if (( ${#NOW_PILEUPS[@]} == 0 )); then
  echo "[MERGE][WARN] No current pileups found in $PILEUP_DIR"
fi

CURR_CHR_PILEUPS=()

for f in "${NOW_PILEUPS[@]}"; do
  base=$(basename "$f" .vcf.gz)   # e.g.: sampleX_sorted_norm_split
  split_vcf="$SPLIT_DIR/${base}.${SAFE_CHR}.vcf.gz"

  if [[ -s "$split_vcf" ]]; then
    echo "[MERGE][SPLIT] SKIP (already exists): $split_vcf"
  else
    echo "[MERGE][SPLIT] Extracting $CHR from $f -> $split_vcf"
    # If CHR is used as CHROM name:
    bcftools view \
      -i "CHROM==\"$CHR\"" \
      -Oz -o "$split_vcf" "$f"
    # If CHR is a genomic region instead, you would use:
    # bcftools view -r "$CHR" -Oz -o "$split_vcf" "$f"
  fi

  # Ensure index
  [[ -s "${split_vcf}.tbi" ]] || tabix -f -p vcf "$split_vcf"

  CURR_CHR_PILEUPS+=("$split_vcf")
done

echo "[MERGE] Number of CURRENT pileups split for $CHR: ${#CURR_CHR_PILEUPS[@]}"

########################################
# 2) SELECT PREVIOUS PILEUPS FOR CHR   #
########################################

PREV_CHR_PILEUPS=()

if [[ "$USE_PREV_PILEUPS" == true ]]; then
  if [[ -f "$PILEUP_TSV" ]]; then
    echo "[MERGE] Reading previous pileups from $PILEUP_TSV..."
    while IFS=$'\t' read -r p; do
      # Skip empty lines or comments
      [[ -z "$p" || "$p" =~ ^# ]] && continue

      # Heuristic: keep only files that appear to belong to this chromosome/region
      if [[ "$p" == *"$SAFE_CHR"* || "$p" == *"$CHR"* ]]; then
        PREV_CHR_PILEUPS+=("$p")
      fi
    done < "$PILEUP_TSV"
  else
    echo "[MERGE][WARN] USE_PREV_PILEUPS=true but PILEUP_TSV does not exist: $PILEUP_TSV"
  fi
fi

echo "[MERGE] Number of PREVIOUS pileups for $CHR: ${#PREV_CHR_PILEUPS[@]}"

########################################
# 3) CONCATENATE LISTS AND CHECK       #
########################################

ALL_PILEUPS_CHR=("${CURR_CHR_PILEUPS[@]}" "${PREV_CHR_PILEUPS[@]}")

if (( ${#ALL_PILEUPS_CHR[@]} == 0 )); then
  echo "[MERGE][ERROR] No VCFs (current or previous) found for $CHR."
  exit 20
fi

# Ensure index for all VCFs
for f in "${ALL_PILEUPS_CHR[@]}"; do
  [[ -s "${f}.tbi" ]] || tabix -f -p vcf "$f"
done

MERGE_LIST_CHR="$MERGE_DIR/merge.${SAFE_CHR}.list"
printf '%s\n' "${ALL_PILEUPS_CHR[@]}" > "$MERGE_LIST_CHR"
echo "[MERGE] Files included in $MERGE_LIST_CHR:"
printf '  %s\n' "${ALL_PILEUPS_CHR[@]}"

########################################
# 4) DECIDE WHETHER TO USE PROBES OR ALL
########################################

OUT_VCF=""
USE_PROBES=false

if [[ -s "$PROBES" ]]; then
  CHR_PROBES_BED="$MERGE_DIR/probes.${SAFE_CHR}.bed"
  echo "[MERGE] Generating per-chromosome BED of probes for $CHR -> $CHR_PROBES_BED"

  # It is expected that the first column of the BED matches the CHR value
  awk -v c="$CHR" 'BEGIN{OFS="\t"} $1 == c {print}' "$PROBES" > "$CHR_PROBES_BED"

  if [[ -s "$CHR_PROBES_BED" ]]; then
    echo "[MERGE] Probes found for $CHR. Merge will be restricted to these probes."
    USE_PROBES=true
    OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.probes.vcf.gz"
  else
    echo "[MERGE][WARN] No probe entries found for $CHR in $PROBES."
    echo "[MERGE] Merge will use the entire content of $CHR (no -R filter)."
    OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.all.vcf.gz"
  fi
else
  echo "[MERGE][WARN] Probe file does not exist or is empty: $PROBES"
  echo "[MERGE] Merge will use the entire content of $CHR (no -R filter)."
  OUT_VCF="$MERGE_DIR/merged.${SAFE_CHR}.all.vcf.gz"
fi

########################################
# 5) RUN BCFTOOLS MERGE                #
########################################

if [[ "$USE_PROBES" == true ]]; then
  echo "[MERGE] Running bcftools merge for $CHR using probes in $CHR_PROBES_BED"
  bcftools merge -Oz --threads "$THREADS" -m none \
    --file-list "$MERGE_LIST_CHR" \
    -R "$CHR_PROBES_BED" \
    -o "$OUT_VCF"
else
  echo "[MERGE] Running FULL bcftools merge for region $CHR (no -R filter)"
  # If needed, you could restrict to CHR with -r:
  # bcftools merge -Oz --threads "$THREADS" -m none \
  #   --file-list "$MERGE_LIST_CHR" \
  #   -r "$CHR" \
  #   -o "$OUT_VCF"
  bcftools merge -Oz --threads "$THREADS" -m none \
    --file-list "$MERGE_LIST_CHR" \
    -o "$OUT_VCF"
fi


tabix -f -p vcf "$OUT_VCF"

merge_end=$(date +%s)
log_timing "merge" "$CHR" "merge_chr" "$merge_start" "$merge_end" "$TIMING_MERGE_TSV"

# ---------------------------------------------------
# Run bcftools call on the merged per-chromosome VCF
# ---------------------------------------------------
call_start=$(date +%s)

CALLED_VCF="$MERGE_DIR/merged.${SAFE_CHR}.called.vcf.gz"
call_bcftools_per_chr "$CHR" "$SAFE_CHR" "$OUT_VCF" "$CALLED_VCF"

call_end=$(date +%s)
log_timing "merge" "$CHR" "bcftools_call_chr" "$call_start" "$call_end" "$TIMING_MERGE_TSV"

echo "[MERGE] Done -> $OUT_VCF"
echo "[CALL]  Done -> $CALLED_VCF"
date



