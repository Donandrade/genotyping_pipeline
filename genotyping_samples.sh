#!/usr/bin/env bash
#SBATCH --job-name=geno_samples
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deandradesilvae@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --output=./logs/geno_samples_%A_%a.txt
#SBATCH --account=munoz
#SBATCH --qos=munoz
# NOTE: --array será passado no sbatch (submit.sh)

set -euo pipefail
pwd; hostname; date

mkdir -p logs

# ===== modules =====
module load trimmomatic
module load bcftools/1.22
module load bwa/0.7.17
module load samtools/1.20
module load bamutil/1.0.15
module load picard/3.2.0

# ===== config =====
source ./genotyping.conf

# agora a variável existe no ambiente
ADAPTER_PE="${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa"
TRIM_OPTS_COMMON="SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50"

# ===== dirs =====
mkdir -p "$OUTDIR" "$TRIM_DIR" "$BAM_TMP_DIR" "$BAM_FINAL_DIR" "$PILEUP_DIR" "$REPORT_DIR"

READCOUNT_TSV="$REPORT_DIR/read_counts.tsv"
FLAGSTAT_TSV="$REPORT_DIR/flagstat_summary.tsv"
TIMING_SAMPLE_TSV="$REPORT_DIR/timing_samples.tsv"

[[ -s "$READCOUNT_TSV" ]]     || echo -e "sample\tr1_raw\tr2_raw\tr1_trimmed_paired\tr2_trimmed_paired" > "$READCOUNT_TSV"
[[ -s "$FLAGSTAT_TSV" ]]      || echo -e "sample\ttotal_reads\tmapped_percent\tpaired_in_sequencing\tproperly_paired_percent" > "$FLAGSTAT_TSV"
[[ -s "$TIMING_SAMPLE_TSV" ]] || echo -e "scope\tid\tstep\tseconds\tstart_epoch\tend_epoch\ttask_id\thost" > "$TIMING_SAMPLE_TSV"

log_timing() {
  local scope="$1" id="$2" step="$3" start_sec="$4" end_sec="$5" outfile="$6"
  local dur=$(( end_sec - start_sec ))
  echo -e "${scope}\t${id}\t${step}\t${dur}\t${start_sec}\t${end_sec}\t${SLURM_ARRAY_TASK_ID:-NA}\t$(hostname)" >> "$outfile"
}

count_fastq_reads () {
  local fq="$1"
  [[ -s "$fq" ]] || { echo 0; return 0; }
  local lines
  if [[ "$fq" == *.gz ]]; then lines=$(zcat "$fq" | wc -l); else lines=$(wc -l < "$fq"); fi
  echo $(( lines / 4 ))
}

trim_pair () {
  local sample="$1" r1="$2" r2="$3"
  local r1p="$TRIM_DIR/${sample}_R1_paired.fq.gz"
  local r1u="$TRIM_DIR/${sample}_R1_unpaired.fq.gz"
  local r2p="$TRIM_DIR/${sample}_R2_paired.fq.gz"
  local r2u="$TRIM_DIR/${sample}_R2_unpaired.fq.gz"

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

  local raw_r1 raw_r2 trim_r1 trim_r2
  raw_r1=$(count_fastq_reads "$r1")
  raw_r2=$(count_fastq_reads "$r2")
  trim_r1=$(count_fastq_reads "$r1p")
  trim_r2=$(count_fastq_reads "$r2p")

  grep -q "^${sample}\b" "$READCOUNT_TSV" || echo -e "${sample}\t${raw_r1}\t${raw_r2}\t${trim_r1}\t${trim_r2}" >> "$READCOUNT_TSV"
}

align_bwa () {
  local sample="$1"
  local r1p="$TRIM_DIR/${sample}_R1_paired.fq.gz"
  local r2p="$TRIM_DIR/${sample}_R2_paired.fq.gz"
  local bam_sorted="$BAM_TMP_DIR/${sample}.sorted.group.bam"

  [[ -s "$bam_sorted" ]] && { echo "SKIP align: $bam_sorted exists"; return 0; }

  bwa mem -t "$THREADS" -M \
    "$REF" "$r1p" "$r2p" \
    -R "@RG\tID:${sample}\tLB:lib1\tPL:ILLUMINA\tPU:10K\tSM:${sample}" \
    | samtools view -hbS - \
    | samtools sort -@ "$THREADS" -o "$bam_sorted" -

  samtools index "$bam_sorted"
}

qc_bam () {
  local sample="$1"
  local bam="$BAM_TMP_DIR/${sample}.sorted.group.bam"
  local flagstat_txt="$BAM_TMP_DIR/${sample}.flagstat.txt"

  samtools flagstat "$bam" > "$flagstat_txt"
  samtools stats "$bam"    > "$BAM_TMP_DIR/${sample}.stats.txt"
  samtools idxstats "$bam" > "$BAM_TMP_DIR/${sample}.idxstats.txt"
  bam validate --in "$bam" --so_coord --verbose > "$BAM_TMP_DIR/${sample}.bamvalidate.txt" || true

  local total_reads mapped_pct paired_in_seq properly_paired_pct
  total_reads=$(grep ' in total' "$flagstat_txt" | head -n1 | awk '{print $1}' || echo 0)
  mapped_pct=$(grep ' mapped (' "$flagstat_txt" | head -n1 | sed -E 's/.*\(([^[:space:]]+)%:.*/\1/' || echo 0)
  paired_in_seq=$(grep ' paired in sequencing' "$flagstat_txt" | head -n1 | awk '{print $1}' || echo 0)
  properly_paired_pct=$(grep ' properly paired (' "$flagstat_txt" | head -n1 | sed -E 's/.*\(([^[:space:]]+)%:.*/\1/' || echo 0)

  grep -q "^${sample}\b" "$FLAGSTAT_TSV" || echo -e "${sample}\t${total_reads}\t${mapped_pct}\t${paired_in_seq}\t${properly_paired_pct}" >> "$FLAGSTAT_TSV"
}

call_variants () {
  local sample="$1"
  local bam="$BAM_TMP_DIR/${sample}.sorted.group.bam"

  local raw_vcf="$PILEUP_DIR/${sample}.vcf.gz"
  local norm_vcf="$PILEUP_DIR/${sample}_sorted_norm_split.vcf.gz"

  local mpileup_opts=()
  if [[ -n "${PROBES:-}" && -s "$PROBES" ]]; then
    mpileup_opts=(-T "$PROBES")
  fi

  bcftools mpileup \
    -f "$REF" --annotate FORMAT/AD,FORMAT/DP --min-MQ 20 \
    "${mpileup_opts[@]}" \
    "$bam" -Oz -o "$raw_vcf"
  tabix -f -p vcf "$raw_vcf"

  bcftools sort "$raw_vcf" \
    | bcftools norm -O u --atomize -f "$REF" \
    | bcftools norm --multiallelics -any -f "$REF" -O z -o "$norm_vcf"
  tabix -f -p vcf "$norm_vcf"
}

finalize_sample () {
  local sample="$1"
  local bam="$BAM_TMP_DIR/${sample}.sorted.group.bam"
  local bai="$BAM_TMP_DIR/${sample}.sorted.group.bam.bai"

  [[ -s "$bam" ]] && mv -f "$bam" "$BAM_FINAL_DIR/"
  [[ -s "$bai" ]] && mv -f "$bai" "$BAM_FINAL_DIR/"
}

# ===== detect header + slice samples for this task =====
HAS_HEADER=0
head -n1 "$SAMPLES_TSV" | grep -q $'\tr1\t' && HAS_HEADER=1

TOTAL_LINES=$(wc -l < "$SAMPLES_TSV")
NSAMPLES=$TOTAL_LINES
(( HAS_HEADER )) && NSAMPLES=$(( TOTAL_LINES - 1 ))

START_NUM=$(( (SLURM_ARRAY_TASK_ID - 1) * PER_TASK + 1 ))
END_NUM=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
(( END_NUM > NSAMPLES )) && END_NUM=$NSAMPLES

if (( START_NUM > NSAMPLES )); then
  echo "[SAMPLES] No samples for task $SLURM_ARRAY_TASK_ID (START_NUM=$START_NUM > NSAMPLES=$NSAMPLES)"
  exit 0
fi

echo "[SAMPLES] Task $SLURM_ARRAY_TASK_ID processes samples $START_NUM..$END_NUM (NSAMPLES=$NSAMPLES)"

if (( HAS_HEADER )); then
  SRC_STREAM=$(awk -v s="$START_NUM" -v e="$END_NUM" 'NR>1 {i=NR-1} (i>=s && i<=e)' "$SAMPLES_TSV")
else
  SRC_STREAM=$(awk -v s="$START_NUM" -v e="$END_NUM" 'NR>=s && NR<=e' "$SAMPLES_TSV")
fi

while IFS=$'\t' read -r sample r1 r2 rgid rglb rgpl rgpu; do
  [[ -z "${sample:-}" || -z "${r1:-}" || -z "${r2:-}" ]] && { echo "Malformed TSV line"; exit 1; }
  [[ -f "$r1" && -f "$r2" ]] || { echo "Missing files: $r1 / $r2"; exit 2; }

  echo "[SAMPLES] sample=$sample"

  t0=$(date +%s); trim_pair "$sample" "$r1" "$r2"; t1=$(date +%s); log_timing "sample" "$sample" "trim_pair" "$t0" "$t1" "$TIMING_SAMPLE_TSV"
  t0=$(date +%s); align_bwa "$sample";                 t1=$(date +%s); log_timing "sample" "$sample" "align_bwa" "$t0" "$t1" "$TIMING_SAMPLE_TSV"
  t0=$(date +%s); qc_bam "$sample";                    t1=$(date +%s); log_timing "sample" "$sample" "qc_bam" "$t0" "$t1" "$TIMING_SAMPLE_TSV"
  t0=$(date +%s); call_variants "$sample";             t1=$(date +%s); log_timing "sample" "$sample" "call_variants" "$t0" "$t1" "$TIMING_SAMPLE_TSV"
  t0=$(date +%s); finalize_sample "$sample";           t1=$(date +%s); log_timing "sample" "$sample" "finalize_sample" "$t0" "$t1" "$TIMING_SAMPLE_TSV"

done <<< "$SRC_STREAM"

echo "[SAMPLES] Done."

