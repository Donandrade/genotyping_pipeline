## 1. Project Overview

This repository contains a **SLURM-oriented genotyping pipeline** designed to efficiently process large sequencing datasets by splitting the workflow into two independent parallel phases:

**Per-sample processing** (FASTQ → BAM → per-sample pileup)

**Per-chromosome** (or region) merging and genotyping

The pipeline performs:

- Read trimming

- Alignment with BWA-MEM

- Read-group assignment

- BAM post-processing and QC

- Per-sample bcftools mpileup

- Per-chromosome / per-region merging

- Optional probe-based restriction using BED files

- Variant calling and summary reports

---


## 2. Workflow Summary

**Phase 1** – Per-sample (SLURM array over samples)

For each sample:

- Trimming of raw `FASTQ` files

- Alignment with `BWA-MEM`

- Read-group assignment

- BAM sorting, indexing, and QC

- `bcftools mpileup` (genome-wide or probe-restricted)

**Phase 2** – Per-chromosome / region (SLURM array over regions)

For each chromosome, either the entire chromosome or the probe-defined region::

- Collection of all per-sample pileups

- Optional probe-based restriction

- `bcftools merge`

- Variant calling

- Per-chromosome summary statistics


## 3. Repository Organization

The workflow is intentionally divided into **independent scripts**, each with a clear responsibility:

```bash
genotyping_pipeline/
├── genotyping_samples.sh   # Phase 1: per-sample processing
├── genotyping_merge.sh     # Phase 2: per-chromosome / region merge
├── genotyping.conf         # Central configuration file
├── submit.sh               # SLURM submission wrapper
├── samples.tsv             # Sample table (input)
├── probes.bed              # Optional probe regions
├── chrom_size.txt          # Chromosome/scaffold sizes
├── old_pileup.list         # Path for the old pileups splited by chromosomes
├── reference/              # Exemple of Reference genome + indexes (Just for test)
├── fastq/                  # Exemple of FASTQ file (Just for test)
└── README.md
```

### Main workflow files

- `genotyping_samples.sh:`
Implements **Phase 1** (per-sample processing).
This script performs read **trimming**, **alignment**, **read-group assignment**, **BAM post-processing** and **QC**, and generates per-sample **pileups**.
It is designed to run as a SLURM array over samples, controlled by PER_TASK.

- `genotyping_merge.sh:`
Implements **Phase 2** (per-chromosome or per-region processing).
This script collects pileups from all samples (and optionally from previous runs), performs chromosome-wise **merging**, and runs **variant calling**.
It is designed to run as a SLURM array over chromosomes or regions, as defined in `CHR_LIST`.

- `genotyping.conf:`
Central configuration file sourced by all pipeline scripts.
It defines input files, reference paths, resource usage, output directories, probe behavior, chromosome lists, and optional reuse of previous pileups.

- `submit.sh:`
SLURM submission wrapper that orchestrates the workflow execution.
It submits the per-sample phase first, followed by the per-chromosome merge phase, ensuring the correct execution order and array configuration.

## 4. Getting started

Clone this repository into your SLURM-based HPC account before running the pipeline.

```bash
git clone https://github.com/Donandrade/genotyping_pipeline.git

cd genotyping_pipeline
```

## 5. Required Input Files and directories.

See all general description of each file and it configuration in the next topc (5. Configuration in `genotyping.sh`).

- `samples.tsv`

```bash
sample_id   r1                          r2
sample001   fastq/sample001_R1.fq.gz    fastq/sample001_R2.fq.gz
sample002   fastq/sample002_R1.fq.gz    fastq/sample002_R2.fq.gz
...
```

- `probes.bed (optional)`

```bash
VaccDscaff1	0	300	VaccDscaff1_probe001
VaccDscaff1	1000	1300	VaccDscaff1_probe002
VaccDscaff1	2000	2300	VaccDscaff1_probe003
VaccDscaff1	3000	3300	VaccDscaff1_probe004
...
```

- `chrom_size.tsv`

```bash
VaccDscaff1   42640288
VaccDscaff2   37844821
...
```

- Reference fasta + index files

```bash
reference/subgenome_blue.multi.fa
reference/subgenome_blue.multi.fa.fai
(reference must also be indexed for BWA)
```

```bash
reference/subgenome_blue.multi.fa
```

- FASTQ directory


```bash
fastq/
├── sample001_R1.fq.gz
├── sample001_R2.fq.gz
└── ...
```

## 6. Configuration (`genotyping.conf`)

All pipeline parameters are centralized in genotyping.conf. `genotyping.sh`

**Core inputs**

```bash
# ===== threads / batching =====
THREADS="${SLURM_CPUS_PER_TASK:-10}"  # Number of CPU threads per SLURM task
PER_TASK=256                          # Number of samples processed per task (Phase 1)

# ===== inputs =====
SAMPLES_TSV="samples.tsv"
REF="reference/subgenome_blue.multi.fa"
CHROM_SIZE="chrom_size.txt"
PROBES="probes.bed"   # optional; empty or unset => genome-wide processing

# ===== reuse previous pileups (merge phase) =====
USE_PREV_PILEUPS=true
PILEUP_TSV="old_pileup.list"   # Per-chromosome list of existing pileups

# ===== outputs =====
OUTDIR="./out"
TRIM_DIR="${OUTDIR}/trimmomatic"
BAM_TMP_DIR="${OUTDIR}/bam_tmp"
BAM_FINAL_DIR="${OUTDIR}/bam"
PILEUP_DIR="${OUTDIR}/pileup"
SPLIT_DIR="${PILEUP_DIR}/split_chr"
MERGE_DIR="${OUTDIR}/merge"
REPORT_DIR="reports"

```

**Input definition**

- `SAMPLES_TSV`
TSV file containing sample IDs and paired-end FASTQ paths.

- `REF`
Reference genome FASTA used for alignment and variant calling.
Must be indexed for BWA and samtools (.fai).

- `CHROM_SIZE`
Two-column TSV file (chromosome/scaffold, length) used for chromosome-aware processing and reporting.

- `PROBES` (optional)
BED file defining probe regions.
When set, both bcftools mpileup and bcftools merge are restricted to these regions.
If empty or unset, the pipeline runs genome-wide.

### Configuring the Array and `PER_TASK`

The number of array tasks (`#SBATCH --array`) and the value of `PER_TASK` must be adjusted according to the size of your dataset.

- `#SBATCH --array` defines how many tasks will run in parallel.
- `PER_TASK` defines how many samples each task will process.

Ensure that: (number of array tasks) × PER_TASK >= total number of samples.  Adjust these values according to your dataset size and available computational resources.

### Trimming parameters

```bash
ADAPTER_PE="${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa"
TRIM_OPTS_COMMON="SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:50"
```

- `ADAPTER_PE`: Adapter file used by Trimmomatic

- `TRIM_OPTS_COMMON`: Common trimming options applied to all samples 


### Reusing previous pileups (optional)

To reuse pileups generated in previous runs, set the following options in `genotyping.conf`:

```bash
USE_PREV_PILEUPS=true
PILEUP_TSV="old_pileup.list"
```

When enabled, previously generated pileups listed in `PILEUP_TSV` are included in the merge step

Useful for incremental runs or combining datasets

The file must be organized by chromosome/region

## 6. Running

```bash
bash submit.sh
```

## 7. Outputs Generated by the Pipeline

The pipeline produces the following groups of files:

### **1. Trimmed FASTQ files**
- `out/<sample>_R1_paired.fq.gz`
- `out/<sample>_R2_paired.fq.gz`

### **2. Aligned BAM files + index**
- `out/bam/<sample>.sorted.group.bam`
- `out/bam/<sample>.sorted.group.bam.bai`

### **3. QC reports (per sample)**
Located in `out/bam_tmp/`:
- `<sample>.flagstat.txt`
- `<sample>.stats.txt`
- `<sample>.idxstats.txt`
- `<sample>.bamvalidate.txt`

### **4. Per-sample VCFs**
- `out/pileup/<sample>_sorted_norm_split.vcf.gz`

### **5. Per-chromosome split pileups**
- `out/pileup/split_chr/<sample>.<CHR>.vcf.gz`

### **6. Merged VCFs (per chromosome/region)**
- `out/merge/merged.<CHR>.all.vcf.gz`  
  or  
- `out/merge/merged.<CHR>.probes.vcf.gz` *(if probes are used)*

### **7. Genotype-called VCF after merge**
- `out/merge/merged.<CHR>.called.vcf.gz`

### **8. Summary reports**
Located in `reports/`:
- `read_counts.tsv` — raw vs trimmed read counts  
- `flagstat_summary.tsv` — QC summary per sample  
- `timing_samples.tsv` and `timing_merge.tsv` — runtime tracking  
- `table_snps_count_last_by_scaffold.tsv` — SNP summary per chromosome (optional)

---

## Directory Structure Created Automatically

```bash
out/
├── trimmomatic/
│   ├── <sample>_R1_paired.fq.gz
│   ├── <sample>_R1_unpaired.fq.gz
│   ├── <sample>_R2_paired.fq.gz
│   └── <sample>_R2_unpaired.fq.gz
│
├── bam/
│   ├── <sample>.sorted.group.bam
│   └── <sample>.sorted.group.bam.bai
│
├── bam_tmp/
│   ├── <sample>.flagstat.txt
│   ├── <sample>.stats.txt
│   ├── <sample>.idxstats.txt
│   └── <sample>.bamvalidate.txt
│
├── pileup/
│   ├── <sample>_sorted_norm_split.vcf.gz
│   ├── <sample>_sorted_norm_split.vcf.gz.tbi
│   └── split_chr/
│       ├── <sample>.<CHR>.vcf.gz
│       └── <sample>.<CHR>.vcf.gz.tbi
│
└── merge/
    ├── merged.<CHR>.all.vcf.gz
    ├── merged.<CHR>.probes.vcf.gz
    ├── merged.<CHR>.called.vcf.gz
    └── *.tbi

```
