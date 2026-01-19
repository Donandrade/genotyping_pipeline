#!/bin/bash
#SBATCH --job-name=nf_blueberry
#SBATCH --output=logs/nf_manager_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=$MY_EMAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=96:00:00
#SBATCH --account=munoz
#SBATCH --qos=munoz-b

# Run this sbatch script as shown below to pass the email via parameter.
## sbatch --mail-user=seu-email@ufl.edu run_pipeline.sh

# 1. Preparação
mkdir -p logs
module purge
module load nextflow

# 2. Otimização para o HiPerGator (Uso do disco local /tmp)
export NXF_TEMP=$SLURM_TMPDIR
export NXF_WORK="./work" 
export NXF_OPTS="-Xms2g -Xmx4g"

echo "========================================================"
echo "Iniciando Pipeline Nextflow"
echo "Data: $(date)"
echo "Work Dir: $PWD"
echo "========================================================"

# 3. Execução Simplificada
# O Nextflow vai ler 'params.ref' e 'params.samples' do nextflow.config

# Use "--probes false" to run the pileup for the intire chromosome 
# Use "--past_calls" to set the path for the old pileups
#    --past_calls 2025_calling/output/04_final_calls/chunks/ \
nextflow run main.nf \
    --samples samples/sample_2_last20.tsv \
    --probes probes.bed \
    --past_calls 2026_calling/output/04_final_calls/chunks/ \
    --chunk_size 10000 \
    -resume \
    -with-report logs/report.html \
    -with-timeline logs/timeline.html

echo "Finalizado: $(date)"
