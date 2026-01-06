#!/bin/bash
#SBATCH --job-name=nf_blueberry
#SBATCH --output=logs/nf_manager_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=youremail@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=96:00:00
#SBATCH --account=munoz
#SBATCH --qos=munoz-b

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
nextflow run main.nf \
     \
    -resume \
    -with-report logs/report.html \
    -with-timeline logs/timeline.html

echo "Finalizado: $(date)"
