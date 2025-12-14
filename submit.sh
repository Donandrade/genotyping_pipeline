#!/usr/bin/env bash
set -euo pipefail

source ./genotyping.conf

# detect header
HAS_HEADER=0
head -n1 "$SAMPLES_TSV" | grep -q $'\tr1\t' && HAS_HEADER=1

TOTAL_LINES=$(wc -l < "$SAMPLES_TSV")
NSAMPLES=$TOTAL_LINES
(( HAS_HEADER )) && NSAMPLES=$(( TOTAL_LINES - 1 ))

NTASKS=$(( (NSAMPLES + PER_TASK - 1) / PER_TASK ))

echo "NSAMPLES=$NSAMPLES PER_TASK=$PER_TASK => NTASKS=$NTASKS"

jid=$(sbatch --array=1-"$NTASKS"%10 genotyping_samples.sh | awk '{print $4}')
echo "Submitted samples job: $jid"

jid2=$(sbatch --dependency=afterok:"$jid" genotyping_merge.sh | awk '{print $4}')
echo "Submitted merge job (afterok): $jid2"

