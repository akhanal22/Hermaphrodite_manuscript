#!/bin/bash
#SBATCH --job-name=hisat_map
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --array=1-35
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load HISAT2/2.2.1-gompi-2023a
module load SAMtools/1.23.1-GCC-13.3.0

THREADS=24

READDIR=cleaned_fastq_ak
OUT=hisat2_mapping_hap2
mkdir -p ${OUT}

NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_name.list | tr -d '\r' | xargs)

R1=${READDIR}/${NAME}_R1.fastq.gz
R2=${READDIR}/${NAME}_R2.fastq.gz

INDEX=hisat2_indexes/hap2_index
# For hap2 instead, use:
# INDEX=hisat2_indexes/hap2_index

echo "Mapping sample: ${NAME}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Index: ${INDEX}"

hisat2 \
  -p ${THREADS} \
  -x ${INDEX} \
  -1 ${R1} \
  -2 ${R2} \
  --dta \
  2> ${OUT}/${NAME}.hisat2.log | \
samtools sort -@ ${THREADS} -o ${OUT}/${NAME}.sorted.bam

samtools index ${OUT}/${NAME}.sorted.bam
