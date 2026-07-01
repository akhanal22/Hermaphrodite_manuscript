#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

module load HISAT2/2.2.1-gompi-2023a

mkdir -p ./hisat2_indexes

hisat2-build -p 24 hap1_softmask_customlib/hap1.softmasked.fa ./hisat2_indexes/hap1_index
hisat2-build -p 24 hap2_softmask_customlib/hap2.softmasked.fa ./hisat2_indexes/hap2_index
