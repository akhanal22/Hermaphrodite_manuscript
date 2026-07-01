#!/bin/bash
#SBATCH -J RM_h1
#SBATCH -p highmem_p
#SBATCH -N 1
#SBATCH --mem=400G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=96:00:00
#SBATCH --array=1-2
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu

set -euo pipefail

# ---- load the modules
module load GenomeTools/1.6.5-GCC-12.3.0
module load LTR_retriever/2.9.0-foss-2022a

# sample name
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst | tr -d '\r' | xargs)
[[ -n "$file" ]] || { echo "ERROR: empty sample name for task $SLURM_ARRAY_TASK_ID"; exit 1; }


RepeatModeler \
-database "$file" \
-threads 16 \
-LTRStruct
