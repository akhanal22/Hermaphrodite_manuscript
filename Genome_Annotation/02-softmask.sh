#!/bin/bash
#SBATCH --job-name=hap_repeatmasker
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=250G
#SBATCH --array=1-2
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/ENTER/etc/profile.d/conda.sh
conda activate repeatmasker_env

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst | tr -d '\r' | xargs)

[[ -n "$sample" ]] || { echo "ERROR: empty sample name"; exit 1; }

ref="${sample}.fa"
lib="${sample}-families.fa"
threads="${SLURM_CPUS_PER_TASK:-8}"

outdir="${sample}_softmask_customlib"
mkdir -p "$outdir"

if [[ ! -s "$ref" ]]; then
  echo "ERROR: reference FASTA not found or empty: $ref" >&2
  exit 1
fi

if [[ ! -s "$lib" ]]; then
  echo "ERROR: RepeatModeler library not found or empty: $lib" >&2
  exit 1
fi

base=$(basename "$ref")

echo "[$(date)] Starting RepeatMasker on $ref"
echo "Using custom library: $lib"
echo "Threads: $threads"

RepeatMasker \
  -pa "$threads" \
  -xsmall \
  -lib "$lib" \
  -gff \
  -dir "$outdir" \
  "$ref"

if [[ -s "$outdir/${base}.masked" ]]; then
  mv -f "$outdir/${base}.masked" "$outdir/${sample}.softmasked.fa"
fi

echo "[$(date)] Done."
