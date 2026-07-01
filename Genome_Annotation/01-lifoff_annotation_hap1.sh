#!/bin/bash
#SBATCH -J liftoff
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64

liftoff -g ./purpurea/annotation/Spurpurea_519_v5.1.gene_exons.gff3 \
        -dir ./genome_final/hap1/genome_final/ \
        -o ./genome_final/hap1/genome_final/hap1.liftoff.gff3 \
        ./genome_final/hap1/genome_final/final_genome.fasta \
        ./purpurea/Spurpurea_519_v5.0.fasta
