#!/bin/bash
#SBATCH -J pacbio_nopostjoining
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128



hifiasm/hifiasm -o raw_assembly_again_nopostjoining/HG002.asm --h1 nextpolish/hap1_round3_re/HiCasm_next_400/fastq/HiC_527M_R1.fastq --h2 nextpolish/hap1_round3_re/HiCasm_next_400/fastq/HiC_527M_R2.fastq pacbio_sn527/hifi.fastq.gz -t 128 -u 0
