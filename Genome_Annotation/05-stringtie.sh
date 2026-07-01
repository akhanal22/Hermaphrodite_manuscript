#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=250G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL



#Load the modules
module load StringTie/3.0.0-GCC-13.3.0
module load SAMtools/1.23.1-GCC-13.3.0

#Define the variables
#BAMDIR=./hisat2_mapping_hap2

#Merge all the bam files
#BAM_LIST=$(ls $BAMDIR/*.sorted.bam)

#Merge all the bam files together
#samtools merge -o ./Stringtie_Ab10_hap2/hap2_AllTissues.bam $BAM_LIST
#Sort the merged bam file
#samtools sort ./Stringtie_Ab10_hap2/hap2_AllTissues.bam -o ./Stringtie_Ab10_hap2/hap2_AllTissues.s.bam

#Run the genome guided transcriptome assembly
stringtie -o Stringtie_Ab10_hap2/hap2_Stringtie_Assembly.gtf Stringtie_Ab10_hap2/hap2_AllTissues.s.bam
