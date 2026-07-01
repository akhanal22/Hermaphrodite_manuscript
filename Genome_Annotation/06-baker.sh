#!/bin/bash
#SBATCH --job-name=braker
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=250G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL


module purge

#load the modules
module load BRAKER/3.0.8-foss-2023a
module load SciPy-bundle/2024.05-gfbf-2024a
module load Biopython/1.84-gfbf-2024a
module load GeneMark-ETP/1.0.0-GCC-12.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load gffread/0.12.7-GCCcore-12.3.0
module load StringTie/3.0.0-GCC-13.3.0

#Define the variables

REF=/scratch/ak18787/hifi_project/hap2_softmask_customlib/hap2.softmasked.fa
BAM=/scratch/ak18787/hifi_project/hisat2_mapping_hap2

#Make the output directory 
#cd $OUT
#mkdir $OUT/Braker
cd Braker_hap2

#List the sorted bam files from the Hisat2 directory
LIST=$(ls $BAM/*.sorted.bam | tr '\n' ,)

#Download the orthodb protein file for plants
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
gunzip Viridiplantae.fa.gz

#Run Braker to annotate the genome.
braker.pl \
--genome $REF \
--bam $LIST \
--prot_seq Viridiplantae.fa \
--threads 40 \
--workingdir Braker/ \
--nocleanup
