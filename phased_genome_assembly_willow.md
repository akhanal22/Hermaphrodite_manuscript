We start with raw assembly of hifi reads. For this, I used hifiasm and used -u 0 for no postjoining step. 
01-pacbio_assembly.sh
```
#!/bin/bash
#SBATCH -J pacbio_nopostjoining
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128



hifiasm/hifiasm -o raw_assembly_again_nopostjoining/HG002.asm --h1 nextpolish/hap1_round3_re/HiCasm_next_400/fastq/HiC_527M_R1.fastq --h2 nextpolish/hap1_round3_
re/HiCasm_next_400/fastq/HiC_527M_R2.fastq pacbio_sn527/hifi.fastq.gz -t 128 -u 0
```
After the hifiasm finished run, I converted the raw gfa to fa using awk: awk '/^S\t/ {print ">"$2"\n"$3}' input.gfa (use whatever your input file is, it my case i did it 
for both haplotypes) > output.fa (I named it as hap1.fa and hap2.fa for each haplotypes) for both haplotypes. Now, for hap1.fa and hap2.fa, I polished these genomes with short reads using
nextDenovo (02-nextpolish_hap1.sh and 03-02-nextpolish_hap2.sh). I only show one haplotype code below because it is the same code with just different haplotype. The codes are also
separately saved in this repository. 
```
#!/bin/bash
#SBATCH -J polish1
#SBATCH -o log/%x.o%j
#SBATCH -e log/%x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128


nextPolish nextpolish_hap1.cfg
```
Ater the raw genome was polished, I map HiC reads against these genome to get the aligned merge_nodups.txt file. We use juicer for that. For Juicer, download it from aiden lab simply using git clone. The directory structure for juicer should be HiCassembly (main directory) and within that we have sub directory: scripts (scripts to run juicer), restriction_sites (from nanopore assembly with generate_site_positions.py),fastq (paired end hic fastq files), references (your genome polished with nextpolish), splits (fastq split files) and a file chrom.sizes that has contig name and sizes. 
For retsriction_sites folder, generate file including restrcition enzyme cutting sites on our raw assembled genome. The python script for thst is under [juicer_directory]/misc named generate_size_positions.py. 
```
#!/bin/bash
#SBATCH -J site_position
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

python ../../juicer/misc/generate_site_positions.py MboI SANI ../references/genome.nextpolish.fasta
```
Here, MboI is restrcition enzyme used to digest the genome. For chrom.sizes, I used Nan's code (seqlength.py) which I have uploaded here as well
python seqlength.py -f polished nanopore assemble fasta > [output file name]. 
For splits folder, I ran the following submission scripts to generate fastq split files:
```
#!/bin/bash
#SBATCH -J splits
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128


split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/HiC_527M_R1.fastq
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/HiC_527M_R2.fastq
```
Finally, after having all files and folder ready, go inside directory, creata submission script named submit.sh and run the juicer.sh Oops also, make sure to index your genome using "bwa index"
```
#!/bin/bash
#SBATCH -J hap1
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=250G


bash juicer.sh -d ../ -z ../references/genome.nextpolish.fasta  -y ../restriction_sites/SANI_MboI.txt -t 64 -p ../chrom.sizes  -a bwa --assembly -S early
```
I used -S early, because i kept error earlier and since I only need merged_nodups.txt file for running 3d-dna, it worked perfectly fine for me. 









































