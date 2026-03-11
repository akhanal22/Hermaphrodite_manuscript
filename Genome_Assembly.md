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
