
# Genome Annotation

This folder contains the scripts used for genome annotation.

1. Use Liftoff to transfer homologous gene annotations from *Salix purpurea*.

2. Use -[RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) to build a custom repeat library.

3. Use [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) with the custom repeat library to mask the genome.

4. Use [Hisat](https://github.com/DaehwanKimLab/hisat2.git) to index the genome and map the rna sequencing reads against the assembly

5. Use [Braker3](https://github.com/Gaius-Augustus/BRAKER.git) to annotate the genome

6. Merge braker and liftoff results to get final gff file. If braker and liftoff annotate the same gene, its gets gene name based on Salix purpurea gene names. Novel genes annotated by braker are named as novel genes. 
