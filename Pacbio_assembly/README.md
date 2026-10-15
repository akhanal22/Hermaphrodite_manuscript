# PacBio HiFi Genome Assembly

This directory documents haplotype-resolved willow genome assembly using PacBio HiFi and Hi-C sequencing data. The workflow includes assembly with hifiasm, short-read polishing with NextPolish, Hi-C scaffolding with Juicer and 3D-DNA, and manual review in Juicebox. Reference alignments and dotplots support chromosome ordering and assembly assessment.

Scripts, polishing configurations, and raw assembly statistics are provided. See `phased_genome_assembly_willow.md` for the detailed workflow. The scripts used are inside scripts directory. Update file paths and SLURM resource settings before running.
