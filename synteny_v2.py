# synteny_pipeline.py

# Step 1: Extract proteins from GFF and FASTA using gffread
# Run in terminal:
# gffread nig_chr7.gff -g nig_chr7.fa -y nig_chr7.protein.fa
# gffread pur_chr7.gff -g pur_chr7.fa -y pur_chr7.protein.fa

# Step 2: Clean protein files using awk
# awk '/^>/ {h=$0; next} /^[ACDEFGHIKLMNPQRSTVWY]+$/ {print h; print; h=""}' input.fa > cleaned.fa

# Step 3: Run BLASTP
# makeblastdb -in pur_chr7.cleaned.fa -dbtype prot
# blastp -query nig_chr7.cleaned.fa -db pur_chr7.cleaned.fa -evalue 1e-5 -outfmt 6 -num_threads 8 -out nig_vs_pur.blast

#Now, let's create a gene coordinate file. For that, the bed file is created using Python. We could probably use bedtools, but I used Python for this.

python3 -m jcvi.formats.gff bed --type=gene nig_chr7.gff3 -o nig_chr7.genes.bed
#Use awk to generae tsv files from bed files.

awk 'BEGIN{OFS="\t"} {mid = int(($2 + $3) / 2); print $4, $2+1, $3, mid}' nig_chr7.genes.bed > genes.tsv

#Now this is the entire script you can create naming synteny.py or anything you would want to name. make sure to replace your file name with mine and import required modules.

import pandas as pd

# Step 1: Load gene coordinate tables
nig_genes = pd.read_csv("nig_chr7_genes.tsv", sep="\t")
dun_genes = pd.read_csv("dun_chr7_genes.tsv", sep="\t")

# Step 2: Load full BLAST results (all columns)
blast = pd.read_csv("nig_vs_dun.blast", sep="\t", header=None)
blast.columns = [
    "gene_nig_blast", "gene_dun_blast", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

# Step 3: Keep only the top (best) match per query based on highest bit score
blast = blast.sort_values("bitscore", ascending=False).drop_duplicates("gene_nig_blast", keep="first")

# Step 4: Extract core gene IDs (e.g., Sadunf07G0052100)
blast["core_nig"] = blast["gene_nig_blast"].str.extract(r"(Sadunf07G\d+)")
blast["core_dun"] = blast["gene_dun_blast"].str.extract(r"(Sadunf07G\d+)")
nig_genes["core_nig"] = nig_genes["gene"].str.extract(r"(Sadunf07G\d+)")
dun_genes["core_dun"] = dun_genes["gene"].str.extract(r"(Sadunf07G\d+)")

# Step 5: Merge BLAST with coordinate tables
merged = blast.merge(nig_genes, on="core_nig", suffixes=("", "_nig"))
merged = merged.merge(dun_genes, on="core_dun", suffixes=("", "_dun"))

# Step 6: Convert midpoints to Mb
merged["mid_nig"] = merged["mid"] / 1e6               # From nig_genes
merged["mid_dun"] = merged["mid_dun"] / 1e6           # From dun_genes

# Step 7: Filter for nig_chr7 regions of interest
focus = merged[((merged["mid_nig"] >= 1) & (merged["mid_nig"] <= 2)) |
               ((merged["mid_nig"] >= 14.7) & (merged["mid_nig"] <= 16))].copy()

# Step 8: Prepare final output
output_df = focus[["gene_nig_blast", "gene_dun_blast", "mid_nig", "mid_dun"]].copy()
output_df.columns = ["gene_nig", "gene_dun", "mid_nig", "mid_pur"]

# Step 9: Save output
output_df.to_csv("ggplot_chr07_fullcopy.tsv", sep="\t", index=False)

#This is done separately. So, the "ggplot_chr07_fullcopy.tsv" file still contains multiple matches for each gene. but i wanted gene by gene match to see where does the gene lie in different genome for different positions so i used this code and extracted the information by keeping header intact. 

(head -n 1 ggplot_chr07_fullcopy.tsv && tail -n +2 ggplot_chr07_fullcopy.tsv | awk -F'\t' '$1 == $2') > matched_gene_pairs.tsv



#After the above file is generated in python use the code below to create plot in R using ggplot. Now you can also do this in python, but i prefer ggplot so i used that
#Syntenyplot
setwd("C:/Users/khana/OneDrive - Texas Tech University/Documents/Texas_Tech/research/third_chapter/dunni_nig_alignment")
library(tidyverse)

# Load the synteny data
df <- read_tsv("matched_gene_pairs.tsv")

head(df)
# Convert positions to Mb
df <- df %>%
  mutate(mid_nig = mid_nig,
         mid_dun = mid_dun)

# Plot
ggplot(df) +
  geom_segment(aes(x = mid_nig, xend = mid_dun,
                   y = 1, yend = 0),
               color = "cornflowerblue", alpha = 0.8) +
  geom_segment(aes(x = 0, xend = max(mid_nig)), y = 1, yend = 1,
               color = "steelblue", size = 2) +
  geom_segment(aes(x = 0, xend = max(c(mid_nig, mid_dun))), y = 0, yend = 0,
               color = "indianred", size = 2) +
  scale_x_continuous(breaks = seq(1, 15, by = 1)) + 
  scale_y_continuous(breaks = c(0, 1), labels = c("dun_chr7", "nig_chr7")) +
  labs(x = "Genomic Position (Mb)", y = NULL,
       title = "Stacked Synteny Plot (Full-Copy Genes: 1–2 Mb & 14.78–16 Mb of nig_chr7)") +
  theme_minimal()


