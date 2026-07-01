#!/bin/bash
#SBATCH --job-name=merge_liftoff_braker
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load BEDTools/2.31.1-GCC-13.3.0

OUTDIR=./merged_liftoff_braker_hap2
mkdir -p "$OUTDIR"

LIFTOFF=./hap2.liftoff.gff3
BRAKER=./Braker_hap2/Braker/braker.gtf

echo "Checking input files..."
ls -lh "$LIFTOFF"
ls -lh "$BRAKER"

echo "Cleaning Liftoff GFF3..."

awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {next}
  NF==9 {print}
  NF!=9 {print "BAD_LIFTOFF_LINE:", NR, $0 > "bad_liftoff_lines.txt"}
' "$LIFTOFF" \
| sort -k1,1V -k4,4n -k5,5n \
> "$OUTDIR/hap2.liftoff.sorted.clean.gff3"

echo "Converting and cleaning BRAKER GTF/GFF..."

awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {next}
  NF!=9 {
    print "BAD_BRAKER_LINE:", NR, $0 > "bad_braker_lines.txt"
    next
  }

  $3=="gene" {
    gene=$9
    $9="ID="gene";Name="gene";source=BRAKER"
    print
    next
  }

  $3=="transcript" {
    tx=$9
    gene=tx
    sub(/\.t[0-9]+$/, "", gene)
    $3="mRNA"
    $9="ID="tx";Parent="gene";Name="tx";source=BRAKER"
    print
    next
  }

  {
    gene="NA"
    tx="NA"

    if (match($9, /gene_id "([^"]+)"/, a)) gene=a[1]
    if (match($9, /transcript_id "([^"]+)"/, b)) tx=b[1]

    if (tx!="NA") {
      $9="ID="tx"."$3"."$4"."$5";Parent="tx";gene_id="gene";source=BRAKER"
      print
    }
  }
' "$BRAKER" \
| sort -k1,1V -k4,4n -k5,5n \
> "$OUTDIR/hap2.braker.sorted.clean.gff3"

echo "Extracting gene features..."

awk -F'\t' '$3=="gene"' "$OUTDIR/hap2.liftoff.sorted.clean.gff3" \
> "$OUTDIR/hap2.liftoff.genes.gff3"

awk -F'\t' '$3=="gene"' "$OUTDIR/hap2.braker.sorted.clean.gff3" \
> "$OUTDIR/hap2.braker.genes.gff3"

echo "Gene counts:"
echo "Liftoff genes:"
wc -l "$OUTDIR/hap2.liftoff.genes.gff3"

echo "BRAKER genes:"
wc -l "$OUTDIR/hap2.braker.genes.gff3"

echo "Finding BRAKER genes not overlapping Liftoff genes..."

bedtools intersect \
  -a "$OUTDIR/hap2.braker.genes.gff3" \
  -b "$OUTDIR/hap2.liftoff.genes.gff3" \
  -v \
> "$OUTDIR/hap2.braker.novel_genes.gff3"

echo "Getting IDs of novel BRAKER genes..."

awk -F'\t' '
{
  if (match($9, /ID=([^;]+)/, a)) {
    print a[1]
  }
}
' "$OUTDIR/hap2.braker.novel_genes.gff3" \
> "$OUTDIR/novel_braker_gene_ids.txt"

echo "Extracting full BRAKER gene models for novel genes..."

awk -F'\t' '
BEGIN{
  OFS="\t"
  while ((getline id < "'$OUTDIR'/novel_braker_gene_ids.txt") > 0) {
    novel[id]=1
  }
}

{
  keep=0

  if (match($9, /ID=([^;]+)/, a)) {
    id=a[1]
    if (id in novel) keep=1
  }

  if (match($9, /Parent=([^;]+)/, b)) {
    parent=b[1]
    split(parent, p, ",")
    for (i in p) {
      tx=p[i]
      gene=tx
      sub(/\.t[0-9]+$/, "", gene)

      if (tx in novel) keep=1
      if (gene in novel) keep=1
    }
  }

  if (match($9, /gene_id=([^;]+)/, c)) {
    geneid=c[1]
    if (geneid in novel) keep=1
  }

  if (keep==1) print
}
' "$OUTDIR/hap2.braker.sorted.clean.gff3" \
> "$OUTDIR/hap2.braker.novel_full_models.gff3"

echo "Renaming novel BRAKER features..."

awk -F'\t' '
BEGIN{OFS="\t"}

{
  old_attr=$9

  old_gene="NA"
  old_tx="NA"

  if (match($9, /ID=([^;]+)/, a)) old_id=a[1]
  else old_id="NA"

  if (match($9, /Parent=([^;]+)/, b)) old_parent=b[1]
  else old_parent="NA"

  if ($3=="gene") {
    gene_count++
    new_gene=sprintf("novel_hap2_BRAKER_%06d", gene_count)
    gene_map[old_id]=new_gene
    $9="ID="new_gene";Name="new_gene";source=BRAKER_novel;original_id="old_id
    print
    next
  }

  if ($3=="mRNA") {
    gene=old_parent
    new_gene=gene_map[gene]
    new_tx=new_gene".t1"
    tx_map[old_id]=new_tx
    $9="ID="new_tx";Parent="new_gene";Name="new_tx";source=BRAKER_novel;original_id="old_id
    print
    next
  }

  parent=old_parent
  new_parent=tx_map[parent]

  if (new_parent=="") {
    new_parent=parent
  }

  feature_count[new_parent,$3]++
  new_id=new_parent"."$3"."feature_count[new_parent,$3]

  $9="ID="new_id";Parent="new_parent";source=BRAKER_novel;original_attr="old_attr
  print
}
' "$OUTDIR/hap2.braker.novel_full_models.gff3" \
> "$OUTDIR/hap2.braker.novel_full_models.renamed.gff3"

echo "Combining Liftoff + novel BRAKER models..."

cat "$OUTDIR/hap2.liftoff.sorted.clean.gff3" \
    "$OUTDIR/hap2.braker.novel_full_models.renamed.gff3" \
| sort -k1,1V -k4,4n -k5,5n \
> "$OUTDIR/hap2.liftoff_plus_braker_novel.final.gff3"

echo "Done."
echo "Final file:"
echo "$OUTDIR/hap2.liftoff_plus_braker_novel.final.gff3"

echo "Counts:"
echo "Liftoff clean:"
wc -l "$OUTDIR/hap2.liftoff.sorted.clean.gff3"

echo "BRAKER clean:"
wc -l "$OUTDIR/hap2.braker.sorted.clean.gff3"

echo "Liftoff genes:"
wc -l "$OUTDIR/hap2.liftoff.genes.gff3"

echo "BRAKER genes:"
wc -l "$OUTDIR/hap2.braker.genes.gff3"

echo "BRAKER novel genes:"
wc -l "$OUTDIR/hap2.braker.novel_genes.gff3"

echo "BRAKER novel gene IDs:"
wc -l "$OUTDIR/novel_braker_gene_ids.txt"

echo "BRAKER novel full models:"
wc -l "$OUTDIR/hap2.braker.novel_full_models.gff3"

echo "BRAKER novel renamed full models:"
wc -l "$OUTDIR/hap2.braker.novel_full_models.renamed.gff3"

echo "Final merged:"
wc -l "$OUTDIR/hap2.liftoff_plus_braker_novel.final.gff3"
