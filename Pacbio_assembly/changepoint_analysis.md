1. Build gene-level 1:1 orthologs using Sapur IDs

Because both hap1 and hap2 annotation are based on S. purpurea, the Sapur gene ID is your natural “bridge”.

You want a table like:

sapur_gene_id          hap1_transcript_id      hap2_transcript_id
Sapur.007G032500       Sapur.007G032500.3.v5.1 Sapur.007G032500.5.v5.1
Sapur.007G065700       ...
...


and only for Chr07.

(A) Extract gene + transcript info for Chr07 from each haplotype GFF

Make a small tabular file from each GFF. Example (bash-ish):

# Hap1 – Chr07 only, mRNA features
awk '$1=="Chr07" && $3=="mRNA"' hap1.gff3 \
  | awk '{
      split($9, a, ";");
      id=""; parent="";
      for (i in a) {
        if (a[i] ~ /^ID=/) { sub("ID=", "", a[i]); id=a[i]; }
        if (a[i] ~ /^Parent=/) { sub("Parent=", "", a[i]); parent=a[i]; }
      }
      print id, parent, $1, $4, $5;
    }' OFS='\t' > hap1_chr7_mrna.tsv

# Hap2 – same idea
awk '$1=="Chr07" && $3=="mRNA"' hap2.gff3 \
  | awk '... same ...' > hap2_chr7_mrna.tsv


Now you have, for each hap:

transcript_id   gene_id      chr   start   end
Sapur.007G032500.1.v5.1  Sapur.007G032500  Chr07  ... ...
...

(B) Pick a single isoform per gene per haplotype

You can do this in Python or R by, e.g., picking the longest CDS or just the lowest isoform number (.1), depending how careful you want to be.

Simple approach: pick the first transcript per gene:

import pandas as pd

hap1 = pd.read_csv("hap1_chr7_mrna.tsv", sep="\t",
                   names=["tx_id", "gene_id", "chr", "start", "end"])
hap2 = pd.read_csv("hap2_chr7_mrna.tsv", sep="\t",
                   names=["tx_id", "gene_id", "chr", "start", "end"])

# Choose one transcript per gene (e.g. first; you can change to longest later)
hap1_one = hap1.sort_values("tx_id").groupby("gene_id", as_index=False).first()
hap2_one = hap2.sort_values("tx_id").groupby("gene_id", as_index=False).first()

# Keep only genes present in BOTH haplotypes
pairs = hap1_one.merge(hap2_one, on="gene_id", suffixes=("_hap1", "_hap2"))

pairs.to_csv("chr7_gene_pairs.tsv", sep="\t", index=False)


Now chr7_gene_pairs.tsv is effectively your 1:1 orthologue list for Chr07.

You actually don’t need OrthoFinder at this point, because the orthology comes from the fact that both haplotypes share the same Sapur gene ID.

2. Extract CDS for those gene pairs (Hap1 and Hap2)

If you haven’t already:

gffread hap1.gff3 -g hap1_genome.fasta -x hap1_cds.fa
gffread hap2.gff3 -g hap2_genome.fasta -x hap2_cds.fa


Now subset CDS to just the transcripts in chr7_gene_pairs.tsv.

Make a list of hap1 transcript IDs and hap2 transcript IDs:

pairs = pd.read_csv("chr7_gene_pairs.tsv", sep="\t")
pairs["tx_id_hap1"].to_csv("hap1_chr7_tx_ids.txt", index=False, header=False)
pairs["tx_id_hap2"].to_csv("hap2_chr7_tx_ids.txt", index=False, header=False)


Then:

seqkit grep -f hap1_chr7_tx_ids.txt hap1_cds.fa > hap1_chr7_pairs_cds.fa
seqkit grep -f hap2_chr7_tx_ids.txt hap2_cds.fa > hap2_chr7_pairs_cds.fa


(If you don’t have seqkit, we can do a little Python script instead.)

3. For each gene pair, make codon alignments and run Ka/Ks Calculator

Loop over each row in chr7_gene_pairs.tsv, extract its CDS from those two FASTAs, build codon alignments, and feed to Ka/Ks Calculator.

Sketch of the workflow (you can turn this into a small bash/Python pipeline):

mkdir kaks_chr7
cd kaks_chr7

# Pseudo-code-ish loop
while read gene tx1 tx2; do
    # Extract CDS for each transcript
    seqkit grep -p "$tx1" ../hap1_chr7_pairs_cds.fa > ${gene}_hap1.cds.fa
    seqkit grep -p "$tx2" ../hap2_chr7_pairs_cds.fa > ${gene}_hap2.cds.fa

    cat ${gene}_hap1.cds.fa ${gene}_hap2.cds.fa > ${gene}.cds.fa

    # Translate to protein (e.g. with transeq or Biopython)
    transeq -sequence ${gene}_hap1.cds.fa -outseq ${gene}_hap1.prot.fa
    transeq -sequence ${gene}_hap2.cds.fa -outseq ${gene}_hap2.prot.fa
    cat ${gene}_hap1.prot.fa ${gene}_hap2.prot.fa > ${gene}.prot.fa

    # Align proteins
    mafft --auto ${gene}.prot.fa > ${gene}.prot.aln.fa

    # Convert to codon alignment (PAL2NAL, PAML format)
    pal2nal.pl ${gene}.prot.aln.fa ${gene}.cds.fa -output paml > ${gene}.paml

done < <(awk 'NR>1{print $1, $2, $6}' ../chr7_gene_pairs.tsv)

Now lets run this python code to get clean file after removing stop codon for getting ds. I used python3 clean_paml.py for running this

```

python <<'EOF'
from pathlib import Path

stops = {"TAA", "TAG", "TGA"}

def read_paml(path):
    lines = [x.strip() for x in open(path) if x.strip()]
    ns, aln_len = map(int, lines[0].split()[:2])

    records = []
    i = 1

    for r in range(ns):
        name = lines[i]
        i += 1

        seq = ""
        while i < len(lines) and len(seq) < aln_len:
            seq += lines[i].replace(" ", "").upper()
            i += 1

        records.append((name, seq))

    return ns, records

for f in Path(".").glob("*.paml"):
    ns, records = read_paml(f)

    cleaned = []
    for name, seq in records:
        if seq[-3:] in stops:
            seq = seq[:-3]
        cleaned.append((name, seq))

    lengths = [len(seq) for name, seq in cleaned]

    if len(set(lengths)) != 1:
        print("BAD_LENGTH", f, lengths)
        continue

    new_len = lengths[0]

    if new_len % 3 != 0:
        print("BAD_FRAME", f, new_len)
        continue

    out = Path("paml_clean") / f.name
    with open(out, "w") as o:
        o.write(f"{ns} {new_len}\n")
        for name, seq in cleaned:
            o.write(name + "\n")
            o.write(seq + "\n")
```

Now lets run yn00 for generating ds value. I create the template for .ctl file and next step would be running this in loop.
```
cat > yn00_template.ctl << EOF
seqfile = INPUT.paml
outfile = OUTPUT.yn00
verbose = 0
icode = 0
weighting = 0
commonf3x4 = 0
```
```
for f in *.paml; do
    gene=${f%.paml}
    sed "s|INPUT.paml|$f|; s|OUTPUT.yn00|${gene}.yn00|" yn00_template.ctl > ${gene}.yn00.ctl
    yn00 ${gene}.yn00.ctl
done
```
I chcked the genes succeeded 
grep -l "Yang & Nielsen" *.yn00 | wc -l (gave me 831 total count)
For the genes failed: grep -L "Yang & Nielsen" *.yn00 > failed_yn00.txt
wc -l failed_yn00.txt (47)
