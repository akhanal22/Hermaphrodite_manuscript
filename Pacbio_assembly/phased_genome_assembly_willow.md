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
02-polish.sh
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
03-submit.sh
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

Now comes 3d-dna where it used the dups.txt file from juicer and polished fasta file to generate scaffold level assembly.
04-scaffolding.sh
```
#!/bin/bash
#SBATCH -J h2_3ddnar0
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128

bash ../3d-dna/run-asm-pipeline.sh -r 0 ../../HiC_assembly/references/genome.nextpolish.fasta ../../HiC_assembly/aligned/old_merged_nodups.txt
```
Now comes the tricky part, after the 3d-dna run, you will get hic and assembly file along with the genome.nextpolish.FINAL.fasta file. You would need use hic and assembly file to manually look over each scaffold in juicebox. After arranging it, run the post assembly review step of 3d-dna.
05-postreview.sh
``
#!/bin/bash
#SBATCH -J hap2_rev
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 34

bash ../3d-dna/run-asm-pipeline-post-review.sh -r ./genome.nextpolish.FINAL.review.assembly ../../HiC_assembly/references/genome.nextpolish.fasta ../../HiC_assembly/aligned/old_merged_nodups.txt
```
Also, look at the alignment with the closest and well assembled reference genome. Here I used Spurpurea version 5 reference genome and generate the maf file using lastz. Ideally you could also do it with paf through minimap but since it is not huge genome,maf would work well as well to generate dotplot.
06-lastz.sh
```
#!/bin/bash
#SBATCH -J l_hap2
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

lastz ../../../nextpolish/hap1_round3_re/Spurpurea_519_v5.0.fasta[multiple] ./genome.nextpolish.FINAL.fasta --format=MAF --chain --gapped --transition --maxwordcount=4 --exact=100 --step=20 > RESULTS_r0.maf
```
Finally I used the following python script to generate dotplots
dotplot.py
```
#!/usr/bin/env python3

import argparse
import re
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Make dotplot from MAF: X = all chromosomes, Y = HiC_scaffold_1..HiC_scaffold_40"
    )
    p.add_argument("maf", help="Input MAF file")
    p.add_argument("-o", "--output", default="dotplot_scaffold1_40.png", help="Output plot file")
    p.add_argument("--figwidth", type=float, default=16, help="Figure width")
    p.add_argument("--figheight", type=float, default=14, help="Figure height")
    p.add_argument("--dpi", type=int, default=300, help="DPI")
    p.add_argument("--linewidth", type=float, default=0.5, help="Line width")
    p.add_argument("--min_block", type=int, default=0, help="Minimum alignment block size")
    p.add_argument("--scaffold_start", type=int, default=1, help="First HiC scaffold number")
    p.add_argument("--scaffold_end", type=int, default=40, help="Last HiC scaffold number")
    return p.parse_args()


def normalize_scaffold_name(name):
    """
    Extract HiC_scaffold_N from names like:
      HiC_scaffold_1
      sample.HiC_scaffold_1
      HiC_scaffold_1:::fragment_2
    """
    m = re.search(r'(HiC_scaffold_\d+)', str(name))
    return m.group(1) if m else str(name)


def normalize_chr_name(name):
    """
    Extract chromosome names like:
      Chr01
      Chr12
      Chr15W
    from possibly decorated names.
    """
    m = re.search(r'(Chr\d+[A-Za-z]*)', str(name))
    return m.group(1) if m else str(name)


def is_scaffold(name):
    return re.search(r'HiC_scaffold_\d+', str(name)) is not None


def is_chr(name):
    return re.search(r'Chr\d+[A-Za-z]*', str(name)) is not None


def chrom_sort_key(name):
    s = str(name)

    # Chr01, Chr02, Chr12
    m = re.match(r'^Chr(\d+)$', s)
    if m:
        return (0, int(m.group(1)), "")

    # Chr15W etc
    m = re.match(r'^Chr(\d+)([A-Za-z]+)$', s)
    if m:
        return (1, int(m.group(1)), m.group(2))

    # fallback
    m = re.search(r'(\d+)', s)
    if m:
        return (2, int(m.group(1)), s)

    return (3, 10**9, s)


def scaffold_sort_key(name):
    s = normalize_scaffold_name(name)
    m = re.match(r'^HiC_scaffold_(\d+)$', s)
    if m:
        return int(m.group(1))
    return 10**9


def to_forward_coords(start, size, strand, src_len):
    if strand == "-":
        new_start = src_len - (start + size)
        new_end = src_len - start
    else:
        new_start = start
        new_end = start + size
    return new_start, new_end


def parse_maf(maf_file):
    """
    Parse pairwise MAF blocks.
    Automatically detects which s-line is chromosome and which is HiC scaffold.
    """
    alignments = []
    seq_lengths = {}
    block = []

    def process_block(lines):
        s_lines = [x for x in lines if x.startswith("s ")]
        if len(s_lines) < 2:
            return

        s1 = s_lines[0].split()
        s2 = s_lines[1].split()

        # MAF: s src start size strand srcSize sequence
        name1 = s1[1]
        start1 = int(s1[2])
        size1 = int(s1[3])
        strand1 = s1[4]
        len1 = int(s1[5])

        name2 = s2[1]
        start2 = int(s2[2])
        size2 = int(s2[3])
        strand2 = s2[4]
        len2 = int(s2[5])

        # detect which is Y scaffold and which is X chromosome
        if is_scaffold(name1) and is_chr(name2):
            y_name_raw, y_start, y_size, y_strand, y_len = name1, start1, size1, strand1, len1
            x_name_raw, x_start, x_size, x_strand, x_len = name2, start2, size2, strand2, len2
        elif is_scaffold(name2) and is_chr(name1):
            y_name_raw, y_start, y_size, y_strand, y_len = name2, start2, size2, strand2, len2
            x_name_raw, x_start, x_size, x_strand, x_len = name1, start1, size1, strand1, len1
        else:
            return

        y_name = normalize_scaffold_name(y_name_raw)
        x_name = normalize_chr_name(x_name_raw)

        y_plot_start, y_plot_end = to_forward_coords(y_start, y_size, y_strand, y_len)
        x_plot_start, x_plot_end = to_forward_coords(x_start, x_size, x_strand, x_len)

        same_orientation = (y_strand == x_strand)

        seq_lengths[y_name] = max(y_len, seq_lengths.get(y_name, 0))
        seq_lengths[x_name] = max(x_len, seq_lengths.get(x_name, 0))

        alignments.append({
            "y_name": y_name,
            "y_start": y_plot_start,
            "y_end": y_plot_end,
            "y_len": y_len,
            "x_name": x_name,
            "x_start": x_plot_start,
            "x_end": x_plot_end,
            "x_len": x_len,
            "size": y_size,
            "same_orientation": same_orientation
        })

    with open(maf_file) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                if block:
                    process_block(block)
                    block = []
                continue

            if line.startswith("a "):
                if block:
                    process_block(block)
                    block = []
                block = [line]
            else:
                block.append(line)

        if block:
            process_block(block)

    return alignments, seq_lengths


def build_offsets(names, lengths, gap):
    offsets = {}
    centers = {}
    current = 0

    for n in names:
        L = lengths.get(n, 1)
        offsets[n] = current
        centers[n] = current + L / 2
        current += L + gap

    return offsets, centers


def main():
    args = parse_args()

    alignments, seq_lengths = parse_maf(args.maf)

    if not alignments:
        sys.exit(
            "ERROR: No usable chromosome-vs-scaffold alignments were found in the MAF.\n"
            "Check that your MAF contains names like Chr01 and HiC_scaffold_1."
        )

    y_order = [f"HiC_scaffold_{i}" for i in range(args.scaffold_start, args.scaffold_end + 1)]

    # keep requested scaffolds
    alignments = [a for a in alignments if a["y_name"] in y_order]

    if args.min_block > 0:
        alignments = [a for a in alignments if a["size"] >= args.min_block]

    if not alignments:
        sys.exit(
            "ERROR: No alignments remain after filtering.\n"
            "Likely causes:\n"
            "  1) your MAF scaffold names are not matching HiC_scaffold_N\n"
            "  2) your chromosome names are not matching Chr01/Chr12 format\n"
            "  3) the chosen scaffold range has no alignments"
        )

    # use all chromosomes present in filtered data
    x_order = sorted({a["x_name"] for a in alignments}, key=chrom_sort_key)

    if not x_order:
        sys.exit("ERROR: No chromosomes found for X axis.")

    # lengths for displayed sequences
    x_lengths = {}
    y_lengths = {}

    for x in x_order:
        vals = [a["x_len"] for a in alignments if a["x_name"] == x]
        x_lengths[x] = max(vals) if vals else seq_lengths.get(x, 1)

    for y in y_order:
        vals = [a["y_len"] for a in alignments if a["y_name"] == y]
        y_lengths[y] = max(vals) if vals else seq_lengths.get(y, 1)

    x_gap = max(x_lengths.values()) * 0.03 if x_lengths else 1
    y_gap = max(y_lengths.values()) * 0.03 if y_lengths else 1

    x_offsets, x_centers = build_offsets(x_order, x_lengths, x_gap)
    y_offsets, y_centers = build_offsets(y_order, y_lengths, y_gap)

    fig, ax = plt.subplots(figsize=(args.figwidth, args.figheight))

    for a in alignments:
        x1 = x_offsets[a["x_name"]] + a["x_start"]
        x2 = x_offsets[a["x_name"]] + a["x_end"]

        if a["same_orientation"]:
            y1 = y_offsets[a["y_name"]] + a["y_start"]
            y2 = y_offsets[a["y_name"]] + a["y_end"]
            color = "red"
        else:
            y1 = y_offsets[a["y_name"]] + a["y_end"]
            y2 = y_offsets[a["y_name"]] + a["y_start"]
            color = "blue"

        ax.plot([x1, x2], [y1, y2], color=color, linewidth=args.linewidth)

    # vertical boundaries for chromosomes
    for x in x_order:
        ax.axvline(x_offsets[x], color="black", linewidth=0.4)
        ax.axvline(x_offsets[x] + x_lengths[x], color="black", linewidth=0.4)

    # horizontal boundaries for scaffolds
    for y in y_order:
        ax.axhline(y_offsets[y], color="black", linewidth=0.4)
        ax.axhline(y_offsets[y] + y_lengths[y], color="black", linewidth=0.4)

    ax.set_xticks([x_centers[x] for x in x_order])
    ax.set_xticklabels(x_order, rotation=45, ha="right")

    ax.set_yticks([y_centers[y] for y in y_order])
    ax.set_yticklabels(y_order)

    ax.set_xlabel("Target chromosomes")
    ax.set_ylabel("Query scaffolds")
    ax.set_title("Dotplot from MAF")

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved: {args.output}")

    # helpful summary
    print(f"Number of plotted alignments: {len(alignments)}")
    print("Chromosomes on X axis:")
    print(", ".join(x_order))


if __name__ == "__main__":
    main()
```
Based on your alignment, you could worked around the juicebox and keep repeating the assembly review step to get the final scaffolded genome.
