from Bio import SeqIO
from Bio.Seq import Seq

# Input and output file paths
input_fasta = "renamed_genome.fasta"
output_fasta = "final_genome.fasta"

# Define renaming and reverse complement instructions
rename_map = {
    "sc2":  ("Chr01", True),
    "sc14": ("Chr02", False),
    "sc30": ("Chr03", False),
    "sc19": ("Chr04", True),
    "sc11": ("Chr05", False),
    "sc3":  ("Chr06", True),
    "sc18": ("Chr07", True),
    "sc9":  ("Chr08", False),
    "sc4":  ("Chr09", True),
    "sc15": ("Chr10", False),
    "sc1":  ("Chr11", True),
    "sc10": ("Chr12", True),
    "sc13": ("Chr13", True),
    "sc8":  ("Chr14", False),
    "sc5":  ("Chr15", False),
    "sc6":  ("Chr16", False),
    "sc12": ("Chr17", True),
    "sc17": ("Chr18", False)
}

# Load all sequences from FASTA
records = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

final_records = []
used_scaffolds = set()

# Apply rename + reverse complement rules
for old_name, (new_name, revcomp) in rename_map.items():
    if old_name in records:
        seq = records[old_name].seq.reverse_complement() if revcomp else records[old_name].seq
        record = records[old_name]
        record.id = new_name
        record.description = ""
        record.seq = seq
        final_records.append(record)
        used_scaffolds.add(old_name)
    else:
        print(f"⚠️ Warning: {old_name} not found in FASTA.")

# Merge sc7 and sc16 → reverse complement → name as Chr19
if "sc7" in records and "sc16" in records:
    merged_seq = (records["sc7"].seq + records["sc16"].seq).reverse_complement()
    merged_record = records["sc7"]
    merged_record.id = "Chr19"
    merged_record.description = ""
    merged_record.seq = merged_seq
    final_records.append(merged_record)
    used_scaffolds.update(["sc7", "sc16"])
else:
    print("⚠️ Warning: sc7 or sc16 not found — Chr19 not created.")

# Add all remaining scaffolds that weren't renamed or merged
for name, rec in records.items():
    if name not in used_scaffolds:
        rec.id = name  # keep original scaffold name
        rec.description = ""
        final_records.append(rec)

# Write all sequences to output FASTA
SeqIO.write(final_records, output_fasta, "fasta")

print("✅ Final genome written to:", output_fasta)

