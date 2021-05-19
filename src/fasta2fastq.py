from Bio import SeqIO

for r in SeqIO.parse("../data/artic_sars-cov-2_primers.fasta", "fasta"):
    r.letter_annotations["solexa_quality"] = [40] * len(r)
    print(r.format("fastq"), end='')
