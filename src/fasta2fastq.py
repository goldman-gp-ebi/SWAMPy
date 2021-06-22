# All this file does is add dummy quality scores to a fasta file to produce a fastq. 
# It's not used in the main code - but it can be handy if you want to use a different set of primers from ARTIC V3. 

from Bio import SeqIO

for r in SeqIO.parse("../data/artic_sars-cov-2_primers.fasta", "fasta"):
    r.letter_annotations["solexa_quality"] = [40] * len(r)
    print(r.format("fastq"), end='')
