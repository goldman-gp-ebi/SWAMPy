# All this file does is take a bed file of primers (tab separated with columns 'genome', 'start', 'end', 'name', 'pool', 'sense')
# and print out a fastq of these primer with dummy quality scores. 
# It's not used in the main code - but it can be handy if you want to use a different set of primers from ARTIC V3. 
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SeqRecord


if __name__ == "__main__":

    genome_file = "../data/MN908947.3.fasta"
    bed_file = "../data/protocols/artic_v4_SARS-CoV-2.primer.bed"

    amplicons = pd.read_csv(bed_file, sep="\t", names=['genome', 'start', 'end', 'name'], usecols=[i for i in range(4)])

    genome = SeqIO.read(genome_file, format="fasta")

    for _, r in amplicons.iterrows():
        start = int(r[1])
        end = int(r[2])
        primer_id = r[3]
        if "_RIGHT" in primer_id:
            primer = SeqRecord(genome.seq[start:end].reverse_complement())
        else:
            primer = SeqRecord(genome.seq[start:end])
        primer.description = ""
        primer.id = primer_id
        primer.letter_annotations["solexa_quality"] = [40] * len(primer)
        print(primer.format("fastq"), end='')