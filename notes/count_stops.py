import pandas as pd
from Bio import SeqIO
import glob
import numpy as np


# look for all mutations to stop codons
# in the contingency tables we are only interested in ORF1ab and S genes
# we are interested in their position and whether or not they are pcr errors
# the final output (per sample) is a bunch of csv files:
#     stops.csv (position, from_base, to_base, frequency, depth)
#     stops_table.csv (contingency table)

reference = SeqIO.read("../data/MN908947.3.fasta", format="fasta")

orfs = pd.read_csv("../data/sars-cov-2-orfs.tsv", delim_whitespace=True, header=None, names=["name", "start", "end"])
orfs["start_mod_3"] = orfs["start"] % 3
orfs["end_mod_3"] = orfs["end"] % 3
orfs["length"] = (orfs["end"] - orfs["start"] + 1) // 3
orfs["start_seq"] = orfs.apply(lambda x: str(reference.seq[x.start-1:x.start+19]), axis=1)
orfs["end_seq"] = orfs.apply(lambda x: str(reference.seq[x.end-20:x.end]), axis=1)

orf_starts = [x.start for x in orfs.itertuples() if x.name != "N"]
assert(orf_starts == sorted(orf_starts))
orf_ends = [x.end for x in orfs.itertuples() if x.name != "N"]
assert(orf_ends == sorted(orf_ends))
orf_frames = [x.start_mod_3 for x in orfs.itertuples() if x.name != "N"]

def align_orfs(new_reference, orf_dataframe):
    """
    Use the orf_dataframe which is the dataframe for orf starts/end positions for the wuhan reference (above).
    Give it a new reference genome. 
    Get out the corresponding dataframe for the new reference. 
    """
    new_dataframe = pd.DataFrame(orf_dataframe["name"])
    starts, ends, start_seqs, end_seqs = [], [], [], []
    
    l = len(new_reference.seq)
    
    for r in orf_dataframe.itertuples():
        left_end = max(0, r.start-500)
        right_end = min(r.start+500, l)
        alignments = pairwise2.align.localxs(new_reference.seq[left_end:right_end], r.start_seq, -6., -1.)
        a = alignments[0]
        s = a.start
        start = left_end + s + 1
        start_seq = a.seqA[s:s+20]
        
        left_end = max(0, r.end-500)
        right_end = min(r.end+500, l)
        alignments = pairwise2.align.localxs(new_reference.seq[left_end:right_end], r.end_seq, -6., -1.)
        a = alignments[0]
        s = a.start
        end = left_end + s + 20
        end_seq = a.seqA[s:s+20]
        
        starts.append(start)
        ends.append(end)
        start_seqs.append(start_seq)
        end_seqs.append(end_seq)
        
    new_dataframe["start"] = starts
    new_dataframe["end"] = ends
    new_dataframe["start_mod_3"] = new_dataframe["start"] % 3
    new_dataframe["end_mod_3"] = new_dataframe["end"] % 3
    new_dataframe["length"] = (new_dataframe["end"] - new_dataframe["start"] + 1) // 3
    new_dataframe["start_seq"] = start_seqs
    new_dataframe["end_seq"] = end_seqs
    
    return new_dataframe


def frame_lookup(position, orf_starts, orf_ends, orf_frames):
    """
    Use the lists of orf starts, ends, and frames to look up which frame any particular position is in. 
    """
    
    if position < orf_starts[0]:
        return 1
    
    for i, start in enumerate(orf_starts):
        if start > position:
            break
    
    if position <= orf_ends[i-1]:
        return -1 * (int(position - orf_frames[i-1]) % 3)
    
    else:
        return 1

def is_stop_codon(position, mutation, reference, orf_starts, orf_ends, orf_frames):
    
    frame = frame_lookup(position, orf_starts, orf_ends, orf_frames)
    
    # if the position is non-coding then it's not a stop codon
    if frame == 1:
        return False
    
    codon = list(reference.seq[position+frame-1:position+frame+2])
    codon[-1 * frame] = mutation
    codon = "".join(codon)
    return (codon in ["TAG", "TAA", "TGA"])


