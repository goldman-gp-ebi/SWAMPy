import os
import numpy as np
import random
from Bio import SeqIO,pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
import pandas as pd
import glob
import sys
from collections import Counter


os.chdir("/nfs/research/goldman/rabia/wastewater/sars-cov-2-metagenomic-simulator/src")

REF=SeqIO.read("../../data/ref/MN908947.3.fasta", format="fasta")

SUBS_RATE=0.01
INS_RATE=0.0001
DEL_RATE=0.0002

DEL_LENGTH_GEOMETRIC_PARAMETER=0.3
INS_MAX_LENGTH=10

SUBS_VAF_DIRICLET_PARAMETER=[0.8,0.2]
INS_VAF_DIRICLET_PARAMETER=[0.8,0.2]
DEL_VAF_DIRICLET_PARAMETER=[0.8,0.2]
VAF=dict(SUBS=SUBS_VAF_DIRICLET_PARAMETER,INS=INS_VAF_DIRICLET_PARAMETER,DEL=DEL_VAF_DIRICLET_PARAMETER)

SUBS_COUNT=int(SUBS_RATE*len(REF.seq))
INS_COUNT=int(INS_RATE*len(REF.seq))
DEL_COUNT=int(DEL_RATE*len(REF.seq))

PRIMER_BED="../example/articV3_no_alt.bed"

#we are picking position in genome without replacement


#Produce an alternative allele for a given error.
def alts(ref,type,len=0):
    if type=="SUBS":
        nucs=["A","C","G","T"]
        nucs.remove(ref)
        substitute=random.choice(nucs)
        return substitute
    
    elif type=="DEL":
        return "-"

    else: #insertion
        insert=random.choices(["A","C","G","T"],k=len)
        return ref+ "".join(insert)


def amplicon_lookup(PRIMER_BED,position):

    df=pd.read_csv(PRIMER_BED,sep="\t", names=["Chr","Start","End","Amplicon","Pool","Strand"])
    df["Handedness"]=df.apply(lambda x: x.Amplicon.split("_")[-1],axis=1)
    df["Amp_no"]=df.apply(lambda x: x.Amplicon.split("_")[-2],axis=1)
    df.drop(["Chr","Amplicon","Pool","Strand"],inplace=True,axis=1)

    df = pd.merge(
          df.loc[df["Handedness"] == "LEFT"], 
          df.loc[df["Handedness"] == "RIGHT"], 
          on=["Amp_no"]
      )
    corresponding_amplicons=list(df.loc[(df["Start_x"]<=position)&(position<=df["End_y"]),"Amp_no"])
    return corresponding_amplicons
 
def choose_reads_to_mutate(df_row):
    reads=[]
    for amp in df_row['amplicons']:
        amp # type: str
        to_mut=df_row.amplicon_source_genome_read_count_to_mutate[amp]
        for genome,n_mut in to_mut.items():
            if n_mut == 0 or df_row.amplicon_source_genome_read_count[amp][genome]==0:
                pass
            else:
                indices=random.sample(list(range(df_row.amplicon_source_genome_read_count[amp][genome]+1)),
                                    k=n_mut)

                indices=[(x+1)*2 for x in indices] #read names are 1 based and even numbers
                for idx in indices:
                    reads.append(f'{genome}_amplicon_{amp}-{idx}/1')
                    reads.append(f'{genome}_amplicon_{amp}-{idx}/2')
    return reads

amp_n_reads=pd.read_csv("../simulation_output/example_amplicon_abundances_summary.tsv",sep="\t")
del amp_n_reads['Unnamed: 0']
amp_n_reads['amplicon_number']=amp_n_reads['amplicon_number'].apply(str)
amp_n_reads_grouped=amp_n_reads.groupby(['amplicon_number']).sum()
amp_n_reads_grouped=amp_n_reads_grouped['n_reads']

errors=pd.DataFrame(dict(errortype=["SUBS"]*SUBS_COUNT + ["INS"]*INS_COUNT + ["DEL"]*DEL_COUNT))
errors['mut_index']=errors.index
errors["length"]=[1]*SUBS_COUNT + random.choices(list(range(1,INS_MAX_LENGTH+1)),k=INS_COUNT) + list(np.random.geometric(p=DEL_LENGTH_GEOMETRIC_PARAMETER, size=DEL_COUNT))
errors["pos"]= random.sample(list(range(len(REF.seq))),k=SUBS_COUNT+INS_COUNT+DEL_COUNT)
errors["ref"]=errors.apply(lambda x: REF.seq[x.pos], axis=1)
errors["alt"]=errors.apply(lambda x: alts(x.ref,x.errortype,x.length), axis=1)
errors["VAF"]=errors.apply(lambda x: np.random.dirichlet(VAF[x.errortype], size=None)[1],axis=1)
errors["amplicons"]=errors.apply(lambda x: amplicon_lookup(PRIMER_BED,x.pos),axis=1)
errors["amplicon_total_reads"]=errors.apply(lambda x:{amp:int(amp_n_reads_grouped[amp]) if amp in amp_n_reads_grouped.index else 0 for amp in x.amplicons} ,axis=1)
errors["read_count_to_mutate_float"]=errors.apply(lambda x: {amp:int(amp_n_reads_grouped[amp]) * x.VAF if amp in amp_n_reads_grouped.index else 0 for amp in x.amplicons},axis=1)
errors["read_count_to_mutate"]=errors.apply(lambda x: {amp:round(int(amp_n_reads_grouped[amp]) * x.VAF) if amp in amp_n_reads_grouped.index else 0 for amp in x.amplicons},axis=1)
errors = errors.loc[errors['VAF']!=0,]
errors["amplicon_source_genome_read_count"]=errors.apply(lambda x:{amp_of_genome:dict(amp_n_reads.loc[amp_n_reads['amplicon_number']==amp_of_genome,["ref","n_reads"]].to_records(index=False)) if amp_of_genome in list(amp_n_reads['amplicon_number']) else {} for amp_of_genome in x.amplicons } ,axis=1)
errors["amplicon_source_genome_read_count_to_mutate"]=errors.apply(lambda x:{amp_of_genome:Counter(random.sample(list(amp_n_reads.loc[amp_n_reads['amplicon_number']==amp_of_genome,"ref"]),counts=list(amp_n_reads.loc[amp_n_reads['amplicon_number']==amp_of_genome,"n_reads"]),k=x.read_count_to_mutate[amp_of_genome])) if amp_of_genome in list(amp_n_reads['amplicon_number']) else {} for amp_of_genome in x.amplicons } ,axis=1)
errors["amplicon_read_to_mutate"]=errors.apply(lambda x:choose_reads_to_mutate(x),axis=1)


reads=[]
for i in errors['amplicon_read_to_mutate'].tolist():
    reads.extend(i)


indices=[[errors["mut_index"][a]]*2*sum(errors.loc[a,"read_count_to_mutate"].values()) for a in range(errors.shape[0])]
indices2=[]
for i in indices:
    indices2.extend(i)

read_df=pd.DataFrame(zip(reads,indices2),columns=["read","mut_index"])
read_df2=read_df.groupby("read").sum()

r1_list=list()
for all_r1 in SeqIO.parse("tmp.sms.all_files_unshuffled1.fastq", format="fastq"):
    if all_r1.id in reads:
        seq_list=list(all_r1.seq)
        qual_list=all_r1.letter_annotations['phred_quality']
        muts=read_df.loc[read_df["read"]==all_r1.id,"mut_index"].tolist()
        align=pairwise2.align.localms(REF.seq, all_r1.seq,5, -1, -3, -0.5)
        print(format_alignment(*align[0],full_sequences=True))

    else:
        pass#write back to the new file

all_r1.id



for idx,err in enumerate(errors.amplicons):
    for j in err:
        print(glob.glob(f"*amplicon_{j}.1.fq"))
    break












lens=[]
for genome in SeqIO.parse("tmp.sms.hCoV-19&Scotland&QEUH-15AFBCE&2021|EPI_ISL_2354066|2021-05-22_amplicon_40.1.fq", format="fastq"):
    gen1=genome
    break


seq1=genome.seq
seq2
set(lens)

seq_list=list(genome.seq)
seq_list[0:10]
seq_list[0]="A"
genome.seq=Seq("".join((seq_list)))
seq3=list(seq1)
seq3=seq3[:50]+seq3[53:]
seq3=Seq("".join(seq3))
seq4=list(seq1)
seq4=seq4[:50]+["A","C","G"]+seq4[50:]
seq4=Seq("".join(seq4))

ref=REF.seq
alignments = pairwise2.align.localms(seq1, seq2.reverse_complement(),10, -10, -3, -0.5)
print(format_alignment(*alignments[0],full_sequences=True))
alignments[0:5]



align1=pairwise2.align.localms(ref,seq1,5, -2, -3, -0.5)
align2=pairwise2.align.localms(ref, seq2.reverse_complement(),10, -10, -10, -0.5)
align3=pairwise2.align.localms(ref, seq3,5, -2, -3, -0.5)
align4=pairwise2.align.localms(ref, seq4,5, -1, -3, -0.5)
print(format_alignment(*align1[0]))
print(format_alignment(*align2[0]))
print(format_alignment(*align3[0]))
print(format_alignment(*align4[0]))
align1[0].start
align2[0].start
ref[12060:12065]
seq1[0:5]


align3[0].seqB[align3[0].start:align3[0].end]
align4[0].seqB[align4[0].start:align4[0].end]
align4[0].seqA[align4[0].start:align4[0].end]

os.getcwd()