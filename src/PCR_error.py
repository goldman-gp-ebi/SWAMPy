import os
import numpy as np
import random
from Bio import SeqIO
import pandas as pd
import subprocess
from io import StringIO
import re


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

#find amplicons corresponding to a given position
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

#create a dataframe of errors that we want to introduce
errors=pd.DataFrame(dict(errortype=["SUBS"]*SUBS_COUNT + ["DEL"]*DEL_COUNT + ["INS"]*INS_COUNT ))
errors["mut_indices"]=errors.index
errors["length"]=[1]*SUBS_COUNT + list(np.random.geometric(p=DEL_LENGTH_GEOMETRIC_PARAMETER, size=DEL_COUNT)) + random.choices(list(range(1,INS_MAX_LENGTH+1)),k=INS_COUNT) 
errors["pos"]= random.sample(list(range(len(REF.seq))),k=SUBS_COUNT+INS_COUNT+DEL_COUNT)
errors["ref"]=errors.apply(lambda x: REF.seq[x.pos], axis=1)
errors["alt"]=errors.apply(lambda x: alts(x.ref,x.errortype,x.length), axis=1)
errors["VAF"]=errors.apply(lambda x: np.random.dirichlet(VAF[x.errortype], size=None)[1],axis=1)
errors["amplicons"]=errors.apply(lambda x: amplicon_lookup(PRIMER_BED,x.pos),axis=1)
errors = errors.loc[errors['VAF']!=0,]

#all amplicons to be mutated
error_amplicons=[]
for i in errors['amplicons'].tolist():
    error_amplicons.extend(i)

#correasponding error index for that amplicon
indices=[[errors.index[a]]*len(errors.loc[a,"amplicons"]) for a in range(errors.shape[0])]
indices2=[]
for i in indices:
    indices2.extend(i)


amplicons=[]
n_reads=[]

for i in df_amplicons.itertuples():
    
    #Sometimes more than 1 error are introduced to the same amplicon.
    #Find all errors that will be introduced to that amplicon
    mut_indices=[indices2[idx] for idx,a in enumerate(error_amplicons) if a==str(i.amplicon_number)]
    
    if mut_indices==[]:

        amplicons.append(i.amplicon_filepath)
        n_reads.append(i.n_reads)
        
    else:
        # | and & are problematic for bash (subprocess), escape those.
        amp=i.amplicon_filepath.replace("&","\\&").replace("|","\\|")

        #align the original amplicon to the reference becuase there could be real indels in the source genome.
        alignment = subprocess.run(
            ["bowtie2", 
            "-x", "../../data/ref/MN908947.3", 
            "-f",f"../example/amplicons/{amp}"], capture_output=True)

        alignment = StringIO(alignment.stdout.decode("UTF-8"))

        # read alignment data as a dataframe
        df = pd.read_csv(alignment, sep="\t", skiprows=[0, 1, 2], header=None, names=[i for i in range(19)])
        df = pd.DataFrame(df[[0,3,5,9]])
        df = df.rename(columns={0:"name", 3:"start", 5:"CIGAR", 9: "seq"})

        #if the amplicon contains too many Ns it will not align, skip introducing PCR error to those
        if df.CIGAR[0]=="*":
            amplicons.append(i.amplicon_filepath)
            n_reads.append(i.n_reads)
        
        else:

            #Record the CIGAR as a long string. i.e. "MMMII" instead of "3M2I"
            CIGAR=""
            for idx,cigar in enumerate(re.split("(M|D|I)", df.CIGAR[0])[:-1]):
                if idx%2==0:
                    prev=cigar[:]
                else:
                    CIGAR+= cigar*int(prev)

            #Start position and the sequence of the alignment(amplicon)
            start_p=df.start[0]-1
            seq=df.seq[0]

            #Create an empty dataframe to hold how many reads each error combination will produce at the end.
            reads_df=pd.DataFrame()

            #For each error that will be introduced to this amplicon, 
            #find the indices for slicing wrt to amplicon seq left end.
            seq_pos=[]
            for mut_idx in mut_indices:

                #we aim for the error's position wrt reference genome
                aim=errors.loc[mut_idx,"pos"]

                #amplicon slicing index starts from 0 (left end)
                seq_idx=0

                #reference index starts from the position where amplicon alignment starts.
                ref_idx=start_p

                for c_idx,c in enumerate(CIGAR):
                    if ref_idx==aim:
                        seq_pos.append(seq_idx)
                        break
                    if c=="M":
                        seq_idx += 1
                        ref_idx +=1
                    elif c=="D":
                        ref_idx+=1
                    elif c=="I":
                        seq_idx+=1
                    if c_idx==len(CIGAR)-1 and ref_idx==aim:#final letter
                        seq_pos.append(seq_idx)


                #How many reads this specific error will have
                mut_reads=int(i.n_reads * errors.loc[mut_idx,"VAF"])

                #if number of reads and/or VAF are small, this can be 0
                if mut_reads == 0:
                    pass
                
                else:
                    #Take that many samples from the total of the imaginary reads of that amplicon
                    reads=sorted(random.sample(range(i.n_reads),k=mut_reads))
                    reads=[str(a) for a in reads]
                    reads=[a+"," for a in reads] #, will be used for grouping

                    muts=[str(mut_idx)]*mut_reads #keep track of the mutation 
                    muts=[a+"," for a in muts]

                    read_df=pd.DataFrame(dict(reads=reads,muts=muts))
                    reads_df=reads_df.append(read_df,ignore_index=True)

            #if there are no errors with a non-zero count, skip introducing errors to that amplicon
            if reads_df.empty:
                amplicons.append(i.amplicon_filepath)
                n_reads.append(i.n_reads)

            else:
                #group wrt imaginary reads to see which ones ended up with wich errors
                reads_df=reads_df.groupby("reads").sum()
                reads_df["count"]=1 
                #group by different combinations of errors to count how many read each combination will produce.
                reads_df=reads_df.groupby("muts",as_index=False).sum()
                #remove the , and turn tham into a list
                reads_df["muts"]=reads_df.apply(lambda x: x.muts.split(",")[:-1] ,axis=1)

                #create a dataframe of all errors of the amplicon. Contains pos, mut_index, errortype, length, alt
                seq_pos_df=pd.DataFrame(dict(seq_pos=seq_pos,mut_indices=mut_indices))
                seq_pos_df=seq_pos_df.merge(errors[["errortype","mut_indices","length","alt"]],on="mut_indices")

                # amplicon's number of reads - total count of all error combination versions is the count of non-mutated (old) version.
                amplicons.append(i.amplicon_filepath)
                n_reads.append(i.n_reads - sum(reads_df["count"]))

                #for all error combination versions of the amplicon
                for idx,pcr_error in enumerate(reads_df.itertuples()):
                    #if a specific combination has 0 reads, pass
                    if pcr_error.count==0:
                        pass
                    else:
                        #create a final df that contains all the errors in that specific combination
                        final_df=seq_pos_df.loc[seq_pos_df["mut_indices"].isin([int(a) for a in pcr_error.muts]),]
                        final_df=final_df.sort_values('seq_pos')
                        
                        #introduce those errors one by one.
                        new_seq=""
                        for indx,final in enumerate(final_df.itertuples()):
                            
                            #the part up to the first error is the same
                            if indx==0:
                                new_seq=new_seq + seq[0:final.seq_pos]

                            #if an error is substition or indel, take the part up to and excluding the error position
                            #add alternative instead of the ref at error pos.
                            if final.errortype=="SUBS" or final.errortype=="INS":
                                new_seq=new_seq+final.alt
                                #then add the part up to the next error
                                try:
                                    new_seq=new_seq+ seq[final.seq_pos+1:final_df.loc[indx+1,"seq_pos"]]
                                #if it is the last error, add all the remaining sequence
                                except KeyError: 
                                    new_seq=new_seq+ seq[final.seq_pos+1:]

                            elif final.errortype=="DEL":
                                #if it is a deletion add the next section of the sequence but leave out the first n bases of it
                                try:
                                    new_seq=new_seq+ seq[final.seq_pos+1:final_df.loc[indx+1,"seq_pos"]][final.length-1:]
                                except KeyError:
                                    new_seq=new_seq+ seq[final.seq_pos+1:][final.length-1:]

                        #add the new amplicon to the list
                        new_path=i.amplicon_filepath[:-6] + "_p" +str(idx+1) + ".fasta"
                        amplicons.append(new_path)
                        n_reads.append(pcr_error.count)
    
                        #write the fasta file of the new amplicon. 
                        #Name all the PCR error combinations as _p1, _p2 and etc.
                        with open(f"../example/amplicons/{new_path}","w") as new_a:
                            new_a.write(f">{new_path[:-6]}\n")
                            new_a.write(new_seq + "\n\n")
    