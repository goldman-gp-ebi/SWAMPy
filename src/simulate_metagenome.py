import argparse
from genericpath import exists
from os.path import dirname, join, abspath, basename
import os
import glob
import logging
from Bio import SeqIO
import subprocess
import pandas as pd
import numpy as np
import random
import shutil

from art_runner import art_illumina
from create_amplicons import build_index, align_primers, write_amplicon
from read_model import get_amplicon_reads_sampler
from PCR_error import add_PCR_errors


# All these caps variables are set once (by user inputs, with default values) but then never touched again.
BASE_DIR = join(dirname(dirname(abspath(__file__))), "example")
TEMP_FOLDER= join(BASE_DIR, "temp")
GENOMES_FILE = join(BASE_DIR, "genomes.fasta")
ABUNDANCES_FILE = join(BASE_DIR, "abundances.tsv")
PRIMER_SET="a1"
PRIMER_SET_FOLDER=join(dirname(dirname(abspath(__file__))), "primer_sets")
OUTPUT_FOLDER = os.getcwd()
OUTPUT_FILENAME_PREFIX = "example"
N_READS = 100000
READ_LENGTH = 250
SEQ_SYS = "MSv3"
SEED = np.random.randint(1000000000)
AMPLICON_DISTRIBUTION = "DIRICHLET_1"
AMPLICON_PSEUDOCOUNTS = 200


##PCR-error related variables:
WUHAN_REF = join(dirname(dirname(abspath(__file__))), "ref","MN908947.3")
U_SUBS_RATE = 0.002485
U_INS_RATE = 0.00002
U_DEL_RATE = 0.000115
R_SUBS_RATE = 0.003357
R_INS_RATE = 0.00002
R_DEL_RATE = 0

DEL_LENGTH_GEOMETRIC_PARAMETER = 0.69
INS_MAX_LENGTH = 14

SUBS_VAF_DIRICLET_PARAMETER = "0.29,1.89"
INS_VAF_DIRICLET_PARAMETER = "0.33,0.45"
DEL_VAF_DIRICLET_PARAMETER = "0.59,0.41"

R_SUBS_VAF_DIRICLET_PARAMETER = SUBS_VAF_DIRICLET_PARAMETER
R_INS_VAF_DIRICLET_PARAMETER = INS_VAF_DIRICLET_PARAMETER
R_DEL_VAF_DIRICLET_PARAMETER = DEL_VAF_DIRICLET_PARAMETER

def setup_parser():
    parser = argparse.ArgumentParser(description="Run SARS-CoV-2 metagenome simulation.")
    parser.add_argument("--genomes_file", metavar='', help="File containing all of the genomes that might be used", default=GENOMES_FILE)
    parser.add_argument("--temp_folder", "-t", metavar='', help="A path for a temporary output folder to store intemediate files. Including FASTA files of genomes, amplicons, and their bowtie2 indices", default=TEMP_FOLDER)
    parser.add_argument("--genome_abundances", "-ab", metavar='', help="TSV of genome abundances.", default=ABUNDANCES_FILE)
    parser.add_argument("--primer_set", "-ps", metavar='', help="Primer set can be either a1 for Artic v1, a4 for Artic v4 and n2 for Nimagen v2, Default is a1.", default="a1",choices=["a1","a4","n2"])
    parser.add_argument("--output_folder", "-o", metavar='', help="A path for a folder where the output fastq files will be stored. Default is working directory", default=OUTPUT_FOLDER)
    parser.add_argument("--output_filename_prefix", "-x", metavar='', help="Name of the fastq files name1.fastq, name2.fastq", default=OUTPUT_FILENAME_PREFIX)
    parser.add_argument("--seqSys", metavar='', help="Name of the sequencing system, options to use are given by the art_illumina help text, and are:" + 
    """GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
           HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
           HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
           MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)""", default="MSv3")
    parser.add_argument("--n_reads", "-n", metavar='', help="Approximate number of reads in fastq file (subject to sampling stochasticity).", default=N_READS)
    parser.add_argument("--read_length", "-l", metavar='', help="Length of reads taken from the sequencing machine.", default=READ_LENGTH)
    parser.add_argument("--seed", "-s", metavar='', help="Random seed", default=SEED)
    parser.add_argument("--quiet", "-q", help="Add this flag to supress verbose output." ,action='store_true')
    parser.add_argument("--amplicon_distribution",help= "Default is DIRICHLET1", metavar='', default=AMPLICON_DISTRIBUTION)
    parser.add_argument("--amplicon_pseudocounts","-c", metavar='', default=AMPLICON_PSEUDOCOUNTS)
    parser.add_argument("--autoremove", action='store_true',help="Delete temproray files after execution.")
    parser.add_argument("--no_pcr_errors", action='store_true',help="Turn off PCR errors. The output will contain only sequencing errors. Other PCR-error related options will be ignored")
    parser.add_argument("--unique_insertion_rate","-ins", metavar='', help="PCR insertion error rate. Unique to one source genome in the mixture Default is 0.00002", default=U_INS_RATE)
    parser.add_argument("--unique_deletion_rate","-del", metavar='', help="PCR deletion error rate. Unique to one source genome in the mixture Default is 0.000115", default=U_DEL_RATE)
    parser.add_argument("--unique_substitution_rate","-subs", metavar='', help="PCR substitution error rate. Unique to one source genome in the mixture Default is 0.002485", default=U_SUBS_RATE)
    parser.add_argument("--recurrent_insertion_rate","-rins", metavar='', help="PCR insertion error rate. Recurs across source genomes. Default is 0.00002", default=R_INS_RATE)
    parser.add_argument("--recurrent_deletion_rate","-rdel", metavar='', help="PCR deletion error rate. Recurs across source genomes. Default is 0", default=R_DEL_RATE)
    parser.add_argument("--recurrent_substitution_rate","-rsubs", metavar='', help="PCR substitution error rate. Recurs across source genomes. Default is 0.003357", default=R_SUBS_RATE)
    parser.add_argument("--deletion_length_p","-dl", metavar='', help="Geometric distribution parameter, p, for PCR deletion length. Default is 0.69", default=DEL_LENGTH_GEOMETRIC_PARAMETER)
    parser.add_argument("--max_insertion_length","-il", metavar='', help="Maximum PCR insertion length in bases (uniform distribution boundry). Default is 14", default=INS_MAX_LENGTH)
    parser.add_argument("--subs_VAF_alpha","-sv", metavar='', help="alpha1,alpha2 of the Dirichlet distribution for VAF of the unique PCR error. Default is 0.29,1.89", default=SUBS_VAF_DIRICLET_PARAMETER)
    parser.add_argument("--del_VAF_alpha","-dv", metavar='', help="alpha1,alpha2 of the Dirichlet distribution for VAF of the unique PCR error. Default is 0.59,0.41", default=DEL_VAF_DIRICLET_PARAMETER)
    parser.add_argument("--ins_VAF_alpha","-iv", metavar='', help="alpha1,alpha2 of the Dirichlet distribution for VAF of the unique PCR error. Default is 0.33,0.45", default=INS_VAF_DIRICLET_PARAMETER)
    parser.add_argument("--r_subs_VAF_alpha","-rsv", metavar='', help="alpha1,alpha2 of the Dirichlet distribution for VAF of the recurrent PCR error. Default is equal to unique erros", default=SUBS_VAF_DIRICLET_PARAMETER)
    parser.add_argument("--r_del_VAF_alpha","-rdv", metavar='', help="alpha1,alpha2 of the Dirichlet distribution for VAF of the recurrent PCR error. Default is equal to unique erros", default=DEL_VAF_DIRICLET_PARAMETER)
    parser.add_argument("--r_ins_VAF_alpha","-riv", metavar='', help="alpha1,alpha2 of the Dirichlet distribution for VAF of the recurrent PCR error. Default is equal to unique erros", default=INS_VAF_DIRICLET_PARAMETER)
    return parser

def load_command_line_args():
    parser = setup_parser()
    args = parser.parse_args()

    global TEMP_FOLDER
    TEMP_FOLDER = args.temp_folder
    if not os.path.exists(TEMP_FOLDER):
        os.makedirs(TEMP_FOLDER,exist_ok=True)

    global GENOMES_FOLDER
    GENOMES_FOLDER = join(TEMP_FOLDER, "genomes")
    if not os.path.exists(GENOMES_FOLDER):
        os.makedirs(GENOMES_FOLDER,exist_ok=True)

    global GENOMES_FILE
    GENOMES_FILE = args.genomes_file
    global GENOMES_FILE2
    GENOMES_FILE2 =join(GENOMES_FOLDER,basename(GENOMES_FILE))

    global AMPLICONS_FOLDER
    AMPLICONS_FOLDER = join(TEMP_FOLDER, "amplicons")
    if not os.path.exists(AMPLICONS_FOLDER):
        os.makedirs(AMPLICONS_FOLDER,exist_ok=True)
    
    global INDICES_FOLDER
    INDICES_FOLDER = join(TEMP_FOLDER, "indices")
    if not os.path.exists(INDICES_FOLDER):
        os.mkdir(INDICES_FOLDER)

    global ABUNDANCES_FILE 
    ABUNDANCES_FILE = args.genome_abundances
    global ABUNDANCES_FILE2
    ABUNDANCES_FILE2=join(GENOMES_FOLDER,basename(ABUNDANCES_FILE))

    global OUTPUT_FOLDER
    OUTPUT_FOLDER = args.output_folder
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER,exist_ok=True)

    global OUTPUT_FILENAME_PREFIX
    OUTPUT_FILENAME_PREFIX = args.output_filename_prefix
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(join(OUTPUT_FOLDER, f"{OUTPUT_FILENAME_PREFIX}.log")),
            logging.StreamHandler()
        ]
    )

    global PRIMER_SET
    PRIMER_SET=args.primer_set

    global PRIMERS_FILE
    if PRIMER_SET=="a1":
        PRIMERS_FILE = join(PRIMER_SET_FOLDER,"artic_v3_primers_no_alts.fastq")
        logging.info(f"Primer set: Artic v1")
    elif PRIMER_SET=="a4":
        PRIMERS_FILE = join(PRIMER_SET_FOLDER,"artic_v4_primers.fastq")
        logging.info(f"Primer set: Artic v4")
    elif PRIMER_SET=="n2":
        PRIMERS_FILE = join(PRIMER_SET_FOLDER,"nimagen_v2_primers.fastq")
        logging.info(f"Primer set: Nimagen v2")

    global N_READS
    N_READS = int(args.n_reads)
    logging.info(f"Number of reads: {N_READS}")

    global READ_LENGTH
    READ_LENGTH = int(args.read_length)

    global SEED
    SEED = args.seed
    np.random.seed(int(SEED))
    random.seed(int(SEED))
    logging.info(f"Random seed: {SEED}")

    global VERBOSE
    VERBOSE = not args.quiet

    global AMPLICON_DISTRIBUTION
    AMPLICON_DISTRIBUTION = args.amplicon_distribution

    global AMPLICON_DISTRIBUTION_FILE
    if PRIMER_SET=="a1":
        AMPLICON_DISTRIBUTION_FILE = join(PRIMER_SET_FOLDER, "artic_v3_amplicon_distribution.tsv")
    elif PRIMER_SET=="a4":
        AMPLICON_DISTRIBUTION_FILE = join(PRIMER_SET_FOLDER, "artic_v4_amplicon_distribution.tsv")
    elif PRIMER_SET=="n2":
        AMPLICON_DISTRIBUTION_FILE = join(PRIMER_SET_FOLDER, "nimagen_v2_amplicon_distribution.tsv")

    global AMPLICON_PSEUDOCOUNTS
    AMPLICON_PSEUDOCOUNTS = int(args.amplicon_pseudocounts)
    logging.info(f"Amplicon pseudocounts/ i.e. quality parameter: {AMPLICON_PSEUDOCOUNTS}")

    global AUTOREMOVE
    AUTOREMOVE = args.autoremove

    ##PCR error arguments

    global NO_PCR_ERRORS
    NO_PCR_ERRORS = args.no_pcr_errors

    global PRIMER_BED 
    if PRIMER_SET=="a1":
        PRIMER_BED = join(PRIMER_SET_FOLDER,"articV3_no_alt.bed")
    elif PRIMER_SET=="a4":
        PRIMER_BED = join(PRIMER_SET_FOLDER,"articV4.bed")
    elif PRIMER_SET=="n2":
        PRIMER_BED = join(PRIMER_SET_FOLDER,"nimagenV2.bed")

    global U_SUBS_RATE
    U_SUBS_RATE = float(args.unique_substitution_rate)

    global U_INS_RATE
    U_INS_RATE = float(args.unique_insertion_rate)

    global U_DEL_RATE
    U_DEL_RATE = float(args.unique_deletion_rate)

    global R_SUBS_RATE
    R_SUBS_RATE = float(args.recurrent_substitution_rate)

    global R_INS_RATE
    R_INS_RATE = float(args.recurrent_insertion_rate)

    global R_DEL_RATE
    R_DEL_RATE = float(args.recurrent_deletion_rate)

    global DEL_LENGTH_GEOMETRIC_PARAMETER
    DEL_LENGTH_GEOMETRIC_PARAMETER = float(args.deletion_length_p)

    global INS_MAX_LENGTH
    INS_MAX_LENGTH = int(args.max_insertion_length)

    global SUBS_VAF_DIRICLET_PARAMETER
    SUBS_VAF_DIRICLET_PARAMETER = args.subs_VAF_alpha.split(",")
    if len(SUBS_VAF_DIRICLET_PARAMETER)!=2:
        logging.error(f"subs_VAF_alpha argument must be a list of 2 values seperated by comma. Example: 0.5,0.4. You entered {args.subs_VAF_alpha} ")
        exit(1)
    else:
        SUBS_VAF_DIRICLET_PARAMETER=[float(a) for a in SUBS_VAF_DIRICLET_PARAMETER]
    
    global INS_VAF_DIRICLET_PARAMETER
    INS_VAF_DIRICLET_PARAMETER = args.ins_VAF_alpha.split(",")
    if len(INS_VAF_DIRICLET_PARAMETER)!=2:
        logging.error(f"ins_VAF_alpha argument must be a list of 2 values seperated by comma. Example: 0.5,0.4. You entered {args.ins_VAF_alpha} ")
        exit(1)
    else:
        INS_VAF_DIRICLET_PARAMETER=[float(a) for a in INS_VAF_DIRICLET_PARAMETER]

    global DEL_VAF_DIRICLET_PARAMETER
    DEL_VAF_DIRICLET_PARAMETER = args.del_VAF_alpha.split(",")
    if len(DEL_VAF_DIRICLET_PARAMETER)!=2:
        logging.error(f"del_VAF_alpha argument must be a list of 2 values seperated by comma. Example: 0.5,0.4. You entered {args.del_VAF_alpha} ")
        exit(1)
    else:
        DEL_VAF_DIRICLET_PARAMETER=[float(a) for a in DEL_VAF_DIRICLET_PARAMETER]

    global R_SUBS_VAF_DIRICLET_PARAMETER
    R_SUBS_VAF_DIRICLET_PARAMETER = args.r_subs_VAF_alpha.split(",")
    if len(R_SUBS_VAF_DIRICLET_PARAMETER)!=2:
        logging.error(f"r_subs_VAF_alpha argument must be a list of 2 values seperated by comma. Example: 0.5,0.4. You entered {args.r_subs_VAF_alpha} ")
        exit(1)
    else:
        R_SUBS_VAF_DIRICLET_PARAMETER=[float(a) for a in R_SUBS_VAF_DIRICLET_PARAMETER]
    
    global R_INS_VAF_DIRICLET_PARAMETER
    R_INS_VAF_DIRICLET_PARAMETER = args.r_ins_VAF_alpha.split(",")
    if len(R_INS_VAF_DIRICLET_PARAMETER)!=2:
        logging.error(f"r_ins_VAF_alpha argument must be a list of 2 values seperated by comma. Example: 0.5,0.4. You entered {args.r_ins_VAF_alpha} ")
        exit(1)
    else:
        R_INS_VAF_DIRICLET_PARAMETER=[float(a) for a in R_INS_VAF_DIRICLET_PARAMETER]

    global R_DEL_VAF_DIRICLET_PARAMETER
    R_DEL_VAF_DIRICLET_PARAMETER = args.r_del_VAF_alpha.split(",")
    if len(R_DEL_VAF_DIRICLET_PARAMETER)!=2:
        logging.error(f"r_del_VAF_alpha argument must be a list of 2 values seperated by comma. Example: 0.5,0.4. You entered {args.r_del_VAF_alpha} ")
        exit(1)
    else:
        R_DEL_VAF_DIRICLET_PARAMETER=[float(a) for a in R_DEL_VAF_DIRICLET_PARAMETER]


if __name__ == "__main__":

    # STEP 0: Read command line arguments
    load_command_line_args()


    # Change spaces with "_" in genomes fasta file and record as a different file.
    os.system(f"sed 's/ /_/g' {GENOMES_FILE} > {GENOMES_FILE2}")
    os.system(f"sed 's/ /_/g' {ABUNDANCES_FILE} > {ABUNDANCES_FILE2}")
    if VERBOSE:
        logging.info("Spaces in the genomes and abundances files are processed as '_' characters if exist")
    
    # STEP 1: Simulate Viral Population

    # Read genome abundances csv file
    genome_abundances = {}
    df_amplicons = pd.DataFrame()

    with open(ABUNDANCES_FILE2) as ab_file:
        for line in ab_file:
            name, relative_abundance = tuple(line.split("\t"))
            genome_abundances[name] = float(relative_abundance)

    if abs(sum(genome_abundances.values()) - 1) > 0.000000001:
        total = sum(genome_abundances.values())
        if total <= 0:
            logging.info(f"The total genome abundance is set to {total}, which is impossible.")
            exit(1)

        logging.info(f"Total of relative abundance values is {total}, not 1.")
        logging.info("Continuing, normalising total of genome abundances to 1.")

        for k in genome_abundances.keys():
            genome_abundances[k] /= total

    n_genomes = len(genome_abundances)

    # Split genome file into multiple separate files
    for genome in SeqIO.parse(GENOMES_FILE2, format="fasta"):
        filepath = genome.description.replace(" ", "&").replace("/", "&").replace(",", "&")
        filepath += ".fasta"
        SeqIO.write(genome, join(GENOMES_FOLDER, filepath), format="fasta")


    # STEP 2: Simulate Amplicon Population
    genome_counter = 0
    for genome_path in genome_abundances:
        genome_counter += 1
        genome_path = genome_path.replace(" ", "_").replace("/", "&").replace(",", "&") + ".fasta"
        genome_path = join(GENOMES_FOLDER, genome_path)
        genome_filename_short = ".".join(basename(genome_path).split(".")[:-1])
        reference = SeqIO.read(genome_path, format="fasta")

        # use bowtie2 to create a dataframe with positions of each primer pair aligned to the genome
        if VERBOSE:
            logging.info(f"Working on genome {genome_counter} of {n_genomes}")
            logging.info(f"Using bowtie2 to align primers to genome {reference.description}")

        build_index(genome_path, genome_filename_short, INDICES_FOLDER)
        df = align_primers(genome_path, genome_filename_short, INDICES_FOLDER, PRIMERS_FILE, False)        
        df["abundance"] = genome_abundances[df["ref"][0]]

        # write the amplicon to a file
        write_amplicon(df, reference, genome_filename_short, AMPLICONS_FOLDER)


        df_amplicons = pd.concat([df_amplicons, df])


    # pick total numbers of reads for each amplicon
    genome_count_sampler, amplicon_hyperparameter_sampler, amplicon_probability_sampler, amplicon_reads_sampler = get_amplicon_reads_sampler(
                                AMPLICON_DISTRIBUTION, 
                                AMPLICON_DISTRIBUTION_FILE, 
                                AMPLICON_PSEUDOCOUNTS, 
                                genome_abundances,
                                N_READS)

    df_amplicons["total_n_reads"] = N_READS

    # for each amplicon, look up what the dirichlet hyperparameter should be (parameter \alpha)
    df_amplicons["hyperparameter"] = df_amplicons.apply(amplicon_hyperparameter_sampler, axis=1)

    # for each genome, sample a total number of reads that should be shared between all of its amplicons
    # N_genome = Multinomial(N_reads, p_genomes)
    df_amplicons["genome_n_reads"] = df_amplicons.apply(genome_count_sampler, axis=1)

    # sample a p_amplicon vector from the dirichlet distribution - p_amplicon = Dir(\alpha)
    df_amplicons["amplicon_prob"] = df_amplicons.apply(amplicon_probability_sampler, axis=1)

    # sample a number of reads for the amplicons of each genome: Multinomial(N_genome, p_amplicon)
    df_amplicons["n_reads"] = df_amplicons.apply(amplicon_reads_sampler, axis=1)

    # write a summary csv
    df_amplicons[
        ["ref", 
        "amplicon_number", 
        "is_alt", 
        "total_n_reads", 
        "abundance",
        "genome_n_reads",
        "hyperparameter",
        "amplicon_prob",
        "n_reads"]].to_csv(join(OUTPUT_FOLDER, f"{OUTPUT_FILENAME_PREFIX}_amplicon_abundances_summary.tsv"), sep="\t")

    df_amplicons.reset_index(drop=True,inplace=True)

    if VERBOSE:
        logging.info(f"Total number of reads was {sum(df_amplicons['n_reads'])}, when {N_READS} was expected.")


    # STEP 3: Library Prep - PCR Amplification of Amplicons

    if NO_PCR_ERRORS:    
        amplicons = list(df_amplicons["amplicon_filepath"])
        n_reads = list(df_amplicons["n_reads"])
    else:
        if VERBOSE:
            logging.info(f"Introducing PCR errors")

        amplicons,n_reads,vcf_errordf=add_PCR_errors(df_amplicons,genome_abundances,PRIMER_BED,WUHAN_REF,AMPLICONS_FOLDER,
                                            U_SUBS_RATE,U_INS_RATE,U_DEL_RATE,R_SUBS_RATE,R_INS_RATE,R_DEL_RATE,DEL_LENGTH_GEOMETRIC_PARAMETER,INS_MAX_LENGTH,
                                            SUBS_VAF_DIRICLET_PARAMETER,INS_VAF_DIRICLET_PARAMETER,DEL_VAF_DIRICLET_PARAMETER,
                                            R_SUBS_VAF_DIRICLET_PARAMETER,R_INS_VAF_DIRICLET_PARAMETER,R_DEL_VAF_DIRICLET_PARAMETER)
               
        if amplicons=="No":
            if VERBOSE:
                logging.info(f"No PCR error was introduced! Possible reason: too low error rates.")
            
            amplicons = list(df_amplicons["amplicon_filepath"])
            n_reads = list(df_amplicons["n_reads"])
        else:
            with open(f"{OUTPUT_FOLDER}/{OUTPUT_FILENAME_PREFIX}_PCR_errors.vcf","w") as o:
                o.write("##fileformat=VCFv4.3\n")
                o.write("##reference=MN908947.3\n")
                o.write('##contig=<ID=MN908947.3,length=29903>\n')
                o.write('##INFO=<ID=VAF,Number=A,Type=Float,Description="Variant Allele Frequency">\n')
                o.write('##INFO=<ID=REC,Number=A,Type=String,Description="Recurrence state across source genomes. R: recurrent; U: unique to genome">\n')
                o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

            vcf_errordf.to_csv(f"{OUTPUT_FOLDER}/{OUTPUT_FILENAME_PREFIX}_PCR_errors.vcf",
                                mode="a",header=False,index=False,sep="\t", float_format='%.5f')
            if VERBOSE:
                logging.info(f'All aimed PCR errros are written to "{OUTPUT_FOLDER}/{OUTPUT_FILENAME_PREFIX}_PCR_errors.vcf"')
    
    amplicons = [join(AMPLICONS_FOLDER, a) for a in amplicons]

    #merge art_illumina runs which have the same read count to optimize

    merged_n_reads=list(set(n_reads))
    if merged_n_reads[0]==0:
        merged_n_reads=merged_n_reads[1:]
        
    merged_amplicons=[join(AMPLICONS_FOLDER,f"merged_amplicon_rcount_{a}.fasta") for a in merged_n_reads]
    
    for readcount,m_amplicon in zip(merged_n_reads,merged_amplicons):
        with open(m_amplicon, "w") as merged:
            for amp in [amplicons[idx] for idx,r in enumerate(n_reads) if r==readcount]:
                with open(amp, "r") as amp_file:
                    shutil.copyfileobj(amp_file,merged)
    
    # STEP 4: Simulate Reads
    logging.info("Generating reads using art_illumina, cycling through all genomes and remaining amplicons.")
    with art_illumina(OUTPUT_FOLDER, OUTPUT_FILENAME_PREFIX, READ_LENGTH, SEQ_SYS,VERBOSE,TEMP_FOLDER,N_READS) as art:
        art.run(merged_amplicons, merged_n_reads)

    # STEP 5: Clean up all of the temp. directories
    for directory in [GENOMES_FOLDER, AMPLICONS_FOLDER, INDICES_FOLDER]:
        logging.info(f"Removing all files in {directory}")
        i = "y"

        if not AUTOREMOVE:
            logging.info(f"Press y and enter if you are ok with all files in the directory {directory}" +
            " being deleted (use flag --autoremove to stop showing this message).")
            i = input()

        if i.lower() == "y":
            files = glob.glob(join(directory, "*"))
            for f in files:
                os.remove(f)
