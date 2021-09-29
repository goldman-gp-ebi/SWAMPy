import argparse
from os.path import dirname, join, abspath, basename
import os
import glob
import logging
from Bio import SeqIO
import subprocess
import pandas as pd
import numpy as np

from art_runner import art_illumina
from create_amplicons import build_index, align_primers, write_amplicon
from read_model import get_amplicon_reads_sampler


# All these caps variables are set once (by user inputs, with default values) but then never touched again.
BASE_DIR = join(dirname(dirname(abspath(__file__))), "example")
GENOMES_FILE = join(BASE_DIR, "genomes.fasta")
GENOMES_FOLDER = join(BASE_DIR, "genomes")
AMPLICONS_FOLDER = join(BASE_DIR, "amplicons")
INDICES_FOLDER = join(BASE_DIR, "indices")
ABUNDANCES_FILE = join(BASE_DIR, "abundances.tsv")
PRIMERS_FILE = join(BASE_DIR, "artic_v3_primers_no_alts.fastq")
OUTPUT_FOLDER = os.getcwd()
OUTPUT_FILENAME_PREFIX = "example"
N_READS = 100000
READ_LENGTH = 250
SEQ_SYS = "MSv3"
SEED = np.random.randint(1000000000)
VERBOSE = True
AMPLICON_DISTRIBUTION = "DIRICHLET_1"
AMPLICON_DISTRIBUTION_FILE = join(BASE_DIR, "artic_v3_amplicon_distribution.tsv")
AMPLICON_PSEUDOCOUNTS = 10000
AUTOREMOVE = False

def setup_parser():
    parser = argparse.ArgumentParser(description="Run SARS-CoV-2 metagenome simulation.")
    parser.add_argument("--genomes_file", help="File containing all of the genomes that might be used", default=GENOMES_FILE)
    parser.add_argument("--genomes_folder", "-g", help="A temporary folder containing fasta files of genomes used in the simulation.", default=GENOMES_FOLDER)
    parser.add_argument("--amplicons_folder", "-am", help="A temporary folder that will contain amplicons of all the genomes.", default=AMPLICONS_FOLDER)
    parser.add_argument("--indices_folder", "-i", help="A temporary folder where bowtie2 indices are created and stored.", default=INDICES_FOLDER)
    parser.add_argument("--genome_abundances", "-ab", help="TSV of genome abundances.", default=ABUNDANCES_FILE)
    parser.add_argument("--primers_file", "-p", help="Path to fastq file of primers. Default ARTIC V1 primers.", default=PRIMERS_FILE)
    parser.add_argument("--output_folder", "-o", help="Folder where the output fastq files will be stored,", default=OUTPUT_FOLDER)
    parser.add_argument("--output_filename_prefix", "-x", help="Name of the fastq files name1.fastq, name2.fastq", default=OUTPUT_FILENAME_PREFIX)
    parser.add_argument("--seqSys", help="Name of the sequencing system, options to use are given by the art_illumina help text, and are:" + 
    """GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
           HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
           HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
           MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)""", default="MSv3")
    parser.add_argument("--n_reads", "-n", help="Approximate number of reads in fastq file (subject to sampling stochasticity).", default=N_READS)
    parser.add_argument("--read_length", "-l", help="Length of reads taken from the sequencing machine.", default=READ_LENGTH)
    parser.add_argument("--seed", "-s", help="Random seed", default=SEED)
    parser.add_argument("--verbose", "-vb", help="Verbose output", default=VERBOSE)
    parser.add_argument("--amplicon_distribution", default=AMPLICON_DISTRIBUTION)
    parser.add_argument("--amplicon_distribution_file", default=AMPLICON_DISTRIBUTION_FILE)
    parser.add_argument("--amplicon_pseudocounts","-c", default=AMPLICON_PSEUDOCOUNTS)
    parser.add_argument("--autoremove", default=AUTOREMOVE)
    
    return parser

def load_command_line_args():
    parser = setup_parser()
    args = parser.parse_args()

    global GENOMES_FILE
    GENOMES_FILE = args.genomes_file

    global GENOMES_FOLDER
    GENOMES_FOLDER = args.genomes_folder

    global AMPLICONS_FOLDER
    AMPLICONS_FOLDER = args.amplicons_folder

    global INDICES_FOLDER
    INDICES_FOLDER = args.indices_folder

    global ABUNDANCES_FILE 
    ABUNDANCES_FILE = args.genome_abundances

    global PRIMERS_FILE
    PRIMERS_FILE = args.primers_file

    global OUTPUT_FOLDER
    OUTPUT_FOLDER = args.output_folder

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

    global N_READS
    N_READS = int(args.n_reads)
    logging.info(f"Number of reads: {N_READS}")

    global READ_LENGTH
    READ_LENGTH = int(args.read_length)

    global SEED
    SEED = args.seed
    np.random.seed(int(SEED))
    logging.info(f"Random seed: {SEED}")

    global VERBOSE
    VERBOSE = args.verbose

    global AMPLICON_DISTRIBUTION
    AMPLICON_DISTRIBUTION = args.amplicon_distribution

    global AMPLICON_DISTRIBUTION_FILE
    AMPLICON_DISTRIBUTION_FILE = args.amplicon_distribution_file
    
    global AMPLICON_PSEUDOCOUNTS
    AMPLICON_PSEUDOCOUNTS = int(args.amplicon_pseudocounts)
    logging.info(f"Amplicon pseudocounts/ i.e. quality parameter: {AMPLICON_PSEUDOCOUNTS}")

    global AUTOREMOVE
    AUTOREMOVE = args.autoremove
    


if __name__ == "__main__":

    # STEP 0: Read command line arguments
    load_command_line_args()

    # STEP 1: Simulate Viral Population

    # Read genome abundances csv file
    genome_abundances = {}
    df_amplicons = pd.DataFrame()

    with open(ABUNDANCES_FILE) as ab_file:
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
    for genome in SeqIO.parse(GENOMES_FILE, format="fasta"):
        filepath = genome.description.replace(" ", "&").replace("/", "&")
        filepath += ".fasta"
        SeqIO.write(genome, join(GENOMES_FOLDER, filepath), format="fasta")


    # STEP 2: Simulate Amplicon Population
    genome_counter = 0
    for genome_path in genome_abundances:
        genome_counter += 1
        genome_path = genome_path.replace(" ", "&").replace("/", "&") + ".fasta"
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
    

        # STEP 3: Library Prep - PCR Amplification of Amplicons
        # Ignore this for now - skip to next step!


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

    if VERBOSE:
        logging.info(f"Total number of reads was {sum(df_amplicons['n_reads'])}, when {N_READS} was expected.")
    
    # STEP 4: Simulate Reads
    amplicons = [join(AMPLICONS_FOLDER, a) for a in df_amplicons["amplicon_filepath"]]
    n_reads = list(df_amplicons["n_reads"])

    logging.info("Generating reads using art_illumina, cycling through all genomes and remaining amplicons.")
    with art_illumina(OUTPUT_FOLDER, OUTPUT_FILENAME_PREFIX, READ_LENGTH, SEQ_SYS) as art:
        art.run(amplicons, n_reads)

    # STEP 5: Clean up all of the temp. directories
    for directory in [GENOMES_FOLDER, AMPLICONS_FOLDER, INDICES_FOLDER]:
        logging.info(f"Removing all files in {directory}")
        i = "y"
        
        if not AUTOREMOVE:
            logging.info(f"Press y and enter if you are ok with all files in the directory {directory}" +
            " being deleted (use option --autoremove True to stop showing this message).")
            i = input()

        if i.lower() == "y":
            files = glob.glob(join(directory, "*"))
            for f in files:
                os.remove(f)