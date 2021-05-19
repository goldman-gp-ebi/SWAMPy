import argparse
from os.path import dirname, join, abspath, basename
import os
import glob
import logging
from Bio import SeqIO
import subprocess
import pandas as pd
import shutil
from io import StringIO
from types import SimpleNamespace
import numpy as np
#0, 1, 5, 10, 20, 50, 80, 90, 95, 99, 100
from art_runner import art_illumina
from create_amplicons import build_index, align_primers, write_amplicon
from read_model import get_amplicon_reads_sampler

BASE_DIR = join(dirname(dirname(abspath(__file__))), "example")
GENOMES_FOLDER = join(BASE_DIR, "genomes")
AMPLICONS_FOLDER = join(BASE_DIR, "amplicons")
INDICES_FOLDER = join(BASE_DIR, "indices")
ABUNDANCES_FILE = join(BASE_DIR, "abundances.csv")
PRIMERS_FILE = join(BASE_DIR, "artic_sars-cov-2_primers_no_alts.fastq")
OUTPUT_FOLDER = os.getcwd()
OUTPUT_FILENAME_PREFIX = "S100_c10000_rep1_R"
N_READS = 100000
SEED = 21
VERBOSE = True
AMPLICON_DISTRIBUTION = "DIRICHLET_1"
AMPLICON_DISTRIBUTION_FILE = join(BASE_DIR, "amplicon_distribution.csv")
AMPLICON_PSEUDOCOUNTS = 10000

np.random.seed = SEED


log = logging.getLogger()

def setup_parser():
    parser = argparse.ArgumentParser(description="Run SARS-CoV-2 metagenome simulation.")
    parser.add_argument("--genomes_folder", "-g", help="Folder containing fasta files of genomes used in the simulation.", default=GENOMES_FOLDER)
    parser.add_argument("--amplicons_folder", "-am", help="Folder that will contain amplicons of all the genomes.", default=AMPLICONS_FOLDER)
    parser.add_argument("--indices_folder", "-i", help="Folder where bowtie2 indices are created and stored.", default=INDICES_FOLDER)
    parser.add_argument("--genome_abundances", "-ab", help="CSV of genome abundances.", default=ABUNDANCES_FILE)
    parser.add_argument("--primers_file", "-p", help="Path to fastq file of primers. Default ARTIC V3 primers.", default=PRIMERS_FILE)
    parser.add_argument("--output_folder", "-o", help="Folder where the output fastq files will be stored,", default=OUTPUT_FOLDER)
    parser.add_argument("--output_filename_prefix", "-x", help="Name of the fastq files name1.fastq, name2.fastq", default=OUTPUT_FILENAME_PREFIX)
    parser.add_argument("--nreads", "-n", help="Approximate number of reads in fastq file (subject to sampling stochasticity).", default=N_READS)
    parser.add_argument("--seed", "-s", help="Random seed", default=SEED)
    parser.add_argument("--verbose", "-vb", help="Verbose output", default=VERBOSE)
    parser.add_argument("--amplicon_distribution", default=AMPLICON_DISTRIBUTION)
    parser.add_argument("--amplicon_distribution_file", default=AMPLICON_DISTRIBUTION_FILE)
    parser.add_argument("--amplicon_pseudocounts","-c", default=AMPLICON_PSEUDOCOUNTS)
    
    return parser

def load_command_line_args():
    parser = setup_parser()
    args = parser.parse_args()
    args = parser.parse_args()

    global BASE_DIR
    BASE_DIR = args.workspace_folder

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

    global N_READS
    N_READS = args.nreads

    global SEED
    SEED = args.seed

    global VERBOSE
    VERBOSE = args.verbose

    global AMPLICON_DISTRIBUTION
    AMPLICON_DISTRIBUTION = args.amplicon_distribution

    global AMPLICON_DISTRIBUTION_FILE
    AMPLICON_DISTRIBUTION_FILE = args.amplicon_distribution_file
    
    global AMPLICON_PSEUDOCOUNTS
    AMPLICON_PSEUDOCOUNTS = args.amplicon_pseudocounts
    


if __name__ == "__main__":

    # STEP 0: Read command line arguments
    load_command_line_args()

    # STEP 1: Simulate Viral Population

    # Read genome abundances csv file
    genome_abundances = {}
    df_amplicons = pd.DataFrame()

    with open(ABUNDANCES_FILE) as ab_file:
        for line in ab_file:
            name, relative_abundance = tuple(line.split(","))
            genome_abundances[name] = float(relative_abundance)

    if abs(sum(genome_abundances.values()) - 1) > 0.000000001:
        total = sum(genome_abundances.values())
        if total <= 0:
            print(f"The total genome abundance is set to {total}, which is impossible.")
            exit(1)

        print(f"Total of relative abundance values is {total}, not 1.")
        print("Continuing, normalising total of genome abundances to 1.")
        
        for k in genome_abundances.keys():
            genome_abundances[k] /= total


    # STEP 2: Simulate Amplicon Population

    for genome_path in glob.iglob(join(GENOMES_FOLDER, "*")):

        genome_filename_short = ".".join(basename(genome_path).split(".")[:-1])
        reference = SeqIO.read(genome_path, format="fasta")

        # use bowtie2 to create a dataframe with positions of each primer pair aligned to the genome
        build_index(genome_path, genome_filename_short, INDICES_FOLDER)
        df = align_primers(genome_path, genome_filename_short, INDICES_FOLDER, PRIMERS_FILE, VERBOSE)        
        df["abundance"] = genome_abundances[df["ref"][0]]
        
        # write the amplicon to a file
        write_amplicon(df, reference, genome_filename_short, AMPLICONS_FOLDER)
    

        # STEP 3: Library Prep - PCR Amplification of Amplicons
        # Ignore this for now - skip to next step!


        df_amplicons = pd.concat([df_amplicons, df])


    # pick total numbers of reads for each amplicon
    amplicon_hyperparameter_sampler, amplicon_probability_sampler, amplicon_reads_sampler = get_amplicon_reads_sampler(
                                AMPLICON_DISTRIBUTION, 
                                AMPLICON_DISTRIBUTION_FILE, 
                                AMPLICON_PSEUDOCOUNTS, 
                                N_READS)

    df_amplicons["total_n_reads"] = N_READS
    
    df_amplicons["hyperparameter"] = df_amplicons.apply(amplicon_hyperparameter_sampler, axis=1)
    df_amplicons["amplicon_prob"] = df_amplicons.apply(amplicon_probability_sampler, axis=1)
    df_amplicons["n_reads"] = df_amplicons.apply(amplicon_reads_sampler, axis=1)
    
    if VERBOSE:
        print(df_amplicons.head())

    # write a summary csv
    df_amplicons[
        ["ref", 
        "amplicon_number", 
        "is_alt", 
        "amplicon_filepath",
        "total_n_reads", 
        "abundance",
        "hyperparameter",
        "amplicon_prob",
        "n_reads"]].to_csv(join(OUTPUT_FOLDER, f"{OUTPUT_FILENAME_PREFIX}_amplicon_abundances_summary.csv"))

    if VERBOSE:
        print(f"Total number of reads was {sum(df_amplicons['n_reads'])}, when {N_READS} was expected.")
    
    # STEP 4: Simulate Reads
    amplicons = [join(AMPLICONS_FOLDER, a) for a in df_amplicons["amplicon_filepath"]]
    n_reads = list(df_amplicons["n_reads"])

    with art_illumina(SEED, OUTPUT_FOLDER, OUTPUT_FILENAME_PREFIX) as art:
        art.run(amplicons, n_reads)