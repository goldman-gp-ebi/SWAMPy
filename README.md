# SARS-CoV-2 Metagenomics Simulator

This project is intended to simulate Sars-CoV-2 metagenomes taken from wastewater samples. 
Synthetic mixtures of amplicons are produced, based on proportions of viral genomes that are
supplied by the user. 

## Installation

Clone this repository and install the python dependancies `pandas` and `biopython`. 

In addition, the software uses the `bowtie2` aligner to extract amplicons from the viral genomes, based on the primers specified using
the `--primers` command line option (defaults to using ARTIC v3 primers). It also uses the `art_illumina` command line tool to simulate
reads from an Illumina sequencing device. 
You need to ensure that both of these commands (`bowtie2` and `art_illumina`) are available from your command line 
(i.e. both of the binaries of these tools are available from your `$PATH` environment variable). 
The simplest way to do this on a Debian-based system, with python and pip already installed, is below: 

```
pip install pandas, biopython
sudo apt-get install art_illumina, bowtie2
```

## Quickstart

```

# You only need to make these directories once.
cd src
mkdir ../example/amplicons
mkdir ../example/indices
mkdir ../simulation_output

# Run the example simulation.
# This creates a synthetic metagenome from the fasta files in the example/genomes folder, using relative genome proportions 
# from the example/abundances.csv, with a random seed of 10. 
# The primers used are from the ARTIC protocol V3, though the primers marked alt have been removed. 

python simulate_metagenome.py 

# Run the same example simulation again with all the parameters explicitly defined. 

python simulate_metagenome.py \
    --genomes_folder ../example/genomes \
    --amplicons_folder ../example/amplicons \
    --indices_folder  ../example/indices \
    --genome_abundances  ../example/abundances.csv \
    --primers_file  ../example/artic_sars-cov-2_primers_no_alts.fastq \
    --output_folder  ../simulation_output \
    --output_filename_prefix my_example \
    --n_reads  100000 \
    --seed 10 \
    --verbose  True \
    --amplicon_distribution  dirichlet_1 \
    --amplicon_distribution_file ../example/amplicon_distribution.csv \
    --amplicon_pseudocounts 10000 

# my_example1.fastq and my_example2.fastq should appear in simulation_output folder, as well as my_example_amplicon_abundances_summary.csv
# While running, some tmp files might appear in your working directory, but they will get cleaned up when the program terminates
# (even if it exits with an error).

```


## How it works

1. Read input files

    - genomes folder - one genome per fasta file

    - abundances.csv - in the format name,relative_abundance (see example/abundance.csv)
        These are relative abundances of viruses, not amplicons. 

    - primers - a fastq file containing all the primers used, the format of each read should be:
        @name_{amplicon number}_{LEFT|RIGHT}
        {primer sequence}
        +
        {dummy scores (ignored)} 
        (see example/artic_sars-cov-2_primers_no_alts.fastq)
        
    - amplicon_distribution.csv - in the format amplicon_number,hyperparameter
        (see example/amplicon_distribution.csv)
        Note that for the Dirichlet distributions, it is expected that the sum of the hyperparameter column is 1. 

    - folders for amplicons and indices:
        The software uses bowtie2 which creates a small index for each genome, these are stored in the folder marked --indices_folder. 
        Each ampilcon is extracted from each genome and stored a fasta file - these are collectively stored in the folder marked --amplicons_folder. 

    - other options (e.g. nreads) are command-line options

2. Create amplicon population

    - the set of primers are aligned with each genome, using bowtie2
    - the output sam file is parsed, and the matching positions and primer lengths, and primer names, are stored in a pandas dataframe. 
    - The above dataframe is separated out into two subframes, LEFT and RIGHT. These subframes are then joined based on amplicon number. 
    - At this stage, any primers that couldn't be aligned to the genome are dropped from the dataframe. 
    - The remaining amplicons have good matches for both their left and right primers - these are extracted from the genomes and stored in fasta files. 
    - The dataframe is appeneded to a 'main' table containing summary information about all amplicons from all genomes processed so far.

3. Simulating numbers of reads per amplicon
    - The amplicon_distribution.csv file is read, the hyperparameter vector (vector a) is stored. 
    - Vector a is expected at this stage to satisfy sum(a) == 1. 
    - Next, a is scaled by the amplicon_pseudocounts, c. I.e. a = a*c
    - Sample an amplicon abundance profile using p = Dir(a), Dirichlet distribution.
    - Each amplicon in the amplicon population dataframe is then given a number of reads according to the formula N_i = Binomial(N_genome, p_i), 
    where i is an index running over the amplicons (the amplicon number) and N_genome = int(NREADS * genome_proportion).
    - All of the above data is adjoined to the amplicon dataframe from step 2. 

4. For each row in the amplicon dataframe, run art_illumina with the required numbers of reads and with the following parameters:

    ```
        art_illumina \
            --amplicon \
            --quiet \
            --paired \
            --rndSeed random_number \
            --noALN \
            --maskN 0  \
            --seqSys", MSv3 \
            --in amplicon_file.fasta \
            --len 250 \
            --rcount n_reads \
            --out out_prefix 

    ```

    All of these output files are then stitched together in alphabetical order for each set of reads 1 and 2. 

## TODOs

- We should be able to pass the software one very big fasta file full of lots of genomes, and make a mixture using possibly a subset - at the moment you need to input one genome per fasta file and put them all in a folder, so this needs to be changed

- DONE Ensure repeatability with same random seeds, including in Nick's shuffling method (therefore might need to change it a little bit)

- Deal with alts - produce all combinations of alts and at a relatively lower proportion (1/2^n)

- Tune bowtie2 so that amplicon dropout occurs in line with experiment

- Tune bowtie2 to allow ACGT -> N substitutions in primer sites, at the moment only single N insertions are allowed

- What to do with genome ends (amplicon 1 and 98)? ATM they are dropped fairly consistantly because the leftmost and rightmost primer sites don't exist on most genomes that we look at except for the Wuhan reference. 

- Simulate other PCR products (how? Chimeras -> Simera + Point mutations -> with a script? /ignore Chimeras for now.)

- genome_abundances.csv should be able to cope with multiple rows / abundances pointing to the same genome. This is actually going to be a bit annoying to fix because the way it's done now, you'd end up with multiple rows with exactly the same id in the fastq. 

- Test the 'custom' amplicon distribution method (Nicola's method) make sure it works as intended - this is already implemented but needs testing. Additionally on this point, Nicola asked for multinomial sampling from each set of reads and also for the ability to read from a SAM or BAM file, this is currently not supported. 

- Generate a log file which saves the input parameters.

- Try to run this program on EBI resources rather than on my laptop and uploading.

- Give Nick a demo of the working software.

- Fix the other sampling forms, at the moment just DIRICHLET_1 works (and possibly Nicola's custom distribution). 

- Cruddiness parameter should be able to be different for each genome, and a place for that information is in the tsv file containing the
genome abundances. If there's a command line value, use that for all the genomes, if there isn't then 

- Random seed, if not given as a command line input make sure it's recorded somewhere in a log file. 

- Command line switch to delete the command line stuff. 

## Aknowledgements

We are using [bowtie2](bowtie-bio.sourceforge.net/bowtie2/index.shtml), [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) and [biopython](https://biopython.org/), as well as other python packages (pandas and numpy)
and are using a primer set coming from [ARTIC v3](https://artic.network/). 
The genomes in the example are courtesy of [Twist Bioscience](https://www.twistbioscience.com/), these genomes are in turn copies of surveillance sequences from [COG UK](https://www.cogconsortium.uk/). 