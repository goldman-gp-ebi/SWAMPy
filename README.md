# SARS-CoV-2 Metagenomics Simulator

This project is intended to simulate Sars-CoV-2 metagenomes taken from wastewater samples. 
Synthetic mixtures of amplicons are produced, based on proportions of viral genomes that are
supplied by the user. 

## Installation

Clone this repository and install the python dependancies `pandas` and `biopython`. 

In addition, the software uses the `bowtie2` aligner to extract amplicons from the viral genomes, based on the primers specified using
the `--primers` command line option (defaults to using ARTIC v1 primers). It also uses the `art_illumina` command line tool to simulate
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
mkdir ../example/genomes
mkdir ../example/amplicons
mkdir ../example/indices
mkdir ../simulation_output

# Run the example simulation.
# This creates a synthetic metagenome from the fasta files in the example/genomes folder, using relative genome proportions 
# from the example/abundances.tsv.
# The primers used are from the ARTIC protocol v1; compared to v3, the primers marked as 'alt' are removed. 

python simulate_metagenome.py 

# Run the same example simulation again with all the parameters explicitly defined, and with the random seed set. 

python simulate_metagenome.py \
    --genomes_file ../example/genomes.fasta \
    --genomes_folder ../example/genomes \
    --amplicons_folder ../example/amplicons \
    --indices_folder  ../example/indices \
    --genome_abundances  ../example/abundances.tsv \
    --primers_file  ../example/artic_sars-cov-2_primers_no_alts.fastq \
    --output_folder  ../simulation_output \
    --output_filename_prefix example \
    --n_reads  100000 \
    --seed 10 \
    --verbose  True \
    --amplicon_distribution  dirichlet_1 \
    --amplicon_distribution_file ../example/amplicon_distribution.tsv \
    --amplicon_pseudocounts 10000 

# my_example1.fastq and my_example2.fastq should appear in simulation_output folder, as well as my_example_amplicon_abundances_summary.tsv
# While running, some tmp files might appear in your working directory, but they will get cleaned up when the program terminates
# (even if it exits with an error).

```

What these CLI arguments mean:

- genomes_file: this is the multi-fasta file containing all of the genomes in the pool of genomes that might be sampled from. 
The fasta id strings in this should not have any tabs in them, and ideally shouldn't have any special characters, though /, |, and spaces are ok. 

- genomes_folder, amplicons_folder, indices_folder: these are folders used to store temporary files (cleaned up before the program exits). 
WARNING: DO NOT STORE ANYTHING IMPORTANT IN THESE FOLDERS!

- genome_abundances: a tab delimited file with one line per genome to be sampled from; each line has the format:
genome_id   relative_genome_abundance

- primers_file: a fastq of all the primers, formatted so that left and right primers can be matched by name after aligning. The default is artic v1. 

- output_folder:  folder where the 4 output files will go

- output_filename_prefix: name for the output files, which could produce e.g. my_file_1.fastq, my_file_2.fastq, my_file.log,
and my_file_amplicon_abundances_summary.tsv (a tsv output summarising statistics about amplicon counts and hyperparameters per genome).

- n_reads: a target number of reads - the true number of reads may be lower if some amplicons drop out

- seed: seed for the random number generator (optional; if not set, a random seed is set and recorded)

- verbose: Verbose mode for both the log file and the terminal output. 

- amplicon_distribution: can be one of 'custom', 'dirichlet_1' (the default), or 'dirichlet_2'. Only use dirichlet_1 at present. 

- amplicon_distribution_file: a file containing a vector for the Dirchlet distribution governing the amplicon distribution. 

- amplicon_pseudocounts: scaling factor for the above vector - higher corresponds to higher 'quality', lower will have more variation in amplicon
distributions between different runs. A reasonable value for this seems to be 10000 (the default). 




## Extra options and potential bugs

Things to watch out for:

- First and last amplicons get dropped for basically every genome except the Wuhan reference.
This is because the leftmost primer 1 and rightmost primer 98 basically never match in a genome. 
Also watch out for long runs of N's in the primer sites, if these are there then that amplicon will drop out. 

- Special characters in the fasta genome ids could potentially cause a problem. In the code, the characters "/" and " " are dealt with,
but other whitespace characters or special characters could cause a bug. 

- Note that the reads in the fastq files are both shuffled by a randomly chosen permutation. 

- The option 'dirichlet_1' relates to how numbers of reads for each genome and for each amplicon are drawn. At the moment
this is done using multiple Multinomial(N_genome, p) draws, where p is drawn once from a Dirichlet(\alpha * c) distribution. 
Here c is a scalar, the amplicon_pseudocounts number, and \alpha is a vector, scaled to have length 1, goverened by the amplicon_distribution file. N_genome is the number of reads for each genome, this is taken from the genome_abundances file. In the future, slightly different methods
of taking these draws may be possible, so this method is called 'dirichlet_1'. 


## How it works

1. Read input files

    - genome_file - a multi-fasta containing all of the genomes that may be sampled. 

    - abundances.tsv - in the format genome_fasta_id \t relative_abundance (see example/abundance.tsv)
        These are relative abundances of viruses, not amplicons. 

    - primers - a fastq file containing all the primers used, the format of each read should be:
        @name_{amplicon number}_{LEFT|RIGHT}
        {primer sequence}
        +
        {dummy quality scores (ignored)} 
        (see example/artic_sars-cov-2_primers_no_alts.fastq)
        
    - amplicon_distribution.tsv - in the format amplicon_number \t hyperparameter
        (see example/amplicon_distribution.csv)
        Note that for the Dirichlet distributions, it is expected that the sum of the hyperparameter column is 1. 

    - folders for amplicons and indices:
        The software uses bowtie2 which creates a small index for each genome, these are stored in the folder marked --indices_folder. 
        Each ampilcon is extracted from each genome and stored as a fasta file - these are collectively stored in the folder marked --amplicons_folder. 

    - other options (e.g. nreads) are command-line options

2. Create amplicon population

    - the set of primers are aligned with each genome, using bowtie2
    - the output sam file is parsed, and the matching positions and primer lengths, and primer names, are stored in a pandas dataframe. 
    - The above dataframe is separated out into two subframes, LEFT and RIGHT. These subframes are then joined based on amplicon number. 
    - At this stage, any primers that couldn't be aligned to the genome are dropped from the dataframe. 
    - The remaining amplicons have good matches for both their left and right primers - these are extracted from the genomes and stored in fasta files. 
    - The dataframe is appended to a 'main' table containing summary information about all amplicons from all genomes processed so far.

3. Simulating numbers of reads per amplicon
    - The amplicon_distribution.csv file is read, the hyperparameter vector (vector a) is stored. 
    - Vector a is expected at this stage to satisfy sum(a) == 1. 
    - Next, a is scaled by the amplicon_pseudocounts, c. I.e. a = a*c
    - Sample an amplicon abundance profile using p = Dir(a), Dirichlet distribution.
    - Each amplicon in the amplicon population dataframe is then given a number of reads according to the formula N_i = Binomial(N_genome, p_i), 
    where i is an index running over the amplicons (the amplicon number) and N_genome = int(NREADS * genome_proportion).
    - All of the above data is adjoined to the amplicon dataframe from step 2. 

4. For each row in the amplicon dataframe, run art_illumina with the required numbers of reads. At present we simulate with parameters appropriate for Illumina MiSeq v3, paired-end, 250bp, using the following parameters:

    ```
        art_illumina \
            --amplicon \
            --quiet \
            --paired \
            --rndSeed random_number \
            --noALN \
            --maskN 0  \
            --seqSys MSv3 \
            --in amplicon_file.fasta \
            --len 250 \
            --rcount n_reads \
            --out out_prefix 

    ```

    All of these output files are then stitched together and shuffled into a random order by a randomly generated permutation. 

5. All of the temporary files (genomes, amplicons, bowtie indices, and art output files) are deleted. 


## TODOs

- Deal with alts - produce all combinations of alts and at a relatively lower proportion (1/2^n)

- Tune bowtie2 so that amplicon dropout occurs in line with experiment

- Tune bowtie2 to allow ACGT -> N substitutions in primer sites, at the moment only single N insertions are allowed

- What to do with genome ends (amplicon 1 and 98)? ATM they are dropped fairly consistantly because the leftmost and rightmost primer sites don't exist on most genomes that we look at except for the Wuhan reference. 

- Simulate other PCR products (how? Chimeras -> Simera + Point mutations -> with a script? /ignore Chimeras for now.)

- genome_abundances.tsv should be able to cope with multiple rows / abundances pointing to the same genome. This is actually going to be a bit annoying to fix because the way it's done now, you'd end up with multiple rows with exactly the same id in the fastq. 

- Test the 'exact' amplicon distribution method (Nicola's method) make sure it works as intended - this is already implemented but needs testing. Additionally on this point, Nicola asked for multinomial sampling from each set of reads and also for the ability to read from a SAM or BAM file, this is currently not supported. 

- Try to run this program on EBI resources rather than on my laptop and uploading.

- Fix the other sampling forms, at the moment just DIRICHLET_1 works (and possibly Nicola's custom distribution). 

- Cruddiness parameter should be able to be different for each genome, and a place for that information is in the tsv file containing the
genome abundances. If there's a command line value, use that for all the genomes, if there isn't then look in the genome abundances, otherwise use a default. 

- Command line switch to delete the temporary folders. 

- Keep track of progress on the bowtie2 stuff so that you can see your progress e.g. "genome 1/3 completed."

- For the loop of reads, each time it's on the outer loop (of genomes) keep track. 

- In the logging remove the dataframe head parts. 

## Aknowledgements

We are using [bowtie2](bowtie-bio.sourceforge.net/bowtie2/index.shtml), [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) and [biopython](https://biopython.org/), as well as other python packages (pandas and numpy)
and are using a primer set coming from [ARTIC v1](https://artic.network/). 
The genomes in the example are courtesy of [Twist Bioscience](https://www.twistbioscience.com/), these genomes are in turn copies of surveillance sequences from [COG UK](https://www.cogconsortium.uk/). 