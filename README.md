# SWAMPy

This project is intended to simulate Sars-CoV-2 metagenomes taken from wastewater samples. 
Synthetic mixtures of amplicons are produced, based on proportions of viral genomes that are
supplied by the user. 

## Installation

Clone this repository and install the python dependencies `pandas` and `biopython`. 

```
git clone https://github.com/goldman-gp-ebi/sars-cov-2-metagenomic-simulator.git
```

In addition, the software uses the `bowtie2` aligner to extract amplicons from the viral genomes, based on the primers specified using
the `--primer_set` command line option (defaults to using ARTIC v1 primers). It also uses the `art_illumina` command line tool to simulate
reads from an Illumina sequencing device. 
You need to ensure that both of these commands (`bowtie2` and `art_illumina`) are available from your command line 
(i.e. both of the binaries of these tools are available from your `$PATH` environment variable). 
The simplest way to do this on a Debian-based system, with python and pip already installed, is below: 

```
pip install pandas, biopython
sudo apt-get install art_illumina, bowtie2
```

Or you can also use conda as follows:

```
conda create -c bioconda -c anaconda -n SWAMPy biopython art bowtie2 pandas 
conda activate SWAMPy
```

You will also need to place a multi-fasta 'genomes.fasta' file in the example folder - we can send this file to you
but since the sequences are from GISAID we have not uploaded these to github. 

## Quickstart

```
# You only need to make these directories once.

cd src

mkdir ../example/genomes
mkdir ../example/amplicons
mkdir ../example/indices
mkdir ../simulation_output
```
### Run the example simulation.

This creates a synthetic metagenome from the fasta files in the example/genomes folder, using relative genome proportions from the example/abundances.tsv. The primers used are from the ARTIC protocol v1; compared to v3, the primers marked as 'alt' are removed. For details, please see the potential bugs section below.

### Run the simulator with default parameters

```
python simulate_metagenome.py 
```

### See the help page 
```
python simulate_metagenome.py --help
```

### Run the same example simulation again with all the parameters explicitly defined, and with the random seed set. 
```
python simulate_metagenome.py \
    --genomes_file ../example/genomes.fasta \
    --genomes_folder ../example/genomes \
    --amplicons_folder ../example/amplicons \
    --indices_folder ../example/indices \
    --genome_abundances ../example/abundances.tsv \
    --primer_set a1 \
    --output_folder ../simulation_output \
    --output_filename_prefix example \
    --n_reads  100000 \
    --seqSys MSv3 \
    --read_length 250 \
    --seed 10 \
    --amplicon_distribution  dirichlet_1 \
    --amplicon_pseudocounts 10000 \
    --unique_insertion_rate 0.00002 \
    --unique_deletion_rate 0.000115 \
    --unique_substitution_rate 0.002485 \
    --recurrent_insertion_rate 0.00002 \
    --recurrent_deletion_rate 0 \
    --recurrent_substitution_rate 0.003357 \
    --max_insertion_length 14 \
    --subs_VAF_alpha 0.36,0.27 \
    --del_VAF_alpha 0.59,0.42 \
    --ins_VAF_alpha 0.33,0.45 \
    --r_subs_VAF_alpha 0.36,0.27 \
    --r_del_VAF_alpha 0.59,0.42 \
    --r_ins_VAF_alpha 0.33,0.45 
```
example1.fastq and example2.fastq should appear in simulation_output folder, as well as example_amplicon_abundances_summary.tsv, example.log and example_PCR_erros.vcf While running, some tmp files might appear in your working directory, but they will get cleaned up when the program terminates (even if it exits with an error).


## What these CLI arguments mean:

- **genomes_file**: this is the multi-fasta file containing all of the genomes in the pool of genomes that might be sampled from. 
The fasta id strings in this should not have any tabs in them, and ideally shouldn't have any special characters, though /, |, and spaces are ok. 

- **genomes_folder**, amplicons_folder, indices_folder: these are folders used to store temporary files (cleaned up before the program exits). 
WARNING: DO NOT STORE ANYTHING IMPORTANT IN THESE FOLDERS!

- **genome_abundances**: a tab delimited file with one line per genome to be sampled from; each line has the format:
genome_id   relative_genome_abundance

- **primer_set**: Can be one of a1(Artic version 1 or Artic version 3 without alts), a4 (Artic version4) or n2 (Nimagen version 2)

- **output_folder**:  folder where the 5 output files will go

- **output_filename_prefix**: name for the output files, which could produce e.g. my_file_1.fastq, my_file_2.fastq, my_file.log,
and my_file_amplicon_abundances_summary.tsv (a tsv output summarising statistics about amplicon counts and hyperparameters per genome).

- **n_reads**: a target number of reads - the true number of reads may be lower if some amplicons drop out

- **seqSys**: name of the sequencing system, options to use are given by the art_illumina help text, and are:
    GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
    HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
    HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
    MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp). Default is MSv3

- **read_length**: lengths of the reads coming off the sequencer (default is 250). 

- **seed**: seed for the random number generator (optional; if not set, a random seed is set and recorded in the log file)

- **quiet**: Supress verbose output for both the log file and the terminal stdout. 

- **autoremove**: Delete temproray files under genomes folder, amplicon folder and indices folder after simulation is completed.

- **amplicon_distribution**: can be one of 'custom', 'dirichlet_1' (the default), or 'dirichlet_2'. Only use dirichlet_1 at present. 

- **amplicon_pseudocounts**: scaling factor for the above vector - higher corresponds to higher 'quality', lower will have more variation in amplicon
distributions between different runs. A reasonable value for this seems to be 10000 (the default). 

- **unique_insertion_rate**:Genome-wide unique high-frequency unique error rate. See PCR error section for details about error types. Default is 0.00002. 

- **unique_deletion_rate**:Genome-wide unique high-frequency unique error rate. See PCR error section for details about error types. Default is 0.000115. 

- **unique_substitution_rate**:Genome-wide unique high-frequency unique error rate. See PCR error section for details about error types. Default is 0.002485. 

- **recurrent_insertion_rate**:Genome-wide high-frequency recurrent error rate. See PCR error section for details about error types. Default is 0.00002.

- **recurrent_deletion_rate**:Genome-wide high-frequency recurrent error rate. See PCR error section for details about error types. Default is 0.

- **recurrent_substitution_rate**:Genome-wide high-frequency recurrent error rate. See PCR error section for details about error types. Default is 0.003357.

- **max_insertion_length**: Maximum length of a high- frequency insertion. Default is 14 .

- **subs_VAF_alpha**:Comma-seperated alpha1 and alpha2 values for the Dirichlet distribution to sample allele frequencies of unique high-frequency errors. (The first value of the sample is for the alternative alelle). Default is 0.36,0.27. 

- **del_VAF_alpha**:Comma-seperated alpha1 and alpha2 values for the Dirichlet distribution to sample allele frequencies of unique high-frequency errors. (The first value of the sample is for the alternative alelle). Default is 0.59,0.42. 

- **ins_VAF_alpha**:Comma-seperated alpha1 and alpha2 values for the Dirichlet distribution to sample allele frequencies of unique high-frequency errors. (The first value of the sample is for the alternative alelle). Default is 0.33,0.45. 

- **r_subs_VAF_alpha**:Comma-seperated alpha1 and alpha2 values for the Dirichlet distribution to sample allele frequencies of recurrent high-frequency errors. (The first value of the sample is for the alternative alelle). Default is equal to unique errors. 

- **r_del_VAF_alpha**:Comma-seperated alpha1 and alpha2 values for the Dirichlet distribution to sample allele frequencies of recurrent high-frequency errors. (The first value of the sample is for the alternative alelle). Default is equal to unique errors. 

- **r_ins_VAF_alpha**:Comma-seperated alpha1 and alpha2 values for the Dirichlet distribution to sample allele frequencies of recurrent high-frequency errors. (The first value of the sample is for the alternative alelle). Default is equal to unique errors .


## Extra options and potential bugs

Things to watch out for:

- First and last amplicons get dropped for basically every genome except the Wuhan reference.
This is because the leftmost primer 1 and rightmost primer 98 basically never match in a genome. 
Also watch out for long runs of N's in the primer sites, if these are there then that amplicon will drop out. 

- Special characters in the fasta genome ids could potentially cause a problem. In the code, the characters "/" and " " are dealt with,
but other whitespace characters or special characters such as "&" could cause a bug. 

- Note that the reads in the fastq files are both shuffled by a randomly chosen permutation. 

- The option 'dirichlet_1' relates to how numbers of reads for each genome and for each amplicon are drawn. At the moment
this is done using multiple Multinomial(N_genome, p) draws, where p is drawn once from a Dirichlet(\alpha * c) distribution. 
Here c is a scalar, the amplicon_pseudocounts number, and \alpha is a vector, scaled to have length 1, goverened by the amplicon_distribution file. N_genome is the number of reads for each genome, this is taken from the genome_abundances file. In the future, slightly different methods
of taking these draws may be possible, so this method is called 'dirichlet_1'. 

- Alt primers are only really relevant to **Artic v3**. Artic v3 has extra primers which are tagged as ALT - these are 'alternative' primer ends for certain amplicons where the original amplicons were dropping out a lot. I believe that these 'alt' primers are added to the primer pool along with the original primers. Therefore for these amplicons (where alt primers exist), there are a total of 4 possible versions of the amplicon that could be made, depending on whether the left/right end was used with the alt primer or not. In the original version of this code, we would only simulate 2 of these possibilities - cases where both ends used the original primer and cases where both ends used the alt primer. However since then, because Artic v4 (which we assume will end up superseding artic v3) was released, and has no alts (and also most other sequencing protocols have no alts either) we decided to ignore alt primers altogether, and our version of artic does not include any simulation of alt primers. If you want this feature please contact us


## How it works

1. Read input files
    **User-given**
    - genome_file - a multi-fasta containing all of the genomes that may be sampled. 

    - abundances.tsv - in the format genome_fasta_id \t relative_abundance (see example/abundance.tsv)
        These are relative abundances of viruses, not amplicons. 

    - folders for amplicons and indices:
        The software uses bowtie2 which creates a small index for each genome, these are stored in the folder marked --indices_folder. 
        Each ampilcon is extracted from each genome and stored as a fasta file - these are collectively stored in the folder marked --amplicons_folder. 

    **primer-set specific input files. Already provided within SWAMPy**

    -   primers - a fastq file containing all the primers used, the format of each read should be:
        @name_{amplicon number}_{LEFT|RIGHT}
        {primer sequence}
        +
        {dummy quality scores (ignored)} 
        (see example/artic_v3_primers_no_alts.fastq)

    - Also a BED file of the primers, where positions are with respect to the Wuhan reference.
        
    - amplicon_distribution.tsv - in the format amplicon_number \t hyperparameter
        (see example/amplicon_distribution.csv)
        Note that for the Dirichlet distributions, it is expected that the sum of the hyperparameter column is 1. 

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
    - Each amplicon in the amplicon population dataframe is then given a number of reads according by drawing from Multinomial(N_genome, p), 
    where N_genome = int(NREADS * genome_proportion).
    - All of the above data is adjoined to the amplicon dataframe from step 2. 

4. High-Frequency error introduction
    - Counts of each type of errors to be introduced are drawn from a Poisson distribution using error rate x genome length as the mean.
    - VAF of each error is drawn from a Dirichlet distribution. (alt,ref)=Dir(alpha1,alpha2)
    - A dataframe which contains error type, position, corresponding amplicon number, alternative allele, recurrency state, VAF and length information of each error is created.
    - Amplicons from step 3 is aligned to the Wuhan reference.
    - Mutants of the amplicons are created using different combinations of errors that are on the same amplicon.
    - Read count of the mutants are determined with a binomial distribution, whose p is VAF of the error.
    - Amplicon dataframe is updated so that it now contains the PCR mutants and corresponding read counts.

5. Grouping amplicons that have the same number of readcounts and merging their FASTA files.

6. For each group, run art_illumina with the required numbers of reads. At present we simulate with parameters appropriate for Illumina MiSeq v3, paired-end, 250bp, using the following parameters:

    ```
        art_illumina \
            --amplicon \
            --quiet \
            --paired \
            --rndSeed random_number \
            --noALN \
            --maskN 0  \
            --seqSys seqSys(default MSv3) \
            --in amplicon_file.fasta \
            --len 250 \
            --rcount n_reads \
            --out out_prefix    
    ```

    All of these output files are then stitched together and shuffled into a random order by a randomly generated permutation. 

7. All of the temporary files (genomes, amplicons, bowtie indices, and art output files) are deleted. 


## TODOs

- Deal with alts - produce all combinations of alts and at a relatively lower proportion (1/2^n)

- Tune bowtie2 so that amplicon dropout occurs in line with experiment

- Tune bowtie2 to allow ACGT -> N substitutions in primer sites, at the moment only single N insertions are allowed

- What to do with genome ends (amplicon 1 and 98)? ATM they are dropped fairly consistantly because the leftmost and rightmost primer sites don't exist on most genomes that we look at except for the Wuhan reference. 

- Simulate other PCR products (how? Chimeras -> Simera + Point mutations -> with a script? /ignore Chimeras for now.)

- Test the 'exact' amplicon distribution method (Nicola's method) make sure it works as intended - this is already implemented but needs testing. Additionally on this point, Nicola asked for multinomial sampling from each set of reads and also for the ability to read from a SAM or BAM file, this is currently not supported. 

- Cruddiness parameter should be able to be different for each genome, and a place for that information is in the tsv file containing the
genome abundances. If there's a command line value, use that for all the genomes, if there isn't then look in the genome abundances, otherwise use a default.  

## Aknowledgements

We are using [bowtie2](bowtie-bio.sourceforge.net/bowtie2/index.shtml), [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) and [biopython](https://biopython.org/), as well as other python packages (pandas and numpy)
and are using a primer set coming from [ARTIC v1](https://artic.network/). 
The genomes in the example are courtesy of [Twist Bioscience](https://www.twistbioscience.com/), these genomes are in turn copies of surveillance sequences from [COG UK](https://www.cogconsortium.uk/). 
