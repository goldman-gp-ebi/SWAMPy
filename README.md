# SWAMPy

This project is intended to simulate SARS-CoV-2 metagenomes taken from wastewater samples. Synthetic mixtures of amplicons are produced, based on proportions of viral genomes that are supplied by the user and a supported primer set of choice. See [our method](https://github.com/goldman-gp-ebi/sars-cov-2-metagenomic-simulator/wiki/SWAMPy-method) to learn about how it works.


<img width="2560" alt="image" src="https://github.com/goldman-gp-ebi/SWAMPy/assets/57137586/5d3f67c2-c697-465d-8978-5ca45c4c0d32">




## Installation

1. Clone this repository. 

```
git clone https://github.com/goldman-gp-ebi/SWAMPy.git

cd SWAMPy
```
2. Install dependencies `pandas` and `biopython`.
You need to ensure `bowtie2` and `art_illumina` are available from your command line (i.e. both of the binaries of these tools are available from your `$PATH` environment variable). The program also requires Python version 3.9.7.

A way to do this on Ubuntu 20.04, with Python 3.9.7 and pip already installed, is below: 

```
pip install pandas, biopython
sudo apt-get install art-nextgen-simulation-tools, bowtie2
```

Or you can also use conda as follows to create a new environment with the dependencies and correct versions (RECOMMENDED):

```
conda env create -f SWAMPy-env.yaml
conda activate SWAMPy
```

**MAC USERS** also need:
```
brew install coreutils
```

**If you have a Mac with an M1 apple chip, the procedure is more annoying, sorry!**
In this circumstance, bioconda doesn't seem to work. Instead, the following process 
seems to work instead (though not everything ended up inside a conda environment...):

```
conda env create -f SWAMPy-env.yaml
# This should successfully install numpy, biopython, pandas, and gsl (2.7).
# If not then download them using pip. 
brew install coreutils
brew install bowtie2
# Install these two things system-wide(!).
```

Next, manually install [art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) from source. 
Download  ART-src-MountRainier-2016.06.05-MacOS.tgz (this is the current latest version as of December 2022, and has been since 2016), unzip the archive and run

```
./configure && make
```

and then move the `art_illumina` binary into wherever the conda environment lives (`conda env list` prints the path to the folder, inside which is a `bin` directory). 

Optionally, run `alias swampy='python /path/to/SWAMPy/src/simulate_metagenome.py'` to use `swampy` as a command rather than `python simulate_metagenome.py`. 


## Docker Installation

The latest SWAMPy Docker image can be found on Docker Hub at https://hub.docker.com/r/wboulton12/swampy, and pulled through the command 

```
docker pull wboulton12/swampy:latest
```

To run SWAMPy through Docker, mount the volumes from the input files (genomes and relative abundances) into the `/swampy_in` directory in the virtual filesystem, and mount the directory where output files are expected to `/swampy_out`. The command below is an example, where the output files are sent to the present working directory:

```
docker run -u $(id -u) -v $(pwd):/swampy_out -v /home/will/Desktop/SWAMPy/docker_example:/swampy_in wboulton12/swampy:latest --genomes_file /swampy_in/my_genomes.fasta --genome_abundances /swampy_in/my_abundances.tsv

# docker run (not as root, as your user), 
# (mount current directory to swampy_out, mount the folder containing genomes and abundance files to swampy_in),
# name of container,
# command line arguments passed to swampy. 
# Note that since the input and output folders are already mapped, you should not pass a --output_folder flag to swampy. 
```

## Singularity Installation

Build the singularity image file after pulling SWAMPy from the Docker Hub:

```
singularity pull docker://wboulton12/swampy:latest
```

Then execute the sif file (use exec rather than run):

```
singularity exec swampy_latest.sif python3 /home/src/simulate_metagenome.py \
    --genomes_file /path/to/genomes.fasta \
    --genome_abundances /path/to/abundances.tsv \
    --temp_folder /some/temporary/folder/you/can/write/to \
    --output_folder $(pwd) 
```

## Quickstart

Example [input files](https://github.com/goldman-gp-ebi/SWAMPy/wiki/SWAMPy-method#1-read-input-files) i.e. genomes.fasta and abudances.tsv are already included in the example directory.

### Run the example simulation

This creates a synthetic metagenome from the fasta files in the example/genomes folder, using relative genome proportions from the example/abundances.tsv. The primers used are from the ARTIC protocol v1; compared to v3, the primers marked as 'alt' are removed. For details, please see the potential bugs section below.

### Run the simulator with default parameters
```
cd src
```

```
python simulate_metagenome.py --primer_set a1 --genomes_file ../example/genomes.fasta --genome_abundances ../example/abundances.tsv --output_folder ../simulation_output --temp_folder ../example/temp/ --autoremove
```

### See the help page 
```
python simulate_metagenome.py --help
```

### Run the same example simulation again with all the parameters explicitly defined, and with the random seed set. 
```
python simulate_metagenome.py \
    --genomes_file ../example/genomes.fasta \
    --temp_folder ../example/temp/
    --genome_abundances ../example/abundances.tsv \
    --primer_set a1 \
    --output_folder ../simulation_output \
    --output_filename_prefix example \
    --n_reads  100000 \
    --seqSys MSv3 \
    --read_length 250 \
    --seed 10 \
    --amplicon_distribution  dirichlet_1 \
    --amplicon_pseudocounts 200 \
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
example1.fastq and example2.fastq should appear in simulation_output folder, as well as example_amplicon_abundances_summary.tsv, example.log and example_PCR_erros.vcf. While running, some tmp files might appear in your working directory, but they will get cleaned up when the program terminates (even if it exits with an error).

## Output files:
- example_R1.fastq & example_R2.fastq: simulated reads.
- example_amplicon_abundances_summary.tsv: a table summarising the amplicons.
- example_PCR_erros.vcf: all the intended high-frequency errors. Observed VAFs may be different from those in the VCF file due to randomness and recurrence.
- example.log: The log file.

## What these CLI arguments mean:

See the [CLI arguments](https://github.com/goldman-gp-ebi/SWAMPy/wiki/CLI-arguments) page.

## Extra options and potential bugs

Things to watch out for:

- If you want to run multiple instances of SWAMPy simultaneously, make sure to use a unique `--TEMP_FOLDER` for each run. Otherwise, they will interfere with each other. Submitting multiple SLURM or LSF jobs or using SWAMPy in a Snakemake rule are examples of situations where this warning may apply. 

- First and last amplicons get dropped for basically every genome except the Wuhan reference.
This is because the leftmost primer 1 and rightmost primer 98 basically never match in a genome. 
Also watch out for long runs of N's in the primer sites, if these are there then that amplicon will drop out. 

- If a source genome amplicon is has deletions at the end, PCR errors that were supposed to be on those loci will be skipped only for that souce genome. A warning is printed in this case.

- Special characters in the fasta genome ids could potentially cause a problem. In the code, the characters "/", "," and " " are dealt with,
but other whitespace characters or special characters such as "&" could cause a bug. 

- Note that the reads in the fastq files are both shuffled by a randomly chosen permutation. 

- The option 'dirichlet_1' relates to how numbers of reads for each genome and for each amplicon are drawn. At the moment
this is done using multiple Multinomial(N_genome, p) draws, where p is drawn once from a Dirichlet(\alpha * c) distribution. 
Here c is a scalar, the amplicon_pseudocounts number, and \alpha is a vector, scaled to have length 1, goverened by the amplicon_distribution file. N_genome is the number of reads for each genome, this is taken from the genome_abundances file. In the future, slightly different methods
of taking these draws may be possible, so this method is called 'dirichlet_1'. 

- Alt primers are only really relevant to **ARTIC v3**. ARTIC v3 has extra primers which are tagged as ALT - these are 'alternative' primer ends for certain amplicons where the original amplicons were dropping out a lot. These 'alt' primers are added to the primer pool along with the original primers. Therefore for these amplicons (where alt primers exist), there are a total of four possible versions of the amplicon that could be made, depending on whether the left/right end was used with the alt primer or not. In the original version of this code, we would only simulate two of these possibilities --- cases where both ends used the original primer and cases where both ends used the alt primer. However since then, because ARTIC v4 (which we assume will end up superceding ARTIC v3) was released, and has no alts (and also most other sequencing protocols have no alts either) we decided to ignore alt primers altogether, and our version of ARTIC does not include any simulation of alt primers. If you want this feature please contact us.

## Citation
If you used SWAMPy, please cite the publication 
[SWAMPy: Simulating SARS-CoV-2 Wastewater Amplicon Metagenomes with Python](https://doi.org/10.1101/2022.12.10.519890).

## Change log
- Allow to pass --qshift options to ART to turn off read errors / shift quality scores. 
- Allow for amplicon fragmentation with mean fragment length (and standard deviation) given by applying the options such as --fragment_amplicons, --fragment_len_mean=150, --fragment_len_sd=66. 
- As of 07365c7, --primer_set option is required from the user and has no default to emphasise the simulation's dependence on the specific primer set used.
- As of e7ae4ec, --disallowed_positions has been added to stop High-frequency error introduction from specific loci.
- As of 0cf5f16, SWAMPy accepts files for custom primer sets (see CLI options) though we still encourage users to add primer set panels as in issue #2 and send us a PR for the benefit of the community.
- As of ce23577, Artic V5 is supported thanks to @j3551ca!

## Acknowledgements

We are using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) and [biopython](https://biopython.org/), as well as other python packages (pandas and numpy)
and are using a primer set coming from [ARTIC v1](https://artic.network/). 
The genomes in the example are courtesy of [Twist Bioscience](https://www.twistbioscience.com/); these genomes are in turn copies of surveillance sequences from [COG UK](https://www.cogconsortium.uk/). 


---
associated repos:
[Will's notes](https://github.com/goldman-gp-ebi/wastewater_Will-s_notes.git)
