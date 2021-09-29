There are two scripts in here:
 - pull_data.sh is used to get all of the data from here: http://www.cgr.liv.ac.uk/illum/WWdump_e26afa198e791cb2/random_real_data/
 and get it onto my laptop. I am only trying to get the R1 and R2 fastq files out of the backup directory. 
 This script is much easier than trying to do it all by hand. 

 - do_alignments.sh is used to process each of these datasets one by one and produces a few sets of csvs - these are then used in the
 random_data_error_stats.ipynb notebook in the notes folder. These are error rates for 3 different things:
    - PCR error rates 
    - "High frequency" error rates (both PCR errors and SNPs)
    - sequencing errors (anything else)
    - I also make pandas dataframes which contain all the information from the vcfs, and save them too. 
    - After the bioinformatics stuff, the notebook which gets run is frequency_spectrum_cleaned.ipynb, 
    which has some notes in it. 