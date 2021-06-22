import shutil
import subprocess
import glob
import os
from os.path import basename, join
from io import StringIO
import numpy as np
from contextlib import contextmanager
import string
import logging


class ArtIllumina:

    def __init__(self, outpath, output_filename_prefix):
        self.outpath = outpath
        self.output_filename_prefix = output_filename_prefix

    def run_once(self, infile, n_reads, out_prefix, rnd_seed):

        subprocess.run([
            "art_illumina", 
            "--amplicon",
            "--quiet",
            "--paired",
            "--rndSeed", str(rnd_seed),
            "--noALN",
            "--maskN", str(0), 
            "--seqSys", "MSv3",
            "--in", infile,
            "--len", "250",
            "--rcount", str(n_reads),
            "--out", out_prefix 
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def run(self, amplicons, n_reads):

        params = zip(amplicons, n_reads)

        for a, n in params:
            
            short_name = ".".join(basename(a).split(".")[:-1])
            logging.info(f"Starting on file {short_name}.fasta with {n} reads")
            self.run_once(a, n, "tmp.sms."+short_name+".", np.random.randint(2 ** 63))


        all_r1_files = sorted([x for x in glob.glob("./tmp.sms.*") if x[-4:] == "1.fq"])
        all_r2_files = sorted([x for x in glob.glob("./tmp.sms.*") if x[-4:] == "2.fq"])

        with open("./tmp.sms.all_files_unshuffled1.fastq", "w") as all_r1:
            for r1 in all_r1_files:
                with open(r1, "r") as r1fh:
                    shutil.copyfileobj(r1fh, all_r1)
        
        with open("./tmp.sms.all_files_unshuffled2.fastq", "w") as all_r2:
            for r2 in all_r2_files:
                with open(r2, "r") as r2fh:
                    shutil.copyfileobj(r2fh, all_r2)

        logging.info("Creating random data for shuffle.")
        with open("./tmp.sms.random_data", "w") as random_data:
            ALPHABET = np.array(list(string.ascii_lowercase))
            random_data.write("".join(np.random.choice(ALPHABET, size=5000000)))

        # shuffle the fastq's so that the reads are in a random order. 
        shuffle_fastq_file("./tmp.sms.all_files_unshuffled1.fastq", join(self.outpath, f"{self.output_filename_prefix}_R1.fastq"), "./tmp.sms.random_data")
        shuffle_fastq_file("./tmp.sms.all_files_unshuffled2.fastq", join(self.outpath, f"{self.output_filename_prefix}_R2.fastq"), "./tmp.sms.random_data")


def shuffle_fastq_file(input_filename, output_filename, random_seed):
    logging.info(f"Shuffling {output_filename}")
    # additionally, this changes all '&' characters back to '/' characters. 
    os.system(f"paste -s -d '\t\t\t\n' {input_filename} | shuf --random-source={random_seed} | tr '\t&' '\n/' > {output_filename}")

@contextmanager
def art_illumina(outpath, output_filename_prefix):
    
    try:
        yield ArtIllumina(outpath, output_filename_prefix)
    
    finally:
        logging.info("Exiting sars-cov-2 metagenome simulator - tidying up.")

        for temp in glob.glob("./tmp.sms.*"):
            os.remove(temp)