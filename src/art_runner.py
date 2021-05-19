import subprocess
import shutil
import glob
import os
from os.path import basename, join
import numpy as np
from contextlib import contextmanager

class ArtIllumina:

    def __init__(self, rnd_seed, outpath, output_filename_prefix):
        self.outpath = outpath
        self.output_filename_prefix = output_filename_prefix
        np.random.seed = rnd_seed

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
            print(f"Starting on file {short_name}.fasta with {n} reads")
            self.run_once(a, n, "tmp.sms."+short_name+".", np.random.randint(2 ** 63))


        all_r1_files = sorted([x for x in glob.glob("./tmp.sms.*") if x[-4:] == "1.fq"])
        all_r2_files = sorted([x for x in glob.glob("./tmp.sms.*") if x[-4:] == "2.fq"])

        with open(f"./{self.output_filename_prefix}1.fastq", "w") as all_r1:
            for r1 in all_r1_files:
                with open(r1, "r") as r1fh:
                    shutil.copyfileobj(r1fh, all_r1)
        
        with open(f"./{self.output_filename_prefix}2.fastq", "w") as all_r2:
            for r2 in all_r2_files:
                with open(r2, "r") as r2fh:
                    shutil.copyfileobj(r2fh, all_r2)

@contextmanager
def art_illumina(rnd_seed, outpath, output_filename_prefix):
    
    try:
        yield ArtIllumina(rnd_seed, outpath, output_filename_prefix)
    
    finally:
        print("Exiting sars-cov-2 metagenome simulator - tidying up.")

        for temp in glob.glob("./tmp.sms.*"):
            os.remove(temp)