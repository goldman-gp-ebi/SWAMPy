import shutil
import subprocess
import glob
import os
from os.path import basename, join
import numpy as np
from contextlib import contextmanager
import string
import logging


class ArtIllumina:

    def __init__(
            self, 
            outpath, 
            output_filename_prefix, 
            read_length, 
            seq_sys, 
            qprof1, 
            qprof2,
            verbose,temp, 
            nreads, 
            fragment_amplicons, 
            fragment_len_mean, 
            fragment_len_sd, 
            qshift):
        
        self.outpath = outpath
        self.output_filename_prefix = output_filename_prefix
        self.read_length = read_length
        self.seq_sys = seq_sys
        self.qprof1 = qprof1 
        self.qprof2 = qprof2
        self.verbose = verbose
        self.temp = temp
        self.nreads = nreads
        self.fragment_amplicons = fragment_amplicons
        self.fragment_len_mean = fragment_len_mean
        self.fragment_len_sd = fragment_len_sd
        self.qshift = qshift


    def run_once(self, infile, n_reads, out_prefix, rnd_seed):
        art_options = [
            "art_illumina", 
        ]
        if self.fragment_amplicons:
            art_options += ["--mflen", str(self.fragment_len_mean), "--sdev", str(self.fragment_len_sd)]
        else:
            art_options += ["--amplicon"]

        if self.seq_sys.lower() == "custom":
            if (not self.qprof1) or (not self.qprof2):
                logging.error("If you supply --seqSys custom then you must supply --qprof1 and --qprof2 files")
                exit("If you supply --seqSys custom then you must supply --qprof1 and --qprof2 files. Exiting.")

            art_options += [
                "--qprof1", self.qprof1,
                "--qprof2", self.qprof2
            ]
        
        else:
            art_options += [ "--seqSys", self.seq_sys ]

            
        art_options += [
            "--paired",
            "--rndSeed", str(rnd_seed),
            "--noALN",
            "--maskN", str(0), 
            "--in", infile,
            "--len", str(self.read_length),
            "--rcount", str(n_reads),
            "--out", out_prefix 
        ]

        if self.qshift != 0:
            art_options += ["--qShift", str(self.qshift), "--qShift2", str(self.qshift)]

        op = subprocess.run(art_options, capture_output=True)

        message_lines = op.stdout.decode("ASCII").split("\n")[-4:-2]
        warning = op.stderr.decode("ASCII")
        if self.verbose:
            for line in message_lines:
                logging.info("art_illumina: " + line)
        if warning != "Warning: your simulation will not output any ALN or SAM file with your parameter settings!\n":
            logging.warning(warning)

        

    def run(self, amplicons, n_reads):

        params = zip(amplicons, n_reads)

        for a, n in params:
            
            short_name = ".".join(basename(a).split(".")[:-1])
            
            if self.verbose:
                logging.info(f"Starting on file {short_name}.fasta with {n} reads")

            self.run_once(a, n, join(self.temp,"tmp.sms.")+short_name+".", np.random.randint(2 ** 63))


        all_r1_files = sorted([x for x in glob.glob(join(self.temp,"tmp.sms.*")) if x[-4:] == "1.fq"])
        all_r2_files = sorted([x for x in glob.glob(join(self.temp,"tmp.sms.*")) if x[-4:] == "2.fq"])

        with open(join(self.temp,"tmp.sms.all_files_unshuffled1.fastq"), "w") as all_r1:
            for r1 in all_r1_files:
                with open(r1, "r") as r1fh:
                    shutil.copyfileobj(r1fh, all_r1)
        
        with open(join(self.temp,"tmp.sms.all_files_unshuffled2.fastq"), "w") as all_r2:
            for r2 in all_r2_files:
                with open(r2, "r") as r2fh:
                    shutil.copyfileobj(r2fh, all_r2)

        logging.info("Creating random data for shuffle.")
        with open(join(self.temp,"tmp.sms.random_data"), "w") as random_data:
            ALPHABET = np.array(list(string.ascii_lowercase))
            random_data.write("".join(np.random.choice(ALPHABET, size=max(5000000, int(2.5*self.nreads)))))

        # shuffle the fastq's so that the reads are in a random order. 
        shuffle_fastq_file(join(self.temp,"tmp.sms.all_files_unshuffled1.fastq"), join(self.outpath, f"{self.output_filename_prefix}_R1.fastq"), join(self.temp,"tmp.sms.random_data"))
        shuffle_fastq_file(join(self.temp,"tmp.sms.all_files_unshuffled2.fastq"), join(self.outpath, f"{self.output_filename_prefix}_R2.fastq"), join(self.temp,"tmp.sms.random_data"))


def shuffle_fastq_file(input_filename, output_filename, random_seed):
    logging.info(f"Shuffling {output_filename}")
    # additionally, this changes all '&' characters back to '/' characters. 
    os.system(f"paste -s -d '\t\t\t\n' {input_filename} | shuf --random-source={random_seed} | tr '\t&' '\n/' > {output_filename}")

@contextmanager
def art_illumina(
    outpath, 
    output_filename_prefix, 
    read_length, 
    seq_sys, 
    qprof1,
    qprof2,
    verbose, 
    temp, 
    nreads, 
    fragment_amplicons, 
    fragment_len_mean, 
    fragment_len_sd, 
    qshift):
    
    try:
        yield ArtIllumina(
            outpath, 
            output_filename_prefix, 
            read_length, 
            seq_sys, 
            qprof1,
            qprof2,
            verbose, 
            temp, 
            nreads, 
            fragment_amplicons, 
            fragment_len_mean, 
            fragment_len_sd, 
            qshift
        )
    
    finally:
        logging.info("Exiting sars-cov-2 metagenome simulator - tidying up.")

        for tem in glob.glob(join(temp,"tmp.sms.*")):
            os.remove(tem)
