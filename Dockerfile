FROM ubuntu:20.04

MAINTAINER goldman_lab

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
        art-nextgen-simulation-tools \
        bowtie2 \
        python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN pip install pandas==1.3.3 biopython==1.79

COPY . /home

RUN mkdir swampy_out && mkdir swampy_in

ENTRYPOINT [\
    "python3",\
    "/home/src/simulate_metagenome.py",\
    "--autoremove",\
    "--output_folder",\
    "swampy_out"\
]
