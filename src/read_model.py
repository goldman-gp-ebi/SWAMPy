from numpy.random import dirichlet, binomial, multinomial
import numpy as np
import pandas as pd
import logging


def get_amplicon_reads_sampler(amplicon_distribution, amplicon_distribution_file, amplicon_pseudocounts_c, genome_abundances, total_n_reads):
    
    if amplicon_distribution.upper() == "EXACT":

        amplicon_distribution_dict = {}
        with open(amplicon_distribution_file) as adf:
            genome_names = adf.readline().split("\t")
            amplicon_distribution_dict = {g: [] for g in genome_names}

        def hyperparam_sampler(dataframe_row):
            return -1

        def genome_count_sampler(dataframe_row):
            return -1

        def prob_sampler(dataframe_row):
            return -1

        def reads_sampler(dataframe_row):
            nonlocal amplicon_distribution_dict
            return amplicon_distribution_dict[dataframe_row.ref][dataframe_row.amplicon_number - 1]


    elif amplicon_distribution.upper() == "DIRICHLET_1":


        genome_counts = multinomial(total_n_reads, [genome_abundances[i] for i in sorted(genome_abundances.keys())])
        genome_counts = {k:genome_counts[i] for i,k in enumerate(sorted(genome_abundances.keys()))}
        
        hyperparams = pd.read_csv(amplicon_distribution_file, sep="\t")
        hyperparams = {t.amplicon_number:t.hyperparameter for t in hyperparams.itertuples()}
        probs = dirichlet(np.array([hyperparams[i] for i in sorted(hyperparams.keys())], dtype=float) * float(amplicon_pseudocounts_c))
        
        amplicon_counts = {ref:multinomial(genome_counts[ref], probs) for ref in genome_abundances.keys()}


        def hyperparam_sampler(dataframe_row):
            nonlocal hyperparams
            return hyperparams[dataframe_row.amplicon_number]

        def genome_count_sampler(dataframe_row):
            nonlocal genome_counts
            return genome_counts[dataframe_row.ref]
        
        def prob_sampler(dataframe_row):
            nonlocal probs
            p = probs[dataframe_row.amplicon_number - 1]
            return p

        def reads_sampler(dataframe_row):
            d = dataframe_row
            return amplicon_counts[d.ref][d.amplicon_number - 1]

    elif amplicon_distribution.upper() == "DIRICHLET_2":
        
        genomes_list = {}
        genome_counts = multinomial(total_n_reads, [genome_abundances[i] for i in sorted(genome_abundances.keys())])
        genome_counts = {k:genome_counts[i] for i, k in enumerate(sorted(genome_abundances.keys()))}

        hyperparams = pd.read_csv(amplicon_distribution_file, sep="\t")
        hyperparams = {t.amplicon_number:t.hyperparameter for t in hyperparams.itertuples()}
        hyperparams = np.array([hyperparams[i] for i in sorted(hyperparams.keys())], dtype=float)
        
        def hyperparam_sampler(dataframe_row):
            nonlocal hyperparams
            return hyperparams[dataframe_row.amplicon_number - 1]

        def genome_count_sampler(dataframe_row):
            nonlocal genome_counts
            return genome_counts[dataframe_row.ref]
        
        def prob_sampler(dataframe_row):
            nonlocal genomes_list
            nonlocal hyperparams
            nonlocal amplicon_pseudocounts_c
            if type(amplicon_pseudocounts_c) == float:
                amplicon_pseudocounts_c = [amplicon_pseudocounts_c for i in range(len(genomes_list))]
            if not dataframe_row.ref in genomes_list:
                genomes_list[dataframe_row.ref] = dirichlet(amplicon_pseudocounts_c * hyperparams)

            p = genomes_list[dataframe_row.ref][dataframe_row.amplicon_number - 1]
            return p
        
        def reads_sampler(dataframe_row):
            d = dataframe_row
            return binomial(round(d.total_n_reads * d.abundance), d.amplicon_prob)

    else:
        logging.info("Amplicon distribution not recognised, pick one of EXACT, DIRICHLET_1, DIRICHLET_2.")
        exit(1)

    return  genome_count_sampler, hyperparam_sampler, prob_sampler, reads_sampler