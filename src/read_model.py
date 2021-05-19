from numpy.random import dirichlet, binomial
import numpy as np
import pandas as pd


def get_amplicon_reads_sampler(amplicon_distribution, amplicon_distribution_file, amplicon_pseudocounts_c, total_n_reads, seed):
    
    np.random.seed = seed 
    if amplicon_distribution.upper() == "EXACT":

        amplicon_distribution_dict = {}
        with open(amplicon_distribution_file) as adf:
            genome_names = adf.readline().split(",")
            amplicon_distribution_dict = {g: [] for g in genome_names}

        def hyperparam_sampler(dataframe_row):
            return -1

        def prob_sampler(dataframe_row):
            return -1

        def reads_sampler(dataframe_row):
            nonlocal amplicon_distribution_dict
            return amplicon_distribution_dict[dataframe_row.ref][dataframe_row.amplicon_number - 1]


    elif amplicon_distribution.upper() == "DIRICHLET_1":

        hyperparams = pd.read_csv(amplicon_distribution_file)
        hyperparams = {t.amplicon_number:t.hyperparameter for t in hyperparams.itertuples()}
        probs = dirichlet(np.array([hyperparams[i] for i in sorted(hyperparams.keys())], dtype=float) * float(amplicon_pseudocounts_c))

        def hyperparam_sampler(dataframe_row):
            nonlocal hyperparams
            return hyperparams[dataframe_row.amplicon_number]
        
        def prob_sampler(dataframe_row):
            nonlocal probs
            p = probs[dataframe_row.amplicon_number - 1]
            return p

        def reads_sampler(dataframe_row):
            d = dataframe_row
            return binomial(round(d.total_n_reads * d.abundance), d.amplicon_prob)

    elif amplicon_distribution.upper() == "DIRICHLET_2":
        
        genomes_list = {}
        hyperparams = pd.read_csv(amplicon_distribution_file)
        hyperparams = {t.amplicon_number:t.hyperparameter for t in hyperparams.itertuples()}
        hyperparams = np.array([hyperparams[i] for i in sorted(hyperparams.keys())], dtype=float)
        
        def hyperparam_sampler(dataframe_row):
            nonlocal hyperparams
            return hyperparams[dataframe_row.amplicon_number - 1]
        
        def prob_samplier(dataframe_row):
            nonlocal genomes_list
            nonlocal hyperparams
            nonlocal amplicon_pseudocounts_c
            if not dataframe_row.name in genomes_list:
                genomes_list[dataframe_row.name] = dirichlet(amplicon_pseudocounts_c * hyperparams)

            p = genomes_list[dataframe_row.name][dataframe_row.amplicon_number - 1]
            return p
        
        def reads_sampler(dataframe_row):
            d = dataframe_row
            return binomial(round(d.total_n_reads * d.abundance), d.amplicon_prob)

    else:
        print("Amplicon distribution not recognised, pick one of EXACT, DIRICHLET_1, DIRICHLET_2.")
        exit(1)

    return  hyperparam_sampler, prob_sampler, reads_sampler