from numpy.random import dirichlet, binomial, multinomial
import numpy as np
import pandas as pd
import logging


def exact_sampler(
    amplicon_df: pd.DataFrame,
    amplicon_distribution_file: str,
) -> pd.DataFrame:
    """
    Determines number of reads per amplicon as listed in `amplicon_distribution_file`:
    | amplicon_number | Genome_1 | Genome_2 | Genome_3 |
    | --------------: | -------: | -------: | -------: |
    | 1               | 100      | 300      | 400      |
    | 2               | 240      | 450      | 300      |
    | 3               | 240      | 360      | 100      |
    | 3               | 230      | 300      | 0        |
    | 3               | 240      | 0        | 100      |
    | 3               | 0        | 300      | 10       |
    """

    amp_dist_df = pd.read_csv(amplicon_distribution_file, sep="\t")

    amp_dist_df.columns = [
        # process the lineage names the same as before
        col.replace(" ", "_")
        if col != "amplicon_number"
        else col for col in amp_dist_df.columns
    ]
    lineages = [col for col in amp_dist_df.columns if col != "amplicon_number"]

    amp_dist_df_long = amp_dist_df.melt(
        id_vars="amplicon_number", var_name="ref",
        value_vars=lineages, value_name="n_reads"
    )

    amplicon_df = amplicon_df.merge(
        amp_dist_df_long,
        on=["ref", "amplicon_number"],
        how="left"
    )

    amplicon_df["total_n_reads"] = amp_dist_df[lineages].sum().sum()
    return amplicon_df


def get_amplicon_reads_sampler(amplicon_distribution, amplicon_distribution_file, amplicon_pseudocounts_c, genome_abundances, total_n_reads):
    
    if amplicon_distribution.upper() == "EXACT":

        def hyperparam_sampler(dataframe_row):
            return -1

        def genome_count_sampler(dataframe_row):
            return -1

        def prob_sampler(dataframe_row):
            return -1

        def reads_sampler(dataframe_row):
            return -1


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