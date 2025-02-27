import os
import datetime
import logging
import pickle
import pyreadr
import seaborn as sns
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from scipy.stats import spearmanr, norm

np.random.seed(42)

# Configure the logger
current_date = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
logging.basicConfig(
    level=logging.DEBUG,  # Set the minimum log level
    format="%(asctime)s - %(levelname)s - %(message)s",  # Define the log format
    datefmt="%Y/%m/%d %H:%M:%S",  # Customize the date format
    handlers=[
        logging.FileHandler(f"run_{current_date}.log", mode='w'),
        logging.StreamHandler()]
)

# Create a logger instance
logger = logging.getLogger(__name__)

def load_gene_set(filename):
    """
    Load the GMT file and return a dictionary with gene set name as key and gene list as value.
    """
    logger.info(f"Loading gene set: {filename}")
    gene_sets = {}
    with open(filename, 'r') as f:
        for line in f:
            val = [ x for x in line.strip().split('\t') if len(x) > 1 ]
            gene_sets[val[0]] = val[1:]
    return gene_sets

def load_result(filename):
    """
    Load the result file and return a dictionary with file name as key and dataframe as value.
    Each row represents a gene set and each column represents a sample.
    """
    logger.info(f"Loading result file: {filename}")
    if filename.endswith('.csv'):
        if "gseapy" in filename:
            # columns: Name,Term,ES,NES
            df = pd.read_csv(filename)
            # Reshaping to create dataframe1 (using 'ES' values)
            df_ES = df.pivot(index='Term', columns='Name', values='ES')
            # Reshaping to create dataframe2 (using 'NES' values)
            df_NES = df.pivot(index='Term', columns='Name', values='NES')
            return {"gseapy:ES": df_ES, "gseapy:NES": df_NES}
        else:
            df = pd.read_csv(filename, index_col=0)
            
    elif filename.endswith('.rds'):
        df = pyreadr.read_r(filename)[None]

    name = filename.rsplit('_',1)[-1][:-4]
    return {name: df}


def visualize_dist(scores=list(), cutoff=None, filename="spearman.png"):
    
    plt.figure(figsize=(10, 6))

    sns.histplot(scores[0], kde=True, color='red', label='In GS', stat='density', bins=20)
    sns.histplot(scores[1], kde=True, color='blue', label='Not in GS', stat='density', bins=20)

    # Customize plot
    plt.title("Score Distribution for Two Conditions")
    plt.xlabel("Score")
    plt.ylabel("Density")
    plt.legend()
    plt.grid(True)

    # Add cutoff if necessary
    if cutoff is not None:
        plt.axvline(cutoff, color='green', linestyle='--', label=f'Optimal Cutoff: {cutoff:.2f}')

    # Show plot
    plt.savefig(filename)


def find_cutoff(scores=list()):
    # 2. Find the optimal cutoff by minimizing the overlap
    # Fit normal distributions to both conditions
    mu1, std1 = norm.fit(scores[0])
    mu2, std2 = norm.fit(scores[1])
    # logger.info("mean1={}, std1={}".format(mu1, std1))
    # logger.info("mean2={}, std2={}".format(mu2, std2))

    # Define the range of possible cutoff values (e.g., within the range of the data)
    cutoff_range = np.linspace(min(np.min(scores[0]), np.min(scores[1])), 
                            max(np.max(scores[0]), np.max(scores[1])), 
                            1000)

    # Calculate the probability density for both distributions at each point
    pdf1 = norm.pdf(cutoff_range, mu1, std1)  # PDF of condition 1
    pdf2 = norm.pdf(cutoff_range, mu2, std2)  # PDF of condition 2

    # Find the points where density_A is approximately equal to density_B
    # This is done by finding where the difference is close to zero
    tolerance = 1e-3  # You can adjust the tolerance based on your data
    equal_indices = np.where(pdf1 - pdf2 < tolerance)[0]

    # Now, if there are such points, select the one with the maximum density
    if len(equal_indices) > 0:
        # Find the index with the maximum density value
        max_density_index = equal_indices[np.argmax(pdf1[equal_indices])]
        cutoff = cutoff_range[max_density_index]
    else:
        logger.warning("No point where pdf1 and pdf2 are approximately equal within the given tolerance: {}".format(tolerance))
        cutoff = min(np.min(scores[0]), np.min(scores[1]))
    return cutoff


def visualize_cutoff_boxplot(cutoffs: dict, filename="cutoff_boxplot.png"):
    plt.figure(figsize=(10, 6))

    values_to_plot = {}
    for tool, df in cutoffs.items():
        values_to_plot[tool] = df["cutoff"].values
    
    sns.boxplot(values_to_plot, orient='v')
    
    # Customize plot
    plt.title("Correlation Cutoff Distribution Per Gene-Set")
    plt.xlabel("Tool")
    plt.ylabel("Cutoff Score")
    plt.legend()
    plt.grid(True)

    # Show plot
    plt.savefig(filename, dpi=100)


def compare(n_sample, dirname, gene_sets, fig_dir, return_abs=False):
    # Load Gene Expression Profile
    # Each row represents a gene and each column represents a sample.
    mat_exp = pyreadr.read_r(f"{dirname}/subsetted_data_{n_sample}.rds")[None]
    
    # Load Results
    res = {}
    for tool in ["gseapy", "reset", "ssgsea"]:
        filename = os.path.join(
            dirname, "test_out", f"subsetted_data_{n_sample}.rds_{tool}.csv")
        if not os.path.exists(filename):
            filename = filename[:-3] + "rds"
            if not os.path.exists(filename):
                logger.error(f"File not found: {filename}")
                continue
        res_dict = load_result(filename)
        for key, val in res_dict.items():
            res[key] = val

    # Calculate Spearman correlation between gene expression profile and each result
    #logger.info(f"***n_sample, tool, gene_set, spearman_r, is_random")
    rscores = {}
    cutoffs = {}

    dict_random_genes = defaultdict(list)

    for tool, mat_score in res.items():
        cutoff_per_gs = []
        logger.info(f"Calculating scores for {tool}...")
        rscores[tool] = {}
        rscores[tool]["default"] = []
        rscores[tool]["random"] = []

        sample_list = mat_score.columns.tolist()
        logger.info(f"Calculating Spearman correlation for {tool}")
        constant_genes = mat_exp.loc[mat_exp.nunique(axis=1) == 1].index.tolist()

        for gene_set, gene_list in gene_sets.items():
            # Filter out constant genes and those already in gene_list
            gene_list = [x for x in gene_list if x in mat_exp.index and x not in constant_genes]
            # Randomly select the same number of genes not included in the gene set as there are in the gene set
            if len(dict_random_genes[gene_set]) == 0:
                available_genes = [x for x in mat_exp.index.tolist() if x not in gene_list and x not in constant_genes]
                dict_random_genes[gene_set] = np.random.choice(
                    available_genes, size=len(gene_list), replace=False)

            r_arr = []
            for gene in gene_list:
                arr_exp = mat_exp.loc[gene, sample_list].values.flatten()
                arr_score = mat_score.loc[gene_set, :].values.flatten()
                spearmanr_value = spearmanr(arr_exp, arr_score)[0]
                #logger.info(f"***{n_sample}, {tool}, {gene_set}, {spearmanr_value}, false")
                if return_abs:
                    spearmanr_value = abs(spearmanr_value)
                r_arr.append(spearmanr_value)
            rscores[tool]["default"].append(r_arr)

            r_arr_rand = []
            for gene in dict_random_genes[gene_set]:
                arr_exp_random = mat_exp.loc[gene, sample_list].values.flatten()
                spearmanr_value_random = spearmanr(arr_exp_random, arr_score)[0]
                #logger.info(f"***{n_sample}, {tool}, {gene_set}, {spearmanr_value_random}, true")
                if return_abs:
                    spearmanr_value = abs(spearmanr_value)
                r_arr_rand.append(spearmanr_value_random)
            rscores[tool]["random"].append(r_arr_rand)
            
            # Find the best cutoff
            scores = [np.array(rscores[tool]["default"][-1]), np.array(rscores[tool]["random"][-1])]
            best_cutoff = find_cutoff(scores)
            cutoff_per_gs.append((gene_set, best_cutoff))
            
        
        cutoffs[tool] = pd.DataFrame.from_records(
                            cutoff_per_gs, columns=['gene_set', 'cutoff'])

    visualize_cutoff_boxplot(cutoffs, f"{fig_dir}/cutoff_boxplot_{n_sample}.png")
    return rscores, cutoffs, dict_random_genes


if __name__ == '__main__':
    dirname = "/home/people/chihs/Gene_Deconvolution/data/"
    gene_sets = load_gene_set("/home/people/chihs/Gene_Deconvolution/GSEApy/tests/extdata/enrichr.KEGG_2016.gmt")

    # Create output directory if not exists
    if not os.path.exists(os.path.join(dirname, "test_out", "cor_png")):
        os.makedirs(os.path.join(dirname, "test_out", "cor_png"))
    fig_dir = os.path.join(dirname, "test_out", "cor_png")

    res = {}
    return_abs = True
    for n_sample in [10000]:
        logger.info(f"Calculating {n_sample} samples ...")
        res[n_sample] = {}
        res[n_sample]["rscores"], res[n_sample]["cutoffs"], res[n_sample]["random_genes"] = compare(n_sample, dirname, gene_sets, fig_dir, return_abs=return_abs)
    
    if return_abs == True:
        outname = f"{dirname}/test_out/corr_return_abs.pickle"
    else:
        outname = f"{dirname}/test_out/corr.pickle"
    with open(outname, "wb") as f:
        pickle.dump(res, f)