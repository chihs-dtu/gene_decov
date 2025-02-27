import gseapy as gp
import argparse
import os
import pyreadr
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Add arguments
    parser.add_argument('-d', '--data', type=str, help='Path to gene expression data.')
    parser.add_argument('--gmt', type=str, help='Path to gene set.')
    parser.add_argument('-o', '--outdir', type=str, help='Path to outdir.')
    parser.add_argument('-n', '--method', type=str, help='Options: gsva | ssgsea')

    # Parse the arguments
    args = parser.parse_args()

    # Read file
    if args.data.endswith(".rds"):
        data = pyreadr.read_r(args.data)[None]
    else:
        data = pd.read_csv(args.data, index_col=0)
    basename = os.path.basename(args.data)

    # Run
    if args.method == "gsva":
        es = gp.gsva(data=data, gene_sets=args.gmt, outdir=args.outdir, max_size=3000, min_size=3)
        os.rename("{}/gseapy.gene_set.gsva.report.csv".format(args.outdir), 
                    "{}/{}_gsva.csv".format(args.outdir, basename))

    elif args.method == "gseapy":
        ss = gp.ssgsea(data=data, gene_sets=args.gmt, outdir=args.outdir,
               sample_norm_method='rank', no_plot=True,
               max_size=3000, min_size=3)
        os.rename("{}/gseapy.gene_set.ssgsea.report.csv".format(args.outdir), 
                    "{}/{}_gseapy.csv".format(args.outdir, basename))

    else:
        print("Invalid method. Please choose either 'gsva' or'ssgsea'.")
        exit(1)
