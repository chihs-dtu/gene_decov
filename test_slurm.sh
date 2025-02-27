#!/bin/bash
# Job name:
#SBATCH --job-name=gene_deconv
#
# Partition:
#SBATCH --partition=cpu
#
# Request one node:
#SBATCH --nodes=1
#
# Specify memory for the job (example):
#SBATCH --mem=150G
# Exclude nodes
#SBATCH --exclude=node11
#
# Processors per task:
#SBATCH --cpus-per-task=4
#
# File for output, use file or /dev/null
#SBATCH -o log/%J.out
#SBATCH -e log/%J.err


GMT_PATH="/home/people/chihs/Gene_Deconvolution/GSEApy/tests/extdata/enrichr.KEGG_2016.gmt"
OUT_PATH="/home/people/chihs/Gene_Deconvolution/data/test_out"
mkdir -p $OUT_PATH

CUR_PATH=$(pwd)


# Activate the virtual invironment
# conda activate r_env

for n_sample in 100 1000 10000 100000; do
	echo "-------------------------------------------------------"
	echo "Testing with ${n_sample} samples..."
    echo "-------------------------------------------------------"
	MAT_PATH="/home/people/chihs/Gene_Deconvolution/data/subsetted_data_${n_sample}.rds"
	for TOOL in ssgsea reset; do
		echo "    Testing ${TOOL}..."
		#/usr/bin/time -v \
		/usr/bin/time -lp \
			Rscript ${CUR_PATH}/test.R ${MAT_PATH} ${GMT_PATH} ${TOOL} ${OUT_PATH}
		echo ""
	done
	for TOOL in gseapy gsva; do
		echo "    Testing ${TOOL}..."
		#/usr/bin/time -v \
		/usr/bin/time -lp \
			python3 run_gseapy.py -d ${MAT_PATH} --gmt ${GMT_PATH} \
								-n ${TOOL} -o ${OUT_PATH}
		echo ""
	done
	echo "-------------------------------------------------------"
done


