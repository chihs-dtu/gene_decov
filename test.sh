#!/bin/bash

set +e


GMT_PATH="/Users/chiaoyuhsieh//Downloads/data/enrichr.KEGG_2016.gmt"
OUT_PATH="/Users/chiaoyuhsieh//Downloads/data/test_out"
mkdir -p $OUT_PATH

# Install packages
Rscript -e 'source("install_packages.R", encoding = "UTF-8")'
pip3 install gseapy pyreadr


# Test R scripts
for n_sample in 100 1000 10000 100000; do
	echo "-------------------------------------------------------"
	echo "Testing with ${n_sample} samples..."
    echo "-------------------------------------------------------"
	MAT_PATH="/Users/chiaoyuhsieh//Downloads/data/subsetted_data_${n_sample}.rds"
	for TOOL in ssgsea reset; do
		echo "    Testing ${TOOL}..."
		#/usr/bin/time -v \
		/usr/bin/time -lp \
			Rscript test.R ${MAT_PATH} ${GMT_PATH} ${TOOL} ${OUT_PATH}
		echo ""
	done
	for TOOL in gseapy; do
		echo "    Testing ${TOOL}..."
		#/usr/bin/time -v \
		/usr/bin/time -lp \
			python3 run_gseapy.py -d ${MAT_PATH} --gmt ${GMT_PATH} \
								-n ${TOOL} -o ${OUT_PATH}
		echo ""
	done
	echo "-------------------------------------------------------"
done
