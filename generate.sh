#!/bin/bash

while getopts i:d:n:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        d) dataset=${OPTARG};;
        n) cell_no=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

module load gcc/8.2.0 intel/perflibs/2020_update4 tools R/4.0.0 scgan/20201002
pwd


setsid Rscript /ngc/people/viklav/bulk2sc/Bulk2SC/deconvolution.R $input $dataset $output
setsid Rscript /ngc/people/viklav/bulk2sc/Bulk2SC/cell_count_maker.R $cell_no $output

module purge
module load scgan/20201002 cuda90/toolkit/9.0.176 cudnn/7.4.2

ROOT_PATH="/ngc/people/viklav/bulk2sc/Bulk2SC/Datasets"
END_PATH="scGAN"
DATE=`date '+%Y-%m-%d'`


SAMPLE=1
while IFS="" read -r p || [ -n "$p" ]
do
  VALUE="$p"
  python /ngc/people/viklav/bulk2sc/Bulk2SC/scGAN/main.py --param $ROOT_PATH/$dataset/$END_PATH/parameters_original.json --generate --cells_no $VALUE --model_path $ROOT_PATH/$dataset/$END_PATH --save_path /ngc/people/viklav/bulk2sc/Bulk2SC/Results/$DATE-$output/generated_sample_$SAMPLE.h5ad
  ((SAMPLE=SAMPLE+1))
done < cell_counts.txt
rm cell_counts.txt

