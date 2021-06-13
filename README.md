# Simulation of realistic single-cell mRNA-seq lookalike data from bulk samples
Code implementations of the thesis work of Viktoria Lavro 

Folder structure:
```
.
├── Raw_data
│   ├── Baron
│   └── PBMC
├── Datasets
│   ├── Baron
|   |    ├── scGAN
|   |    └── Deconvolution
|   |           └──references
│   └── PBMC
|   |    ├── scGAN
|   |    └── Deconvolution
|   |           └──references
├── scGAN
│   ├── preprocessing
|   ...
│   └── main.py
├── cell_count_maker.R
├── deconvolution.R
├── generate.sh
├── h5ad_to_mtx.py
├── helper_functions.R
├── make_reference.R
├── mtx_to_rds.R
└── result_plotter_cells.py


## Usage
### Training
Instructions for training are identical to the original implementation of cscGAN:
See: https://github.com/imsb-uke/scGAN


