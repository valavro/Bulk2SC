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
|        ├── scGAN
|        └── Deconvolution
|              └──references
├── scGAN
│   ├── preprocessing
|   ...
│   └── main.py
├── Results
│   ├── 2021-03-14-Result_1
|   ...
│   └── 2021-06-12-Result_n
├── cell_count_maker.R
├── deconvolution.R
├── generate.sh
├── h5ad_to_mtx.py
├── helper_functions.R
├── make_reference.R
├── mtx_to_rds.R
└── result_plotter_cells.py
```

## Usage
### Training
Instructions for training are identical to the original implementation of cscGAN:
See: https://github.com/imsb-uke/scGAN

### Reference matrices for deconvolution
When the training is done, the resulting `.h5ad` file is used to create the reference gene expression profile matrix C. First the `.h5ad` file must be converted to `.RDS` by first using the command
'''
python 
'''




