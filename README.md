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
Folders to create in advance:
 - Datasets
 - Dataset folders within Datasets
 - Deconvolution folders within each dataset folder
 - Results

## Usage
### Training
Instructions for training are identical to the original implementation of cscGAN:
See: https://github.com/imsb-uke/scGAN

### Reference matrices for deconvolution
When the training is done, the resulting `.h5ad` file is used to create the reference gene expression profile matrix C. First the `.h5ad` file must be converted to `.RDS` by first using the command
```
python h5ad_to_mtx.py --dataset "Dataset"
```
where "Dataset" is the corresponding folder name within the Datasets folder (e.g. PBMC, Baron, etc.). The same guidelines are used for the next step in
```
Rscript mtx_to_rds.R "Dataset"
```
Finally, all reference files are created using the next command. In case you want to create bulk pseudomixtures, you may use an integer as a second argument, to specify the number of mixtures to generate. These will be created from 1000 cells each.
```
Rscript make_reference.R "Dataset" 20
```

### Generation of simulated cells
To simulate single-cell lookalike data, you can use the `generate.sh` script with the options `-i` as the location of the bulk sample file, `-d` as the name of the dataset (same way as in previous section), `-n` as the number of cells to simulate per sample and `-o` as the name to give to the folder, where the results should be saved (the current date will be automatically included in it)
```
bash generate.sh -i /path/to/bulk/samples/file -d "Dataset" -n 1000 -o Chosen_folder_name
```
You will find the generate `.h5ad` files at `.../Results/yyy-mm-dd-Chosen_folder_name`

### Plotting cells simulated from pseudomixtures
If you used the bulk pseudomixtures, you can use `result_plotter_cells.py` to plot UMAP and barplots collected in a PDF and heatmaps as images.

![image](https://user-images.githubusercontent.com/54985154/121823963-f0b68e80-cca8-11eb-85c2-26ec69d57827.png)
![image](https://user-images.githubusercontent.com/54985154/121823971-04fa8b80-cca9-11eb-8276-71b0489b2d73.png)



