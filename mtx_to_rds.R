library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
dataset = args[1]

### Load expression matrix, cell IDs (barcode.names) and gene names (feature.names)
path = paste0("Datasets/", dataset, "/Deconvolution/raw_matrix")

print("Reading expression matrix:")
exp_mat <- readMM(paste0(path, "/exp_matrix.mtx"))


barcode.names = read.delim(paste0(path, "/barcodes.tsv"), 
                           header = FALSE,
                           stringsAsFactors = FALSE)


feature.names = read.delim(paste0(path, "/gene_names.tsv"), 
                           header = FALSE,
                           stringsAsFactors = FALSE)

### Set column and row names
colnames(exp_mat) <- barcode.names$V1
rownames(exp_mat) <- feature.names$V1
print(paste("Number of genes:", dim(exp_mat)[1]))
print(paste("Number of cells:", dim(exp_mat)[2]))

### Save sparse matrix in proper format
saveRDS(exp_mat, paste0("Datasets/", dataset, "/Deconvolution/", dataset, ".rds"))
print(paste("Saved as", paste0(dataset, ".rds")))