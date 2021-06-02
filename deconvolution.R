source("helper_functions.R")
args <- commandArgs(trailingOnly=TRUE)
#args <- c("whatevs", "PBMC", "test")


input_data <- args[1]
dataset <- args[2]
output_name <- args[3]


save_path= paste0("Results/", Sys.Date(), "-", output_name)
dir.create(file.path(save_path))

### Read reference matrix and markers for deconvolution

reference_path = paste("Datasets", dataset, "Deconvolution/reference/", sep="/")
C = read.table(paste0(reference_path, "C"))
markers = readRDS(paste0(reference_path, "markers.rds"))


### Read input bulk mRNA_seq data (.txt for now)
T <- read.delim(input_data)
#T= read.delim(paste0(reference_path, "mixes.txt"))#[, 1]
#T1 = readRDS(paste0(reference_path, "mixes.rds"))
#print(typeof(T))
refProfiles.var = read.table(paste0(reference_path, "Profiles.var"))


to_remove = "none"

# marker selection (on training data)
marker_distrib = marker_strategies(markers, "all", C)
  
RESULTS = Deconvolution(T = T, C = C, method = "nnls", elem = to_remove, marker_distrib = marker_distrib, refProfiles.var = refProfiles.var) 

RESULTS_wide = reshape(RESULTS, idvar="tissue", timevar="CT", direction="wide")
write.table(RESULTS_wide, file = paste(save_path,"cell_proportions.txt",sep="/"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")