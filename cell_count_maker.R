args <- commandArgs(trailingOnly=TRUE)

cell_no <- as.numeric(args[1])
output_name <- args[2]


cellCountMaker <- function(prop, total = 40000){
  counts = round(total*prop)
  return(counts)
}

save_path= paste0("Results/", Sys.Date(), "-", output_name)
props <- read.table(paste(save_path, "cell_proportions.txt",   sep="/"), header=TRUE, row.names=1)
cell_counts <- cellCountMaker(props, cell_no)
write.table(cell_counts, file ="cell_counts.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")
