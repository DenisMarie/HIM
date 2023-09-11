#!/usr/bin/env Rscript --vanilla
# Author: Ressom Lab
# Description: Run the HIMr package, an implementation of the HIM method created by Marie Denis
# Parameters:
## --path_to_data - path to where the relevant files are located
## --ncpu - number of cores/threads on which to run (can only set >1 for *nix machines, e.g. Mac, Linux)
##


# Install and load packages
package_list <- c(
  "optparse",
  "parallel",
  "tictoc",
  "pbapply",
  "Matrix",
  "huge",
  "igraph",
  "coda",
  "glmnet",
  "bayesreg",
  "mvtnorm",
  "ggplot2",
  "glmnetUtils",
  "BoomSpikeSlab",
  "truncnorm",
  "extraDistr",
  "igraph"
)

for(package in package_list){
  if(!require(package, character.only = T, quietly = T)){
    install.packages(package, quiet = T, ask = F)
    require(package, character.only = T)
  }
}

# Parse and verify params:
opt <- list()
opt$path_to_data <- "data/"
opt$path_to_functions <- "R/functions_final.R"
opt$threads = 4
opt$delim = ";"
opt$verbose = T
opt$subset_size = 200

# option_list = list(
#   optparse::make_option(c("-p", "--path_to_data"), type="character", default="data/",
#               help="Path to directory where relevant files (CSV files for miRNA, mRNA, and Score_mat) are located. Defaults to data/", metavar="character"),
#   optparse::make_option(c("-t", "--threads"), type="integer", default=1,
#               help="output file name [default= %default]", metavar="integer"),
#   optparse::make_option(c("-d", "--delim"), type = "character", default = ";",
#                         help = "Single character delimiter specifying spearator of all CSV files. For legacy reasons, defaults to semicolon (';')"),
#   optparse::make_option(c("-v", "--verbose"), type = "logical", action = "store_true",
#                         help = "Output messages and progress meters.")
#
# );

#opt_parser = optparse::OptionParser(option_list=option_list);
#opt = optparse::parse_args(opt_parser);

# Specify verbosity
verbose = isTRUE(opt$verbose)

# Check directory and contents
if(verbose) message("Working directory: ", getwd())
stopifnot("ERROR: Path to data directory is not correctly specified, or does not exist." = dir.exists(opt$path_to_data))
wd <- opt$path_to_data
miRNA.csv <- dir(wd, pattern = glob2rx("miRNA*csv"), full.names = T)
mRNA.csv <- dir(wd, pattern = glob2rx("mRNA*csv"), full.names = T)
Score_mat.csv <- dir(wd, pattern = glob2rx("Score_mat*csv"), full.names = T)
sample_labels.csv <- dir(wd, pattern = glob2rx("sample*csv"), full.names = T)
stopifnot("ERROR: No CSV with miRNA information found. Ensure that miRNA CSV file contains 'miRNA' in the file name." = length(miRNA.csv)>0)
stopifnot("ERROR: No CSV with mRNA information found. Ensure that mRNA CSV file contains 'mRNA' in the file name." = any(grepl(glob2rx("mRNA*csv"), ignore.case = T, dir(wd))))
stopifnot("ERROR: No CSV with Score_mat information found. Ensure that Score_mat CSV file contains 'Score_mat' in the file name." = any(grepl(glob2rx("Score_mat*csv"), ignore.case = T, dir(wd))))
stopifnot("ERROR: No CSV with sample information found. Ensure that sample CSV file contains 'sample' in the file name." = any(grepl(glob2rx("sample*csv"), ignore.case = T, dir(wd))))


# Check if delimiter is correct
delim = opt$delim
if(!any(grepl(delim, readLines(con = mRNA.csv, n = 1)))){
  stop(paste0("Delimiter [",delim,"] is not found in first line of mRNA CSV. Please specify correct delimiter with --delim."))
}

# Check number of threads specified
ncpu <- parallel::detectCores()-1
if(verbose) message("Threads detected: ", ncpu, ". Threads specified: ", opt$threads, ".")
ncpu <- min(ncpu, opt$threads)
if(verbose) message("Using ", ncpu, " thread(s).")

### START ###
# Check that function and data files exist in wd
stopifnot("ERROR: Path to functions R script does not exist." = file.exists(opt$path_to_functions))
source(opt$path_to_functions)



# Read CSV files
if(verbose) message("Reading in data files...")
nor_miRNA <- read.csv("data/miRNA_Enzy.csv", sep = opt$delim)
nor_mRNA <- read.csv("data/mRNA_Enzy.csv", sep = opt$delim)
score_mat <- read.csv("data/Score_mat_Enzy.csv", sep = opt$delim)
sample_labels <- as.numeric(read.csv("data/sample_labels.csv", sep = opt$delim)[,1])
XY.spvas <- read.csv("data/df_xy_miRNA_full.csv", sep = opt$delim)
# Initialize score matrix:
if(verbose) message("Cleaning score matrix...")
cor_mat <- cor(nor_mRNA, nor_miRNA)
score_mat[abs(cor_mat) < 0.4] <- 0

# for (r in 1: nrow(score_mat)){
#   pbapply::setpb(pb = pb, value = r)
#   for (c in 1:ncol(score_mat)){
#     if (score_mat[r,c] !=0){
#       if (
#         abs(cor(nor_mRNA[,r], nor_miRNA[,c])) < 0.4
#       )
#         score_mat[r,c] <- 0
#     }
#   }
# }
# pbapply::closepb(pb)


if(opt$subset_size > 5){
  subset <- 1:opt$subset_size
}else{
  if(verbose) message("Subset size not specified or too small. Using all available data.")
  subset <- 1:ncol(nor_mRNA)
}
system.time({
  him <- HIM(Y = sample_labels,
             X1=as.matrix(nor_mRNA),
             X2=as.matrix(nor_miRNA),
             S=score_mat,
             niter=10000, burn=1000, burnin=100, nmc=1000, thres=0.2, ggm = F,
             n_threads = ncpu)#, quiet = !verbose)
})


if(verbose) message("Saving output...")
save(list="him", file = paste0("Result_",Sys.Date(),".RData"))
if(verbose) message("Done!")
# str(him)

# represent the adjusted graph
# myg.mRNA <- graph_from_adjacency_matrix(as.matrix(him$G), mode="undirected")
# gsize(myg.mRNA)
# deg <- degree(myg.mRNA, mode="all")
# V(myg.mRNA)$size <- 5
# V(myg.mRNA)$color <- as.factor(info_mRNA[,2])
# V(myg.mRNA)$label <- colnames(nor_mRNA)
# V(myg.mRNA)$label.cex <- 0.8
# lol <- cluster_fast_greedy(myg.mRNA)
# plot.igraph(myg.mRNA,layout=layout.fruchterman.reingold)

# selecting the variables
# dim(him$gg)
# gg.hat <- colMeans(him$gg)
# plot(gg.hat, main = '106 Genes')
# thres <- 0.2
# colnames(him$gg)[gg.hat > thres]









