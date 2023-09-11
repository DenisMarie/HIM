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



#opt_parser = optparse::OptionParser(option_list=option_list);
#opt = optparse::parse_args(opt_parser);

# Specify verbosity
verbose = isTRUE(opt$verbose)

# Check number of threads specified
ncpu <- parallel::detectCores()-1
if(verbose) message("Threads detected: ", ncpu, ". Threads specified: ", opt$threads, ".")
ncpu <- min(ncpu, opt$threads)
if(verbose) message("Using ", ncpu, " thread(s).")

### START ###
# Check that function and data files exist in wd
stopifnot("ERROR: Path to functions R script does not exist." = file.exists(opt$path_to_functions))
source(opt$path_to_functions)


# Simulated datasets
n<- 60; p <- 100;  q <-50
X1 <- matrix(abs(rnorm(n *p , mean =20, sd = 5)), ncol = p)
colnames(X1) <- paste0("G",1:p)
X2 <- matrix(abs(rnorm(n *q , mean =10, sd = 4)), ncol = q)
colnames(X2)  <- paste0("miRNA",1:q)
score_mat <- matrix(rbinom(p * q , size = 1, prob = 0.2), nrow = p, ncol = q)
Y <- rep(c(0,1), each = n/2)


system.time({
  him <- HIM(Y = Y,
             X1=X1,
             X2=X2,
             S=score_mat,
             niter=100, burn=10, burnin=10, nmc=100, thres=0.2, ggm = F,
             n_threads = ncpu)
})

if(verbose) message("Saving output...")
save(list="him", file = paste0("Result_",Sys.Date(),".RData"))
if(verbose) message("Done!")

