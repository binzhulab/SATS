#!usr/bin/env/ Rscript

# Get input and output files
args    <- commandArgs(trailingOnly=TRUE)

# Check args
nargs <- length(args)
if (nargs < 2) {
  stop("ERROR: an input file (arg[1]) and output file (arg[2]) must be specified")
}

infile  <- args[1]
outfile <- args[2]

# Check that input file exists and can be loaded
if (!file.exists(infile)) {
  msg <- paste0("ERROR file not found: ", infile)
  stop(msg)
}
tmp <- try(load(infile), silent=FALSE)
if ("try-error" %in% class(tmp)) {
  stop("ERROR: input file (arg[1]) must be an .rda file to be loaded with the load() function")
}
obj <- c("genomic_information", "Types")
if (!(obj[1] %in% tmp)) stop("ERROR: input file must contain the object 'genomic_information'")
if (!(obj[2] %in% tmp)) stop("ERROR: input file must contain the object 'Types'")

# Install devtools package 
if (!require("devtools")) install.packages("devtools", repos="https://cloud.r-project.org")
devtools::install_github("binzhulab/SATS", subdir="source", upgrade="never", lib="./")

.libPaths(c("./", .libPaths()))

# Load SATS and run the function
library(SATS)

ret <- SATS::GeneratePanelSize(genomic_information, Types)
save(ret, file=outfile)
