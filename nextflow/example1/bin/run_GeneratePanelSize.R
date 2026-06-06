#!/usr/bin/env Rscript

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
obj <- c("genomic_information", "Class", "SBS_order", "ref.genome")
if (!(obj[1] %in% tmp)) stop("ERROR: input file must contain the object 'genomic_information'")
if (!(obj[2] %in% tmp)) Class      <- "SBS"
if (!(obj[3] %in% tmp)) SBS_order  <- "COSMIC"
if (!(obj[4] %in% tmp)) ref.genome <- "hg19"

# Use an installed SATS package when available. Otherwise, install the current
# GitHub source into a local library in the Nextflow work directory.
local_lib <- file.path(getwd(), "r_libs")
dir.create(local_lib, showWarnings=FALSE, recursive=TRUE)
.libPaths(c(local_lib, .libPaths()))

if (!requireNamespace("SATS", quietly=TRUE)) {
  if (!requireNamespace("remotes", quietly=TRUE)) {
    install.packages("remotes", repos="https://cloud.r-project.org", lib=local_lib)
  }
  remotes::install_github("binzhulab/SATS", subdir="source", upgrade="never", lib=local_lib)
}

# Load SATS and run the function
library(SATS)

ret <- SATS::GeneratePanelSize(genomic_information, Class=Class, SBS_order=SBS_order,
                               ref.genome=ref.genome)
save(ret, file=outfile)


