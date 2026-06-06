if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github("binzhulab/SATS", subdir = "source", upgrade = "never")
