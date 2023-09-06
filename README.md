# SATS (Signature Analyzer for Targeted Sequencing)
<br/>

### Introduction
This package is created to perform mutational signature 
  analysis for targeted sequenced tumors. 
  Unlike the canonical analysis of mutational signatures, 
  SATS factorizes the mutation counts matrix into a panel 
  context matrix (measuring the size of the targeted sequenced genome 
  for each tumor in the unit of million base pairs (Mb)), 
  a signature profile matrix, and a signature activity matrix. 
  SATS also calculates the expected number of mutations attributed 
  by a signature, namely signature expectancy, 
  for each targeted sequenced tumor.

For more information please refer to the [user guide](https://github.com/binzhulab/SATS/blob/main/User_Guide_SATS.pdf).
<br/>

### Installation
To install from Github, use the devtools R package:
```r
if (!requireNamespace("devtools", quietly = TRUE))  
	install.packages("devtools")
devtools::install_github("binzhulab/SATS/source")
```
Alternatively, download the package and follow the steps below. Download SATS_0.0.7.tar.gz (for Unix) or SATS_0.0.7.zip (for Windows, R version >= 4.1). To install SATS on Unix, enter the command from a Unix prompt:
```
R CMD INSTALL SATS_0.0.7.tar.gz -l path_to_install_package
```
Alternatively, SATS_0.0.7.tar.gz (for Unix) or SATS_0.0.7.zip (for Windows, R version >= 4.1) from the [Github page](https://github.com/binzhulab/SATS) are available and one may use the following commands:
```
install.packages("./SATS_0.0.7.tar.gz", repos = NULL, type = "source")
install.packages("./SATS_0.0.7.zip", repos = NULL, type = "win.binary")
```
Once the installation is successful, it can be loaded in **R** by calling 
```
library(SATS)
```
<br/>
