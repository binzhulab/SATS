
\name{SATS-package}
\alias{SATS}
\docType{package}
\title{
    SATS (Signature Analyzer for Targeted Sequencing)
}
\description{
  This package performs mutational signature analysis 
  for tumors profiled by targeted sequencing. 
  Unlike the canonical analysis of mutational signatures, 
  SATS factorizes the mutation counts matrix into a panel 
  context matrix (measuring mutation opportunities per million base pairs
  for each tumor), 
  a signature profile matrix, and a signature activity matrix. 
  SATS also calculates the expected number of mutations attributed 
  by a signature, namely signature burden, 
  for each tumor profiled by targeted sequencing.
}
\details{
This package includes a novel algorithm, SATS, 
to perform mutational signature analysis for targeted 
sequencing data. The algorithm first applies the signeR 
algorithm to extract profiles of de novo mutational signatures 
by appropriately adjusting for various panel sizes. 
Next, the profiles of identified de novo mutational signatures 
are mapped to the profiles of catalog signatures of tumor mutation 
burden (TMB), measured as the number of mutations per million base pairs, 
using penalized non-negative least squares. 
Then, given the panel sizes and profiles of mapped TMB catalog signatures, 
signature activities are estimated for all samples simultaneously 
through the Expectation-Maximization (EM) algorithm. 
Finally, the expected number of mutations attributed by a signature, 
namely signature burden, is calculated for each tumor profiled by targeted sequencing.

The main functions in this package are \code{\link{GenerateVMatrix}},
\code{\link{GenerateLMatrix}}, \code{\link{ValidateSATSInputs}},\cr
\code{\link{EstimateSigActivity}}, \code{\link{CalculateSignatureBurdens}},
and \code{\link{MappingSignature}}.

The accompanying manuscript applies SATS to 111,711 tumors profiled by targeted sequencing
from AACR Project GENIE to construct a real-world, panel-calibrated pan-cancer
catalogue of mutational signatures.
}
\references{
Lee, D., Hua, M., Wang, D., Song, L., Yu, K., Yang, X., Shi, J., Landi, M., 
Zhu, B. A real-world pan-cancer catalogue of mutational signatures from
111,711 tumors. Submitted.
}
\author{Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>}
\keyword{package}
