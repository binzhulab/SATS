
\name{MappingSignature}
\alias{MappingSignature}
\title{Find a subset of TMB-based catalog SBS signatures}
\description{
   This function finds a subset of TMB-based catalog SBS signatures 
   whose linear combinations approximate
   \code{de novo} SBS signatures detected by signeR.
}
\usage{ MappingSignature(W_hat, W_ref=NULL, niter=100, cutoff.I2=0.1, min.repeats=80, 
                         COSMICv="v3.4")
}
\arguments{
    \item{W_hat}{Matrix or data frame of de novo signatures from signeR}
    \item{W_ref}{NULL or a matrix or data frame of TMB-based catalog signatures.
                 If NULL, then it will default to \code{RefTMB$TMB_SBS_v3.2} or 
                 \code{RefTMB$TMB_SBS_v3.4}
                 depending on the value of \code{COSMICv} below
                 (see \code{\link{RefTMB}}).}
    \item{niter}{Number of iterations. The default is 100.}
    \item{cutoff.I2}{Coefficient cutoff used to select reference signatures.
                     The default is 0.1, requiring a retained reference
                     signature to contribute more than 10 percent to a
                     \code{de novo} profile in a repeated fit.}
    \item{min.repeats}{Minimum number of repeated fits in which a reference
                       signature must satisfy \code{coeff > cutoff.I2}. The
                       default is 80.}
    \item{COSMICv}{Version of the TMB-based COSMIC signatures ("v3.2" or "v3.4"). 
                   This option is ignored if \code{W_ref} is not \code{NULL}. 
                   The default is "v3.4".}
}
\details{ 
MappingSignature() applies penalized non-negative least squares 
(pNNLS) for selecting the TMB-based catalog signatures. 
Specifically, it repeats pNNLS 100 times (niter) to reduce 
the randomness of cross-validation involved in pNNLS. 
Then TMB-based catalog signatures are selected with a coefficient
greater than \code{cutoff.I2} in at least \code{min.repeats} repeats.
With the default settings, a retained signature must have a coefficient
greater than 0.1 in at least 80 of 100 repeated fits.
}
\value{ 
A data frame with one row per retained reference signature. The
\code{Reference} column contains the corresponding column names of \code{W_ref}
(COSMIC SBS names when the default COSMIC reference matrix is used), and
\code{freq} gives the number of repeated fits in which the signature coefficient
exceeded \code{cutoff.I2}.
}
%\references{}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
%\examples{}
%\keyword{ mutational signatures } 
%\seealso{ \code{\link{} }
