
\name{CalculateSigExpectancy}
\alias{CalculateSigExpectancy}
\title{CalculateSigExpectancy}
\description{
   Estimation of the expected number of mutations attributed 
   by TMB-based catalog signatures (signature expectancy) given 
   the panel size matrix, the catalog signature profile matrix 
   and the signature activities matrix.
}
\usage{ CalculateSigExpectancy(L, W, H)
}
\arguments{
    \item{L}{Panel size matrix or data frame with samples in columns}
    \item{W}{Catalog signature profiles matrix or data frame with signatures in columns}
    \item{H}{Activity matrix or data frame with samples in columns}
}
\details{ 
The panel size matrix L is of size P (the mutation context) by N 
(the sample size). The catalog signature profile matrix has dimension 
of P by K (the number of signatures) and the activity matrix H is of
 size K by N. For single base substitutions (SBS), P is 96. 
If \code{K} is the number of signatures and \code{N} is
  the number of samples, then \code{H} must be of dimension \code{K X N},
  \code{ncol(L) = N}, and \code{ncol(W) = K}.
}
\value{ 
  A matrix of dimension \code{K X N}, where \code{K} is the number
  of signatures and \code{N} is the number of samples.
}
%\references{}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
    data(SimData, package="SATS")

    CalculateSigExpectancy(SimData$L, SimData$TrueW_TMB, SimData$TrueH)
}
\keyword{ mutational signatures } 
\seealso{ \code{\link{EstimateSigActivity}} }
