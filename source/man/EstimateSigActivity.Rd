
\name{EstimateSigActivity}
\alias{EstimateSigActivity}
\title{Estimate signature activity}
\description{
   Estimate signature activities from a mutation-count matrix, a
   sample-level panel-context matrix, and a catalog signature-profile matrix.
}
\usage{ EstimateSigActivity(V, L, W, n.start=50, iter.max=5000, eps=1e-5)
}
\arguments{
  \item{V}{Mutation-count matrix or data frame with samples in columns}
  \item{L}{Sample-level panel-context matrix or data frame with samples in columns, see
           \code{\link{GenerateLMatrix}}}
  \item{W}{Catalog signature-profile matrix or data frame with signatures in columns}
  \item{n.start}{Number of initializations. The default is 50.}
  \item{iter.max}{Maximum number of iterations in the EM algorithm. The default is 5000.}
  \item{eps}{Stopping tolerance in the EM algorithm. The default is 1e-5.}
}
\details{ 
The sample-level panel-context matrix \code{L} and mutation-count matrix
\code{V} are of size P
(the mutation context) by N (the sample size). 
The catalog signature-profile matrix \code{W} has dimension P by K
(the number of signatures). For single base substitutions (SBS), 
P is 96. For the objects V, L, and W, we must have dim(V) = dim(L) 
and ncol(W) = K, where K is the number of signatures.
The input matrices must contain numeric, finite and non-negative values, and
\code{V} must contain integer-like mutation counts. When row names and column
names are present, \code{EstimateSigActivity()} checks that \code{V}, \code{L}
and \code{W} contain matched mutation contexts and matched sample IDs. If the
same names are present in a different order, \code{L} and \code{W} are reordered
to match \code{V}. If one object has names on an axis used for alignment and
the paired object does not, or if the named entries differ between matrices, the
function stops with an informative error. Users can run
\code{\link{ValidateSATSInputs}} directly to check and align inputs before
calling this function.
	\code{EstimateSigActivity()} uses an EM algorithm to estimate the
	signature activity matrix. The arguments \code{n.start}, \code{iter.max}
	and \code{eps} control the number of random starts, maximum number of EM
	iterations and stopping tolerance, respectively. Because convergence to a
	local optimum can occur, the default uses multiple initial values
	(\code{n.start = 50}). For each initial value, the default maximum number
	of EM iterations is \code{iter.max = 5000}, and the stopping tolerance is
	\code{eps = 1e-5}. For the catalog
	signature profile matrix \code{W}, reference SBS TMB signature
	profiles in \code{data(SimData)} can be used.
}
\value{ 
  A list containing the estimated activity matrix \code{H}, the
  log-likelihood \code{loglike}, and the logical value \code{converged}.
}
%\references{}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
    data(SimData, package="SATS")
  \donttest{
    EstimateSigActivity(SimData$V, SimData$L, SimData$TrueW_TMB)
  } 
  # For more detailed usage, please refer to README and the user manual
  # in https://github.com/binzhulab/SATS/tree/main.
}
\keyword{ mutational signatures } 
\seealso{ \code{\link{CalculateSignatureBurdens}}, \code{\link{ValidateSATSInputs}} }
