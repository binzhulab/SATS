\name{ValidateSATSInputs}
\alias{ValidateSATSInputs}
\title{Validate SATS input matrices}
\description{
  Validate and align SATS mutation-count, panel-context, signature-profile and
  signature-activity matrices before running refitting or burden calculation.
}
\usage{
ValidateSATSInputs(V=NULL, L=NULL, W=NULL, H=NULL)
}
\arguments{
  \item{V}{Optional mutation-count matrix or data frame with mutation contexts
           in rows and samples in columns. Values must be non-negative,
           finite and integer-like.}
  \item{L}{Optional panel-context matrix or data frame with mutation contexts
           in rows and samples in columns. Values must be non-negative and
           finite.}
  \item{W}{Optional catalog signature-profile matrix or data frame with
           mutation contexts in rows and signatures in columns. Values must be
           non-negative and finite.}
  \item{H}{Optional signature-activity matrix or data frame with signatures in
           rows and samples in columns. Values must be non-negative and finite.}
}
\details{
\code{ValidateSATSInputs()} performs the same matrix-level checks used
internally by \code{\link{EstimateSigActivity}} and
\code{\link{CalculateSignatureBurdens}}. It checks object type, matrix
dimensions, numeric values, finite values, non-negative values, duplicated or
missing row and column identifiers, mutation-context alignment, sample-ID
alignment and signature-name alignment.

When compatible names are present in a different order, the returned matrices
are reordered to a common order. For example, \code{L} and \code{W} are
reordered to match \code{V} for activity estimation, and \code{L} and \code{W}
are reordered to match \code{H} for burden calculation. If one object has names
on an axis used for alignment and the paired object does not, or if the named
entries differ, the function stops with an informative error. If paired axes
are both unnamed, validation proceeds by dimension.
}
\value{
  A named list containing the supplied matrices after compatible named axes have
  been aligned. Objects that were not supplied are omitted from the returned
  list.
}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(SimData, package="SATS")

  obj <- ValidateSATSInputs(
      V = SimData$V[, 1:3],
      L = SimData$L[, c(3, 1, 2)],
      W = SimData$TrueW_TMB
  )
  identical(colnames(obj$L), colnames(obj$V))
}
\keyword{ mutational signatures }
\seealso{
  \code{\link{EstimateSigActivity}},
  \code{\link{CalculateSignatureBurdens}}
}
