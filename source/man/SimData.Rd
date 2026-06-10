% --- Source file: SimData.Rd ---
\name{SimData}
\alias{SimData}

\title{
    Simulated example data
}
\description{
  Simulated data used in SATS examples and tests.
}
\details{
  This file consists of the list \code{SimData} with the following objects: 
  \itemize{
    \item{V} : Simulated mutation-count matrix of size 96 by 10027
	            generated from a Poisson distribution with mean
				corresponding to each element of L*(WH), where * denotes
				elementwise multiplication, L is the sample-level panel-context
				matrix below, W is the tumor mutation burden-normalized
				signature-profile matrix (TrueW_TMB below), and H is the
				signature-activity matrix (TrueH below).
	
    \item{L} : Sample-level panel-context matrix of size 96 by 10027.
    \item{TrueH} : Simulated signature-activity matrix of size 6 by 10027 used to generate V.
    \item{TrueW_TMB} : Tumor mutation burden-normalized signature-profile matrix of size 96 by 6
                        for 6 SBS signatures (SBS1, SBS2/13, SBS4, SBS5, SBS40, SBS89) 
						used to generate V.
    \item{SingleTumorEx} : A single simulated tumor (column \code{singleV}) 
                    and its sequencing context (column \code{singleL}).
    \item{PanelEx}: An example sequencing panel information with Chromosome, Start_Position,
                    End_Position, SEQ_ASSAY_ID.
    \item{PatientInfo}: An example patient data set containing patient IDs associated with
                        SEQ_ASSAY_ID.
  }
 }
%\examples{}
\keyword{data}
