
\name{GeneratePanelSize}
\alias{GeneratePanelSize}
\title{Generate a panel-context matrix}
\description{
   Generate a panel-context matrix from panel-coordinate information.
}
\usage{ GeneratePanelSize(genomic_information, Class = c("SBS", "DBS"), 
          SBS_order = c("COSMIC", "signeR"), ref.genome="hg19")
}
\arguments{
    \item{genomic_information}{Data frame of panel information (see details).}
    \item{Class}{A character string specifying the mutation class. Must be either
                 "SBS" (for single base substitutions) or "DBS" (for double base
                 substitutions). If "DBS" is selected, the \code{SBS_order} 
                 argument is ignored. }
    \item{SBS_order}{Mutation-type order, either \code{"COSMIC"} or
    \code{"signeR"} (see details).}
    \item{ref.genome}{The reference genome, either \code{"hg19"} or
    \code{"hg38"}.}
}
\details{ 
The first argument \code{genomic_information} must contain the columns
\code{Chromosome}, \code{Start_Position}, \code{End_Position} and
\code{SEQ_ASSAY_ID}. The \code{Chromosome} column identifies the chromosome,
\code{Start_Position} and \code{End_Position} define the targeted interval, and
\code{SEQ_ASSAY_ID} identifies the sequencing assay or panel. These column names
must be provided exactly.
BED target-region files can first be converted into this format with
\code{\link{ReadBEDAsPanelInfo}}.

For SBS analyses, \code{SBS_order = "COSMIC"} returns the COSMIC-style SBS96
mutation-channel order and \code{SBS_order = "signeR"} returns the order used by
the \code{signeR} package. This argument controls the row order of the returned
mutation-opportunity matrix only; the COSMIC reference-signature version used for
mapping is controlled separately by \code{\link{MappingSignature}}.

\bold{Note}: The result of \code{GeneratePanelSize()} may not be an \code{L}
matrix.
For the standard preprocessing workflow, users can call
\code{\link{GenerateLMatrix}} directly with panel-coordinate and sample-panel
annotation data. \code{GeneratePanelSize()} is retained as a lower-level helper
for users who want to inspect or reuse panel-level context counts.
The resulting panel-level matrix can be used as the opportunity matrix for
\code{signeR()} or expanded to a sample-level \code{L} matrix with
\code{\link{GenerateLMatrix}}. The mutation-channel row order must match the row
order of the mutation-count matrix \code{V}.

}
\value{
 A data frame of P by S, where P is the number of mutation channels
 (96 for SBS or 78 for DBS) and S is the number of sequencing panels
 (\code{SEQ_ASSAY_ID}). Entries denote mutation opportunities per million
 base pairs for each panel.
}
%\references{}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(SimData, package="SATS")
  keep <- match(SimData$PatientInfo$SEQ_ASSAY_ID,
                unique(SimData$PanelEx$SEQ_ASSAY_ID), nomatch=0) > 0
  PatientInfo <- SimData$PatientInfo[keep, ]

  ## For the single base substitution (SBS) case:
  Panel_context1 <- GeneratePanelSize(genomic_information = SimData$PanelEx, Class = "SBS",  
                                      SBS_order = "COSMIC", ref.genome="hg19")
  L_mat1 <- GenerateLMatrix(Panel_context1, PatientInfo)
  Panel_context2 <- GeneratePanelSize(genomic_information = SimData$PanelEx, Class = "SBS",  
                                      SBS_order = "signeR", ref.genome="hg19")
  L_mat2 <- GenerateLMatrix(Panel_context2, PatientInfo)

  ## For the double base substitution (DBS) case:
  Panel_context_DBS <- GeneratePanelSize(genomic_information = SimData$PanelEx, Class = "DBS", 
                                         ref.genome = "hg19")
  L_mat_DBS  <- GenerateLMatrix(Panel_context_DBS, PatientInfo)

  # For more detailed usage, please refer to README and the user manual
  # in https://github.com/binzhulab/SATS/tree/main 
}
%\keyword{ mutational signatures } 
\seealso{ \code{\link{GenerateLMatrix}}, \code{\link{GenerateVMatrix}},
\code{\link{ReadBEDAsPanelInfo}}, \code{\link{EstimateSigActivity}} }
