
\name{GenerateLMatrix}
\alias{GenerateLMatrix}
\title{Generate a sample-level panel-context matrix}
\description{
   Generate the sample-level \code{L} matrix from panel-context or
   panel-coordinate information and sample-panel annotations.
}
\usage{ GenerateLMatrix(Panel_context, Patient_Info, Class = c("SBS", "DBS"),
          SBS_order = c("COSMIC", "signeR"), ref.genome = "hg19")
}
\arguments{
    \item{Panel_context}{Either a panel-context matrix returned from
    \code{\link{GeneratePanelSize}} or a panel-coordinate data frame containing
    \code{Chromosome}, \code{Start_Position}, \code{End_Position} and
    \code{SEQ_ASSAY_ID}.}
    \item{Patient_Info}{A sample-panel annotation data frame containing
    \code{SEQ_ASSAY_ID} and either \code{PATIENT_ID} or \code{SAMPLE_ID}.}
    \item{Class}{A character string specifying the mutation class. Used when
    \code{Panel_context} is a panel-coordinate data frame.}
    \item{SBS_order}{Mutation-type order, either \code{"COSMIC"} or
    \code{"signeR"}. Used for SBS analyses when \code{Panel_context} is a
    panel-coordinate data frame.}
    \item{ref.genome}{The reference genome, either \code{"hg19"} or
    \code{"hg38"}. Used when \code{Panel_context} is a panel-coordinate data
    frame.}
}
\details{\code{GenerateLMatrix()} links each sample to its specific
sequencing assay, identified by \code{SEQ_ASSAY_ID}. For the standard
preprocessing workflow, provide panel-coordinate information as the first
argument and sample-panel annotation as the second argument. The function then
calculates panel-level mutation-context opportunity counts internally and
expands them to a sample-level \code{L} matrix.

For backward compatibility, the first argument can also be a panel-context
matrix returned by \code{\link{GeneratePanelSize}}. In that lower-level workflow,
\code{GenerateLMatrix()} only expands the provided panel-context matrix to the
sample level.

BED target-region files can first be converted into SATS-compatible
panel-coordinate tables with \code{\link{ReadBEDAsPanelInfo}}.
}

\value{ 
 A data frame (L matrix) of P by N, where P is the number of mutation
 channels (96 for SBS or 78 for DBS) and N is the number of tumors.
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

L_mat <- GenerateLMatrix(SimData$PanelEx, PatientInfo, Class="SBS",
                         SBS_order="COSMIC", ref.genome="hg19")

Panel_context <- GeneratePanelSize(genomic_information=SimData$PanelEx,
                                   Class="SBS", SBS_order="COSMIC")
L_mat2 <- GenerateLMatrix(Panel_context, PatientInfo)

stopifnot(identical(colnames(L_mat), colnames(L_mat2)))
stopifnot(identical(rownames(L_mat), rownames(L_mat2)))
# For more detailed usage, please refer to README and the user manual
# in https://github.com/binzhulab/SATS/tree/main.
}
%\keyword{ mutational signatures } 
\seealso{ \code{\link{GeneratePanelSize}}, \code{\link{GenerateVMatrix}},
\code{\link{ReadBEDAsPanelInfo}} }
