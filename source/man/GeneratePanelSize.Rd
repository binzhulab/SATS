
\name{GeneratePanelSize}
\alias{GeneratePanelSize}
\title{Generate Panel Size Matrix}
\description{
   Generation of the panel size matrix given the panel information.
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
    \item{SBS_order}{Mutation type order either one of "COSMIC" or "signeR" (see details).}
    \item{ref.genome}{The reference genome "hg19" or "hg38".}
}
\details{ 
The first argument 'genomic_information' should contain columns 'Chromosome', 
'Start_Position', 'End_Position', 'SEQ_ASSAY_ID'. 
The column 'Chromosome' contains chromosome number where 'Start_Position' and
'End_Position' columns are start and end positions of the targeted panel. 
The last column 'SEQ_ASSAY_ID' distinguishes different panels consisting of the result. 
Please note that the column names of 'genomic_information' identical to 'Chromosome', 
'Start_Position', 'End_Position', 'SEQ_ASSAY_ID'.
BED target-region files can first be converted into this format with
\code{\link{ReadBEDAsPanelInfo}}.
The second argument specifies mutation type order as either one of "COSMIC" or 
"signeR" where "COSMIC" corresponds to the 96 SBS mutation-type order used by COSMIC v3.2 and 
"signeR" corresponds to the order from the signeR package.
\bold{Note}: The result of 'GeneratePanelSize()' may NOT be an 'L' matrix.
For the standard preprocessing workflow, users can call
\code{\link{GenerateLMatrix}} directly with panel-coordinate and sample-panel
annotation data. \code{GeneratePanelSize()} is retained as a lower-level helper
for users who want to inspect or reuse panel-level context counts.
The resulting augmented matrix can be used as the opportunity matrix for 'signeR()' function, 
'L' matrix for 'EstimateSigActivity()' and 'CalculateSignatureBurdens()' functions. 
Therefore, it is important the mutation type order (row names) should be the same as 
input matrix (Mutation type matrix 'V'). 
We highly recommend confirming that both 'V' and 'L' matrices have the same mutation type 
order. The SBS_order argument controls mutation-type ordering only; the COSMIC
reference-signature version used for mapping is controlled separately by
\code{\link{MappingSignature}}.

}
\value{ 
 A data frame of 96 by 'S' (the number of panels, 'SEQ_ASSAY_ID') where entries 
 denote the number of trinucleotides per million base pairs.
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
