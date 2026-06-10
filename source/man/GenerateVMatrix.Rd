\name{GenerateVMatrix}
\alias{GenerateVMatrix}
\title{Generate a Mutation Count Matrix from Mutation Records}
\description{
Generate an SBS96 or DBS78 mutation count matrix from a MAF-like mutation
record table.
}
\usage{
GenerateVMatrix(mutation_record, Class = c("SBS", "DBS"),
                ref.genome = "hg19", mutation_order = NULL)
}
\arguments{
  \item{mutation_record}{A MAF-like data frame containing mutation records.}
  \item{Class}{A character string specifying the mutation class. Must be either
  \code{"SBS"} for single base substitutions or \code{"DBS"} for double base
  substitutions.}
  \item{ref.genome}{The reference genome, either \code{"hg19"} or \code{"hg38"}.}
  \item{mutation_order}{Optional character vector specifying the row order of
  the returned mutation count matrix. If \code{NULL}, SATS uses the row order of
  \code{RefTMB$TMB_SBS_v3.4} for SBS or \code{RefTMB$TMB_DBS_v3.4} for DBS.}
}
\details{
The \code{mutation_record} data frame must contain the columns
\code{Chromosome}, \code{Start_Position}, \code{End_Position},
\code{Variant_Type}, \code{Reference_Allele}, \code{Tumor_Seq_Allele2} and
\code{Tumor_Sample_Barcode}. SBS analyses use records with
\code{Variant_Type == "SNP"} and generate COSMIC-style 96-channel mutation
contexts. DBS analyses use records with \code{Variant_Type == "DNP"} and
generate 78 double-base substitution contexts.

This function is intended for MAF-like tabular mutation records. Simple
single-sample VCF files can first be converted into this tabular format with
\code{\link{ReadVCFAsMutationRecord}}.
}
\value{
A numeric matrix with mutation contexts in rows and tumor samples in columns.
}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
dir <- system.file("extdata", "refitting_examples", package = "SATS")
sbs_file <- file.path(dir, "SBS_MAF_two_samples.txt")
sbs_mut <- read.table(sbs_file, header = TRUE, sep = "\t", quote = "",
                      stringsAsFactors = FALSE)

V_sbs <- GenerateVMatrix(sbs_mut, Class = "SBS", ref.genome = "hg19")
dim(V_sbs)
}
\seealso{
\code{\link{GenerateLMatrix}}, \code{\link{GeneratePanelSize}},
\code{\link{ReadVCFAsMutationRecord}}
}
