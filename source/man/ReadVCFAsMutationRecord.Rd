\name{ReadVCFAsMutationRecord}
\alias{ReadVCFAsMutationRecord}
\title{Convert a Simple VCF File to a SATS Mutation Record Table}
\description{
Convert a simple single-sample Variant Call Format (VCF) file into the MAF-like mutation
record table used by \code{\link{GenerateVMatrix}}.
}
\usage{
ReadVCFAsMutationRecord(vcf_file, sample_id = NULL,
                        keep_filter = c("PASS", "."))
}
\arguments{
  \item{vcf_file}{Path to a VCF file. Files ending in \code{.gz} are read with
  \code{gzfile}.}
  \item{sample_id}{Optional sample identifier for the returned
  \code{Tumor_Sample_Barcode} column. If \code{NULL}, SATS uses the single VCF
  sample column when exactly one is present, or derives a sample identifier from
  the file name when the VCF has no sample column. Multi-sample VCF genotype
  parsing is not supported by this lightweight converter.}
  \item{keep_filter}{Character vector of FILTER values to retain. Use
  \code{NULL} to retain all records.}
}
\details{
The returned data frame contains the columns required by
\code{GenerateVMatrix}: \code{Chromosome}, \code{Start_Position},
\code{End_Position}, \code{Variant_Type}, \code{Reference_Allele},
\code{Tumor_Seq_Allele2} and \code{Tumor_Sample_Barcode}.

This function is a lightweight input-preparation utility for standard
single-sample targeted-panel VCF files. It does not replace a full clinical
variant-normalization pipeline. Multiallelic records are split by ALT allele,
but complex normalization, multi-sample genotype parsing, tumor-normal
interpretation, phasing and decomposition of complex events should be performed
upstream when needed. \code{GenerateVMatrix} uses
\code{Variant_Type == "SNP"} for SBS analyses and \code{Variant_Type == "DNP"}
for DBS analyses.
}
\value{
A MAF-like data frame suitable for \code{\link{GenerateVMatrix}}.
}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
dir <- system.file("extdata", "refitting_examples", package = "SATS")
vcf_file <- file.path(dir, "SBS_two_variants.vcf")

mutation_record <- ReadVCFAsMutationRecord(vcf_file)
V <- GenerateVMatrix(mutation_record, Class = "SBS", ref.genome = "hg19")
dim(V)
}
\seealso{
\code{\link{GenerateVMatrix}}, \code{\link{ReadBEDAsPanelInfo}}
}
