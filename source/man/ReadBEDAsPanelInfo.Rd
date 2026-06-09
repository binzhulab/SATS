\name{ReadBEDAsPanelInfo}
\alias{ReadBEDAsPanelInfo}
\title{Convert BED Target Regions to SATS Panel Information}
\description{
Convert a Browser Extensible Data (BED) target-region file into the
SATS-compatible panel-coordinate table used by \code{\link{GeneratePanelSize}}
and \code{\link{GenerateLMatrix}}.
}
\usage{
ReadBEDAsPanelInfo(bed_file, seq_assay_id = NULL, has_header = FALSE,
                   chromosome_col = 1, start_col = 2, end_col = 3,
                   seq_assay_col = NULL, name_col = NULL)
}
\arguments{
  \item{bed_file}{Path to a BED file. Files ending in \code{.gz} are read with
  \code{gzfile}.}
  \item{seq_assay_id}{Sequencing assay identifier to assign to every interval
  when \code{seq_assay_col} is not supplied.}
  \item{has_header}{Logical value indicating whether the BED file has a header
  row.}
  \item{chromosome_col}{Column index or name for chromosome.}
  \item{start_col}{Column index or name for the BED start coordinate.}
  \item{end_col}{Column index or name for the BED end coordinate.}
  \item{seq_assay_col}{Optional column index or name containing
  \code{SEQ_ASSAY_ID}.}
  \item{name_col}{Optional column index or name to preserve as
  \code{Hugo_Symbol}.}
}
\details{
BED coordinates are 0-based and half-open. SATS panel-coordinate tables use
1-based inclusive coordinates. Therefore this function returns
\code{Start_Position = BED_start + 1} and \code{End_Position = BED_end}.

The returned data frame contains at least \code{Chromosome},
\code{Start_Position}, \code{End_Position} and \code{SEQ_ASSAY_ID}; these
columns can be passed directly to \code{GeneratePanelSize} or
\code{GenerateLMatrix}.
}
\value{
A data frame of panel target intervals in SATS coordinate format.
}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
dir <- system.file("extdata", "refitting_examples", package = "SATS")
bed_file <- file.path(dir, "SATS_example_panel.bed")

panel_info <- ReadBEDAsPanelInfo(
    bed_file = bed_file,
    seq_assay_id = "SATS_EXAMPLE_PANEL",
    name_col = 4
)
head(panel_info)
}
\seealso{
\code{\link{GeneratePanelSize}}, \code{\link{GenerateLMatrix}},
\code{\link{ReadVCFAsMutationRecord}}
}
