context("Input conversion")

dir <- system.file("extdata", "refitting_examples", package = "SATS",
                   mustWork = TRUE)
vcf_file <- file.path(dir, "SBS_two_variants.vcf")
bed_file <- file.path(dir, "SATS_example_panel.bed")

test_that("ReadVCFAsMutationRecord returns SATS mutation-record columns", {
  mut <- SATS::ReadVCFAsMutationRecord(vcf_file)

  expect_equal(nrow(mut), 2)
  expect_true(all(c("Chromosome", "Start_Position", "End_Position",
                    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2",
                    "Tumor_Sample_Barcode") %in% colnames(mut)))
  expect_equal(unique(mut$Tumor_Sample_Barcode), "SATS_VCF_SAMPLE")
  expect_equal(mut$Variant_Type, c("SNP", "SNP"))
})

test_that("ReadVCFAsMutationRecord rejects multi-sample VCF files", {
  multi_vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", "sample_1", "sample_2", sep = "\t"),
    paste("1", "100", ".", "C", "T", ".", "PASS", ".", "GT", "0/1", "0/1",
          sep = "\t")
  ), multi_vcf)

  expect_error(
    SATS::ReadVCFAsMutationRecord(multi_vcf),
    "single-sample VCF",
    fixed = TRUE
  )
})

test_that("ReadBEDAsPanelInfo converts BED coordinates to SATS coordinates", {
  panel <- SATS::ReadBEDAsPanelInfo(
    bed_file = bed_file,
    seq_assay_id = "SATS_EXAMPLE_PANEL",
    name_col = 4
  )

  expect_equal(nrow(panel), 2)
  expect_equal(panel$Start_Position, c(25398271L, 28248241L))
  expect_equal(panel$End_Position, c(25398300L, 28248270L))
  expect_equal(unique(panel$SEQ_ASSAY_ID), "SATS_EXAMPLE_PANEL")
  expect_equal(panel$Hugo_Symbol, c("KRAS", "ALK"))
})

test_that("VCF and BED converters feed GenerateVMatrix and GenerateLMatrix", {
  mut <- SATS::ReadVCFAsMutationRecord(vcf_file)
  panel <- SATS::ReadBEDAsPanelInfo(
    bed_file = bed_file,
    seq_assay_id = "SATS_EXAMPLE_PANEL",
    name_col = 4
  )
  clinical_sample <- data.frame(
    SAMPLE_ID = unique(mut$Tumor_Sample_Barcode),
    SEQ_ASSAY_ID = "SATS_EXAMPLE_PANEL",
    stringsAsFactors = FALSE
  )

  V <- SATS::GenerateVMatrix(mut, Class = "SBS", ref.genome = "hg19")
  L <- SATS::GenerateLMatrix(panel, clinical_sample, Class = "SBS",
                             ref.genome = "hg19")

  expect_equal(dim(V), c(96, 1))
  expect_equal(dim(L), c(96, 1))
  expect_equal(colSums(V), c("SATS_VCF_SAMPLE" = 2))
  expect_true(identical(rownames(V), rownames(L)))
  expect_true(identical(colnames(V), colnames(L)))
})
