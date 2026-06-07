context("GenerateVMatrix")

dir <- system.file("extdata", package = "SATS", mustWork = TRUE)
exdir <- file.path(dir, "refitting_examples")

sbs_file <- file.path(exdir, "SBS_MAF_two_samples.txt")
dbs_file <- file.path(exdir, "DBS_MAF_two_samples.txt")

sbs_mut <- read.table(sbs_file, header = TRUE, sep = "\t", quote = "",
                      stringsAsFactors = FALSE)
dbs_mut <- read.table(dbs_file, header = TRUE, sep = "\t", quote = "",
                      stringsAsFactors = FALSE)

data(RefTMB, package = "SATS")
data(SimData, package = "SATS")

test_that("GenerateVMatrix returns COSMIC-ordered SBS counts", {
  V <- SATS::GenerateVMatrix(sbs_mut, Class = "SBS", ref.genome = "hg19")

  expect_equal(dim(V), c(96, 2))
  expect_true(identical(rownames(V), rownames(RefTMB$TMB_SBS_v3.4)))
  expect_equal(colSums(V), c("GENIE-DFCI-050984-218969" = 2,
                             "GENIE-DFCI-109295-436549" = 5))
})

test_that("GenerateVMatrix returns DBS78-ordered DBS counts", {
  V <- SATS::GenerateVMatrix(dbs_mut, Class = "DBS", ref.genome = "hg19")

  expect_equal(dim(V), c(78, 2))
  expect_true(identical(rownames(V), rownames(RefTMB$TMB_DBS_v3.4)))
  expect_equal(colSums(V), c("GENIE-DFCI-000735-437254" = 2,
                             "GENIE-MSK-P-0071479-T02-IH4" = 2))
})

test_that("GenerateVMatrix and GenerateLMatrix align SBS V and L matrices", {
  clinical_sample <- data.frame(
    SAMPLE_ID = unique(sbs_mut$Tumor_Sample_Barcode),
    SEQ_ASSAY_ID = SimData$PatientInfo$SEQ_ASSAY_ID[1],
    stringsAsFactors = FALSE
  )

  V <- SATS::GenerateVMatrix(sbs_mut, Class = "SBS", ref.genome = "hg19")
  L <- SATS::GenerateLMatrix(
    Panel_context = SimData$PanelEx,
    Patient_Info = clinical_sample,
    Class = "SBS",
    ref.genome = "hg19"
  )

  expect_true(identical(rownames(V), rownames(L)))
  expect_true(identical(colnames(V), colnames(L)))
  expect_equal(dim(V), dim(L))
})

test_that("GenerateVMatrix and GenerateLMatrix align DBS V and L matrices", {
  clinical_sample <- data.frame(
    SAMPLE_ID = unique(dbs_mut$Tumor_Sample_Barcode),
    SEQ_ASSAY_ID = SimData$PatientInfo$SEQ_ASSAY_ID[1],
    stringsAsFactors = FALSE
  )

  V <- SATS::GenerateVMatrix(dbs_mut, Class = "DBS", ref.genome = "hg19")
  L <- SATS::GenerateLMatrix(
    Panel_context = SimData$PanelEx,
    Patient_Info = clinical_sample,
    Class = "DBS",
    ref.genome = "hg19"
  )

  expect_true(identical(rownames(V), rownames(L)))
  expect_true(identical(colnames(V), colnames(L)))
  expect_equal(dim(V), dim(L))
})

test_that("GenerateLMatrix direct genomic-information input matches low-level workflow", {
  Panel_context <- SATS::GeneratePanelSize(SimData$PanelEx, Class = "SBS",
                                           SBS_order = "COSMIC",
                                           ref.genome = "hg19")
  L_low_level <- expect_warning(
    SATS::GenerateLMatrix(Panel_context, SimData$PatientInfo),
    "SEQ_ASSAY_ID values not present"
  )
  L_direct <- expect_warning(
    SATS::GenerateLMatrix(SimData$PanelEx, SimData$PatientInfo, Class = "SBS",
                          SBS_order = "COSMIC", ref.genome = "hg19"),
    "SEQ_ASSAY_ID values not present"
  )

  expect_equal(L_direct, L_low_level)
  expect_equal(ncol(L_direct), 20)
  expect_equal(setdiff(unique(SimData$PatientInfo$SEQ_ASSAY_ID),
                       colnames(Panel_context)), "VICC-01-T7")
  expect_false(ncol(L_direct) == ncol(SimData$L))
})
