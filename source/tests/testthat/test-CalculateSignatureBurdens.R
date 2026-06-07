
context("CalculateSignatureBurdens")

dir <- system.file("extdata", package="SATS", mustWork=TRUE)

# Load the baseline object obj0
rdafile <- file.path(dir, "test_objects", "test_CalculateSignatureBurdens.rda")
load(rdafile)

# Load simulated data
data(SimData, package="SATS")

# Call the function to test
obj <- SATS::CalculateSignatureBurdens(SimData$L, SimData$TrueW_TMB, SimData$TrueH)


# Compare result to the baseline. 
test_that("CalculateSignatureBurdens",
{
  expect_equal(obj0, obj, tolerance=1e-4)
})

test_that("CalculateSignatureBurdens aligns named L, W and H matrices",
{
  H_small <- SimData$TrueH[, 1:3, drop=FALSE]
  L_small <- SimData$L[, c(3, 1, 2), drop=FALSE]
  W_small <- SimData$TrueW_TMB[, rev(seq_len(ncol(SimData$TrueW_TMB))), drop=FALSE]

  obj_aligned <- SATS::CalculateSignatureBurdens(L_small, W_small, H_small)

  expect_identical(colnames(obj_aligned), colnames(H_small))
  expect_identical(rownames(obj_aligned), rownames(H_small))
})

test_that("CalculateSignatureBurdens rejects mismatched L and H sample IDs",
{
  H_small <- SimData$TrueH[, 1:3, drop=FALSE]
  L_bad <- SimData$L[, 1:3, drop=FALSE]
  colnames(L_bad)[1] <- "sample_not_in_H"

  expect_error(
    SATS::CalculateSignatureBurdens(L_bad, SimData$TrueW_TMB, H_small),
    "sample IDs"
  )
})

