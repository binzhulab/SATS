
context("EstimateSigActivity")

dir <- system.file("extdata", package="SATS", mustWork=TRUE)

# Load the baseline objects obj0 and matrices L, V, W
rdafile <- file.path(dir, "test_objects", "test_EstimateSigActivity.rda")
load(rdafile)

# Load simulated data
data(SimData, package="SATS")

set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion") 

# Call the function to test
obj <- SATS::EstimateSigActivity(V, L, W)

# Compare result to the baseline. 
test_that("EstimateSigActivity",
{
  expect_equal(obj0$H, obj$H, tolerance=1e-4)
  expect_equal(obj0$loglike, obj$loglike, tolerance=1e-4)
  expect_equal(obj0$converged, obj$converged)
})

test_that("EstimateSigActivity aligns named V and L sample columns",
{
  V_small <- SimData$V[, 1:3, drop=FALSE]
  L_small <- SimData$L[, c(3, 1, 2), drop=FALSE]
  W_small <- SimData$TrueW_TMB

  obj_aligned <- SATS::EstimateSigActivity(V_small, L_small, W_small,
                                           n.start=1, iter.max=50)
  expect_identical(colnames(obj_aligned$H), colnames(V_small))
})

test_that("EstimateSigActivity rejects mismatched V and L sample IDs",
{
  V_small <- SimData$V[, 1:3, drop=FALSE]
  L_bad <- SimData$L[, 1:3, drop=FALSE]
  colnames(L_bad)[1] <- "sample_not_in_V"

  expect_error(
    SATS::EstimateSigActivity(V_small, L_bad, SimData$TrueW_TMB,
                              n.start=1, iter.max=50),
    "sample IDs"
  )
})

test_that("EstimateSigActivity ignores counts at zero-opportunity entries",
{
  V_zero <- matrix(c(5, 2), nrow=2,
                   dimnames=list(c("ctx1", "ctx2"), "sample1"))
  L_zero <- matrix(c(0, 1), nrow=2,
                   dimnames=list(c("ctx1", "ctx2"), "sample1"))
  W_zero <- matrix(c(0.5, 0.5), nrow=2,
                   dimnames=list(c("ctx1", "ctx2"), "sig1"))

  V_expected <- V_zero
  L_expected <- L_zero
  zero_opportunity <- L_expected == 0
  V_expected[zero_opportunity] <- 0
  L_expected[zero_opportunity] <- 1

  set.seed(2026)
  obj_zero <- SATS::EstimateSigActivity(V_zero, L_zero, W_zero,
                                        n.start=2, iter.max=50)
  set.seed(2026)
  obj_expected <- SATS::EstimateSigActivity(V_expected, L_expected, W_zero,
                                            n.start=2, iter.max=50)

  expect_equal(obj_zero$H, obj_expected$H, tolerance=1e-8)
  expect_true(is.finite(obj_zero$loglike))
})
