
context("ValidateSATSInputs")

data(SimData, package="SATS")

V_small <- SimData$V[, 1:3, drop=FALSE]
L_small <- SimData$L[, 1:3, drop=FALSE]
W_small <- SimData$TrueW_TMB
H_small <- SimData$TrueH[, 1:3, drop=FALSE]

test_that("ValidateSATSInputs aligns compatible named SATS matrices",
{
  L_reordered <- L_small[, c(3, 1, 2), drop=FALSE]
  W_reordered <- W_small[, rev(seq_len(ncol(W_small))), drop=FALSE]

  obj <- SATS::ValidateSATSInputs(V=V_small, L=L_reordered, W=W_reordered)

  expect_identical(colnames(obj$L), colnames(V_small))
  expect_identical(rownames(obj$L), rownames(V_small))
  expect_identical(rownames(obj$W), rownames(V_small))
})

test_that("ValidateSATSInputs aligns L, W and H for burden calculation",
{
  L_reordered <- L_small[, c(2, 3, 1), drop=FALSE]
  W_reordered <- W_small[, rev(seq_len(ncol(W_small))), drop=FALSE]

  obj <- SATS::ValidateSATSInputs(L=L_reordered, W=W_reordered, H=H_small)

  expect_identical(colnames(obj$L), colnames(H_small))
  expect_identical(colnames(obj$W), rownames(H_small))
})

test_that("ValidateSATSInputs rejects malformed numeric inputs",
{
  V_negative <- V_small
  V_negative[1, 1] <- -1
  expect_error(SATS::ValidateSATSInputs(V=V_negative, L=L_small, W=W_small),
               "non-negative")

  V_decimal <- V_small
  V_decimal[1, 1] <- V_decimal[1, 1] + 0.25
  expect_error(SATS::ValidateSATSInputs(V=V_decimal, L=L_small, W=W_small),
               "integer-like")

  W_character <- W_small
  W_character[1, 1] <- NA
  storage.mode(W_character) <- "character"
  expect_error(SATS::ValidateSATSInputs(V=V_small, L=L_small, W=W_character),
               "numeric")
})

test_that("ValidateSATSInputs rejects duplicate or missing identifiers",
{
  V_duplicate <- V_small
  colnames(V_duplicate)[2] <- colnames(V_duplicate)[1]
  expect_error(SATS::ValidateSATSInputs(V=V_duplicate, L=L_small, W=W_small),
               "duplicate")

  L_missing_name <- L_small
  colnames(L_missing_name)[1] <- ""
  expect_error(SATS::ValidateSATSInputs(V=V_small, L=L_missing_name, W=W_small),
               "missing")

  W_duplicate <- W_small
  colnames(W_duplicate)[2] <- colnames(W_duplicate)[1]
  expect_error(SATS::ValidateSATSInputs(L=L_small, W=W_duplicate, H=H_small),
               "duplicate")
})

test_that("ValidateSATSInputs rejects mismatched named axes",
{
  L_bad_context <- L_small
  rownames(L_bad_context)[1] <- "not_a_context"
  expect_error(SATS::ValidateSATSInputs(V=V_small, L=L_bad_context, W=W_small),
               "mutation context names")

  H_bad_signature <- H_small
  rownames(H_bad_signature)[1] <- "not_a_signature"
  expect_error(SATS::ValidateSATSInputs(L=L_small, W=W_small, H=H_bad_signature),
               "signature names")
})

test_that("Core outputs preserve names, dimensions and non-negative values",
{
  H_hat <- SATS::EstimateSigActivity(V_small, L_small, W_small,
                                     n.start=1, iter.max=50)
  expect_true(all(c("H", "loglike", "converged") %in% names(H_hat)))
  expect_identical(dim(H_hat$H), c(ncol(W_small), ncol(V_small)))
  expect_identical(rownames(H_hat$H), colnames(W_small))
  expect_identical(colnames(H_hat$H), colnames(V_small))
  expect_true(all(is.finite(H_hat$H)))
  expect_true(all(H_hat$H >= 0))

  burdens <- SATS::CalculateSignatureBurdens(L_small, W_small, H_hat$H)
  expect_identical(dim(burdens), c(ncol(W_small), ncol(L_small)))
  expect_identical(rownames(burdens), colnames(W_small))
  expect_identical(colnames(burdens), colnames(L_small))
  expect_true(all(is.finite(burdens)))
  expect_true(all(burdens >= 0))
})
