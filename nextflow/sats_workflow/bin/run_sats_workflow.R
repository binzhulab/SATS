#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

n_samples <- as.integer(get_arg("--n-samples", "25"))
out_prefix <- get_arg("--out-prefix", "sats")

if (is.na(n_samples) || n_samples < 1) {
  stop("ERROR: --n-samples must be a positive integer")
}

suppressPackageStartupMessages(library(SATS))

data(SimData, package = "SATS")
data(RefTMB, package = "SATS")

V_all <- SimData$V
L_all <- SimData$L

validated_all <- ValidateSATSInputs(V = V_all, L = L_all)
V_all <- validated_all$V
L_all <- validated_all$L

n_use <- min(n_samples, ncol(V_all))
sample_id <- colnames(V_all)[seq_len(n_use)]
V_mat <- V_all[, sample_id, drop = FALSE]
L_mat <- L_all[, sample_id, drop = FALSE]

validated_subset <- ValidateSATSInputs(V = V_mat, L = L_mat)
V_mat <- validated_subset$V
L_mat <- validated_subset$L

reference_names <- intersect(c("SBS1", "SBS5"), colnames(RefTMB$TMB_SBS_v3.4))
if (length(reference_names) < 2L) {
  reference_names <- colnames(RefTMB$TMB_SBS_v3.4)[seq_len(2L)]
}

# This compact workflow uses known reference profiles as stand-ins for de novo
# profiles so the containerized example remains fast and deterministic. Full
# analyses should replace W_hat with de novo profiles from a cohort-level
# extraction step.
W_hat <- as.matrix(RefTMB$TMB_SBS_v3.4[, reference_names, drop = FALSE])

set.seed(1)
MappedSig <- MappingSignature(
  W_hat = W_hat,
  W_ref = RefTMB$TMB_SBS_v3.4,
  niter = 10,
  cutoff.I2 = 0.1,
  min.repeats = 8
)

if (is.null(MappedSig) || nrow(MappedSig) == 0L) {
  stop("ERROR: MappingSignature did not select any reference signatures")
}

SBS.list <- unique(MappedSig$Reference)
W_star <- as.matrix(RefTMB$TMB_SBS_v3.4[, SBS.list, drop = FALSE])

validated_refit <- ValidateSATSInputs(V = V_mat, L = L_mat, W = W_star)
V_mat <- validated_refit$V
L_mat <- validated_refit$L
W_star <- validated_refit$W

set.seed(1)
H_hat <- EstimateSigActivity(
  V = V_mat,
  L = L_mat,
  W = W_star,
  n.start = 5,
  iter.max = 1000,
  eps = 1e-5
)

SigBdn <- CalculateSignatureBurdens(
  L = L_mat,
  W = W_star,
  H = H_hat$H
)

validated_burden <- ValidateSATSInputs(L = L_mat, W = W_star, H = H_hat$H)
L_mat <- validated_burden$L
W_star <- validated_burden$W

write.csv(MappedSig, paste0(out_prefix, "_mapping_results.csv"),
          row.names = FALSE)
write.csv(H_hat$H, paste0(out_prefix, "_activity_matrix.csv"))
write.csv(SigBdn, paste0(out_prefix, "_signature_burdens.csv"))

saveRDS(
  list(
    n_samples = n_use,
    mapped_signatures = MappedSig,
    activity = H_hat$H,
    signature_burdens = SigBdn,
    loglike = H_hat$loglike,
    converged = H_hat$converged
  ),
  file = paste0(out_prefix, "_workflow_outputs.rds")
)

writeLines(
  c(
    "SATS containerized workflow example completed.",
    paste0("Samples analyzed: ", n_use),
    paste0("Mapped signatures: ", paste(SBS.list, collapse = ", ")),
    paste0("EM converged: ", H_hat$converged),
    paste0("Log-likelihood: ", signif(H_hat$loglike, 6))
  ),
  con = paste0(out_prefix, "_workflow_summary.txt")
)
