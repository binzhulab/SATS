check_L_W_H <- function(L, W, H) {

  check_mat_df(L, "L")
  check_mat_df(W, "W")
  check_mat_df(H, "H")
  if (nrow(H) != ncol(W)) stop("ERROR: nrow(H) != ncol(W)")
  if (ncol(L) != ncol(H)) stop("ERROR: ncol(H) != ncol(L)")
  if (nrow(L) != nrow(W)) stop("ERROR: nrow(L) != nrow(W)")

  NULL
}

check_L_W_V <- function(L, W, V) {

  check_mat_df(L, "L")
  check_mat_df(W, "W")
  check_mat_df(V, "V")
  if (any(dim(L) != dim(V))) stop("ERROR: dim(L) != dim(V)")
  if (nrow(L) != nrow(W)) stop("ERROR: nrow(L) != nrow(W)")

  NULL
}


check_mat_df <- function(x, nm) {

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(paste0("ERROR: ", nm, " must be a matrix or data frame"))
  }
  if (!nrow(x)) stop(paste0("ERROR: ", nm, " has no rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", nm, " has no columns"))
  if (any(!is.finite(as.matrix(x)))) {
    stop(paste0("ERROR: ", nm, " contains non-finite values"))
  }

  NULL
}

check_number <- function(x, nm, min=NULL, pos=FALSE, max=NULL) {

  if (length(x) != 1) stop(paste0("ERROR: ", nm, " must be a single numeric value"))
  if (!is.null(min) && (x < min)) {
    stop(paste0("ERROR: ", nm, " must be greater than or equal to ", min))
  }
  if (pos && (x <= 0)) {
    stop(paste0("ERROR: ", nm, " must be positive"))
  }
  if (!is.null(max) && (x > max)) {
    stop(paste0("ERROR: ", nm, " must be less than or equal to ", max))
  }

  NULL

}