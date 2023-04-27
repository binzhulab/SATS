CalculateSigExpectancy <- function(L, W, H) {

  check_L_W_H(L, W, H) 

  N   <- ncol(H)
  ret <- sigExp_setReturn(L, W, H)

  if (!is.matrix(L)) L <- as.matrix(L)
  if (!is.matrix(W)) W <- as.matrix(W)
  H <- t(as.matrix(H))

  for(n in 1:N){
    tmp      <- L[, n, drop=FALSE] %*% H[n, , drop=FALSE]
    ret[, n] <- colSums(tmp*W)
  }

  ret

}

sigExp_setReturn <- function(L, W, H) {

  N   <- ncol(H)
  K   <- ncol(W)
  ret <- matrix(data=NA, nrow=K, ncol=N)

  nms <- colnames(L)
  if (is.null(nms)) nms <- colnames(H) 
  if (is.null(nms)) nms <- paste0("Sample", 1:N)
  colnames(ret) <- nms

  nms <- colnames(W)
  if (is.null(nms)) nms <- paste0("deNovo", LETTERS[1:K])
  rownames(ret) <- nms
  
  ret

}

EstimateSigActivity <- function(V, L, W, n.start=50, iter.max=5000, eps=1e-5) {

  check_L_W_V(L, W, V)
  check_number(n.start, "n.start", min=1)
  check_number(iter.max, "iter.max", min=1)
  check_number(eps, "eps", pos=TRUE)
  op  <- list(n.start=n.start, iter.max=iter.max, eps=eps, print=0)
  ret <- estSigAct_main(V, L, W, op) 
  ret
}

estSigAct_main <- function(V, L, W, op) {

  DEBUG <- 0
  if (DEBUG) {
    print(paste0("nrow(V)=", nrow(V), ", ncol(V)=", ncol(V)))
    print(paste0("nrow(L)=", nrow(L), ", ncol(L)=", ncol(L)))
    print(paste0("nrow(W)=", nrow(W), ", ncol(W)=", ncol(W)))
  }
  n     <- ncol(V)
  k     <- ncol(W)
  p     <- nrow(V)
  iargs <- c(n, k, p, op$n.start, op$iter.max, op$print, DEBUG)

  L[L == 0] <- 1
  V[L == 0] <- 0
  lower     <- 1e-6
  upper     <- sum(V)/sum(L)
  dargs     <- c(op$eps, lower, upper)

  ret_H     <- rep(-9999.0e200, k*n)
  ret_ll    <- -9999.0
  ret_conv  <- 0

  tmp <- .C("C_call_salmon", as.integer(iargs), as.numeric(dargs), 
                       as.numeric(t(V)), as.numeric(t(L)), as.numeric(t(W)), 
                       ret_H=as.numeric(ret_H), ret_ll=as.numeric(ret_ll), 
                       ret_conv=as.integer(ret_conv),
                       PACKAGE="SATS")

  ret_H  <- matrix(tmp$ret_H, nrow=k, ncol=n, byrow=TRUE)
  ret_H  <- setNames_H(ret_H, V, L, W)
  ret_ll <- tmp$ret_ll
  conv   <- tmp$ret_conv
  list(H=ret_H, loglike=ret_ll, converged=conv)

}

setNames_H <- function(H, V, L, W) {

  nms <- colnames(V)
  if (is.null(nms)) nms <- colnames(L)
  if (is.null(nms)) nms <- paste0("Sample", 1:ncol(V))
  colnames(H) <- nms
 
  nms <- colnames(W)
  if (is.null(nms)) nms <- paste0("deNovo", LETTERS[1:ncol(W)])
  rownames(H) <- nms
 
  H
}

EstimateSigActivity_R <- function(V, L, W, n.start=50, iter.max=5000, eps=1e-5) {

  check_L_W_V(L, W, V)
  check_number(n.start, "n.start", min=1)
  check_number(iter.max, "iter.max", min=1)
  check_number(eps, "eps", pos=TRUE)
  ret <- SALMON_H(n.start, V, L, W, iter.max=iter.max, EPS=eps)
  ret
}

SALMON_H <- function(nstart, V, L, W0, iter.max=5000, EPS=1e-5){
  
  ret       <- NULL
  maxll     <- -Inf
  V         <- as.matrix(V)
  L         <- as.matrix(L)
  W0        <- as.matrix(W0)
  L[L == 0] <- 1
  V[L == 0] <- 0
  upper     <- sum(V)/sum(L)

  for(i in 1:nstart)
  {
    tmp <- try(SALMON_H_1(V, L, W0, upper, iter.max=iter.max, EPS=EPS))
    if (!("try-error" %in% class(tmp))) {
      ll <- tmp$loglike
      if (ll > maxll) {
        ret   <- tmp
        maxll <- ll
      }
    }
  }
  
  ret
}

SALMON_H_1 <- function(V, L, W0, upper, iter.max=5000, EPS=1e-5){
  
  n  <- ncol(V)
  p  <- nrow(V) 
  r  <- ncol(W0)
  np <- n*p
  H0 <- matrix(runif(r*n,min=1e-6,max = upper),nrow = r, ncol=n, byrow=TRUE)

  #update the expectation
  WH.hat <- W0 %*% H0
  V.hat  <- WH.hat * L
  
  #calculate observed likelihood
  Lik_mat <- V*log(V.hat)-V.hat 
  lik0    <- sum(Lik_mat)/np
  if (!is.finite(lik0)) stop("Initial likelihood is not finite")
  #print(paste0("ll0=", lik0))
  
  iter     <- 0
  conv     <- FALSE
  trans.W0 <- t(W0)
  tWL      <- trans.W0 %*% L
  
  while(1){
    iter <- iter + 1
    #cat("Iteration:",iter, "\n")
    
    ##M step: update H 
    H1 <- H0*( trans.W0 %*% (V/(WH.hat)))/tWL

    #E step: update the expectation
    WH.hat <- W0 %*% H1
    V.hat  <- WH.hat * L

    #calculate observed likelihood
    Lik_mat <- V*log(V.hat)-V.hat 
    lik1    <- sum(Lik_mat)/np #normalized by number of obs
    #print(paste0("ll=", lik1))
    if (!is.finite(lik1)) return(NULL)
    
    #cat(lik1, "\n")
    if (abs((lik0-lik1)/lik0) < EPS) {
      conv <- TRUE
      break
    }
    if (iter >= iter.max) break
    H0   <- H1
    lik0 <- lik1
  }
  #print(paste0("final loglike = ", lik1))
  return(list(H=H1, loglike=lik1, converged=conv))
}


