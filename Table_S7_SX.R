library(gtools) 
library(DEoptim)
library(stats)
library(mgcv)
library(psych)
library(parallel)
library(foreach)
library(doParallel)
library(utils)
library(Matrix)
###########################example################################
Design <- matrix(c(4,0,2,3,1, 0,2,3,1,4, 2,3,1,4,0, 3,1,4,0,2, 1,4,0,2,3,
                   4,2,1,0,3, 2,1,0,3,4, 1,0,3,4,2, 0,3,4,2,1, 3,4,2,1,0),
                 nrow=10, byrow=TRUE)
Design <- Design + 1
###########################example################################
getPolyD <- function(Design, type = "F", weights = NULL, crit = "D",
                     returnModel = FALSE, returnWholeModel = FALSE) {
  # Here Design is supplied as a component matrix / OofA matrix.
  # Convert each row to the corresponding PofC / position matrix row.
  m <- ncol(Design)
  stopifnot(is.matrix(Design) || is.data.frame(Design))
  Design <- as.matrix(Design)
  if (!all(Design %in% seq_len(m))) {
    stop("Design must contain component labels 1,...,m.")
  }
  if (!all(apply(Design, 1, function(x) identical(sort(as.integer(x)), seq_len(m))))) {
    stop("Each row of Design must be a permutation of 1:m.")
  }
  for (i in seq_len(nrow(Design))) {
    Design[i, ] <- order(Design[i, ])
  }

  cont <- contr.poly(m) * sqrt(m)  # matches the paper's polynomial normalization

  trace_mat <- function(M) sum(diag(M))
  build_info <- function(Model, weights = NULL) {
    if (is.null(weights)) {
      return(crossprod(Model) / nrow(Model))
    }
    weights <- as.numeric(weights)
    if (length(weights) != nrow(Model)) {
      stop("weights must have length equal to nrow(Design).")
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be finite and nonnegative.")
    }
    if (sum(weights) <= 0) {
      stop("weights must sum to a positive value.")
    }
    weights <- weights / sum(weights)
    Info <- matrix(0, nrow = ncol(Model), ncol = ncol(Model))
    for (i in seq_len(nrow(Model))) {
      rowi <- Model[i, , drop = FALSE]
      Info <- Info + weights[i] * crossprod(rowi)
    }
    Info
  }

  # Intercept + linear terms.
  if (!returnWholeModel) {
    Model <- cbind(1, matrix(cont[Design[, -ncol(Design), drop = FALSE], 1], ncol = m - 1))
  } else {
    Model <- cbind(1, matrix(cont[Design, 1], ncol = m))
  }

  if (type == "F") {
    if (returnModel) return(Model[, -1, drop = FALSE])
    Info <- build_info(Model, weights)
    if (crit == "D") return(det(Info))
    if (crit == "A") return(trace_mat(solve(Info)))
    stop("crit must be 'D' or 'A'.")
  }

  # Quadratic terms.
  for (j in seq_len(m - 1)) {
    Model <- cbind(Model, cont[Design[, j], 2])
  }
  if (returnWholeModel) {
    Model <- cbind(Model, cont[Design[, m], 2])
  }

  if (type == "Q") {
    if (returnModel) return(Model[, -1, drop = FALSE])
    Info <- build_info(Model, weights)
    if (crit == "D") return(det(Info))
    if (crit == "A") return(trace_mat(solve(Info)))
    stop("crit must be 'D' or 'A'.")
  }

  # Second-order model.
  if (!returnWholeModel) {
    Model <- Model[, -ncol(Model), drop = FALSE]  # drop last quadratic term; not estimable
    for (k in seq_len(m - 2)) {
      for (l in (k + 1):(m - 1)) {
        Model <- cbind(Model, cont[Design[, k], 1] * cont[Design[, l], 1])
      }
    }
  } else {
    for (k in seq_len(m - 1)) {
      for (l in (k + 1):m) {
        Model <- cbind(Model, cont[Design[, k], 1] * cont[Design[, l], 1])
      }
    }
  }

  if (type == "S" || is.null(type)) {
    if (returnModel) return(Model[, -1, drop = FALSE])
    Info <- build_info(Model, weights)
    if (crit == "D") return(det(Info))
    if (crit == "A") return(trace_mat(solve(Info)))
    stop("crit must be 'D' or 'A'.")
  }

  stop("type must be one of 'F', 'Q', or 'S'.")
}


getCCP1D_parallel <- function(Design, q) {
  n <- nrow(Design)  
  m <- ncol(Design)  
  s <- max(q)
  
  M <- matrix(0, ncol = (m-1)*(s-1) + 1, nrow = (m-1)*(s-1) + 1)
  cl <- makeCluster(max(1L, min(25L, parallel::detectCores())))
  on.exit(stopCluster(cl), add = TRUE)
  registerDoParallel(cl)
  result <- foreach(i = 1:n, .combine = "+") %dopar% {
    g_xi <- matrix(0, ncol = (m-1)*(s-1) + 1, nrow = 1)
    g_xi[1, 1] <- 1  # 
    count <- 2  # 
    for(k in 1:(s-1)){
      for(j in 1:(m-1)){
        if(q[Design[i,j]] == k){
          g_xi[1, count] <- 1
        }
        count <- count + 1
      }
    }
    (t(g_xi) %*% g_xi) / n
  }
  M <- result
  eigenvalues <- eigen(M, symmetric = TRUE)$values  
  if (all(eigenvalues > 0)) {
    logD <- sum(log(eigenvalues)) 
    return(logD)
  } else {
    warning("The matrix M is not positive definite.")
    return(NULL)
  }
}

m <- ncol(Design)
FullDesign <- permutations(m, m, 1:m)
PolyD_F <- getPolyD(Design, type = "F")
print(PolyD_F)
PolyD_F_full <- getPolyD(FullDesign, type = "F")
print(PolyD_F_full)
PolyD_Q <- getPolyD(Design, type = "Q")
print(PolyD_Q)
PolyD_Q_full <- getPolyD(FullDesign, type = "Q")
print(PolyD_Q_full)
q <- c(1, 1, 2, 2, 3)
s <- max(q)
CCP1D <- getCCP1D_parallel(Design, q)
if (!is.numeric(CCP1D)) CCP1D <- NA
CCP1D_full <- getCCP1D_parallel(FullDesign, q)

Deff_PolyD_F <- NA
Deff_PolyD_F <- (PolyD_F / PolyD_F_full)^(1 / m)
print(Deff_PolyD_F)
Deff_PolyD_Q <- NA
Deff_PolyD_Q <- (PolyD_Q / PolyD_Q_full)^(1 / (2 * m - 1))
print(Deff_PolyD_Q)
Deff_CCP1D <- NA
if (!is.na(CCP1D)) Deff_CCP1D <- (exp(CCP1D - CCP1D_full))^(1/((m-1)*(s-1)+1))
print(Deff_CCP1D)