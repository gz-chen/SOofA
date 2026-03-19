library(DEoptim)
library(stats)
library(mgcv)
library(psych)
library(parallel)
library(utils)
library(Matrix)
library(combinat)
library(gtools) # for permn

stratum1 <- c(1, 2) 
stratum2 <- c(3, 4) 
stratum3 <- c(5, 6)
stratum4 <- c(7, 8) 

perms1 <- list(c(1, 2), c(2, 1)) 
perms2 <- list(c(3, 4), c(4, 3)) 
perms3 <- list(c(5, 6), c(6, 5)) 
perms4 <- list(c(7, 8), c(8, 7)) 

Design <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8,  
                   3, 4, 1, 2, 7, 8, 5, 6,
                   5, 6, 7, 8, 1, 2, 3, 4,
                   7, 8, 5, 6, 3, 4, 1, 2,
                   1, 3, 5, 7, 6, 8, 2, 4,
                   3, 1, 7, 5, 8, 6, 4, 2,
                   5, 7, 1, 3, 2, 4, 6, 8,
                   7, 5, 3, 1, 4, 2, 8, 6,
                   1, 4, 7, 6, 2, 3, 8, 5,
                   3, 2, 5, 8, 4, 1, 6, 7,
                   5, 8, 3, 2, 6, 7, 4, 1,
                   7, 6, 1, 4, 8, 5, 2, 3,
                   1, 5, 6, 2, 8, 4, 3, 7,
                   3, 7, 8, 4, 6, 2, 1, 5,
                   5, 1, 2, 6, 4, 8, 7, 3,
                   7, 3, 4, 8, 2, 6, 5, 1,
                   1, 6, 8, 3, 4, 7, 5, 2,
                   3, 8, 6, 1, 2, 5, 7, 4,
                   5, 2, 4, 7, 8, 3, 1, 6,
                   7, 4, 2, 5, 6, 1, 3, 8,
                   1, 7, 2, 8, 3, 5, 4, 6,
                   3, 5, 4, 6, 1, 7, 2, 8,
                   5, 3, 6, 4, 7, 1, 8, 2,
                   7, 1, 8, 2, 5, 3, 6, 4,
                   1, 8, 4, 5, 7, 2, 6, 3,
                   3, 6, 2, 7, 5, 4, 8, 1,
                   5, 4, 8, 1, 3, 6, 2, 7,
                   7, 2, 6, 3, 1, 8, 4, 5),
                 nrow = 28, byrow = TRUE)


n_rows <- nrow(Design) # 28
n_cols <- ncol(Design) # 8
print(paste("Number of rows (N):", n_rows))
print(paste("Number of columns (m):", n_cols))

allowable_perms <- list()
idx <- 1
for (p1 in perms1) {
  for (p2 in perms2) {
    for (p3 in perms3) {
      for (p4 in perms4) {
        pi_map <- 1:n_cols 
        pi_map[stratum1] <- p1
        pi_map[stratum2] <- p2
        pi_map[stratum3] <- p3
        pi_map[stratum4] <- p4
        allowable_perms[[idx]] <- pi_map
        idx <- idx + 1
      }
    }
  }
}

print(allowable_perms) 
print(paste("Number of allowable permutations:", length(allowable_perms))) 

getPolyD <- function(Design, type = "F", weights = NULL, crit = "D",
                     returnModel = FALSE, returnWholeModel = FALSE) {
  # Design is assumed to be PofC / position matrix with entries 1,...,m.
  m <- ncol(Design)
  stopifnot(is.matrix(Design) || is.data.frame(Design))
  Design <- as.matrix(Design)
  if (!all(Design %in% seq_len(m))) {
    stop("Design must contain positions 1,...,m.")
  }
  if (!all(apply(Design, 1, function(x) identical(sort(as.integer(x)), seq_len(m))))) {
    stop("Each row of Design must be a permutation of 1:m (PofC / position matrix).")
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


num_perms <- length(allowable_perms)
vals_PolyD_F <- numeric(num_perms)
vals_PolyD_Q <- numeric(num_perms)
n_rows <- nrow(Design)
n_cols <- ncol(Design)
m <- n_cols
FullDesign <- permutations(m, m, 1:m)
apply_mapping_to_rows <- function(matrix_data) {
  m <- ncol(matrix_data)
  result_matrix <- matrix(nrow = nrow(matrix_data), ncol = m)
  for (i in 1:nrow(matrix_data)) {
    z <- matrix_data[i, ]  # OofA: z[i] = component at position i
    x <- numeric(m)        # PofC: x[j] = position of component j
    for (pos in 1:m) {
      comp <- z[pos]       # component at position pos
      x[comp] <- pos       # position of component comp
    }
    result_matrix[i, ] <- x
  }
  return(result_matrix)
}
FullDesign <- apply_mapping_to_rows(FullDesign)

stopifnot(all(apply(Design, 1, function(x) identical(sort(as.integer(x)), seq_len(m)))))

PolyD_F_full <- getPolyD(FullDesign, type = "F")
PolyD_Q_full <- getPolyD(FullDesign, type = "Q")

for (i in seq_along(allowable_perms)) {
  pi <- allowable_perms[[i]]
  Design_integer <- matrix(pi[Design], nrow = n_rows, ncol = n_cols)
  m <- n_cols
  PolyD_F <- getPolyD(Design_integer, type = "F")
  print(PolyD_F)
  print(PolyD_F_full)
  if (!is.na(PolyD_F) && !is.na(PolyD_F_full))
    vals_PolyD_F[i] <- (PolyD_F / PolyD_F_full)^(1 / m)
  PolyD_Q <- getPolyD(Design_integer, type = "Q")
  print(PolyD_Q)
  print(PolyD_Q_full)
  if (!is.na(PolyD_Q) && !is.na(PolyD_Q_full))
    vals_PolyD_Q[i] <- (PolyD_Q / PolyD_Q_full)^(1 / (2 * m - 1))
}

summary_stats <- data.frame(
  Criterion = c("PolyD_F", "PolyD_Q"),
  Mean = c(mean(vals_PolyD_F, na.rm = TRUE), mean(vals_PolyD_Q, na.rm = TRUE)),
  Variance = c(var(vals_PolyD_F, na.rm = TRUE), var(vals_PolyD_Q, na.rm = TRUE)),
  SD = c(sd(vals_PolyD_F, na.rm = TRUE), sd(vals_PolyD_Q, na.rm = TRUE))
)
print(summary_stats)