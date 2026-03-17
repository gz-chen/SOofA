library(DEoptim)
library(stats)
library(mgcv)
library(psych)
library(parallel)
library(foreach)
library(doParallel)
library(utils)
library(Matrix)
library(combinat)
library(gtools) 
###################################example############################################
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

getCCP1D_parallel <- function(Design, q) {
  n <- nrow(Design)  
  m <- ncol(Design)  
  s <- max(q)
  M <- matrix(0, ncol = (m-1)*(s-1) + 1, nrow = (m-1)*(s-1) + 1)
  cl <- makeCluster(25)
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
  M <- Reduce(`+`, result)
  stopCluster(cl)
  eigenvalues <- eigen(M, symmetric = TRUE)$values  
  if (all(eigenvalues > 0)) {
    logD <- sum(log(eigenvalues)) 
    return(logD)
  } else {
    warning("The matrix M is not positive definite.")
    return(NULL)
  }
}

num_perms <- length(allowable_perms)
vals_Deff <- numeric(num_perms)
sigma <- c(1,1,2,2,3,3,4,4)  # 1-based
s <- max(sigma)
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

for (i in seq_along(allowable_perms)) {
  pi <- allowable_perms[[i]]
  Design_integer <- matrix(pi[Design], nrow = n_rows, ncol = n_cols) 
  m <- ncol(Design_integer)
  CCP1D <- getCCP1D_parallel(Design_integer, sigma)
  print(CCP1D)
  CCP1D_full <- getCCP1D_parallel(FullDesign, sigma)
  print(CCP1D_full)
  if (!is.na(CCP1D))
    vals_Deff[i] <- (exp(CCP1D - CCP1D_full))^(1/((m-1)*(s-1)+1))
}

summary_stats <- data.frame(
  Criterion = c("Deff"),
  Mean = c(mean(vals_Deff)),
  Variance = c(var(vals_Deff)),
  SD = c(sd(vals_Deff))
)

print(summary_stats)