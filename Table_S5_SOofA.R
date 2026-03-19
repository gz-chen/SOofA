library(combinat)
library(gtools) 

stratum1 <- c(0, 1)    # {0, 1}
stratum2 <- c(2, 3)    # {2, 3}
stratum3 <- c(4, 5)    # {4, 5}
stratum4 <- c(6, 7)    # {6, 7}

perms1 <- list(c(0, 1), c(1, 0))  # {0, 1} 
perms2 <- list(c(2, 3), c(3, 2))  # {2, 3} 
perms3 <- list(c(4, 5), c(5, 4))  # {4, 5} 
perms4 <- list(c(6, 7), c(7, 6))  # {6, 7} 

Design <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7,
                   2, 3, 0, 1, 6, 7, 4, 5,
                   4, 5, 6, 7, 0, 1, 2, 3,
                   6, 7, 4, 5, 2, 3, 0, 1,
                   0, 2, 4, 6, 5, 7, 1, 3,
                   2, 0, 6, 4, 7, 5, 3, 1,
                   4, 6, 0, 2, 1, 3, 5, 7,
                   6, 4, 2, 0, 3, 1, 7, 5,
                   0, 3, 6, 5, 1, 2, 7, 4,
                   2, 1, 4, 7, 3, 0, 5, 6,
                   4, 7, 2, 1, 5, 6, 3, 0,
                   6, 5, 0, 3, 7, 4, 1, 2,
                   0, 4, 5, 1, 7, 3, 2, 6,
                   2, 6, 7, 3, 5, 1, 0, 4,
                   4, 0, 1, 5, 3, 7, 6, 2,
                   6, 2, 3, 7, 1, 5, 4, 0,
                   0, 5, 7, 2, 3, 6, 4, 1,
                   2, 7, 5, 0, 1, 4, 6, 3,
                   4, 1, 3, 6, 7, 2, 0, 5,
                   6, 3, 1, 4, 5, 0, 2, 7,
                   0, 6, 1, 7, 2, 4, 3, 5,
                   2, 4, 3, 5, 0, 6, 1, 7,
                   4, 2, 5, 3, 6, 0, 7, 1,
                   6, 0, 7, 1, 4, 2, 5, 3,
                   0, 7, 3, 4, 6, 1, 5, 2,
                   2, 5, 1, 6, 4, 3, 7, 0,
                   4, 3, 7, 0, 2, 5, 1, 6,
                   6, 1, 5, 2, 0, 7, 3, 4),
                 nrow = 28, byrow = TRUE)

n_rows <- nrow(Design)  # 28
n_cols <- ncol(Design)  # 8
print(paste("Number of rows (N):", n_rows))
print(paste("Number of columns (m):", n_cols))

allowable_perms <- list()
idx <- 1
for (p1 in perms1) {
  for (p2 in perms2) {
    for (p3 in perms3) {
      for (p4 in perms4) {
        pi_map <- 0:(n_cols - 1)  
        pi_map[stratum1 + 1] <- p1
        pi_map[stratum2 + 1] <- p2
        pi_map[stratum3 + 1] <- p3
        pi_map[stratum4 + 1] <- p4
        allowable_perms[[idx]] <- pi_map
        idx <- idx + 1
      }
    }
  }
}

print(allowable_perms)  
print(paste("Number of allowable permutations:", length(allowable_perms))) 

calculate_DD <- function(Design) {
  n <- nrow(Design)
  s <- ncol(Design)
  a <- 2
  b <- 1
  qj <- apply(Design, 2, function(x) max(x) + 1)
  DD <- 1
  for (j in 1:s) {
    DD <- DD * (1/s + (qj[j] - 1)/(s * qj[j]))
  }
  sum_term <- 0
  for (ii in 1:n) {
    for (k in 1:n) {
      product <- 1
      for (j in 1:s) {
        delta <- as.integer(Design[ii,j] == Design[k,j])
        product <- product * (a^delta) * (b^(1 - delta))
      }
      sum_term <- sum_term + product
    }
  }
  DD <- sqrt(DD + (1/n^2) * sum_term)
  return(DD)
}

normalize_matrix <- function(Design) {
  min_val <- min(Design)
  max_val <- max(Design)
  normalized_matrix <- (Design - min_val) / (max_val - min_val)
  return(normalized_matrix)
}

calculate_separation_distance <- function(Design) {
  distance_matrix <- as.matrix(dist(Design))
  separation_distance <- min(distance_matrix[upper.tri(distance_matrix)])
  return(separation_distance)
}

calculate_MD <- function(Design) {
  n <- nrow(Design)
  s <- ncol(Design)
  term1 <- (4/3)^s
  sum_prod_3_xij2 <- 0
  for (ii in 1:n) {
    prod_3_xij2 <- 1
    for (j in 1:s) {
      prod_3_xij2 <- prod_3_xij2 * (3 - Design[ii, j]^2)
    }
    sum_prod_3_xij2 <- sum_prod_3_xij2 + prod_3_xij2
  }
  term2 <- (2^(1 - s) / n) * sum_prod_3_xij2
  sum_prod_2_max <- 0
  for (ii in 1:n) {
    for (l in 1:n) {
      prod_2_max <- 1
      for (j in 1:s) {
        prod_2_max <- prod_2_max * (2 - max(Design[ii, j], Design[l, j]))^2
      }
      sum_prod_2_max <- sum_prod_2_max + prod_2_max
    }
  }
  term3 <- (1/n^2) * sum_prod_2_max
  MD <- sqrt(term1 - term2 + term3)
  return(MD)
}

calculate_CD <- function(Design) {
  n <- nrow(Design)
  s <- ncol(Design)
  term1 <- (13/12)^s
  sum_prod_1 <- 0
  for (ii in 1:n) {
    prod_term_1 <- 1
    for (j in 1:s) {
      prod_term_1 <- prod_term_1 * (1 + 0.5 * abs(Design[ii, j] - 0.5) - 0.5 * (Design[ii, j] - 0.5)^2)
    }
    sum_prod_1 <- sum_prod_1 + prod_term_1
  }
  term2 <- (2/n) * sum_prod_1
  sum_prod_2 <- 0
  for (ii in 1:n) {
    for (k in 1:n) {
      prod_term_2 <- 1
      for (j in 1:s) {
        prod_term_2 <- prod_term_2 * (1 + 0.5 * abs(Design[ii, j] - 0.5) + 0.5 * abs(Design[k, j] - 0.5) - 0.5 * abs(Design[ii, j] - Design[k, j]))
      }
      sum_prod_2 <- sum_prod_2 + prod_term_2
    }
  }
  term3 <- (1/n^2) * sum_prod_2
  CD <- sqrt(term1 - term2 + term3)
  return(CD)
}

calculate_WD <- function(Design) {
  n <- nrow(Design)
  s <- ncol(Design)
  term1 <- (4/3)^s
  term2 <- (1/n) * (3/2)^s
  sum_prod_3 <- 0
  for (ii in 1:(n-1)) {
    for (k in (ii+1):n) {
      prod_term_3 <- 1
      for (j in 1:s) {
        prod_term_3 <- prod_term_3 * ((3/2) - abs(Design[ii, j] - Design[k, j]) + abs(Design[ii, j] - Design[k, j])^2)
      }
      sum_prod_3 <- sum_prod_3 + prod_term_3
    }
  }
  term3 <- (2/n^2) * sum_prod_3
  WD <- sqrt(term1 - term2 + term3)
  return(WD)
}

calculate_L2_distance <- function(row, permutation) {
  return(sqrt(sum((row - permutation)^2)))
}

calculate_d_fill <- function(Design) {
  n <- nrow(Design)
  m <- ncol(Design)
  permutations <- permn(0:(m-1)) 
  d_fill <- -Inf
  for (pi in permutations) {
    pi_norm <- pi / (m - 1)  
    min_distance <- Inf
    for (i in 1:n) {
      distance <- calculate_L2_distance(Design[i, ], pi_norm)
      min_distance <- min(min_distance, distance)
    }
    d_fill <- max(d_fill, min_distance)
  }
  return(d_fill)
}

num_perms <- length(allowable_perms)
vals_DD <- numeric(num_perms)
vals_separation <- numeric(num_perms)
vals_MD <- numeric(num_perms)
vals_CD1 <- numeric(num_perms)
vals_WD <- numeric(num_perms)
vals_d_fill <- numeric(num_perms)  


for (i in seq_along(allowable_perms)) {
  pi <- allowable_perms[[i]]
  
  # relabel: π o X 
  Design_integer <- matrix(pi[Design + 1], nrow = n_rows, ncol = n_cols) 
  
  #  DD 
  vals_DD[i] <- calculate_DD(Design_integer)
  
  Design_norm <- normalize_matrix(Design_integer)
  
  # separation_distance
  vals_separation[i] <- calculate_separation_distance(Design_norm)
  
  # MD
  vals_MD[i] <- calculate_MD(Design_norm)
  
  # CD1
  vals_CD1[i] <- calculate_CD(Design_norm)
  
  # WD
  vals_WD[i] <- calculate_WD(Design_norm)
  
  # d_fill 
  vals_d_fill[i] <- calculate_d_fill(Design_norm)
}

summary_stats <- data.frame(
  Criterion = c("DD", "separation_distance", "MD", "CD1", "WD", "d_fill"),
  Mean = c(mean(vals_DD), mean(vals_separation), mean(vals_MD), mean(vals_CD1), mean(vals_WD), mean(vals_d_fill)),
  Variance = c(var(vals_DD), var(vals_separation), var(vals_MD), var(vals_CD1), var(vals_WD), var(vals_d_fill)),
  SD = c(sd(vals_DD), sd(vals_separation), sd(vals_MD), sd(vals_CD1), sd(vals_WD), sd(vals_d_fill))
)
print(summary_stats)