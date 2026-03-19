library(gtools)
Design <- matrix(c(
  0,1,2,3,4,5,6,7,8,
  1,2,0,4,5,3,7,8,6,
  2,0,1,5,3,4,8,6,7,
  3,4,5,6,7,8,0,1,2,
  4,5,3,7,8,6,1,2,0,
  5,3,4,8,6,7,2,0,1,
  6,7,8,0,1,2,3,4,5,
  7,8,6,1,2,0,4,5,3,
  8,6,7,2,0,1,5,3,4,
  0,2,1,6,8,7,3,5,4,
  1,0,2,7,6,8,4,3,5,
  2,1,0,8,7,6,5,4,3,
  3,5,4,0,2,1,6,8,7,
  4,3,5,1,0,2,7,6,8,
  5,4,3,2,1,0,8,7,6,
  6,8,7,3,5,4,0,2,1,
  7,6,8,4,3,5,1,0,2,
  8,7,6,5,4,3,2,1,0,
  0,3,6,4,7,1,8,2,5,
  1,4,7,5,8,2,6,0,3,
  2,5,8,3,6,0,7,1,4,
  3,6,0,7,1,4,2,5,8,
  4,7,1,8,2,5,0,3,6,
  5,8,2,6,0,3,1,4,7,
  6,0,3,1,4,7,5,8,2,
  7,1,4,2,5,8,3,6,0,
  8,2,5,0,3,6,4,7,1,
  0,4,8,7,2,3,5,6,1,
  1,5,6,8,0,4,3,7,2,
  2,3,7,6,1,5,4,8,0,
  3,7,2,1,5,6,8,0,4,
  4,8,0,2,3,7,6,1,5,
  5,6,1,0,4,8,7,2,3,
  6,1,5,4,8,0,2,3,7,
  7,2,3,5,6,1,0,4,8,
  8,0,4,3,7,2,1,5,6,
  0,5,7,1,3,8,2,4,6,
  1,3,8,2,4,6,0,5,7,
  2,4,6,0,5,7,1,3,8,
  3,8,1,4,6,2,5,7,0,
  4,6,2,5,7,0,3,8,1,
  5,7,0,3,8,1,4,6,2,
  6,2,4,7,0,5,8,1,3,
  7,0,5,8,1,3,6,2,4,
  8,1,3,6,2,4,7,0,5,
  0,6,3,8,5,2,4,1,7,
  1,7,4,6,3,0,5,2,8,
  2,8,5,7,4,1,3,0,6,
  3,0,6,2,8,5,7,4,1,
  4,1,7,0,6,3,8,5,2,
  5,2,8,1,7,4,6,3,0,
  6,3,0,5,2,8,1,7,4,
  7,4,1,3,0,6,2,8,5,
  8,5,2,4,1,7,0,6,3,
  0,7,5,2,6,4,1,8,3,
  1,8,3,0,7,5,2,6,4,
  2,6,4,1,8,3,0,7,5,
  3,1,8,5,0,7,4,2,6,
  4,2,6,3,1,8,5,0,7,
  5,0,7,4,2,6,3,1,8,
  6,4,2,8,3,1,7,5,0,
  7,5,0,6,4,2,8,3,1,
  8,3,1,7,5,0,6,4,2,
  0,8,4,5,1,6,7,3,2,
  1,6,5,3,2,7,8,4,0,
  2,7,3,4,0,8,6,5,1,
  3,2,7,8,4,0,1,6,5,
  4,0,8,6,5,1,2,7,3,
  5,1,6,7,3,2,0,8,4,
  6,5,1,2,7,3,4,0,8,
  7,3,2,0,8,4,5,1,6,
  8,4,0,1,6,5,3,2,7
),nrow = 72, byrow = TRUE) #example

# Function to calculate the encounter numbers (lambda values) for all pairs of runs
calculate_lambda <- function(Design) {
  encounter_numbers <- numeric()
  for (i in 1:(nrow(Design) - 1)) {
    for (j in (i + 1):nrow(Design)) {
      encounter_numbers <- c(encounter_numbers, sum(Design[i, ] == Design[j, ]))
    }
  }
  encounter_numbers
}

# Updated function to calculate a_k values
calculate_a_k <- function(k, m) {
  if (k < 4) {
    # Return the value for a_1, a_2, a_3 directly
    return(c(0, 1, 2)[k])
  } else {
    # Calculate a_k based on the summation formula for k >= 4
    a_k <- sum(sapply(4:k, function(j) (-1)^j * factorial(k) / (factorial(j) * factorial(k - j)))) + factorial(k) / 3
    return(a_k)
  }
}

# Function to compute delta for the discrepancy calculation
calculate_delta <- function(m, a, b) {
  delta <- a^m
  for (i in 1:(m-2)) {  
    k <- m - i
    a_k_value <- calculate_a_k(k, m)  # Remember to update this function with the correct logic for a_k
    #n_i <- a_k_value * choose(m, k)
    n_i <- a_k_value * choose(m, i)  # Updated to use m and i based on the new understanding
    delta <- delta + n_i * a^i * b^(m-i)
  }
  return(delta)
}

# Function to calculate discrepancy based on the lambda values and the delta
calculate_discrepancy <- function(Design, a, b) {
  m <- ncol(Design)
  n <- nrow(Design)
  lambda_values <- calculate_lambda(Design)
  delta <- calculate_delta(m, a, b)
  discrepancy <- (-delta / factorial(m)) + (a^m / n) + (2*b^m / (n^2)) * sum((a/b)^lambda_values)
  return(discrepancy)
}

# Parameters for discrepancy calculation
a <- 2
b <- 1

#计算所有偏差除离散偏差外都需要标准化设计矩阵
# 计算DD偏差
calculate_DD <- function(Design) {
  n <- nrow(Design)
  s <- ncol(Design)
  a <- 2  
  b <- 1
  qj <- apply(Design, 2, function(x) max(x) + 1)
  #  qj <- rep(s, s)
  DD <- 1
  # 计算DD的第一部分
  for (j in 1:s) {
    DD <- DD * (1/s + (qj[j] - 1)/(s * qj[j]))
  }
  
  # 计算DD的第二部分
  sum_term <- 0
  for (i in 1:n) {
    for (k in 1:n) {
      product <- 1
      for (j in 1:s) {
        delta <- as.integer(Design[i,j] == Design[k,j])
        product <- product * (a^delta) * (b^(1 - delta))
      }
      sum_term <- sum_term + product
    }
  }
  # 结合两部分计算最终的DD值
  DD <- sqrt(DD + (1/n^2) * sum_term)
  return(DD)
}

# 计算DD
DD_result <- calculate_DD(Design)
print(DD_result)

# 标准化设计矩阵
normalize_matrix <- function(Design) {
  # 获取矩阵中的最小值和最大值
  min_val <- min(Design)
  max_val <- max(Design)
  # 标准化矩阵到[0, 1]
  normalized_matrix <- (Design - min_val) / (max_val - min_val)
  return(normalized_matrix)
}
Design <- normalize_matrix(Design)

# 计算两点之间的L2距离
calculate_L2_distance <- function(row, permutation) {
  return(sqrt(sum((row - permutation)^2)))
}

# 修正 d_fill
calculate_d_fill <- function(Design) {
  n <- nrow(Design)
  m <- ncol(Design)
  permutations <- permn(0:(m-1))  # 修正为 0:(m-1)
  d_fill <- -Inf
  for (pi in permutations) {
    pi_norm <- pi / (m - 1)  # 标准化 pi 到 [0,1]
    min_distance <- Inf
    for (i in 1:n) {
      distance <- calculate_L2_distance(Design[i, ], pi_norm)
      min_distance <- min(min_distance, distance)
    }
    d_fill <- max(d_fill, min_distance)
  }
  return(d_fill)
}

result <- calculate_d_fill(Design)
print(result)  # 输出单个d_fill值

# 计算所有行之间的欧几里得距离
distance_matrix <- as.matrix(dist(Design))
# 找到非对角线元素中的最小距离（分离距离）
separation_distance <- min(distance_matrix[upper.tri(distance_matrix)])
# 打印分离距离
print(separation_distance)

calculate_separation_distance <- function(Design) {
  # 计算所有行之间的欧几里得距离
  distance_matrix <- as.matrix(dist(Design))
  # 找到非对角线元素中的最小距离（分离距离）
  separation_distance <- min(distance_matrix[upper.tri(distance_matrix)])
  return(separation_distance)
}
separation_distance <- calculate_separation_distance(Design)
print(separation_distance)

library(combinat)

calculate_MD <- function(Design) {
  n <- nrow(Design) 
  s <- ncol(Design) 
  # 第一部分
  term1 <- (4/3)^s
  # 第二部分
  sum_prod_3_xij2 <- 0
  for (i in 1:n) {
    prod_3_xij2 <- 1
    for (j in 1:s) {
      prod_3_xij2 <- prod_3_xij2 * (3 - Design[i, j]^2)
    }
    sum_prod_3_xij2 <- sum_prod_3_xij2 + prod_3_xij2
  }
  term2 <- (2^(1 - s) / n) * sum_prod_3_xij2
  # 第三部分
  sum_prod_2_max <- 0
  for (i in 1:n) {
    for (l in 1:n) {
      prod_2_max <- 1
      for (j in 1:s) {
        prod_2_max <- prod_2_max * (2 - max(Design[i, j], Design[l, j]))^2
      }
      sum_prod_2_max <- sum_prod_2_max + prod_2_max
    }
  }
  term3 <- (1/n^2) * sum_prod_2_max
  # 计算MD偏差
  MD <- sqrt(term1 - term2 + term3)
  return(MD)
}

# 计算MD偏差
MD_result <- calculate_MD(Design)
print(MD_result)

calculate_CD <- function(Design) {
  n <- nrow(Design) 
  s <- ncol(Design) 
  
  # 第一部分
  term1 <- (13/12)^s
  # 第二部分
  sum_prod_1 <- 0
  for (i in 1:n) {
    prod_term_1 <- 1
    for (j in 1:s) {
      prod_term_1 <- prod_term_1 * (1 + 0.5 * abs(Design[i, j] - 0.5) - 0.5 * (Design[i, j] - 0.5)^2)
    }
    sum_prod_1 <- sum_prod_1 + prod_term_1
  }
  term2 <- (2/n) * sum_prod_1
  # 第三部分
  sum_prod_2 <- 0
  for (i in 1:n) {
    for (k in 1:n) {
      prod_term_2 <- 1
      for (j in 1:s) {
        prod_term_2 <- prod_term_2 * (1 + 0.5 * abs(Design[i, j] - 0.5) + 0.5 * abs(Design[k, j] - 0.5) - 0.5 * abs(Design[i, j] - Design[k, j]))
      }
      sum_prod_2 <- sum_prod_2 + prod_term_2
    }
  }
  term3 <- (1/n^2) * sum_prod_2
  
  # 计算CD偏差
  CD <- sqrt(term1 - term2 + term3)
  return(CD)
}

# 计算CD偏差
CD_result <- calculate_CD(Design)
print(CD_result)

calculate_WD <- function(Design) {
  n <- nrow(Design) 
  s <- ncol(Design) 
  
  # 第一部分
  term1 <- (4/3)^s
  # 第二部分
  term2 <- (1/n) * (3/2)^s
  # 第三部分
  sum_prod_3 <- 0
  for (i in 1:(n-1)) {
    for (k in (i+1):n) {
      prod_term_3 <- 1
      for (j in 1:s) {
        prod_term_3 <- prod_term_3 * ((3/2) - abs(Design[i, j] - Design[k, j]) + abs(Design[i, j] - Design[k, j])^2)
      }
      sum_prod_3 <- sum_prod_3 + prod_term_3
    }
  }
  term3 <- (2/n^2) * sum_prod_3
  
  # 计算WD偏差
  WD <- sqrt(term1 - term2 + term3)
  return(WD)
}
# 计算WD偏差
WD_result <- calculate_WD(Design)
print(WD_result)

D <- Design
calculate_CD <- function(D, s) {
  n <- nrow(D)  
  m <- ncol(D)  
  s <- ncol(D)  
  Z <- (2 * D - s + 1) / (2 * s)
  sum1 <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      prod <- 1
      for (k in 1:m) {
        prod <- prod * (1 + 0.5 * abs(Z[i,k] + Z[j,k]) - 0.5 * abs(Z[i,k] - Z[j,k]))
      }
      sum1 <- sum1 + prod
    }
  }
  sum2 <- 0
  for (i in 1:n) {
    prod <- 1
    for (k in 1:m) {
      prod <- prod * (1 + 0.5 * abs(Z[i,k]) - 0.5 * abs(Z[i,k])^2)
    }
    sum2 <- sum2 + prod
  }
  CD_D <- (1/n^2) * sum1 - (2/n) * sum2 + (13/12)^m
  return(sqrt(CD_D))
}

calculate_phi <- function(D) {
  n <- nrow(D)  
  m <- ncol(D)  
  CD_values <- sapply(combn(m, 2), function(cols) {
    Du <- D[, cols, drop = FALSE]  
    #   print(paste(paste(cols, collapse = ",")))
    #   print(Du)
    
    CD_value <- calculate_CD(Du) 
    return(CD_value)  
  })
  
  # 计算φ(D)
  phi_D <- (2 / (m * (m - 1))) * sum(CD_values)
  return(phi_D)
}

# 调用函数计算CD
CD_result <- calculate_CD(D, s)
print(CD_result)
phi_result <- calculate_phi(D)
print(phi_result)
