library(magrittr)
library(tidyverse)

## Parameters
iterations_per_alpha <- 100
N                    <- 10000
sigma                <- 0.1
bin_partitions       <- 10
K                    <- 1000
alpha_levels         <- 10^seq(-3, 3, length.out = 21)

result <- numeric(length(alpha_levels))
for(i in 1:length(alpha_levels)) {
  cat("Running with alpha = ", alpha_levels[i], "\n")
  proportion <- 0
  for(j in 1:iterations_per_alpha) {
    cat("Iteration ", j, "/", iterations_per_alpha)
    
    delta <- bin_partitions^2 / (N * exp(alpha_levels[i] / K) + bin_partitions^2 - N)
    alpha <- K * log((1 - delta) * bin_partitions^2 / (N * delta) + 1)
    cat(" ; delta = ", delta, ", alpha =", alpha, "\n")
    
    ## Generate mock data
    X1  <- runif(N)
    eps <- rnorm(N, sd = sigma, mean = 0)
    X2 <- X1 + eps
    X <- data.frame(X1, X2)
    while((X %>% 
           filter(X1 <= 1, X1 >= 0, 
                  X2 <= 1, X2 >= 0) %>% 
           nrow) < N) {
      
      temp <- 
        X %>% 
        filter(X1 > 1 | X1 < 0 | X2 > 1 | X2 < 0) %>% 
        nrow
      X1  <- runif(temp)
      eps <- rnorm(temp, sd = sigma, mean = 0)
      X2  <- X1 + eps
      
      temp <- data.frame(X1, X2)
      X %<>% 
        filter(X1 <= 1, X1 >= 0, 
               X2 <= 1, X2 >= 0) %>%
        rbind(temp)
    }
    
    ## Generate bins
    X1 <- numeric(bin_partitions^2)
    Y1 <- numeric(bin_partitions^2)
    X2 <- numeric(bin_partitions^2)
    Y2 <- numeric(bin_partitions^2)
    
    # Partition the intervals
    binx     <- seq(from = 0, to = 1, length.out = bin_partitions + 1)
    biny     <- seq(from = 0, to = 1, length.out = bin_partitions + 1)
    
    num <- 1
    for(k in 1:(length(biny) - 1)) {
      for(l in 1:(length(binx) - 1)) {
        X1[num] <- binx[l]
        Y1[num] <- biny[k]
        X2[num] <- binx[l + 1]
        Y2[num] <- biny[k + 1]
        num <- num + 1
      }
    }
    
    bins <- data.frame(X1, Y1, X2, Y2) 
    
    ## Generate the histogram
    C <- matrix(nrow = bin_partitions, ncol = bin_partitions)
    
    rownum <- 1
    for(k in 0:(nrow(bins) - 1)) {
      C[rownum, (k %% bin_partitions) + 1] <-
        X %>%
        filter(X1 >= bins[k + 1, 1], X1 < bins[k + 1, 3],
               X2 >= bins[k + 1, 2], X2 < bins[k + 1, 4]) %>% 
        nrow
      
      if(!((k + 1) %% bin_partitions)) {
        rownum <- rownum + 1 
      }
    }
    
    p <- C / nrow(X)
    
    # Empirical histogram
    histogram <- function(z) {
      if(z[1] < 0 | z[1] > 1 | z[2] < 0 | z[2] > 1) {
        return(0) 
      }
      
      rownum <- 1
      sum <- 0
      for(k in 0:(nrow(bins) - 1)) {
        sum <- sum + p[rownum, (k %% bin_partitions) + 1] / (1 / bin_partitions)^2 * 
          (z[1] >= bins[k + 1, 1] & z[1] < bins[k + 1, 3] & z[2] >= bins[k + 1, 2] & z[2] < bins[k + 1, 4])
        
        if(!((k + 1) %% bin_partitions)) {
          rownum <- rownum + 1 
        }
      }
      
      return(sum)
    }
    
    # Smoothed histogram
    smoothed_histogram <- function(z) {
      if(z[1] < 0 | z[1] > 1 | z[2] < 0 | z[2] > 1) {
        return(0) 
      }
      
      return((1 - delta) * histogram(z) + delta) 
    }
    
    bin_distributions <- rmultinom(1, K, (1 - delta) * c(p) / (1 / bin_partitions)^2 + delta)
    Z <- matrix(nrow = K, ncol = 2)
    X1 <- c()
    X2 <- c()
    for(k in 1:length(bin_distributions)) {
      X1 <- c(X1, runif(bin_distributions[k], bins[k,1], bins[k,3]))
      X2 <- c(X2, runif(bin_distributions[k], bins[k,2], bins[k,4]))
    }
    Z <- data.frame(X1, X2)
    temp <- cor.test(Z[,1], Z[,2], method = "spearman")
    if(temp$p.value < 0.05) {
      proportion <- proportion + 1 
    }
  }
  result[i] <- proportion / iterations_per_alpha
}

par(mfrow=c(1,1))
std_deviation <- sqrt(result * (1 - result) / iterations_per_alpha)
plot(x = log(alpha_levels) / log(10), y = result, type = "l", panel.first = grid(),
     xlab = "log10(alpha)", ylab = "Power",
     main = "Power of Spearman rank test for correlation")
lines(x = log(alpha_levels)/log(10), y = result + std_deviation, lty = 2, col = "red")
lines(x = log(alpha_levels)/log(10), y = result - std_deviation, lty = 2, col = "red")
legend("topleft", legend = c("Power", "Standard deviation"),
       col = c("black", "red"), lty = 1:2, cex = 1)
