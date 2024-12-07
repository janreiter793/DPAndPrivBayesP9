library(magrittr)
library(tidyverse)

## Parameters
N              <- 10000
sigma          <- 0.1
bin_partitions <- 4
delta          <- 0.01
K              <- 5000

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
for(i in 1:(length(biny) - 1)) {
  for(j in 1:(length(binx) - 1)) {
    X1[num] <- binx[j]
    Y1[num] <- biny[i]
    X2[num] <- binx[j + 1]
    Y2[num] <- biny[i + 1]
    num <- num + 1
  }
}

bins <- data.frame(X1, Y1, X2, Y2) 

## Generate the histogram
C <- matrix(nrow = bin_partitions, ncol = bin_partitions)

rownum <- 1
for(i in 0:(nrow(bins) - 1)) {
  C[rownum, (i %% bin_partitions) + 1] <-
    X %>%
    filter(X1 >= bins[i + 1, 1], X1 < bins[i + 1, 3],
           X2 >= bins[i + 1, 2], X2 < bins[i + 1, 4]) %>% 
    nrow
  
  if(!((i + 1) %% bin_partitions)) {
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
  for(i in 0:(nrow(bins) - 1)) {
    sum <- sum + p[rownum, (i %% bin_partitions) + 1] / (1 / bin_partitions)^2 * 
      (z[1] >= bins[i + 1, 1] & z[1] < bins[i + 1, 3] & z[2] >= bins[i + 1, 2] & z[2] < bins[i + 1, 4])
  
    if(!((i + 1) %% bin_partitions)) {
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
for(i in 1:length(bin_distributions)) {
  X1 <- c(X1, runif(bin_distributions[i], bins[i,1], bins[i,3]))
  X2 <- c(X2, runif(bin_distributions[i], bins[i,2], bins[i,4]))
}
Z <- data.frame(X1, X2)

alpha <- K * log((1 - delta) * bin_partitions^2 / (N * delta) + 1)

par(mfrow=c(1,2))
plot(X, panel.first = grid(),
     main = "Original dataset X")
plot(Z, panel.first = grid(),
     xlab = "X1", ylab = "X2", main = "Sampled dataset Z")

