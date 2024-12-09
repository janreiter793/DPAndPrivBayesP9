library(dplyr)
library(tidyverse)
library(MASS)
library(entropy)
library(magrittr)
library(VGAM)
library(future.apply)
library(future)
library(progressr)

# Parameters
eps   <- c(0.1, 1, 10)
betas <- seq(from = 0.1, to = 0.9, length.out = 100)
iterations <- 500
k <- 3       # Degree of the Bayesian network
M <- 1000    # Number synthetic samples
n <- 1000    # Number of rows in original dataset
d <- 5       # Number of attributes (Will be binary)

# Correlation matrix
A <- matrix(
  c(1.0, 0.6, 0.5, 0.4, 0.3,
    0.6, 1.0, 0.6, 0.5, 0.4,
    0.5, 0.6, 1.0, 0.6, 0.5,
    0.4, 0.5, 0.6, 1.0, 0.6,
    0.3, 0.4, 0.5, 0.6, 1.0),
  nrow = d
)

# Test if matrix is symmetric and positive definite
all(A == t(A))
chol(A)

# The sensitivity of I according to Lemma 4.1 of PrivBayes-article
Delta_I <- 1 / n * log2(n) + (n - 1) / n * log2(n / (n - 1))

######################## Functions #############################################
# Takes a vector of elements v, and a size parameter s. Then it finds all
# subsets of v with at most s elements in them. Does not generate an empty
# subset, all nodes except first one must have parents
generate_subsets <- function(v, s) {
  subsets <- list()
  for(i in 1:min(s, length(v))) {
    # Concatenate list of subsets with length i
    subsets <- c(subsets, 
                 combn(v, i, simplify = FALSE)) 
  }
  return(subsets)
}

# Takes a vector of unnormalized density values and normalizes it
normalize <- function(vec) {
  return(vec / sum(vec)) 
}

# Takes a probability table representing a joint distribution. Marginalizes first
# variable out, and returns probability table representing joint distribution of
# the rest of the variables
marginalizeFirstVar <- function(tab) {
  temp <- array(dim = dim(tab)[-1])
  
  temp_entries <-
    temp %>% 
    dim %>% 
    lapply(seq_len) %>% 
    expand.grid
  
  for(i in 1:nrow(temp_entries)) {
    temp_coord <- temp_entries[i, ] %>% as.numeric
    temp <-
      do.call("[<-",
              c(list(temp), 
                as.list(temp_coord),
                list(  
                  do.call("[", c(list(tab), c(list(TRUE), temp_coord))) %>% 
                    sum
                ))
      )
  }
  
  return(temp)
}

# Takes a dataframe of bernoulli columns, then returns a matrix with 1's and 0's
# where a 1 represents an accepted null-hypothesis and 0 represents a rejected
# null-hypothesis
chisquare_test <- function(df) {
  temp <- matrix(nrow = d, ncol = d)
  for(i in 1:d) {
    for(j in 1:d) {
      tab <- table(df[,i], df[,j])
      temp[i, j] <- chisq.test(tab)$p.value >= 0.05
    }
  }
  
  return(temp + 0)
}

PrivBayes <- function(epsilon, beta) {
  alpha_1 <- beta * epsilon       # Level of diff-privacy for Bayesian network-generation
  alpha_2 <- (1 - beta) * epsilon # Level of diff-privacy for noise injected in conditional props
  
  # Generate data using copula-technique
  X <-
    mvrnorm(n = n,
            mu = numeric(d),
            Sigma = A) %>% 
    pnorm %>%
    {. >= 0.5} %>% 
    as.numeric %>% 
    matrix(nrow = n,
           ncol = d)
  
  # Takes a node number and a set of potential parent-nodes, and then returns the
  # estimated mutual information.
  I <- function(node, parents) {
    subset <- X[, c(node, parents)]
    
    # Estimate the joint distribution by counting cases
    joint_entropy <-
      subset %>% 
      as.data.frame %>% 
      table %>% 
      prop.table %>% 
      entropy(method = "ML")
    
    # Distribution of Y
    Y_entropy <-
      subset[, -1] %>% 
      as.data.frame %>% 
      table %>% 
      prop.table %>% 
      entropy(method = "ML")
    
    # Distribution of the attribute
    X_entropy <-
      subset[, 1] %>% 
      as.data.frame %>% 
      table %>% 
      prop.table %>% 
      entropy(method = "ML")
    
    # I(X, Y) = H(X) + H(Y) - H((X, Y))
    return(X_entropy + Y_entropy - joint_entropy)
  }
  
  # Generate a Bayesian network differentially private using the exponential
  # mechanism, and the mutual information as the scoring function.
  N <- list()                              # Init the Bayesian network.
  N[[1]] <- list(node = 1, parents = NULL) # First node has no ancestors.
  # NOTE: Program should try different
  # permutations, instead of first
  # column never having ancestors.
  
  # This part generates the Bayesian network
  gamma <- (d - 1) * Delta_I / alpha_1
  for(i in 2:d) {
    subsets <- generate_subsets(1:(i - 1), k)
    g_x <- subsets %>% length %>% numeric
    for(j in 1:length(subsets)) {
      g_x[j] <- exp(I(i, subsets[[j]]) / (2 * gamma)) # Should try to avoid the exponential
      # function using log
    }
    
    g_x %<>% normalize
    N[[i]] <- list(node = i, 
                   parents = subsets[[
                     sample(1:length(g_x), size = 1, prob = g_x)
                   ]])
  }

  # Estimate the noisy conditionals based on the structure of the Bayesian net-
  # work obtained in the former section. Injects noise using the Laplace-mechanism
  conditional_props <- list()
  conditional_props[[1]] <- table(X[,1]) / n
  
  # Constructs arrays that resemble the conditional distributions of the attri-
  # butes
  for(i in 2:d) {
    conditional_props[[i]] <- table(as.data.frame(X[, c(N[[i]]$node, 
                                                        N[[i]]$parents)])) / n
    marg_density <- marginalizeFirstVar(conditional_props[[i]])
    indices_1 <- c(list(1), rep(TRUE, length(dim(conditional_props[[i]])) - 1))
    indices_2 <- c(list(2), rep(TRUE, length(dim(conditional_props[[i]])) - 1))
    
    conditional_props[[i]] <-
      do.call("[<-",
              c(list(conditional_props[[i]]),
                indices_1,
                list(
                  do.call("[",
                          c(list(conditional_props[[i]]),
                            indices_1)) /
                    marg_density
                )))
    
    conditional_props[[i]] <-
      do.call("[<-",
              c(list(conditional_props[[i]]),
                indices_2,
                list(
                  do.call("[",
                          c(list(conditional_props[[i]]),
                            indices_2)) /
                    marg_density
                )))
  }
  
  # Inject Laplace-noise into the arrays
  for(i in 1:d) {
    if(is.numeric(conditional_props[[i]])) {
      conditional_props[[i]] <-
        conditional_props[[i]] +
        rlaplace(length(conditional_props[[i]]),
                 location = 0,
                 scale = d / n * 2 / alpha_2)
    } else {
      laplace <- 
        dim(conditional_props[[i]]) %>% 
        prod %>% 
        rlaplace(scale = d / n * 2^length(dim(conditional_props[[i]])) / 
                   alpha_2) %>% 
        array(dim = dim(conditional_props[[i]]))
      
      conditional_props[[i]] <- 
        conditional_props[[i]] + laplace
    }
    
    # If estimate is below zero set to zero, and normalize
    conditional_props[[i]][conditional_props[[i]] < 0] <- 0
    indices_1 <- c(list(1), rep(TRUE, length(dim(conditional_props[[i]])) - 1))
    indices_2 <- c(list(2), rep(TRUE, length(dim(conditional_props[[i]])) - 1))
    normalizing_const <- 1 / (do.call("[", c(list(conditional_props[[i]]), indices_1)) +
                                do.call("[", c(list(conditional_props[[i]]), indices_2)))
    conditional_props[[i]] <- do.call("[", c(list(conditional_props[[i]]), indices_2)) *
      normalizing_const
    conditional_props[[i]][is.infinite(conditional_props[[i]])] <- 1
    conditional_props[[i]][is.nan(conditional_props[[i]])] <- 0
  }
  
  # This function samples a single sample from conditional_props
  sampleObservation <- function(epsilon, beta) {
    
    synth_sample <- numeric(d)
    
    # The first node has no parents and is therefore sampled directly from its dis-
    # tribution
    synth_sample[1] <- rbinom(1, 1, prob = conditional_props[[1]])
    
    # For each node sample based on conditional_props
    for(i in 2:d) {
      parent_samples <- as.list(synth_sample[N[[i]]$parents] + 1)
      prob <- do.call("[", c(list(conditional_props[[i]]), c(parent_samples)))
      synth_sample[i] <- rbinom(1, 1, prob)
    }
    
    return(synth_sample)
  }
  
  Z <- matrix(ncol = d, nrow = M)
  for(i in 1:M) {
    Z[i,] <- sampleObservation()
  }
  
  return(chisquare_test(Z))
}

# Call this function for each thread
threads <- function(N) {
  res <- list()
  for(e in eps) {
    res[[as.character(e)]] <- list()
    for(beta in betas[(10 * N - 9):(10 * N)]) {
      z <- matrix(0, nrow = d, ncol = 5)
      #cat("Running", iterations, "iterations with eps =",
      #    e, "and beta =", beta, "\n")
      iter_star <- iterations
      for(iter in 1:iterations) {
        #cat(iter, "/", iterations, "\n")
        privbayes_res <- try({  
          PrivBayes(epsilon = e, beta = beta)
        }, silent = TRUE)
        
        if(inherits(privbayes_res, "try-error")) {
          iter_star <- iter_star - 1
          next
        }
        z <- z + privbayes_res
      }
      res[[as.character(e)]][[as.character(beta)]] <- z / iter_star
    }
  }
  
  return(res)
}

plan(multisession(workers = 10))
handlers("progress")  # Enable progress bar
res <- future_lapply(1:10, threads, future.seed = TRUE)
plan(sequential)

res <-
  mapply(c, 
         res[[1]], res[[2]], res[[3]], res[[4]], res[[5]],
         res[[6]], res[[7]], res[[8]], res[[9]], res[[10]],
         SIMPLIFY = FALSE)


# Extract the chi-square power from the result
plotPowers <- function(privacy) {
  privacy <- as.character(privacy)
  attr_1_and_2 <- betas %>% length %>% numeric
  attr_1_and_3 <- betas %>% length %>% numeric
  attr_1_and_4 <- betas %>% length %>% numeric
  attr_1_and_5 <- betas %>% length %>% numeric
  for(i in 1:length(betas)) {
    name <- names(res[[privacy]])[i]
    attr_1_and_2[i] <- res[[privacy]][[name]][2,1]
    attr_1_and_3[i] <- res[[privacy]][[name]][3,1] 
    attr_1_and_4[i] <- res[[privacy]][[name]][4,1]
    attr_1_and_5[i] <- res[[privacy]][[name]][5,1]
  }
  matplot(x = betas, cbind(1 - attr_1_and_2,
                           1 - attr_1_and_3,
                           1 - attr_1_and_4,
                           1 - attr_1_and_5), 
          type = "l", lty = 1,
          ylab = "Power", xlab = "beta",
          main = paste0("Estimated power of Chi-squared test on PrivBayes data (eps = ",
                        privacy, ")"),
          panel.first = grid(),
          col = c("black", "red", "blue", "green"))
  legend("bottomright", legend = c("Z1 and Z2",
                                   "Z1 and Z3",
                                   "Z1 and Z4",
                                   "Z1 and Z5"),
         col = c("black", "red", "blue", "green"),
         lty = 1)
}

png(file = "~/Sem9/Projekt/power1.png", width = 850, height = 550)
  plotPowers(eps[1])
dev.off()
