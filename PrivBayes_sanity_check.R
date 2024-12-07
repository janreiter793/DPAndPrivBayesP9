library(tidyverse)
library(MASS)
library(entropy)
library(magrittr)
library(VGAM)
set.seed(100)

# Parameters
epsilon <- c(0.1, 1, 10, 100)
beta <- 0.5
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
    entropy
  
  # Distribution of Y
  Y_entropy <-
    subset[, -1] %>% 
    as.data.frame %>% 
    table %>% 
    prop.table %>% 
    entropy
  
  # Distribution of the attribute
  X_entropy <-
    subset[, 1] %>% 
    as.data.frame %>% 
    table %>% 
    prop.table %>% 
    entropy
  
  # I(X, Y) = H(X) + H(Y) - H((X, Y))
  return(X_entropy + Y_entropy - joint_entropy)
}

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

# This function samples a single sample from conditional_props
sampleObservation <- function() {
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

# Takes a privacy budget [eps], runs PrivBayes and returns a list of two fitted
# Bayesian networks. N_priv is the differentially private generated one, and
# N_greedy is the greedily picked one. They should somewhat resemble each other
# for high eps, and differ a bit for low eps
PrivBayes_sanity <- function(eps) {
  # Distribute the privacy budget
  alpha_1 <- beta * eps
  alpha_2 <- (1 - beta) * eps
  
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
  
  # Generate a Bayesian network differentially private using the exponential
  # mechanism, and the mutual information as the scoring function.
  
  N <- list()                              # Init the Bayesian network.
  N[[1]] <- list(node = 1, parents = NULL) # First node has no ancestors.
  
  # This part generates the Bayesian network
  gamma <- (d - 1) * Delta_I / alpha_1
  for(i in 2:d) {
    subsets <- generate_subsets(1:(i - 1), k)
    g_x <- subsets %>% length %>% numeric
    for(j in 1:length(subsets)) {
      g_x[j] <- exp(I(i, subsets[[j]]) / (2 * gamma))
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
  conditional_props[[1]] <- table(X[,1])[2] / n
  
  # Constructs arrays that resemble the conditional distributions of the attri-
  # butes
  for(i in 2:length(N)) {
    temp <- table(as.data.frame(X[, c(N[[i]]$node, 
                                      N[[i]]$parents)])) / n
    indices <- c(list(2), rep(TRUE, length(dim(temp)) - 1))
    conditional_props[[i]] <-
      do.call("[", c(list(temp), indices))
  }
  
  # Inject Laplace-noise into the arrays
  for(i in 1:length(N)) {
    if(is.numeric(conditional_props[[i]])) {
      conditional_props[[i]] <-
        conditional_props[[i]] +
        rlaplace(length(conditional_props[[i]]),
                 location = 0,
                 scale = d / n * 2 * length(conditional_props[[i]]) / alpha_2)
    } else {
      laplace <- 
        dim(conditional_props[[i]]) %>% 
        prod %>% 
        rlaplace(scale = d / n * 2^(1 + length(dim(conditional_props[[i]]))) / 
                   alpha_2) %>% 
        array(dim = dim(conditional_props[[i]]))
      
      conditional_props[[i]] <- 
        conditional_props[[i]] + laplace
    }
    
    # If estimate is below zero set to zero, and if above one set to one
    conditional_props[[i]][conditional_props[[i]] < 0] <- 0
    conditional_props[[i]][conditional_props[[i]] > 1] <- 1
  }
  
  Z <- matrix(ncol = d, nrow = M)
  for(i in 1:M) {
    Z[i,] <- sampleObservation()
  }
  
  # Fit a Bayesian network using a greedy algorithm
  N_greedy <- list()                              # Init the Bayesian network.
  N_greedy[[1]] <- list(node = 1, parents = NULL) # First node has no ancestors.
  
  for(i in 2:d) {
    subsets <- generate_subsets(1:(i - 1), k)
    mutual_informations <- subsets %>% length %>% numeric
    for(j in 1:length(subsets)) {
      mutual_informations[j] <- I(i, subsets[[j]])
    }
    
    N_greedy[[i]] <- list(node = i, 
                          parents = subsets[[
                            which.max(mutual_informations)
                          ]])
  }
  
  return(list(N_priv = N, N_greedy = N_greedy))
}

res_0.1  <- PrivBayes_sanity(epsilon[1])
res_1    <- PrivBayes_sanity(epsilon[2])
res_10   <- PrivBayes_sanity(epsilon[3])
res_100  <- PrivBayes_sanity(epsilon[4])
