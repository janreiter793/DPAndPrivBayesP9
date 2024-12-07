library(tidyverse)
library(MASS)
library(entropy)
library(magrittr)
library(VGAM)

# Parameters
alpha_1 <- 5 # Level of diff-privacy for Bayesian network-generation
alpha_2 <- 5 # Level of diff-privacy for noise injected in conditional props
k <- 3       # Degree of the Bayesian network
M <- 1000    # Number synthetic samples
n <- 1000    # Number of rows in original dataset
d <- 5       # Number of attributes (Will be binary)

# Correlation matrix
A <- matrix(
  c(1.0, 0.5, 0.4, 0.3, 0.2,
    0.5, 1.0, 0.5, 0.4, 0.3,
    0.4, 0.5, 1.0, 0.5, 0.4,
    0.3, 0.4, 0.5, 1.0, 0.5,
    0.2, 0.3, 0.4, 0.5, 1.0),
  nrow = d
)

# Test if matrix is symmetric and positive definite
all(A == t(A))
chol(A)

######################## Generate Data #########################################
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

######################## Generate a Bayesian Network ###########################
# Generate a Bayesian network differentially private using the exponential
# mechanism, and the mutual information as the scoring function.

N <- list()                              # Init the Bayesian network.
N[[1]] <- list(node = 1, parents = NULL) # First node has no ancestors.
                                         # NOTE: Program should try different
                                         # permutations, instead of first
                                         # column never having ancestors.

# Takes a node number and a set of potential parent-nodes, and then returns the
# estimated mutual information.
I <- function(node, parents) {
  subset <- X[, c(node, parents)]
  
  # Estimate the joint distribution by counting cases
  joint_distribution <-
    subset %>% 
    as.data.frame %>% 
    table %>% 
    prop.table
  
  # Estimate the entropy for each variable individually
  marginal_entropies <- 
    subset %>% 
    ncol %>% 
    numeric
  for(i in 1:ncol(subset)) {
    marginal_entropies[i] <-
      subset[,i] %>% 
      table %>% 
      prop.table %>% 
      entropy
  }
  
  # Joint entropy
  joint_entropy <- 
    joint_distribution %>% 
    entropy
  
  # I(X1, X2, ..., Xi) = H(X1) + H(X2) + ... + H(Xi) - H(X1, X2, ..., Xi)
  # needs proof
  return(sum(marginal_entropies) - joint_entropy)
}

# The sensitivity of I according to Lemma 4.1 of PrivBayes-article
Delta_I <- 1 / n * log2(n) + (n - 1) / n * log2(n / (n - 1))

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

######################## Estimate the noisy conditionals #######################
# Estimate the noisy conditionals based on the structure of the Bayesian net-
# work obtained in the former section. Injects noise using the Laplace-mechanism
conditional_props <- list()
conditional_props[[1]] <- table(X[,1])[2] / n

# Constructs arrays that resemble the conditional distributions of the attri-
# butes
for(i in 2:length(N)) {
  node <- N[[i]]$node
  parents <- N[[i]]$parents
  temp <- table(as.data.frame(X[, c(node, parents)]))
  
  total_marginals <- apply(temp, seq_along(dim(temp))[-1], sum)
  indices <- c(list(1), rep(TRUE, 
                            length(dim(temp)) - 1)
               )
  conditional_props[[i]] <- 
    do.call("[", c(list(temp), indices)) /
    total_marginals
}

# Inject Laplace-noise into the arrays
for(i in 1:length(N)) {
  if(is.numeric(conditional_props[[i]])) {
    conditional_props[[i]] <-
      conditional_props[[i]] +
      rlaplace(length(conditional_props[[i]]),
               location = 0,
               scale = 1 / alpha_2)
  } else {
    laplace <- 
      dim(conditional_props[[i]]) %>% 
      prod %>% 
      rlaplace(scale = 1 / alpha_2) %>% 
      array(dim = dim(conditional_props[[i]]))
    
    conditional_props[[i]] <- 
      conditional_props[[i]] + laplace
  }
  
  # If estimate is below zero set to zero, and if above one set to one
  conditional_props[[i]][conditional_props[[i]] < 0] <- 0
  conditional_props[[i]][conditional_props[[i]] > 1] <- 1
}

######################## Sample a synthetic dataset ############################
# This function samples a single sample from conditional_props
sampleObservation <- function() {
    synth_sample <- numeric(d)
    
    # The first node has no parents and is therefore sampled directly from its dis-
    # tribution
    synth_sample[1] <- rbinom(1, 1, prob = conditional_props[[1]])
    
    # For each node sample based on conditional_props
    for(i in 2:d) {
      parent_samples <- synth_sample[N[[i]]$parents] + 1
      prob <- conditional_props[[i]][parent_samples]
      synth_sample[i] <- rbinom(1, 1, prob)
    }
    
    return(synth_sample)
}

Z <- matrix(ncol = d, nrow = M)
for(i in 1:M) {
  Z[i,] <- sampleObservation()
}

# Compare covariance matrices
A
Z %>% data.frame %>% cor(method = "spearman")
X %>% data.frame %>% cor(method = "spearman") # Compare with estimated cor-
                                              # relation matrix of the original
                                              # dataset
