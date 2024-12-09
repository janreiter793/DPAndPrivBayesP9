library(tidyverse)
library(MASS)
library(entropy)
library(magrittr)
library(VGAM)

# Parameters
eps <- 0.1
beta <- 0.5
alpha_1 <- beta * eps # Level of diff-privacy for Bayesian network-generation
alpha_2 <- (1 - beta) * eps # Level of diff-privacy for noise injected in conditional props
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
    g_x[j] <- exp(I(i, subsets[[j]]) / (2 * gamma)) # Should try to avoid the exponential
                                                    # function using log
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
conditional_props[[1]] <- table(X[,1]) / n

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

######################## Sample a synthetic dataset ############################
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

Z <- matrix(ncol = d, nrow = M)
for(i in 1:M) {
  Z[i,] <- sampleObservation()
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

chisquare_test(Z)
chisquare_test(X)

######################## Greedy Bayes ##########################################
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
