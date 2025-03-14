
#    EM ALGORITHM
#---------------------

# k: Number of classes (j)
# D: Dimensions 
# N: Number of observations (i)
# max_iter: maximum number of iterations
# tol: convergence tolerance

em_gaussian_mixture <- function(X, K, max_iter = 100, tol = 1e-6) {
  
N <- nrow(X)  # Number of data points
D <- ncol(X)  # Number of dimensions

# Set initial parameters
pi <-  rep(1/K, K)    # Mixing coefficients

# mu
idx <- sample(1:N, K)
mu <- X[idx, , drop = FALSE]

# Covariance Matrix
sigma <- list()
for (j in 1:K) {
  sigma[[j]] <- diag(1, D)
}

# Initialize log-likelihood
prev_log_lik <- -Inf
log_likelihood <- c()

# Main EM loop
for (iter in 1:max_iter) {
  cat("Iteration:", iter, "\n")
  
#----------------------------
#   EXPECTATION (E-STEP)
#----------------------------
cat("  E-step: Computing responsibilities\n")

# Initialize prob matrix (N x K)
gamma <- matrix(0, nrow = N, ncol = K)


# Compute responsibilities for each data point and cluster
for (i in 1:N) {
  for (j in 1:K) {
    # Calculate multivariate Gaussian density
    gamma[i, j] <- pi[j] * dmvnorm(X[i, ], mean = mu[j, ], sigma = sigma[[j]])
  }
  # Normalize responsibilities for data point i
  if (sum(gamma[i, ]) > 0) {
    gamma[i, ] <- gamma[i, ] / sum(gamma[i, ])
  }
}

#----------------------------
#   MAXIMIZATION (M-STEP)
#----------------------------
cat("  M-step: Updating parameters\n")

# Calculate effective number of points in each class
N_k <- colSums(gamma) # weighted sum (prob of being in each class)

# Update mixing coefficients
pi <- N_k / N

# Update means
for (j in 1:K) {
  if (N_k[j] > 0) {
    # Weighted average of data points
    mu[j, ] <- t(gamma[, j]) %*% X / N_k[j]
  }
}

# Update covariance matrices
for (j in 1:K) {
  sigma[[j]] <- matrix(0, nrow = D, ncol = D)
  if (N_k[j] > 0) {
    for (i in 1:N) {
      diff <- X[i, ] - mu[j, ]
      sigma[[j]] <- sigma[[j]] + gamma[i, j] * (diff %*% t(diff))
    }
    sigma[[j]] <- sigma[[j]] / N_k[j]
    
    # Add small regularization to prevent singularity
    sigma[[j]] <- sigma[[j]] + diag(1e-6, D)
    }
  }

#--------------------------
# Check convergence
#--------------------------
# Compute log-likelihood
current_log_lik <- 0
for (i in 1:N) {
  # Sum weighted probability over all clusters
  p_xi <- 0
  for (j in 1:K) {
    p_xi <- p_xi + pi[j] * dmvnorm(X[i, ], mean = mu[j, ], sigma = sigma[[j]])
  }
  current_log_lik <- current_log_lik + log(p_xi)
}

log_likelihood <- c(log_likelihood, current_log_lik)

# Check for convergence
improvement <- abs(current_log_lik - prev_log_lik)
cat("  Log-likelihood:", current_log_lik, "Improvement:", improvement, "\n")

if (iter > 1 && improvement < tol) {
  cat("Converged after", iter, "iterations\n")
  break
}

prev_log_lik <- current_log_lik
}

# Return results
return(list(
  pi = pi,                   # Mixing coefficients
  mu = mu,                   # Means
  sigma = sigma,             # Covariance matrices
  gamma = gamma,             # Responsibilities
  log_likelihood = log_likelihood,  # Log-likelihood history
  K = K,                     # Number of clusters
  D = D                      # Number of dimensions
))
}

#----------------
#   PREDICTION
#----------------
# Function to predict cluster assignments
predict_cluster <- function(model, X_new) {
  N_new <- nrow(X_new)
  K <- model$K
  
  # Calculate probabilities for each point and cluster
  probs <- matrix(0, nrow = N_new, ncol = K)
  for (i in 1:N_new) {
    for (j in 1:K) {
      probs[i, j] <- model$pi[j] * dmvnorm(X_new[i, ], model$mu[j, ], model$sigma[[j]])
    }
  }
  
  # Return most likely cluster
  return(max.col(probs))
}

#------------------------
#   plot   convergence
#------------------------
# Function to plot convergence
plot_convergence <- function(model) {
  plot(model$log_likelihood, type = "o", 
       xlab = "Iteration", ylab = "Log-likelihood",
       main = "EM Algorithm Convergence")
}


