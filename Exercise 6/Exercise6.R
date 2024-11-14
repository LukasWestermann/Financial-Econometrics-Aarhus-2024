####Computational part####

# Load required libraries
library(rugarch)

# Gaussian DCC model function
dcc_gaussian <- function(data, garch_spec) {
  # Step 1: Estimate individual GARCH models
  n_assets <- ncol(data)
  garch_models <- lapply(1:n_assets, function(i) ugarchfit(spec = garch_spec, data = data[, i]))
  
  # Step 2: Extract standardized residuals and volatilities
  std_residuals <- sapply(garch_models, function(model) residuals(model, standardize = TRUE))
  volatilities <- sapply(garch_models, function(model) sigma(model))
  
  # Step 3: Initialize parameters for DCC
  n <- nrow(data)
  dcc_params <- c(0.05, 0.94)  # Initial guesses for DCC parameters (a and b)
  
  # Define the DCC log-likelihood function
  dcc_loglik <- function(params) {
    a <- params[1]
    b <- params[2]
    Q_bar <- cov(std_residuals)  # Initial unconditional correlation matrix
    Q_t <- array(0, dim = c(n_assets, n_assets, n))
    R_t <- array(0, dim = c(n_assets, n_assets, n))
    log_likelihood <- 0
    
    # DCC recursive updates
    for (t in 2:n) {
      Q_t[, , t] <- (1 - a - b) * Q_bar + a * (std_residuals[t - 1, ] %*% t(std_residuals[t - 1, ])) + b * Q_t[, , t - 1]
      R_t[, , t] <- diag(1 / sqrt(diag(Q_t[, , t]))) %*% Q_t[, , t] %*% diag(1 / sqrt(diag(Q_t[, , t])))
      log_likelihood <- log_likelihood - 0.5 * (log(det(R_t[, , t])) + t(std_residuals[t, ]) %*% solve(R_t[, , t]) %*% std_residuals[t, ])
    }
    
    return(-log_likelihood)  # Return negative log-likelihood
  }
  
  # Optimize DCC parameters
  optim_result <- optim(dcc_params, dcc_loglik, method = "L-BFGS-B", lower = c(0.0001, 0.0001), upper = c(0.99, 0.99))
  dcc_estimates <- optim_result$par
  
  # Filtered correlations and covariances
  list(log_likelihood = -optim_result$value, filtered_correlations = R_t, estimated_params = dcc_estimates)
}

# Example usage with simulated data (data should be your actual asset return matrix)
garch_spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0)), distribution.model = "norm")
data <- matrix(rnorm(1000), ncol = 2)  # Replace with actual data
result_gaussian_dcc <- dcc_gaussian(data, garch_spec)






# Multivariate Student’s t DCC model function
dcc_student_t <- function(data, garch_spec, df_init = 10) {
  n_assets <- ncol(data)
  garch_models <- lapply(1:n_assets, function(i) ugarchfit(spec = garch_spec, data = data[, i]))
  
  std_residuals <- sapply(garch_models, function(model) residuals(model, standardize = TRUE))
  volatilities <- sapply(garch_models, function(model) sigma(model))
  
  Q_bar <- cov(std_residuals)
  dcc_params <- c(0.05, 0.94)
  df <- df_init
  
  # Define the Student’s t DCC log-likelihood
  dcc_loglik_student <- function(params) {
    a <- params[1]
    b <- params[2]
    df <- params[3]
    Q_t <- array(0, dim = c(n_assets, n_assets, n))
    R_t <- array(0, dim = c(n_assets, n_assets, n))
    log_likelihood <- 0
    
    for (t in 2:n) {
      Q_t[, , t] <- (1 - a - b) * Q_bar + a * (std_residuals[t - 1, ] %*% t(std_residuals[t - 1, ])) + b * Q_t[, , t - 1]
      R_t[, , t] <- diag(1 / sqrt(diag(Q_t[, , t]))) %*% Q_t[, , t] %*% diag(1 / sqrt(diag(Q_t[, , t])))
      term1 <- lgamma((df + n_assets) / 2) - lgamma(df / 2)
      term2 <- -0.5 * log(det(R_t[, , t])) - n_assets / 2 * log(df * pi)
      term3 <- -(df + n_assets) / 2 * log(1 + (t(std_residuals[t, ]) %*% solve(R_t[, , t]) %*% std_residuals[t, ]) / df)
      log_likelihood <- log_likelihood + (term1 + term2 + term3)
    }
    
    return(-log_likelihood)
  }
  
  # Optimize DCC and degrees of freedom parameters
  optim_result <- optim(c(dcc_params, df), dcc_loglik_student, method = "L-BFGS-B", lower = c(0.0001, 0.0001, 2.1), upper = c(0.99, 0.99, 30))
  dcc_estimates <- optim_result$par
  
  list(log_likelihood = -optim_result$value, filtered_correlations = R_t, estimated_params = dcc_estimates)
}

# Example usage
result_student_dcc <- dcc_student_t(data, garch_spec)
