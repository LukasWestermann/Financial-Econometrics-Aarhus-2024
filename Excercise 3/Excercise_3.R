####Computational part####
##a##
# Define the Kalman Filter function
kalman_filter = function(y,Z,T,D,H,Q,R,a0,P0){
  #number of observations
  n = length(y)
  
  #Initalize matrices to store results
  at = array(NA,dim = c(length(a0),n)) #State estimates
  Pt = array(NA,dim = c(length(a0), length(a0),n)) #State covariances
  vt = numeric(n) #Prediction errors
  Ft = numeric(n) #Variance of the predicted errors
  
  #Inital states
  at[,1] = a0
  Pt[, , 1] = P0
   
  #Inital log likelihood
  ll = 0
  
  #Kalman Filter recursions
  for(t in 1:n){
    #Prediction for observation
    vt[t] = y[t] - Z %*% at[,t]
    Ft[t] = Z %*% P[,,t] %*% t(Z)
  }
  
}

# Define quasi-likelihood function for optimization
quasi_likelihood <- function(params, y) {
  # Transform log-scale parameters to original scale
  sigma <- exp(params[1])
  phi <- params[2]
  sigma_eta <- exp(params[3])
  
  # Define matrices based on parameters
  Z <- matrix(1, nrow = 1)
  T <- matrix(phi, nrow = 1)
  D <- matrix(0, nrow = 1)
  H <- matrix(1, nrow = 1)
  Q <- matrix(sigma_eta^2, nrow = 1)
  R <- sigma^2
  
  # Initial state estimates
  a0 <- 0
  P0 <- matrix(sigma_eta^2 / (1 - phi^2), nrow = 1)
  
  # Run Kalman filter
  kf_results <- kalman_filter_with_log_likelihood(y, Z, T, D, H, Q, R, a0, P0)
  
  # Return negative log-likelihood
  return(-kf_results$log_likelihood)
}

# Simulate data for testing (from Task 2-2 and 2-3)
set.seed(123)
sigma <- 1
phi <- 0.9
sigma_eta <- sqrt(0.25)
T <- 1000
y <- numeric(T)
h <- numeric(T)

h[1] <- rnorm(1, mean = 0, sd = sigma_eta / sqrt(1 - phi^2))
for (t in 2:T) {
  h[t] <- phi * h[t-1] + rnorm(1, mean = 0, sd = sigma_eta)
}
y <- exp(h / 2) * rnorm(T, mean = 0, sd = sigma)

# Define initial parameter guesses for log(sigma), phi, and log(sigma_eta)
initial_params <- c(log(1), 0.5, log(sqrt(0.25)))

# Optimize the quasi-likelihood function using optim
fit <- optim(
  par = initial_params,
  fn = quasi_likelihood,
  y = y,
  method = "BFGS",
  control = list(maxit = 1000, reltol = 1e-8)
)

# Transform estimated parameters back to original scale
estimated_sigma <- exp(fit$par[1])
estimated_phi <- fit$par[2]
estimated_sigma_eta <- exp(fit$par[3])

# Print results
estimated_params <- c(sigma = estimated_sigma, phi = estimated_phi, sigma_eta = estimated_sigma_eta)
print(estimated_params)






kalman_filter <- function(y, Z, T, D, H, Q, R, a0, P0) {
  
  
  
  
  # Kalman Filter recursions
  for (t in 1:n) {
    # Prediction for observation
    vt[t] <- y[t] - Z %*% at[, t]
    Ft[t] <- Z %*% Pt[, , t] %*% t(Z) + R  # Ft = ZPtZ' + σ²εDD'
    
    # Kalman Gain
    Kt <- T %*% Pt[, , t] %*% t(Z) %*% solve(Ft[t])  # Kt = TPtZ'Ft⁻¹
    
    # Update state estimate and covariance
    at[, t + 1] <- T %*% at[, t] + Kt * vt[t]
    Lt <- T - Kt %*% Z  # Lt = T - KtZ
    Pt[, , t + 1] <- T %*% Pt[, , t] %*% t(Lt) + H %*% Q %*% t(H)  # Pt+1 = TPtLt' + HQH'
    
    ll <- ll - 0.5 * (log(Ft[t]) + (vt[t]^2) / Ft[t] + log(2 * pi))
  }
  
  # Return results
  return(list(at = at, Pt = Pt, vt = vt, Ft = Ft, ll = ll))
}
