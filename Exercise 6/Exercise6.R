# Load necessary packages
library(rugarch)
library(rmgarch)
source("help_fn.R")

#Function estimate a Gaussian DDC model using MLE (two step approch)

estimate_gaussian_dcc = function(returns){
  
  #Step 1: Fit univarative GRCH(1,1) models using ugarch to estimate conditional variance
  n_assets = ncol(returns)
  garch_specs = vector(mode = "list", length = n_assets)
  garch_fits = vector(mode = "list", length = n_assets)
  
  for(i in 1:n_assets){
    garch_specs[[i]] = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = 
                                    list(model = "sGARCH"), distribution.model = "norm")
    garch_fits[[i]] = ugarchfit(spec = garch_specs[[i]], data = returns[,i])
  }
  
  #Step 2 extract standardized residuals and estimated variances
  residuals = sapply(garch_fits, function(fit) residuals(fit, standardize = TRUE))
  sigma_t = sapply(garch_fits, function(fit) sigma(fit))

  #Step 3 Estimate DCC parameters using MLE
  initial_params = c(0.05,0.93)
  mle_fit = optim(par = initial_params, fn = log_likelihood_dcc_normal, residuals = residuals, method = "L-BFGS-B",
                  lower = c(0.01,0.01), upper = c(1 - 0.01, 1 - 0.01)) # to ensure p.d. of Qt matrix
  alpha = mle_fit$par[1]
  beta = mle_fit$par[2]
  log_lik = mle_fit$value

  #Compute initial correlation matrix
  Q_bar = cor(residuals)
  Q_t = Q_bar

  #Time varying correlation matrix
  T = nrow(returns)
  R_t_list = vector(mode = "list", length = T)

  #caluclate correlation matrix with optimal parameters obatained via MLE

  for(t in 2:T){
    Q_t = (1-alpha -beta) * Q_bar + alpha * (residuals[t-1, ] %*% t(residuals[t-1, ])) + beta * Q_t
    D_t = diag(sqrt(diag(Q_t)))
    R_t = solve(D_t) %*% Q_t %*% solve(D_t)
    R_t_list[[t]] = R_t
  }

  return(list(correlations = R_t_list, variances = sigma_t,log_lik = log_lik, parameters = list(alpha = alpha, beta = beta)))

}

#Function to estimate a multivartaive student's t DCC model

estimate_student_t_dcc = function(returns){
  
  #Step 1: Fit univarative GRCH(1,1) models using ugarch to estimate conditional variance, same as in Gaussain case
  n_assets = ncol(returns)
  garch_specs = vector(mode = "list", length = n_assets)
  garch_fits = vector(mode = "list", length = n_assets)
  
  for(i in 1:n_assets){
    garch_specs[[i]] = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = 
                                    list(model = "sGARCH"), distribution.model = "norm")
    garch_fits[[i]] = ugarchfit(spec = garch_specs[[i]], data = returns[,i])
  }
  
  #Step 2 extract standardized residuals and estimated variances
  residuals = sapply(garch_fits, function(fit) residuals(fit, standardize = TRUE))
  sigma_t = sapply(garch_fits, function(fit) sigma(fit))
  
  #Step 3 Estimate DCC parameters using MLE, add degress of freedom paramter
  initial_params = c(0.05,0.9,10)
  mle_fit = optim(par = initial_params, fn = log_likelihood_dcc_t, residuals = residuals, method = "L-BFGS-B",
                  lower = c(0.01,0.01,3), upper = c(1 - 0.01, 1 - 0.01,100)) # to ensure p.d. of Qt matrix
  alpha = mle_fit$par[1]
  beta = mle_fit$par[2]
  nu = mle_fit$par[3]
  
  #Compute initial correlation matrix
  Q_bar = cor(residuals)
  Q_t = Q_bar
  
  #Time varying correlation matrix
  T = nrow(returns)
  R_t_list = vector(mode = "list", length = T)
  
  #caluclate correlation matrix with optimal parameters obatained via MLE
  
  for(t in 2:T){
    Q_t = (1-alpha -beta) * Q_bar + alpha * (residuals[t-1, ] %*% t(residuals[t-1, ])) + beta * Q_t
    D_t = diag(sqrt(diag(Q_t)))
    R_t = solve(D_t) %*% Q_t %*% solve(D_t)
    R_t_list[[t]] = R_t
  }
  
  return(list(correlations = R_t_list, variances = sigma_t, parameters = list(alpha = alpha, beta = beta, nu = nu)))
  
}





# Example usage for Gaussian and Student's t DCC estimation
# Load data (e.g., dji30ret dataset from rugarch package)
data(dji30ret)
returns <- dji30ret[, 1:2]  # Select two assets for bivariate estimation

results = estimate_gaussian_dcc(returns)
results$log_lik

# Specify a DCC-GARCH model with Gaussian distribution
spec_garch <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                         variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                         distribution.model = "norm")
spec_dcc_gaussian <- dccspec(uspec = multispec(replicate(2, spec_garch)),
                             dccOrder = c(1, 1),
                             distribution = "mvnorm")

# Fit the Gaussian DCC model
dcc_fit_gaussian <- dccfit(spec = spec_dcc_gaussian, data = returns)

# Specify a DCC-GARCH model with Student's t distribution
spec_garch_t <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                           variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                           distribution.model = "std")
spec_dcc_t <- dccspec(uspec = multispec(replicate(2, spec_garch_t)),
                      dccOrder = c(1, 1),
                      distribution = "mvt")

# Fit the Student's t DCC model
dcc_fit_t <- dccfit(spec = spec_dcc_t, data = returns)

# Extract filtered correlations for both models
correlations_gaussian <- rcor(dcc_fit_gaussian)
correlations_t <- rcor(dcc_fit_t)

# Compare filtered correlations (example for plotting correlations over time)
par(mfrow = c(2, 1))
plot(ts(correlations_gaussian[1, 2, ]), main = "Filtered Correlations (Gaussian DCC)", ylab = "Correlation", xlab = "Time")
plot(ts(correlations_t[1, 2, ]), main = "Filtered Correlations (Student's t DCC)", ylab = "Correlation", xlab = "Time")

# Compare models using BIC
bic_gaussian <- infocriteria(dcc_fit_gaussian)[2]
bic_t <- infocriteria(dcc_fit_t)[2]

cat("Gaussian DCC BIC:", bic_gaussian, "\n")
cat("Student's t DCC BIC:", bic_t, "\n")

# Determine which model is selected based on BIC
if (bic_gaussian < bic_t) {
  cat("The Gaussian DCC model is preferred based on BIC.\n")
} else {
  cat("The Student's t DCC model is preferred based on BIC.\n")
}




