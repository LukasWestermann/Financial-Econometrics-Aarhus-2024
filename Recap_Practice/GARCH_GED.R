library(Rsolnp) #For (constrained) nonlinear optimisation
library(fGarch) # For GED 
source("utils.R") # For simualting GARCH process for testing


# GED Log-Likelihood Function for GARCH(1,1)
garch11_ged_loglik <- function(vPar, dY) {
  # Parameters
  dOmega <- vPar[1]
  dAlpha <- vPar[2]
  dBeta  <- vPar[3]
  dNu    <- vPar[4]
  
  iT = length(dY)
  
  ## initialize the variance
  vSigma2 = numeric(iT)
  
  ## set the variance at time t = 1 at its unconditional value
  vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta)
  
  
  # Calculate lambda (scaling factor for GED)
  lambda <- (2^(-2/dNu) * gamma(1/dNu) / gamma(3/dNu))^(0.5)
  
  # Iteratively compute sigma^2 (conditional variance)
  for (t in 2:iT) {
    vSigma2[t] <- dOmega + dAlpha * dY[t-1]^2 + dBeta * vSigma2[t-1]
  
  }
  
  # Log-likelihood ## see slide 54 lecture 4 for general likelihood form of GARCH models
  dLLK <-  T * log(dNu / (lambda * gamma(1/dNu))) - 
    sum((abs(dY / sqrt(vSigma2)))^dNu / lambda^dNu) - 0.5 * sum(log(vSigma2))
  
  lOut = list()
  lOut[["dLLK"]] = dLLK
  lOut[["vSigma2"]] = vSigma2
  return(lOut)  # Return log likelihood and estimated varainces
}

# Estimation Function (Optimization)
estimate_garch11_ged <- function(dY) {
  # Initial values for parameters
  vPar <- c(dOmega = 0.01, dAlpha = 0.05, dBeta = 0.9, dNu = 5)
  
  

  
  optimizer = solnp(vPar, fun = function(vPar, dY) {
    
    
    dnLLK = -as.numeric(garch11_ged_loglik(vPar, dY)$dLLK)
    
    return(dnLLK)
    
  }, ineqfun = function(vPar, ...) {
    vPar[2] + vPar[3]
  }, ineqLB = 1e-4, ineqUB = 0.999, # to ensure stationarity
  LB = c(1e-4, 1e-4, 1e-4, 2.01), UB = c(0.999, 0.999, 0.999, 30),
  dY = dY)
  
  # extract optimal parameters
  vPar <- optimizer$pars
  
  # Compute conditional variances with optimal paramas obtained via MLE
  vSigma2 = garch11_ged_loglik(vPar, dY)$vSigma2
  
  # extract likelihood at its maximum
  dLLK = - optimizer$values[length(optimizer$values)]
  
  #Compute average BIC
  iT = length(dY)
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT
  
  #return list with estimated conditional variances, parameters, likelihood value and BIC
  lOut = list(vSigma2 = vSigma2, vPar = vPar, dLLK = dLLK, BIC = BIC)
  
  return(lOut)
}

# h-step Ahead Forecast for GARCH(1,1) as in lecture 4
forecast_garch11 <- function(iH, Sigma2_t, return_t, dOmega, dAlpha, dBeta) {
  
  vSigma2_pred = numeric(iH)
  
  vSigma2_pred[1] = dOmega + dAlpha * return_t^2 + dBeta * Sigma2_t #Compute forecast at h=1
  
  if (iH > 1) {
    
    dSigma2_unconditional = dOmega/(1.0 - dAlpha - dBeta) #Compute unconditional variance
    
    #main loop for calcualting h>=2 forecasts
    for (h in 2:iH) {
      
      vSigma2_pred[h] = dSigma2_unconditional + (dAlpha + dBeta)^(h-1) * ( vSigma2_pred[1] - dSigma2_unconditional)
    }
  }
  return(vSigma2_pred)
}





# Example Usage
set.seed(123)
data <- sim_garch11(T = 1000, dOmega = 0.02, dAlpha = 0.05, dBeta = 0.9, dNu = 2)
fit <- estimate_garch11_ged(data$dY)

print(fit$vPar)  # Estimated Parameters
print(fit$dLLK) # Log-Likelihood Value
print(fit$BIC)  #BIC


plot(fit$vSigma2, type = 'l', main = 'Estimated vs Actual Variance', ylab = 'Variance', ylim = c(0.2,1))
lines(data$variance, col = 'red')  # Plot actual variance for comparison
legend('topright', legend = c('Estimated', 'Actual'), col = c('black', 'red'), lty = 1)
#something with the scaling does not work here with the ged function and my implementation
#probably different parameteriaztions

h_forecast <- forecast_garch11(iH = 20, tail(fit$vSigma2,1), tail(data$dY,1), fit$vPar[1], fit$vPar[2], fit$vPar[3])  
print(h_forecast)
plot(h_forecast, type = 'l', main = 'h-step Ahead Forecast', ylab = 'Forecasted Variance', cex.main = 1.5, cex.lab = 1.2)


#With package
# spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), variance.model = 
#               list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "ged")
# fit_pac = ugarchfit(spec, data$dY)
# coef(fit_pac)
# ugarchforecast(fit_pac, n.ahead = 20)
# infocriteria(fit_pac)[2]

