# Simulate Data for Testing (GARCH(1,1) Process)
sim_garch11 <- function(T, dOmega, dAlpha, dBeta, dNu) {
  
  z <- rged(T, mean = 0, sd = 1, nu = dNu)  # GED innovationvs
  lambda <- (2^(-2/dNu) * gamma(1/dNu) / gamma(3/dNu))^(0.5)
  z <- z * lambda
  
  sigma2 <- rep(0, T)
  dY <- rep(0, T)
  sigma2[1] <- dOmega / (1 - dAlpha - dBeta)
  
  for (t in 2:T) {
    sigma2[t] <- dOmega + dAlpha * dY[t-1]^2 + dBeta * sigma2[t-1]
    dY[t] <- z[t] * sqrt(sigma2[t])
  }
  return(list(dY = dY, variance = sigma2))
}