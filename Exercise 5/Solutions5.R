

# Ex2

tGAS_LLK <- function(vY, dOmega, dAlpha, dBeta, dNu) {

  iT = length(vY)

  vPhi = numeric(iT)
  vPhi_tilde = numeric(iT)

  #initialize at its unconditional value
  vPhi_tilde[1] = dOmega/(1.0 - dBeta)
  vPhi[1] = exp(vPhi_tilde[1])

  # compute the log density at time 1
  dLLK = dt(vY[1]/vPhi[1], dNu, log = TRUE) - log(vPhi[1]) #add the log term since the dt funtcinion works without a scale parameter, so we also have to set x/scale to accoutn for that!!!

  for (t in 2:iT) {

    #update
    vPhi_tilde[t] = dOmega + dAlpha * (((dNu + 1.0) * (vY[t - 1]/vPhi[t - 1])^2)/(dNu + (vY[t - 1]/vPhi[t - 1])^2) - 1.0) +
      dBeta * vPhi_tilde[t - 1]
    #map
    vPhi[t] = exp(vPhi_tilde[t])
    #density at time t
    dLLK = dLLK + dt(vY[t]/vPhi[t], dNu, log = TRUE) - log(vPhi[t])
  }
  #output the results
  lOut = list(vPhi = vPhi,
              dLLK = dLLK)

  return(lOut)

}

Estimate_tGAS <- function(vY) {

  ##Starting values
  vPar = c(omega = log(sd(vY)*sqrt(3/5)) * 0.9,
           alpha = 0.05,
           phi = 0.9,
           nu = 5)

  #optimizer
    optimizer = optim(vPar, function(vPar, vY) {

      dnLLK = -tGAS_LLK(vY, vPar[1], vPar[2], vPar[3], vPar[4])$dLLK

      if (!is.finite(dnLLK)) {
        dnLLK = 1000
      }

      return(dnLLK)

    }, lower = c(-3.0, 0.0001, 0.00001, 2.01), upper = c(3.0, 2.0, 0.9999, 50),
    method = "L-BFGS-B", vY = vY)

    # estimated parameters
    vPar = optimizer$par

    #compute phi (with optimal parameters)
    vPhi = tGAS_LLK(vY, vPar[1], vPar[2], vPar[3], vPar[4])$vPhi

    #compute sigma
    vSigma = vPhi * sqrt(vPar[4]/(vPar[4] - 2.0))

    #output the results
    lOut = list(vSigma = vSigma,
                vPar = vPar,
                vPhi = vPhi)

    return(lOut)

}

# Ex3

library(quantmod)

mPrice = getSymbols("^GSPC", from = "2005-01-01", to = "2018-01-01", auto.assign = FALSE)

vR = as.numeric((diff(log(mPrice[, 6]))*100)[-1])

vR[vR == 0] = mean(vR)

#Fit the tGAS model
Fit_tGAS = Estimate_tGAS(vR)

#Estimated parameters
Fit_tGAS$vPar

# Compare with SV
plot.ts(Fit_tGAS$vSigma)

## you can compare this with the volatility obtained from the previous exercise set.



