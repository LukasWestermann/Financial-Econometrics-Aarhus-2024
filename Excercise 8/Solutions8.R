####################################################################
############################ Exercise 2 ############################
####################################################################

MinimumVariancePortfolio <- function(mSigma) {
  #compute the inverse of the covariance matrix
  mSigmaInv = solve(mSigma)
  #vector of ones
  vOnes = matrix(1, nrow = ncol(mSigma))
  #compute weights
  vOmega = mSigmaInv %*% vOnes
  #normalize weights
  vOmega = vOmega/sum(vOmega)
  return(vOmega)
}


# This function computes the vector of weights
# that minimize the porfolio variance subject to
# a fixed return. (It is slide 4 of lecture 15).
EfficientSet <- function(mSigma, vMu, dK) {

  mSigma_i = solve(mSigma)

  vI = matrix(1, nrow = 2)

  dA = as.numeric(t(vI) %*% mSigma_i %*% vMu)
  dB = as.numeric(t(vMu) %*% mSigma_i %*% vMu)
  dC = as.numeric(t(vI) %*% mSigma_i %*% vI)

  dD = dB*dC - dA^2

  vG = (dB * (mSigma_i %*% vI) - dA*(mSigma_i  %*% vMu))/dD
  vH = (dC*(mSigma_i %*% vMu) - dA * (mSigma_i  %*% vI))/dD

  vOmega = vG + vH*dK
  dSigma2_p = as.numeric(t(vOmega) %*% mSigma %*% vOmega)

  return(list("weights" = vOmega,
              "variance" = dSigma2_p,
              "mean" = dK))
}

###

# Example: Draws the efficient frontier
# Assume these mean and variance
mSigma = matrix(c(5, -4,
                  -4, 15),2, byrow = TRUE)
vMu = c(2, 5)

# different levels of expected returns
vK = seq(0.01, 10.0, 0.01)

mFrontier = t(sapply(vK, function(dK) {
  set = EfficientSet(mSigma, vMu, dK)
  c(set$mean, set$variance)
}))

# plot of the Frontier
plot(mFrontier[,2], mFrontier[,1], type = "l", ylab = "variance", xlab = "mean")


###


EstimateCCC <- function(vY1, vY2) {
  require(rugarch)

  #Model Specification
  ModelSpec = ugarchspec(mean.model = list(armaOrder = c(1, 0)))

  #Model estimation
  Fit_1 = ugarchfit(ModelSpec, vY1)
  Fit_2 = ugarchfit(ModelSpec, vY2)

  #standardized residuas
  vZ_1 = residuals(Fit_1, standardize = TRUE)
  vZ_2 = residuals(Fit_2, standardize = TRUE)

  #unconditional correlation
  mR = cor(cbind(vZ_1, vZ_2))

  ##prediction
  Forc_1 = ugarchforecast(Fit_1, n.ahead = 1)
  Forc_2 = ugarchforecast(Fit_2, n.ahead = 1)

  #one step ahead standard deviations
  vSigma1_tp1 = as.numeric(sigma(Forc_1))
  vSigma2_tp1 = as.numeric(sigma(Forc_2))
  #one step ahead mean
  vMu1_tp1 = as.numeric(fitted(Forc_1))
  vMu2_tp1 = as.numeric(fitted(Forc_2))

  #one step ahead covariance matrix
  mSigma_tp1 = diag(c(vSigma1_tp1, vSigma2_tp1)) %*% mR %*% diag(c(vSigma1_tp1, vSigma2_tp1))
  #one step ahead mean vector
  vMu_tp1 = c(vMu1_tp1, vMu2_tp1)

  #output
  lOut = list()

  lOut[["Fit_1"]] = Fit_1
  lOut[["Fit_2"]] = Fit_2
  lOut[["mR"]] = mR
  lOut[["mSigma_tp1"]] = mSigma_tp1
  lOut[["vMu_tp1"]] = vMu_tp1

  return(lOut)

}

# Estimate the DCC model of Engle (2002) and the Patton's (2006) model

# Filter for the dynamic correlation of the DCC model of Engle (2002)
# mEta is a iN x iT matrix of standardized residuals
# dA and dB are the two parameters of the DCC model
# mQ is the unconditional correlation matrix compute as cor(mEta)
# The function also returns the one step ahead prediction of the correlation matrix
DCCFilter <- function(mEta, dA, dB, mQ) {

  iN = ncol(mEta)
  iT = nrow(mEta)

  # initialize the array for the correlations
  aCor = array(0, dim = c(iN, iN, iT + 1))
  # initialize the array for the Q matrices
  aQ = array(0, dim = c(iN, iN, iT + 1))

  ## initialization at the unconditional cor
  aCor[,, 1] = mQ
  aQ[,,1] = mQ

  #Compute the first likelihood contribution
  dLLK = mEta[1, , drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) -
    mEta[1, , drop = FALSE]%*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1]))

  #main loop
  for (t in 2:(iT + 1)) {
    #update the Q matrix
    aQ[,, t] = mQ * (1 - dA - dB) + dA * t(mEta[t - 1, , drop = FALSE]) %*% mEta[t - 1, , drop = FALSE] +
      dB * aQ[,,t - 1]

    ## Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2}
    aCor[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t])))

    if (t <= iT) {
      #augment the likelihood value
      dLLK = dLLK + mEta[t, , drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t, , drop = FALSE]) -
        mEta[t, , drop = FALSE] %*% t(mEta[t, , drop = FALSE]) + log(det(aCor[,, t]))
    }
  }

  lOut = list()
  #remember to include the -1/2 term in the likelihood evaluation
  #see the equations in the corresponding lecture
  lOut[["dLLK"]] = -0.5 * dLLK
  lOut[["aCor"]] = aCor

  return(lOut)
}

#DCC estimation
# StartingDCCPar is a vector of starting parameters for the DCC estimation
# we will use previous estimates of the DCC parameters as starting value during
# the for loop in the empirical part. This will speed up the estimation.
EstimateDCC <- function(vY1, vY2) {
  require(rugarch)

  #Model Specification
  ModelSpec = ugarchspec(mean.model = list(armaOrder = c(1, 0)))

  #Model estimation -- univariate
  Fit_1 = ugarchfit(ModelSpec, vY1)
  Fit_2 = ugarchfit(ModelSpec, vY2)

  #standardized residuas
  vZ_1 = residuals(Fit_1, standardize = TRUE)
  vZ_2 = residuals(Fit_2, standardize = TRUE)

  #unconditional correlation
  mR = cor(cbind(vZ_1, vZ_2))

  #Model estimation -- multivariate

  ## maximization of the DCC likelihood
  vPar = c(0.04, 0.9)

  #maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ) {

    Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)

  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-4, ineqUB = 0.999,
  LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
  mEta = cbind(vZ_1, vZ_2), mQ = mR)

  vPar = optimizer$pars

  ##prediction
  Forc_1 = ugarchforecast(Fit_1, n.ahead = 1)
  Forc_2 = ugarchforecast(Fit_2, n.ahead = 1)

  #one step ahead standard deviations
  vSigma1_tp1 = as.numeric(sigma(Forc_1))
  vSigma2_tp1 = as.numeric(sigma(Forc_2))
  #one step ahead mean
  vMu1_tp1 = as.numeric(fitted(Forc_1))
  vMu2_tp1 = as.numeric(fitted(Forc_2))

  #Filter DCC
  Filter_DCC = DCCFilter(mEta = cbind(vZ_1, vZ_2), vPar[1], vPar[2], mR)

  #prediction DCC
  mR_tp1 = Filter_DCC$aCor[,, length(vY1) + 1]

  #one step ahead covariance matrix
  mSigma_tp1 = diag(c(vSigma1_tp1, vSigma2_tp1)) %*% mR_tp1 %*% diag(c(vSigma1_tp1, vSigma2_tp1))
  #one step ahead mean vector
  vMu_tp1 = c(vMu1_tp1, vMu2_tp1)

  #output
  lOut = list()

  lOut[["Fit_1"]] = Fit_1
  lOut[["Fit_2"]] = Fit_2
  lOut[["Filter_DCC"]] = Filter_DCC
  lOut[["optimizer"]] = optimizer
  lOut[["vDCCPar"]] = vPar
  lOut[["mR"]] = mR
  lOut[["mSigma_tp1"]] = mSigma_tp1
  lOut[["vMu_tp1"]] = vMu_tp1

  return(lOut)

}

####################################################################
############################ Exercise 3 ############################
####################################################################
library(rugarch)
library(Rsolnp)
data("dji30ret")

#full sample size
iT = 2000

vY1 = tail(dji30ret[, "HPQ"] * 100, iT)
vY2 = tail(dji30ret[, "PG"] * 100, iT)

## we (arbitrarily) set K such that we have an annualized expected return of 7%
dK = 7/225

#length of the out of sample period
# if IF = 1000 you get the answer of the Exercise set.
# here we use IF = 100 because it is faster to run.
iF = 100

aWeights = array(NA, dim = c(iF, 2, 2, 2),
                 dimnames = list(NULL, c("CCC", "DCC"), c("MVP", "FixMean"), c("omega1", "omega2")))

## it is computationally expensive !!! try with a fixed t before.

for (t in (iT - iF):(iT - 1)) {
  #Estimate models and make prediction for time t + 1
  Fit_CCC = EstimateCCC(vY1[(t - iT + iF + 1):t], vY2[(t - iT + iF + 1):t])
  Fit_DCC = EstimateDCC(vY1[(t - iT + iF + 1):t], vY2[(t - iT + iF + 1):t])

  #Compute Tangency Portfolio
  aWeights[t - iT + iF + 1, "CCC", "FixMean", ] = EfficientSet(mSigma = Fit_CCC$mSigma_tp1, vMu = Fit_CCC$vMu_tp1, dK)$weights
  aWeights[t - iT + iF + 1, "DCC", "FixMean", ] = EfficientSet(mSigma = Fit_DCC$mSigma_tp1, vMu = Fit_DCC$vMu_tp1, dK)$weights

  #Compute MVP
  aWeights[t - iT + iF + 1, "CCC", "MVP", ] = MinimumVariancePortfolio(mSigma = Fit_CCC$mSigma_tp1)
  aWeights[t - iT + iF + 1, "DCC", "MVP", ] = MinimumVariancePortfolio(mSigma = Fit_DCC$mSigma_tp1)

  cat(paste(t, "\n"))

}

#array of portfolio returns computed according to different strategies
#and models.
mPortfolioReturns = matrix(NA, iF, 4, dimnames = list(NULL, c("CCC-MVP", "CCC-FixMean", "DCC-MVP", "DCC-FixMean")))

mPortfolioReturns[, "CCC-MVP"] = aWeights[, "CCC", "MVP", "omega1"] * tail(vY1, iF) +
  aWeights[, "CCC", "MVP", "omega2"] * tail(vY2, iF)

mPortfolioReturns[, "DCC-MVP"] = aWeights[, "DCC", "MVP", "omega1"] * tail(vY1, iF) +
  aWeights[, "DCC", "MVP", "omega2"] * tail(vY2, iF)

mPortfolioReturns[, "CCC-FixMean"] = aWeights[, "CCC", "FixMean", "omega1"] * tail(vY1, iF) +
  aWeights[, "CCC", "FixMean", "omega2"] * tail(vY2, iF)

mPortfolioReturns[, "DCC-FixMean"] = aWeights[, "DCC", "FixMean", "omega1"] * tail(vY1, iF) +
  aWeights[, "DCC", "FixMean", "omega2"] * tail(vY2, iF)


# Porfolio statistics
mStat = matrix(NA, 5, 4, dimnames = list(c("Mean", "SD", "SR", "Kurtosis", "Skewness"),
                                         c("CCC-MVP", "CCC-TP", "DCC-MVP", "DCC-TP")))

mStat["Mean", ] = colMeans(mPortfolioReturns)
mStat["SD", ] = apply(mPortfolioReturns, 2, sd)
mStat["SR", ] = mStat["Mean", ]/mStat["SD", ]
mStat["Kurtosis", ] = apply(mPortfolioReturns, 2, function(x) mean((x - mean(x))^4)/(sd(x)^4))
mStat["Skewness", ] = apply(mPortfolioReturns, 2, function(x) mean((x - mean(x))^3)/(sd(x)^3))


#Graphical representation of the portfolio weights
par(mfrow = c(2, 1))

par(mar = c(2,3,3,2))

plot.ts(aWeights[,, "MVP", "omega1"], plot.type = "single", col = 1:2,
        ylab = "", xlab = "", main = "Omega1 for CCC-MVP and DCC-MVP")

par(mar = c(2,3,3,2))

legend("topright", legend = c("CCC-MVP", "DCC-MVP"), col = 1:2, lty = 1)

plot.ts(aWeights[,, "TP", "omega1"], plot.type = "single", col = 1:2,
        ylab = "", xlab = "", main = "Omega1 for CCC-FixMean and DCC-FixMean")

legend("topright", legend = c("CCC-FixMean", "DCC-FixMean"), col = 1:2, lty = 1)
