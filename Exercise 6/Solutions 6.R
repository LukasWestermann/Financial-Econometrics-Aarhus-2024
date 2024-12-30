# Estimate the DCC model of Engle (2002) and the Patton's (2006) model

# Filter for the dynamic correlation of the DCC model of Engle (2002)
# mEta is a iN x iT matrix of standardized residuals
# dA and dB are the two parameters of the DCC model
# mQ is the unconditional correlation matrix compute as cor(mEta)
DCCFilter <- function(mEta, dA, dB, mQ) {

  iN = ncol(mEta)
  iT = nrow(mEta)

  # initialize the array for the correlations
  aCor = array(0, dim = c(iN, iN, iT))
  # initialize the array for the Q matrices
  aQ = array(0, dim = c(iN, iN, iT))

  ## initialization at the unconditional cor
  aCor[,, 1] = mQ
  aQ[,,1] = mQ

  #Compute the first likelihood contribution of the correlation component)
  dLLK = mEta[1, , drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) -
    mEta[1, , drop = FALSE]%*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1]))

  #main loop
  for (t in 2:iT) {
    #update the Q matrix
    aQ[,, t] = mQ * (1 - dA - dB) + dA * t(mEta[t - 1, , drop = FALSE]) %*% mEta[t - 1, , drop = FALSE] +
      dB * aQ[,,t - 1]

    ## Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2}
    aCor[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t])))

    #augment the likelihood value
    dLLK = dLLK + mEta[t, , drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t, , drop = FALSE]) -
      mEta[t, , drop = FALSE] %*% t(mEta[t, , drop = FALSE]) + log(det(aCor[,, t]))
  }

  lOut = list()
  #remember to include the -1/2 term in the likelihood evaluation
  #see the equations in the corresponding lecture
  lOut[["dLLK"]] = -0.5 * dLLK
  lOut[["aCor"]] = aCor

  return(lOut)
}

# Function to estimate the DCC model
# y_t = Sigma_t^{1/2}z_t
# Sigma_t = D_t^{1/2} R_t D_t^{1/2}
# where D_t is a diagonal matrix with
# typical element D_iit = sigma_it^2.
# we define with eta_t = D_t^{-1/2} y_t
Estimate_DCC <- function(mY) {

  ## estimate the marginal GARCH models
  require(rugarch) #load the rugarch package
  require(Rsolnp)

  #####################################################
  # The following part of GARCH estimation can be written
  # in a nicer way by using the multifit function in the rugarch
  # package. See help(multifit)
  #############################################################

  #Marginal garch specifications
  SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0, 0)))

  #list where marginal models are stored
  lFit_univariate = list()

  #estimate the univariate GARCH models
  for(n in 1:ncol(mY)) {
    lFit_univariate[[n]] = ugarchfit(SpecGARCH, mY[, n])
  }

  #Compute the residuals
  mEta = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(residuals(Fit, standardize = TRUE))
  }))

  #####################################################

  ## maximization of the DCC likelihood

  #initial parameters
  vPar = c(0.04, 0.9)

  #unconditional correlation
  mQ = cor(mEta)

  #maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ) {

    Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)

  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-4, ineqUB = 0.999,
  LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
  mEta = mEta, mQ = mQ)

  #Extract the estimated parameters
  vPar = optimizer$pars

  #Filter the dynamic correlation using the estimated parameters
  Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)

  #extract univariate volatilities
  mSigma = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(sigma(Fit))
  }))

  #extract univariate estimated parameters
  mCoef = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(coef(Fit))
  }))

  #compute the likelihood of the volatility  part
  dLLK_V = do.call(sum, lapply(lFit_univariate, function(Fit) {
    as.numeric(likelihood(Fit))
  }))

  #compute the total likelihood
  ## this can be computed in the following way or
  ## or as in slides 35 of lecture 10
  aCor = Filter[["aCor"]]
  aCov = array(NA, dim = dim(aCor))
  iT = dim(aCor)[3]
  iN = dim(aCor)[2]
  dkern = 0
  for (t in 1:iT) {
    aCov[,,t] = diag(mSigma[t, ]) %*% aCor[,,t] %*% diag(mSigma[t, ])
    dkern = dkern + as.numeric(determinant(aCov[,,t], logarithm = TRUE)$modulus + t(t(mY[t, ])) %*% solve(aCov[,,t]) %*% t(mY[t, ]))
  }
  dLLK = -0.5*(iT * iN * log(2*pi) + dkern)

  ## Compute z_t

  iT = nrow(mY)

  mZ = matrix(0, iT, ncol(mY))

  for (t in 1:iT) {
    mZ[t, ] = diag(1/mSigma[t, ]) %*% solve(chol(aCor[,,t])) %*% as.numeric(mY[t, ])
  }

  BIC = log(iT) * 8 - 2 * dLLK

  lOut = list()

  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["aCor"]] = aCor
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["BIC"]] = BIC

  return(lOut)

}

# Two step estimation of DCC model with Student's t
# distributed shocks

Estimate_DCC_t <- function(mY) {

  # QML estimation of Sigma_t

  Fit_QML = Estimate_DCC(mY)

  # extract standardized residuals
  mZ = Fit_QML$mZ

  iT = nrow(mZ)
  iP = ncol(mZ)

  # Estimate a Multivariate Student's t
  # on the standardized residuals

  optimizer = optim(5, function(dNu, mZ, iT, iP) {

    dKern = 0.0
    for (t in 1:iT) {
      dKern = dKern + log(1.0 + c(t(mZ[t, ]) %*% (mZ[t, ]))/(dNu - 2.0))
    }
    dLLK = iT * (lgamma(0.5*(dNu + iP)) - lgamma(0.5*dNu) - 0.5 * iP * log(dNu - 2.0)) - 0.5 * (dNu + iP) * dKern

    return(-dLLK)

  }, method = "L-BFGS-B", lower = 2.01, upper = 50, mZ = mZ, iT = iT, iP = iP)

  dNu = optimizer$par

  lOut = list()
  lOut[["Fit_QML"]] = Fit_QML
  lOut[["dNu"]] = dNu

  return(lOut)

}

##############################################################
#                   Estimation part                          #
##############################################################

library(rugarch)
data("dji30ret")
mY = dji30ret[1:1000, 1:2]

#estimate Gaussian DCC
Fit_DCC = Estimate_DCC(mY)

#estimate Student's t DCC in two step with QML for the first step
Fit_DCC_t = Estimate_DCC_t(mY)

# The filtered correlation of the Gaussian and Student's t DCC models are the same
plot.ts(Fit_DCC$aCor[1,2,])

# The two are very similar

# Compare the two DCC and copula models using BIC. Which model is selected?
# This point is subtle because we haven't computed the Likelihood for the DCC model with
# Student's t shocks. Remeber that we have estimated the model by QML !
# We first compute the likelihood of the DCC t model

mY = as.matrix(mY)
iT = nrow(mY)
iP = ncol(mY)

mSigma_DCC_t = Fit_DCC_t$Fit_QML$mSigma
aCor_DCC_t   = Fit_DCC_t$Fit_QML$aCor
dNu = Fit_DCC_t$dNu

dLLK_DCC_t = 0

for (t in 1:iT) {

  mSigma = diag(mSigma_DCC_t[t, ]) %*% aCor_DCC_t[,,t] %*% diag(mSigma_DCC_t[t, ])

  dLLK_DCC_t = dLLK_DCC_t +
    lgamma(0.5*(dNu + iP)) - lgamma(0.5*dNu) - 0.5 * iP * log(dNu - 2.0) - 0.5 * log(det(mSigma)) +
    -0.5 * (dNu + iP) * log(1.0 + c(mY[t, ] %*% solve(mSigma) %*% mY[t, ]) /(dNu - 2.0))
}

vBIC = c(
  BIC_DCC_t        = log(iT) * 9 - 2 * dLLK_DCC_t,
  BIC_DCC_Gauss    = Fit_DCC$BIC
)

sort(vBIC)

# The model selected by BIC is DCC_t

