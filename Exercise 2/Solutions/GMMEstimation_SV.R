
## This function computes E[sigma_t^p]
E.sigma_p <- function(p, dAlpha,
                      dBeta2) {
  dE = exp(p * dAlpha/2 + p^2*dBeta2/8)
  return(dE)
}

## This function computes E[sigma_t^p sigma_{t-j}^s]
E.sigma_p_sigmaj_s <- function(p, s, j, dAlpha, dPhi,
                               dBeta2) {

  Esigma_p = E.sigma_p(p, dAlpha, dBeta2)
  Esigma_s = E.sigma_p(s, dAlpha, dBeta2)

  dEsigma_p_sigmaj_s = Esigma_p * Esigma_s *
    exp(p * s * dPhi^j * dBeta2/4)

  return(dEsigma_p_sigmaj_s)

}


# This function computes the theoretical moments of the SV model
# for a given set of parameters in vPar. It returns a vector of length
# 24 with the theoretical moments ordered as in slide 21 of lecture 7
f.theo <- function(vPar) {

  dOmega = vPar[1]
  dPhi = vPar[2]
  dSigma2 = vPar[3]

  dAlpha = dOmega/(1 - dPhi)
  dBeta2 = dSigma2/(1 - dPhi^2)

  vM = numeric(24)

  # E[|r_t|]
  vM[1] = sqrt(2/pi) * E.sigma_p(1, dAlpha, dBeta2)
  # E[r_t^2]
  vM[2] = E.sigma_p(2, dAlpha, dBeta2)
  # E[|r_t|^3]
  vM[3] = 2 * sqrt(2/pi) * E.sigma_p(3, dAlpha, dBeta2)
  # E[r_t^4]
  vM[4] = 3 * E.sigma_p(4, dAlpha, dBeta2)
  # E[|r_t r_{t-j}|]
  for (j in 1:10) {
    vM[4 + j] = (2/pi) * E.sigma_p_sigmaj_s(1, 1, j, dAlpha, dPhi,
                                            dBeta2)
  }
  # E[r_t^2 r_{t-j}^2]
  for (j in 1:10) {
    vM[14 + j] = E.sigma_p_sigmaj_s(2, 2, j, dAlpha, dPhi,
                                    dBeta2)
  }

  return(vM)
}


# This function computes the empirical moments of a vector of observations
# It returns a matrix of length T - 10 x 24. Note that 10 observations are lost
# when computing y_t^2 * y_{t-10}^2.
f.EmpiricalMoments <- function(vY) {

  iT = length(vY)

  mM = matrix(NA, iT, 24)

  mM[, 1] = abs(vY)
  mM[, 2] = vY^2
  mM[, 3] = abs(vY)^3
  mM[, 4] = vY^4

  for (j in 1:10) {
    mM[-(1:j), 4 + j] = abs(vY[-(1:j)]) * abs(vY[-((iT - j + 1):iT)])
  }

  for (j in 1:10) {
    mM[-(1:j), 14 + j] = vY[-(1:j)]^2 * vY[-((iT - j + 1):iT)]^2
  }

  mM = mM[-(1:10), ]

  return(mM)

}

# This function computes the S_HAC matrix for efficient GMM estimation
# of the SV model as in slide 17 of lecture 7. It accepts a iT - 10 x 24
# matrix, mG, for which the generic row is given by the difference between
# the theoretical and empirical moments (the g_t(w_t, theta) vector).
# The function uses a triangular kernel k(x) = 1 - |x| for |x| < 1. Other choices
# can be implemented, see for example help(kernel). As in slide 17 we average
# all lag-j cross moments of the g_t(w_t, theta) variable using the weights obtained by the kernel.
# The lag-j matrix of cross moments is obtained via the acf() function (this is the last equation
# of slide 17). In acf() we set type = "covariance" and demean = FALSE to indicate that
# we want to compute the empirical counterpart of E[g_t(w_t, theta) g_{t-j}(w_{t-j}, theta)] for
# j = 1, 2, ..., lag.max. See help(acf). The threshold q(T) is set to round(iT^(1/3)), that is,
# we compute autocovariances up to round(iT^(1/3)) of the sample.
f.Make_S_HAC <- function(mG) {

  iT = nrow(mG)

  # this is the bandwith used by
  # Newey, W. K., & West, K. D. (1987). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation. Econometrica, 55(3), 703-708.
  Q_t = round(iT^(1/3))

  ## This function gives you the autocovariance function. The argument demean = FALSE imply that E[X_t X_{t-j}'] is computed
  ## and not E[(X_t - E[X_t])(X_{t-j} - E[X_{t-j}])'] as in the usual covariance formula
  mACF = acf(mG, lag.max = Q_t + 1, plot = FALSE, demean = FALSE, type = "covariance")$acf

  mS = matrix(0, ncol(mG), ncol(mG))

  # we use the triangular kernel k(x) = 1 - |x| for |x| < 1

  for (j in 1:Q_t) {
    mS = mS + (1 - j/(Q_t + 1)) * (mACF[j + 1,,] + t(mACF[j + 1,,]))
  }

  mS =  mS + mACF[1,,]

  return(mS)

}

# This function estimate the SV model using the efficient GMM estimator.
GMM_Estimator <- function(vY) {

  iT = length(vY)
  
  ## vector of initial parameters
  ## we use MM in slide 10 of Lecture 10 
  ## to initialize omega and sigma
  ## we set phi = 0.95
  kappa_hat = mean(abs(vY))
  mu_2_hat = mean(vY^2)
  alpha_hat = log(pi^2*kappa_hat^4/(4*mu_2_hat))
  beta2_hat = log(16*mu_2_hat^4/(pi^4*kappa_hat^8))
  
  vPar = c(alpha_hat*(1-0.95), 0.95, beta2_hat*(1-0.95^2))

  ## moments (note we need to compute this just once)
  mM = f.EmpiricalMoments(vY)

  # here we minimize the quadratic form T g_T(theta)' W g_T(theta),
  # where W is the inverse of S_HAC. See slides 14 and 16 of lecture 7.
  # Note that we compute the empirical moments outside the optimizer which
  # leads to computational advantages.
  optimizer = optim(vPar, function(vPar, mM, iT) {

    dLoss = 1e7
    
    try({
    
    vTheo = f.theo(vPar)

    # this is the matrix for which the generic row is given by the difference between
    # the theoretical and empirical moments (the g_t(w_t, theta) vector).
    mG = t(t(mM) - vTheo)

    # # Efficient weight matrix
    mS_HAC = f.Make_S_HAC(mG)
    # ## ginv() computes the generalized inverse (slower but more robust.
    # ## the usual inverse is computed by solve()
    mW = MASS::ginv(mS_HAC)
    # mW = diag(length(vTheo))
    
    # the g_T vector at slide 14 of lecture 7
    vAvgMoments = colMeans(mG)

    # the quadratic form J(theta, W) at slide 14 of lecture 7
    dLoss = (vAvgMoments %*% mW %*% vAvgMoments) * iT
    }, silent = TRUE)

    return(dLoss)

  }, method = "L-BFGS-B", lower = c(-2, 0.0001, 0.001), upper = c(1.0, 0.999, 1.0),
  mM = mM, iT = iT)

  return(optimizer)

}
