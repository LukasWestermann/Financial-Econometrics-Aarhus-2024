###########################
# (1): Computational part #
###########################

LSE <- function(vX) {
  
  dM = max(vX)
  dOut = dM + log(sum(exp(vX-dM)))
  return(dOut)
  
}

### Exercise 1)

BootstrapFilter <- function(vY, dOmega, dPhi, dTau, iN = 10000) {

  iT = length(vY)

  ## matrix of particles
  mAlpha = matrix(NA, iT, iN)
  ## vector with filtered volatility values
  vVol = numeric(iT)
  ## importance weight at iteration t
  vW = numeric(iN)

  ##initialize at time t = 1
  #sample from the unconditional dist as an initial proposal distribution
  mAlpha[1, ] = rnorm(iN, dOmega/(1.0 - dPhi), sqrt(dTau^2/(1 - dPhi^2)))
  #compute the importance weights
  vW = dnorm(vY[1], 0, exp(mAlpha[1, ]/2.0))
  #normalized importance weights
  vW = vW/sum(vW)
  #approximate E[exp(alpha_1) | y_{1}]
  vVol[1] = sum(exp(mAlpha[1, ]/2) * vW)
  #resample
  mAlpha[1, ] = sample(mAlpha[1, ], iN, replace = TRUE, prob = vW)

  for (t in 2:iT) {
    #sample from the proposal distribution p(alpha_t|alpha_{t-1})
    mAlpha[t, ] = dOmega + dPhi * mAlpha[t - 1, ] + rnorm(iN, 0, dTau)
    #compute the importance weights
    vW = dnorm(vY[t], 0, exp(mAlpha[t, ]/2.0))
    #normalized importance weights
    vW = vW/sum(vW)
    #approximate E[exp(alpha_t) | y_{1:t}]
    vVol[t] = sum(exp(mAlpha[t, ]/2) * vW)
    #resample
    mAlpha[t, ] = sample(mAlpha[t, ], iN, replace = TRUE, prob = vW)

    ### NOTE THAT SINCE WE RESAMPLE AT EACH ITERATION WE AVOID THE STEP vW = 1
    ### CONSEQUENTLY WHEN WE COMPUTE THE IMPORTANCE WEIGHTS AT TIME t WE NEGLECT
    ### THE EFFECT OF THE IMPORTANCE WEIGHTS AT TIME t-1, SEE LECTURE 10
  }

  return(vVol)

}


BootstrapFilter_LSE <- function(vY, dOmega, dPhi, dTau, iN = 10000) {
  
  iT = length(vY)
  
  ## matrix of particles
  mAlpha = matrix(NA, iT, iN)
  ## vector with filtered volatility values
  vVol = numeric(iT)
  ## importance weight at iteration t
  vLW = numeric(iN)
  
  ##initialize at time t = 1
  #sample from the unconditional dist as an initial proposal distribution
  mAlpha[1, ] = rnorm(iN, dOmega/(1.0 - dPhi), sqrt(dTau^2/(1 - dPhi^2)))
  #compute the importance weights
  vLW = dnorm(vY[1], 0, exp(mAlpha[1, ]/2.0), log = TRUE)
  #normalized importance weights
  vW = exp(vLW - LSE(vLW))
  #approximate E[exp(alpha_1) | y_{1}]
  vVol[1] = sum(exp(mAlpha[1, ]/2) * vW)
  #resample
  mAlpha[1, ] = sample(mAlpha[1, ], iN, replace = TRUE, prob = vW)
  
  for (t in 2:iT) {
    #sample from the proposal distribution p(alpha_t|alpha_{t-1})
    mAlpha[t, ] = dOmega + dPhi * mAlpha[t - 1, ] + rnorm(iN, 0, dTau)
    #compute the importance weights
    vLW = dnorm(vY[t], 0, exp(mAlpha[t, ]/2.0), log = TRUE)
    #normalized importance weights
    vW = exp(vLW - LSE(vLW))
    #approximate E[exp(alpha_t) | y_{1:t}]
    vVol[t] = sum(exp(mAlpha[t, ]/2) * vW)
    #resample
    mAlpha[t, ] = sample(mAlpha[t, ], iN, replace = TRUE, prob = vW)
    
    ### NOTE THAT SINCE WE RESAMPLE AT EACH ITERATION WE AVOID THE STEP vW = 1
    ### CONSEQUENTLY WHEN WE COMPUTE THE IMPORTANCE WEIGHTS AT TIME t WE NEGLECT
    ### THE EFFECT OF THE IMPORTANCE WEIGHTS AT TIME t-1, SEE LECTURE 10
  }
  
  return(vVol)
  
}


### Exercise 2)
BootstrapFilter_EES <- function(vY, dOmega, dPhi, dTau, iN = 10000, dG = 0.75) {

  iT = length(vY)

  ## matrix of particles
  mAlpha = matrix(NA, iT, iN)
  ## vector with filtered volatility values
  vVol = numeric(iT)
  ## matrix importance weight at iteration t
  ## A more efficient version of the code can be written by
  ## keeping track only of two consecutive vectors of particles
  mW_tilde = matrix(NA, iT, iN)
  ## vector of normalized importance weights
  vW = numeric(iN)

  ##initialize at time t = 1
  #sample from the unconditional dist as an initial proposal distribution
  mAlpha[1, ] = rnorm(iN, dOmega/(1.0 - dPhi), sqrt(dTau^2/(1 - dPhi^2)))
  #compute the importance weights
  mW_tilde[1, ] = dnorm(vY[1], 0, exp(mAlpha[1, ]/2.0))
  #normalized importance weights
  vW = mW_tilde[1, ]/sum(mW_tilde[1, ])
  #approximate E[exp(alpha_1) | y_{1}]
  vVol[1] = sum(exp(mAlpha[1, ]/2) * vW)
  #resample
  dESS = 1.0/sum(vW^2)
  if (dESS < dG * iN) {
    mAlpha[1, ] = sample(mAlpha[1, ], iN, replace = TRUE, prob = vW)
    mW_tilde[1, ] = 1.0
  }
  for (t in 2:iT) {
    #sample from the proposal distribution p(alpha_t|alpha_{t-1})
    mAlpha[t, ] = dOmega + dPhi * mAlpha[t - 1, ] + rnorm(iN, 0, dTau)
    #compute the importance weights
    mW_tilde[t, ] = mW_tilde[t-1, ] * dnorm(vY[t], 0, exp(mAlpha[t, ]/2.0))
    #normalized importance weights
    vW =  mW_tilde[t, ]/sum(mW_tilde[t, ])
    #approximate E[exp(alpha_t) | y_{1:t}]
    vVol[t] = sum(exp(mAlpha[t, ]/2) * vW)
    #resample
    dESS = 1.0/sum(vW^2)
    if (dESS < dG * iN) {
      mAlpha[t, ] = sample(mAlpha[t, ], iN, replace = TRUE, prob = vW)
      mW_tilde[t, ] = 1.0
    }
  }

  return(vVol)

}


### Exercise 3)
# simulate the SVAR model with different parameterization
# y_t = exp(alpha_t/2)*u_t
# alpha_{t+1} = omega + phi * alpha_t + tau * eta_t
SV_sim <- function(iT, dOmega, dPhi, dTau) {

  vEps = rnorm(iT)
  vAlpha = numeric(iT)

  vAlpha[1] = rnorm(1, dOmega/(1.0 - dPhi), sqrt(dTau^2/(1 - dPhi^2)))
  if (iT > 1) {
    for (t in 2:iT) {
      vAlpha[t] = dOmega + dPhi * vAlpha[t - 1] + dTau * rnorm(1)
    }
  }

  vY = exp(vAlpha/2) * vEps

  return(list(vY = vY,
              vAlpha = vAlpha,
              vSigma = exp(vAlpha/2)))


}

# set the seed
set.seed(123)

# model parameters
iT = 500
dOmega = 0.0
dPhi = 0.9
dTau = 0.5

# simulate
lSim = SV_sim(iT, dOmega, dPhi, dTau)

# extract the realizations
vY = lSim$vY

plot.ts(vY)
plot.ts(lSim$vSigma)
### Exercise 4)
# number of particles
iN = 10000

#compute the filtered volatility using N = 10000 particles
set.seed(123)
vVol_Boot = BootstrapFilter(vY, dOmega, dPhi, dTau, iN)
vVol_Boot_LSE = BootstrapFilter_LSE(vY, dOmega, dPhi, dTau, iN)

lines(vVol_Boot, col = "red")
lines(vVol_Boot_LSE, col = "blue")


vY = c(1.3, -1)*40

vW_tilde = dnorm(vY)

vW_tilde/sum(vW_tilde)

vlogW_tilde = dnorm(vY, log = TRUE)

exp(vlogW_tilde - LSE(vlogW_tilde))



r## Note vVol_Boot is the same as vVol_Boot_g1 if we re set the seed
set.seed(123)
vVol_Boot_g1 = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 1.0)

all.equal(vVol_Boot, vVol_Boot_g1) #TRUE

## we re set the seed each time in order to have the same draw. In this way,
## the only differences between the filtered volatilities are due to the
## g parameter. Of course this is not required in a general setting.
set.seed(123)
vVol_Boot_g075 = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 0.75)
set.seed(123)
vVol_Boot_g05 = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 0.5)

plot.ts(lSim$vSigma)
lines(vVol_Boot_g1, col = "red")
lines(vVol_Boot_g075, col = "blue")
lines(vVol_Boot_g05, col = "purple")

### There is no much difference! The reason is that we have employed a large number of
### particles. Let's repeat with only N = 50

# number of particles
iN = 50

#compute the filtered volatility using N = 10000 particles
set.seed(123)
vVol_Boot = BootstrapFilter(vY, dOmega, dPhi, dTau, iN)

## Note vVol_Boot is the same as vVol_Boot_g1 if we re set the seed
set.seed(123)
vVol_Boot_g1 = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 1.0)

all.equal(vVol_Boot, vVol_Boot_g1) #TRUE

## we re set the seed each time in order to have the same draw. This is not
## necessary
set.seed(123)
vVol_Boot_g075 = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 0.75)
set.seed(123)
vVol_Boot_g05 = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 0.5)

plot.ts(lSim$vSigma)
lines(vVol_Boot_g1, col = "red")
lines(vVol_Boot_g075, col = "blue")
lines(vVol_Boot_g05, col = "purple")

## Now differences are more pronounced.

# Which is closer to the true one?
mean((lSim$vSigma-vVol_Boot_g1)^2)
mean((lSim$vSigma-vVol_Boot_g075)^2)
mean((lSim$vSigma-vVol_Boot_g05)^2)

## The one with g = 0.75 reports the lowest mean
## squared error. The reason is that, reseampling at
## each step increases the Monte Carlo variance,
## which might translate in filtered estimates with higher
## variance!


### Exercise 5)
############################################
#   Below is the code from Exercise set 3  #
############################################
# Kalman filter and smoother for the state space:
# Y_t         = Z * alpha_t + eps_t, eps_t ~ N(0, S)
# alpha_{t+1} = T * alpha_t + H*eta_t, eta_t ~N(0, Q)
#
# Y_t is p x 1
# Z   is p x m
# S   is p x p
# T   is m x m
# H   is m x l
# Q   is l x l
# a1 is the initialization for a
# P1 is the initialization for P
#
KF_R <- function(mY, mZ, mS, mT, mH, mQ, a1, P1, Smoothing = TRUE) {
  n = ncol(mY)
  p = nrow(mY)
  r_int = ncol(mH)
  m = length(a1)
  
  v = matrix(0, p,n);
  F = array(0, dim = c(p, p, n));
  K = array(0, dim = c(m, p, n));
  a_filt = matrix(0, m,n);
  a_pred = matrix(0, m,n);
  P_filt = array(0, dim = c(m, m, n));
  P_pred = array(0, dim = c(m,m,n));
  
  r = matrix(0, m,n);
  N = array(0, dim = c(m,m,n));
  a_smoot = matrix(0, m, n);
  V = array(0, dim = c(m, m, n));
  L = array(0, dim = c(m,m,n));
  
  eps_smoot = matrix(0, p,n);
  eta_smoot = matrix(0, r_int, n);
  vLLK = numeric(n);
  
  # //initialise
  v[, 1] = mY[, 1];
  a_filt[, 1] = a1;
  a_pred[, 1] = a1;
  P_filt[,,1] = P1;
  P_pred[,,1] = P1;
  
  HQH = mH %*% mQ %*% t(mH);
  dC = -0.5 * (n * p * 1.0) * log(pi * 2.0);
  
  # //filtering recursion
  for(t in 1:n) {
    #prediction error
    v[, t] = mY[, t] - mZ %*% a_pred[, t];
    # variance of the prediction error
    F[,,t] = mZ %*% P_pred[,,t] %*% t(mZ) + mS;
    # filtered state mean E[alpha_t |Y_1:t]
    a_filt[, t] = a_pred[,t] + P_pred[,,t] %*% t(mZ) %*% solve(F[,,t]) %*% v[,t];
    # filtered state variance Var[alpha_t |Y_{1:t}]
    P_filt[,,t] = P_pred[,,t] - P_pred[,,t] %*% t(mZ) %*% solve(F[,,t]) %*% mZ %*% P_pred[,,t];
    # kalman gain
    K[,,t]      = mT %*% P_pred[,,t] %*% t(mZ) %*% solve(F[,,t]);
    # likelihood contribution
    vLLK[t] = (log(det(as.matrix(F[,,t]))) + c(v[,t] %*% solve(F[,,t]) * v[,t]));
    if(t < n){
      # predicted state mean E[alpha_{t+1}|Y_{1:t}]
      a_pred[,t + 1] = mT %*% a_pred[,t] + K[,,t] %*% v[, t];
      # predicted state variance Var[alpha_{t+1}|Y_{1:t}]
      P_pred[,,t + 1] = mT %*% P_pred[,,t] %*% t((mT - K[,,t] %*% mZ)) + HQH;
    }
  }
  # //Smoothing recursion
  if(Smoothing){
    for(t in n:2) {
      L[,,t]= mT - K[,,t] %*% mZ;
      r[, t - 1] = t(mZ) %*% solve(F[,,t]) %*% v[, t] + t(L[,,t]) %*% r[, t];
      N[,,t - 1] = t(mZ) %*% solve(F[,,t]) %*% mZ + t(L[,,t]) %*% N[,,t] %*% L[,,t];
      # smoothed state mean E[alpha_t | Y_{1:n}]
      a_smoot[, t] = a_pred[, t] + P_pred[,,t] %*% r[, t - 1];
      # smoothed state variance Var[alpha_t | Y_{1:n}]
      V[,,t] = P_pred[,,t] - P_pred[,,t] %*% N[,,t - 1] %*% P_pred[,,t];
      # //error smoothing
      eps_smoot[, t] = mS %*% (solve(F[,,t]) %*% v[, t] - t(K[,,t]) %*% r[, t]);
      eta_smoot[, t] = mQ %*% t(mH) %*% r[, t];
    }
  }
  
  KF = list();
  
  KF[["v"]] = v;
  KF[["a_filt"]] = a_filt;
  KF[["a_pred"]] = a_pred;
  KF[["P_filt"]] = P_filt;
  KF[["P_pred"]] = P_pred;
  KF[["F"]] = F;
  KF[["K"]] = K;
  KF[["N"]] = N;
  KF[["a_smoot"]] = a_smoot;
  KF[["V"]] = V;
  KF[["L"]] = L;
  KF[["eps_smoot"]] = eps_smoot;
  KF[["eta_smoot"]] = eta_smoot;
  KF[["vLLK"]] = vLLK;
  KF[["dLLK"]] = -dC - 0.5 * sum(vLLK);
  
  return(KF);
  
}

# THIS FUNCTION WILL BE USEFUL LATER!
# Kalman filter and smoother for the model
# y_t = alpha_t + sigma * eps_t
# alpha_{t+1} = phi * alpha_t + eta * xi_t
# vPar is the vector of parameters vPar = (phi, sigma, eta)
# note that the vPar contains the standard deviations (sigma and eta)
# and not the variances (sigma^2, eta^2).
# vY is the vector of observations
KalmanFilter_AR1plusNoise <- function(vY, vPar, Smoothing = FALSE) {
  
  dPhi = vPar[1]
  dSigma = vPar[2]
  dEta = vPar[3]
  
  mY = matrix(vY, nrow = 1)
  mZ = matrix(1, 1, 1)
  mS = matrix(dSigma^2, 1, 1)
  mT = matrix(dPhi, 1, 1)
  mH = matrix(1, 1, 1)
  mQ = matrix(dEta^2, 1, 1)
  a1 = 0
  mP1 = matrix(dEta^2/(1-dPhi^2), 1, 1)
  
  KF_R(mY, mZ, mS, mT, mH, mQ, a1, mP1, Smoothing)
  
}

# This function maximizes the Quasi Log Likelihood
# computed by the Kalman Filter for a log-linearized
# stochastic volatility model.
# The SV model is
# y_t = sigma*exp(w_t/2)*z_t, z_t ~ iid N(0,1)
# w_t = rho * w_{t-1} + eta_t, eta_t ~ iid N(0, sigma_eta^2)
# It's log linearized version is
# zeta_t = mu + w_t + eps_t,
# w_t = rho * w_{t-1} + eta_t, eta_t ~ iid N(0, sigma_eta^2)
# We have that E[eps_t|F_{t-1}] = 0, and E[eps_t^2|F_{t-1}] = pi^2/2
# and mu = log(sigma^2) + E[log(z_t^2)]
# where E[log(z_t^2)] = -1.270376
## vY are the returns
QML_SV <- function(vY) {

  ## Be sure that there are no zeros
  if(any(vY==0)) {
    vY[vY==0] = mean(vY)
  }

  # 1) transform the data
  vXi = log(vY^2)

  # Estimate sigma
  dSigma = exp((mean(vXi) + 1.270376)/2)

  # 2) Demean vXi
  vXi = vXi - mean(vXi)

  iT = length(vXi)

  #just two parameters rho and sigma2_tau, we constrain sigma_eps^2 = pi^2/2
  # starting values
  vPar = c(cor(vXi[-1], vXi[-iT]), var(vXi) * 0.1)

  # maximize the likelihood
  optimizer = optim(vPar, function(vPar, vXi, iT) {

    dRho = vPar[1]
    dSigma2_eta = vPar[2]

    vPar_tot = c(dRho, sqrt(pi^2/2), sqrt(dSigma2_eta))

    KF = KalmanFilter_AR1plusNoise(vXi, vPar_tot, Smoothing = FALSE)

    dLLK = KF$dLLK

    return(-dLLK)

  }, method = "L-BFGS-B", lower = c(0.001, 0.001), upper = c(0.99, 1.0),
  vXi = vXi, iT = iT)

  #extract estimated parameters
  vPar = optimizer$par
  dRho = vPar[1]
  dSigma2_eta = vPar[2]

  #compute the vector of total parameters
  vPar_tot = c(dRho, sqrt(pi^2/2), sqrt(dSigma2_eta))

  # filtering and smoothing
  Filter = KalmanFilter_AR1plusNoise(vY, vPar_tot, Smoothing = TRUE)

  vModelPar = c(sigma = dSigma,
                rho = dRho,
                sigma2eta = dSigma2_eta)

  lOut = list(vPar = vModelPar,
              Filter = Filter)

  return(lOut)

}


### Estimate the model by QML
Fit = QML_SV(vY)

## The estimated parameters are for the formulation of
## slide 31 of lecture 7, i.e. (sigma, rho, sigma2_eta)
Fit$vPar
## Let's use the results from the results of exercise set 3
## in order to map these estimates into the parameterization we
## consider in this code.

dOmega_hat = log(Fit$vPar["sigma"]^2)*(1.0 - Fit$vPar["rho"])
dPhi_hat = Fit$vPar["rho"]
dTau_hat = sqrt(Fit$vPar["sigma2eta"])

## The estimated omega is -0.005, the true one is 0
## The estimated phi is  0.925, the true one is 0.9
## The estimated tau is 0.391, the true one is 0.5

## Let's perform filtering usinf the estimated parameters. We use g = 1

### Exercise 4)
# number of particles
iN = 10000

#Filtering with the true parameters and g = 0.75
set.seed(123)
vVol_Boot_TRUE = BootstrapFilter_EES(vY, dOmega, dPhi, dTau, iN, dG = 0.75)
#Filtering with the estimated parameters and g = 0.75
set.seed(123)
vVol_Boot_Est = BootstrapFilter_EES(vY, dOmega_hat, dPhi_hat, dTau_hat, iN, dG = 0.75)

# Let's compare the results in a figure

plot.ts(lSim$vSigma)
lines(vVol_Boot_TRUE, col = "red")
lines(vVol_Boot_Est, col = "blue")

## Results are similar, and estimator error does not play a crucial role in our
## example. We are then satisfied with the results.
## Note that, this example shows us that we will never be able to recover
## the "true" volatility recursion.

###########################
# (1): Real Data          #
###########################

library(quantmod)

### Exercise 1)
## download the series
mPrice = getSymbols("^GSPC", from = "2005-01-01", to = "2018-01-01", auto.assign = FALSE)

## compute log returns
vY = as.numeric((diff(log(mPrice[, 6]))*100)[-1])

## replace zeros with the empirical mean
vY[vY == 0] = mean(vY)

### Exercise 2)
## estimate the model by QML (may take a bit of time)
Fit_QML = QML_SV(vY)
## Map the parameters in our parametrization
dOmega_hat = log(Fit_QML$vPar["sigma"]^2)*(1.0 - Fit_QML$vPar["rho"])
dPhi_hat = Fit_QML$vPar["rho"]
dTau_hat = sqrt(Fit_QML$vPar["sigma2eta"])

### Exercise 3)
vSigma_SV = BootstrapFilter_EES(vY, dOmega_hat, dPhi_hat, dTau_hat, iN = 10000, dG = 1.0)

### Exercise 4)
library(rugarch)

# specify a GARCH(1,1) model with Gaussian shocks. The model is
# y_t = sigma_t eps_t, with eps_t iid N(0, 1)
# sigma_t^2 = omega + alpha y_{t-1}^2 + beta sigma_{t-1}^2
spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                  distribution.model = "norm")

Fit_GARCH = ugarchfit(spec, vY)

## extract the estimated volatility from GARCH
# if you don't remember the sigma() function do
# help(ugarchfit), scroll down until you reach the "Value"
# part of the help, and click on uGARCHfit. Here you have
# all the functions to extract results from an estimated
# garch model.

## here we use as.numeric because sigma() returns an xts object
vSigmaGarch = as.numeric(sigma(Fit_GARCH))

### Exercise 5)
plot(vSigma_SV, type = "l", las = 1, ylab = "volatility", xlab = "time")
lines(vSigmaGarch, col = "red")
legend("topright", legend = c("SV", "GARCH"), col = c("black", "red"),
       lty = 1)

## results are very similar!
