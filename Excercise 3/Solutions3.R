
#####################################
#     SOLUTIONS TO EXERCISE (2)     #
#####################################

##### POINT 1)

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


##### POINT 2)
# Function to simulate from the model:
# y_t = sigma*exp(w_t/2)*z_t, z_t ~ iid N(0,1)
# w_t = rho * w_{t-1} + eta_t, eta_t ~ iid N(0, sigma_eta^2)
#
# iT is the length of sample to simulate
# dSigma is the constant in the measurement equation
# dRho is the autoregressive parameter of the log volatility
# dSigma2_eta is the variance of the volatility shocks
SV_sim <- function(iT, dSigma, dRho, dSigma2_eta) {

  vW = numeric(iT)
  vEta = rnorm(iT, mean = 0, sd = sqrt(dSigma2_eta))
  vZeta = rnorm(iT)
  ## initialize w sampling from the unconditional distribution
  # also setting vW[1] = 0 (i.e. vW[1] = E[w_t]) is fine
  vW[1] = rnorm(1, 0, sqrt(dSigma2_eta/(1-dRho^2)))

  ## recursion for the volatility
  for (t in 2:iT) {
    vW[t] = dRho*vW[t-1] + vEta[t]
  }

  # compute the returns
  vY = dSigma*exp(vW/2)*vZeta

  #output the returns and the volatility
  lOut = list(vY = vY,
              vW = vW)

  return(lOut)
}


### POINT 3
set.seed(123)

iT = 5000
dSigma = 1.0
dRho = 0.9 #(in the question it is phi)
dSigma2_eta = 0.25

lSim = SV_sim(iT, dSigma, dRho, dSigma2_eta)

### POINT 4

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

### POINT 5

vY = lSim$vY
Fit = QML_SV(vY)

# estimated parameters are
# sigma = 1.0671994
# rho   = 0.9073020
# sigma2eta = 0.1423963
# We see that estimated for sigma and rho
# are reasonable but sigma2eta is heavily
# underestimated. If you increase T = 5000
# using the same seed you obtain
# sigma = 0.9915461
# rho   = 0.8900296
# sigma2eta = 0.2628427
# such that, as expected, the precision of the estimates
# increases with T
Fit$vPar

plot.ts(lSim$vW)
lines(c(Fit$Filter$a_filt), col = "red")

plot.ts(exp(lSim$vW)[1:1000])
lines(exp(c(Fit$Filter$a_filt)), col = "red")

#####################################
#     SOLUTIONS TO EXERCISE (3)     #
#####################################

library(quantmod)

## download the series
mPrice = getSymbols("^GSPC", from = "2005-01-01", to = "2018-01-01", auto.assign = FALSE)

## compute log returns
vY = as.numeric((diff(log(mPrice[, 6]))*100)[-1])

## replace zeros with the empirical mean
vY[vY == 0] = mean(vY)

## estimate the model by QML (may take a bit of time)
Fit_QML = QML_SV(vY)

### From Exercise set 2
# setwd("G:/Dropbox/Teaching/Aarhus_Financial Econometrics/2021/ExerciseSets/2/Solutions/")
setwd("C:/Users/au588008/Dropbox/Teaching/Aarhus_Financial Econometrics/2021/ExerciseSets/2/Solutions/")
source("GMMEstimation_SV.R")

## estimate the model by GMM (may take a bit of time)
Fit_GMM = GMM_Estimator(vY)

## see estimated parameters
Fit_GMM$par

## see estimated parameters
Fit_QML$vPar

# map between the GMM and QML estimate
# sigma = exp(omega/(2*(1-rho)))
exp(Fit_GMM$par[1]/(2*(1.0 - Fit_GMM$par[2]))) # 0.6991863

# should be similar to
Fit_QML$vPar[1] # 0.7453121

# We see that estimates of the
# the autoregressive parameter and
# the variance of the volatility shock are quite different
# between the GMM and QML estimator. Estimated from the QML
# estimator seems much more plausible and should be preferred.
# estimates of the intercept are similar between the two
# estimator.

