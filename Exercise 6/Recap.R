## Estimate Gaussian DCC model

#Filter for dynmaic correlations

DCCFilter = fucntion(mEta, dA, dB, mQ) {
  
  iN = col(mEta)
  iT = nrow(mEta)
  
  #initalize array of correlations
  aCor = arrary(0, dim = c(iN,iN, iT))
  #iniatlize array for the correlation matrices
  aQ =  arrary(0, dim = c(iN, iN, iT))
  
  #initalize uncondtional correlatinon
  aCor[,, 1] = mQ
  aQ[,, 1] = mQ
  
  #Compute first likelihood contribution
  dLLK = mEta[1,,drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) -
    mEta[1, , drop = FLASE] %*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1]))
  
  #main loop for filtered correlations
  for(t in 2:iT){
    aQ[,, t ] = mQ * (1- dA - dB) + dA *t(mEta[t-1,,drop = FALSE]) %*% mEta[t-1,,drop = FALSE] +
      dB * aQ[,,t-1]
    
    ##Compute the correlation
    aCor[,,t] = diag(sqrt(1/diag(aQ[,,t]))) %*% aQ[,,t] %*% diag(sqrt(1/diag(aQ[,,t])))
    
    ##Augment Likelihod
    dLLk = dLLK + mEta[t,,drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t, , drop = FALSE]) -
      mEta[t, , drop = FLASE] %*% t(mEta[t, , drop = FALSE]) + log(det(aCor[,, t]))
  }
  
  lOut = list()
  lOut[["dLLK"]] = -0.5 * dLKK
  lOut[["aCor"]] = aCor
  
}

##Fucntion to estimate the DCC model
# y_t = Sigma_t^{1/2}z_t
# Sigma_t = D_t^{1/2} R_t D_t^{1/2}
# where D_t is a diagonal matrix with
# typical element D_iit = sigma_it^2.
# we define with eta_t = D_t^{-1/2} y_t
Estimate_DCC = function(mY){
  
  ##estimate marginal GARCH models##
  require(rugarch)
  require(Rsolnp)
  
  SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = 
                           list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "norm")
  
  mspec = multispec(replicate(ncol(mY), SpecGARCH))
  
  fitlist = multifit(multispec = mspec, data = mY)
  
  #Compute residuals
  mEta = as.numeric(residuals(fitlist, standardize = TRUE))
  
  ##maximize DCC likelihood##
  
  #inital parameters
  vPar = c(0.04, 0.9)
  
  #unconditional correlation
  
  mQ = cor(mEta)
  
  #maximize DCC likelihood
  
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ){
    
    Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)
  }, ineqfun = function(vPar, ...){
    sum(vPar)
  }, ineqLB =  1e-14, ineqUB = 0.999,
  LB = c(1e-14,1e-14), UB = c(0.999, 0.999),
  mEta = mEta, mQ = mQ)
  
  #Extract the estimated parameters#
  vPar = optimizer$pars
  
  #Filter the dynmaic correlations using the estimated parameters
  Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
  
  #extract univarative volatilities
  mSigma = 


  
  
}






