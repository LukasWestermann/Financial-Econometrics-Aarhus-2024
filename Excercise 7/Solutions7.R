# Estimate the Patton's (2006) model

#Filter for the Patton (2006) model 
# mU is a iT x iN matrix with estimated PIT
# CopType is either "t" or "norm" for the t and normal copula, respectively.
# dOmeha, dAlpha, dBeta are the parameters driving the 
# copula correlation dynamic. dNu is the degree of freedom patameter
# in the t copula case. 
PattonFilter <- function(mU, CopType, dOmega, dAlpha, dBeta, dNu) {
  
  #lambda function
  LambdaFun <- function(x) {
    (1 - exp(-x))/(1 + exp(-x))
  }
  
  ##transforming the PIT
  if (CopType == "norm") {
    mU_Tr = qnorm(mU)
  }
  if (CopType == "t") {
    mU_Tr = qt(mU, dNu)
  }
  
  iT = nrow(mU)
  
  #initialize the correlation dynamic
  vCor = numeric(iT)
  vCor[1] = cor(mU_Tr[, 1], mU_Tr[, 2]) #empirical one
  
  #compute the first likelihood contribution
  if (CopType == "norm") {
    norm.cop <- normalCopula(vCor[1])
    dLLK = dCopula(mU[1, ], norm.cop, log = TRUE)
  }
  if (CopType == "t") {
    t.cop = tCopula(vCor[1], df = dNu)
    dLLK = dCopula(mU[1, ], t.cop, log = TRUE)
  }
  
  # main loop
  for (t in 2:iT) {
    
    #update the correlation using the Patton's (2006) recursion
    vCor[t] = LambdaFun(dOmega + dBeta * vCor[t - 1] + dAlpha * mU_Tr[t - 1, 1] * mU_Tr[t - 1, 2])
    
    #compute the lileihood contribution at time t
    if (CopType == "norm") {
      norm.cop <- normalCopula(vCor[t])
      dLLK = dLLK + dCopula(mU[t, ], norm.cop, log = TRUE)
    }
    if (CopType == "t") {
      t.cop = tCopula(vCor[t], df = dNu)
      dLLK = dLLK + dCopula(mU[t, ], t.cop, log = TRUE)
    }
    
  }
  
  #output the result
  lOut = list()
  lOut[["dLLK"]] = dLLK
  lOut[["vCor"]] = vCor
  
  return(lOut)
}

# Main function to estimate the Patton model with t or normal copula
# mY are the observations. The univariate models are GARCH with t distributed errors
Estimate_Patton <- function(mY, CopType) {
  
  ## estimate the marginal GARCH models
  require(rugarch) #load the rugarch package
  require(copula)
  require(Rsolnp)
  
  #specify the univariate GARCH models
  SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0, 0)), distribution.model = "std")
  
  #Estimate the univariate GARCH models
  lFit_univariate = list()
  
  for(n in 1:ncol(mY)) {
    lFit_univariate[[n]] = ugarchfit(SpecGARCH, mY[, n])
  }
  
  # Extract the PIT
  mU = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(pit(Fit))
  }))
  
  ## robustify
  mU[mU > 0.999] = 0.999
  mU[mU < 0.001] = 0.001
  
  ## maximization of the Copula likelihood in the t and norm copula cases
  # here we impose the constraint alpha + beta < 1 to avoid explosive patterns
  # in the correlation parameter
  if (CopType == "norm") {
    # approximated unc cor
    dCor_app = cor(mU[, 1], mU[, 2]) * 0.16
    # approximated unmapped unc cor
    dOmega_starting = log(dCor_app + 1) - log(1 - dCor_app) 
    vPar = c(dOmega_starting, 0.04, 0.8)
    
    optimizer = solnp(vPar, fun = function(vPar, mU) {
      
      Filter = PattonFilter(mU, CopType = "norm", vPar[1], vPar[2], vPar[3], dNu = NA)
      dNLLK = -as.numeric(Filter$dLLK)
      
      if (!is.finite(dNLLK)) {
        dNLLK = 1e4
      }
      
      if (!is.numeric(dNLLK)) {
        dNLLK = 1e4
      }
      
      return(dNLLK)
      
    }, 
    LB = c(-3, -0.999, 1e-4), UB = c(3, 0.999, 0.9999), 
    mU = mU)
  }
  
  if (CopType == "t") {
    ##unmap initial value
    # approximated unc cor
    dCor_app = cor(mU[, 1], mU[, 2]) * 0.16
    # approximated unmapped unc cor
    dOmega_starting = log(dCor_app + 1) - log(1 - dCor_app) 
    vPar = c(dOmega_starting, 0.04, 0.8, 5)
    
    optimizer = solnp(vPar, fun = function(vPar, mU) {
      
      Filter = PattonFilter(mU, CopType = "t", vPar[1], vPar[2], vPar[3], dNu = vPar[4])
      dNLLK = -as.numeric(Filter$dLLK)
      
      if (!is.finite(dNLLK)) {
        dNLLK = 1e4
      }
      
      if (!is.numeric(dNLLK)) {
        dNLLK = 1e4
      }
      
      return(dNLLK)
      
    },  
    LB = c(-3, -0.999, 1e-4, 2.01), UB = c(3, 0.999, 0.9999, 30), 
    mU = mU)
  }
  
  vPar = optimizer$pars
  dLLK_C = -tail(optimizer$values, 1)
  
  # compute the filtered correlation parameter
  if (CopType == "norm") {
    Filter = PattonFilter(mU, CopType = "norm", vPar[1], vPar[2], vPar[3], dNu = NA)
  }
  if (CopType == "t") {
    Filter = PattonFilter(mU, CopType = "t", vPar[1], vPar[2], vPar[3], dNu = vPar[4])
  }
  
  #extract univariate volatilities
  mSigma = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(sigma(Fit))
  }))
  
  #extract univariate estimated parameters
  mCoef = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(coef(Fit))
  }))
  
  #compute the likelihood od the univariate models
  dLLK_V = do.call(sum, lapply(lFit_univariate, function(Fit) {
    as.numeric(likelihood(Fit))
  }))
  
  #compute the total likelihood of the model
  dLLK = dLLK_V + dLLK_C
  
  if (CopType == "norm") {
    iK = 9 # 3 pars for each marginal + 3 for the copula
  }
  if (CopType == "t") {
    iK = 10 # 3 pars for each marginal + 4 for the copula
  }
  
  iT = nrow(mY)
  
  BIC = log(iT) * iK - 2 * dLLK
  
  lOut = list()
  
  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["vCor"]] = Filter[["vCor"]]
  lOut[["BIC"]] = BIC
  
  return(lOut)
  
}

##############################################################
#                   Estimation part                          #
##############################################################

library(rugarch)
library(copula)
data("dji30ret")
mY = dji30ret[1:1000, 1:2]

#it is time consuming
Fit_Patton_t = Estimate_Patton(mY, CopType = "t")
#it is time consuming
Fit_Patton_norm = Estimate_Patton(mY, CopType = "norm")

plot.ts(Fit_Patton_t$vCor)
lines(Fit_Patton_norm$vCor, col = "red")
# The two are very similar

vBIC = c(
  BIC_Patton_Gauss = Fit_Patton_norm$BIC,
  BIC_Patton_t     = Fit_Patton_t$BIC
)

sort(vBIC)

# The model selected by BIC is Student t model

