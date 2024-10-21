## Load the quantmod package
library(quantmod)

## set the working directory
setwd("/Users/au588008/Dropbox/Teaching/Aarhus_Financial Econometrics/2024/ExerciseSet/1/Solutions/")

## load the csv file
mTicker = read.csv("DJITicker.csv", header = TRUE, sep = ";", 
                   stringsAsFactors = FALSE)

# DWDP and UTX are no more valid tickers
mTicker = mTicker[-c(10, 29), ]

## create the list where prices will be stored
lPrices = list()

## loop over the elements of the first column of mTicker
for (sTicker in mTicker[, "Symbol"]) {
  
  ## download the price series of the ticker 
  ## Symbols = sTicker
  ## start = from "2022-01-01" 
  ## end = to "2024-01-01"
  ## auto.assign = FALSE indicates that the result should be put in lPrices
  lPrices[[sTicker]] = getSymbols(Symbols = sTicker,
                                  from = "2022-01-01",
                                  to = "2024-01-01",
                                  auto.assign = FALSE)
  
}

## create the list of where we store percentage log returns
lRet = list()

## loop over the elements of the first column of mTicker
for (sTicker in mTicker[, "Symbol"]) {
  
  ## compute the percentage returns:
  ## na.omit removes NAs
  ## log does the logarithm of the prices
  ## diff returns the first difference
  ## [-1] is used in oder to remove the NA resulting from diff
  ## since we don't have the price at time 0 diff return an NA at the 
  ## first position
  ## the output is moltiplied by 100 in order to have percentage returns
  lRet[[sTicker]] = diff(log(na.omit(lPrices[[sTicker]][, 6])))[-1] * 100
  
}

## create the matrix where to store descriptive statistics
DescStat = matrix(NA, 
                  nrow = nrow(mTicker),
                  ncol = 7,
                  ## the names of the rows and columns are defined with the dimnames argument
                  dimnames = list(mTicker[, "Symbol"],
                                  c("mean", "median", "variance", "kurtosis", "skewness", "rho", "rho2"))
                  )

## function to compute the kurtosis coefficient
kurtosis <- function(vY) {
  mean((vY - mean(vY))^4)/(mean((vY - mean(vY))^2)^2)
}

## function to compute the skewness coefficient
skewness <- function(vY) {
  mean((vY - mean(vY))^3)/(mean((vY - mean(vY))^2)^(3/2))
}

## loop over the elements of the first column of mTicker
for (sTicker in mTicker[, "Symbol"]) {
  DescStat[sTicker, "mean"] = mean(lRet[[sTicker]])
  DescStat[sTicker, "median"] = median(lRet[[sTicker]])
  DescStat[sTicker, "variance"] = var(lRet[[sTicker]])
  DescStat[sTicker, "skewness"] = skewness(lRet[[sTicker]])
  DescStat[sTicker, "kurtosis"] = kurtosis(lRet[[sTicker]])
  
  ## here we exploit the acf function in order to compute the first order autocorrelation
  DescStat[sTicker, "rho"] = acf(lRet[[sTicker]], lag.max = 1, plot = FALSE)$acf[2, 1, 1]
  
  ##notice how to extract the autocorrelations: see help(acf)
  DescStat[sTicker, "rho2"] = acf(lRet[[sTicker]]^2, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
}

##function to perform OLS estimation
OLS_Estimator <- function(vY, vX) {
  
  dBeta   = cov(vY, vX)/var(vX)
  dAlpha  = mean(vY) - dBeta * mean(vX) 
  
  dSigma2 = var(vY - dBeta * vX - dAlpha)
  
  vOut = c("beta" = dBeta,
           "alpha" = dAlpha,
           "sigma2" = dSigma2)
  
  return(vOut)
}

## create a matrix of returns
mRet = do.call(cbind, lRet)

## download SP500 returns
getSymbols("^GSPC", from = "2022-01-01", to = "2024-01-01")
## approximate the market index with SP500 returns
vX = as.numeric(diff(log(na.omit(GSPC[, 6])))[-1] * 100)

## create a matrix where to store estimates coefficients
mEst = matrix(NA, nrow = nrow(mTicker), ncol = 3,
              dimnames = list(mTicker[, "Symbol"],
                              c("beta", "alpha", "sigma2")))

## loop over the elements of the first column of mTicker
for (sTicker in mTicker[, "Symbol"]) {
  
  ## insert the output of OLS_Estimator to the row sTicker of mEst
  mEst[sTicker, ] = OLS_Estimator(as.numeric(lRet[[sTicker]]), vX)
  
}

## 5) Simulation from GARCH

GARCHSim <- function(iT, dOmega, dAlpha, dBeta) {
  
  ## initialize the vector of simulated returns and variances
  vY = numeric(iT)
  vSigma2 = numeric(iT)
  
  ## initialize the variance at time t = 1 with its unconditional value
  vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta)
  ## sample the first observations
  vY[1] = rnorm(1, mean = 0, sd = sqrt(vSigma2[1]))
  
  ##loop over iT. We start from t = 2 since t = 1 has already been sampled
  for (t in 2:iT) {
    #update the volatility
    vSigma2[t] = dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1]
    #sample a new observarions
    vY[t] = rnorm(1, mean = 0, sd = sqrt(vSigma2[t]))
  }
  
  ## we return a list with two components: the sampled returns and the volatility
  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vSigma2"]] = vSigma2
  
  ## output lOut
  return(lOut)
} 

## 6) Estimation

## Function to evaluate the likelihood of an ARCH(1) model
ARCHLLK <- function(vY, dOmega, dAlpha) {
  
  ## number of observations
  iT = length(vY)
  
  ## initialize the variance
  vSigma2 = numeric(iT)
  
  ## set the variance at time t = 1 at its unconditional value
  vSigma2[1] = dOmega/(1.0 - dAlpha)
  
  # ## compute the log--likelihood of the first obs
  # dLLK = dnorm(vY[1], 0, sqrt(vSigma2[1]), log = TRUE)
  # 
  # #######################################################################
  # ## The following loop is not the best efficient solution in our case  #
  # ## can you do better?                                                 #
  # #######################################################################
  # 
  # ## loop over iT
  # for (t in 2:iT) {
  #   # update the volatility
  #   vSigma2[t] = dOmega + dAlpha * vY[t-1]^2
  #   # add the log-likelihood contribution at time t
  #   dLLK = dLLK + dnorm(vY[t], 0, sqrt(vSigma2[t]), log = TRUE)
  # }
  
  ## alternative strategy
  vSigma2[2:iT] = dOmega + dAlpha * vY[1:(iT - 1)]^2
  dLLK = sum(dnorm(vY, 0, sqrt(vSigma2), log = TRUE))
  
  # return the likelihood
  return(dLLK)
  
}

## function to estimate an ARCH(1) model

#we first code the objective function i.e. the negative log likelihood
#se specify a vector of coefficients to be estimated, the first coefficient
#is omega, the second is alpha

ObjectiveFunction <- function(vPar, vY) {
  
  dOmega = vPar[1]
  dAlpha = vPar[2]
  dLLK = ARCHLLK(vY, dOmega, dAlpha)
  
  return(-dLLK)
}

## function to estimate the ARCH(1) model
EstimateARCH <- function(vY) {
  
  # We set starting value for alpha equal to 0.8, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # ARCH model
  
  dAlpha = 0.8 # can you find a better starting value for alpha?
  dOmega = var(vY) * (1.0 - dAlpha)
  
  ## vector of starting parameters
  vPar = c(dOmega, dAlpha)
  
  ##optimization step
  optimizer = optim(vPar, fn = ObjectiveFunction, method = "L-BFGS-B", vY = vY,
                    ## note that we set suitable constraints on the model parameters
                    lower = c(0.00001, 0.0001), upper = c(10.0, 0.999)) 
  
  ## extract estimated parameters
  vPar = optimizer$par
  
  ## extract the likelihood computed at its maximum
  dLLK = -optimizer$value
  
  ## Compute the Average BIC
  iT = length(vY)
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT
  
  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC
  
  return(lOut)
}

## 7) Monte Carlo

#set the seed in order to get the same results
set.seed(123)

## number of replicates
iB = 500
## different sample sizes
vT = c(200, 500, 1000)

## We all estimated results in an array
aCoef = array(NA, dim = c(iB, length(vT), 2, 2),
              dimnames = list(NULL, vT, c("omega", "alpha"), c("Correct", "Misspecified")))

## correctly specified model
for (iT in vT) { #loop over the sample sized
  for (b in 1:iB) { #loop over the replicates
    
    #simulate an ARCH model
    lSim = GARCHSim(iT, dOmega = 0.3, dAlpha = 0.7, dBeta = 0.0)
    #estimate an ARCH model
    Fit = EstimateARCH(lSim[["vY"]])
    #collect the estimated coefficients in the array
    aCoef[b, paste(iT), ,"Correct"] = Fit$vPar
    
  }
}

## misspecified model
for (iT in vT) { #loop over the sample size
  for (b in 1:iB) { #loop over the replicates
    
    #simulate an ARCH model
    lSim = GARCHSim(iT, dOmega = 0.3, dAlpha = 0.1, dBeta = 0.8)
    #estimate an ARCH model
    Fit = EstimateARCH(lSim[["vY"]])
    #collect the estimated coefficients in the array
    aCoef[b, paste(iT), ,"Misspecified"] = Fit$vPar
    
  }
}

## Plot the density of the estimated parameters
## we organize the plots in a 2x2 grid

## create a 2 x 2 grid of plot
par(mfrow = c(2, 2))

## plot omega in the correctly specified case
plot(density(aCoef[, "200", "omega", "Correct"]), col = "blue", main = "Omega - Correct Spec",
     xlab = "", ylab = "", lwd = 2, ylim = c(0, 18))

lines(density(aCoef[, "500", "omega", "Correct"]), col = "red", lwd = 2)
lines(density(aCoef[, "1000", "omega", "Correct"]), col = "purple", lwd = 2)
## add the true value
abline(v = 0.3)

## plot alpha in the correctly specified case
plot(density(aCoef[, "200", "alpha", "Correct"]), col = "blue", main = "Alpha - Correct Spec",
     xlab = "", ylab = "", lwd = 2, ylim = c(0, 6))

lines(density(aCoef[, "500", "alpha", "Correct"]), col = "red", lwd = 2)
lines(density(aCoef[, "1000", "alpha", "Correct"]), col = "purple", lwd = 2)
## add the true value
abline(v = 0.7)

## plot omega in the misspecified case
plot(density(aCoef[, "200", "omega", "Misspecified"]), col = "blue", main = "Omega - Misspecified Spec",
     xlab = "", ylab = "", lwd = 2, ylim = c(0, 2.0), xlim = c())

lines(density(aCoef[, "500", "omega", "Misspecified"]), col = "red", lwd = 2)
lines(density(aCoef[, "1000", "omega", "Misspecified"]), col = "purple", lwd = 2)
## add the true value
abline(v = 0.3)

## plot alpha in the misspecified specified case
plot(density(aCoef[, "200", "alpha", "Misspecified"]), col = "blue", main = "Alpha - Misspecified Spec",
     xlab = "", ylab = "", lwd = 2, ylim = c(0, 7.5))

lines(density(aCoef[, "500", "alpha", "Misspecified"]), col = "red", lwd = 2)
lines(density(aCoef[, "1000", "alpha", "Misspecified"]), col = "purple", lwd = 2)
## add the true value
abline(v = 0.1)

## 8) Real Data

## load quantmod
library(quantmod)

## load the csv file
mTicker = read.csv("DJITicker.csv", header = TRUE, sep = ";", 
                   stringsAsFactors = FALSE)

# DWDP and UTX are no more valid tickers
mTicker = mTicker[-c(10, 29), ]

## create the list where prices will be stored
lPrices = list()

## loop over the elements of the first column of mTicker
for (sTicker in mTicker[, "Symbol"]) {
  
  ## download the price series of the ticker 
  ## Symbols = sTicker
  ## start from "2017-01-01" end today
  ## auto.assign = FALSE indicates that the result should be put in lPrices
  lPrices[[sTicker]] = getSymbols(Symbols = sTicker,
                                  from = "2017-01-01",
                                  to = "2020-01-01",
                                  auto.assign = FALSE)
  
}

## create the list of where we store percentage log returns
lRet = list()

## loop over the elements of the first column of mTicker
for (sTicker in mTicker[, "Symbol"]) {
  
  ## compute the percentage returns:
  ## na.omit removes NAs
  ## log does the logarithm of the prices
  ## diff returns the first difference
  ## [-1] is used in oder to remove the NA resulting from diff
  ## since we don't have the price at time 0 diff return an NA at the 
  ## first position
  ## the output is moltiplied by 100 in order to have percentage returns
  lRet[[sTicker]] = as.numeric(diff(log(na.omit(lPrices[[sTicker]][, 6])))[-1] * 100)
  
  ## here we use as.numeric since we don't care about the time index
}

## We store all the BIC in a martrix
mBIC = matrix(NA, length(lRet), 4, 
              dimnames = list(names(lRet), c("ARCH", "GARCH", "EGARCH", "GJRGARCH")))


## Estimate ARCH and store the BIC
for (sTicker in mTicker[, "Symbol"]) {
  
  ## Fir ARCH
  Fit = EstimateARCH(lRet[[sTicker]])
  ## store the BIC
  mBIC[sTicker, "ARCH"] = Fit$BIC
  
}

## install the rugarch package
# install.packages("rugarch")

#load the rugarch package
library(rugarch)

## Specify the three GARCH models
#GARCH
GARCHSpec = ugarchspec(variance.model = list(model = "sGARCH"), 
                       mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
#EGARCH
EGARCHSpec = ugarchspec(variance.model = list(model = "eGARCH"), 
                        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
#GJRGARCH
GJRGARCHSpec = ugarchspec(variance.model = list(model = "gjrGARCH"), 
                          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))

#### Store the estimated models

lFit = list()
lFit[["GARCH"]] = list()
lFit[["EGARCH"]] = list()
lFit[["GJRGARCH"]] = list()

## Estimate the models and store the BIC
for (sTicker in mTicker[, "Symbol"]) {
  
  ## GARCH
  Fit_GARCH = ugarchfit(GARCHSpec, lRet[[sTicker]])
  ## EGARCH
  Fit_EGARCH = ugarchfit(EGARCHSpec, lRet[[sTicker]])
  ## GJRGARCH
  Fit_GJRGARCH = ugarchfit(GJRGARCHSpec, lRet[[sTicker]])
  
  ## extract the BIC
  mBIC[sTicker, "GARCH"] = infocriteria(Fit_GARCH)[2]
  mBIC[sTicker, "EGARCH"] = infocriteria(Fit_EGARCH)[2]
  mBIC[sTicker, "GJRGARCH"] = infocriteria(Fit_GJRGARCH)[2]
  
  ## GARCH
  lFit[["GARCH"]][[sTicker]] = Fit_GARCH
  ## EGARCH
  lFit[["EGARCH"]][[sTicker]] = Fit_EGARCH
  ## GJRGARCH
  lFit[["GJRGARCH"]][[sTicker]] = Fit_GJRGARCH
  
}

## Select the best model for each asset

vBest = apply(mBIC, 1, which.min)

vBest[] = colnames(mBIC)[vBest]

Selection = as.data.frame(vBest)

Selection
