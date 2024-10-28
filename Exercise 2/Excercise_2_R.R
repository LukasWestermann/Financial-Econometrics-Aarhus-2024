######Computational part#####

#####libaries####
library(gmm) #for gmm estimation of SV model parameters
library(quantmod) # for S&P data
source("help_functions.R")


#####Simulation of SV model####

set.seed(123)

y = simulate_SV(1000,0, 0.9, 0.25)

#Plot the simulated data
plot(y,type = "l",main = "Simulated returns", xlab = "Time", ylab = "Returns")


Fit_GMM = GMM_Estimator(y)

Fit_GMM$par

########################################################
#                         EX 3                         #
########################################################

library(quantmod)

## download the series
mPrice = getSymbols("^GSPC", from = "2005-01-01", to = "2018-01-01", auto.assign = FALSE)

## compute log returns
vR = as.numeric((diff(log(mPrice[, 6]))*100)[-1])

## replace zeros with the empirical mean
vR[vR == 0] = mean(vR)

## estimate the model by GMM (may take a bit of time)
Fit_GMM = GMM_Estimator(vR)

## see estimated parameters
Fit_GMM$par

# #Plot real vs fitted values

fitted_returns = simulate_SV(length(vR),Fit_GMM$par[1],Fit_GMM$par[2],Fit_GMM$par[3])
plot(as.numeric(vR), type = "l", col = "blue")
lines(fitted_returns,type = "l", col = "red")
legend("topleft",legend = c("Actual Returns", "Estimated Returns"), fill = c("blue","red"),cex = 0.1,lwd = 3)







# # Define moment conditions
# 
# #Initial parameter guess
# theta_start = c(0, 0.5, 0.1)
# 
# gmm_model = gmm(moment_conditions_SV, y, t0 = theta_start)
# summary(gmm_model)
# 
# #Plot actual against fitted values
# fitted_returns = simulate_SV(1000,coef(gmm_model)[1],coef(gmm_model)[2],coef(gmm_model)[3])
# plot(y, type = "l", col = "blue")
# lines(fitted_returns,type = "l", col = "red")
# 
# ####Real Data Application####
# 
# #Download S&P 500 data
# getSymbols("^GSPC", from = "2005-01-01", to = "2018-01-01")
# sp500 = Cl(GSPC) #closing returns
# 
# #Compute log returns
# log_returns = diff(log(sp500))[-1] * 100
# #Replace zero returns with their empirical mean (to ensure that we can take exp)
# log_returns[log_returns == 0] = mean(log_returns)
# 
# #Estimate SV model using GMM
# gmm_model_real = gmm(moment_conditions_SV, log_returns, t0 = c(0,0.9,0.25),optfct = "optim")
# summary(gmm_model_real)
# 
# #Plot real vs fitted values
# 
# fitted_returns = simulate_SV(length(log_returns),coef(gmm_model_real)[1],coef(gmm_model_real)[2],coef(gmm_model_real)[3])
# plot(as.numeric(log_returns), type = "l", col = "blue")
# lines(fitted_returns,type = "l", col = "red")
# legend("topleft",legend = c("Actual Returns", "Estimated Returns"), fill = c("blue","red"),cex = 0.25)
# 
# 
# 
# 
# 
