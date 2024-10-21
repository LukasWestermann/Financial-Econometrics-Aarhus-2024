##libraries##
library(quantmod)
library(moments)
library(tseries) 
library(rugarch)
#Import helper functions
source("help_functions.R")

##1##
DJITicker = read.csv("DJITicker.csv", sep = ";")

##2##

# List of tickers (assuming 'dji_tickers$Ticker' contains the tickers)
ticker <- DJITicker$Symbol

# Use lapply to get the log returns for each ticker
lRet <- lapply(ticker, get_log_returns)

valid_ticker <- ticker[!sapply(lRet, is.null)]
lRet <- lRet[!sapply(lRet, is.null)]

# Name the list elements with the valid tickers
names(lRet) <- valid_ticker



##3##

DescStat = matrix(nrow = 28, ncol = 7)
rownames(DescStat) = valid_ticker
colnames(DescStat) = c("mean", "median", "variance", "kurtosis", "skewness", "rho", "rho2")


#Fill the rows

for(i in 1:length(valid_ticker)){
  returns = lRet[[valid_ticker[i]]]
  DescStat[i, "mean"] = mean(returns, na.rm = TRUE)
  DescStat[i, "median"] = median(returns, na.rm = TRUE)
  DescStat[i, "variance"] = var(returns, na.rm = TRUE)
  DescStat[i, "kurtosis"] = kurtosis(returns, na.rm = TRUE)
  DescStat[i, "skewness"] = skewness(returns, na.rm = TRUE)
  DescStat[i, "rho"] = acf(returns, lag.max = 1, 
                            plot = FALSE, na.action = na.pass)[[1]][2]
  DescStat[i, "rho2"] = acf(returns^2, lag.max = 1, 
                            plot = FALSE, na.action = na.pass)[[1]][2]
}

print(DescStat)


##4##

# Download the S&P 500 index data
getSymbols("^GSPC", scr = "yahoo", from = "2023-01-01", to = "2023-12-31")

# Compute S&P 500 returns
sp500_returns <- diff(log(Ad(GSPC))) * 100
sp500_returns = sp500_returns[-(1:2)] #to fit length of ticker data

#Initalize CAPM matrix
CAPM_results = matrix(NA,nrow = length(valid_ticker), ncol = 3)
rownames(CAPM_results) = valid_ticker
colnames(CAPM_results) = c("alpha","beta","sigma2")

#Estimate CAPM for each asset

for(i in 1:length(valid_ticker)){
  asset_returns = lRet[[i]][-1]
  capm_model = lm(asset_returns~sp500_returns)
  
  # Extract coefficients and residual variance
  CAPM_results[i, "alpha"] = coef(capm_model)[1]
  
  CAPM_results[i, "alpha"] = coef(capm_model)[1]
  CAPM_results[i, "beta"] = coef(capm_model)[2]
  CAPM_results[i, "sigma2"] = var(residuals(capm_model))
}

##5 Simulation##

#see help_functions.R

##6 Estimation##

#see help_functions.R

#Test
garch_simulation <- simulate_garch(100, 0.1, 0.3, 0.6)
arch_estimate <- arch_mle(garch_simulation$y, c(0.1,0.1))
print(arch_estimate)

##7 Monte Carlo Simulations##


## Monte Carlo for T = 200, 500, 1000 for correctly specified model
B <- 500
omega_true <- 0.3
alpha_true <- 0.7
beta_true = 0

mc_results <- array(NA, dim = c(B, 2, 3))
dimnames(mc_results) <- list(NULL, c("omega", "alpha"), c("200", "500", "1000"))

# Run Monte Carlo simulations for different values of T
set.seed(123)
mc_results[, , 1] <- peform_montecarlo(200, omega_true, alpha_true,beta_true, B)
mc_results[, , 2] <- peform_montecarlo(500, omega_true, alpha_true,beta_true, B)
mc_results[, , 3] <- peform_montecarlo(1000, omega_true, alpha_true,beta_true, B)

A = peform_montecarlo(1000, omega_true, alpha_true,beta_true, B)
# Example: display the first few results for T = 200
print(head(mc_results[, , 1]))  # Results for T = 200

# Plot densities of omega and alpha for T = 200, 500, and 1000
par(mfrow = c(2, 3))

# Densities for ω
plot(density(mc_results[, "omega", 1]), main = "Density of ω (T = 200)", xlab = "ω")
abline(v = 0.3)
plot(density(mc_results[, "omega", 2]), main = "Density of ω (T = 500)", xlab = "ω")
abline(v = 0.3)
plot(density(mc_results[, "omega", 3]), main = "Density of ω (T = 1000)", xlab = "ω")
abline(v = 0.3)

mean(mc_results[, "omega", 1]-0.3)
mean(mc_results[, "omega", 2]-0.3)
mean(mc_results[, "omega", 3]-0.3)

mean((mc_results[, "omega", 1]-0.3)^2)
mean((mc_results[, "omega", 2]-0.3)^2)
mean((mc_results[, "omega", 3]-0.3)^2)


# Densities for α
plot(density(mc_results[, "alpha", 1]), main = "Density of α (T = 200)", xlab = "α")
abline(v = 0.7)
plot(density(mc_results[, "alpha", 2]), main = "Density of α (T = 500)", xlab = "α")
abline(v = 0.7)
plot(density(mc_results[, "alpha", 3]), main = "Density of α (T = 1000)", xlab = "α")
abline(v = 0.7)

mean(mc_results[, "alpha", 1]-0.7)
mean(mc_results[, "alpha", 2]-0.7)
mean(mc_results[, "alpha", 3]-0.7)

mean((mc_results[, "alpha", 1]-0.3)^2)
mean((mc_results[, "alpha", 2]-0.3)^2)
mean((mc_results[, "alpha", 3]-0.3)^2)


#Monte Carlo Simulation for misspecified model
B <- 500
omega_true <- 0.3
alpha_true <- 0.1
beta_true = 0.8

mc_results_miss <- array(NA, dim = c(B, 2, 3))
dimnames(mc_results_miss) <- list(NULL, c("omega", "alpha"), c("200", "500", "1000"))

# Run Monte Carlo simulations for different values of T
set.seed(123)
mc_results_miss[, , 1] <- peform_montecarlo(200, omega_true, alpha_true,beta_true, B)
mc_results_miss[, , 2] <- peform_montecarlo(500, omega_true, alpha_true,beta_true, B)
mc_results_miss[, , 3] <- peform_montecarlo(1000, omega_true, alpha_true,beta_true, B)

# Example: display the first few results for T = 200
print(head(mc_results_miss[, , 1]))  # Results for T = 200

# Plot densities of omega and alpha for T = 200, 500, and 1000
par(mfrow = c(2, 3))

# Densities for ω
plot(density(mc_results_miss[, "omega", 1]), main = "Density of ω (T = 200)", xlab = "ω")
abline(v = 0.3)
plot(density(mc_results_miss[, "omega", 2]), main = "Density of ω (T = 500)", xlab = "ω")
abline(v = 0.3)
plot(density(mc_results_miss[, "omega", 3]), main = "Density of ω (T = 1000)", xlab = "ω")
abline(v = 0.3)

# Densities for α
plot(density(mc_results_miss[, "alpha", 1]), main = "Density of α (T = 200)", xlab = "α")
abline(v = 0.1)
plot(density(mc_results_miss[, "alpha", 2]), main = "Density of α (T = 500)", xlab = "α")
abline(v = 0.1)
plot(density(mc_results_miss[, "alpha", 3]), main = "Density of α (T = 1000)", xlab = "α")
abline(v = 0.1)




##8 Real Data Analysis##
# Install the rugarch package (if needed)


# Function to estimate ARCH(1), GARCH(1,1), E-GARCH(1,1), and GJR-GARCH(1,1)
estimate_garch_models <- function(returns) {
  # Define model specifications
  spec_arch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 0)),
                          mean.model = list(armaOrder = c(0, 0)))
  
  spec_garch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(0, 0)))
  
  spec_egarch <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                            mean.model = list(armaOrder = c(0, 0)))
  
  spec_gjr <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)))
  
  # Fit models
  fit_arch <- ugarchfit(spec_arch, returns)
  fit_garch <- ugarchfit(spec_garch, returns)
  fit_egarch <- ugarchfit(spec_egarch, returns)
  fit_gjr <- ugarchfit(spec_gjr, returns)
  
  # Collect BIC values
  bic_values <- c(infocriteria(fit_arch)[2], 
                  infocriteria(fit_garch)[2],
                  infocriteria(fit_egarch)[2],
                  infocriteria(fit_gjr)[2])
  
  names(bic_values) <- c("ARCH(1)", "GARCH(1,1)", "E-GARCH(1,1)", "GJR-GARCH(1,1)")
  return(bic_values)
}

# Apply the function to each series in the DJIA
all_bic <- lapply(lRet, function(returns) {
  estimate_garch_models(na.omit(returns))
})

# Identify the best model for each series
best_models <- sapply(all_bic, function(bic_values) {
  return(names(which.min(bic_values)))
})

# Display the best models for each series
print(best_models)

