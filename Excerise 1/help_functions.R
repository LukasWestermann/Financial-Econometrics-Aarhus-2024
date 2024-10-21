
#get log returns
get_log_returns = function(ticker){
  tryCatch({
    getSymbols(ticker, scr = "yahoo", from = "2023-01-01", to = "2023-12-31")
    prices = Ad(get(ticker)) #prices
    log_returns = diff(log(na.omit(prices)))[-1] * 100 #Log returns in percentage
    return(log_returns)
  },error = function(e) {
    # If there is an error (e.g., the ticker is delisted), return NULL
    message(paste("Data for", ticker, "not found or delisted."))
    return(NULL)
  })
}

#Simulate GARCH

simulate_garch = function(T,omega,alpha,beta){
  
  sigma2 = rep(NA,T)
  y = rep(NA,T)
  
  #Set initial values
  sigma2[1] = omega/(1-alpha-beta)
  y[1] = rnorm(1,sd = sqrt(sigma2[1]))
  
  #For loop over T
  for(t in 2:T){
    sigma2[t] = omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
    y[t] = rnorm(1,sd = sqrt(sigma2[t]))
  }
  
  return(list(sigma2 = sigma2,y= y))
  
}

#Function to calculate ARCH(1) log-likelihood function, assuming normality

arch_ll = function(params, y){
  #initialize parameters
  omega = params[1]
  alpha = params[2]
  T = length(y)
  
  sigma2 = rep(NA,T)
  sigma2[1] = omega/(1-alpha) # set to unconditional value
  
  sigma2[2:T] = omega + alpha * y[1:(T-1)]^2
  LLK = sum(dnorm(y,0,sqrt(sigma2), log = TRUE))
  return(-LLK)
}

#Function to estimate ARCH(1) by MLE

arch_mle = function(y,start_params){ #start values for omega and alpha
  opt = optim(start_params,arch_ll,method = "L-BFGS-B", 
              lower = c(0.000001,0.000001), upper = c(10,0.99999), y = y)# bounds to ensure stationary
  #use L-BFGS-B because we have a constraint optimization problem
  #can add a bunch of other outputs,e.g LLK, BIC etc.
  return(opt$par)
  #good guess for the intitail parameters
  #omega = var(vY,)
}



#Function to perform Monte carlo simulation

peform_montecarlo = function(T,omega_true, alpha_true, beta_true,B){
  #matrix to store estimates
  estimates = matrix(NA,nrow = B, ncol = 2)
  
  for(b in 1:B){
    #simulate data for ARCH(1) model with true parameters
    simulation = simulate_garch(T, omega_true, alpha_true, beta_true)
    estimates[b,] = arch_mle(simulation$y, c(omega_true, alpha_true))
  }
  return(estimates)
}