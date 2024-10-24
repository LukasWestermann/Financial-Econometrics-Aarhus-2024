####function to simulate data from a SV model (from lecture 6)####
simulate_SV = function(T,omega,phi,sigma2_eta){
  set.seed(123)
  sigma2_eta = sqrt(sigma2_eta)
  h = numeric(T)
  #Initial values
  h[1] = omega/(1-phi) #intial value, set to expetced value of h
  
  #Loop over T
  for(t in 2:T){
    h[t] = omega + phi * h[t - 1] + rnorm(1,mean = 0,sd = sigma_eta)
  }
  y = exp(h/2) * rnorm(T)
  return(y)
}
  
  

T = 1000
omega = 0
phi = 0.9
sigma2_eta = 0.25
sigma_eta = sqrt(sigma2_eta)

#Simulate latent volatilty process

h = numeric(T)
h[1] = omega/(1-phi) #intial value, set to expetced value of h

#Loop over recursions
for(t in 2:T){
  h[t] = omega + phi * h[t - 1] + rnorm(1,mean = 0,sd = sigma_eta)
}


####Function to calculate the moment conditions of the SV model of lecture 6####

moment_conditions_SV = function(theta, y){
  #fucntion to comppute the (empirical) moment conditions
  #theta vector of parameters
  #y data input
  
  omega <- theta[1]
  phi <- theta[2]
  sigma2_eta <- theta[3]
  
  # Reparameterization
  alpha <- omega / (1 - phi)
  beta_2 <- sigma2_eta / (1 - phi^2)
  
  
  #calculate moments
  m1 = sqrt(2/pi) * exp(alpha/2 + beta_2/8)
  m2 = exp(alpha  + beta_2/2)
  m3 = 2 * sqrt(2/pi) * exp( 3 * alpha/2 + 9 * beta_2/8)
  m4 = 3 * exp( 2 * alpha + 2* beta_2)
  
  
  #calculate sample moments and take difference
  sample_moments = c(mean(abs(y)), mean(y^2), mean(abs(y)^3), mean(y^4))
  moment_diff = sample_moments - c(m1, m2, m3, m4)
  
  #go over the moment condtions from 1:10 for the abs avule of E[y_t, y_t-1] and calculate emptical means
  function_diff_1 = function(j){
    mean(abs(y*lag(y, k = j)), na.rm = TRUE) - (sqrt(2/pi) * 2*exp(alpha/2 + beta_2/8) * 
                                                  exp(phi^j*beta_2/4))
  }
  autocov_abs_diff = sapply(1:10, function_diff_1)
  
  #go over the moment condtions from 1:10 for E[y_t^2, y_t-1^^2] and calculate emprical means
  function_diff_2 = function(j){
    mean(y^2*lag(y, k = j)^2 , na.rm = TRUE) - (exp(alpha  + beta_2/2)^2 * exp(phi^j*beta_2))
  }
  
  autocov_sq_diff = sapply(1:10, function_diff_2)
  
  #Combine all 24 moment conditions
  return(c(moment_diff, autocov_abs_diff, autocov_sq_diff))
  
}