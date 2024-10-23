######Computational part#####

#####Simulation of SV model####

set.seed(123)

T = 1000
omega = 0
phi = 0.9
sigma2_eta = 0.25
sigma_eta = sqrt(0.25)

#Simulate latent volatilty process

h = rep(0,T)
h[1] = omega/(1-phi) #intial value, set to expetced value of h

#Loop over recursions
for(t in 2:T){
  h[t] = omega + phi * h[t-1] + rnorm(1,0,sigma_eta)
}

#Simulate returns y_t
y = exp(h/2) * rnorm(T)

#Plot the simulated data
plot(y,type = "l",main = "Simulated returns", xlab = "Time", ylab = "Returns")


