#Log-Likelihood fucntion of DCC model (for second step,correlation component)

log_likelihood_dcc_normal = function(params, residuals){
  alpha = params[1]
  beta = params[2]
  T = nrow(residuals)
  Q_bar = cor(residuals) #uncondtional correlation
  Q_t = Q_bar #intialize correlation matrix
  log_likelihhod = 0 #intialize log likelihood
  
  for(t in 2:T){
    Q_t = (1-alpha-beta) * Q_bar + alpha * residuals[t-1, ] %*% t(residuals[t-1, ]) + beta * Q_t
    D_t = diag(sqrt(diag(Q_t)))
    R_t = solve(D_t) %*% Q_t %*% solve(D_t)
    log_likelihhod = log_likelihhod - 0.5 *( log(det(R_t)) + t(residuals[t, ]) %*% solve(R_t) %*% residuals[t, ] 
                                             - t(residuals[t, ]) %*% residuals[t, ] )
  }# use log fucntio of det to calcuate the log of a dtermindent it is more robust!!!
  return(-log_likelihhod) #return negative log-likehood for minimzatioon
}


log_likelihood_dcc_t = function(params,residuals){
  alpha = params[1]
  beta = params[2]
  nu = params[3]
  T = nrow(residuals)
  Q_bar = cor(residuals) #uncondtional correlation
  Q_t = Q_bar #intialize correlation matrix
  log_likelihhod = 0 #intialize log likelihood
  
  for(t in 2:T){
    Q_t = (1-alpha-beta) * Q_bar + alpha * residuals[t-1, ] %*% t(residuals[t-1, ]) + beta * Q_t
    D_t = diag(sqrt(diag(Q_t)))
    R_t = solve(D_t) %*% Q_t %*% solve(D_t)
    term_1 = lgamma((nu + ncol(residuals)) / 2) - lgamma(nu / 2) - (ncol(residuals) / 2) * log(pi * (nu - 2))
      - 0.5 * log(det(R_t))
    term_2 = -(nu + ncol(residuals)) / 2 * log(1 + t(residuals[t, ]) %*% solve(R_t) %*% residuals[t, ])
    log_likelihhod = log_likelihhod + term_1 + term_2
  }
  return(-log_likelihhod)
}
