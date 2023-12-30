# Author: Jannik Steenbergen

####################### Implementation of TVP-LPs ##############################
KF_tvplp <- function(params, y, X, robust){
  # Description
  # a function defining the tvp-lp in state space representation and running the Kalman Filter on it
  # Arguments
  # y: nx1 vector of dependent variable obs
  # X: nxm matrix of regressors
  # params: (m+1)x1 vector of parameters ordered as c(sigma2_h, sigma2_eta_1, ..., sigma2_eta_m)
  # robust: If robust == TRUE, autocorrelation correction by Bryson and Henrikson (1968) is applied. 
  # If robust == TRUE, params becomes (m+2)x1 with c(sigma2_h, sigma2_eta_1, ..., sigma2_eta_m, phi)
  
  m <- ncol(X)    # number of state variables
  n <- nrow(X)    # sample size
  n_params <- length(params) # number of parameters
  
  # define parameters
  
  if(robust == TRUE){
    sigma2_h <- params[1]
    sigma2_eta_vec <- params[2:(n_params-1)]
    phi <- params[n_params]
  } else {
    sigma2_h <- params[1]
    sigma2_eta_vec <- params[-1]
  }
  
  # Define state space parameters
  if (robust == TRUE){
    y = y[-1] - phi* y[-n]
    Z = X[-1,] - phi*X[-n,]
  } else {
    y = y
    Z = X
  }
  T = diag(m)
  R = diag(m)
  d = 0
  c = 0
  H = sigma2_h
  Q = diag(sigma2_eta_vec)
  
  # initialization
  initialization <- vector("list", m)
  initialization[1:m] <- "diffuse"
  
  # aim 
  aim <- "p"
  
  KF_output <- KF(y=y, Z=Z, T=T, R=R, d=d, c=c, H=H, Q=Q, initialization=initialization, aim=aim)
  
  return(KF_output)
}

KFS_tvplp <- function(params, y, X, robust){
  # Description
  # a function defining the tvp-lp in state space representation and running the Kalman Filter Smoother on it
  # Arguments
  # y: nx1 vector of dependent variable obs
  # X: nxm matrix of regressors
  # params: (m+1)x1 vector of parameters ordered as c(sigma2_h, sigma2_eta_1, ..., sigma2_eta_m)
  # robust: If robust == TRUE, autocorrelation correction by Bryson and Henrikson (1968) is applied. 
  # If robust == TRUE, params becomes (m+2)x1 with c(sigma2_h, sigma2_eta_1, ..., sigma2_eta_m, phi)
  
  m <- ncol(X)    # number of state variables
  n <- nrow(X)    # sample size
  n_params <- length(params) # number of parameters
  
  # define parameters
  if(robust == TRUE){
    sigma2_h <- params[1]
    sigma2_eta_vec <- params[2:(n_params-1)]
    phi <- params[n_params]
  } else {
    sigma2_h <- params[1]
    sigma2_eta_vec <- params[-1]
  }
  
  # Define state space parameters
  if (robust == TRUE){
    y = y[-1] - phi* y[-n]
    Z = X[-1,] - phi*X[-n,]
  } else {
    y = y
    Z = X
  }
  T = diag(m)
  R = diag(m)
  d = 0
  c = 0
  H = sigma2_h
  Q = diag(sigma2_eta_vec)
  
  # initialization
  initialization <- vector("list", m)
  initialization[1:m] <- "diffuse"
  
  # aim 
  aim <- "p"
  
  # run Kalman filter smoother
  KFS_output <- KFS(y=y, Z=Z, T=T, R=R, d=d, c=c, H=H, Q=Q, initialization=initialization)
  
  return(KFS_output)
}

ML_tvplp <- function(y, X, robust){
  # Description: 
  # Function performing ML estimation of a tvp-lp model. The model specification is defined in function KF_tvplp
  # Arguments
  # y: nx1 vector of dependent variable obs
  # X: nxm vector of regressors
  # robust: If robust == TRUE, autocorrelation correction by Bryson and Henrikson (1968) is applied. 
  
  m <- ncol(X)        # number of state variables 
  
  if (robust == TRUE){
    n_params <- m+2     # number of params to be estimated
    
    ## 1: initial parameter estimates
    #initial_params <- rep(0.002, (n_params-1))  # initial variance estimates     # old
    #initial_params <- c(initial_params, 0.5)    # append initial phi-estimate
    initial_params <- rep(0.002, (n_params-2))  # initial variance estimates
    initial_params <- c(0.1, initial_params, 0.5)    # append initial phi-estimate
    
    ## 2: estimation
    
    # bounds on parameters
    lower <- rep(0.0000001, n_params)
    upper <- rep(10, (n_params-1))
    upper <- c(upper, 0.999)    # append upper bound on phi
    
  } else {
    n_params <- m+1     # number of params to be estimated
    
    ## 1: initial parameter estimates
    #initial_params <- rep(0.002, n_params) #old 
    initial_params <- rep(0.002, (n_params-1))
    initial_params <- c(0.1, initial_params)
    
    ## 2: estimation
    
    # bounds on parameters
    lower <- rep(0.0000001, n_params)
    upper <- rep(10, n_params)
  }
  
  # obtain optimization results
  result = optim(par=initial_params, fn=function(params, y, X, robust){    # for optim fct
    # result = solnp(initial_params, fun = function(params, y, X) {     # for solnp fct
    
    ## 2.1 obtain total log-likelihood
    
    # run Kalman Filter
    KF_output <- KF_tvplp(params=params, y=y, X=X, robust=robust)
    
    # extract v and F
    v <- c(KF_output$v)[-1]  # equivalent to letting likelihood contr start at t=2
    F <- c(KF_output$F)[-1]
    
    # Define sample size
    n <- length(v)
    
    # compute log-likelihood
    ll_tot <- -(n/2) * log(2*pi) - (1/2) * ( sum(log(abs(F))) + sum(v * F^(-1) * v) )
    
    
    ## 2.2 obtain loss 
    loss <- -(ll_tot)/n 
    
    return(loss)
    
  }, method="L-BFGS-B", lower=lower, upper=upper, y=y, X=X, robust=robust)            # for optim fct
  #}, LB = lower, UB = upper,  y=y, X=X)     # for solnp fct
  
  # ensure convergence
  if(result$convergence != 0){
    print(result$message)
  }
  
  ML_estimate <- result$par
  
  return(ML_estimate)
}

tvplp_estimator <- function(y, x, W, dates=NULL, robust){
  # Description
  # Estimates a tvp-lp model and returns filtered and smoothed tvp-lp path estimates
  # Arguments:
  # y: dependent variable data; (Tx1) vector
  # x: regressor of interest data; (Tx1) vector
  # W: matrix of further regressors (incl. controls and constant vector)
  # dates: a character vector of dates (dates correspond to dates of y) (optional)
  # robust: If robust == TRUE, autocorrelation correction by Bryson and Henrikson (1968) is applied. 
  
  ## Step 1: arrange input data
  X <- cbind(x, W)
  
  ## Step 2: parameter estimation 
  ML_estimates <- ML_tvplp(y=y, X=X, robust=robust)
  
  ## Step 3: obtain filtered path estimate
  KF_output <- KF_tvplp(params=ML_estimates, y=y, X=X, robust=robust)
  filtered_path <- KF_output$a[1,]          # variable nr. 1 is the regressor of interest
  
  ## Step 4: obtain smoothed path estimate
  KFS_output <- KFS_tvplp(params=ML_estimates, y=y, X=X, robust=robust)
  smoothed_path <- KFS_output$alpha_hat[1,] # variable nr. 1 is the regressor of interest
  
  ## Step 5: obtain confidence intervals for filtered path
  P <- KF_output$P[1,1,]                    # variable nr 1 is the regressor of interest
  P[P <= 0] <- NA                           # replace negative variances with NA (should be adressed at later point)
  ci_f_l <- filtered_path - qnorm(0.95)*sqrt(P)    
  ci_f_u <- filtered_path + qnorm(0.95)*sqrt(P)
  
  ## Step 6: obtain confidence intervals for smoothed path
  V <- KFS_output$V[1,1,]                   # variable nr. 1 is the regressor of interest
  V[V <= 0] <- NA                           # replace negative variances with NA (should be adressed at later point)
  ci_s_l <- smoothed_path - qnorm(0.95)*sqrt(V)    
  ci_s_u <- smoothed_path + qnorm(0.95)*sqrt(V)
  
  ## Step 7: combine output in matrix and add dates 
  
  paths_output <- cbind(filtered_path, smoothed_path, ci_f_l, ci_f_u, ci_s_l, ci_s_u)
  
  if (!is.null(TRUE)) rownames(paths_output) <- dates
  
  ## Step 7: return output 
  lOut <- list("ML_estimates"=ML_estimates, "paths_output"=paths_output)
  
  return(lOut)
}
