# Author: Jannik Steenbergen

####Functions used for implementation of simulation study and empirical study###

# simulation function 
simulate_dgp <- function(T, phi){
  # Description:
    # Simulates from a model of the form 
      # y_t = mu_t + beta_t*x_t + gamma_t*w_{t-1} + v_t
  # Arguments: 
    # phi: if obs eq. errors v_t should be uncorrelated phi=0; 
    # else phi=c(phi1) or phi=c(phi1, phi2, phi3)
  
  ## parameter choices
  sigma2_v <- 0.16
  sigma2_u_mu <- 0.001
  sigma2_u_beta <- 0.001
  sigma2_u_gamma <- 0.001
  
  ## define constant (needed to later estimate mu_t)
  ones <- rep(1, T)
  
  ## draw x_t (x) and w_{t-1} (w_lag1)
  x <- rnorm(n=T, mean=0, sd=1)
  w_lag1 <- rnorm(n=T, mean=0, sd=1)
  
  ## draw errors (v)
  
  if (length(phi) == 1) if(phi == 0) v <- rnorm(n=T, mean=0, sd=sqrt(sigma2_v))
  
  if (length(phi) == 1){
      
    iid_errors <- rnorm(n=(T+1), mean=0, sd=sqrt(sigma2_v)) # +1 since 1 pre-obs needed
      
    v <- c()
    v[1] <- iid_errors[1]   # pre-observation
      
    for (t in 2:(T+1)) v[t] <- phi*v[t-1] + iid_errors[t]
    
    v <- v[-1] # discard pre-observation
  }
  
  if (length(phi) == 3){
     
    iid_errors <- rnorm(n=(T+3), mean=0, sd=sqrt(sigma2_v)) # +3 since 3 pre-obs needed
    
    v <- c()
    v[1:3] <- iid_errors[1:3]   # pre-observations
    
    for (t in 4:(T+3)) v[t] <- phi[1]*v[t-1] + phi[2]*v[t-2] + phi[3]*v[t-3] + iid_errors[t]
    
    v <- v[-c(1:3)] # discard pre-observations
  }
  
  ## simulate parameter paths
  
  # draw errors
  u_mu <- rnorm(n=T, mean=0, sd=sqrt(sigma2_u_mu))
  u_beta <- rnorm(n=T, mean=0, sd=sqrt(sigma2_u_beta))
  u_gamma <- rnorm(n=T, mean=0, sd=sqrt(sigma2_u_gamma))
  
  # initialize paths at 0
  mu <- c(0)
  beta <- c(0)
  gamma <- c(0)
  
  for (t in 2:T){
    
    mu[t] <- mu[t-1] + u_mu[t]
    beta[t] <- beta[t-1] + u_beta[t]
    gamma[t] <- gamma[t-1] + u_gamma[t]
    
  }
  
  ## obtain y
  y <- c()
  
  for (t in 1:T){
    
    y[t] <- mu[t] + beta[t]*x[t] + gamma[t]*w_lag1[t] + v[t]
    
  }
  
  output <- cbind(y, ones, x, w_lag1, v, mu, beta, gamma)
  
  return(output)
}



# empirical study function 
empirical_study <- function(y, x, W, H_min, H_max, dates, smoothed, robust){
  # Description:
    # empirical study of tvplp estimator
  # Arguments:
    # y: dependent variable data; (Tx1) vector
    # x: regressor of interest data; (Tx1) vector
    # W: other regressors (incl. controls); matrix 
    # dates: a character vector of dates (dates correspond to dates of y)
   # smoothed: if TRUE, returns smoothed path estimates in surface matrix, else filtered paths 
   # robust: If robust == TRUE, autocorrelation correction by Bryson and Henrikson (1968) is applied. 
  
  T <- length(y)  # sample size
  n_h <- H_max-H_min+1  # number of horizons 
  
  
  paths_output <- list()
  ML_estimates <- list()
  surface_matrix <- matrix(NA, nrow=(T-H_max), ncol=n_h)
  if (robust == TRUE) surface_matrix <- matrix(NA, nrow=(T-H_max-1), ncol=n_h)  # correction: robust loses 1 state est
  
  benchmark_results <- matrix(NA, nrow=n_h, ncol=3, dimnames = list(c(), c("beta_lp", "ci_l_lp", "ci_u_lp"))) # 3 col for estimate and conf intervals
  
  for (i in 1:n_h){
    
    ## Step 1: obtain results 
    
    if (H_min==0) h <- i-1
    if (H_min==1) h <- i
    if (H_min!=0 & H_min!=1) print("ERROR: H_min VALUE NOT ALLOWED")
    
    print(paste("h =", h, sep=" "))
    
    # produce data for h-step-ahead direct regression
    y_h <- y[(1+h):T]
    x_h <- x[1:(T-h)]
    W_h <- W[1:(T-h),]
    
    dates_h <- dates[(1+h):T]
    if (robust == TRUE) dates_h <- dates[(1+h):(T-1)] # correction since robust loses last state est
    
    # run TVP-LP estimation
    result <- tvplp_estimator(y=y_h, x=x_h, W=W_h, dates=dates_h, robust=robust)
    paths_output[[i]] <- result[["paths_output"]]
    ML_estimates[[i]] <- result[["ML_estimates"]]
    
    # run LP estimation
    data_lp <- as.data.frame(cbind(y_h,x_h,W_h[, !(colnames(W_h) %in% "ones")]))  # excluding constant
    lp_model <- lm(y_h ~ ., data=data_lp)
    beta_lp <- as.numeric(lp_model$coefficients[2]) # note: element 2 since constant is nr 1
    cov_matrix_nw <- NeweyWest(lp_model)
    se_beta <- sqrt(cov_matrix_nw[2,2]) # note: element 2 since constant is nr 1
    ci_l_lp <- beta_lp - qnorm(0.95)*se_beta
    ci_u_lp <- beta_lp + qnorm(0.95)*se_beta
    
    benchmark_results[i, ] <- c(beta_lp, ci_l_lp, ci_u_lp)
    
    ## Step 2: obtain surface plot results
    if (smoothed == TRUE) {
      beta_h <- paths_output[[i]][,"smoothed_path"]
    } else {
      beta_h <- paths_output[[i]][,"filtered_path"]
    }
    beta_h <- beta_h[(H_max+1-h):length(beta_h)]  # ensures all paths have same length
    surface_matrix[, i] <- beta_h
    
  }
  
  # add dates to surface_matrix 
  if (robust != TRUE) rownames(surface_matrix) <- dates[(H_max+1):T]
  if (robust == TRUE) rownames(surface_matrix) <- dates[(H_max+1):(T-1)]  # quick fix 
  
  return(list("paths_output" = paths_output, "ML_estimates" = ML_estimates, "surface_matrix" = surface_matrix, 
              "benchmark_results" = benchmark_results))
}














































