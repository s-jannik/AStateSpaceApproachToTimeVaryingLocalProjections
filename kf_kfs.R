# Author: Jannik Steenbergen

############################### Kalman Filter ##################################


KF <- function(y, Z, T, R, d, c, H, Q, initialization, aim){
  # Description: 
    # This functions performs Kalman filtering on a state space model with a scalar dependent variable and state vector (i.e. alpha is 1xm, where m is the number of state variables)
  # Arguments:
    # y: nx1  vector
    # Z: nxm  matrix (where m is number of state variables)
    # T: mxm  matrix 
    # R: mxq  matrix (where q is the number of state error processes)
    # d: nx1  vector (vector of time-varying obs. eq. means. Currently has to be a vector (can also be a vector of constants))
    # c: mxn  matrix (matrix of time-varying state eq. means. Currently has to be a matrix  (can also be a matrix of constants))
    # H: scalar
    # Q: qxq  matrix 
    # initialization: an m-element list. The m-th element that holds a character "diffuse" if the m-th state variable should be initialized with diffuse initialization. The m-th element holds a vector in form (a1, P1) if the m-th state variable should be initialized at a specific set of starting values.
  
  # define N
  N <- length(y)
  
  # define m
  if (is.matrix(T)){
    m <- dim(T)[1] 
  } 
  else{
    m <- 1
  }
   
  # initializations
  v <- c()    # will be (nx1)
  F <- c()    # will be (nx1)
  K <- matrix(NA, nrow=m, ncol=N) # will be (mxn)
  
  # initialization of a and P
  a <- matrix(NA, nrow=m, ncol=N+1)   # +1 because forecasting
  P <- array(0, dim = c(m, m, N+1))   # note: 0 instead of NA to ensure no dependence between states (i.e. all off-diagonals are 0)
  
  a_f <- matrix(NA, nrow=m, ncol=N)     # filtering (nowcasting)
  P_f <- array(NA, dim = c(m, m, N))
  
  # initial conditions for a and P
  for (z in 1:m){
    
    if (initialization[z] == "diffuse"){        # note: possibly consider double brackets [[z]] (appears to not make a difference currently)
      
      a[z, 1] <- 0          # initialized at 0
      P[z, z, 1] <- 10^7    # initialized at 10^7
      
    } 
    else {
      
      a[z, 1] <- initialization[[z]][1]
      P[z, z, 1] <- initialization[[z]][2]
      
    }
  }
  
  # recursions
  
  for (t in 1:N){
    
    # initial definitions (for brevity and ensuring all elements have right dimensions)
    
    if (is.matrix(Z)){
      Zt <- Z[t, , drop=FALSE]
    } else {
      Zt <- Z
    }
    
    if(is.matrix(c)){
      ct <- c[,t]
    } else {
      ct <- c
    }
    
    if(length(d)>1){      # i.e. if d is time-varying
      dt <- d[t]
    } else {
      dt <- d
    }
    
    at <- a[ , t, drop=FALSE]
    Pt <- P[, , t]
    
    # the recursions
    v[t] <- y[t] - Zt %*% at - dt
    F[t] <- Zt %*% Pt %*% t(Zt) + H
    K[, t] <- T %*% Pt %*% t(Zt) * F[t]^(-1)    # important that K[, t] has right dimensions
    Kt <- K[, t]
    
    # filtering step
    a_f[,t] <- at + Pt %*% t(Zt) * F[t] ^ (-1) * v[t]
    P_f[, , t] <- Pt - Pt %*% t(Zt) %*% Zt %*% Pt * F[t]^(-1)
    
    # prediction step
    a[, (t + 1)] <- T %*% at + Kt * v[t] + ct
    P[, , (t + 1)] <- T %*% Pt %*% t(T) + R %*% Q %*% t(R)  - Kt %*% t(Kt) * F[t]
    
  }
  
  # Get rid of last element (N+1) of a and P
  a<-a[,-(N+1), drop=FALSE]
  P<-P[, , -(N+1), drop=FALSE]
  
  # Gather all output in list and return list
  
  if (aim == "f") {
    # output for filtering
    output<-list()
    output[[1]]<-a_f
    output[[2]]<-P_f
    output[[3]]<-v
    output[[4]]<-F
    output[[5]]<-K
    names(output)<-c("a_f", "P_f", "v", "F", "K")
    return(output)
    
    
  } else if (aim == "p") {
    # output for prediction 
    output<-list()
    output[[1]]<-a
    output[[2]]<-P
    output[[3]]<-v
    output[[4]]<-F
    output[[5]]<-K
    names(output)<-c("a", "P", "v", "F", "K")
    return(output)
  }
  
}


############################ Kalman Filter Smoother ############################

KFS <- function(y, Z, T, R, d, c, H, Q, initialization){
  # Description:
    # Function running the Kalman Filter Smoother
  # Arguments:
    # see function "KF"
  
  # run Kalman filter
  KF_output <- KF(y=y, Z=Z, T=T, R=R, d=d, c=c, H=H, Q=Q, initialization=initialization, aim="p")   # note: for KFS we are only interested in a_t, not a_(t|t)
  
  a <- KF_output$a
  P <- KF_output$P
  v <- KF_output$v
  F <- KF_output$F
  K <- KF_output$K
  
  # define N
  n <- length(y)
  
  # define m
  if (is.matrix(T)){
    m <- dim(T)[1] 
  } 
  else{
    m <- 1
  }
  
  # define empty matrices/arrays to store KFS output
  L <- array(NA, dim=c(m, m, (n+1)))        # ("+1" because we also need r_0 and N_0)
  r <- matrix(NA, nrow=m, ncol=(n+1))
  N <- array(NA, dim=c(m, m, (n+1)))
  alpha_hat <- matrix(NA, nrow=m, ncol=(n+1))
  V <- array(NA, dim=c(m, m, (n+1)))
  
  # initializations 
  r[, (n+1)] <- 0
  N[, , (n+1)] <- 0 
  
  # Kalman filter smoother recursions:
  for (t in n:1){
    
    t_a=t+1     # for practical purposes relating to store of results
    
    # initial definitions (for brevity)
    Kt <- K[, t, drop=FALSE]
    
    if (is.matrix(Z)){
      Zt <- Z[t, , drop=FALSE]
    } else {
      Zt <- Z
    }
    
    at <- a[ , t, drop=FALSE]
    Pt <- P[, , t]
    
    
    # the recursions
    L[, , t_a] <- T - Kt %*% Zt
    Lt <- L[, , t_a]      # the Lt vs t_a can be confusing
    
    rt <- r[ , t_a]
    r[ , (t_a-1)] <- t(Zt) %*% F[t]^(-1) * v[t] + t(Lt) %*% rt
    rt_l1 <- r[ , (t_a-1)]
    
    Nt <- N[, , t_a]
    N[, , (t_a-1)] <- t(Zt) %*% F[t]^(-1) %*% Zt + t(Lt) %*% Nt %*% Lt
    Nt_l1 <- N[, , (t_a-1)]
    
    alpha_hat[, t_a] <- at + Pt %*% rt_l1
    V[, ,(t_a)] <- Pt - Pt %*% Nt_l1 %*% Pt
    
    
  }
  
  # discard first "element" of each object (containing NA for all other than r and N)
  L <- L[, , -1, drop=FALSE]
  r <- r[, -1, drop=FALSE]
  N <- N[, , -1, drop=FALSE]
  alpha_hat <- alpha_hat[, -1, drop=FALSE]
  V <- V[, , -1, drop=FALSE]
  
  # combine output
  KFS_output <- list(L, r, N, alpha_hat, V)
  names(KFS_output) <- c("L", "r", "N", "alpha_hat", "V")
  
  return(KFS_output)
  
}

################################################################################
########Implementations for replicating figures from Durbin-Koopman book########
################################################################################

#################### Quasi Maximum Likelihood Estimator ########################

initial_estimate_func_unconstrained<-function(y, model){
  
  if (model == "LLM"){
    
    sigma2_eps <- 15000     # semi-arbitrarily chosen based on DK-book estimates (DK-book estimates: sigma2_eps = 15099, sigma2_eta = 1469.1)
    sigma2_eta <- 1450
    
    initial_params <- c(sigma2_eps, sigma2_eta)
    
  }
  
  if (model == "SVM"){
      
    phi <- 0.9950 # Given in DK book Chapter 14.5 page 320
    
    omega <- (1-phi)*(mean(y)+1.27)
    sigma2_eta <- (1-phi^2)*(var(y)-((pi^2)/2))
    
    # Define bounds
    lb_phi=0; ub_phi=1; lb_sigma2_eta=0
    
    # transform phi and sigma2_eta to unconstrained form 
    z=(phi-lb_phi)/(ub_phi-lb_phi)
    phi_unconstrained=log(z/(1-z))
    
    sigma2_eta_unconstrained=log(sigma2_eta - lb_sigma2_eta)
    
    initial_params <- c(omega, phi_unconstrained, sigma2_eta_unconstrained)
    
  }
  return(initial_params)
}

get_constrained_params<-function(unconstrained_params, model){
  
  if (model == "LLM"){
    
    constrained_params <- unconstrained_params
  }
  
  if (model == "SVM"){
  
    # Define bounds
    lb_phi=0; ub_phi=1; lb_sigma2_eta=0     # currently being defined also in another fct.: suboptimal
    
    # transform to constrained form
    constrained_params <- unconstrained_params
    constrained_params[2] <- lb_phi+(ub_phi - lb_phi) / (1.0 + exp(-unconstrained_params[2]))
    constrained_params[3] <- exp(unconstrained_params[3])+lb_sigma2_eta
  }
  
  return(constrained_params)
}

ll_func <- function(params, y, X = NULL, model){
  # Description:
    # Function returning the total log likelihood
  # Arguments:
    # params: a vector of parameters to be estimated
      # y: (nx1) vector of observations
      # X: (nxd) matrix of exogeneous regressors (optional)
      # model: a character taking one of the following options:
        # "LLM"     (local-level model)
        # "SVM"     (stochastic volatility model)
        # "TVP-LP"
  
  
  # Define model-specific state space form parameters
  
  if (model == "LLM"){
    
    # define parameters to be estimated
    sigma2_eps <- params[1]
    sigma2_eta <- params[2]
    
    # define state space form parameters
    Z=1; T=1; R=1; d=0; c=0; H=sigma2_eps; Q=sigma2_eta
    
    # define initialization and aim 
    initialization = "diffuse"
    aim = "p"
      
  }
  
  if (model == "SVM"){
    
    # define parameters to be estimated 
    omega <- params[1]
    phi <- params[2]
    sigma2_eta <- params[3]
    
    Z=1; T=phi; R=1; d=-1.27; c=omega; H=(pi^2)/2; Q=sigma2_eta
    
    # define initialization and aim
    a1 <- omega/(1-phi)
    P1 <- sigma2_eta/(1-(phi^2))
    initialization <- list(c(a1, P1))
    
    aim <- "p"
    
  }
  
  # Run Kalman Filter
  KF_output <- KF(y=y, Z=Z, T=T, R=R, d=d, c=c, H=H, Q=Q, initialization=initialization, aim=aim)
  
  # Extract v and F
  v <- c(KF_output$v)
  F <- c(KF_output$F)
  
  # Define sample size
  n <- length(v)
  
  # compute log-likelihood
  ll_tot <- -(n/2) * log(2*pi) - (1/2) * ( sum(log(abs(F))) + sum(v * F^(-1) * v) )
    # t=1 is currently included in the likelihood. I should consider whether to remove it. If i do, remember to change the sample size used to compute avg. likelihood in loss_func
  
  return(ll_tot)
}

loss_func <- function(params, y, X = NULL, model){
  
  unconstrained_params <- params
  
  constrained_params <- get_constrained_params(unconstrained_params, model=model)
  
  ll_tot <- ll_func(params=constrained_params, y=y, X=X, model=model)
  
  n <- length(y)
  
  loss = -(ll_tot)/n  # note: division should be corrected if likelihood in ll_func starts at t=2
  
  return(loss)
}

# Maximum likelihood estimation function
MLE_func<-function(y, X = NULL, model){
  
  # estimates
  initial_params <- initial_estimate_func_unconstrained(y, model=model)
  
  # result = optim(par=initial_params, fn=loss_func, y=y, method=c("L-BFGS-B"), lower = lb, upper = ub)
  result <- optim(par=initial_params, fn=loss_func, y=y, X=X, model=model, method=c("Nelder-Mead"), control=list(maxit=100000))
  
  # ensure convergence 
  if(result$convergence != 0){
    print(result$message)
  }
  
  unconstrained_estimate=result$par
  ML_estimate=get_constrained_params(unconstrained_estimate, model=model)
  
  return(ML_estimate)
}

