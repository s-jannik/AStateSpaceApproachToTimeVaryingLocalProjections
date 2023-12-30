# Author: Jannik Steenbergen

# remove all plots currently loaded:
while (!is.null(dev.list()))  dev.off()

# clear command window
cat("\014") 

# remove all variables from history:
rm(list = ls())

# set working directory
setwd("~/AU master/Topics")

# libraries
library(readxl)
library(plotly)
library(sandwich)
source("kf_kfs.R")
source("support_functions.R")
source("studies.R")
source("tvplp.R")

fred_qd <- read.csv("fred_qd.csv")

############################ data preparation ##################################

## (original) data to be used
o_data <- fred_qd
o_data <- preprocess_data(data=o_data, "fred_qd")

## plot untransformated series
par(mfrow=c(2,2))
for (i in 1:ncol(o_data)) ts.plot(o_data[,i], ylab="", main=paste(colnames(o_data)[i], "(untransformed)", sep=" "))


## stationarity transformations
tcodes <- c("CPI"=9, "IR"=2, "UR"=2)      # variable transformation codes
data <- transform_data_tcode(data=o_data, tcode=tcodes, use_which_data="all_available", as_ts_object=FALSE)
data <- na.omit(data)   # removes first 4 rows containing NAs

## plot transformated series
par(mfrow=c(2,2))
for (i in 1:ncol(data)) ts.plot(data[,i], ylab="", main=paste(colnames(data)[i], "(transformed)", sep=" "))

## nice plot for paper 
par(mfrow=c(3,2))
for (i in c(1,3)) {
  ts1 <- ts(o_data[,i], start=c(1959, 1), frequency=4)
  ts2 <- ts(data[,i], start=c(1960, 1), frequency=4)
  ts.plot(ts1, ylab="", main=paste(colnames(o_data)[i], "(untransformed)", sep=" "), xlab="")
  ts.plot(ts2, ylab="", main=paste(colnames(data)[i], "(transformed)", sep=" "), xlab="")
}

ts1 <- ts(o_data[,2], start=c(1959, 1), frequency=4)
ts.plot(ts1, ylab="", main=paste(colnames(o_data)[2], "(untransformed)", sep=" "), xlab="")


########################## Example simulation study ############################

## Example simulated dataset
# simulate dataset
T <- 500    # length of simulate series
sim_data <- simulate_dgp(T=T, phi=c(0.6, -0.2, 0.1))
y <- sim_data[,"y"]
x <- sim_data[,"x"]
W <- sim_data[,c("ones", "w_lag1")]

## plot simulated dataset
par(mfrow=c(3,2))
for(i in 1:ncol(sim_data)) ts.plot(sim_data[,i], ylab="", main=colnames(sim_data)[i])

# Example implementation LP Jorda (2005)
A <- as.data.frame(cbind(y,x,W))

lp_model <- lm(y ~ -1 + x + ones + w_lag1, data=A)
lp_model <- lm(y ~ ., data=A)
beta_lp <- as.numeric(lp_model$coefficients[1])
cov_matrix_nw <- NeweyWest(lp_model)
se_beta <- sqrt(cov_matrix_nw[1,1])

ci_l_lp <- beta_lp - qnorm(0.95)*se_beta
ci_u_lp <- beta_lp + qnorm(0.95)*se_beta

par(mfrow=c(1,1))
plot(NULL, ylab="y label", xlab="x label", main="title",
     xlim=c(0, 10), ylim=c(-0.5, 0.5))
abline(h=c(beta_lp, ci_l_lp, ci_u_lp), lty=c(1, 2, 2), lwd=c(3,3,3))




############################### Simulation study ###############################

par(mfrow=c(3,1))

## Simulation parameters
T <- 550                # length of simulate series
burn_in <- 50           # length of burn-in period
lost_obs <- 1           # number of obs lost from using robust TVP-LP
T_c <- 550-burn_in-lost_obs # corrected T

M <- 1                # number of simulations 
phi = c(0.6, -0.2, 0.1) # autocorrelation structure
n_params = 5

## initialize array to hold all paths
paths <- array(NA, dim=c((T-burn_in-lost_obs), 3, 4), dimnames=list(c(), c("beta_hat", "ci_l", "ci_u"), 
                                                                    c("filtered_tvplp", "filtered_tvplp_robust", "smoothed_tvplp", "smoothed_tvplp_robust")))
## initialize result matrices
colnames <- c("tvplp1", "tvplp2", "lp")
coverage_rates <- matrix(NA, nrow=M, ncol=3, dimnames=list(c(), colnames))
ci_performance <- array(NA, dim=c(M, 3, 4), dimnames=list(c(), colnames, c("coverage_filtered", "coverage_smoothed", "width_filtered", "width_smoothed")))
RMSE <- array(NA, dim=c(M, 3, 2), dimnames=list(c(), colnames, c("filtered", "smoothed")))
MAE <- array(NA, dim=c(M, 3, 2), dimnames=list(c(), colnames, c("filtered", "smoothed")))
param_estimates <- array(NA, di=c(n_params, 3, M), dimnames = list(c(), colnames, c()))


## Run simulation study 
for (m in 1:M){
  print(m)
  
  # simulate dataset
  sim_data <- simulate_dgp(T=T, phi=phi)
  y <- sim_data[,"y"]
  x <- sim_data[,"x"]
  W <- sim_data[,c("ones", "w_lag1")]
  
  ## plot simulated dataset
  #par(mfrow=c(3,2))
  #for(i in 1:ncol(sim_data)) ts.plot(sim_data[,i], ylab="", main=colnames(sim_data)[i])
  
  ## Obtain true beta
  beta_true <- sim_data[,"beta"]
  beta_true <- beta_true[-T]                  # correction for lost_obs
  beta_true <- beta_true[-c(1:burn_in)]       # discarding burn-in period
  
  ## Estimation for robust=FALSE
  result1 <- tvplp_estimator(y=y, x=x, W=W, dates=NULL, robust=FALSE)
  paths_output1 <- result1$paths_output
  
  # Correct for lost_obs
  paths_output1 <- paths_output1[-T,] 
  
  # Discard burn-in period
  paths_output1 <- paths_output1[-c(1:burn_in),]
  
  # Define paths 
  beta_filtered1 <- paths_output1[,"filtered_path"]
  beta_smoothed1 <- paths_output1[,"smoothed_path"]
  ci_f_l1 <- paths_output1[,"ci_f_l"] # filtered intervals
  ci_f_u1 <- paths_output1[,"ci_f_u"]
  ci_s_l1 <- paths_output1[,"ci_s_l"] # smoother intervals
  ci_s_u1 <- paths_output1[,"ci_s_u"]
  
  ## Estimation for robust=TRUE
  result2 <- tvplp_estimator(y=y, x=x, W=W, dates=NULL, robust=TRUE)
  paths_output2 <- result2$paths_output
  
  # Discard burn-in period
  paths_output2 <- paths_output2[-c(1:burn_in),]
  
  # Define paths 
  beta_filtered2 <- paths_output2[,"filtered_path"]
  beta_smoothed2 <- paths_output2[,"smoothed_path"]
  ci_f_l2 <- paths_output2[,"ci_f_l"] # filtered intervals
  ci_f_u2 <- paths_output2[,"ci_f_u"]
  ci_s_l2 <- paths_output2[,"ci_s_l"] # smoother intervals
  ci_s_u2 <- paths_output2[,"ci_s_u"]
  
  ## Estimation for benchmark LP
  data_lp <- as.data.frame(cbind(y,x,W))
  
  lp_model <- lm(y ~ -1 + x + ones + w_lag1, data=data_lp)
  beta_lp <- as.numeric(lp_model$coefficients[1])
  cov_matrix_nw <- NeweyWest(lp_model)
  se_beta <- sqrt(cov_matrix_nw[1,1])
  
  ci_l_lp <- beta_lp - qnorm(0.95)*se_beta
  ci_u_lp <- beta_lp + qnorm(0.95)*se_beta
  
  
  ## Plot
  ylim_lower <- min(ci_s_l1[!is.na(ci_s_l1)], beta_true)
  ylim_upper <- max(ci_s_u1[!is.na(ci_s_u1)], beta_true)
  ts.plot(beta_true, type="l", col="black", ylim=c(ylim_lower,ylim_upper), lwd=1.5, ylab="")
  #lines(beta_smoothed1, col="blue")
  #lines(beta_smoothed2, col="red")
  
  lines(beta_filtered1, col="blue", lty=1, lwd=1.5)
  lines(beta_filtered2, col="red", lty=1, lwd=1.5)
  
  #lines(ci_s_l1, col="blue", lty=3, lwd=1.5)
  #lines(ci_s_u1, col="blue", lty=3, lwd=1.5)
  lines(ci_f_l1, col="grey", lty=3, lwd=1.5)
  lines(ci_f_u1, col="grey", lty=3, lwd=1.5)
  
  #lines(ci_s_l2, col="red", lty=3, lwd=1.5)
  #lines(ci_s_u2, col="red", lty=3, lwd=1.5)
  #lines(ci_f_l2, col="grey", lty=3, lwd=1.5)
  #lines(ci_f_u2, col="grey", lty=3, lwd=1.5)
  
  abline(h=beta_lp, col="green", lty=1, lwd=1.5)
  abline(h=ci_l_lp, col="green", lty=2, lwd=1.5)
  abline(h=ci_u_lp, col="green", lty=2, lwd=1.5)
         
  # store RMSE 
  RMSE[m, "tvplp1", "filtered"] <- sqrt(mean((beta_true - beta_filtered1)^2))
  RMSE[m, "tvplp1", "smoothed"] <- sqrt(mean((beta_true - beta_smoothed1)^2))
  RMSE[m, "tvplp2", "filtered"] <- sqrt(mean((beta_true - beta_filtered2)^2))
  RMSE[m, "tvplp2", "smoothed"] <- sqrt(mean((beta_true - beta_smoothed2)^2))
  RMSE[m, "lp", "filtered"] <- sqrt(mean((beta_true - beta_lp)^2))
  RMSE[m, "lp", "smoothed"] <- sqrt(mean((beta_true - beta_lp)^2))

  # store MAE
  MAE[m, "tvplp1", "filtered"] <- mean(abs(beta_true - beta_filtered1))
  MAE[m, "tvplp1", "smoothed"] <- mean(abs(beta_true - beta_smoothed1))
  MAE[m, "tvplp2", "filtered"] <- mean(abs(beta_true - beta_filtered2))
  MAE[m, "tvplp2", "smoothed"] <- mean(abs(beta_true - beta_smoothed2))
  MAE[m, "lp", "filtered"] <- mean(abs(beta_true - beta_lp))
  MAE[m, "lp", "smoothed"] <- mean(abs(beta_true - beta_lp))
  
  # store coverage rate
  not_na <- !is.na(ci_f_u1) # indicator of not being na
  ci_performance[m, "tvplp1", "coverage_filtered"] <- sum(beta_true[not_na] < ci_f_u1[not_na] & beta_true[not_na] > ci_f_l1[not_na])/T_c
  
  not_na <- !is.na(ci_s_u1)
  ci_performance[m, "tvplp1", "coverage_smoothed"] <- sum(beta_true[not_na ] < ci_s_u1[not_na ] & beta_true[not_na ] > ci_s_l1[not_na ])/T_c
  
  not_na <- !is.na(ci_f_u2)
  ci_performance[m, "tvplp2", "coverage_filtered"] <- sum(beta_true[not_na] < ci_f_u2[not_na] & beta_true[not_na] > ci_f_l2[not_na])/T_c
  
  not_na <- !is.na(ci_s_u2)
  ci_performance[m, "tvplp2", "coverage_smoothed"] <- sum(beta_true[not_na] < ci_s_u2[not_na] & beta_true[not_na] > ci_s_l2[not_na])/T_c
  
  not_na <- !is.na(ci_u_lp)
  ci_performance[m, "lp", "coverage_filtered"] <- sum(beta_true[not_na] < ci_u_lp[not_na] & beta_true[not_na] > ci_l_lp[not_na])/T_c
  
  # store interval width
  ci_performance[m, "tvplp1", "width_filtered"] <- mean((ci_f_u1-ci_f_l1)[!is.na(ci_f_u1)])
  
  ci_performance[m, "tvplp1", "width_smoothed"] <- mean((ci_s_u1-ci_s_l1)[!is.na(ci_s_u1)])
  
  ci_performance[m, "tvplp2", "width_filtered"] <- mean((ci_f_u2-ci_f_l2)[!is.na(ci_f_u2)])
  
  ci_performance[m, "tvplp2", "width_smoothed"] <- mean((ci_s_u2-ci_s_l2)[!is.na(ci_s_u2)])
  
  ci_performance[m, "lp", "width_filtered"] <- mean((ci_u_lp-ci_l_lp)[!is.na(ci_u_lp)])
  
  # store parameter estimates
  param_estimates[1:(n_params-1),"tvplp1",m] <- result1$ML_estimates
  param_estimates[1:n_params,"tvplp2",m] <- result2$ML_estimates
  #param_estimates[1:(n_params-1),"lp",m] <- lp_estimates
  
  
  #print(result1$ML_estimates)
  #print(result2$ML_estimates)
}

## aRMSE ratio
aRMSE <- matrix(NA, nrow=2, ncol=3, dimnames=list(c("filtered", "smoothed"), c("tvplp1", "tvplp2", "lp")))
aRMSE_ratio <- matrix(NA, nrow=2, ncol=3, dimnames=list(c("filtered", "smoothed"), c("tvplp1", "tvplp2", "lp")))
aMAE_ratio <- matrix(NA, nrow=2, ncol=3, dimnames=list(c("filtered", "smoothed"), c("tvplp1", "tvplp2", "lp")))

for (i in 1:dim(RMSE)[2]){
  for (j in 1:dim(RMSE)[3]){
    
    aRMSE[j, i] <- mean(RMSE[,i,j])
    aRMSE_ratio[j, i] <- mean(RMSE[,i,j] / RMSE[,"lp",j])
    aMAE_ratio[j, i] <- mean(MAE[,i,j] / MAE[,"lp",j])
    
  }
}

## CI performance
a_ci_performance <- matrix(NA, nrow=4, ncol=3, dimnames=list(c("coverage_filtered", "coverage_smoothed",
                                                               "width_filtered", "width_smoothed"), c("tvplp1", "tvplp2", "lp")))

for (i in 1:dim(ci_performance)[2]){
  for (j in 1:dim(ci_performance)[3]){
    
    a_ci_performance[j, i] <- mean(ci_performance[,i,j])
    
  }
}

## MSE params
MSE_params <- matrix(NA, nrow=5, ncol=3)
true_params <- c(0.16, 0.001, 0.001, 0.001, NA)

for (i in 1:dim(param_estimates)[1]){
  for (j in 1:dim(param_estimates)[2]){
    
    MSE_params[i, j] <- mean((true_params[i] - param_estimates[i,j, ])^2)
    
  }
}

# print all results
aRMSE
aRMSE_ratio
aMAE_ratio
a_ci_performance
MSE_params





############################## Empirical study #################################

# Settings
H_min <- 0
H_max <- 8
smoothed = FALSE
robust = FALSE
maxlags <- c("CPI"=4, "IR"=4, "UR"=4)

# define variables for h=0 model
data_reg <- as.matrix(get_dataset(data=data, maxlags=maxlags))
y <- data_reg[,"CPI"]
x <- data_reg[,"IR"]
W <- data_reg[,!(colnames(data_reg) %in% c("CPI", "IR"))]
dates <- rownames(data_reg)
T <- nrow(data_reg)

# run empirical study
results <- empirical_study(y=y, x=x, W=W, H_min=H_min, H_max=H_max, dates=dates, smoothed=smoothed, robust=robust)
paths_output <- results$paths_output
ML_estimates <- results$ML_estimates
surface_matrix <- results$surface_matrix
benchmark_results <- results$benchmark_results


################### individual paths 

# plot individual paths 
par(mfrow=c(3,3))
for (i in 1:ncol(surface_matrix)) ts.plot(surface_matrix[,i], ylab="", main=paste("h =", (i-1), sep=" "))

# plot individual paths with benchmark
par(mfrow=c(3,3))
for (i in 1:ncol(surface_matrix)) {
  
  # tvplp estimates
  beta_tvplp <- surface_matrix[,i]
  if (smoothed == TRUE){
    ci_l_tvplp <- paths_output[[i]][,"ci_s_l"]
    ci_u_tvplp <- paths_output[[i]][,"ci_s_u"]
  } else {
    ci_l_tvplp <- paths_output[[i]][,"ci_f_l"]
    ci_u_tvplp <- paths_output[[i]][,"ci_f_u"]
  }
  
  # lp estimates
  beta_lp <- benchmark_results[i, "beta_lp"]
  ci_l_lp <- benchmark_results[i, "ci_l_lp"]
  ci_u_lp <- benchmark_results[i, "ci_u_lp"]
  
  # plot
  highest_obs <- max(ci_u_tvplp[-(1:25)], ci_u_lp)
  lowest_obs <- min(ci_l_tvplp[-(1:25)], ci_l_lp)
  ylim=c(lowest_obs, highest_obs)
  
  length(seq(from=1962, to=2021, by=1)) 
  dates <- c(1962, 1962, 1962, floor(seq(from=1963, to=2022.75, by=0.25)), 2023, 2023, 2023)
  
  # make sure surface_matrix and CI's have same dimensions (quick fix)
  if (nrow(surface_matrix) != length(ci_l_tvplp)){
    ci_l_tvplp <- ci_l_tvplp[-c(1:(length(ci_l_tvplp)-nrow(surface_matrix)))]
    ci_u_tvplp <- ci_u_tvplp[-c(1:(length(ci_u_tvplp)-nrow(surface_matrix)))]
  } 
  
  # define ts matrix object
  ts.series <- ts(cbind(surface_matrix[,i], ci_l_tvplp, ci_u_tvplp), start=c(1962, 2), frequency=4)
  
  ts.plot(ts.series, ylim=ylim, ylab="", main=paste("h =", (i-1), sep=" "), lwd=c(1.2, 1.2, 1.2), lty=c(1,2,2), xlab="")
  abline(h=beta_lp, col="red", lwd=1.2)
  abline(h=ci_l_lp, col="red", lty=2, lwd=1.2)
  abline(h=ci_u_lp, col="red", lty=2, lwd=1.2)
  
}

################### accumulated paths 

# plot accumulated path
par(mfrow=c(1,1))
ts.plot(apply(surface_matrix, 1, sum))

# plot accumulated path with benchmark
par(mfrow=c(1,1))
ts.plot(apply(surface_matrix, 1, sum))
abline(h=sum(benchmark_results[,1]), col="red")


################# IRF function 
horizons <- 0:8
plot(horizons, benchmark_results[,1], type="l", xlab="Horizon", ylab="", lwd=1.5)

write.csv(surface_matrix, "surface_matrix_cpi4_ir0_ur4.csv", row.names=FALSE)

surface_matrix[72, , drop="FALSE"]
surface_matrix[69, , drop="FALSE"]
