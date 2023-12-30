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
source("kf_kfs.R")

################################################################################

# nile data for LLM (DK-book)
nile_data <- read.table("Nile.csv", header = TRUE, sep = ",")
y_nile <- nile_data$Nile

# returns data for SVM (DK-book)
returns <- readxl::read_excel("sv.xlsx")
y_returns <- returns$GBPUSD/100

# perform the SV model data transformation (see e.g. DK-book)
x <- log((y_returns-mean(y_returns))^2)

################# Testing the KF function on the Nile data #####################

# Kalman Filter output for fig 2.1
KF_output=KF(y=y_nile, Z=1, T=1, R=1, d=0, c=0, H=15099, Q=1469.1, initialization=list("diffuse"), aim="p")

# Extract all KF output separately and combine in matrix 
a <- c(KF_output$a)
P <- c(KF_output$P)
v <- c(KF_output$v)
F <- c(KF_output$F)
K <- c(KF_output$K)

o_2.1 <- cbind(a, P, v, F, K)
colnames(o_2.1) <- c("a", "P", "v", "F", "K")

# Print output
print(o_2.1)

# Construct confidence intervals for a (new section)
a_lb <- o_2.1[, "a"] - qnorm(0.95)*sqrt(o_2.1[,"P"])
a_ub <- o_2.1[, "a"] + qnorm(0.95)*sqrt(o_2.1[,"P"])

## Plots

par(mfrow=c(2,2))

# (i)
plot(y_nile[-1], main="(i) y_t, a_t", xlab="", ylab="")   # note: leave out initialization
lines(o_2.1[-1, "a"])
lines(a_lb[-1], col="grey")
lines(a_ub[-1], col="grey")

# (ii)
plot(o_2.1[-1, "P"], main="(ii) P_t", type="l", xlab="", ylab="")

# (iii)
plot(o_2.1[-1, "v"], main="(iii) v_t", type="l", xlab="", ylab="")
abline(h = 0, col = "grey")

# (iv)
plot(o_2.1[-1, "F"], main="(iv) F_t", type="l", xlab="", ylab="")

par(mfrow=c(1,1))


################## Testing the KFS function on the Nile data ###################

# Kalman Filter smoother output for fig 2.2
KFS_output <- KFS(y=y_nile, Z=1, T=1, R=1, d=0, c=0, H=15099, Q=1469.1, initialization=list("diffuse"))

# Extract all KFS output separately and combine in matrix
L <- c(KFS_output$L)
r <- c(KFS_output$r)
N <- c(KFS_output$N)
alpha_hat <- c(KFS_output$alpha_hat)
V <- c(KFS_output$V)

o_2.2 <- cbind(L, r, N, alpha_hat, V)
colnames(o_2.2) <- c("L", "r", "N", "alpha_hat", "V")


# Print output
print(o_2.2)

# Construct confidence intervals for a (new section)
alpha_hat_lb <- o_2.2[, "alpha_hat"] - qnorm(0.95)*sqrt(o_2.2[,"V"])
alpha_hat_ub <- o_2.2[, "alpha_hat"] + qnorm(0.95)*sqrt(o_2.2[,"V"])


## Plots

par(mfrow=c(2,2))

# (i)
plot(y_nile[-1], main="(i) y_t, alpha_hat_t", xlab="", ylab="")   # note: leave out initialization
lines(o_2.2[-1,"alpha_hat"])
lines(alpha_hat_lb[-1], col="grey")
lines(alpha_hat_ub[-1], col="grey")

# (ii)
plot(o_2.2[-1,"V"], main="(ii) V_t", type="l", xlab="", ylab="")

# (iii)
plot(o_2.2[,"r"], main="(iii) r_t", type="l", xlab="", ylab="")
abline(h = 0, col = "grey")

# (iv)
plot(o_2.2[,"N"], main="(iv) N_t", type="l", xlab="", ylab="")

par(mfrow=c(1,1))


####################### Testing how a_f relates to a ###########################

output_p <- KF(y=y_nile, Z=1, T=1, R=1, d=0, c=0, H=15099, Q=1469.1, initialization=list("diffuse"), aim="p")
output_f <- KF(y=y_nile, Z=1, T=1, R=1, d=0, c=0, H=15099, Q=1469.1, initialization=list("diffuse"), aim="f")

a <- output_p$a
a_f <- output_f$a_f

# (i)
plot(y_nile[-1], main="(i) y_t, a_t, a_f_t", xlab="")   # note: leave out initialization
lines(a[-1])
lines(a_f[-1], col="blue")

# difference in values
#a[-1]-a_f[-length(a_f)]


################## Testing Quasi Maximum Likelihood function ###################

## LLM

# obtain estimates
LLM_estimates <- MLE_func(y=y_nile, model="LLM")

# print estimates
cat("sigma2_eps:", LLM_estimates[1],
    "\nsigma2_eta:", LLM_estimates[2],
    "\nDK-book sigma2_eps:", 15099,
    "\nDK-book sigma2_eta:", 1469.1,
    "\n")


## SVM

# obtain estimates 
SVM_estimates <- MLE_func(y=x, model="SVM")

# print estimates
cat("omega (Nelder-Mead):", SVM_estimates[1],
    "\nphi (Nelder-Mead):", SVM_estimates[2],
    "\nsigma2_eta (Nelder-Mead):", SVM_estimates[3],
    "\nomega (BFGS):", -0.077092759,
    "\nphi (BFGS):", 0.992279226,
    "\nsigma2_eta (BFGS):", 0.006362378,
    "\n")



## replicate SVM plots

# estimates
omega <- SVM_estimates[1] 
phi <- SVM_estimates[2]   
sigma2_eta <- SVM_estimates[3]  
ksi <- omega/(1-phi)         

# unconditional moments for initialization:
a1 <- omega/(1-phi)
P1 <- sigma2_eta/(1-phi^2)

# Apply Kalman filter
KF_output_h <- KF(y=x, Z=1, R=1, T=phi, d=-1.27, c=omega, H=(pi^(2)/2), Q=sigma2_eta, initialization=list(c(a1, P1)), aim="f")

# define filtered h as h_f
h_f <- c(KF_output_h$a_f)

# Apply Kalman filter smoother
KFS_output_h <- KFS(y=x, Z=1, R=1, T=phi, d=-1.27, c=omega, H=(pi^(2)/2), Q=sigma2_eta, initialization=list(c(a1, P1)))

# define smoothed h as h_s
h_s <- c(KFS_output_h$alpha_hat)


# Nonlinear model: h_tilde
h_f_tilde=h_f-ksi
h_s_tilde=h_s-ksi



### Plots
par(mfrow=c(1,1))

# Plot 1:
plot(x, xlab="", ylab="", pch=20, cex=0.5)
lines(h_f[-length(h_f)], col="blue")
lines(h_s[-length(h_s)], col="red")

# Plot 2:
plot(h_f_tilde[-length(h_f_tilde)], col="blue", type="l", xlab="", ylab="") 
lines(h_s_tilde[-length(h_s_tilde)], col="red")



