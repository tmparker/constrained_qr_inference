# Inference process for a simple hypothesis.
# Maintained hypothesis: F dominates N(0, 1)
# Test F = N(0,1) against F FOSD N(0,1) subject to FOSD constraint.
# Makes this a type B problem.

library(quantreg)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../const_inf_utils.cpp", cacheDir="../cpp")
set.seed(8675309)

erep <- 1000 # number of repetitions of data generation & estimation
rrep <- 1000 # number of reps used for reference dbn. simulation
N <- c(100, 400, 1000)
Tau <- 5:95 / 100

# Objective function depends on residuals u and quantile level tau.
# u is a vector of residuals
# tau is the (scalar) quantile level
Vtau <- function(u, tau) { sum(u * (tau - (u < 0))) }

# est() estimates a quantile regression model and three test statistics for
# testing the hypothesis that the quantile is equal to the corresponding
# standard normal quantile, when the estimation is constrained so that the
# conditional quantiles stochastically dominate the normal distribution.
# This experiment uses a simple regression model:
# Q is a scalar quantile level (stitched together later to make test statistic
# processes)
# y is the data vector (since this is a one-sample location problem, X is coded
# in as a vector of ones)
# This function returns a likelihood ratio, Wald and score statistic.
est <- function(Q, y) {
  n <- length(y)
  X <- matrix(1, n, 1)
  ureg <- rq.fit.br(X, y, tau = Q)
  cco <- max(ureg$coef, qnorm(Q))
  creg <- list("coef" = cco, "resid" = y - cco)
  Vhat <- Vtau(ureg$resid, Q)
  Vtil <- Vtau(creg$resid, Q)
  rhostat <- 2 * dnorm(qnorm(Q)) * (Vtil - Vhat) / (Q * (1 - Q))
  waldstat <- n * pmax(qnorm(Q) - ureg$coef, 0)^2 *
                    dnorm(qnorm(Q))^2 / (Q * (1 - Q))
  if (abs(creg$coef - qnorm(Q)) < 1e-06) {
    Sn <- sum(Q - (creg$resid <= 0))
  }
  else {
    bn <- Q - (creg$resid <= 0)
    h <- which.min(abs(creg$resid))
    bn[h] <- -sum(Q - (creg$resid[-h] <= 0))
    Sn <- sum(bn)
  }
  scorestat <- Sn^2 / (n * Q * (1 - Q))
  ans <- c(rhostat, waldstat, scorestat)
  ans
}

# Performs one data generating and estimation step.  
# i is a dummy for lapply later
# n is the sample size 
# tau is a vector of quantile levels
one_est <- function(i, ssize, tau) {
  y <- rnorm(ssize)
  ans <- suppressWarnings(t(sapply(tau, est, y = y)))
  dimnames(ans)[[2]] <- c("QLR", "Wald", "score")
  ans
}

# This creates a list of lists.  Top level: sample sizes.  In each sample size
# there are ereps repetitions (different lists) of nt x 3 matrices with rho,
# wald and score processes in each column in that order.
sim_run <- lapply(N, function(nn) {lapply(1:erep, one_est, ssize = nn, 
                                          tau = Tau)})

#####
# Simulate a reference distribution for the process.
# This is a least favorable distribution for a type B problem
# The function makebridges() is in const_inf_utils.cpp.
BB <- makebridges(Tau, rrep, 1)
refproc <- pmin(BB, 0)^2 / (Tau * (1 - Tau))
if (Tau[1] == 0) {
  lapply(refproc, function(a) {a[1] <- 0})
}
if (Tau[length(Tau)] == 1) {
  lapply(refproc, function(a) {a[length(Tau)] <- 0})
}

# Save all this data!
dat <- list("estimations" = sim_run, "reference_dist" = refproc)
save(dat, file = "normal_sim.rda")

