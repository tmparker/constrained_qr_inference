# Simulation for inference processes for a treatment effect.
# Maintained hypotheis \beta_p \geq 0
# Test \beta_p = 0 against \beta_p > 0.
# Makes this a type A problem.

library(quantreg)
#library(coneproj)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../const_inf_utils.cpp", cacheDir="../cpp")
set.seed(8675309)

erep <- 1000
rrep <- 1000
N <- c(100, 400, 1000)
Tau <- seq(0.1, 0.9, by = 0.05)

p <- 4
pm1 <- p - 1
# H0: the pth coefficient coordinate is zero.
R <- matrix(c(rep(0, pm1), 1), 1, p)
r <- 0

# Sum of checked residuals.
Vtau <- function(u, tau) { sum(u * (tau - (u < 0))) }

n <- N[1]
x <- sapply(1:pm1, function(i) rnorm(n))
y <- 1 + .rowSums(x, n, pm1-1) + rnorm(n) # last regressor irrelevant
X <- cbind(1, x)
tau <- Tau

# This est() is constructed for this particular R and r, even though they are
# arguments to the function.
est <- function(tau, x, y, R, r) {
  n <- length(y)
  nt <- length(tau)
  X <- cbind(1, x)
  p <- dim(X)[[2]]
  X0 <- X[, -p] # tailored to the hypotheses.
  XR <- rbind(X, R)

  # Like the summary.rq() default but scaled by n.
  Dn <- crossprod(X) / n
  xxinv <- diag(p)
  xxinv <- backsolve(qr(X)$qr[1:p, 1:p, drop = FALSE], xxinv)
  Dinv <- n * tcrossprod(xxinv)

  #####
  # creg is inequality constrained and ereg is equality constrained.
  ereg <- rq(y ~ x[, 1:(pm1 - 1)], tau = Tau)
  creg <- rq(y ~ x, method = "fnc", R = R, r = r, tau = Tau)
  esum <- summary(ereg, se = "iid", covariance = TRUE)
  csum <- summary(creg, se = "iid", covariance = TRUE)

  # Pieces taken from rq() and summary.rq()
  eco <- cco <- matrix(0, p, nt)
  Vbar <- Vtil <- rep(0, nt)
  ubar <- util <- matrix(0, n, nt)
  spar <- double(nt)
  elam <- matrix(0, n, nt)
  clam <- matrix(0, n + 1, nt)
  for (i in 1:nt) {
    ereg <- rq.fit.br(X[, -p], y, tau = tau[i])
    creg <- rq.fit.fnc(X, y, R, r, tau = tau[i])
    eco[1:pm1, i] <- ereg$coef
    cco[, i] <- creg$coef
    ubar[, i] <- ereg$resid
    util[, i] <- creg$resid
    Vbar[i] <- Vtau(ereg$resid, tau[i])
    Vtil[i] <- Vtau(creg$resid, tau[i])

    # Suggestion from Koenker & Machado (1999) - almost option 'se = "nid"' in
    # summary.rq().  A little easier because in this example the data is assumed
    # to be iid, not nid for this experiment (to demonstrate all the inference
    # processes together, and the QLR statistics demand iid data).
    xbar <- colMeans(X)
    h <- bandwidth.rq(tau[i], n, hs = TRUE)
    if (tau[i] + h > 1) 
      stop("tau + h > 1:  error in summary.rq")
    if (tau[i] - h < 0) 
      stop("tau - h < 0:  error in summary.rq")
    bhi <- rq.fit.fnc(X, y, R = R, r = r, tau = tau[i] + h)$coef
    blo <- rq.fit.fnc(X, y, R = R, r = r, tau = tau[i] - h)$coef
    spar[i] <- xbar %*% (bhi - blo) / (2 * h)

    # Rankscore construction
    elam[, i] <- tau[i] - (ereg$resid < 0)
    he <- smallind(abs(ereg$resid), pm1) # C++ code
    elam[he, i] <- -tcrossprod(solve(t(X0[he, ])), X0[-he, ]) %*% elam[-he, i]
    clam[1:n, i] <- tau[i] - (creg$resid < 0)
    if (abs(creg$coef[p]) < 1e-06) { # Relies on this particular null
      hc <- c(smallind(abs(creg$resid), pm1), n + 1) # C++ code
    } else {
      hc <- smallind(abs(creg$resid), p) # C++ code
    }
    clam[hc, i] <- -tcrossprod(solve(t(XR[hc, ])), XR[-hc, ]) %*% clam[-hc, i]
  }

  Sigma <- rep(tau * (1 - tau) * spar^2, each = p^2) * replicate(nt, Dinv)
  bbar <- elam + 1 - rep(Tau, each = n)
  btil <- clam[1:n, ] + 1 - rep(Tau, each = n)
  Sndiff <- c(Dinv[p, ] %*% crossprod(X, bbar - btil) / sqrt(n))

  #####
  # Inference processes
  Ln <- 2 * (Vbar - Vtil) / (tau * (1 - tau) * spar)
  Lamn <- 2 * Vtil * log(Vbar / Vtil) / (tau * (1 - tau) * spar)
  Wn <- n * cco[p, ]^2 / Sigma[p, p, ]
  Tn <- Sndiff^2 / (tau * (1 - tau) * Dn[p, p])
  testprocs <- cbind(Ln, Lamn, Wn, Tn)
  return(list(proc = testprocs, cfun = Sigma))
}

# Performs one data generating and estimation step.  i is a dummy for lapply
# later, n is the sample size and Tau is a vector of quantile levels.
one_est <- function(i, ssize, tau) {
  x <- sapply(1:pm1, function(i) rnorm(ssize))
  y <- 1 + .rowSums(x, ssize, pm1-1) + rnorm(ssize) # last regressor irrelevant
  ans <- suppressWarnings(est(tau, x, y, R, r))
  ans
}
#tst <- one_est(1, N[1], Tau)

# This creates a list of lists.  Top level: sample sizes.  In each sample size
# there are ereps repetitions (different lists) of nt x 4 matrices with rho (two
# kinds), wald and rankscore processes in each column in that order.
sim_run <- lapply(N, function(nn) {lapply(1:erep, one_est, ssize = nn, 
                                          tau = Tau)})
# The reference distributions will depend on some Brownian bridges.
BB <- makebridges(Tau, rrep, 1)

# Save all this data!
dat <- list("estimates" = sim_run, "bridges" = BB)
save(dat, file="oneregressor_treatment_sim.rda")

