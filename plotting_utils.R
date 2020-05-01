# Functions used to make p-value plots for both simulation experiments.

erep <- 1000 # number of repetitions of data generation & estimation
N <- c(100, 400, 1000)
Tau <- 5:95 / 100

# Construct simulated reference distributions of statistics:
L2fun <- function(a) {
  b <- a^2
  trap <- sqrt(0.5 * sum((b[-1] + b[-nT]) * diff(Tau)))
  return(trap)
}
pvsim <- function(a, ref) {mean(ref > a)}
L2fun <- function(a) {sum(a^2) * diff(Tau)[1]} # assumes even spacing

# p-value plots
# A re-used plotting parameter
colvec <- 1:length(N) + 1 #c("red", "blue", "green")
#labvec <- c("rho process", "Wald process", "rankscore process")
#labvec <- c("rho", "Wald", "rankscore") # differs from printed figures.

# Plot the empirical CDF of p-values against the uniform CDF
# pvl is a list of p-value matrices
# lab is a vector of labels
# num is which test to use (among those in the label vector)
pfun <- function(pvl, lab, num) {
  plot(1, 1, type = "n", xlab = "Uniform p-value", ylab = "Empirical p-values", 
       main = paste("p-value CDFs,", lab[num]), xlim = 0:1, ylim = 0:1)
  for (i in 1:length(N)) {
    pv <- pvl[[i]][, num]
    u <- c(0, unique(sort(pv)))
    yu <- cumsum(tabulate(match(pv, u))) / erep
    lines(u, yu, type = "s", col = colvec[i], lty = i+1)
    abline(0, 1, lwd = 0.25)
  }
  legend("bottomright", lty = 1:length(N) + 1, col = colvec, 
          legend = paste("n =", N), bty = "n")
}

# Plot the difference between the empirical CDF of p-values and the uniform
# Arguments are the same as for pfun
pdiffun <- function(pvl, lab, num) {
  plot(1, 1, type = "n", xlab = "Uniform p-value", ylab = "Emp. minus theo.",
        #ylim = c(min(yu[-length(u)] - u[-1]), max(yu - u)), 
       main = paste("Difference from unif,", lab[num]), xlim = 0:1, 
       ylim = c(-0.05, 0.05))
  for (i in 1:length(N)) {
    pv <- pvl[[i]][, num]
    u <- c(0, unique(sort(pv)))
    yu <- cumsum(tabulate(match(pv, u))) / erep
    lines(u, yu - u, col = colvec[i], lty = i+1)
    abline(h = 0, lwd = 0.25)
  }
  legend("bottomright", lty = 1:length(N) + 1, col = colvec, 
          legend = paste("n =", N), bty = "n")
}

