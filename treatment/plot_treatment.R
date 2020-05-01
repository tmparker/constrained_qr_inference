# Create p-value plots for the second simulation experiment.
source("../plotting_utils.R")
load("./oneregressor_treatment_sim.rda") # loads an object called dat
erep <- 1000
N <- c(100, 400, 1000)
Tau <- seq(0.1, 0.9, by = 0.05)

# Simulate a reference distribution for the process.
refproc <- pmax(dat$bridges, 0)^2 / (Tau * (1 - Tau))
if (Tau[1] == 0) {
  lapply(refproc, function(a) {a[1] <- 0})
}
if (Tau[length(Tau)] == 1) {
  lapply(refproc, function(a) {a[length(Tau)] <- 0})
}
refsup <- apply(refproc, 3, max)
refL2 <- apply(refproc, 3, L2fun) # defined in plotting_utils.R

# Collect statistics ans p-values
sim_pro <- lapply(dat$estimates, function(x) lapply(x, "[[", 1))
sups <- lapply(sim_pro, function(a) {do.call(rbind, lapply(a, 
                                          function(b) {apply(b, 2, max)}))})
L2s <- lapply(sim_pro, function(a) {do.call(rbind, lapply(a, 
                                          function(b) {apply(b, 2, L2fun)}))})
pv_sup <- lapply(sups, function(a) {matrix(sapply(a, pvsim, ref = refsup), 
                                    nrow = erep)})
pv_L2 <- lapply(L2s, function(a) {matrix(sapply(a, pvsim, ref = refL2), 
                                    nrow = erep)})
names(pv_sup) <- names(pv_L2) <- paste("n =", N)
nlist <- colnames(dat$est[[1]][[1]]$proc)
lapply(pv_sup, function(M) {colnames(M) <- nlist; M})
lapply(pv_L2, function(M) {colnames(M) <- nlist; M})

# Save p-value plots
labvec <- c("L_n", "Lambda_n", "W_n", "T_n")

pdf(file = "pvplots_supnorm_sim2.pdf")
par(mfcol = c(4, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
sapply(1:4, pfun, pvl = pv_sup, lab = labvec)
sapply(1:4, pdiffun, pvl = pv_sup, lab = labvec)
dev.off()

pdf(file = "pvplots_l2_sim2.pdf")
par(mfcol = c(4, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
sapply(1:4, pfun, pvl = pv_L2, lab = labvec)
sapply(1:4, pdiffun, pvl = pv_L2, lab = labvec)
dev.off()

