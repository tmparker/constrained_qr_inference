# Create p-value plots for the first simulation experiment.
source("../plotting_utils.R")
load("./normal_sim.rda") # loads an object called dat
erep <- 1000 # number of repetitions of data generation & estimation
N <- c(100, 400, 1000)
refsup <- apply(dat$reference_dist, 3, max)
refL2 <- apply(dat$reference_dist, 3, L2fun) # defined in plotting_utils.R

# Collect statistics and p-values
sups <- lapply(dat$est, function(a) {do.call(rbind, lapply(a, 
                                          function(b) {apply(b, 2, max)}))})
L2s <- lapply(dat$est, function(a) {do.call(rbind, lapply(a, 
                                          function(b) {apply(b, 2, L2fun)}))})
pv_sup <- lapply(sups, function(a) {matrix(sapply(a, pvsim, ref = refsup), 
                                    nrow = erep)})
pv_L2 <- lapply(L2s, function(a) {matrix(sapply(a, pvsim, ref = refL2), 
                                    nrow = erep)})
names(pv_sup) <- names(pv_L2) <- paste("n =", N)
lapply(pv_sup, function(M) {colnames(M) <- c("QLR", "W", "RS"); M})
lapply(pv_L2, function(M) {colnames(M) <- c("QLR", "W", "RS"); M})

# Save p-value plots
labvec <- c("rho", "Wald", "rankscore") # differs from printed figures.
pdf(file = "pvplots_supnorm_normal.pdf")
par(mfcol = c(3, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
sapply(1:3, pfun, pvl = pv_sup, lab = labvec)
sapply(1:3, pdiffun, pvl = pv_sup, lab = labvec)
dev.off()

pdf(file = "pvplots_l2_normal.pdf")
par(mfcol = c(3, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
sapply(1:3, pfun, pvl = pv_L2, lab = labvec)
sapply(1:3, pdiffun, pvl = pv_L2, lab = labvec)
dev.off()

