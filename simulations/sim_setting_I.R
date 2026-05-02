## Simulation experiments: Setting I

# Functions and libraries ----

source("functions/sim_fit_inhomogeneousHMM.r")
source("functions/ar_approach.r")
source("functions/bb_approach.r")
source("functions/dir_reg.r")
source("functions/sim_study.r")

colour = c("orange", "skyblue", "seagreen")

# Parameters ----
N <- 3 
n <- 2000 
num_simulations <- 200 # number of simulated datasets 
rho <- 0.95 # highly persistent covariate
epsilon <- sqrt(1 - rho^2)
periodic <- FALSE


mu <- c(6, 15, 20)
sig = c(3, 1.5, 1.5)

beta <- matrix(c( 1, 2,   # 1-> 2
                 1, -2.5,  # 1-> 3
                 -1, 0.5,   # 2-> 1
                 -1, -2.5, # 2-> 3
                 1, -1,  # 3-> 1
                 1, 1), # 3-> 2
              nrow = 6, byrow = TRUE)

par = c(beta = c(rep(-1, 6), rep(0,6)),
        logitdelta = c(0,1), 
        mu = mu, 
        sig = sig)


# Simulation
sim <-  simCovHMM(n=n, rho=rho, mu=mu, sig=sig, beta=beta, periodic=periodic)
par = c(beta = c(rep(-1, 6), rep(0,6)),
        logitdelta = c(0,1), 
        mu = c(6, 15, 20), 
        sig = c(log(3),log(3),log(3)))
fit <- fitCovHMM(par=par, x=sim$x, Z=matrix(sim$z))

Gamma <- tpm_g(sim$z, beta)
Delta <- matrix(NA, n, N)
Delta[1,] <- rep(1/N, N)
for (t in 2:n) Delta[t,] <- Delta[t-1,] %*% Gamma[,,t]


### GROUNDTRUTH (T=10,000,000) ----
n_z = 10000000
zGT2 = numeric(n_z)
zGT2[1] = rnorm(1)

for(t in 2:n_z){ zGT2[t] = rho * zGT2[t - 1] + rnorm(1, sd = sqrt(1 - rho^2))}
zGT2=matrix(zGT2)

GammaGT = tpm_g(zGT2, beta)
delta = rep(1/3, 3) 

DeltaGT = matrix(NA, n_z, 3)
DeltaGT[1,] = rep(1/3, 3)
for(t in 2:n_z) DeltaGT[t,] = DeltaGT[t-1,] %*% GammaGT[,,t]

z_binsGT2 <- seq(min(zGT2), max(zGT2), length.out = 100)
bin_midpointsGT2 <- (z_binsGT2[-1] + z_binsGT2[-length(z_binsGT2)]) / 2

bin_ids2 <- cut(zGT2, breaks = z_binsGT2, include.lowest = TRUE, labels = FALSE)
bin_ids2_factor <- factor(bin_ids2, levels = 1:(length(z_binsGT2) - 1))  

mean_stateprobsGT2 <- matrix(NA, nrow = length(z_binsGT2) - 1, ncol = 3)

counts2 <- as.vector(table(bin_ids2_factor))  

for (state in 1:3) {
  sums <- rowsum(DeltaGT[, state], bin_ids2_factor)
  mean_stateprobsGT2[, state] <- NA  
  nonzero_bins <- counts2 > 0
  mean_stateprobsGT2[nonzero_bins, state] <- sums[nonzero_bins] / counts2[nonzero_bins]
}

mean_stateprobsGT2 <- matrix(NA, nrow = length(z_binsGT2) - 1, ncol = 3)
q2.5_stateprobsGT2 <- matrix(NA, nrow = length(z_binsGT2) - 1, ncol = 3)
q97.5_stateprobsGT2 <- matrix(NA, nrow = length(z_binsGT2) - 1, ncol = 3)

for (state in 1:3) {
  vals_split <- split(DeltaGT[, state], bin_ids2_factor)
  for (b in seq_along(vals_split)) {
    vals <- vals_split[[b]]
    if (length(vals) > 0 && !all(is.na(vals))) {
      mean_stateprobsGT2[b, state]  <- mean(vals, na.rm = TRUE)
      q2.5_stateprobsGT2[b, state]  <- quantile(vals, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT2[b, state] <- quantile(vals, 0.975, na.rm = TRUE)
    }
  }
}

mean_stateprobsGT2 <- na.approx(mean_stateprobsGT2, x = bin_midpointsGT2, na.rm = FALSE)
zseqGT2 = seq(min(zGT2), max(zGT2), length = 200)
GammaseqGT2 = tpm_g(zseqGT2, beta)
DeltaseqGT2 = matrix(NA, length(zseqGT2), 3)
for(t in 1:length(zseqGT2)) DeltaseqGT2[t,] = LaMa::stationary(GammaseqGT2[,,t]) 

# hypothetical stationary distribution (dashed lines) vs. true state occupancy distribution (solid lines)
#pdf("./simulations/figures/true_state_occup_dist_I.pdf", width=7, height=4)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth=1), col = colour[1], lwd = 3, ylab = "Pr(state i | z), i=1,2,3", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n", xlim=c(-5,5.5))
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
lines(ksmooth(zseqGT2, DeltaseqGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT2, DeltaseqGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT2, DeltaseqGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3, lty=2)
legend("right",col = c(colour, "transparent", "black", "black"),
       lwd = 3,bty = "n",lty = c(1,1,1, NA, 1, 3),legend = expression(
       state~1, state~2, state~3, "", delta, rho))
#dev.off()


### AR -----

num_covsim <- 200
num_simulations <- 200

results <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- simOneAr(n = n, rho=rho, mu = mu, sig = sig, beta = beta, periodic = periodic, par=par, num_covsim = num_covsim)
    p(message = sprintf("Done %d/%d", i, num_simulations))
    res
  })
})

curve <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations) {
  curve$State1[[i]] <- results[[i]]$State1
  curve$State2[[i]] <- results[[i]]$State2
  curve$State3[[i]] <- results[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i|z), i=1,2,3", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve$State1[[i]]$x, curve$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve$State2[[i]]$x, curve$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve$State3[[i]]$x, curve$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
lines(ksmooth(zseqGT2, DeltaseqGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT2, DeltaseqGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT2, DeltaseqGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3, lty=2)


legend("right",
       col = c(colour, "transparent", "black", "black"),
       lwd = 3,
       bty = "n",
       lty = c(1,1,1, NA, 1, 3),
       legend = expression(
         state~1, state~2, state~3, 
         "", 
         delta, 
         rho
       ))

#### AR and Hypothetical  -----
num_covsim <- 200
num_simulations <- 200

results2 <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- simOneArHypothetical(
      n = n, rho = rho, mu = mu, sig = sig,
      beta = beta, periodic = periodic, par = par,
      num_covsim = num_covsim
    )
    p(message = sprintf("Done %d/%d", i, num_simulations))
    res
  })
})


curve2 <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical2 <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:num_simulations) {
  curve2$State1[[i]] <- results2[[i]]$resultar$State1
  curve2$State2[[i]] <- results2[[i]]$resultar$State2
  curve2$State3[[i]] <- results2[[i]]$resultar$State3
  
  hypothetical2$State1[[i]] <- results2[[i]]$resulthypothetical$State1
  hypothetical2$State2[[i]] <- results2[[i]]$resulthypothetical$State2
  hypothetical2$State3[[i]] <- results2[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main="Setting (I): AR resampling", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve2$State1[[i]]$x, curve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve2$State2[[i]]$x, curve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve2$State3[[i]]$x, curve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
legend("topleft",col = colour, lwd = 3,bty = "n",legend = expression(state~1, state~2, state~3), cex=1.3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "Setting (I): hypothetical stationary distribution", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(hypothetical2$State1[[i]]$x, hypothetical2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical2$State2[[i]]$x, hypothetical2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical2$State3[[i]]$x, hypothetical2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)



### Flexible Dirichlet regression ----

gam_results2 <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- oneDirGAM(n = n, rho=rho, mu = mu, sig = sig, beta = beta, periodic = periodic, par=par)
    p()
    res
  })
})

gamcurve2 <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations) {
  gamcurve2$State1[[i]] <- gam_results2[[i]]$State1
  gamcurve2$State2[[i]] <- gam_results2[[i]]$State2
  gamcurve2$State3[[i]] <- gam_results2[[i]]$State3
}

trim_to_range <- function(x, y, xmin=min(zGT2), xmax=max(zGT2)) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}

ks1 <- ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1)
cut1 <- trim_to_range(ks1$x, ks1$y)
ks2 <- ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1)
cut2 <- trim_to_range(ks2$x, ks2$y)
ks3 <- ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1)
cut3 <- trim_to_range(ks3$x, ks3$y)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i|z), i=1,2,3", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve2$State1[[i]]$x, gamcurve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve2$State2[[i]]$x, gamcurve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve2$State3[[i]]$x, gamcurve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(cut1$x, cut1$y, col = colour[1], lwd = 3)
lines(cut2$x, cut2$y, col = colour[2], lwd = 3)
lines(cut3$x, cut3$y, col = colour[3], lwd = 3)
