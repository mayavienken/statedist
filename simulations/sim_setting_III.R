## Simulation experiments: Setting III

# Functions and libraries ----

source("functions/sim_fit_inhomogeneousHMM.r")
source("functions/ar_approach.r")
source("functions/bb_approach.r")
source("functions/dir_reg.r")

colour = c("orange", "skyblue", "seagreen")

# Parameters ----
N <- 3
n <- 2000
num_simulations <- 200 # number of simulated datasets 
mu <- c(5, 10, 15)
sig <- c(3, 3, 4)
periodic <- TRUE # periodic covariate
cycle_length = 24
base_amplitude = 4
amplitude_var = 0.6
noise_sd = 1

beta = matrix(c( -1, -1.5,   # 1-> 2
                 -1, 1.5,  # 1-> 3
                 -1, 0.5,   # 2-> 1
                 -1, 1.5, # 2-> 3
                 -1, 1.5,  # 3-> 1
                 -1, -1), # 3-> 2
              nrow = 6, byrow = TRUE)

par = c(beta = c(rep(-1, 6), rep(0,6)),
        logitdelta = c(0,0), 
        mu = mu, 
        sig = log(sig)) 

sim <- simCovHMM(n = n, mu = mu, sig = sig, beta = beta, periodic = TRUE)

# Time series of simulated covariate and states
par(mfrow=c(1,1))
plot(ts(sim$z[1:300]), ylab = "Covariate z", xlab = "Time", main = "Simulated periodic covariate (first 300 time steps)", bty="n")
points(ts(sim$z[1:300]), col = colour[sim$s], pch = 20)
legend("topright", legend = c("state 1", "state 2", "state 3"), col = colour[1:3], pch = 20, bty="n")

Gamma = tpm_g(sim$z, beta)

Delta = matrix(NA, n, 3)
Delta[1,] = rep(1/3, 3)
for(t in 2:n) Delta[t,] = Delta[t-1,] %*% Gamma[,,t]

xseq = seq(min(sim$z), max(sim$z), length = 200)
Gammaseq = tpm_g(xseq, beta)
Deltaseq = matrix(NA, length(xseq), 3)
for(t in 1:length(xseq)) Deltaseq[t,] = LaMa::stationary(Gammaseq[,,t]) 

# hypothetical stationary distribution (lines) vs. empirical state probs (points)
par(mfrow=c(1,1))
plot(sim$z, Delta[, 1], 
     xlim = c(min(sim$z)-1, max(sim$z)+1), ylim = c(0, 1), 
     pch = 16, col = alpha(colour[1], 0), 
     bty = "n", ylab = expression(paste("State probabilities ", delta[t])), xlab = expression(z[t]), 
     main="") 

for (state in 1:3) {
  points(sim$z, Delta[, state], col = alpha(colour[state], 0.1), pch = 16)
  lines(xseq, Deltaseq[, state], col = colour[state], lwd = 3)
}
legend("left", legend = paste("state", 1:3), col = colour, lwd = 2, bty = "n")


### GT ----
nGT3 <- 20000000
tGT3 <- 0:(nGT3 - 1)
amp <- base_amplitude + rnorm(nGT3, 0, amplitude_var)
zGT3 <- amp * sin(2 * pi * tGT3 / cycle_length) + rnorm(nGT3, 0, noise_sd)

ZGT3=matrix(zGT3)

GammaGT3 = tpm_g(ZGT3, beta)
delta = rep(1/3, 3) 

DeltaGT3 = matrix(NA, nGT3, 3)
DeltaGT3[1,] = rep(1/3, 3)

for(t in 2:nGT3) DeltaGT3[t,] = DeltaGT3[t-1,] %*% GammaGT3[,,t]

z_binsGT3 <- seq(min(zGT3), max(zGT3), length.out = 100)
bin_midpointsGT3 <- (z_binsGT3[-1] + z_binsGT3[-length(z_binsGT3)]) / 2

bin_ids3 <- cut(zGT3, breaks = z_binsGT3, include.lowest = TRUE, labels = FALSE)
bin_ids_factor3 <- factor(bin_ids3, levels = 1:(length(z_binsGT3) - 1))  

mean_stateprobsGT3 <- matrix(NA, nrow = length(z_binsGT3) - 1, ncol = 3)
q2.5_stateprobsGT3 <- matrix(NA, nrow = length(z_binsGT3) - 1, ncol = 3)
q97.5_stateprobsGT3 <- matrix(NA, nrow = length(z_binsGT3) - 1, ncol = 3)

for (state in 1:3) {
  vals_split <- split(DeltaGT3[, state], bin_ids_factor3)
  for (b in seq_along(vals_split)) {
    vals <- vals_split[[b]]
    if (length(vals) > 0 && !all(is.na(vals))) {
      mean_stateprobsGT3[b, state]  <- mean(vals, na.rm = TRUE)
      q2.5_stateprobsGT3[b, state]  <- quantile(vals, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT3[b, state] <- quantile(vals, 0.975, na.rm = TRUE)
    }
  }
}

zseqGT3 = seq(min(zGT3), max(zGT3), length = 200)
GammaseqGT3 = tpm_g(zseqGT3, beta)
DeltaseqGT3 = matrix(NA, length(zseqGT3), 3)
for(t in 1:length(zseqGT3)) DeltaseqGT3[t,] = LaMa::stationary(GammaseqGT3[,,t]) 

# hypothetical stationary distribution (dashed lines) vs. true state occupancy distribution (solid lines)
par(mfrow=c(1,1))
plot(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3, ylab = "Pr(state i|z), i=1,2,3", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n", xlim=c(-10,10))
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
lines(ksmooth(zseqGT3, DeltaseqGT3[, 1], "normal", bandwidth = 1), col = alpha(colour[1]), lwd = 3, lty=2)
lines(ksmooth(zseqGT3, DeltaseqGT3[, 2], "normal", bandwidth = 1), col = alpha(colour[2]), lwd = 3, lty=2)
lines(ksmooth(zseqGT3, DeltaseqGT3[, 3], "normal", bandwidth = 1), col = alpha(colour[3]), lwd = 3, lty=2)


#### AR and Hypothetical -----
num_covsim <- 200
num_simulations <- 20

results3 <- with_progress({
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

curve3 <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical3 <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:num_simulations) {
  curve3$State1[[i]] <- results3[[i]]$resultar$State1
  curve3$State2[[i]] <- results3[[i]]$resultar$State2
  curve3$State3[[i]] <- results3[[i]]$resultar$State3
  
  hypothetical3$State1[[i]] <- results3[[i]]$resulthypothetical$State1
  hypothetical3$State2[[i]] <- results3[[i]]$resulthypothetical$State2
  hypothetical3$State3[[i]] <- results3[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:num_simulations) {
  lines(curve3$State1[[i]]$x, curve3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve3$State2[[i]]$x, curve3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve3$State3[[i]]$x, curve3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)



plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:num_simulations) {
  lines(hypothetical3$State1[[i]]$x, hypothetical3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical3$State2[[i]]$x, hypothetical3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical3$State3[[i]]$x, hypothetical3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

### Flexible Dirichlet regression ----
source("00_GamApproach.R")
num_simulations <- 200

gam_results3 <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- oneDirGAM(n = n, rho=rho, mu = mu, sig = sig, beta = beta, periodic = periodic, par=par)
    p()
    res
  })
})

trim_to_range <- function(x, y, xmin = -8.05, xmax = 8.05) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}
ks1_3 <- ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1)
cut1_3 <- trim_to_range(ks1_3$x, ks1_3$y)
ks2_3 <- ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1)
cut2_3 <- trim_to_range(ks2_3$x, ks2_3$y)
ks3_3 <- ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1)
cut3_3 <- trim_to_range(ks3_3$x, ks3_3$y)


gamcurve3 <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations) {
  gamcurve3$State1[[i]] <- gam_results3[[i]]$State1
  gamcurve3$State2[[i]] <- gam_results3[[i]]$State2
  gamcurve3$State3[[i]] <- gam_results3[[i]]$State3
}


plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:num_simulations) {
  lines(gamcurve3$State1[[i]]$x, gamcurve3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve3$State2[[i]]$x, gamcurve3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve3$State3[[i]]$x, gamcurve3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(cut1_3$x, cut1_3$y, col = colour[1], lwd = 3)
lines(cut2_3$x, cut2_3$y, col = colour[2], lwd = 3)
lines(cut3_3$x, cut3_3$y, col = colour[3], lwd = 3)
legend("right",col = c(colour, "transparent"),lwd = 3,bty = "n",
       lty = c(1,1,1),legend = expression(state~1, state~2, state~3))


lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
legend("right",col = c(colour, "transparent", "black", "black"),lwd = 3,bty = "n",
       lty = c(1,1,1, NA, 1, 3),legend = expression(state~1, state~2, state~3, "", delta, rho))



### AR Approach ----

num_covsim <- 200

num_simulations <- 200

results3AR <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- simOneAr(n = n, rho=rho, mu = mu, sig = sig, beta = beta, periodic = periodic, par=par, num_covsim = num_covsim)
    p(message = sprintf("Done %d/%d", i, num_simulations))
    res
  })
})

curve3AR <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations) {
  curve3AR$State1[[i]] <- results3AR[[i]]$State1
  curve3AR$State2[[i]] <- results3AR[[i]]$State2
  curve3AR$State3[[i]] <- results3AR[[i]]$State3
}


par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i|z), i=1,2,3", bty = "n",
     xlim=c(min(zGT3)+1, max(zGT3)-1), main="AR resampling")

for (i in 1:num_simulations) {
  lines(curve3AR$State1[[i]]$x, curve3AR$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve3AR$State2[[i]]$x, curve3AR$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve3AR$State3[[i]]$x, curve3AR$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 5), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 5), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 5), col = colour[3], lwd = 3)

### BB Approach ----
num_covsim <- 2000
num_simulations <- 200

results_3BB <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- simOneBB(n = n, rho = rho, mu = mu, sig = sig,
                    beta = beta, periodic = periodic, par = par,
                    num_covsim = num_covsim)
    p(message = sprintf("Done %d/%d", i, num_simulations))
    res
  })
})

curve_3BB <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations) {
  curve_3BB$State1[[i]] <- results_3BB[[i]]$State1
  curve_3BB$State2[[i]] <- results_3BB[[i]]$State2
  curve_3BB$State3[[i]] <- results_3BB[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (III): BB resampling", bty = "n",
     xlim=c(min(zGT3)+1, max(zGT3)-1))

for (i in 1:num_simulations) {
  lines(curve_3BB$State1[[i]]$x, curve_3BB$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve_3BB$State2[[i]]$x, curve_3BB$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve_3BB$State3[[i]]$x, curve_3BB$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour), lwd = 2, bty = "n")

