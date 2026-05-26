## Simulation experiments: Setting III
## Large n version with "_l" suffix

# Functions and libraries ----

source("functions/sim_fit_inhomogeneousHMM.r")
source("functions/ar_approach.r")
source("functions/bb_approach.r")
source("functions/dir_reg.r")
source("functions/sim_study.r")

colour_l = c("orange", "skyblue", "seagreen")

# Parameters ----
set.seed(22)
N_l <- 3
n_l <- 10000 # CHANGED TO LARGE N
num_simulations_l <- 200 # number of simulated datasets 
mu_l <- c(5, 10, 15)
sig_l <- c(3, 3, 4)
periodic_l <- TRUE # periodic covariate
cycle_length_l = 24
base_amplitude_l = 4
amplitude_var_l = 0.6
noise_sd_l = 1

beta_l = matrix(c( -1, -1.5,   # 1-> 2
                   -1, 1.5,  # 1-> 3
                   -1, 0.5,   # 2-> 1
                   -1, 1.5, # 2-> 3
                   -1, 1.5,  # 3-> 1
                   -1, -1), # 3-> 2
                nrow = 6, byrow = TRUE)

par_l = list(
  beta = c(rep(-1, 6), rep(0,6)),
  delta = c(0,1), 
  mu = mu_l, 
  sigma = log(sig_l)
)

sim_l <- simCovHMM(n = n_l, mu = mu_l, sig = sig_l, beta = beta_l, periodic = TRUE)

dat_l = list(x = sim_l$x,
             Z = matrix(sim_l$z), 
             N=3)

# Time series of simulated covariate and states
par(mfrow=c(1,1))
plot(ts(sim_l$z[1:300]), ylab = "Covariate z", xlab = "Time", main = "Simulated periodic covariate (first 300 time steps)", bty="n")
points(ts(sim_l$z[1:300]), col = colour_l[sim_l$s], pch = 20)
legend("topright", legend = c("state 1", "state 2", "state 3"), col = colour_l[1:3], pch = 20, bty="n")

Gamma_l = tpm_g(sim_l$z, beta_l)

Delta_l = matrix(NA, n_l, 3)
Delta_l[1,] = rep(1/3, 3)
for(t in 2:n_l) Delta_l[t,] = Delta_l[t-1,] %*% Gamma_l[,,t]

xseq_l = seq(min(sim_l$z), max(sim_l$z), length = 200)
Gammaseq_l = tpm_g(xseq_l, beta_l)
Deltaseq_l = matrix(NA, length(xseq_l), 3)
for(t in 1:length(xseq_l)) Deltaseq_l[t,] = LaMa::stationary(Gammaseq_l[,,t]) 

# hypothetical stationary distribution (lines) vs. empirical state probs (points)
par(mfrow=c(1,1))
plot(sim_l$z, Delta_l[, 1], 
     xlim = c(min(sim_l$z)-1, max(sim_l$z)+1), ylim = c(0, 1), 
     pch = 16, col = alpha(colour_l[1], 0), 
     bty = "n", ylab = expression(paste("State probabilities ", delta[t])), xlab = expression(z[t]), 
     main="") 

for (state in 1:3) {
  points(sim_l$z, Delta_l[, state], col = alpha(colour_l[state], 0.1), pch = 16)
  lines(xseq_l, Deltaseq_l[, state], col = colour_l[state], lwd = 3)
}
legend("left", legend = paste("state", 1:3), col = colour_l, lwd = 2, bty = "n")


### GT ----
nGT3_l <- 20000000
tGT3_l <- 0:(nGT3_l - 1)
amp_l <- base_amplitude_l + rnorm(nGT3_l, 0, amplitude_var_l)
zGT3_l <- amp_l * sin(2 * pi * tGT3_l / cycle_length_l) + rnorm(nGT3_l, 0, noise_sd_l)

ZGT3_l=matrix(zGT3_l)

GammaGT3_l = tpm_g(ZGT3_l, beta_l)
delta_l = rep(1/3, 3) 

DeltaGT3_l = matrix(NA, nGT3_l, 3)
DeltaGT3_l[1,] = rep(1/3, 3)

for(t in 2:nGT3_l) DeltaGT3_l[t,] = DeltaGT3_l[t-1,] %*% GammaGT3_l[,,t]

z_binsGT3_l <- seq(min(zGT3_l), max(zGT3_l), length.out = 100)
bin_midpointsGT3_l <- (z_binsGT3_l[-1] + z_binsGT3_l[-length(z_binsGT3_l)]) / 2

bin_ids3_l <- cut(zGT3_l, breaks = z_binsGT3_l, include.lowest = TRUE, labels = FALSE)
bin_ids_factor3_l <- factor(bin_ids3_l, levels = 1:(length(z_binsGT3_l) - 1))  

mean_stateprobsGT3_l <- matrix(NA, nrow = length(z_binsGT3_l) - 1, ncol = 3)
q2.5_stateprobsGT3_l <- matrix(NA, nrow = length(z_binsGT3_l) - 1, ncol = 3)
q97.5_stateprobsGT3_l <- matrix(NA, nrow = length(z_binsGT3_l) - 1, ncol = 3)

for (state in 1:3) {
  vals_split_l <- split(DeltaGT3_l[, state], bin_ids_factor3_l)
  for (b in seq_along(vals_split_l)) {
    vals_l <- vals_split_l[[b]]
    if (length(vals_l) > 0 && !all(is.na(vals_l))) {
      mean_stateprobsGT3_l[b, state]  <- mean(vals_l, na.rm = TRUE)
      q2.5_stateprobsGT3_l[b, state]  <- quantile(vals_l, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT3_l[b, state] <- quantile(vals_l, 0.975, na.rm = TRUE)
    }
  }
}

zseqGT3_l = seq(min(zGT3_l), max(zGT3_l), length = 200)
GammaseqGT3_l = tpm_g(zseqGT3_l, beta_l)
DeltaseqGT3_l = matrix(NA, length(zseqGT3_l), 3)
for(t in 1:length(zseqGT3_l)) DeltaseqGT3_l[t,] = LaMa::stationary(GammaseqGT3_l[,,t]) 

# hypothetical stationary distribution (dashed lines) vs. true state occupancy distribution (solid lines)
par(mfrow=c(1,1))
plot(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3, ylab = "Pr(state | z)", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n", xlim=c(-10,10))
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
lines(ksmooth(zseqGT3_l, DeltaseqGT3_l[, 1], "normal", bandwidth = 1), col = alpha(colour_l[1]), lwd = 3, lty=2)
lines(ksmooth(zseqGT3_l, DeltaseqGT3_l[, 2], "normal", bandwidth = 1), col = alpha(colour_l[2]), lwd = 3, lty=2)
lines(ksmooth(zseqGT3_l, DeltaseqGT3_l[, 3], "normal", bandwidth = 1), col = alpha(colour_l[3]), lwd = 3, lty=2)


#### AR and Hypothetical -----
num_covsim_l <- 200
num_simulations_l <- 200

results3_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- simOneArHypothetical(
      n = n_l, rho = rho_l, mu = mu_l, sig = sig_l,
      beta = beta_l, periodic = periodic_l, par = par_l, dat = dat_l,
      num_covsim = num_covsim_l
    )
    p(message = sprintf("Done %d/%d", i, num_simulations_l))
    res
  })
})

curve3_l <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical3_l <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:num_simulations_l) {
  curve3_l$State1[[i]] <- results3_l[[i]]$resultar$State1
  curve3_l$State2[[i]] <- results3_l[[i]]$resultar$State2
  curve3_l$State3[[i]] <- results3_l[[i]]$resultar$State3
  
  hypothetical3_l$State1[[i]] <- results3_l[[i]]$resulthypothetical$State1
  hypothetical3_l$State2[[i]] <- results3_l[[i]]$resulthypothetical$State2
  hypothetical3_l$State3[[i]] <- results3_l[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:num_simulations_l) {
  lines(curve3_l$State1[[i]]$x, curve3_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve3_l$State2[[i]]$x, curve3_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve3_l$State3[[i]]$x, curve3_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)



plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:num_simulations_l) {
  lines(hypothetical3_l$State1[[i]]$x, hypothetical3_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(hypothetical3_l$State2[[i]]$x, hypothetical3_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(hypothetical3_l$State3[[i]]$x, hypothetical3_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)

### Flexible Dirichlet regression ----
num_simulations_l <- 200

gam_results3_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- oneDirGAM(n = n_l, rho=rho_l, mu = mu_l, sig = sig_l, beta = beta_l, periodic = periodic_l, 
                     par = par_l, dat = dat_l, min_pred = -8, max_pred = 8)
    p()
    res
  })
})

trim_to_range_l <- function(x, y, xmin = -8.05, xmax = 8.05) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}
ks1_3_l <- ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1)
cut1_3_l <- trim_to_range_l(ks1_3_l$x, ks1_3_l$y)
ks2_3_l <- ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1)
cut2_3_l <- trim_to_range_l(ks2_3_l$x, ks2_3_l$y)
ks3_3_l <- ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1)
cut3_3_l <- trim_to_range_l(ks3_3_l$x, ks3_3_l$y)


gamcurve3_l <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations_l) {
  gamcurve3_l$State1[[i]] <- gam_results3_l[[i]]$State1
  gamcurve3_l$State2[[i]] <- gam_results3_l[[i]]$State2
  gamcurve3_l$State3[[i]] <- gam_results3_l[[i]]$State3
}


plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:num_simulations_l) {
  lines(gamcurve3_l$State1[[i]]$x, gamcurve3_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(gamcurve3_l$State2[[i]]$x, gamcurve3_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(gamcurve3_l$State3[[i]]$x, gamcurve3_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}

lines(cut1_3_l$x, cut1_3_l$y, col = colour_l[1], lwd = 3)
lines(cut2_3_l$x, cut2_3_l$y, col = colour_l[2], lwd = 3)
lines(cut3_3_l$x, cut3_3_l$y, col = colour_l[3], lwd = 3)
legend("right",col = c(colour_l, "transparent"),lwd = 3,bty = "n",
       lty = c(1,1,1),legend = expression(state~1, state~2, state~3))


lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
legend("right",col = c(colour_l, "transparent", "black", "black"),lwd = 3,bty = "n",
       lty = c(1,1,1, NA, 1, 3),legend = expression(state~1, state~2, state~3, "", delta, rho))



### AR Approach ----

num_covsim_l <- 200

num_simulations_l <- 200

results3AR_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- simOneAr(n = n_l, rho=rho_l, mu = mu_l, sig = sig_l, beta = beta_l, periodic = periodic_l, par=par_l, num_covsim = num_covsim_l)
    p(message = sprintf("Done %d/%d", i, num_simulations_l))
    res
  })
})

curve3AR_l <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations_l) {
  curve3AR_l$State1[[i]] <- results3AR_l[[i]]$State1
  curve3AR_l$State2[[i]] <- results3AR_l[[i]]$State2
  curve3AR_l$State3[[i]] <- results3AR_l[[i]]$State3
}


par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", bty = "n",
     xlim=c(min(zGT3_l)+1, max(zGT3_l)-1), main="AR resampling")

for (i in 1:num_simulations_l) {
  lines(curve3AR_l$State1[[i]]$x, curve3AR_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve3AR_l$State2[[i]]$x, curve3AR_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve3AR_l$State3[[i]]$x, curve3AR_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)

### BB Approach ----
num_covsim_l <- 200
num_simulations_l <- 200

results_3BB_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- simOneBB(n = n_l, rho = rho_l, mu = mu_l, sig = sig_l,
                    beta = beta_l, periodic = periodic_l, par = par_l, dat = dat_l,
                    num_covsim = num_covsim_l, blocklength = 24)
    p(message = sprintf("Done %d/%d", i, num_simulations_l))
    res
  })
})

curve_3BB_l <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations_l) {
  curve_3BB_l$State1[[i]] <- results_3BB_l[[i]]$State1
  curve_3BB_l$State2[[i]] <- results_3BB_l[[i]]$State2
  curve_3BB_l$State3[[i]] <- results_3BB_l[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (III): BB resampling", bty = "n",
     xlim=c(min(zGT3_l)+1, max(zGT3_l)-1))

for (i in 1:num_simulations_l) {
  lines(curve_3BB_l$State1[[i]]$x, curve_3BB_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve_3BB_l$State2[[i]]$x, curve_3BB_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve_3BB_l$State3[[i]]$x, curve_3BB_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_l, mean_stateprobsGT3_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour_l), lwd = 2, bty = "n")