## Simulation experiments: Setting III
## Alternative n version with "_s" suffix and error handling

# Functions and libraries ----
set.seed(22)
source("functions/sim_fit_inhomogeneousHMM.r")
source("functions/ar_approach.r")
source("functions/bb_approach.r")
source("functions/dir_reg.r")
source("functions/sim_study.r")

colour_s = c("orange", "skyblue", "seagreen")

# Parameters ----
N_s <- 3
n_s <- 500 # dimension of simulated dataset (length of time series) - CHANGED
num_simulations_s <- 200 # number of simulated datasets 
mu_s <- c(5, 10, 15)
sig_s <- c(3, 3, 4)
periodic_s <- TRUE # periodic covariate
cycle_length_s = 24
base_amplitude_s = 4
amplitude_var_s = 0.6
noise_sd_s = 1

beta_s = matrix(c( -1, -1.5,   # 1-> 2
                   -1, 1.5,  # 1-> 3
                   -1, 0.5,   # 2-> 1
                   -1, 1.5, # 2-> 3
                   -1, 1.5,  # 3-> 1
                   -1, -1), # 3-> 2
                nrow = 6, byrow = TRUE)

par_s = list(
  beta = c(rep(-1, 6), rep(0,6)),
  delta = c(0,1), 
  mu = c(6, 15, 20), 
  sigma = log(sig_s)
)

# Simulation
sim_s <-  simCovHMM(n=n_s, rho=rho_s, mu=mu_s, sig=sig_s, beta=beta_s, periodic=periodic_s)

dat_s = list(
  x = sim_s$x,
  Z = matrix(sim_s$z), 
  N=3
)

fit_s <- fitCovHMM(par=par_s, dat= dat_s)
# Time series of simulated covariate and states
par(mfrow=c(1,1))
plot(ts(sim_s$z[1:300]), ylab = "Covariate z", xlab = "Time", main = "Simulated periodic covariate (first 300 time steps)", bty="n")
points(ts(sim_s$z[1:300]), col = colour_s[sim_s$s], pch = 20)
legend("topright", legend = c("state 1", "state 2", "state 3"), col = colour_s[1:3], pch = 20, bty="n")

Gamma_s = tpm_g(sim_s$z, beta_s)

Delta_s = matrix(NA, n_s, 3)
Delta_s[1,] = rep(1/3, 3)
for(t in 2:n_s) Delta_s[t,] = Delta_s[t-1,] %*% Gamma_s[,,t]

xseq_s = seq(min(sim_s$z), max(sim_s$z), length = 200)
Gammaseq_s = tpm_g(xseq_s, beta_s)
Deltaseq_s = matrix(NA, length(xseq_s), 3)
for(t in 1:length(xseq_s)) Deltaseq_s[t,] = LaMa::stationary(Gammaseq_s[,,t]) 

# hypothetical stationary distribution (lines) vs. empirical state probs (points)
par(mfrow=c(1,1))
plot(sim_s$z, Delta_s[, 1], 
     xlim = c(min(sim_s$z)-1, max(sim_s$z)+1), ylim = c(0, 1), 
     pch = 16, col = alpha(colour_s[1], 0), 
     bty = "n", ylab = expression(paste("State probabilities ", delta[t])), xlab = expression(z[t]), 
     main="") 

for (state in 1:3) {
  points(sim_s$z, Delta_s[, state], col = alpha(colour_s[state], 0.1), pch = 16)
  lines(xseq_s, Deltaseq_s[, state], col = colour_s[state], lwd = 3)
}
legend("left", legend = paste("state", 1:3), col = colour_s, lwd = 2, bty = "n")


### GT ----
nGT3_s <- 20000000
tGT3_s <- 0:(nGT3_s - 1)
amp_s <- base_amplitude_s + rnorm(nGT3_s, 0, amplitude_var_s)
zGT3_s <- amp_s * sin(2 * pi * tGT3_s / cycle_length_s) + rnorm(nGT3_s, 0, noise_sd_s)

ZGT3_s=matrix(zGT3_s)

GammaGT3_s = tpm_g(ZGT3_s, beta_s)
delta_s = rep(1/3, 3) 

DeltaGT3_s = matrix(NA, nGT3_s, 3)
DeltaGT3_s[1,] = rep(1/3, 3)

for(t in 2:nGT3_s) DeltaGT3_s[t,] = DeltaGT3_s[t-1,] %*% GammaGT3_s[,,t]

z_binsGT3_s <- seq(min(zGT3_s), max(zGT3_s), length.out = 100)
bin_midpointsGT3_s <- (z_binsGT3_s[-1] + z_binsGT3_s[-length(z_binsGT3_s)]) / 2

bin_ids3_s <- cut(zGT3_s, breaks = z_binsGT3_s, include.lowest = TRUE, labels = FALSE)
bin_ids_factor3_s <- factor(bin_ids3_s, levels = 1:(length(z_binsGT3_s) - 1))  

mean_stateprobsGT3_s <- matrix(NA, nrow = length(z_binsGT3_s) - 1, ncol = 3)
q2.5_stateprobsGT3_s <- matrix(NA, nrow = length(z_binsGT3_s) - 1, ncol = 3)
q97.5_stateprobsGT3_s <- matrix(NA, nrow = length(z_binsGT3_s) - 1, ncol = 3)

for (state in 1:3) {
  vals_split_s <- split(DeltaGT3_s[, state], bin_ids_factor3_s)
  for (b in seq_along(vals_split_s)) {
    vals_s <- vals_split_s[[b]]
    if (length(vals_s) > 0 && !all(is.na(vals_s))) {
      mean_stateprobsGT3_s[b, state]  <- mean(vals_s, na.rm = TRUE)
      q2.5_stateprobsGT3_s[b, state]  <- quantile(vals_s, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT3_s[b, state] <- quantile(vals_s, 0.975, na.rm = TRUE)
    }
  }
}

zseqGT3_s = seq(min(zGT3_s), max(zGT3_s), length = 200)
GammaseqGT3_s = tpm_g(zseqGT3_s, beta_s)
DeltaseqGT3_s = matrix(NA, length(zseqGT3_s), 3)
for(t in 1:length(zseqGT3_s)) DeltaseqGT3_s[t,] = LaMa::stationary(GammaseqGT3_s[,,t]) 

# hypothetical stationary distribution (dashed lines) vs. true state occupancy distribution (solid lines)
par(mfrow=c(1,1))
plot(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3, ylab = "Pr(state i|z), i=1,2,3", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n", xlim=c(-10,10))
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
lines(ksmooth(zseqGT3_s, DeltaseqGT3_s[, 1], "normal", bandwidth = 1), col = alpha(colour_s[1]), lwd = 3, lty=2)
lines(ksmooth(zseqGT3_s, DeltaseqGT3_s[, 2], "normal", bandwidth = 1), col = alpha(colour_s[2]), lwd = 3, lty=2)
lines(ksmooth(zseqGT3_s, DeltaseqGT3_s[, 3], "normal", bandwidth = 1), col = alpha(colour_s[3]), lwd = 3, lty=2)


#### AR and Hypothetical with error handling -----
num_covsim_s <- 200
num_simulations_s <- 200

results3_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneArHypothetical(
        n = n_s, rho = rho_s, mu = mu_s, sig = sig_s,
        beta = beta_s, periodic = periodic_s, par = par_s,dat = dat_s, 
        num_covsim = num_covsim_s, min_z = -8, max_z = 8
      )
    }, error = function(e) {
      message(sprintf("Simulation %d failed: %s", i, e$message))
      NULL
    })
    
    if(!is.null(res)) {
      successful_results[[j]] <- res
      j <- j + 1
    }
    
    p(message = sprintf("Done %d/%d (%d successful)", i, num_simulations_s, j-1))
  }
  
  successful_results
})

message(sprintf("AR+Hypothetical: Completed %d out of %d simulations successfully", 
                length(results3_s), num_simulations_s))

curve3_s <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical3_s <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:length(results3_s)) {
  curve3_s$State1[[i]] <- results3_s[[i]]$resultar$State1
  curve3_s$State2[[i]] <- results3_s[[i]]$resultar$State2
  curve3_s$State3[[i]] <- results3_s[[i]]$resultar$State3
  
  hypothetical3_s$State1[[i]] <- results3_s[[i]]$resulthypothetical$State1
  hypothetical3_s$State2[[i]] <- results3_s[[i]]$resulthypothetical$State2
  hypothetical3_s$State3[[i]] <- results3_s[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:length(results3_s)) {
  lines(curve3_s$State1[[i]]$x, curve3_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve3_s$State2[[i]]$x, curve3_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve3_s$State3[[i]]$x, curve3_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)



plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:length(results3_s)) {
  lines(hypothetical3_s$State1[[i]]$x, hypothetical3_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(hypothetical3_s$State2[[i]]$x, hypothetical3_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(hypothetical3_s$State3[[i]]$x, hypothetical3_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)

### Flexible Dirichlet regression with error handling ----
num_simulations_s <- 200

gam_results3_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      oneDirGAM(n = n_s, rho=rho_s, mu = mu_s, sig = sig_s, beta = beta_s, 
                periodic = periodic_s, par=par_s, dat = dat_s, min_pred=-8, max_pred=8)
    }, error = function(e) {
      message(sprintf("Simulation %d failed: %s", i, e$message))
      NULL
    })
    
    if(!is.null(res)) {
      successful_results[[j]] <- res
      j <- j + 1
    }
    
    p(message = sprintf("Done %d/%d (%d successful)", i, num_simulations_s, j-1))
  }
  
  successful_results
})

message(sprintf("Dirichlet GAM: Completed %d out of %d simulations successfully", 
                length(gam_results3_s), num_simulations_s))

trim_to_range_s <- function(x, y, xmin = -8.05, xmax = 8.05) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}
ks1_3_s <- ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1)
cut1_3_s <- trim_to_range_s(ks1_3_s$x, ks1_3_s$y)
ks2_3_s <- ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1)
cut2_3_s <- trim_to_range_s(ks2_3_s$x, ks2_3_s$y)
ks3_3_s <- ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1)
cut3_3_s <- trim_to_range_s(ks3_3_s$x, ks3_3_s$y)


gamcurve3_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(gam_results3_s)) {
  gamcurve3_s$State1[[i]] <- gam_results3_s[[i]]$State1
  gamcurve3_s$State2[[i]] <- gam_results3_s[[i]]$State2
  gamcurve3_s$State3[[i]] <- gam_results3_s[[i]]$State3
}


plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:length(gam_results3_s)) {
  lines(gamcurve3_s$State1[[i]]$x, gamcurve3_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(gamcurve3_s$State2[[i]]$x, gamcurve3_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(gamcurve3_s$State3[[i]]$x, gamcurve3_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}

lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)
legend("right",col = c(colour_s, "transparent"),lwd = 3,bty = "n",
       lty = c(1,1,1),legend = expression(state~1, state~2, state~3))


lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
legend("right",col = c(colour_s, "transparent", "black", "black"),lwd = 3,bty = "n",
       lty = c(1,1,1, NA, 1, 3),legend = expression(state~1, state~2, state~3, "", delta, rho))



### AR Approach with error handling ----

num_covsim_s <- 200
num_simulations_s <- 200

results3AR_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneAr(n = n_s, rho=rho_s, mu = mu_s, sig = sig_s, beta = beta_s, 
               periodic = periodic_s, par=par_s, dat= dat_s, num_covsim = num_covsim_s)
    }, error = function(e) {
      message(sprintf("Simulation %d failed: %s", i, e$message))
      NULL
    })
    
    if(!is.null(res)) {
      successful_results[[j]] <- res
      j <- j + 1
    }
    
    p(message = sprintf("Done %d/%d (%d successful)", i, num_simulations_s, j-1))
  }
  
  successful_results
})

message(sprintf("AR: Completed %d out of %d simulations successfully", 
                length(results3AR_s), num_simulations_s))

curve3AR_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(results3AR_s)) {
  curve3AR_s$State1[[i]] <- results3AR_s[[i]]$State1
  curve3AR_s$State2[[i]] <- results3AR_s[[i]]$State2
  curve3AR_s$State3[[i]] <- results3AR_s[[i]]$State3
}


par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", bty = "n",
     xlim=c(min(zGT3_s)+1, max(zGT3_s)-1), main="AR resampling")

for (i in 1:length(results3AR_s)) {
  lines(curve3AR_s$State1[[i]]$x, curve3AR_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve3AR_s$State2[[i]]$x, curve3AR_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve3AR_s$State3[[i]]$x, curve3AR_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)

### BB Approach with error handling ----
num_covsim_s <- 2000
num_simulations_s <- 200

results_3BB_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneBB(n = n_s, rho = rho_s, mu = mu_s, sig = sig_s,
               beta = beta_s, periodic = periodic_s, par = par_s, dat= dat_s,
               num_covsim = num_covsim_s, blocklength = 24)
    }, error = function(e) {
      message(sprintf("Simulation %d failed: %s", i, e$message))
      NULL
    })
    
    if(!is.null(res)) {
      successful_results[[j]] <- res
      j <- j + 1
    }
    
    p(message = sprintf("Done %d/%d (%d successful)", i, num_simulations_s, j-1))
  }
  
  successful_results
})

message(sprintf("BB resampling: Completed %d out of %d simulations successfully", 
                length(results_3BB_s), num_simulations_s))

curve_3BB_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(results_3BB_s)) {
  curve_3BB_s$State1[[i]] <- results_3BB_s[[i]]$State1
  curve_3BB_s$State2[[i]] <- results_3BB_s[[i]]$State2
  curve_3BB_s$State3[[i]] <- results_3BB_s[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (III): BB resampling", bty = "n",
     xlim=c(min(zGT3_s)+1, max(zGT3_s)-1))

for (i in 1:length(results_3BB_s)) {
  lines(curve_3BB_s$State1[[i]]$x, curve_3BB_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_3BB_s$State2[[i]]$x, curve_3BB_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_3BB_s$State3[[i]]$x, curve_3BB_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3_s, mean_stateprobsGT3_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour_s), lwd = 2, bty = "n")