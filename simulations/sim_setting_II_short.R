## Simulation experiments: Setting II (moderately persistent covariate following an AR(1) process)
## Alternative n version with "_s" suffix and error handling

# Functions and libraries ----

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
rho_s <- 0.7 # moderately persistent covariate 
epsilon_s <- sqrt(1 - rho_s^2)
periodic_s <- FALSE


mu_s <- c(6, 15, 20)
sig_s = c(3, 1.5, 1.5)

beta_s <- matrix(c( 1, 2,   # 1-> 2
                    1, -2.5,  # 1-> 3
                    -1, 0.5,   # 2-> 1
                    -1, -2.5, # 2-> 3
                    1, -1,  # 3-> 1
                    1, 1), # 3-> 2
                 nrow = 6, byrow = TRUE)

par_s = c(beta = c(rep(-1, 6), rep(0,6)),
          logitdelta = c(0,1), 
          mu = mu_s, 
          sig = sig_s)

sim_s <-  simCovHMM(n=n_s, rho=rho_s, mu=mu_s, sig=sig_s, beta=beta_s, periodic=periodic_s)
par_s = c(beta = c(rep(-1, 6), rep(0,6)),
          logitdelta = c(0,1), 
          mu = c(6, 15, 20), 
          sig = c(log(3),log(3),log(3)))
fit_s <- fitCovHMM(par=par_s, x=sim_s$x, Z=matrix(sim_s$z))


Gamma_s <- tpm_g(sim_s$z, beta_s)
Delta_s <- matrix(NA, n_s, N_s)
Delta_s[1,] <- rep(1/N_s, N_s)
for (t in 2:n_s) Delta_s[t,] <- Delta_s[t-1,] %*% Gamma_s[,,t]



### GROUNDTRUTH (T=10,000,000) ----
# simulate a very long covariate series from the same AR(1) process 
# to obtain a smooth estimate of the true state occupancy distribution 
# as a function of the covariate value (solid lines in figure)

n_z_s = 10000000
zGT2_s = numeric(n_z_s)
zGT2_s[1] = rnorm(1)

for(t in 2:n_z_s){ zGT2_s[t] = rho_s * zGT2_s[t - 1] + rnorm(1, sd = sqrt(1 - rho_s^2))}
zGT2_s=matrix(zGT2_s)

# transition probability matrices (TPMs)
GammaGT_s = tpm_g(zGT2_s, beta_s)
delta_s = rep(1/3, 3) 

# empirical state probabilities at each time point
DeltaGT_s = matrix(NA, n_z_s, 3)
DeltaGT_s[1,] = rep(1/3, 3)
for(t in 2:n_z_s) DeltaGT_s[t,] = DeltaGT_s[t-1,] %*% GammaGT_s[,,t]

# split covariate range into bins and compute mean state probabilities within each bin
z_binsGT2_s <- seq(min(zGT2_s), max(zGT2_s), length.out = 100)
bin_midpointsGT2_s <- (z_binsGT2_s[-1] + z_binsGT2_s[-length(z_binsGT2_s)]) / 2

bin_ids2_s <- cut(zGT2_s, breaks = z_binsGT2_s, include.lowest = TRUE, labels = FALSE)
bin_ids2_factor_s <- factor(bin_ids2_s, levels = 1:(length(z_binsGT2_s) - 1))  

mean_stateprobsGT2_s <- matrix(NA, nrow = length(z_binsGT2_s) - 1, ncol = 3)

counts2_s <- as.vector(table(bin_ids2_factor_s))  

for (state in 1:3) {
  sums_s <- rowsum(DeltaGT_s[, state], bin_ids2_factor_s)
  mean_stateprobsGT2_s[, state] <- NA  
  nonzero_bins_s <- counts2_s > 0
  mean_stateprobsGT2_s[nonzero_bins_s, state] <- sums_s[nonzero_bins_s] / counts2_s[nonzero_bins_s]
}

mean_stateprobsGT2_s <- matrix(NA, nrow = length(z_binsGT2_s) - 1, ncol = 3)
q2.5_stateprobsGT2_s <- matrix(NA, nrow = length(z_binsGT2_s) - 1, ncol = 3)
q97.5_stateprobsGT2_s <- matrix(NA, nrow = length(z_binsGT2_s) - 1, ncol = 3)

for (state in 1:3) {
  vals_split_s <- split(DeltaGT_s[, state], bin_ids2_factor_s)
  for (b in seq_along(vals_split_s)) {
    vals_s <- vals_split_s[[b]]
    if (length(vals_s) > 0 && !all(is.na(vals_s))) {
      mean_stateprobsGT2_s[b, state]  <- mean(vals_s, na.rm = TRUE)
      q2.5_stateprobsGT2_s[b, state]  <- quantile(vals_s, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT2_s[b, state] <- quantile(vals_s, 0.975, na.rm = TRUE)
    }
  }
}

mean_stateprobsGT2_s <- na.approx(mean_stateprobsGT2_s, x = bin_midpointsGT2_s, na.rm = FALSE)
zseqGT2_s = seq(min(zGT2_s), max(zGT2_s), length = 200)
GammaseqGT2_s = tpm_g(zseqGT2_s, beta_s)
DeltaseqGT2_s = matrix(NA, length(zseqGT2_s), 3)
for(t in 1:length(zseqGT2_s)) DeltaseqGT2_s[t,] = LaMa::stationary(GammaseqGT2_s[,,t]) 

# hypothetical stationary distribution (dashed lines) vs. true state occupancy distribution (solid lines)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 1], "normal", bandwidth=1), col = colour_s[1], lwd = 3, ylab = "Pr(state i | z), i=1,2,3", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n", xlim=c(-5,5.5))
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
lines(ksmooth(zseqGT2_s, DeltaseqGT2_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_s, DeltaseqGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_s, DeltaseqGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3, lty=2)
legend("right",col = c(colour_s, "transparent", "black", "black"),
       lwd = 3,bty = "n",lty = c(1,1,1, NA, 1, 3),legend = expression(
         state~1, state~2, state~3, "", delta, rho))


### AR with error handling -----

num_covsim_s <- 200
num_simulations_s <- 200

results_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneAr(n = n_s, rho=rho_s, mu = mu_s, sig = sig_s, beta = beta_s, 
               periodic = periodic_s, par=par_s, num_covsim = num_covsim_s)
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
                length(results_s), num_simulations_s))

curve_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(results_s)) {
  curve_s$State1[[i]] <- results_s[[i]]$State1
  curve_s$State2[[i]] <- results_s[[i]]$State2
  curve_s$State3[[i]] <- results_s[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i|z), i=1,2,3", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:length(results_s)) {
  lines(curve_s$State1[[i]]$x, curve_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_s$State2[[i]]$x, curve_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_s$State3[[i]]$x, curve_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
lines(ksmooth(zseqGT2_s, DeltaseqGT2_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_s, DeltaseqGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_s, DeltaseqGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3, lty=2)


legend("right",
       col = c(colour_s, "transparent", "black", "black"),
       lwd = 3,
       bty = "n",
       lty = c(1,1,1, NA, 1, 3),
       legend = expression(
         state~1, state~2, state~3, 
         "", 
         delta, 
         rho
       ))

#### AR and Hypothetical with error handling -----
num_covsim_s <- 200
num_simulations_s <- 200

results2_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneArHypothetical(
        n = n_s, rho = rho_s, mu = mu_s, sig = sig_s,
        beta = beta_s, periodic = periodic_s, par = par_s,
        num_covsim = num_covsim_s
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
                length(results2_s), num_simulations_s))


curve2_s <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical2_s <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:length(results2_s)) {
  curve2_s$State1[[i]] <- results2_s[[i]]$resultar$State1
  curve2_s$State2[[i]] <- results2_s[[i]]$resultar$State2
  curve2_s$State3[[i]] <- results2_s[[i]]$resultar$State3
  
  hypothetical2_s$State1[[i]] <- results2_s[[i]]$resulthypothetical$State1
  hypothetical2_s$State2[[i]] <- results2_s[[i]]$resulthypothetical$State2
  hypothetical2_s$State3[[i]] <- results2_s[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main="Setting (II): AR resampling", bty = "n",
     xlim = c(-4, 4))

for (i in 1:length(results2_s)) {
  lines(curve2_s$State1[[i]]$x, curve2_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve2_s$State2[[i]]$x, curve2_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve2_s$State3[[i]]$x, curve2_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
legend("topleft",col = colour_s, lwd = 3,bty = "n",legend = expression(state~1, state~2, state~3), cex=1.3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "Setting (II): hypothetical stationary distribution", bty = "n",
     xlim = c(-4, 4))

for (i in 1:length(results2_s)) {
  lines(hypothetical2_s$State1[[i]]$x, hypothetical2_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(hypothetical2_s$State2[[i]]$x, hypothetical2_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(hypothetical2_s$State3[[i]]$x, hypothetical2_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)



### Flexible Dirichlet regression with error handling ----

gam_results2_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      oneDirGAM(n = n_s, rho=rho_s, mu = mu_s, sig = sig_s, beta = beta_s, 
                periodic = periodic_s, par=par_s)
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
                length(gam_results2_s), num_simulations_s))

gamcurve2_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(gam_results2_s)) {
  gamcurve2_s$State1[[i]] <- gam_results2_s[[i]]$State1
  gamcurve2_s$State2[[i]] <- gam_results2_s[[i]]$State2
  gamcurve2_s$State3[[i]] <- gam_results2_s[[i]]$State3
}

trim_to_range_s <- function(x, y, xmin=min(zGT2_s), xmax=max(zGT2_s)) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}

ks1_s <- ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 1], "normal", bandwidth = 1)
cut1_s <- trim_to_range_s(ks1_s$x, ks1_s$y)
ks2_s <- ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 2], "normal", bandwidth = 1)
cut2_s <- trim_to_range_s(ks2_s$x, ks2_s$y)
ks3_s <- ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 3], "normal", bandwidth = 1)
cut3_s <- trim_to_range_s(ks3_s$x, ks3_s$y)

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i|z), i=1,2,3", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:length(gam_results2_s)) {
  lines(gamcurve2_s$State1[[i]]$x, gamcurve2_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(gamcurve2_s$State2[[i]]$x, gamcurve2_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(gamcurve2_s$State3[[i]]$x, gamcurve2_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}

lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)



### BB Approach with error handling ----
num_covsim_s <- 2000
num_simulations_s <- 200

results_2BB_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneBB(n = n_s, rho = rho_s, mu = mu_s, sig = sig_s,
               beta = beta_s, periodic = periodic_s, par = par_s,
               num_covsim = num_covsim_s, blocklength=20)
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
                length(results_2BB_s), num_simulations_s))

curve_2BB_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(results_2BB_s)) {
  curve_2BB_s$State1[[i]] <- results_2BB_s[[i]]$State1
  curve_2BB_s$State2[[i]] <- results_2BB_s[[i]]$State2
  curve_2BB_s$State3[[i]] <- results_2BB_s[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (II): BB resampling", bty = "n",
     xlim=c(min(zGT2_s)+1, max(zGT2_s)-1))

for (i in 1:length(results_2BB_s)) {
  lines(curve_2BB_s$State1[[i]]$x, curve_2BB_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_2BB_s$State2[[i]]$x, curve_2BB_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_2BB_s$State3[[i]]$x, curve_2BB_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_s, mean_stateprobsGT2_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour_s), lwd = 2, bty = "n")