## Simulation experiments: Setting I (highly persistent covariate following an AR(1) process)
## Alternative n version with "_s" suffix (n=500)

# Functions and libraries ----
set.seed(22)
source("functions/sim_fit_inhomogeneousHMM.r")
source("functions/ar_approach.r")
source("functions/bb_approach.r")
source("functions/dir_reg.r")
source("functions/sim_study.r")

colour_s = c("orange", "skyblue", "seagreen")

# Parameters ----
N_s <- 3 # number of states
n_s <- 500 # dimension of simulated dataset (length of time series) - CHANGED
num_simulations_s <- 200 # number of simulated datasets 
rho_s <- 0.95 # highly persistent covariate
epsilon_s <- sqrt(1 - rho_s^2)
periodic_s <- FALSE # is the covariate periodic? 


mu_s <- c(6, 15, 20) # state-dependent means
sig_s = c(3, 1.5, 1.5) # state-dependent standard deviations

beta_s <- matrix(c( 1, 2,   # 1-> 2
                    1, -2.5,  # 1-> 3
                    -1, 0.5,   # 2-> 1
                    -1, -2.5, # 2-> 3
                    1, -1,  # 3-> 1
                    1, 1), # 3-> 2
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
Gamma_s <- tpm_g(sim_s$z, beta_s)
Delta_s <- matrix(NA, n_s, N_s)
Delta_s[1,] <- rep(1/N_s, N_s)
for (t in 2:n_s) Delta_s[t,] <- Delta_s[t-1,] %*% Gamma_s[,,t]


### GROUNDTRUTH (T=10,000,000) ----
n_z_s = 20000000
zGT1_s = numeric(n_z_s)
zGT1_s[1] = rnorm(1)

for(t in 2:n_z_s){ zGT1_s[t] = rho_s * zGT1_s[t - 1] + rnorm(1, sd = sqrt(1 - rho_s^2))}
zGT1_s=matrix(zGT1_s)

GammaGT_s = tpm_g(zGT1_s, beta_s)
delta_s = rep(1/3, 3) 

DeltaGT_s = matrix(NA, n_z_s, 3)
DeltaGT_s[1,] = rep(1/3, 3)
for(t in 2:n_z_s) DeltaGT_s[t,] = DeltaGT_s[t-1,] %*% GammaGT_s[,,t]

z_binsGT1_s <- seq(min(zGT1_s), max(zGT1_s), length.out = 100)
bin_midpointsGT1_s <- (z_binsGT1_s[-1] + z_binsGT1_s[-length(z_binsGT1_s)]) / 2

bin_ids1_s <- cut(zGT1_s, breaks = z_binsGT1_s, include.lowest = TRUE, labels = FALSE)
bin_ids1_factor_s <- factor(bin_ids1_s, levels = 1:(length(z_binsGT1_s) - 1))  

mean_stateprobsGT1_s <- matrix(NA, nrow = length(z_binsGT1_s) - 1, ncol = 3)

counts1_s <- as.vector(table(bin_ids1_factor_s))  

for (state in 1:3) {
  sums_s <- rowsum(DeltaGT_s[, state], bin_ids1_factor_s)
  mean_stateprobsGT1_s[, state] <- NA  
  nonzero_bins_s <- counts1_s > 0
  mean_stateprobsGT1_s[nonzero_bins_s, state] <- sums_s[nonzero_bins_s] / counts1_s[nonzero_bins_s]
}

mean_stateprobsGT1_s <- matrix(NA, nrow = length(z_binsGT1_s) - 1, ncol = 3)
q2.5_stateprobsGT1_s <- matrix(NA, nrow = length(z_binsGT1_s) - 1, ncol = 3)
q97.5_stateprobsGT1_s <- matrix(NA, nrow = length(z_binsGT1_s) - 1, ncol = 3)


for (state in 1:3) {
  vals_split_s <- split(DeltaGT_s[, state], bin_ids1_factor_s)
  for (b in seq_along(vals_split_s)) {
    vals_s <- vals_split_s[[b]]
    if (length(vals_s) > 0 && !all(is.na(vals_s))) {
      mean_stateprobsGT1_s[b, state]  <- mean(vals_s, na.rm = TRUE)
      q2.5_stateprobsGT1_s[b, state]  <- quantile(vals_s, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT1_s[b, state] <- quantile(vals_s, 0.975, na.rm = TRUE)
    }
  }
}

mean_stateprobsGT1_s <- na.approx(mean_stateprobsGT1_s, x = bin_midpointsGT1_s, na.rm = FALSE)
zseqGT1_s = seq(min(zGT1_s), max(zGT1_s), length = 200)
GammaseqGT1_s = tpm_g(zseqGT1_s, beta_s)
DeltaseqGT1_s = matrix(NA, length(zseqGT1_s), 3)
for(t in 1:length(zseqGT1_s)) DeltaseqGT1_s[t,] = LaMa::stationary(GammaseqGT1_s[,,t]) 

par(mfrow=c(1,1))
plot(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 1], "normal", bandwidth=1), col = colour_s[1], lwd = 3, ylab = "Pr(state | z)", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n")
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[,1], "normal", bandwidth =1), col=colour_s[1], lwd=3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
lines(ksmooth(zseqGT1_s, DeltaseqGT1_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT1_s, DeltaseqGT1_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT1_s, DeltaseqGT1_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3, lty=2)


### AR -----
# note: here we added error handling to skip simulations that fail to converge
# generally maximum of 1 or 2 out of 200 simulations fail 
system.time({
num_covsim_s <- 500
num_simulations_s <- 200

results_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneAr(n = n_s, rho=rho_s, mu = mu_s, sig = sig_s, beta = beta_s, 
               periodic = periodic_s, par=par_s, dat=dat_s, num_covsim = num_covsim_s)
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

message(sprintf("Completed %d out of %d simulations successfully", length(results_s), num_simulations_s))

curve_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(results_s)) {
  curve_s$State1[[i]] <- results_s[[i]]$State1
  curve_s$State2[[i]] <- results_s[[i]]$State2
  curve_s$State3[[i]] <- results_s[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:length(results_s)) {
  lines(curve_s$State1[[i]]$x, curve_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_s$State2[[i]]$x, curve_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_s$State3[[i]]$x, curve_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)

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
})
## Runtime on Apple M4 with 16GB RAM
# user   system   elapsed 
# 370.358 246.981 146.307 

#### AR and hypothetical -----

system.time({
num_covsim_s <- 500
num_simulations_s <- 200

# AR and hypothetical with error handling
results1_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneArHypothetical(
        n = n_s, rho = rho_s, mu = mu_s, sig = sig_s,
        beta = beta_s, periodic = periodic_s, par = par_s, dat=dat_s,
        num_covsim = num_covsim_s, min_z = -4, max_z = 4
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
                length(results1_s), num_simulations_s))



curve1_s <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical_s <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1: length(results1_s)) {
  curve1_s$State1[[i]] <- results1_s[[i]]$resultar$State1
  curve1_s$State2[[i]] <- results1_s[[i]]$resultar$State2
  curve1_s$State3[[i]] <- results1_s[[i]]$resultar$State3
  
  hypothetical_s$State1[[i]] <- results1_s[[i]]$resulthypothetical$State1
  hypothetical_s$State2[[i]] <- results1_s[[i]]$resulthypothetical$State2
  hypothetical_s$State3[[i]] <- results1_s[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "AR resampling (Setting I)", bty = "n",
     xlim = c(-4, 4))

for (i in 1: length(results1_s)) {
  lines(curve1_s$State1[[i]]$x, curve1_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve1_s$State2[[i]]$x, curve1_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve1_s$State3[[i]]$x, curve1_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "Hypothetical stationary distribution (Setting I)", bty = "n",
     xlim = c(-4, 4))

for (i in 1: length(results1_s)) {
  lines(hypothetical_s$State1[[i]]$x, hypothetical_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(hypothetical_s$State2[[i]]$x, hypothetical_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(hypothetical_s$State3[[i]]$x, hypothetical_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
})

## Runtime: Apple M4 with 16GB RAM
# user   system   elapsed 
# 381.576 252.350 147.605 


### Flexible Dirichlet regression ----

system.time({
num_simulations_s <- 200

# flexible Dirichlet regression with error handling
gam_results1_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      oneDirGAM(n = n_s, rho=rho_s, mu = mu_s, sig = sig_s, 
                beta = beta_s, periodic = periodic_s, par=par_s, 
                dat = dat_s, min_pred = -4, max_pred = 4)
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
                length(gam_results1_s), num_simulations_s))

gamcurve1_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(gam_results1_s)) {
  gamcurve1_s$State1[[i]] <- gam_results1_s[[i]]$State1
  gamcurve1_s$State2[[i]] <- gam_results1_s[[i]]$State2
  gamcurve1_s$State3[[i]] <- gam_results1_s[[i]]$State3
}

trim_to_range_s <- function(x, y, xmin = -4.05, xmax = 4.05) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}

ks1_1_s <- ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 1], "normal", bandwidth = 1)
cut1_1_s <- trim_to_range_s(ks1_1_s$x, ks1_1_s$y)
ks2_1_s <- ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 2], "normal", bandwidth = 1)
cut2_1_s <- trim_to_range_s(ks2_1_s$x, ks2_1_s$y)
ks3_1_s <- ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 3], "normal", bandwidth = 1)
cut3_1_s <- trim_to_range_s(ks3_1_s$x, ks3_1_s$y)


plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (I)", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:length(gam_results1_s)) {
  lines(gamcurve1_s$State1[[i]]$x, gamcurve1_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(gamcurve1_s$State2[[i]]$x, gamcurve1_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(gamcurve1_s$State3[[i]]$x, gamcurve1_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}

lines(cut1_1_s$x, cut1_1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_1_s$x, cut2_1_s$y, col = colour_s[2], lwd = 3)
lines(cut3_1_s$x, cut3_1_s$y, col = colour_s[3], lwd = 3)
legend("topleft",col = c(colour_s, "transparent"),lwd = 3,bty = "n",
       lty = c(1,1,1),legend = expression(state~1, state~2, state~3))
})
## Runtime: Apple M4 with 16GB RAM
# user   system   elapsed 
# 161.771   1.347 165.614 

### BB approach ----
system.time({
num_covsim_s <- 400
num_simulations_s <- 200

# BB approach with error handling
results_1BB_s <- with_progress({
  p <- progressor(steps = num_simulations_s)
  successful_results <- list()
  j <- 1
  
  for(i in 1:num_simulations_s) {
    res <- tryCatch({
      simOneBB(n = n_s, rho = rho_s, mu = mu_s, sig = sig_s,
               beta = beta_s, periodic = periodic_s, par = par_s, dat=dat_s,
               num_covsim = num_covsim_s, blocklength=40)
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
                length(results_1BB_s), num_simulations_s))

curve_1BB_s <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:length(results_1BB_s)) {
  curve_1BB_s$State1[[i]] <- results_1BB_s[[i]]$State1
  curve_1BB_s$State2[[i]] <- results_1BB_s[[i]]$State2
  curve_1BB_s$State3[[i]] <- results_1BB_s[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (I): BB resampling", bty = "n",
     xlim=c(min(zGT1_s)+1, max(zGT1_s)-1))

for (i in 1:length(results_1BB_s)) {
  lines(curve_1BB_s$State1[[i]]$x, curve_1BB_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_1BB_s$State2[[i]]$x, curve_1BB_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_1BB_s$State3[[i]]$x, curve_1BB_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 1], "normal", bandwidth = 1), col = colour_s[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 2], "normal", bandwidth = 1), col = colour_s[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1_s, mean_stateprobsGT1_s[, 3], "normal", bandwidth = 1), col = colour_s[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour_s), lwd = 2, bty = "n")
})
## Runtime on Apple M4 with 16GB RAM
# user   system   elapsed 
# 179.513 156.304 188.564 