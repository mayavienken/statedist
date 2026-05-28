## Simulation experiments: Setting I (highly persistent covariate following an AR(1) process)
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
N_l <- 3 # number of states
n_l <- 10000 # dimension of simulated dataset (length of time series) - CHANGED TO LARGE N
num_simulations_l <- 200 # number of simulated datasets 
rho_l <- 0.95 # highly persistent covariate
epsilon_l <- sqrt(1 - rho_l^2)
periodic_l <- FALSE # is the covariate periodic? 


mu_l <- c(6, 15, 20) # state-dependent means
sig_l = c(3, 1.5, 1.5) # state-dependent standard deviations

beta_l <- matrix(c( 1, 2,   # 1-> 2
                    1, -2.5,  # 1-> 3
                    -1, 0.5,   # 2-> 1
                    -1, -2.5, # 2-> 3
                    1, -1,  # 3-> 1
                    1, 1), # 3-> 2
                 nrow = 6, byrow = TRUE)

# Simulation
sim_l <-  simCovHMM(n=n_l, rho=rho_l, mu=mu_l, sig=sig_l, beta=beta_l, periodic=periodic_l)

par_l = list(
  beta = c(rep(-1, 6), rep(0,6)),
  delta = c(0,1), 
  mu = mu_l, 
  sigma = log(sig_l)
)

dat_l = list(x = sim$x,
           Z = matrix(sim$z), 
           N=3)

fit <- fitCovHMM(par=par_l, dat=dat_l)

Gamma_l <- tpm_g(sim_l$z, beta_l)
Delta_l <- matrix(NA, n_l, N_l)
Delta_l[1,] <- rep(1/N_l, N_l)
for (t in 2:n_l) Delta_l[t,] <- Delta_l[t-1,] %*% Gamma_l[,,t]


### GROUNDTRUTH (T=10,000,000) ----
n_z_l = 10000000
zGT1_l = numeric(n_z_l)
zGT1_l[1] = rnorm(1)

for(t in 2:n_z_l){ zGT1_l[t] = rho_l * zGT1_l[t - 1] + rnorm(1, sd = sqrt(1 - rho_l^2))}
zGT1_l=matrix(zGT1_l)

GammaGT_l = tpm_g(zGT1_l, beta_l)
delta_l = rep(1/3, 3) 

DeltaGT_l = matrix(NA, n_z_l, 3)
DeltaGT_l[1,] = rep(1/3, 3)
for(t in 2:n_z_l) DeltaGT_l[t,] = DeltaGT_l[t-1,] %*% GammaGT_l[,,t]

z_binsGT1_l <- seq(min(zGT1_l), max(zGT1_l), length.out = 100)
bin_midpointsGT1_l <- (z_binsGT1_l[-1] + z_binsGT1_l[-length(z_binsGT1_l)]) / 2

bin_ids1_l <- cut(zGT1_l, breaks = z_binsGT1_l, include.lowest = TRUE, labels = FALSE)
bin_ids1_factor_l <- factor(bin_ids1_l, levels = 1:(length(z_binsGT1_l) - 1))  

mean_stateprobsGT1_l <- matrix(NA, nrow = length(z_binsGT1_l) - 1, ncol = 3)

counts1_l <- as.vector(table(bin_ids1_factor_l))  

for (state in 1:3) {
  sums_l <- rowsum(DeltaGT_l[, state], bin_ids1_factor_l)
  mean_stateprobsGT1_l[, state] <- NA  
  nonzero_bins_l <- counts1_l > 0
  mean_stateprobsGT1_l[nonzero_bins_l, state] <- sums_l[nonzero_bins_l] / counts1_l[nonzero_bins_l]
}

mean_stateprobsGT1_l <- matrix(NA, nrow = length(z_binsGT1_l) - 1, ncol = 3)
q2.5_stateprobsGT1_l <- matrix(NA, nrow = length(z_binsGT1_l) - 1, ncol = 3)
q97.5_stateprobsGT1_l <- matrix(NA, nrow = length(z_binsGT1_l) - 1, ncol = 3)


for (state in 1:3) {
  vals_split_l <- split(DeltaGT_l[, state], bin_ids1_factor_l)
  for (b in seq_along(vals_split_l)) {
    vals_l <- vals_split_l[[b]]
    if (length(vals_l) > 0 && !all(is.na(vals_l))) {
      mean_stateprobsGT1_l[b, state]  <- mean(vals_l, na.rm = TRUE)
      q2.5_stateprobsGT1_l[b, state]  <- quantile(vals_l, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT1_l[b, state] <- quantile(vals_l, 0.975, na.rm = TRUE)
    }
  }
}

mean_stateprobsGT1_l <- na.approx(mean_stateprobsGT1_l, x = bin_midpointsGT1_l, na.rm = FALSE)
zseqGT1_l = seq(min(zGT1_l), max(zGT1_l), length = 200)
GammaseqGT1_l = tpm_g(zseqGT1_l, beta_l)
DeltaseqGT1_l = matrix(NA, length(zseqGT1_l), 3)
for(t in 1:length(zseqGT1_l)) DeltaseqGT1_l[t,] = LaMa::stationary(GammaseqGT1_l[,,t]) 

par(mfrow=c(1,1))
plot(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 1], "normal", bandwidth=1), col = colour_l[1], lwd = 3, ylab = "Pr(state | z)", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n")
lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[,1], "normal", bandwidth =1), col=colour_l[1], lwd=3)
lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
lines(ksmooth(zseqGT1_l, DeltaseqGT1_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT1_l, DeltaseqGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT1_l, DeltaseqGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3, lty=2)


### AR -----

system.time({
  num_covsim_l <- 200
  num_simulations_l <- 200
  
  results_l <- with_progress({
    p <- progressor(steps = num_simulations_l)
    lapply(1:num_simulations_l, function(i) {
      res <- simOneAr(n = n_l, rho=rho_l, mu = mu_l, sig = sig_l, beta = beta_l, 
                      periodic = periodic_l, par=par_l, dat=dat_l, num_covsim = num_covsim_l)
      p(message = sprintf("Done %d/%d", i, num_simulations_l))
      res
    })
  })
  
  curve_l <- list(State1 = list(), State2 = list(), State3 = list())
  for (i in 1:num_simulations_l) {
    curve_l$State1[[i]] <- results_l[[i]]$State1
    curve_l$State2[[i]] <- results_l[[i]]$State2
    curve_l$State3[[i]] <- results_l[[i]]$State3
  }
  
  plot(NULL, ylim = c(0, 1),
       xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n",
       xlim=c(-4, 4))
  
  for (i in 1:num_simulations_l) {
    lines(curve_l$State1[[i]]$x, curve_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
    lines(curve_l$State2[[i]]$x, curve_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
    lines(curve_l$State3[[i]]$x, curve_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
  }
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
  lines(ksmooth(zseqGT1_l, DeltaseqGT1_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3, lty=2)
  lines(ksmooth(zseqGT1_l, DeltaseqGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3, lty=2)
  lines(ksmooth(zseqGT1_l, DeltaseqGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3, lty=2)
  
  legend("right",
         col = c(colour_l, "transparent", "black", "black"),
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

#### AR and hypothetical -----
system.time({
  num_covsim_l <- 200
  num_simulations_l <- 200
  
  results1_l <- with_progress({
    p <- progressor(steps = num_simulations_l)
    lapply(1:num_simulations_l, function(i) {
      res <- simOneArHypothetical(
        n = n_l, rho = rho_l, mu = mu_l, sig = sig_l,
        beta = beta_l, periodic = periodic_l, par = par_l, dat = dat_l,
        num_covsim = num_covsim_l, min_z = -4, max_z = 4
      )
      p(message = sprintf("Done %d/%d", i, num_simulations_l))
      res
    })
  })
  
  
  curve1_l <- list(State1 = list(), State2 = list(), State3 = list())
  hypothetical_l <- list(State1 = list(), State2 = list(), State3 = list())
  
  for (i in 1:num_simulations_l) {
    curve1_l$State1[[i]] <- results1_l[[i]]$resultar$State1
    curve1_l$State2[[i]] <- results1_l[[i]]$resultar$State2
    curve1_l$State3[[i]] <- results1_l[[i]]$resultar$State3
    
    hypothetical_l$State1[[i]] <- results1_l[[i]]$resulthypothetical$State1
    hypothetical_l$State2[[i]] <- results1_l[[i]]$resulthypothetical$State2
    hypothetical_l$State3[[i]] <- results1_l[[i]]$resulthypothetical$State3
  }
  
  par(mfrow=c(1,2))
  
  plot(NULL, ylim = c(0, 1),
       xlab = "covariate value z", ylab = "Pr(state | z)",
       main = "AR resampling (Setting I)", bty = "n",
       xlim = c(-4, 4))
  
  for (i in 1:num_simulations_l) {
    lines(curve1_l$State1[[i]]$x, curve1_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
    lines(curve1_l$State2[[i]]$x, curve1_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
    lines(curve1_l$State3[[i]]$x, curve1_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
  }
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
  
  plot(NULL, ylim = c(0, 1),
       xlab = "covariate value z", ylab = "Pr(state | z)",
       main = "Hypothetical stationary distribution (Setting I)", bty = "n",
       xlim = c(-4, 4))
  
  for (i in 1:num_simulations_l) {
    lines(hypothetical_l$State1[[i]]$x, hypothetical_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
    lines(hypothetical_l$State2[[i]]$x, hypothetical_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
    lines(hypothetical_l$State3[[i]]$x, hypothetical_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
  }
  
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
})

### Flexible Dirichlet regression ----
system.time({
  par(mfrow=c(1,1))  
  num_simulations_l <- 200
  
  gam_results1_l <- with_progress({
    p <- progressor(steps = num_simulations_l)
    lapply(1:num_simulations_l, function(i) {
      res <- oneDirGAM(n = n_l, rho=rho_l, mu = mu_l, sig = sig_l, beta = beta_l, 
                       periodic = periodic_l, par=par_l, dat = dat_l, min_pred = -4.05, max_pred = 4.05)
      p()
      res
    })
  })
  
  gamcurve1_l <- list(State1 = list(), State2 = list(), State3 = list())
  for (i in 1:num_simulations_l) {
    gamcurve1_l$State1[[i]] <- gam_results1_l[[i]]$State1
    gamcurve1_l$State2[[i]] <- gam_results1_l[[i]]$State2
    gamcurve1_l$State3[[i]] <- gam_results1_l[[i]]$State3
  }
  
  trim_to_range_l <- function(x, y, xmin = -4.05, xmax = 4.05) {
    keep <- x >= xmin & x <= xmax
    list(x = x[keep], y = y[keep])
  }
  
  ks1_1_l <- ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 1], "normal", bandwidth = 1)
  cut1_1_l <- trim_to_range_l(ks1_1_l$x, ks1_1_l$y)
  ks2_1_l <- ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 2], "normal", bandwidth = 1)
  cut2_1_l <- trim_to_range_l(ks2_1_l$x, ks2_1_l$y)
  ks3_1_l <- ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 3], "normal", bandwidth = 1)
  cut3_1_l <- trim_to_range_l(ks3_1_l$x, ks3_1_l$y)
  
  
  plot(NULL, ylim = c(0, 1),
       xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (I)", bty = "n", 
       xlim=c(-4, 4))
  
  for (i in 1:num_simulations_l) {
    lines(gamcurve1_l$State1[[i]]$x, gamcurve1_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
    lines(gamcurve1_l$State2[[i]]$x, gamcurve1_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
    lines(gamcurve1_l$State3[[i]]$x, gamcurve1_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
  }
  
  lines(cut1_1_l$x, cut1_1_l$y, col = colour_l[1], lwd = 3)
  lines(cut2_1_l$x, cut2_1_l$y, col = colour_l[2], lwd = 3)
  lines(cut3_1_l$x, cut3_1_l$y, col = colour_l[3], lwd = 3)
  legend("topleft",col = c(colour_l, "transparent"),lwd = 3,bty = "n",
         lty = c(1,1,1),legend = expression(state~1, state~2, state~3))
})


### BB approach ----
system.time({
  num_covsim_l <- 200
  num_simulations_l <- 200
  
  results_1BB_l <- with_progress({
    p <- progressor(steps = num_simulations_l)
    lapply(1:num_simulations_l, function(i) {
      res <- simOneBB(n = n_l, rho = rho_l, mu = mu_l, sig = sig_l,
                      beta = beta_l, periodic = periodic_l, par = par_l, dat = dat_l,
                      num_covsim = num_covsim_l, blocklength=60)
      p(message = sprintf("Done %d/%d", i, num_simulations_l))
      res
    })
  })
  
  curve_1BB_l <- list(State1 = list(), State2 = list(), State3 = list())
  for (i in 1:num_simulations_l) {
    curve_1BB_l$State1[[i]] <- results_1BB_l[[i]]$State1
    curve_1BB_l$State2[[i]] <- results_1BB_l[[i]]$State2
    curve_1BB_l$State3[[i]] <- results_1BB_l[[i]]$State3
  }
  
  plot(NULL, ylim = c(0, 1),
       xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (I): BB resampling", bty = "n",
       xlim=c(min(zGT1_l)+1, max(zGT1_l)-1))
  
  for (i in 1:num_simulations_l) {
    lines(curve_1BB_l$State1[[i]]$x, curve_1BB_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
    lines(curve_1BB_l$State2[[i]]$x, curve_1BB_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
    lines(curve_1BB_l$State3[[i]]$x, curve_1BB_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
  }
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
  lines(ksmooth(bin_midpointsGT1_l, mean_stateprobsGT1_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
  legend("right", legend = c("state 1", "state 2", "state 3"),
         col = c(colour_l), lwd = 2, bty = "n")
})
