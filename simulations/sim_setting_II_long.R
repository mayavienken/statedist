## Simulation experiments: Setting II (moderately persistent covariate following an AR(1) process)
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
n_l <- 10000 # dimension of simulated dataset (length of time series) - CHANGED TO LARGE N
num_simulations_l <- 200 # number of simulated datasets 
rho_l <- 0.7 # moderately persistent covariate 
epsilon_l <- sqrt(1 - rho_l^2)
periodic_l <- FALSE


mu_l <- c(6, 15, 20)
sig_l = c(3, 1.5, 1.5)

beta_l <- matrix(c( 1, 2,   # 1-> 2
                    1, -2.5,  # 1-> 3
                    -1, 0.5,   # 2-> 1
                    -1, -2.5, # 2-> 3
                    1, -1,  # 3-> 1
                    1, 1), # 3-> 2
                 nrow = 6, byrow = TRUE)

par_l = list(
  beta = c(rep(-1, 6), rep(0,6)),
  delta = c(0,1), 
  mu = mu_l, 
  sigma = log(sig_l)
)

sim_l <-  simCovHMM(n=n_l, rho=rho_l, mu=mu_l, sig=sig_l, beta=beta_l, periodic=periodic_l)

dat_l = list(x = sim_l$x,
             Z = matrix(sim_l$z), 
             N=3)

fit_l <- fitCovHMM(par=par_l, dat=dat_l)


Gamma_l <- tpm_g(sim_l$z, beta_l)
Delta_l <- matrix(NA, n_l, N_l)
Delta_l[1,] <- rep(1/N_l, N_l)
for (t in 2:n_l) Delta_l[t,] <- Delta_l[t-1,] %*% Gamma_l[,,t]



### GROUNDTRUTH (T=10,000,000) ----
# simulate a very long covariate series from the same AR(1) process 
# to obtain a smooth estimate of the true state occupancy distribution 
# as a function of the covariate value (solid lines in figure)

n_z_l = 10000000
zGT2_l = numeric(n_z_l)
zGT2_l[1] = rnorm(1)

for(t in 2:n_z_l){ zGT2_l[t] = rho_l * zGT2_l[t - 1] + rnorm(1, sd = sqrt(1 - rho_l^2))}
zGT2_l=matrix(zGT2_l)

# transition probability matrices (TPMs)
GammaGT_l = tpm_g(zGT2_l, beta_l)
delta_l = rep(1/3, 3) 

# empirical state probabilities at each time point
DeltaGT_l = matrix(NA, n_z_l, 3)
DeltaGT_l[1,] = rep(1/3, 3)
for(t in 2:n_z_l) DeltaGT_l[t,] = DeltaGT_l[t-1,] %*% GammaGT_l[,,t]

# split covariate range into bins and compute mean state probabilities within each bin
z_binsGT2_l <- seq(min(zGT2_l), max(zGT2_l), length.out = 100)
bin_midpointsGT2_l <- (z_binsGT2_l[-1] + z_binsGT2_l[-length(z_binsGT2_l)]) / 2

bin_ids2_l <- cut(zGT2_l, breaks = z_binsGT2_l, include.lowest = TRUE, labels = FALSE)
bin_ids2_factor_l <- factor(bin_ids2_l, levels = 1:(length(z_binsGT2_l) - 1))  

mean_stateprobsGT2_l <- matrix(NA, nrow = length(z_binsGT2_l) - 1, ncol = 3)

counts2_l <- as.vector(table(bin_ids2_factor_l))  

for (state in 1:3) {
  sums_l <- rowsum(DeltaGT_l[, state], bin_ids2_factor_l)
  mean_stateprobsGT2_l[, state] <- NA  
  nonzero_bins_l <- counts2_l > 0
  mean_stateprobsGT2_l[nonzero_bins_l, state] <- sums_l[nonzero_bins_l] / counts2_l[nonzero_bins_l]
}

mean_stateprobsGT2_l <- matrix(NA, nrow = length(z_binsGT2_l) - 1, ncol = 3)
q2.5_stateprobsGT2_l <- matrix(NA, nrow = length(z_binsGT2_l) - 1, ncol = 3)
q97.5_stateprobsGT2_l <- matrix(NA, nrow = length(z_binsGT2_l) - 1, ncol = 3)

for (state in 1:3) {
  vals_split_l <- split(DeltaGT_l[, state], bin_ids2_factor_l)
  for (b in seq_along(vals_split_l)) {
    vals_l <- vals_split_l[[b]]
    if (length(vals_l) > 0 && !all(is.na(vals_l))) {
      mean_stateprobsGT2_l[b, state]  <- mean(vals_l, na.rm = TRUE)
      q2.5_stateprobsGT2_l[b, state]  <- quantile(vals_l, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT2_l[b, state] <- quantile(vals_l, 0.975, na.rm = TRUE)
    }
  }
}

mean_stateprobsGT2_l <- na.approx(mean_stateprobsGT2_l, x = bin_midpointsGT2_l, na.rm = FALSE)
zseqGT2_l = seq(min(zGT2_l), max(zGT2_l), length = 200)
GammaseqGT2_l = tpm_g(zseqGT2_l, beta_l)
DeltaseqGT2_l = matrix(NA, length(zseqGT2_l), 3)
for(t in 1:length(zseqGT2_l)) DeltaseqGT2_l[t,] = LaMa::stationary(GammaseqGT2_l[,,t]) 

# hypothetical stationary distribution (dashed lines) vs. true state occupancy distribution (solid lines)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 1], "normal", bandwidth=1), col = colour_l[1], lwd = 3, ylab = "Pr(state | z)", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n", xlim=c(-5,5.5))
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
lines(ksmooth(zseqGT2_l, DeltaseqGT2_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_l, DeltaseqGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_l, DeltaseqGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3, lty=2)
legend("right",col = c(colour_l, "transparent", "black", "black"),
       lwd = 3,bty = "n",lty = c(1,1,1, NA, 1, 3),legend = expression(
         state~1, state~2, state~3, "", delta, rho))


### AR -----

num_covsim_l <- 200
num_simulations_l <- 200

results_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- simOneAr(n = n_l, rho=rho_l, mu = mu_l, sig = sig_l, beta = beta_l, periodic = periodic_l, par=par_l, num_covsim = num_covsim_l, dat=dat_l)
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
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
lines(ksmooth(zseqGT2_l, DeltaseqGT2_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_l, DeltaseqGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT2_l, DeltaseqGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3, lty=2)


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

#### AR and Hypothetical -----
num_covsim_l <- 200
num_simulations_l <- 200

results2_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- simOneArHypothetical(
      n = n_l, rho = rho_l, mu = mu_l, sig = sig_l,
      beta = beta_l, periodic = periodic_l, par = par_l, dat=dat_l,
      num_covsim = num_covsim_l
    )
    p(message = sprintf("Done %d/%d", i, num_simulations_l))
    res
  })
})


curve2_l <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical2_l <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:num_simulations_l) {
  curve2_l$State1[[i]] <- results2_l[[i]]$resultar$State1
  curve2_l$State2[[i]] <- results2_l[[i]]$resultar$State2
  curve2_l$State3[[i]] <- results2_l[[i]]$resultar$State3
  
  hypothetical2_l$State1[[i]] <- results2_l[[i]]$resulthypothetical$State1
  hypothetical2_l$State2[[i]] <- results2_l[[i]]$resulthypothetical$State2
  hypothetical2_l$State3[[i]] <- results2_l[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main="Setting (II): AR resampling", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(curve2_l$State1[[i]]$x, curve2_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve2_l$State2[[i]]$x, curve2_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve2_l$State3[[i]]$x, curve2_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
legend("topleft",col = colour_l, lwd = 3,bty = "n",legend = expression(state~1, state~2, state~3), cex=1.3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "Setting (II): hypothetical stationary distribution", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(hypothetical2_l$State1[[i]]$x, hypothetical2_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(hypothetical2_l$State2[[i]]$x, hypothetical2_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(hypothetical2_l$State3[[i]]$x, hypothetical2_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)



### Flexible Dirichlet regression ----

gam_results2_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- oneDirGAM(n = n_l, rho = rho_l, mu = mu_l, sig = sig_l, beta = beta_l, periodic = periodic_l, 
                     par = par_l, dat = dat_l, min_pred = -4, max_pred = 4)
    p()
    res
  })
})

gamcurve2_l <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations_l) {
  gamcurve2_l$State1[[i]] <- gam_results2_l[[i]]$State1
  gamcurve2_l$State2[[i]] <- gam_results2_l[[i]]$State2
  gamcurve2_l$State3[[i]] <- gam_results2_l[[i]]$State3
}

trim_to_range_l <- function(x, y, xmin = -4.1, xmax = 4.1) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}

ks1_2_l <- ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 1], "normal", bandwidth = 1)
cut1_2_l <- trim_to_range_l(ks1_2_l$x, ks1_2_l$y)
ks2_2_l <- ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 2], "normal", bandwidth = 1)
cut2_2_l <- trim_to_range_l(ks2_2_l$x, ks2_2_l$y)
ks3_2_l <- ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 3], "normal", bandwidth = 1)
cut3_2_l <- trim_to_range_l(ks3_2_l$x, ks3_2_l$y)

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(gamcurve2_l$State1[[i]]$x, gamcurve2_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(gamcurve2_l$State2[[i]]$x, gamcurve2_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(gamcurve2_l$State3[[i]]$x, gamcurve2_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}

lines(cut1_2_l$x, cut1_2_l$y, col = colour_l[1], lwd = 3)
lines(cut2_2_l$x, cut2_2_l$y, col = colour_l[2], lwd = 3)
lines(cut3_2_l$x, cut3_2_l$y, col = colour_l[3], lwd = 3)



### BB Approach ----
num_covsim_l <- 2000
num_simulations_l <- 200

results_2BB_l <- with_progress({
  p <- progressor(steps = num_simulations_l)
  lapply(1:num_simulations_l, function(i) {
    res <- simOneBB(n = n_l, rho = rho_l, mu = mu_l, sig = sig_l,
                    beta = beta_l, periodic = periodic_l, par = par_l, dat=dat_l,
                    num_covsim = num_covsim_l, blocklength=20)
    p(message = sprintf("Done %d/%d", i, num_simulations_l))
    res
  })
})

curve_2BB_l <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations_l) {
  curve_2BB_l$State1[[i]] <- results_2BB_l[[i]]$State1
  curve_2BB_l$State2[[i]] <- results_2BB_l[[i]]$State2
  curve_2BB_l$State3[[i]] <- results_2BB_l[[i]]$State3
}

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "Setting (II): BB resampling", bty = "n",
     xlim=c(min(zGT2_l)+1, max(zGT2_l)-1))

for (i in 1:num_simulations_l) {
  lines(curve_2BB_l$State1[[i]]$x, curve_2BB_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve_2BB_l$State2[[i]]$x, curve_2BB_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve_2BB_l$State3[[i]]$x, curve_2BB_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 1], "normal", bandwidth = 1), col = colour_l[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 2], "normal", bandwidth = 1), col = colour_l[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2_l, mean_stateprobsGT2_l[, 3], "normal", bandwidth = 1), col = colour_l[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour_l), lwd = 2, bty = "n")