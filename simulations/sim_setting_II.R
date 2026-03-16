## Simulation experiments: Setting II

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
rho <- 0.7 # moderately persistent covariate 
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
n_z = 20000000
zGT1 = numeric(n_z)
zGT1[1] = rnorm(1)

for(t in 2:n_z){ zGT1[t] = rho * zGT1[t - 1] + rnorm(1, sd = sqrt(1 - rho^2))}
zGT1=matrix(zGT1)

GammaGT = tpm_g(zGT1, beta)
delta = rep(1/3, 3) 

DeltaGT = matrix(NA, n_z, 3)
DeltaGT[1,] = rep(1/3, 3)
for(t in 2:n_z) DeltaGT[t,] = DeltaGT[t-1,] %*% GammaGT[,,t]

z_binsGT1 <- seq(min(zGT1), max(zGT1), length.out = 100)
bin_midpointsGT1 <- (z_binsGT1[-1] + z_binsGT1[-length(z_binsGT1)]) / 2

bin_ids1 <- cut(zGT1, breaks = z_binsGT1, include.lowest = TRUE, labels = FALSE)
bin_ids1_factor <- factor(bin_ids1, levels = 1:(length(z_binsGT1) - 1))  

mean_stateprobsGT1 <- matrix(NA, nrow = length(z_binsGT1) - 1, ncol = 3)

counts1 <- as.vector(table(bin_ids1_factor))  

for (state in 1:3) {
  sums <- rowsum(DeltaGT[, state], bin_ids1_factor)
  mean_stateprobsGT1[, state] <- NA  
  nonzero_bins <- counts1 > 0
  mean_stateprobsGT1[nonzero_bins, state] <- sums[nonzero_bins] / counts1[nonzero_bins]
}

mean_stateprobsGT1 <- matrix(NA, nrow = length(z_binsGT1) - 1, ncol = 3)
q2.5_stateprobsGT1 <- matrix(NA, nrow = length(z_binsGT1) - 1, ncol = 3)
q97.5_stateprobsGT1 <- matrix(NA, nrow = length(z_binsGT1) - 1, ncol = 3)


for (state in 1:3) {
  vals_split <- split(DeltaGT[, state], bin_ids1_factor)
  for (b in seq_along(vals_split)) {
    vals <- vals_split[[b]]
    if (length(vals) > 0 && !all(is.na(vals))) {
      mean_stateprobsGT1[b, state]  <- mean(vals, na.rm = TRUE)
      q2.5_stateprobsGT1[b, state]  <- quantile(vals, 0.025, na.rm = TRUE)
      q97.5_stateprobsGT1[b, state] <- quantile(vals, 0.975, na.rm = TRUE)
    }
  }
}

mean_stateprobsGT1 <- na.approx(mean_stateprobsGT1, x = bin_midpointsGT1, na.rm = FALSE)
zseqGT1 = seq(min(zGT1), max(zGT1), length = 200)
GammaseqGT1 = tpm_g(zseqGT1, beta)
DeltaseqGT1 = matrix(NA, length(zseqGT1), 3)
for(t in 1:length(zseqGT1)) DeltaseqGT1[t,] = LaMa::stationary(GammaseqGT1[,,t]) 

par(mfrow=c(1,1))
plot(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth=1), col = colour[1], lwd = 3, ylab = "Pr(state i|z), i=1,2,3", 
     xlab = "covariate value z", main = "", ylim=c(0,1), type="l", bty="n")
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[,1], "normal", bandwidth =1), col=colour[1], lwd=3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
lines(ksmooth(zseqGT1, DeltaseqGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3, lty=2)
lines(ksmooth(zseqGT1, DeltaseqGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3, lty=2)
lines(ksmooth(zseqGT1, DeltaseqGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3, lty=2)


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

#### AR and Hypothetical -----
num_covsim <- 200
num_simulations <- 200

results1 <- with_progress({
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


curve1 <- list(State1 = list(), State2 = list(), State3 = list())
hypothetical <- list(State1 = list(), State2 = list(), State3 = list())

for (i in 1:num_simulations) {
  curve1$State1[[i]] <- results1[[i]]$resultar$State1
  curve1$State2[[i]] <- results1[[i]]$resultar$State2
  curve1$State3[[i]] <- results1[[i]]$resultar$State3
  
  hypothetical$State1[[i]] <- results1[[i]]$resulthypothetical$State1
  hypothetical$State2[[i]] <- results1[[i]]$resulthypothetical$State2
  hypothetical$State3[[i]] <- results1[[i]]$resulthypothetical$State3
}

par(mfrow=c(1,2))

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "AR resampling (Setting II)", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve1$State1[[i]]$x, curve1$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve1$State2[[i]]$x, curve1$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve1$State3[[i]]$x, curve1$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "Hypothetical stationary distribution (Setting II)", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(hypothetical$State1[[i]]$x, hypothetical$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical$State2[[i]]$x, hypothetical$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical$State3[[i]]$x, hypothetical$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)


### Flexible Dirichlet Regression ----

num_simulations <- 200

gam_results1 <- with_progress({
  p <- progressor(steps = num_simulations)
  lapply(1:num_simulations, function(i) {
    res <- oneDirGAM(n = n, rho=rho, mu = mu, sig = sig, beta = beta, periodic = periodic, par=par)
    p()
    res
  })
})

gamcurve1 <- list(State1 = list(), State2 = list(), State3 = list())
for (i in 1:num_simulations) {
  gamcurve1$State1[[i]] <- gam_results1[[i]]$State1
  gamcurve1$State2[[i]] <- gam_results1[[i]]$State2
  gamcurve1$State3[[i]] <- gam_results1[[i]]$State3
}

trim_to_range <- function(x, y, xmin = -4.05, xmax = 4.05) {
  keep <- x >= xmin & x <= xmax
  list(x = x[keep], y = y[keep])
}

ks1_1 <- ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1)
cut1_1 <- trim_to_range(ks1_1$x, ks1_1$y)
ks2_1 <- ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1)
cut2_1 <- trim_to_range(ks2_1$x, ks2_1$y)
ks3_1 <- ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1)
cut3_1 <- trim_to_range(ks3_1$x, ks3_1$y)


par(mfrow=c(1,3), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3, cex.main=1.3)
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (I)", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve1$State1[[i]]$x, gamcurve1$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve1$State2[[i]]$x, gamcurve1$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve1$State3[[i]]$x, gamcurve1$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(cut1_1$x, cut1_1$y, col = colour[1], lwd = 3)
lines(cut2_1$x, cut2_1$y, col = colour[2], lwd = 3)
lines(cut3_1$x, cut3_1$y, col = colour[3], lwd = 3)
legend("topleft",col = c(colour, "transparent"),lwd = 3,bty = "n",
       lty = c(1,1,1),legend = expression(state~1, state~2, state~3))

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (II)", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve2$State1[[i]]$x, gamcurve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve2$State2[[i]]$x, gamcurve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve2$State3[[i]]$x, gamcurve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}


lines(cut1$x, cut1$y, col = colour[1], lwd = 3)
lines(cut2$x, cut2$y, col = colour[2], lwd = 3)
lines(cut3$x, cut3$y, col = colour[3], lwd = 3)


plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (III)", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:num_simulations) {
  lines(gamcurve3$State1[[i]]$x, gamcurve3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve3$State2[[i]]$x, gamcurve3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve3$State3[[i]]$x, gamcurve3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(cut1_3$x, cut1_3$y, col = colour[1], lwd = 3)
lines(cut2_3$x, cut2_3$y, col = colour[2], lwd = 3)
lines(cut3_3$x, cut3_3$y, col = colour[3], lwd = 3)

