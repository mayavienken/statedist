### Script for creating plots for the paper (and talks) for further illustration
source("functions/sim_fit_inhomogeneousHMM.r")
source("functions/ar_approach.r")

library(scales)
library(LaMa)
library(mvtnorm)

set.seed(22)

# parameters
n <- 2000 
rho_high <- 0.95
rho_low  <- 0.7
N <- 3

mu  <- c(6, 15, 20)
sig <- c(3, 1.5, 1.5)

beta <- matrix(c(
  -1,  2,
  -1, -2.5,
  -1,  0.5,
  -1, -2.5,
  -1, -1,
  -1,  1
), nrow = 6, byrow = TRUE)

periodic <- FALSE

# compute state probabilities for a given covariate sequence
compute_delta <- function(z, beta, N) {
  n <- length(z)
  Gamma <- tpm_g(z, beta)
  
  Delta <- matrix(NA, n, N)
  Delta[1, ] <- rep(1 / N, N)
  
  for (t in 2:n) {
    Delta[t, ] <- Delta[t - 1, ] %*% Gamma[, , t]
  }
  
  list(Gamma = Gamma, Delta = Delta)
}

# hypothetical stationary distribution (stationary state probabilities)
compute_stationary_curve <- function(zseq, beta, N) {
  Gammaseq <- tpm_g(zseq, beta)
  
  Deltaseq <- matrix(NA, length(zseq), N)
  for (t in seq_along(zseq)) {
    Deltaseq[t, ] <- LaMa::stationary(Gammaseq[, , t])
  }
  
  Deltaseq
}

# simulate two differnt data sets
# both use the identical parameters, except for the persistence of the covariate process (rho)
sim_high <- simCovHMM(n=n, rho=rho_high, mu=mu, sig=sig, beta=beta, periodic=periodic)
sim_low  <- simCovHMM(n=n, rho=rho_low, mu=mu, sig=sig, beta=beta, periodic=periodic)

res_high <- compute_delta(sim_high$z, beta, N)
res_low  <- compute_delta(sim_low$z,  beta, N)

# burn-in removal
burn <- 10

z_high <- sim_high$z[burn:n]
z_low  <- sim_low$z[burn:n]

Delta_high <- res_high$Delta[burn:n, ]
Delta_low  <- res_low$Delta[burn:n, ]

# both range somewhat from -4 to 4 so we can use the same grid for the hypothetical stationary distribution
zseq <- seq(-4, 4, length = 200)
Deltaseq <- compute_stationary_curve(zseq, beta, N)

plot_state_probs <- function(z, Delta, zseq, Deltaseq, main_title) {
  plot(0, 0,xlim = c(-4, 4.5), ylim = c(0, 1),
       col = "transparent",xlab = "covariate value z",
       ylab = "Pr(state i | z)",main = main_title,bty = "n")
  
  for (state in 1:ncol(Delta)) {
    points(jitter(z, 0), Delta[, state],
           col = alpha(colour[state], 0.1), pch = 20)
    
    lines(zseq, Deltaseq[, state],
          col = colour[state], lwd = 3, lty = 2)
  }
}

# comparison of just the hypotheticals
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot_state_probs(z_high, Delta_high, zseq, Deltaseq,
                 "Setting I: high persistence")
plot_state_probs(z_low, Delta_low, zseq, Deltaseq,
                 "Setting II: low persistence")

# 2x2 plot with time series and state probabilities + hypotheticals
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(ts(sim_high$z[1:500]), type = "l", xlab = "time", ylab = "covariate value z", main = "Setting I: high persistence (first 500 obs.)", bty="n")
plot(ts(sim_low$z[1:500]), type = "l", xlab = "time", ylab = "covariate value z", main = "Setting II: low persistence", bty="n")
plot_state_probs(z_high, Delta_high, zseq, Deltaseq, "")
plot_state_probs(z_low, Delta_low, zseq, Deltaseq,"")


### Plot showing time series with two highlighted points (for talk)
# pdf("figures/talk_time_series.pdf", width=9, height=4)
par(mfrow=c(1,1),mar=c(3,1,1,4))
plot(ts(sim_low$z[25:100]), type="o", axes=FALSE, ann=FALSE, bty="n", col="black",bg="white", pch=21)

abline(v=6+25,  lty=2, col="darkorange", lwd=2)
abline(v=14+25, lty=2, col="darkorange", lwd=2)

text(6+26,  sim_low$z[55]+0.1,  labels=expression(z[80]),  pos=2.9, col="darkorange", cex=1.2)
points(6+25, sim_low$z[55], pch=16, col="darkorange", cex=1.1)
cols_or <- colorRampPalette(c("#FFFF99", "darkorange"))(15)
cols_wh <- colorRampPalette(c("white", "#FFFF99"))(17)[1:16]
points(17:31, sim_low$z[41:55] , pch = 16, col = cols_or, cex = 1.1)
points(1:16, sim_low$z[25:40] , pch = 16, col = cols_wh, cex = 1.1)
text(14+30, sim_low$z[63], labels=expression(z[88]), pos=2.9, col="darkorange", cex=1.2)
points(14+25, sim_low$z[63], pch=16, col="darkorange", cex=1.1)

arrows(x0 = 1, y0 = min(sim_low$z[25:101]) - 0.05 * diff(range(sim_low$z[25:101])),
       x1 = 77, y1 = min(sim_low$z[25:101]) - 0.05 * diff(range(sim_low$z[25:101])),
       length = 0.12, xpd = TRUE)

text(75, min(sim_low$z[25:101]) - 0.13*diff(range(sim_low$z[25:101])),
     labels="time", xpd=TRUE, cex=1.2)
# dev.off()

# autocovariances based on Yule Walker equations
gen_autocovariances = function(phi, sigma2, max_lag) {
  p = length(phi)
  gamma = numeric(max_lag + 1)
  gamma[1] = sigma2 / (1 - sum(phi * phi))  # Variance at lag 0
  for (h in 1:max_lag) {
    if (h <= p) {
      gamma[h + 1] = sum(phi[1:h] * rev(gamma[1:h]))
    } else {
      gamma[h + 1] = sum(phi * rev(gamma[(h - p + 1):h]))
    }
  }
  return(gamma)
}
# Generate the covariance matrix for AR(p) process
gen_Sigma_ARp = function(sigmasq, phi, k) {
  autocovariances = gen_autocovariances(phi, sigmasq, k - 1)
  Sigma = matrix(0, nrow = k, ncol = k)
  for (i in 1:k) { for (j in 1:k) {
    Sigma[i, j] = autocovariances[abs(i - j) + 1]
  }}
  return(Sigma)
}
# Calculate the conditional distribution for AR(p) process
cond_dist_ARp = function(mu, sigma, phi, k, x) {
  Sigma = gen_Sigma_ARp(sigma^2, phi, k)
  mu_cond = mu + Sigma[-k, k] / Sigma[k, k] * (x - mu)
  Sigma_cond = Sigma[-k, -k] - Sigma[-k, k] %*% t(Sigma[k, -k]) / Sigma[k, k]
  return(list(mu_cond = mu_cond, Sigma_cond = Sigma_cond))
}

# AR fit with automatic order selection
fit = ar(ts(sim_low$z))
mu = 0
(phi_hat = fit$ar)

(sigma_hat = sqrt(fit$var.pred))




### comparing hypothetical stationary observation assumption vs. number of possible observed paths 

# pdf("figures/paper_paths_hypothetical_real.pdf", width=7, height=3)
par(mfrow=c(1,2), mar=c(3,1,1,1))

ylim <- range(c(y, paths))
cols <- colorRampPalette(c("skyblue", "darkblue"))(20)
y <- rep(sim_low$z[100], times = 20)
plot(ts(y),type = "n",axes = FALSE,ann = FALSE,bty = "n",ylim=ylim)
for(i in 1:19){lines(i:(i+1), y[i:(i+1)], col = cols[i], lwd = 2)}
points(1:20, y, pch = 16, col = cols, cex = 1.1)

cols <- colorRampPalette(c("skyblue", "darkblue"))(20)
y <- sim_low$z[81:100]

k <- 20
B <- 100
par <- cond_dist_ARp(mu, sigma_hat, phi_hat, k, x0)
xsim <- mvtnorm::rmvnorm(B, par$mu_cond, par$Sigma_cond)
paths <- cbind(xsim, rep(x0,B))

plot(ts(y), type="n", axes=FALSE, ann=FALSE, bty="n", ylim=ylim)
for(b in 1:10)
  lines(1:k, paths[b,], col="lightgrey")
for(i in 1:(length(y)-1))
  lines(i:(i+1), y[i:(i+1)], col=cols[i], lwd=2)
points(1:length(y), y, pch=16, col=cols, cex=1.1)
# dev.off()


