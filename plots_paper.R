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
       ylab = "Pr(state | z)",main = main_title,bty = "n")
  
  for (state in 1:ncol(Delta)) {
    points(jitter(z, 0), Delta[, state],
           col = alpha(colour[state], 0.05), pch = 20)
    
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
# pdf("figures/talk_samehypo_differentpersistence.pdf", width=8, height=5)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(ts(sim_high$z[1:500]), type = "l", xlab = "time", ylab = "covariate value z", main = "High persistence (first 500 obs.)", bty="n")
points(84, sim_high$z[84], pch=16, col="red", cex=1)
plot(ts(sim_low$z[1:500]), type = "l", xlab = "time", ylab = "covariate value z", main = "Low persistence", bty="n")
plot_state_probs(z_high, Delta_high, zseq, Deltaseq, "")
points(rep(sim_high$z[84],3), Delta_high[84, ], pch=16, col=colour, cex=1)
points(rep(sim_high$z[84],3), Delta_high[84, ], pch=1, col="red", cex=1)
plot_state_probs(z_low, Delta_low, zseq, Deltaseq,"")
# dev.off()

### Plot showing time series with two highlighted points (for talk)
# pdf("figures/talk_time_series.pdf", width=9, height=4)
par(mfrow=c(1,1),mar=c(3,1,1,4))
plot(ts(sim_low$z[25:100]), type="o", axes=FALSE, ann=FALSE, bty="n", col="black",bg="white", pch=21)

abline(v=6+25,  lty=2, col="skyblue2", lwd=2)
abline(v=14+25, lty=2, col="skyblue2", lwd=2)

text(6+26,  sim_low$z[55]+0.1,  labels=expression(z[80]),  pos=2.9, col="skyblue2", cex=1.2)
points(6+25, sim_low$z[55], pch=16, col="skyblue2", cex=1.1)
#cols_or <- colorRampPalette(c("#FFFF99", "darkorange"))(15)
#cols_wh <- colorRampPalette(c("white", "#FFFF99"))(17)[1:16]
cols_or <- alpha("skyblue2", seq(0, 1, length.out=31))
points(1:31, sim_low$z[25:55] , pch = 16, col = cols_or, cex = 1.1)
text(14+30, sim_low$z[63], labels=expression(z[88]), pos=2.9, col="skyblue2", cex=1.2)
points(14+25, sim_low$z[63], pch=16, col="skyblue2", cex=1.1)

arrows(x0 = 1, y0 = min(sim_low$z[25:101]) - 0.05 * diff(range(sim_low$z[25:101])),
       x1 = 77, y1 = min(sim_low$z[25:101]) - 0.05 * diff(range(sim_low$z[25:101])),
       length = 0.12, xpd = TRUE)

text(75, min(sim_low$z[25:101]) - 0.13*diff(range(sim_low$z[25:101])),
     labels="time", xpd=TRUE, cex=1.2)
# dev.off()



## creating realistic paths backsampled for given time point
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
mu1 = 0
(phi_hat = fit$ar)

(sigma_hat = sqrt(fit$var.pred))

### figure comparing hypothetical stationary distribution assumption vs. covariate sequence
# + realistic paths backsampled from fitted AR model, for given time point 


k <- 20 
B <- 100 # number of paths to simulate
y <- sim_low$z[230:249]
x0 <- sim_low$z[249]
par <- cond_dist_ARp(mu1, sigma_hat, phi_hat, k, x0)
xsim <- mvtnorm::rmvnorm(B, par$mu_cond, par$Sigma_cond)
paths <- cbind(xsim, rep(x0,B))
ylim <- range(c(y, paths))
# colours are brightening the further we go in past, to illustrate the decreasing influence of past values
#cols <- colorRampPalette(c("skyblue", "darkblue"))(20)
cols <- alpha("darkblue", seq(0.1, 1, length.out=20))


# pdf("figures/paper_paths_hypothetical_real_fade.pdf", width=7, height=3)
par(mfrow=c(1,2), mar=c(3,1,1,1))
y <- rep(sim_low$z[249], times = 20)
# covariate path under hypothetical stationary distribution assumption 
plot(ts(y),type = "n",axes = FALSE,ann = FALSE,bty = "n",ylim=ylim)
for(i in 1:19){lines(i:(i+1), y[i:(i+1)], col = cols[i], lwd = 2)}
points(1:20, y,pch = 16,col = "white",cex = 1)
points(1:20, y, pch = 16, col = cols, cex = 1)

y <- sim_low$z[230:249]

# actual covariate path + alternative paths backsampled from fitted AR model, for given time point
plot(ts(y), type="n", axes=FALSE, ann=FALSE, bty="n", ylim=ylim)
for(b in 1:10)
  lines(1:k, paths[b,], col="lightgrey")
for(i in 1:(length(y)-1))
  lines(i:(i+1), y[i:(i+1)], col=cols[i], lwd=2)
points(1:20, y,pch = 16,col = "white",cex = 1)
points(1:length(y), y, pch=16, col=cols, cex=1)
# dev.off()

# added formulas to figure with keynote 'create_fig_hyp_path.key'

## figure to explain resampling procedure
par(mfrow=c(1,1))
n_bins <- 50
z_breaks <- seq(min(sim_low$z)-0.2, max(sim_low$z)+0.2, length.out = n_bins + 1)
z_bins <- cut(sim_low$z, breaks = z_breaks, include.lowest = TRUE)
bin_colors <- colorRampPalette(c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"))(n_bins)
point_colors <- bin_colors[as.numeric(z_bins)]

# pdf("figures/bins_orig_cov.pdf", width=6, height=5)
par(mar = c(4, 4, 3, 2), cex.lab = 1.3, cex.axis = 1.1)
plot(1:(n/2), sim_low$z[1:1000], type = "n",xlab = "time",ylab = "covariate value z",
     main = "",cex.main = 1.2,bty = "l", ylim=c(min(sim_low$z)-0.2, max(sim_low$z)+0.2))
abline(h=z_breaks, col="gray90", lty=1)

# the line and points
lines(1:(n/2), sim_low$z[1:1000], col = "gray50", lwd = 1)
points(1:(n/2), sim_low$z[1:1000], 
       col = point_colors, 
       pch = 16, 
       cex = 0.8)
# dev.off()
par = list(beta = c(rep(-1, 6), rep(0,6)),
        delta = c(0,0,1), 
        mu = c(6, 15, 20), 
        sigma = log(sig))

ar_sim <- applyAr(sim_low, par, 1000)


par(mfrow=c(1,1), mar=c(4,4,2,1))
plot_state_probs(z_low, Delta_low, zseq, Deltaseq,
                 "Setting II: low persistence")

ar_sim <- applyAr(sim_low, par, 400)

par(mfrow=c(1,1), mar=c(4,4,2,4), cex.lab = 1.3, cex.axis = 1.1)

# pdf("figures/paper_bins_state_probs.pdf", width=6, height=5)
plot(ar_sim$result$State1$x, ar_sim$result$State1$y, type="n", ylim=c(0,1), 
     xlab="covariate value z", ylab="Pr(state i | z)",bty = "l", cex.main = 1.2)

draw_smooth_gradient_line <- function(x, y, colors) {
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]
  
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  smoothed <- ksmooth(x, y, kernel="normal", bandwidth=1, x.points=x)
  
  for(i in 1:(length(smoothed$x)-1)) {
    lines(smoothed$x[i:(i+1)], smoothed$y[i:(i+1)], col=colors[i], lwd=2)
  }
  
  return(smoothed)
}

smooth1 <- draw_smooth_gradient_line(ar_sim$result$State1$x, ar_sim$result$State1$y, bin_colors)
smooth2 <- draw_smooth_gradient_line(ar_sim$result$State2$x, ar_sim$result$State2$y, bin_colors)
smooth3 <- draw_smooth_gradient_line(ar_sim$result$State3$x, ar_sim$result$State3$y, bin_colors)

text(max(smooth1$x)-1.1, smooth1$y[length(smooth1$y)]+0.03, "State 1", pos=4, cex=1)
text(max(smooth2$x)-1.1, smooth2$y[length(smooth2$y)]-0.03, "State 2", pos=4, cex=1)
text(max(smooth3$x)-1.1, smooth3$y[length(smooth3$y)]+0.03, "State 3", pos=4, cex=1)
# dev.off()

