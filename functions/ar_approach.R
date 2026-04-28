
source("./functions/sim_fit_inhomogeneousHMM.r")

library(LaMa) 
library(scales)
library(parallel)
library(progressr)
library(zoo)

# simulate covariate series using fitted AR model
simAr <- function(n, z) {
  
  fit <- ar(ts(z))            # fit AR model to observed covariate series 
  phi <- fit$ar               # extract AR coefficients
  sigma <- sqrt(fit$var.pred) # innovation standard deviation
  
  z_gen <- numeric(n)
  p <- length(phi)            # AR order
  
  z_gen[1:p] <- rnorm(p, mean(z), sigma) # initialise with random draws
  
  for (t in (p + 1):n) {
    # generate new values from AR recursion
    z_gen[t] <- sum(phi * rev(z_gen[(t - p):(t - 1)])) + rnorm(1, 0, sigma)
  }
  return(z_gen)
}

# compute state probabilities for a given covariate sequence
compStateProbs <- function(z, beta, n) {
  
  Gamma <- tpm_g(z, beta)
  
  Delta <- matrix(NA, n, 3)
  Delta[1, ] <- rep(1/3, 3) # initialise with uniform distribution
  
  for (t in 2:n) {
    # propagate state probabilities forward 
    Delta[t, ] <- Delta[t-1, ] %*% Gamma[, , t]
  }
  return(Delta)
}


simOneAr <- function(n, rho, mu, sig, beta, periodic, par, num_covsim) {
  
  sim <- simCovHMM(n = n, rho = rho, mu = mu, sig = sig, beta = beta, periodic = periodic)
  fit <- fitCovHMM(par = par, x = sim$x, Z = matrix(sim$z))

  sim_z <- mclapply(1:num_covsim, function(i) simAr(n=n, sim$z), mc.cores = max(1, detectCores() - 2))
  sim_z <- do.call(cbind, sim_z)
  
  sim_delta <- mclapply(1:num_covsim, function(i) compStateProbs(sim_z[, i], fit$beta, n=n), mc.cores = max(1, detectCores() - 2))
  sim_delta <- array(unlist(sim_delta), dim = c(n=n, 3, num_covsim))
  
  z_bins <- seq(min(sim_z), max(sim_z), length.out = 30)
  bin_midpoints <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
  
  mean_state_probs <- matrix(NA, nrow = length(bin_midpoints), ncol = 3)
  for (state in 1:3) {
    mean_state_probs[, state] <- sapply(1:(length(z_bins) - 1), function(b) {
      mean(sim_delta[, state, ][sim_z >= z_bins[b] & sim_z < z_bins[b + 1]], na.rm = TRUE)
    })
  }
  
  result <- list(
    State1 = data.frame(x = bin_midpoints, y = mean_state_probs[, 1]),
    State2 = data.frame(x = bin_midpoints, y = mean_state_probs[, 2]),
    State3 = data.frame(x = bin_midpoints, y = mean_state_probs[, 3])
  )
  
  return(result)
}


# compute stationary state probabilities over covariate grid 
# (hypothetical stationary distribution (Patterson et al., 2009))
fithypothetical <-  function(n, rho, mu, sig, beta, periodic, par, num_covsim) {
  
  sim <- simCovHMM(n = n, rho = rho, mu = mu, sig = sig, beta = beta, periodic = periodic)
  fit <- fitCovHMM(par = par, x = sim$x, Z = matrix(sim$z))
  
  # define covariate grid
  zseq <- seq(min(sim$z), max(sim$z), length = 200)
  
  Gammaseq <- tpm_g(zseq, fit$beta)
  Deltaseq <- matrix(NA, length(zseq), N)
  
  # compute stationary distribution at each covariate value
  for (t in 1:length(zseq)) Deltaseq[t,] <- LaMa::stationary(Gammaseq[,,t])
  
  result <- list(
    State1 = data.frame(x = zseq, y = Deltaseq[, 1]),
    State2 = data.frame(x = zseq, y = Deltaseq[, 2]),
    State3 = data.frame(x = zseq, y = Deltaseq[, 3])
  )
  
  return(result)
}
  

simOneArHypothetical <- function(n, rho, mu, sig, beta, periodic, par, num_covsim) {
  
  sim <- simCovHMM(n = n, rho = rho, mu = mu, sig = sig, beta = beta, periodic = periodic)
  fit <- fitCovHMM(par = par, x = sim$x, Z = matrix(sim$z))
  
  sim_z <- mclapply(1:num_covsim, function(i) simAr(n=n, sim$z), mc.cores = max(1, detectCores() - 2))
  sim_z <- do.call(cbind, sim_z)
  
  sim_delta <- mclapply(1:num_covsim, function(i) compStateProbs(sim_z[, i], fit$beta, n=n), mc.cores = max(1, detectCores() - 2))
  sim_delta <- array(unlist(sim_delta), dim = c(n=n, 3, num_covsim))
  
  z_bins <- seq(min(sim_z), max(sim_z), length.out = 30)
  bin_midpoints <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
  
  mean_state_probs <- matrix(NA, nrow = length(bin_midpoints), ncol = 3)
  for (state in 1:3) {
    mean_state_probs[, state] <- sapply(1:(length(z_bins) - 1), function(b) {
      mean(sim_delta[, state, ][sim_z >= z_bins[b] & sim_z < z_bins[b + 1]], na.rm = TRUE)
    })
  }
  
  resultar <- list(
    State1 = data.frame(x = bin_midpoints, y = mean_state_probs[, 1]),
    State2 = data.frame(x = bin_midpoints, y = mean_state_probs[, 2]),
    State3 = data.frame(x = bin_midpoints, y = mean_state_probs[, 3])
  )
  
  zseq <- seq(min(sim$z), max(sim$z), length = 200)
  Gammaseq <- tpm_g(zseq, fit$beta)
  Deltaseq <- matrix(NA, length(zseq), N)
  for (t in 1:length(zseq)) Deltaseq[t,] <- LaMa::stationary(Gammaseq[,,t])
  
  resulthypothetical <- list(
    State1 = data.frame(x = zseq, y = Deltaseq[, 1]),
    State2 = data.frame(x = zseq, y = Deltaseq[, 2]),
    State3 = data.frame(x = zseq, y = Deltaseq[, 3])
  )
  
  return(list(
    resultar = resultar,
    resulthypothetical = resulthypothetical
  ))
}



block_bootstrap <- function(ts, block_size, n) {
  num_blocks <- ceiling(n / block_size)
  block_starts <- sample(1:(length(ts) - block_size + 1), num_blocks, replace = TRUE)
  
  bootstrap_sample <- unlist(lapply(block_starts, function(start) ts[start:(start + block_size - 1)]))
  return(bootstrap_sample[1:n])  # trim to match original length
}


simOneBB <- function(n, rho, mu, sig, beta, periodic, par, num_covsim) {
  
  sim <- simCovHMM(n = n, rho = rho, mu = mu, sig = sig, beta = beta, periodic = periodic)
  fit <- fitCovHMM(par = par, x = sim$x, Z = matrix(sim$z))
  
  sim_z <- mclapply(1:num_covsim, function(i) block_bootstrap(sim$z, 24, n), mc.cores = max(1, detectCores() - 2))
  sim_z <- do.call(cbind, sim_z)
  
  sim_delta <- mclapply(1:num_covsim, function(i) compStateProbs(sim_z[, i], fit$beta, n=n), mc.cores = max(1, detectCores() - 2))
  sim_delta <- array(unlist(sim_delta), dim = c(n=n, 3, num_covsim))
  
  z_bins <- seq(min(sim_z), max(sim_z), length.out = 30)
  bin_midpoints <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
  
  mean_state_probs <- matrix(NA, nrow = length(bin_midpoints), ncol = 3)
  for (state in 1:3) {
    mean_state_probs[, state] <- sapply(1:(length(z_bins) - 1), function(b) {
      mean(sim_delta[, state, ][sim_z >= z_bins[b] & sim_z < z_bins[b + 1]], na.rm = TRUE)
    })
  }
  #mean_state_probs <- na.approx(mean_state_probs, x = bin_midpoints, rule = 2)
  
  result <- list(
    State1 = data.frame(x = bin_midpoints, y = mean_state_probs[, 1]),
    State2 = data.frame(x = bin_midpoints, y = mean_state_probs[, 2]),
    State3 = data.frame(x = bin_midpoints, y = mean_state_probs[, 3])
  )
  
  return(result)
}

applyAr <- function(sim, par, num_covsim) {
  z <- sim$z
  x <- sim$x
  fit <- fitCovHMM(par = par, x = x, Z = matrix(z))
  
  sim_z <- mclapply(1:num_covsim, function(i) simAr(n=n, z), mc.cores = max(1, detectCores() - 2))
  sim_z <- do.call(cbind, sim_z)
  
  sim_delta <- mclapply(1:num_covsim, function(i) compStateProbs(sim_z[, i], fit$beta, n=n), mc.cores = max(1, detectCores() - 2))
  sim_delta <- array(unlist(sim_delta), dim = c(n=n, 3, num_covsim))
  
  z_bins <- seq(min(sim_z), max(sim_z), length.out = 30)
  bin_midpoints <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
  
  mean_state_probs <- matrix(NA, nrow = length(bin_midpoints), ncol = 3)
  for (state in 1:3) {
    mean_state_probs[, state] <- sapply(1:(length(z_bins) - 1), function(b) {
      mean(sim_delta[, state, ][sim_z >= z_bins[b] & sim_z < z_bins[b + 1]], na.rm = TRUE)
    })
  }
  
  result <- list(
    State1 = data.frame(x = bin_midpoints, y = mean_state_probs[, 1]),
    State2 = data.frame(x = bin_midpoints, y = mean_state_probs[, 2]),
    State3 = data.frame(x = bin_midpoints, y = mean_state_probs[, 3])
  )
  
  return(result)
}
