### Functions used for the simulation study 

source("./functions/sim_fit_inhomogeneousHMM.r")
source("./functions/ar_approach.r")
source("./functions/bb_approach.r")
source("./functions/dir_reg.r")

# for simulation study
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


## Dirichlet regression approach
oneDirGAM <- function(n, rho, mu, sig, beta, periodic, par) {
  
  # simulate data + fit HMM
  sim <- simCovHMM(n = n, rho=rho, mu = mu, sig = sig, beta = beta, periodic = periodic)
  fit <- fitCovHMM(par = par, x = sim$x, Z = matrix(sim$z))
  
  # compute state probabilities for the observed covariate series
  z <- sim$z
  Gamma <- tpm_g(z, fit$beta)
  
  Delta <- matrix(NA, n, 3)
  Delta[1, ] <- rep(1/3, 3)
  for (t in 2:n) Delta[t, ] <- Delta[t-1, ] %*% Gamma[,,t]
  
  Y <- Delta[10:n,]
  x <- z[10:n]
  
  # fit Dirichlet regression to the state probabilities as a function of the covariate
  system.time(
    mod <- dir_reg(Y, x, k = 8, bs="tp", lambda0=100)
  )
  
  x_p <- seq(-4, 4, length = 200)
  Mean <- mod$predict(x_p)
  
  list(
    State1 = list(x = x_p, y = Mean[,1]),
    State2 = list(x = x_p, y = Mean[,2]),
    State3 = list(x = x_p, y = Mean[,3])
  )
}
