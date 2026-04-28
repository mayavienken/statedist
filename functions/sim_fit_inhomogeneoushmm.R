library(LaMa)

colour = c("orange", "skyblue", "seagreen")

cycle_length = 24
base_amplitude = 4
amplitude_var = 0.6
noise_sd = 1

# simulate HMM with covariate-driven transition probabilities
simCovHMM <- function(n, mu, sig, beta, periodic = TRUE,
                      rho = 0.95, cycle_length = 24,
                      base_amplitude = 4, amplitude_var = 0.6, noise_sd = 1) {
  z <- numeric(n) # covariate series

  if (periodic) {
    # # generate periodic covariate with noise
    t <- 0:(n - 1)
    amp <- base_amplitude + rnorm(n, 0, amplitude_var)
    z <- amp * sin(2 * pi * t / cycle_length) + rnorm(n, 0, noise_sd)
  } else {
    # generate AR(1) covariate process
    z[1] <- rnorm(1)
    for (i in 2:n) {
      z[i] <- rho * z[i - 1] + rnorm(1, sd = sqrt(1 - rho^2))
    }
  }

  Z <- matrix(z)
  Gamma <- tpm_g(Z[-1, ], beta)
  
  # # simulate hidden states
  delta <- rep(1/3, 3)
  s <- numeric(n)
  s[1] <- sample(1:3, 1, prob = delta)
  for (i in 2:n) {
    s[i] <- sample(1:3, 1, prob = Gamma[s[i - 1], , i - 1])
  }
  
  # simulate observations from state-dependent normal distributions
  x <- rnorm(n, mean = mu[s], sd = sig[s])

  return(list(x = x, s = s, z = z))
}

# fit HMM with covariate-dependent transitions via maximum likelihood
fitCovHMM <- function(par, x, Z) {
  
  # negative log-likelihood using forward algorithm
  nll_cov <- function(par, x, Z) {
    
    beta <- matrix(par[1:12], nrow = 6)
    Gamma <- tpm_g(Z[-1, ], beta)
    delta <- c(1, exp(par[13]), exp(par[14]))
    delta <- delta / sum(delta)
    
    mu <- par[15:17]        # state means
    sig <- exp(par[18:20])  # state standard deviations
    
    # emission probabilities
    allprobs <- matrix(1, length(x), 3)
    for (j in 1:3) allprobs[, j] <- dnorm(x, mu[j], sig[j])
    
    -forward_g(delta, Gamma, allprobs)
  }
  
  # optimise likelihood
  mod <- nlm(nll_cov, par, x = x, Z = Z, hessian = TRUE)
  
  # extract parameter estimates
  beta_hat <- matrix(mod$estimate[1:12], nrow = 6)
  Gamma_hat <- tpm_g(Z[-1, ], beta_hat)
  delta_hat <- c(1, exp(mod$estimate[13]), exp(mod$estimate[14]))
  delta_hat <- delta_hat / sum(delta_hat)
  mu_hat <- mod$estimate[15:17]
  sig_hat <- exp(mod$estimate[18:20])
  
  return(list(
    beta = beta_hat,
    #Gamma = Gamma_hat,
    #delta = delta_hat,
    mu = mu_hat,
    sig = sig_hat,
    logLik = -mod$minimum,
    convergence = mod$code,
    hessian = mod$hessian
  ))
}


