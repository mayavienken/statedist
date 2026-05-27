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
fitCovHMM <- function(par, dat) {
  
  # negative log-likelihood using forward algorithm
  nll_cov <- function(par) {
    getAll(par, dat) # make everything accessible without $
    
    beta <- matrix(beta, nrow = N*(N-1))
    Gamma <- tpm_g(Z[-1, ], beta) # transition probability matrices
    
    delta <- c(1, exp(delta[1]), exp(delta[2]))
    delta <- delta / sum(delta)
    
    mu <- mu        # state-dependent means
    REPORT(mu); ADREPORT(mu)
    sigma <- exp(sigma)  # state-dependent standard deviations
    REPORT(sigma); ADREPORT(sigma)
    
    # emission probabilities
    allprobs <- matrix(1, length(x), N) 
    for (j in 1:3) allprobs[, j] <- dnorm(x, mu[j], sigma[j])
    
    -forward_g(delta, Gamma, allprobs)
  }
  
  # create automatically differentiable objective
  obj = MakeADFun(nll_cov, par, silent = TRUE)
  
  # optimise likelihood
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  # report results
  mod = report(obj)
  # extract parameter estimates
  beta_hat <- matrix(mod$beta, nrow = 6)
  Gamma_hat <- tpm_g(dat$Z[-1, ], beta_hat)
  delta_hat <- c(1, mod$delta)
  delta_hat <- delta_hat / sum(delta_hat)
  mu_hat <- mod$mu
  sig_hat <- mod$sigma
  
  return(list(
    beta = beta_hat,
    #Gamma = Gamma_hat,
    #delta = delta_hat,
    mu = mu_hat,
    sig = sig_hat,
    logLik = -mod$ll
  ))
}
