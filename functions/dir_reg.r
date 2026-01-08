source("./functions/sim_fit_inhomogeneousHMM.r")
library(LaMa) # qREML for penalty selection
library(scales)
library(RTMBdist) # dirichlet distribution
library(splines2) # B-spline design matrix
library(Matrix) 
library(progressr)

dir_reg <- function(y, x, k = 10, bs = "tp", lambda0 = 1e3, alpha = 0.1, silent = 1) {
  form <- eval(substitute(~ s(x, bs = bs, k = k), list(bs = bs, k = k)))
  modmat <- make_matrices(form, data = data.frame(x = x)) #to build design and penalty matrix for models involving penalised splines
  X <- modmat$Z # spline design matrix
  S <- modmat$S[[1]] # penalty matrix
  # Likelihood
  pnll <- function(par) {
    getAll(par, dat, warn = FALSE)
    Alpha <- exp(X %*% t(cbind(beta0, beta))) #linear predictor built (exp converts to positive Dirichlet parameters alpha)
    -sum(ddirichlet(y, Alpha, log = TRUE), na.rm = TRUE) + #negative log-likelihood plus penalty
      penalty(beta, S, lambda) #adds a smoothness penalty proportional to spline roughness
  }
  # Model fitting
  par <- list(
    beta0 = rep(0, ncol(y)),
    beta = matrix(0, ncol(y), ncol(X)-1))
  dat <- list(y = y, X = X, lambda = rep(lambda0, ncol(y)), S = S)
  TapeConfig(matmul = "plain") # internal handling of matrix operations (sparsity)
  mod <- qreml(pnll, par, dat, random = "beta", 
               alpha = alpha, silent = silent)
  mod$predict <- function(x_p) {
    X_p <- pred_matrix(modmat, newdata = data.frame(x = x_p)) # new spline design matrix
    Alpha_p <- exp(X_p %*% t(cbind(mod$par$beta0, mod$par$beta))) # computes predicted concentration parameters
    Mean_p <- Alpha_p / rowSums(Alpha_p) # normalised to obtain Dirichlet mean proportions
    return(Mean_p)
  }
  return(mod)
}

apply_dir_reg <- function(y, x, k = 8, bs = "tp", lambda0 = 100, colour = c("orange", "skyblue", "seagreen"), N, x_range_extension = 0.5, covname="covariate") {
  
  mod <- dir_reg(y, x, k = k, bs = bs, lambda0 = lambda0)
  
  x_p <- seq(min(x) - x_range_extension, max(x) + x_range_extension, length = 200)
  Mean <- mod$predict(x_p)
  
  par(mfrow = c(1,1), mgp = c(1.8, 0.5, 0), mar = c(3,3,1,1), cex.lab = 1.3)
  plot(0, 0,
       xlim = c(min(x), max(x)),
       ylim = c(0, 1),
       pch = 16,
       col = alpha(colour[1], 0),
       bty = "n",
       ylab = paste0("Pr(state i | x), i = 1,...,", N),
       xlab = paste0(covname),
       main = "Flexible Dirichlet regression")
  
  for (state in 1:N) {
    points(jitter(x, factor = 0), y[, state], col = alpha(colour[state], 0.1), pch = 16)
  }
  
  for (state in 1:N) {
    darker_col <- adjustcolor(colour[state], red.f = 0.9, green.f = 0.9, blue.f = 0.9, alpha.f = 1)
    lines(x_p, Mean[, state], col = "white", lwd = 7, lty = 1)
    lines(x_p, Mean[, state], col = darker_col, lwd = 3, lty = 1)
  }
  
  return(list(model = mod, x_p = x_p, predicted_mean = Mean))
}
