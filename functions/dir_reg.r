# Dirichlet regression with spline-based covariate effects

source("./functions/sim_fit_inhomogeneousHMM.r")

library(LaMa).    # qREML for penalty selection
library(scales).  # for colour transparency
library(RTMBdist) # Dirichlet density
library(splines2) # B-spline design matrix
library(Matrix)   
library(progressr)

# Dirichlet regression with spline-based covariate effect
dir_reg <- function(y, x, k = 10, bs = "tp", lambda0 = 1e3, alpha = 0.1, silent = 1) {
  
  # define spline formula for covariate effect 
  form <- eval(substitute(~ s(x, bs = bs, k = k), list(bs = bs, k = k)))
  
  # to construct design and penalty matrices for models involving penalised splines
  modmat <- make_matrices(form, data = data.frame(x = x)) 
  X <- modmat$Z       # spline design matrix
  S <- modmat$S[[1]]  # penalty matrix
  
  # penalised negative log-likelihood
  pnll <- function(par) {
    getAll(par, dat, warn = FALSE)
    
    # compute positive Dirichlet parameters via exponential link
    Alpha <- exp(X %*% t(cbind(beta0, beta))) 
    # combine likelihood with smoothness penalty (proportional to spline roughness)
    -sum(ddirichlet(y, Alpha, log = TRUE), na.rm = TRUE) + 
      penalty(beta, S, lambda) 
  }
  # initialise model parameters 
  par <- list(
    beta0 = rep(0, ncol(y)),
    beta = matrix(0, ncol(y), ncol(X)-1))
  
  dat <- list(y = y, X = X, lambda = rep(lambda0, ncol(y)), S = S)
  
  # internal handling of matrix operations (sparsity)
  TapeConfig(matmul = "plain") 
  
  # fit model using restricted maximum likelihood
  mod <- qreml(pnll, par, dat, random = "beta", 
               alpha = alpha, silent = silent)
  # define prediction method for new covariate values
  mod$predict <- function(x_p) {
    X_p <- pred_matrix(modmat, newdata = data.frame(x = x_p))     # new spline design matrix
    Alpha_p <- exp(X_p %*% t(cbind(mod$par$beta0, mod$par$beta))) # computes predicted concentration parameters
    Mean_p <- Alpha_p / rowSums(Alpha_p)                          # normalised to obtain Dirichlet mean proportions
    return(Mean_p)
  }
  return(mod)
}

# wrapper: fit model and produce diagnostic plot
apply_dir_reg <- function(y, x, k = 8, bs = "tp", lambda0 = 100, colour = c("orange", "skyblue", "seagreen"), N, x_range_extension = 0.5, covname="covariate") {
  
  # fit Dirichlet regression model
  mod <- dir_reg(y, x, k = k, bs = bs, lambda0 = lambda0)
  
  # generate grid of covariate values for smooth curves
  x_p <- seq(min(x) - x_range_extension, max(x) + x_range_extension, length = 200)
  Mean <- mod$predict(x_p)
  
  # initialise plot
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
  
  # plot observed proportions with transparency
  for (state in 1:N) {
    points(jitter(x, factor = 0), y[, state], col = alpha(colour[state], 0.1), pch = 16)
  }
  
  # overlay fitted smooth curves
  for (state in 1:N) {
    darker_col <- adjustcolor(colour[state], red.f = 0.9, green.f = 0.9, blue.f = 0.9, alpha.f = 1)
    lines(x_p, Mean[, state], col = "white", lwd = 7, lty = 1)
    lines(x_p, Mean[, state], col = darker_col, lwd = 3, lty = 1)
  }
  
  # return model and predictions for further use
  return(list(model = mod, x_p = x_p, predicted_mean = Mean))
}
