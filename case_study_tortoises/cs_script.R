## Case study on giant Galápagos tortoises

# Bastille-Rousseau G, Yackulic CB, Gibbs J, Frair JL, Cabrera F, Blake S. 2019. 
# Data from: Migration triggers in a large herbivore: Galápagos giant tortoises navigating resource gradients on volcanoes.
# Movebank Data Repository. https://doi.org/10.5441/001/1.6gr485fk

# Functions and libraries ----

source("./functions/bb_approach.R")
source("./functions/dir_reg.R")

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(future.apply)
library(parallel)
library(LaMa)         
library(RTMBdist)
library(splines2)
library(Matrix)
library(ggmap)
library(moveHMM)
library(zoo)    

# Data ----
load("./case_study_tortoises/data/df_caro.RData")
df.animal <- df.caro # movement data for one individual

# defining colours for the three states 
colour = c("orange", "skyblue", "seagreen")


# Model fitting ----

N = 3 # number of hidden states

# initial values for state-dependent distributions
mu.1 <- c(4.5,12, 40)       # step length mean (in metres)
sigma.1 <- c(2,4, 20)       # step length SD
kappa.1 <- c(1.5, 1.5, 0.5) # turning angle concentration

# parameter list for RTMB (transformations ensuring parameter constraints)
par = list(
  logmu = log(mu.1), 
  logsigma = log(sigma.1), 
  logkappa = log(kappa.1),
  beta = matrix(c(rep(-1, 6), rep(0,6)), nrow=6), 
  logitdelta = c(0, 0))

# data passed to likelihood
dat = list(
  step = df.animal$step, 
  angle = df.animal$angle, 
  Z = matrix(df.animal$temperature), 
  N = N, 
  trackID = df.animal$track_id)

# negative log-likelihood with RTMB
nll.rtmb = function(par) {
  
  getAll(par, dat)  
  Gamma = tpm_g(Z, beta) # transition matrices depending on temperature
  ADREPORT(Gamma)
  
  delta = c(1, logitdelta)
  delta = delta / sum(delta)
  ADREPORT(delta)
  
  # state-dependent parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  
  mu.kappa = c(pi, pi, 0) # angle means
  allprobs = matrix(1, nrow = length(step), ncol = N)
  
  # handle missing step/angle observations
  indStep <- which(!is.na(step))
  indAngle <- which(!is.na(angle))
  indBoth <- intersect(indStep, indAngle)
  indStep <- setdiff(indStep, indBoth)
  indAngle <- setdiff(indAngle, indBoth)
  
  for(j in 1:N) {
    allprobs[indBoth,j] = dgamma2(step[indBoth], mu[j], sigma[j]) *dvm(angle[indBoth], mu.kappa[j], kappa[j])
    allprobs[indStep,j] = dgamma2(step[indStep], mu[j], sigma[j])
    allprobs[indAngle,j] = dvm(angle[indAngle], mu.kappa[j], kappa[j])
    
  }
  
  # forward algorithm likelihood
  -forward_g(delta, Gamma, allprobs, trackID)
  
}

# model estimation
obj.rtmb = MakeADFun(nll.rtmb, par, silent = T) 
opt.rtmb = nlminb(obj.rtmb$par, obj.rtmb$fn, obj.rtmb$gr, control=list(iter.max=700, eval.max=700)) # optimisation

# extract fitted quantities 
mod.rtmb = obj.rtmb$report() 

# most likely state sequence (Viterbi)
mod.rtmb$states = viterbi_g(mod.rtmb$delta, mod.rtmb$Gamma, mod.rtmb$allprobs, df.animal$track_id)

# unpack estimated parameters
mu.hat <- mod.rtmb$mu
sigma.hat <- mod.rtmb$sigma
kappa.hat <- mod.rtmb$kappa
beta.hat <- mod.rtmb$beta
Gamma.hat <- tpm_g(matrix(df.animal$temperature), beta.hat)
delta.hat <- mod.rtmb$delta


# alternative fit using moveHMM (for CI estimation) 
moveHMMmodel <- fitHMM(data=df.animal, nbStates=3, 
                       stepPar0= c(mu.1, sigma.1), 
                       anglePar0= c(pi, pi, 0, kappa.1),
                       formula=~temperature)

# hypothetical stationary distribution with CIs
plotStationary(moveHMMmodel, plotCI=TRUE)

# transition probabilities with CIs over temperature grid
predTPM <- predictTPM(moveHMMmodel, newData = data.frame(temperature=seq(0,45,length=100)), returnCI = TRUE)
vis_z <- seq(0, 45, length=100)

#pdf("./case_study_tortoises/figures/transprobs_ci.pdf", width=8, height=5)
par(mfrow=c(3,3), mar=c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.lab=1.3)
for(i in 1:3){
  for(j in 1:3){
    
    ylab_expr <- bquote(gamma[.(as.numeric(paste0(i, j)))])
    xlab_here <- if(i == 3) "temperature" else ""
    
    plot(vis_z, predTPM$mle[i,j,],
         type="n",
         xlab=xlab_here,
         ylab=ylab_expr,
         main= bquote(Pr(.(i) %->% .(j))),
         ylim=c(0,1),
         xlim=c(0,45),
         bty="n",
         cex.lab=1.5,
         cex.axis=1.2,
         cex.main=1.5,
         xaxt="n")
    
    axis(1, at = c(0,10,20,30,40), labels = c(0,10,20,30,40))
    
    polygon(
      x = c(vis_z, rev(vis_z)),
      y = c(predTPM$lci[i,j,], rev(predTPM$uci[i,j,])),
      col = adjustcolor(colour[j], alpha.f = 0.25),
      border = NA
    )
    
    lines(vis_z, predTPM$mle[i,j,], lwd = 2, col = colour[j])
    
  }
}
#dev.off()

# Plots ----

# GPS track coloured by decoded states 
register_stadiamaps(key = Sys.getenv("STADIAMAPS_KEY"))

map_bg <- get_stadiamap(
  bbox = c(left = -90.48167, bottom = -0.7252189, 
           right = -90.39135, top = -0.658182),
  zoom = 15, 
  maptype = "stamen_terrain_background"
)
#pdf("./case_study_tortoises/figures/track.pdf", width=4, height=3, bg = "transparent")
# plot track with state colouring (note that some points overlap, despite transparency)
ggmap(map_bg) +
  geom_point(aes(x = x, y = y), 
             data = df.animal[which(mod.rtmb$states == 3), ], 
             colour = alpha(colour[3], 0.8), size = 0.02) +
  geom_point(aes(x = x, y = y), 
             data = df.animal[which(mod.rtmb$states == 2), ], 
             colour = alpha(colour[2], 0.5), size = 0.02) +
  geom_point(aes(x = x, y = y), 
             data = df.animal[which(mod.rtmb$states == 1), ], 
             colour = alpha(colour[1], 0.5), size = 0.02) +
  labs( x = "longitude", 
        y = "latitude")+
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title=element_text(size=14))
#dev.off()

# covariate time series
df.animal <- df.animal %>%
  arrange(timestamp) %>%
  mutate(
    # compute time difference between consecutive timestamps (in hours)
    time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "hours")),
    # start a new group if gap > 24 observations (or 24 hours if regularly sampled)
    gap_group = cumsum(ifelse(is.na(time_diff) | time_diff > 24, 1, 0))
  )

# zoomd-in view on temperature time series 
par(mfrow=c(1,1), mar=c(4,4,1,1), cex.lab=1.2, mgp = c(1.8, 0.5, 0) )
ggplot(df.animal[50:350,], aes(x = timestamp, y = temperature, group = gap_group)) +
  geom_line(color = "black") +
  labs(title = "", x = "time", y = "temperature (°C)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title=element_text(size=12))

#pdf("./case_study_tortoises/figures/temp_ts.pdf", width=6, height=4)
par(mfrow=c(1,1), mar=c(4,4,1,1), cex.lab=1.2, mgp = c(1.8, 0.5, 0) )
ggplot(df.animal, aes(x = timestamp, y = temperature, group = gap_group)) +
  geom_line(color = "black") +
  labs(title = "", x = "time", y = "temperature (°C)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title=element_text(size=12))
#dev.off()

# observation processes decoded
#pdf("case_study_tortoises/figures/obs_processes.pdf", width=8, height=6)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1), cex.lab = 1.2, mgp = c(1.8, 0.5, 0) )
plot(df.animal$step[1:500], type="h", col=colour[mod.rtmb$states], ylab="step length", xlab="time", main="", bty = "n")
plot(df.animal$angle[1:500], type="h", col=colour[mod.rtmb$states], ylab="turning angle", xlab="time", main="", bty = "n")
#dev.off()

# plotting the transition probabilities
vis_z <- seq(0, 45, length=100)
Gamma_vis <- tpm_g(vis_z, beta.hat)
#pdf("./case_study_tortoises/figures/transprobs.pdf", width=8, height=5)
par(mfrow=c(3,3), mar=c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.lab=1.3)
for(i in 1:3){
  for(j in 1:3){
    ylab_expr <-bquote(gamma[.(as.numeric(paste0(i, j)))])
    xlab_here <- if(i == 3) "temperature" else ""
    plot(vis_z, Gamma_vis[i,j,],
         type="l",
         lwd=2,
         col=colour[j],
         xlab=xlab_here,
         ylab=ylab_expr,
         main=paste0("Pr(", i, "\u2192", j, ")"),
         ylim=c(0,1),
         xlim=c(0,45),
         bty="n",
         cex.lab=1.5,
         cex.axis=1.2,
         cex.main=1.5,
         xaxt="n")   # remove default x axis
    
    axis(1, at = c(0,10,20,30,40), labels = c(0,10,20,30,40))
    
  }
}
#dev.off()

# state dependent distributions
zseq = seq(min(df.animal$temperature), max(df.animal$temperature), by = 0.01)
Gamma_seq = tpm_g(Z = cbind(zseq), beta.hat)
Prob = matrix(nrow = length(zseq), ncol = 3)
for(i in 1:length(zseq)){ Prob[i,] = LaMa::stationary(Gamma_seq[,,i]) }
prob = apply(Prob, 2, mean)

#pdf("./case_study_tortoises/figures/obs_dist.pdf", width=4.5, height=3)
par(mfrow = c(1, 2), mar = c(4, 3, 2, 1), cex.lab = 1, mgp = c(1.8, 0.5, 0))
hist(df.animal$step, breaks = 60, col = "lightgrey", xlab = "step length (in metres)", ylab="density",  main = "", prob = T, border=0)
for(j in 1:3) curve(prob[j] * dgamma2(x, mu.hat[j], sigma.hat[j]), 
                    lwd = 2, add = T, col = colour[j])
curve(prob[1]*dgamma2(x, mu.hat[1], sigma.hat[1]) + prob[2]*dgamma2(x, mu.hat[2], sigma.hat[2]) + 
        prob[3]*dgamma2(x, mu.hat[3], sigma.hat[3]), lwd = 2, lty = 2, add = T, n=500)

mu.kappa = c(pi, pi, 0)
legend("topright", col=c(colour, "black"), legend=expression(state~1, state~2, state~3, marginal), lwd=2, lty=c(1, 1, 1, 2), bty="n", cex=1)
hist(df.animal$angle, breaks = 60, col = "lightgrey", xlab = "turning angle", ylab="", main = "", prob = T, border=0)
for(j in 1:3) curve(prob[j] * dvm(x, mu.kappa[j], kappa.hat[j]), 
                    lwd = 2, add = T, col = colour[j])
curve(prob[1]*dvm(x, mu.kappa[1], kappa.hat[1]) + prob[2]*dvm(x, mu.kappa[2], kappa.hat[2]) + prob[3]*dvm(x, mu.kappa[3], kappa.hat[3]), 
      lwd = 2, lty = 2, add = T, n=500)
#dev.off()

# hypothetical stationary distribution
n <- dim(df.animal)[1]
Delta <- matrix(NA, n, N)
Delta[1,] <- rep(1/N, N)
for (t in 2:n) Delta[t,] <- Delta[t-1,] %*% Gamma.hat[,,t]

Delta <- Delta[10:n,] # we remove first values as initial distribution is arbitrary  
z <- df.animal$temperature[10:n]

zseq <- seq(min(z)-0.5, max(z)+0.5, length = 200)
Gammaseq <- tpm_g(zseq, beta.hat)
Deltaseq <- matrix(NA, length(zseq), N)
for (t in 1:length(zseq)) Deltaseq[t,] <- LaMa::stationary(Gammaseq[,,t])

#pdf("./case_study_tortoises/figures/hypothetical.pdf", width=6, height=4)
par(mfrow = c(1,1), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3) 
plot(z, Delta[, 1],
     xlim = c(min(z), max(z)), ylim = c(0, 1),
     pch = 16, col = alpha(colour[1], 0),
     bty = "n",
     ylab = "Pr(state i | temp), i=1,2,3",
     xlab = "temperature",
     main = "")
for (state in 1:N) {
  points(jitter(z, factor = 0), Delta[, state], col = alpha(colour[state], 0.1), pch = 16)
}
for (state in 1:N) {
  darker_col <- adjustcolor(colour[state], red.f = 0.9, green.f = 0.9, blue.f = 0.9, alpha.f = 1)
  lines(zseq, Deltaseq[, state], col = "white", lwd = 7, lty = 1)
  lines(zseq, Deltaseq[, state], col = darker_col, lwd = 3, lty = 2)
}

legend("topright",col = c( "orange",  "skyblue", "seagreen"), lty = c(3,3,3),   
       lwd = 3, bty = "n", legend = expression("", rho, ""), y.intersp = 0.3, cex = 1.3)
#dev.off()


# BB approach ----

df.animal$doy <- as.numeric(format(df.animal$timestamp, "%j"))
mod_lm <- lm(temperature ~ sin((2*pi*doy)/365) + cos((2*pi*doy)/365), data=df.animal)

simulated_x <- mclapply(1:100, function(i) block_bootstrap(mod_lm$residuals, 24, n), mc.cores = max(1, detectCores() - 2))
simulated_x <- do.call(cbind, simulated_x)
simulated_x <- simulated_x + mod_lm$fitted.values

sim_delta <- future_lapply(1:100,
  function(i) compStateProbs(simulated_x[, i], mod.rtmb$beta, n=length(simulated_x[,1])))
sim_delta <- array(unlist(sim_delta), dim = c(n, N, 100))

x_bins <- seq(min(simulated_x)-1, max(simulated_x)+1, length.out = 30)
bin_midpoints <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2

mean_state_probs <- matrix(NA, nrow = length(bin_midpoints), ncol = N)
for (state in 1:N) {
  mean_state_probs[, state] <- sapply(1:(length(x_bins) - 1), function(b) {
    mean(sim_delta[, state, ][simulated_x >= x_bins[b] & simulated_x < x_bins[b + 1]], na.rm = TRUE)
  })
}
mean_state_probs <- na.approx(mean_state_probs)

#pdf("./case_study_tortoises/figures/bb_approach.pdf", width=6, height=4)
par(mfrow = c(1,1), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3) 
plot(z, Delta[, 1],
     xlim = c(min(z),max(z)), ylim = c(0, 1),
     pch = 16, col = alpha(colour[1], 0),
     bty = "n", ylab = "Pr(state i | temp), i=1,2,3",
     xlab = "temperature", main="BB resampling")
for (state in 1:N) {
  points(jitter(z, factor = 0), Delta[, state], col = alpha(colour[state], 0.1), pch = 16)
}
for (state in 1:N) {
  darker_col <- adjustcolor(colour[state], red.f = 0.9, green.f = 0.9, blue.f = 0.9, alpha.f = 1)
  lines(ksmooth(bin_midpoints, mean_state_probs[,state], "normal", bandwidth = 5), col = "white", lwd = 7, lty = 1)
  lines(ksmooth(bin_midpoints, mean_state_probs[,state], "normal", bandwidth = 5), 
        col = darker_col, lwd = 3)
}
#legend("topright",col = c(colour, "transparent", "black"),
#       pch = c(16, 16, 16, NA, NA), lty = c(NA, NA, NA, NA, 1),   
#       lwd = 3, bty = "n", legend = expression(state~1, state~2, state~3, "", delta))

#dev.off()

# Dirichlet regression ----

system.time(
  mod <- dir_reg(Delta, z, k = 8, bs="tp", lambda0=100)
)

summary(mod)
x_p <- seq(min(z)-0.5, max(z)+0.5, length = 200)
Mean <- mod$predict(x_p)

par(mfrow = c(1,1), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3) 

#pdf("./case_study_tortoises/figures/dirichlet_regression.pdf", width=6, height=4)
plot(0,0,xlim=c(min(df.animal$temperature), max(df.animal$temperature)), ylim = c(0, 1),
     pch = 16, col = alpha(colour[1], 0),
     bty = "n", ylab = "Pr(state i | temp), i=1,2,3",
     xlab = "temperature", main="Flexible Dirichlet regression")  
for (state in 1:N) {
  points(jitter(z, factor = 0), Delta[, state], col = alpha(colour[state], 0.1), pch = 16)
}
for (state in 1:N) {
  darker_col <- adjustcolor(colour[state], red.f = 0.9, green.f = 0.9, blue.f = 0.9, alpha.f = 1)
  lines(x_p, Mean[,state], col = "white", lwd = 7, lty = 1)
  lines(x_p, Mean[,state], col = darker_col, lwd = 3, lty = 1)
}
#dev.off()

# or by using this function that generates the plot automatically
apply_dir_reg(y=Delta,x=z,N=3, covname="temperature")

# CI BB approach ----
nCI <- 1000
num_covsim <- 50

sds <- solve(obj.rtmb$he()[10:21, 10:21])
betas <- rmvnorm(nCI, c(mod.rtmb$beta), sds)

all_mean_probs_II <- vector("list", length = nCI)

mod_lm <- lm(temperature ~ sin((2*pi*doy)/365) + cos((2*pi*doy)/365), data=df.animal)

for (i in 1:nCI) {
  
  simulated_x <- mclapply(1:num_covsim, function(i) block_bootstrap(mod_lm$residuals, 24, n), mc.cores = max(1, detectCores() - 2))
  simulated_x <- do.call(cbind, simulated_x)
  simulated_x <- simulated_x + mod_lm$fitted.values
  
  simbeta <- matrix(betas[i,], nrow=6)
  
  sim_delta <- future_lapply(
    1:num_covsim,
    function(i) compStateProbs(simulated_x[, i], simbeta, n=length(simulated_x[,1])))  
  sim_delta <- array(unlist(sim_delta), dim = c(length(simulated_x[,1]), N, 50))
  
  x_bins <- seq(min(simulated_x), max(simulated_x)+0.1, length.out = 30)
  bin_midpoints <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
  mean_state_probs <- matrix(NA, nrow = length(bin_midpoints), ncol = N)
  for (state in 1:N) {
    mean_state_probs[, state] <- sapply(1:(length(x_bins) - 1), function(b) {
      mean(sim_delta[, state, ][simulated_x >= x_bins[b] & simulated_x < x_bins[b + 1]], na.rm = TRUE)
    })
  }
  mean_state_probs <- na.approx(mean_state_probs)
  
  all_mean_probs_II[[i]] <- mean_state_probs
  
}

all_probs_array_II <- array(NA, dim = c(nCI,29, 3))
for (i in 1:nCI) {
  all_probs_array_II[i, , ] <- all_mean_probs_II[[i]]
}

par(mfrow = c(1,1), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3) 
#pdf("./case_study_tortoises/figures/bb_approach_ci.pdf", width=6, height=4)
plot(NULL, xlim = c(min(df.animal$temperature), max(df.animal$temperature)), ylim = c(0, 1), bty = "n",
     xlab = "temperature", ylab = "Pr(state i | temp), i=1,2,3", main = "BB resampling")

for (state in 1:3) {
  y <- all_probs_array_II[,,state]
  yCI <- apply(y, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  lines(ksmooth(bin_midpoints, colMeans(y, na.rm = TRUE),"normal", bandwidth=5), col = colour[state], lwd = 3)
  polygon(c(bin_midpoints, rev(bin_midpoints)), 
          c(yCI[1, ], rev(yCI[2, ])), col = alpha(colour[state], 0.2), border = NA)
}
legend("topright", legend = paste("state", 1:3), col = colour, lwd = 3, bty = "n")
#dev.off()

# CI Dirichlet regression ----

nCI <- 1000
N <- 3
n <- nrow(df.animal)

sds <- solve(obj.rtmb$he()[10:21, 10:21])
betas <- MASS::mvrnorm(nCI, mu = c(mod.rtmb$beta), Sigma = sds)

z <- df.animal$temperature[10:n]
zseq <- seq(min(z), max(z), length = 200)

Delta_boot <- array(NA, dim = c(length(zseq), N, nCI))

run_boot <- function(b) {
  beta_b <- betas[b, ]
  Gamma_b <- tpm_g(matrix(df.animal$temperature), matrix(beta_b, nrow=6))
  Delta_b <- matrix(NA, n, N)
  Delta_b[1, ] <- rep(1/N, N)
  for (t in 2:n) Delta_b[t, ] <- Delta_b[t-1, ] %*% Gamma_b[, , t]
  Delta_b <- Delta_b[10:n, , drop = FALSE]
  
  mod <- dir_reg(Delta_b, z, k = 8, bs="tp", lambda0=100, silent=2)
  x_p <- seq(min(z), max(z), length = 200)
  Mean <- mod$predict(x_p)
  Mean
}

Delta_list <- lapply(1:nCI, run_boot)


for (b in 1:nCI) {
  Delta_boot[, , b] <- Delta_list[[b]]
}

Delta_median <- apply(Delta_boot, c(1, 2), median, na.rm = TRUE)
Delta_low <- apply(Delta_boot, c(1, 2), function(x) quantile(x, 0.025, na.rm = TRUE))
Delta_high <- apply(Delta_boot, c(1, 2), function(x) quantile(x, 0.975, na.rm = TRUE))

pdf("./case_study_tortoises/figures/dirichlet_regression_ci.pdf", width=6, height=4)
par(mfrow = c(1,1), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3) 
plot(zseq, Delta_median[, 1], type = "n",
     ylim = c(0, 1), bty="n", xlim=c(min(df.animal$temperature), max(df.animal$temperature)),
     xlab = "temperature", ylab = "Pr(state i | temp), i=1,2,3", main="Flexible Dirichlet regression")  

for (state in 1:N) {
  col_state <- adjustcolor(colour[state], red.f = 0.9, green.f = 0.9, blue.f = 0.9, alpha.f = 1)
  polygon(c(zseq, rev(zseq)),
          c(Delta_low[, state], rev(Delta_high[, state])),
          col = adjustcolor(col_state, alpha.f = 0.2), border = NA)
  lines(x_p, Mean[,state], col = colour[state], lwd = 3, lty = 1)
}
#dev.off()