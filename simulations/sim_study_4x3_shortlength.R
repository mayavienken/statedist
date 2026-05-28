## Plot results of the simulation study of four approaches x all three settings
## Hypothetical stationary distribution (Patterson et al., 2009), AR approach, 
## BB approach, Dirichlet regression
# n = 500 

# before running this script one must run
## sim_setting_I_short.R
## sim_setting_II_short.R
## sim_setting_III_short.R


# pdf("simulations/figures/all_500.pdf")

par(
  mfcol = c(4, 3),
  mar = c(2.5, 2.5, 1, 0.5),
  oma = c(0, 4, 0.2, 0),   # extra outer margin on the left
  mgp = c(1.2, 0.3, 0),
  tcl = -0.2
)


## Setting I -------------------------------------------------------------------

# Hypothetical 
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)",
     main = "Setting (I)", bty = "n",
     xlim = c(-4, 4))

for (i in 1: length(results1_s)) {
  lines(hypothetical_s$State1[[i]]$x, hypothetical_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(hypothetical_s$State2[[i]]$x, hypothetical_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(hypothetical_s$State3[[i]]$x, hypothetical_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_1_s$x, cut1_1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_1_s$x, cut2_1_s$y, col = colour_s[2], lwd = 3)
lines(cut3_1_s$x, cut3_1_s$y, col = colour_s[3], lwd = 3)

mtext("Hypothetical \n stationary distribution", side = 2, line = 4, las = 0, cex=0.9)

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-4, 4))

for (i in 1: length(results1_s)) {
  lines(curve1_s$State1[[i]]$x, curve1_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve1_s$State2[[i]]$x, curve1_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve1_s$State3[[i]]$x, curve1_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_1_s$x, cut1_1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_1_s$x, cut2_1_s$y, col = colour_s[2], lwd = 3)
lines(cut3_1_s$x, cut3_1_s$y, col = colour_s[3], lwd = 3)

mtext("AR approach", side = 2, line = 4, las = 0, cex=0.9)

## BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:length(results_1BB_s)) {
  lines(curve_1BB_s$State1[[i]]$x, curve_1BB_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_1BB_s$State2[[i]]$x, curve_1BB_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_1BB_s$State3[[i]]$x, curve_1BB_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_1_s$x, cut1_1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_1_s$x, cut2_1_s$y, col = colour_s[2], lwd = 3)
lines(cut3_1_s$x, cut3_1_s$y, col = colour_s[3], lwd = 3)

mtext("BB approach", side = 2, line = 4, las = 0, cex=0.9)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:length(gam_results1_s)) {
  lines(gamcurve1_s$State1[[i]]$x, gamcurve1_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(gamcurve1_s$State2[[i]]$x, gamcurve1_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(gamcurve1_s$State3[[i]]$x, gamcurve1_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_1_s$x, cut1_1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_1_s$x, cut2_1_s$y, col = colour_s[2], lwd = 3)
lines(cut3_1_s$x, cut3_1_s$y, col = colour_s[3], lwd = 3)

mtext("Dirichlet regression", side = 2, line = 4, las = 0, cex=0.9)


## Setting II ------------------------------------------------------------------

# Hypothetical
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "",
     main = "Setting (II)", bty = "n",
     xlim = c(-4, 4))

for (i in 1:length(results2_s)) {
  lines(hypothetical2_s$State1[[i]]$x, hypothetical2_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(hypothetical2_s$State2[[i]]$x, hypothetical2_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(hypothetical2_s$State3[[i]]$x, hypothetical2_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:length(results_s)) {
  lines(curve_s$State1[[i]]$x, curve_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_s$State2[[i]]$x, curve_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_s$State3[[i]]$x, curve_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)

# BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:length(results_2BB_s)) {
  lines(curve_2BB_s$State1[[i]]$x, curve_2BB_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_2BB_s$State2[[i]]$x, curve_2BB_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_2BB_s$State3[[i]]$x, curve_2BB_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:length(gam_results2_s)) {
  lines(gamcurve2_s$State1[[i]]$x, gamcurve2_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(gamcurve2_s$State2[[i]]$x, gamcurve2_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(gamcurve2_s$State3[[i]]$x, gamcurve2_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)


## Setting III -----------------------------------------------------------------

# Hypothetical
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "",
     main = "Setting (III)", bty = "n",
     xlim = c(-8, 8))

for (i in 1:length(results3_s)) {
  lines(hypothetical3_s$State1[[i]]$x, hypothetical3_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(hypothetical3_s$State2[[i]]$x, hypothetical3_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(hypothetical3_s$State3[[i]]$x, hypothetical3_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "",
     main = "", bty = "n",
     xlim=c(-8, 8))

for (i in 1:length(results3AR_s)) {
  lines(curve3AR_s$State1[[i]]$x, curve3AR_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve3AR_s$State2[[i]]$x, curve3AR_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve3AR_s$State3[[i]]$x, curve3AR_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)

# BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(-8, 8))

for (i in 1:length(results_3BB_s)) {
  lines(curve_3BB_s$State1[[i]]$x, curve_3BB_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(curve_3BB_s$State2[[i]]$x, curve_3BB_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(curve_3BB_s$State3[[i]]$x, curve_3BB_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "", main = "", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:length(gam_results3_s)) {
  lines(gamcurve3_s$State1[[i]]$x, gamcurve3_s$State1[[i]]$y, col = alpha(colour_s[1], 0.1), lwd = 1)
  lines(gamcurve3_s$State2[[i]]$x, gamcurve3_s$State2[[i]]$y, col = alpha(colour_s[2], 0.1), lwd = 1)
  lines(gamcurve3_s$State3[[i]]$x, gamcurve3_s$State3[[i]]$y, col = alpha(colour_s[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)


# dev.off()
