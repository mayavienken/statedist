## Plot results of the simulation study of four approaches x all three settings
## Hypothetical stationary distribution (Patterson et al., 2009), AR approach, 
## BB approach, Dirichlet regression
# n = 10,000 

# before running this script one must run
## sim_setting_I_long.R
## sim_setting_II_long.R
## sim_setting_III_long.R


# pdf("simulations/figures/all_10000.pdf")

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

for (i in 1:num_simulations_l) {
  lines(hypothetical_l$State1[[i]]$x, hypothetical_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(hypothetical_l$State2[[i]]$x, hypothetical_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(hypothetical_l$State3[[i]]$x, hypothetical_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_1_l$x, cut1_1_l$y, col = colour_l[1], lwd = 3)
lines(cut2_1_l$x, cut2_1_l$y, col = colour_l[2], lwd = 3)
lines(cut3_1_l$x, cut3_1_l$y, col = colour_l[3], lwd = 3)

mtext("Hypothetical \n stationary distribution", side = 2, line = 4, las = 0, cex=0.9)

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(curve1_l$State1[[i]]$x, curve1_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve1_l$State2[[i]]$x, curve1_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve1_l$State3[[i]]$x, curve1_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_1_l$x, cut1_1_l$y, col = colour_l[1], lwd = 3)
lines(cut2_1_l$x, cut2_1_l$y, col = colour_l[2], lwd = 3)
lines(cut3_1_l$x, cut3_1_l$y, col = colour_l[3], lwd = 3)
mtext("AR approach", side = 2, line = 4, las = 0, cex=0.9)

## BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(curve_1BB_l$State1[[i]]$x, curve_1BB_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve_1BB_l$State2[[i]]$x, curve_1BB_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve_1BB_l$State3[[i]]$x, curve_1BB_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_1_l$x, cut1_1_l$y, col = colour_l[1], lwd = 3)
lines(cut2_1_l$x, cut2_1_l$y, col = colour_l[2], lwd = 3)
lines(cut3_1_l$x, cut3_1_l$y, col = colour_l[3], lwd = 3)

mtext("BB approach", side = 2, line = 4, las = 0, cex=0.9)

# Dirichlet

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(gamcurve1_l$State1[[i]]$x, gamcurve1_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(gamcurve1_l$State2[[i]]$x, gamcurve1_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(gamcurve1_l$State3[[i]]$x, gamcurve1_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
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
for (i in 1:num_simulations_l) {
  lines(hypothetical2_l$State1[[i]]$x, hypothetical2_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(hypothetical2_l$State2[[i]]$x, hypothetical2_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(hypothetical2_l$State3[[i]]$x, hypothetical2_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(-4, 4))

for (i in 1:num_simulations_l) {
  lines(curve_l$State1[[i]]$x, curve_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve_l$State2[[i]]$x, curve_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve_l$State3[[i]]$x, curve_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)

# BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(-4, 4))
for (i in 1:num_simulations_l) {
  lines(curve_2BB_l$State1[[i]]$x, curve_2BB_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve_2BB_l$State2[[i]]$x, curve_2BB_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve_2BB_l$State3[[i]]$x, curve_2BB_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_s$x, cut1_s$y, col = colour_s[1], lwd = 3)
lines(cut2_s$x, cut2_s$y, col = colour_s[2], lwd = 3)
lines(cut3_s$x, cut3_s$y, col = colour_s[3], lwd = 3)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "", main = "", bty = "n", 
     xlim=c(-4, 4))
for (i in 1:num_simulations_l) {
  lines(gamcurve2_l$State1[[i]]$x, gamcurve2_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(gamcurve2_l$State2[[i]]$x, gamcurve2_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(gamcurve2_l$State3[[i]]$x, gamcurve2_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
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
for (i in 1:num_simulations_l) {
  lines(hypothetical3_l$State1[[i]]$x, hypothetical3_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(hypothetical3_l$State2[[i]]$x, hypothetical3_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(hypothetical3_l$State3[[i]]$x, hypothetical3_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "",
     main = "", bty = "n",
     xlim=c(-8, 8))
for (i in 1:num_simulations_l) {
  lines(curve3AR_l$State1[[i]]$x, curve3AR_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve3AR_l$State2[[i]]$x, curve3AR_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve3AR_l$State3[[i]]$x, curve3AR_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)

# BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(-8, 8))
for (i in 1:num_simulations_l) {
  lines(curve_3BB_l$State1[[i]]$x, curve_3BB_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(curve_3BB_l$State2[[i]]$x, curve_3BB_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(curve_3BB_l$State3[[i]]$x, curve_3BB_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "", main = "", bty = "n", 
     xlim=c(-8, 8))
for (i in 1:num_simulations_l) {
  lines(gamcurve3_l$State1[[i]]$x, gamcurve3_l$State1[[i]]$y, col = alpha(colour_l[1], 0.1), lwd = 1)
  lines(gamcurve3_l$State2[[i]]$x, gamcurve3_l$State2[[i]]$y, col = alpha(colour_l[2], 0.1), lwd = 1)
  lines(gamcurve3_l$State3[[i]]$x, gamcurve3_l$State3[[i]]$y, col = alpha(colour_l[3], 0.1), lwd = 1)
}
lines(cut1_3_s$x, cut1_3_s$y, col = colour_s[1], lwd = 3)
lines(cut2_3_s$x, cut2_3_s$y, col = colour_s[2], lwd = 3)
lines(cut3_3_s$x, cut3_3_s$y, col = colour_s[3], lwd = 3)


# dev.off()
