## Supplemenary material: Plot 3 Settings x 4 Approaches 
# n = 2000

pdf("simulations/figures/all_2000.pdf")
par(
  mfcol = c(4, 3),
  mar = c(2.5, 2.5, 1, 0.5),
  oma = c(0, 4, 0, 0),   # extra outer margin on the left
  mgp = c(1.2, 0.3, 0),
  tcl = -0.2
)
### Setting I

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)",
     main = "Setting (I)", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve1$State1[[i]]$x, curve1$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve1$State2[[i]]$x, curve1$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve1$State3[[i]]$x, curve1$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

mtext("AR approach", side = 2, line = 4, las = 0, cex=0.9)

## BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)", main = "", bty = "n",
     xlim=c(min(zGT1)+1, max(zGT1)-1))

for (i in 1:num_simulations) {
  lines(curve_1BB$State1[[i]]$x, curve_1BB$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve_1BB$State2[[i]]$x, curve_1BB$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve_1BB$State3[[i]]$x, curve_1BB$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

mtext("BB approach", side = 2, line = 4, las = 0, cex=0.9)

# Dirichlet

plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "Pr(state | z)", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve1$State1[[i]]$x, gamcurve1$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve1$State2[[i]]$x, gamcurve1$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve1$State3[[i]]$x, gamcurve1$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(cut1_1$x, cut1_1$y, col = colour[1], lwd = 3)
lines(cut2_1$x, cut2_1$y, col = colour[2], lwd = 3)
lines(cut3_1$x, cut3_1$y, col = colour[3], lwd = 3)

mtext("Dirichlet regression", side = 2, line = 4, las = 0, cex=0.9)

# Hypothetical 
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state | z)",
     main = "", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(hypothetical$State1[[i]]$x, hypothetical$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical$State2[[i]]$x, hypothetical$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical$State3[[i]]$x, hypothetical$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

mtext("Hypothetical \n stationary distribution", side = 2, line = 4, las = 0, cex=0.9)

## Setting II

# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main="Setting (II)", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve2$State1[[i]]$x, curve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve2$State2[[i]]$x, curve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve2$State3[[i]]$x, curve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

# BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(min(zGT2)+1, max(zGT2)-1))

for (i in 1:num_simulations) {
  lines(curve_2BB$State1[[i]]$x, curve_2BB$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve_2BB$State2[[i]]$x, curve_2BB$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve_2BB$State3[[i]]$x, curve_2BB$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve2$State1[[i]]$x, gamcurve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve2$State2[[i]]$x, gamcurve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve2$State3[[i]]$x, gamcurve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)


# Hypothetical
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "",
     main = "", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(hypothetical2$State1[[i]]$x, hypothetical2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical2$State2[[i]]$x, hypothetical2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical2$State3[[i]]$x, hypothetical2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)


# Setting III
# AR
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "",
     main = "Setting (III)", bty = "n",
     xlim=c(min(zGT3)+1, max(zGT3)-1))

for (i in 1:num_simulations) {
  lines(curve3$State1[[i]]$x, curve3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve3$State2[[i]]$x, curve3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve3$State3[[i]]$x, curve3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(cut1_3$x, cut1_3$y, col = colour[1], lwd = 3)
lines(cut2_3$x, cut2_3$y, col = colour[2], lwd = 3)
lines(cut3_3$x, cut3_3$y, col = colour[3], lwd = 3)

# BB
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n",
     xlim=c(min(zGT3)+1, max(zGT3)-1))

for (i in 1:num_simulations) {
  lines(curve_3BB$State1[[i]]$x, curve_3BB$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve_3BB$State2[[i]]$x, curve_3BB$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve_3BB$State3[[i]]$x, curve_3BB$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(cut1_3$x, cut1_3$y, col = colour[1], lwd = 3)
lines(cut2_3$x, cut2_3$y, col = colour[2], lwd = 3)
lines(cut3_3$x, cut3_3$y, col = colour[3], lwd = 3)

# Dirichlet
plot(NULL, ylim = c(0, 1),
     xlab = "", ylab = "", main = "", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:num_simulations) {
  lines(gamcurve3$State1[[i]]$x, gamcurve3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve3$State2[[i]]$x, gamcurve3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve3$State3[[i]]$x, gamcurve3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(cut1_3$x, cut1_3$y, col = colour[1], lwd = 3)
lines(cut2_3$x, cut2_3$y, col = colour[2], lwd = 3)
lines(cut3_3$x, cut3_3$y, col = colour[3], lwd = 3)


# Hypothetical
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "",
     main = "", bty = "n",
     xlim = c(-10, 10))

for (i in 1:num_simulations) {
  lines(hypothetical3$State1[[i]]$x, hypothetical3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical3$State2[[i]]$x, hypothetical3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical3$State3[[i]]$x, hypothetical3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
dev.off()
