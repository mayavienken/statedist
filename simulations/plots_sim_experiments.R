## Plots related to the simulation experiments that can be found in our paper

## After running the three scripts: 
# "sim_setting_I.R", "sim_setting_II.R" and "sim_setting_III.R"

# 1) Comparing AR resampling and hypothetical stationary distribution ----
# for Setting I and II (covariates follow AR(1) processes)

#pdf("./simulations/figures/ar_hypothetical_I_II.pdf", width=10, height=6)
par(mfrow=c(2,2), cex.main=1.2, cex.lab=1.2, mar=c(4,3,1,1), mgp = c(1.8, 0.5, 0) )
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main="Setting (I): AR resampling", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve2$State1[[i]]$x, curve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve2$State2[[i]]$x, curve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve2$State3[[i]]$x, curve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
legend("topleft",col = colour, lwd = 3,bty = "n",legend = expression(state~1, state~2, state~3), cex=1.3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "Setting (II): AR resampling", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(curve1$State1[[i]]$x, curve1$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve1$State2[[i]]$x, curve1$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve1$State3[[i]]$x, curve1$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "Setting (I): hypothetical stationary distribution", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(hypothetical2$State1[[i]]$x, hypothetical2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical2$State2[[i]]$x, hypothetical2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical2$State3[[i]]$x, hypothetical2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT2, mean_stateprobsGT2[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3",
     main = "Setting (II): hypothetical stationary distribution", bty = "n",
     xlim = c(-4, 4))

for (i in 1:num_simulations) {
  lines(hypothetical$State1[[i]]$x, hypothetical$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(hypothetical$State2[[i]]$x, hypothetical$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(hypothetical$State3[[i]]$x, hypothetical$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}

lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT1, mean_stateprobsGT1[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)

#dev.off()


# 2) Comparing AR and BB resampling for Setting III (periodic covariate) ----

#pdf("./simulations/figures/bb_ar_III.pdf", width = 10, height = 3.5)
par(mfrow=c(1,2), cex.main=1.2, cex.lab=1.2, mar=c(4,3,1,1), mgp = c(1.8, 0.5, 0) )
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", bty = "n",
     xlim=c(min(zGT3)+1, max(zGT3)-1), main="Setting (III): AR resampling")

for (i in 1:num_simulations) {
  lines(curve3AR$State1[[i]]$x, curve3AR$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve3AR$State2[[i]]$x, curve3AR$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve3AR$State3[[i]]$x, curve3AR$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)


plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (III): BB resampling", bty = "n",
     xlim=c(min(zGT3)+1, max(zGT3)-1))

for (i in 1:num_simulations) {
  lines(curve_3BB$State1[[i]]$x, curve_3BB$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(curve_3BB$State2[[i]]$x, curve_3BB$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(curve_3BB$State3[[i]]$x, curve_3BB$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 1], "normal", bandwidth = 1), col = colour[1], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 2], "normal", bandwidth = 1), col = colour[2], lwd = 3)
lines(ksmooth(bin_midpointsGT3, mean_stateprobsGT3[, 3], "normal", bandwidth = 1), col = colour[3], lwd = 3)
legend("right", legend = c("state 1", "state 2", "state 3"),
       col = c(colour), lwd = 2, bty = "n")

#dev.off()


# 3) Dirichlet Regression for all three settings ----

#pdf("./simulations/figures/bb_ar_III.pdf", width=9, height=3)
par(mfrow=c(1,3), mgp = c(1.8, 0.5, 0), mar=c(3, 3, 1,1), cex.lab=1.3, cex.main=1.3)
# Setting I
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (I)", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve1$State1[[i]]$x, gamcurve1$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve1$State2[[i]]$x, gamcurve1$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve1$State3[[i]]$x, gamcurve1$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(cut1_1$x, cut1_1$y, col = colour[1], lwd = 3)
lines(cut2_1$x, cut2_1$y, col = colour[2], lwd = 3)
lines(cut3_1$x, cut3_1$y, col = colour[3], lwd = 3)
legend("topleft",col = c(colour, "transparent"),lwd = 3,bty = "n",
       lty = c(1,1,1),legend = expression(state~1, state~2, state~3))

# Setting II
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (II)", bty = "n", 
     xlim=c(-4, 4))

for (i in 1:num_simulations) {
  lines(gamcurve2$State1[[i]]$x, gamcurve2$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve2$State2[[i]]$x, gamcurve2$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve2$State3[[i]]$x, gamcurve2$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(cut1$x, cut1$y, col = colour[1], lwd = 3)
lines(cut2$x, cut2$y, col = colour[2], lwd = 3)
lines(cut3$x, cut3$y, col = colour[3], lwd = 3)

# Setting III
plot(NULL, ylim = c(0, 1),
     xlab = "covariate value z", ylab = "Pr(state i | z), i=1,2,3", main = "Setting (III)", bty = "n", 
     xlim=c(-8, 8))

for (i in 1:num_simulations) {
  lines(gamcurve3$State1[[i]]$x, gamcurve3$State1[[i]]$y, col = alpha(colour[1], 0.1), lwd = 1)
  lines(gamcurve3$State2[[i]]$x, gamcurve3$State2[[i]]$y, col = alpha(colour[2], 0.1), lwd = 1)
  lines(gamcurve3$State3[[i]]$x, gamcurve3$State3[[i]]$y, col = alpha(colour[3], 0.1), lwd = 1)
}
lines(cut1_3$x, cut1_3$y, col = colour[1], lwd = 3)
lines(cut2_3$x, cut2_3$y, col = colour[2], lwd = 3)
lines(cut3_3$x, cut3_3$y, col = colour[3], lwd = 3)

#dev.off()
