recovery <- read.csv("results/nesi_recovery_summary.csv")

png("results/nesi_recovery_summary.png", width = 1200, height = 750, res = 160)
plot(
  recovery$kappa_difference,
  recovery$mean_absolute_change_time_error,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#2166AC",
  xlab = expression(abs(kappa[2] - kappa[1])),
  ylab = "Mean absolute change-time error",
  main = "Change-time recovery across 4,382 NeSI simulations"
)
grid()
lines(
  recovery$kappa_difference,
  recovery$mean_absolute_change_time_error,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#2166AC"
)
dev.off()
