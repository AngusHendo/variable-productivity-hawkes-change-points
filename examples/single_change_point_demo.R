source("R/VariableProductivityHawkes.R")

set.seed(124)
TIME <- 100
OMEGA <- c(0, 10, 0, 10)
WINDOW <- list(t = c(0, TIME), s = as.owin(OMEGA))
SPLITS <- c(0.5 * TIME, TIME)
PARAMS <- list(mu = 0.75, alpha = 2, beta = 1)
K_PARAMS <- list(K1 = 0.1, K2 = 0.8)

onerealiz <- sim_hawkes_K_fn(
  stationary_params = PARAMS,
  K_params = K_PARAMS,
  K_fn_consts = SPLITS,
  K_fn = K_FN_discrete,
  window = WINDOW
)

t_candidates <- seq(25, 75, by = 2.5)
profile <- data.frame(
  t_star = t_candidates,
  loglik = NA,
  K1 = NA,
  K2 = NA
)

for (i in 1:length(t_candidates)) {
  fit <- fit_hawkes_K_fn(
    K_inits = c(0.4, 0.4),
    K_fn = K_FN_discrete,
    K_fn_consts = c(t_candidates[i], TIME),
    stationary_params = PARAMS,
    realiz = onerealiz,
    window = WINDOW,
    maxit = 200
  )
  profile$loglik[i] <- fit$value
  profile$K1[i] <- fit$par[1]
  profile$K2[i] <- fit$par[2]
}

best <- profile[which.max(profile$loglik), ]
profile$relative_loglik <- profile$loglik - max(profile$loglik)

write.csv(profile, "results/single_change_point_profile.csv", row.names = FALSE)

png("results/single_change_point_demo.png", width = 1400, height = 1000, res = 160)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(
  c(0, SPLITS[1], TIME),
  c(K_PARAMS$K1, K_PARAMS$K2, K_PARAMS$K2),
  type = "s",
  lwd = 3,
  col = "black",
  ylim = c(0, 1),
  xlab = "Time",
  ylab = expression(kappa(t)),
  main = "Simulated and fitted productivity"
)
lines(
  c(0, best$t_star, TIME),
  c(best$K1, best$K2, best$K2),
  type = "s",
  lwd = 3,
  lty = 2,
  col = "#B2182B"
)
rug(onerealiz$t, col = "#2166AC")
legend(
  "topleft",
  legend = c("True productivity", "Profile fit", "Events"),
  col = c("black", "#B2182B", "#2166AC"),
  lty = c(1, 2, 1),
  lwd = c(3, 3, 1),
  bty = "n"
)

plot(
  profile$t_star,
  profile$relative_loglik,
  type = "b",
  pch = 19,
  col = "#2166AC",
  xlab = "Candidate change time",
  ylab = "Relative profile log likelihood",
  main = "Change-point likelihood profile"
)
abline(v = SPLITS[1], lwd = 2)
abline(v = best$t_star, lwd = 2, lty = 2, col = "#B2182B")
legend(
  "bottomleft",
  legend = c("True change time", "Fitted change time"),
  col = c("black", "#B2182B"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)
dev.off()

print(best)
