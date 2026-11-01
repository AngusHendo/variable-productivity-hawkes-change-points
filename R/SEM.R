SEM <- function(
  theta_0,
  realiz,
  window,
  n_iter = 5,
  optim_maxit = 10,
  tol = 0.001,
  n_candidates = 100,
  trace = 0,
  eps_step = 0.1
) {
  theta <- theta_0
  prev_theta <- theta_0
  t_candidates <- sort(runif(n_candidates, window$t[1], window$t[2]))
  sigma <- (window$t[2] - window$t[1]) / 4
  for (i in 1:n_iter) {
    fit <- optim(
      par = theta,
      fn = augmented_loglik,
      method = "Nelder-Mead",
      control = list(fnscale = -1, trace = trace, maxit = optim_maxit),
      t_candidates = t_candidates,
      realiz = realiz,
      window = window
    )
    prev_theta <- theta
    theta <- fit$par
    t_star <- get_best_t_star(theta, t_candidates, realiz, window)
    t_candidates <- sort(runif(
      n_candidates,
      max(t_star - sigma, window$t[1]),
      min(t_star + sigma, window$t[2])
    ))
    sigma <- sigma - eps_step
    cat("Iteration", i, "\tt_star:", t_star, "\ttheta:", theta, "\n\n")
    if (max(abs(theta - prev_theta)) < tol) break
  }
  return(list(theta = theta, t_star = t_star, sigma = sigma * 2))
}

augmented_loglik <- function(
  theta,
  t_candidates,
  realiz,
  window
) {
  theta <- unlist(theta)
  stationary_params <- list(
    mu = theta[1],
    alpha = theta[2],
    beta = theta[3]
  )
  LL_cands <- rep(0, length(t_candidates))
  for (i in 1:length(t_candidates)) {
    LL_cands[i] <- loglik_hawk_K_fn(
      K_params = theta[c(4, 5)],
      K_fn = K_FN_discrete,
      K_fn_consts = c(t_candidates[i], window$t[2]),
      stationary_params = stationary_params,
      realiz = realiz,
      window = window
    )
  }
  w <- LL_weights(LL_cands)
  q <- 1 / length(LL_cands) * sum(w * LL_cands)
  return(q)
}

LL_weights <- function(LL_cands) {
  adjust <- exp(LL_cands / 1000)
  w <- adjust / sum(adjust)
  return(w)
}

get_best_t_star <- function(
  theta,
  t_candidates,
  realiz,
  window
) {
  LL_cands <- rep(0, length(t_candidates))
  stationary_params <- list(
    mu = theta[1],
    alpha = theta[2],
    beta = theta[3]
  )
  for (i in 1:length(t_candidates)) {
    LL_cands[i] <- loglik_hawk_K_fn(
      K_params = theta[c(4, 5)],
      K_fn = K_FN_discrete,
      K_fn_consts = c(t_candidates[i], window$t[2]),
      stationary_params = stationary_params,
      realiz = realiz,
      window = window
    )
  }
  w <- LL_weights(LL_cands)
  t_star <- median(t_candidates[which(w == max(w))])
  return(t_star)
}
