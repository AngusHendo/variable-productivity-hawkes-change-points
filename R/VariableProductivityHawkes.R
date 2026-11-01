library(spatstat)
library(data.table)
library(ggplot2)

###### METHODS FOR FITTING STATIONARY PARAMETERS (MU, ALPHA, BETA, K_OVERALL) #####

#' Calculates the spatio-temporal integral directly
#'
#' Computes the integral of the spatio-temporal intensity function over the given spatial and temporal bounds.
#'
#' @param realiz A dataframe containing `t`, `x`, and `y` values for each event.
#' @param bigT Numeric, the upper bound of the temporal window (lower bound is assumed to be 0).
#' @param bigX Numeric, the upper bound of the x-coordinate spatial window (lower bound is assumed to be 0).
#' @param bigY Numeric, the upper bound of the y-coordinate spatial window (lower bound is assumed to be 0).
#' @param alpha Numeric, the spatial decay parameter for the exponential kernel.
#' @param beta Numeric, the temporal decay parameter for the exponential kernel.
#' @return A numeric value representing the spatio-temporal integral.
get_sti <- function(realiz, bigT, bigX, bigY, alpha, beta) {
  erf_star <- function(S, x, alpha) {
    pnorm(S, mean = x, sd = 1 / sqrt(2 * alpha))
  }
  t_comp <- 1 - exp(-beta * (bigT - realiz$t))
  x_comp <- erf_star(bigX, realiz$x, alpha) - erf_star(0, realiz$x, alpha)
  y_comp <- erf_star(bigY, realiz$y, alpha) - erf_star(0, realiz$y, alpha)
  return(sum(t_comp * x_comp * y_comp))
}

#' Log-likelihood for a stationary Hawkes process
#'
#' Computes the log-likelihood of the data under a stationary Hawkes process with constant productivity (K).
#'
#' @param params A numeric vector containing the parameters \code{mu}, \code{alpha}, \code{beta}, and \code{K}.
#' @param realiz A dataframe of realized events with `t`, `x`, and `y` columns.
#' @param window A list specifying the temporal (`t`) and spatial (`s`) windows. `window$s` should be an `owin` object.
#' @param K_limit Numeric, an upper limit for \code{K} to ensure stability (default: \code{0.99999}).
#' @return The log-likelihood value as a numeric scalar.
loglik_stationary_hawkes<- function(
  params,
  realiz,
  window,
  K_limit = .99999
) {
  bigT <- window$t[2] - window$t[1]
  bigX <- window$s$xrange[2] - window$s$xrange[1]
  bigY <- window$s$yrange[2] - window$s$yrange[1]
  
  if (class(params) == "list") { params <- unlist(params) }
  
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  K <- params[4]
  
  # Restrict to non-negative parameters.
  if (min(mu, K, alpha, beta) < 0) return(-99999)
  
  # K < 1 to avoid explosive growth.
  if (K > K_limit) return(-99999)
  
  # Pre-compute distances
  x_diff <- outer(realiz$x, realiz$x, "-")
  y_diff <- outer(realiz$y, realiz$y, "-")
  t_diff <- outer(realiz$t, realiz$t, "-")
  space_dist <- x_diff^2 + y_diff^2
  time_dist <- t_diff
  exp_decay_mat <- exp(-beta * time_dist - alpha * space_dist)
  
  # For intlam to work according to Facilitated estimation of ETAS.
  mu_star <- mu / (spatstat.geom::area(window$s))
  
  # Calculate the large spatio-temporal integral directly
  spatio_temp_integral <- get_sti(realiz, bigT, bigX, bigY, alpha, beta)
  intlam <- mu * bigT + K * spatio_temp_integral
  
  # K is E # points, alpha and beta from exponential decay.
  const <- K * (alpha / pi) * beta
  
  # Initialize log sum:
  sum_log <- log(mu_star)
  
  for (j in 2:nrow(realiz)) {
    gij <- sum(exp(-beta * time_dist[j, 1:(j - 1)] - alpha * space_dist[j, 1:(j - 1)]))
    lamj <- mu_star + const * gij
    if (is.na(lamj) || lamj < 0) return(-99999)
    sum_log <- sum_log + log(lamj)
  }
  loglik <- sum_log - intlam
  if (loglik == -Inf) return(-99999)
  
  return(loglik)
}

#' Fits a stationary Hawkes process
#'
#' Optimizes the parameters of a stationary Hawkes process with constant productivity (K).
#'
#' @param params_init A numeric vector or list of initial parameter values (\code{mu}, \code{alpha}, \code{beta}, \code{K}).
#' @param realiz A dataframe of realized events.
#' @param window A list specifying the temporal (`t`) and spatial (`s`) windows.
#' @param trace Numeric, the trace level for optimization (default: \code{0}).
#' @param maxit Integer, the maximum number of iterations for optimization (default: \code{100}).
#' @return A list of fitted parameters (\code{mu}, \code{alpha}, \code{beta}, \code{K}).
fit_stationary_hawkes <- function(
  params_init,
  realiz,
  window,
  K_limit = .9999,
  trace = 0,
  maxit = 100
) {
  if (class(params_init) == "list") { params_init <- unlist(params_init) }
  fit <- optim(
    par = params_init,
    fn = loglik_stationary_hawkes,
    method = "Nelder-Mead",
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    realiz = realiz,
    window = window,
    K_limit = K_limit
  )
  
  return(as.list(fit$par))
}

##### ITERATIVE METHODS FOR FITTING PIECEWISE CONSTANT KAPPA #####

#' Calculates the spatio-temporal integral with varied K
#'
#' Computes the spatio-temporal integral for a Hawkes process with time-varying productivity \code{K}.
#'
#' @param realiz A dataframe of `t`, `x`, and `y` values for each event.
#' @param bigT Numeric, the upper bound of the temporal window.
#' @param bigX Numeric, the upper bound of the x-coordinate spatial window.
#' @param bigY Numeric, the upper bound of the y-coordinate spatial window.
#' @param alpha Numeric, the spatial decay parameter.
#' @param beta Numeric, the temporal decay parameter.
#' @param K_per_t A numeric vector of \code{K} values corresponding to time intervals.
#' @return The value of the spatio-temporal integral as a numeric scalar.
get_sti_K_per_t <- function(realiz, bigT, bigX, bigY, alpha, beta, K_per_t) {
  erf_star <- function(S, x, alpha) {
    pnorm(S, mean = x, sd = 1 / sqrt(2 * alpha))
  }
  t_comp <- 1 - exp(-beta * (bigT - realiz$t))
  x_comp <- erf_star(bigX, realiz$x, alpha) - erf_star(0, realiz$x, alpha)
  y_comp <- erf_star(bigY, realiz$y, alpha) - erf_star(0, realiz$y, alpha)
  return(sum(K_per_t * t_comp * x_comp * y_comp))
}

#' Log-likelihood for a Hawkes process with varied K
#'
#' Computes the log-likelihood of a Hawkes process with time-varying productivity \code{K}.
#'
#' @param params A numeric vector of \code{K} values to be optimized.
#' @param previous_K A numeric vector of previously fitted \code{K} values for earlier time intervals.
#' @param splits A numeric vector of time split points defining the intervals for \code{K}.
#' @param mu_hat Numeric, the background intensity parameter.
#' @param alpha Numeric, the spatial decay parameter.
#' @param beta Numeric, the temporal decay parameter.
#' @param realiz A dataframe of realized events.
#' @param window A list specifying the temporal (`t`) and spatial (`s`) windows.
#' @return The log-likelihood value as a numeric scalar.
loglik_hawkes_varied_K <- function(
  params,
  previous_K,
  splits,
  mu_hat,
  alpha,
  beta,
  realiz,
  window
) {
  bigT <- window$t[2] - window$t[1]
  bigX <- window$s$xrange[2] - window$s$xrange[1]
  bigY <- window$s$yrange[2] - window$s$yrange[1]
  mu <- mu_hat
  
  K <- params
  
  all_K <- c(previous_K, K)
  
  # Ensure change points are strictly increasing
  if (!(all(diff(splits) > 0))) {
    stop("Change points are not strictly increasing!")
    return(-99999)
  }
  
  t_intervals <- findInterval(realiz$t[], splits) + 1
  
  K_per_t <- all_K[t_intervals]
  
  # If there is 0/1 points in the first interval
  if (sum(realiz$t < min(splits)) < 2) {
    stop("There is less than 2 observations in the first interval!")
    return(-99999)
  }
  
  # Restrict to non-negative parameters.
  if (min(mu, all_K, alpha, beta) < 0) return(-99999)
  
  # Pre-compute distances
  x_diff <- outer(realiz$x, realiz$x, "-")
  y_diff <- outer(realiz$y, realiz$y, "-")
  t_diff <- outer(realiz$t, realiz$t, "-")
  space_dist <- x_diff^2 + y_diff^2
  time_dist <- t_diff
  
  # For intlam to work according to Facilitated estimation of ETAS.
  mu_star <- mu / (spatstat.geom::area(window$s))
  
  # Calculate the large spatio-temporal integral directly
  spatio_temp_integral <- get_sti_K_per_t(realiz, bigT, bigX, bigY, alpha, beta, K_per_t)
  intlam <- mu * bigT + spatio_temp_integral
  
  const <- (alpha / pi) * beta
  
  # Initialize log sum:
  sum_log <- log(mu_star)
  for (j in 2:nrow(realiz)) {
    gij <- sum(exp(-beta * time_dist[j, 1:(j - 1)] - alpha * space_dist[j, 1:(j - 1)]))
    lamj <- mu_star + K_per_t[j] * const * gij
    if (is.na(lamj) || lamj < 0) return(-99999)
    sum_log <- sum_log + log(lamj)
  }
  
  loglik <- sum_log - intlam
  if (loglik == -Inf) return(-99999)
  
  return(loglik)
}

#' Fits a Hawkes process with varied K
#'
#' Optimizes the \code{K} values for a Hawkes process with time-varying productivity.
#'
#' @param params_init Initial \code{K} values for optimization.
#' @param previous_K Previously fitted \code{K} values.
#' @param splits A numeric vector of time split points defining the intervals for \code{K}.
#' @param mu Numeric, the background intensity parameter.
#' @param alpha Numeric, the spatial decay parameter.
#' @param beta Numeric, the temporal decay parameter.
#' @param realiz A dataframe of realized events.
#' @param window A list specifying the temporal (`t`) and spatial (`s`) windows.
#' @param trace Numeric, the trace level for optimization.
#' @param maxit Integer, the maximum number of iterations for optimization.
#' @return A list containing the optimized \code{K} values.
fit_hawkes_varied_K <- function(params_init,
                                previous_K,
                                splits,
                                mu,
                                alpha,
                                beta,
                                realiz,
                                window,
                                trace = 0,
                                maxit = 1000,
                                K_max = 2) {
  fit <- optim(
    par = params_init,
    fn = loglik_hawkes_varied_K,
    method = "Brent",
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    previous_K = previous_K,
    splits = splits,
    mu_hat = mu,
    alpha = alpha,
    beta = beta,
    realiz = realiz,
    window = window,
    lower = 0,
    upper = K_max # This can be changed to model greater explosive behaviour.
  )
  
  return(fit)
}

#' Simulates a Hawkes process with variable K
#'
#' Generates events for a Hawkes process where the productivity parameter \code{K} varies with time.
#'
#' @param stationary_params A list containing the stationary parameters (\code{mu}, \code{alpha}, \code{beta}).
#' @param K_fn A function returning \code{K} values for given time points.
#' @param K_params A list of parameters for the \code{K_fn} function.
#' @param K_fn_consts Constants required by \code{K_fn}.
#' @param window A list specifying the temporal (`t`) and spatial (`s`) windows.
#' @return A dataframe of simulated events with columns `t`, `x`, and `y`.
fit_var_prod_hawkes <- function(
  inits,
  splits,
  realiz,
  window,
  trace = 0,
  maxit = 1000,
  output = TRUE
) {
  # Fit the overall mu, alpha, beta, K-hat values
  if (output) print("Fitting overall parameters...")
  
  stationary_params <- fit_hawkes_const_K(
    params_init = inits,
    realiz = realiz,
    window = window,
    trace = 0,
    maxit = 1000
  )
  
  if (output) {
    print("Overall fitted parameters:")
    print(params_fitted$par)
  }
  
  K_hat_list <- c()
  
  # Sequentially fit K_j for each interval, holding the previous K's constant.
  for (intrvl in 1:length(splits)) {
    if (output) cat("Fitting interval", intrvl, "...\n")
    K_hat_varied <- fit_hawkes_varied_K(
      params_init = stationary_params$K, 
      previous_K = K_hat_list,
      splits = splits,
      mu = stationary_params$mu,
      alpha = stationary_params$alpha,
      beta = stationary_params$beta,
      realiz = realiz[realiz$t < splits[intrvl],],
      window = list(t = c(0, splits[intrvl]), s = window$s),
      trace = trace,
      maxit = maxit
    )
    
    K_hat_list <- c(K_hat_list, K_hat_varied$par)
  }
  
  names(K_hat_list) <- paste0("K", 1:length(splits))
  output <- c(mu_hat, alpha_hat, beta_hat)
  names(output) <- c("mu", "alpha", "beta")
  output <- append(output, K_hat_list)
  return(output)
}


##### METHODS FOR FITTING KAPPA AS A FUNCTION OF TIME #####

#' Evaluate a Time-Dependent Function K(t) for Discrete Intervals
#'
#' Computes the value of \( K(t) \) based on predefined time splits and corresponding parameter values.
#'
#' @param t Numeric vector of times at which \( K(t) \) is evaluated.
#' @param K_params List or numeric vector of \( K \) values for each interval.
#' @param K_fn_consts Numeric vector defining the time splits (interval boundaries).
#' @return Numeric vector of \( K(t) \) values corresponding to the input times.
#' @examples
#' K_params <- list(1, 2, 3)
#' K_fn_consts <- c(0, 5, 10)
#' t <- c(2, 6, 11)
#' K_FN_discrete(t, K_params, K_fn_consts)
K_FN_discrete <- function(t, K_params, K_fn_consts) {
  splits <- K_fn_consts
  all_K <- unlist(K_params)
  t_intervals <- findInterval(t, splits) + 1
  K_per_t <- all_K[t_intervals]
  return(K_per_t)
}


#' Simulate a Spatio-Temporal Hawkes Process with Variable K
#'
#' Generates a realization of a Hawkes process where the productivity parameter \( K \) varies as a function of time.
#'
#' @param stationary_params List containing background intensity (\code{mu}), spatial decay (\code{alpha}), and temporal decay (\code{beta}).
#' @param K_fn Function that calculates \( K(t) \) based on the provided parameters and constants.
#' @param K_params Parameters for the \( K \) function, typically a list or numeric vector.
#' @param K_fn_consts Constants for the \( K \) function, such as time splits.
#' @param window List specifying the spatio-temporal domain. Should contain:
#'   - \code{t}: Numeric vector of time boundaries (\code{c(start, end)}).
#'   - \code{s}: Spatial window object (e.g., \code{spatstat.geom::owin}).
#' @param ... Additional arguments for \code{K_fn}.
#' @return A dataframe of simulated events with columns:
#'   - \code{t}: Event times.
#'   - \code{x}: Event x-coordinates.
#'   - \code{y}: Event y-coordinates.
#' @examples
#' stationary_params <- list(mu = 1, alpha = 2, beta = 1)
#' K_fn <- K_FN_discrete
#' K_params <- list(1, 2)
#' K_fn_consts <- c(0, 5)
#' window <- list(t = c(0, 10), s = spatstat.geom::owin(c(0, 1), c(0, 1)))
#' sim_hawkes_K_fn(stationary_params, K_fn, K_params, K_fn_consts, window)
sim_hawkes_K_fn <- function(
  stationary_params,
  K_fn,
  K_params,
  K_fn_consts,
  window,
  ...
){
  mu <- stationary_params$mu
  alpha <- stationary_params$alpha
  beta <- stationary_params$beta
  K_vec <- unlist(K_params)
  
  # # Ensure keys and values are of the same length
  # if (length(K_vec) != length(K_splits)) {
  #   cat("Length of K productivities:", length(K_vec), "\nLength of time splits:", length(K_splits), "\n")
  #   print(K_splits)
  #   stop("Number of K productivities must be the same as the number of time splits.")
  # }
  
  # Initialize the output list of events
  events <- list()
  events$t <- numeric()
  events$x <- numeric()
  events$y <- numeric()
  
  # Generate background events
  n_bg <- rpois(1, mu * (window$t[2] - window$t[1]))
  events$t <- sort(runif(n_bg, min=window$t[1], max=window$t[2]))
  locations = spatstat.random::runifpoint(n_bg, win = window$s)
  events$x <- locations$x
  events$y <- locations$y
  
  # Initialize the queue with background events
  event_queue <- data.frame(t = events$t, x = events$x, y = events$y)
  event_queue <- data.table::as.data.table(event_queue)
  
  # Initialize the list to store new events
  new_events_list <- vector("list", dim(event_queue)[1] * (1 / (max(1 - max(unlist(K_params)), 1e-5))) * 2)  # Pre-allocate space for new events
  list_index <- 1
  
  # Iterate through the new events, producing their offspring
  while (nrow(event_queue) > 0) {
    current_event <- event_queue[1, ]
    event_queue <- event_queue[-1, ]  # Remove the processed event
    
    # Get the corresponding value of K, based on the time split
    K_curr <- K_fn(
      t = current_event$t, 
      K_params = K_params, 
      K_fn_consts = K_fn_consts,
      ...
    )
    n_aftershocks <- rpois(1, K_curr)
    
    if (n_aftershocks > 0) {
      delay <- rexp(n_aftershocks, rate = beta)
      dist <- sqrt(rexp(n_aftershocks, rate = alpha))
      angle <- runif(n_aftershocks, min = 0, max = 2 * pi)
      new_x <- current_event$x + dist * cos(angle)
      new_y <- current_event$y + dist * sin(angle)
      new_t <- current_event$t + delay
      
      # Remove points outside window
      valid_events <- spatstat.geom::inside.owin(x = new_x, y = new_y, w = window$s) &
        new_t <= window$t[2]
      valid_events = which(valid_events)
      
      if (length(valid_events) != 0) {
        events$t <- c(events$t, new_t[valid_events])
        events$x <- c(events$x, new_x[valid_events])
        events$y <- c(events$y, new_y[valid_events])
        
        # Add the valid new events to the queue
        new_events <- data.table(t = new_t[valid_events],
                                 x = new_x[valid_events],
                                 y = new_y[valid_events])
        new_events_list[[list_index]] <- new_events
        list_index <- list_index + 1
      }
    }
    # Concatenate new events to event_queue only if we have only one event left to go
    if (dim(event_queue)[1] <= 1) {
      event_queue <- rbindlist(list(event_queue, rbindlist(new_events_list, use.names = TRUE)))
      new_events_list <- vector("list", length(new_events_list))  # Reset the list
      list_index <- 1
    }
  }
  events <- data.frame(events)
  return(events[order(events$t), ])
}


#' Log-Likelihood for a Spatio-Temporal Hawkes Process with Variable K
#'
#' Computes the log-likelihood of a given realization of a Hawkes process with time-varying \( K(t) \).
#'
#' @param K_params Parameters for the \( K \) function, typically a list or numeric vector.
#' @param K_fn Function that calculates \( K(t) \) based on the provided parameters and constants.
#' @param K_fn_consts Constants for the \( K \) function, such as time splits.
#' @param stationary_params List containing background intensity (\code{mu}), spatial decay (\code{alpha}), and temporal decay (\code{beta}).
#' @param realiz Dataframe of observed events with columns:
#'   - \code{t}: Event times.
#'   - \code{x}: Event x-coordinates.
#'   - \code{y}: Event y-coordinates.
#' @param window List specifying the spatio-temporal domain. Should contain:
#'   - \code{t}: Numeric vector of time boundaries (\code{c(start, end)}).
#'   - \code{s}: Spatial window object (e.g., \code{spatstat.geom::owin}).
#' @param ... Additional arguments for \code{K_fn}.
#' @return Numeric value of the log-likelihood for the observed data given the model parameters.
#' @examples
#' stationary_params <- list(mu = 1, alpha = 2, beta = 1)
#' realiz <- data.frame(t = c(1, 2), x = c(0.5, 0.6), y = c(0.4, 0.3))
#' window <- list(t = c(0, 10), s = spatstat.geom::owin(c(0, 1), c(0, 1)))
#' loglik_hawk_K_fn(list(1, 2), K_FN_discrete, c(0, 5), stationary_params, realiz, window)
loglik_hawk_K_fn <- function(
  K_params,
  K_fn,
  K_fn_consts,
  stationary_params,
  realiz,
  window,
  ...
) {
  mu <- stationary_params$mu
  alpha <- stationary_params$alpha
  beta <- stationary_params$beta
  
  bigT <- window$t[2] - window$t[1]
  bigX <- window$s$xrange[2] - window$s$xrange[1]
  bigY <- window$s$yrange[2] - window$s$yrange[1]
  
  # Ensure change points are strictly increasing
  # if (!(all(diff(K_fn_consts) > 0))) return(-99999)
  
  # Get K values for each time t in realization
  K_per_t <- K_fn(
    t = realiz$t, 
    K_params = K_params, 
    K_fn_consts = K_fn_consts, 
    ...
  )
  # K_tilde_hat <- mean(K_per_t) 
  
  # if (any(K_per_t < 1e-5)) return(-99999)
  
  # Restrict to non-negative parameters.
  if (min(mu, alpha, beta) < 0) return(-99999)
  
  # Pre-compute distances
  x_diff <- outer(realiz$x, realiz$x, "-")
  y_diff <- outer(realiz$y, realiz$y, "-")
  t_diff <- outer(realiz$t, realiz$t, "-")
  space_dist <- x_diff^2 + y_diff^2
  time_dist <- t_diff
  
  # For intlam to work according to Facilitated estimation of ETAS.
  mu_star <- mu / (spatstat.geom::area(window$s))
  
  # Calculate the large spatio-temporal integral directly
  spatio_temp_integral <- get_sti_K_per_t(realiz, bigT, bigX, bigY, alpha, beta, K_per_t)
  intlam <- mu * bigT + spatio_temp_integral
  
  # const <- K_weighted * (alpha / pi) * beta
  const <- (alpha / pi) * beta
  
  # Initialize log sum:
  sum_log <- log(mu_star)
  for (j in 2:nrow(realiz)) {
    gij <- sum(exp(-beta * time_dist[j, 1:(j - 1)] - alpha * space_dist[j, 1:(j - 1)]))
    lamj <- mu_star + K_per_t[j] * const * gij
    if (is.na(lamj) || lamj < 0) return(-99999)
    sum_log <- sum_log + log(lamj)
  }
  
  loglik <- sum_log - intlam
  if (loglik == -Inf) return(-99999)
  
  return(loglik)
}

#' Fit a Hawkes Process with Variable Productivity K
#'
#' Performs maximum likelihood estimation for a Hawkes process with a time-dependent \( K(t) \) function.
#'
#' @param K_inits Initial values for the \( K \) parameters to be optimized.
#' @param K_fn Function that calculates \( K(t) \) based on the provided parameters and constants.
#' @param K_fn_consts Constants for the \( K \) function, such as time splits.
#' @param stationary_params List containing background intensity (\code{mu}), spatial decay (\code{alpha}), and temporal decay (\code{beta}).
#' @param realiz Dataframe of observed events with columns:
#'   - \code{t}: Event times.
#'   - \code{x}: Event x-coordinates.
#'   - \code{y}: Event y-coordinates.
#' @param window List specifying the spatio-temporal domain. Should contain:
#'   - \code{t}: Numeric vector of time boundaries (\code{c(start, end)}).
#'   - \code{s}: Spatial window object (e.g., \code{spatstat.geom::owin}).
#' @param trace Integer specifying the verbosity level of the optimization process. Default is \code{0} (no output).
#' @param maxit Maximum number of iterations for the optimization algorithm. Default is \code{100}.
#' @param ... Additional arguments for \code{loglik_hawk_K_fn}.
#' @return A list containing the optimization results, including:
#'   - \code{par}: The fitted \( K \) parameters.
#'   - \code{value}: The log-likelihood at the optimum.
#'   - Other outputs from \code{optim}.
#' @examples
#' stationary_params <- list(mu = 1, alpha = 2, beta = 1)
#' K_inits <- c(1, 2)
#' K_fn <- K_FN_discrete
#' K_fn_consts <- c(0, 5)
#' realiz <- data.frame(t = c(1, 2), x = c(0.5, 0.6), y = c(0.4, 0.3))
#' window <- list(t = c(0, 10), s = spatstat.geom::owin(c(0, 1), c(0, 1)))
#' fit_hawkes_K_fn(K_inits, K_fn, K_fn_consts, stationary_params, realiz, window)
fit_hawkes_K_fn <- function(
  K_inits,
  K_fn,
  K_fn_consts,
  stationary_params,
  realiz,
  window,
  trace = 0,
  maxit = 100,
  ...
) {
  
  fit <- optim(
    par = K_inits,
    fn = loglik_hawk_K_fn,
    method = "Nelder-Mead",
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    K_fn = K_fn,
    K_fn_consts = K_fn_consts,
    stationary_params = stationary_params,
    realiz = realiz,
    window = window,
    ...
  )
  
  return(fit)
}

#' Fit a Hawkes Process with Variable Productivity K Across Intervals
#'
#' Fits a Hawkes process with spatio-temporal variation in the productivity parameter \( K(t) \), using a two-step approach:
#' 1. Estimates overall parameters (\( \mu \), \( \alpha \), \( \beta \)).
#' 2. Optimizes the \( K \) values across time intervals.
#'
#' @param overall_inits Initial values for the overall parameters (\( \mu \), \( \alpha \), \( \beta \)).
#' @param n_K Number of \( K \) parameters to estimate, corresponding to the number of time intervals.
#' @param K_fn Function that calculates \( K(t) \) based on the provided parameters and constants.
#' @param K_fn_consts Constants for the \( K \) function, such as time splits.
#' @param realiz Dataframe of observed events with columns:
#'   - \code{t}: Event times.
#'   - \code{x}: Event x-coordinates.
#'   - \code{y}: Event y-coordinates.
#' @param window List specifying the spatio-temporal domain. Should contain:
#'   - \code{t}: Numeric vector of time boundaries (\code{c(start, end)}).
#'   - \code{s}: Spatial window object (e.g., \code{spatstat.geom::owin}).
#' @param trace Integer specifying the verbosity level of the optimization process. Default is \code{0} (no output).
#' @param maxit Maximum number of iterations for the optimization algorithm. Default is \code{1000}.
#' @param output Logical flag. If \code{TRUE}, prints fitting progress and results. Default is \code{TRUE}.
#' @param ... Additional arguments for the optimization process.
#' @return A list containing:
#'   - \code{overall_params}: Fitted overall parameters (\( \mu \), \( \alpha \), \( \beta \)).
#'   - \code{K_params}: Fitted \( K \) parameters for each time interval.
#' @examples
#' overall_inits <- c(1, 2, 1)  # Initial guesses for mu, alpha, beta
#' n_K <- 3
#' K_fn <- K_FN_discrete
#' K_fn_consts <- c(0, 5, 10)
#' realiz <- data.frame(t = c(1, 2, 8), x = c(0.5, 0.6, 0.7), y = c(0.4, 0.3, 0.5))
#' window <- list(t = c(0, 10), s = spatstat.geom::owin(c(0, 1), c(0, 1)))
#' fit_var_prod_hawkes_K_fn(overall_inits, n_K, K_fn, K_fn_consts, realiz, window)
fit_var_prod_hawkes_K_fn <- function(
  overall_inits,
  n_K,
  K_fn,
  K_fn_consts,
  realiz,
  window,
  trace = 0,
  maxit = 1000,
  output = TRUE,
  ...
) {
  # Fit the overall mu, alpha, beta, K-hat values
  if (output) print("Fitting overall parameters...")
  stationary_params <- fit_hawkes_const_K(
    params_init = overall_inits,
    realiz = realiz,
    window = window,
    trace = trace,
    maxit = maxit
  )
  
  if (output) cat("Overall fitted parameters:\n", params_fitted$par, "\n")
  
  if (output) print("Fitting K parameters...")
  K_hat_varied <- fit_hawkes_K_fn(
    params_init = rep(overall_K_hat, n_K), 
    K_fn = K_fn,
    K_fn_consts = K_fn_consts,
    stationary_params = stationary_params,
    realiz = realiz,
    window = window,
    trace = trace,
    maxit = maxit,
    ...
  )
  
  return(K_hat_varied)
}

##### KAPPA FUNCTIONS #####

#' Function to calculate K values over time intervals
#'
#' Generates a time-varying parameter K over specified intervals, defined by change points.
#' This function handles both constant and linearly interpolated parameter values between
#' change points.
#'
#' @param t Numeric vector of time points.
#' @param k_params Numeric vector of K values at specified intervals.
#' @param splits Numeric vector of change points defining intervals.
#' @return Numeric vector of K values for each time point in \code{t}.
#' @export
K_FN_spline <- function(t, K_params, K_fn_consts) {
  all_K <- unlist(K_params)
  # Adjustments for midpoint splits instead of splits per K
  # all_splits <- c(0, splits, TIME)
  # midpoints <- all_splits[-length(all_splits)] + diff(all_splits) / 2
  midpoints <- K_fn_consts
  
  K_per_t <- numeric(length(t))
  
  # Loop through each t value
  for (i in seq_along(t)) {
    t_val <- t[i]
    # Check if t_val is before the first or after the last change point
    if (t_val < midpoints[1]) {
      K_per_t[i] <- all_K[1]  # Constant value before first change point
    } else if (t_val > midpoints[length(midpoints)]) {
      K_per_t[i] <- all_K[length(all_K)]  # Constant value after last change point
    } else {
      # Find the interval where t_val falls
      for (j in 1:(length(midpoints) - 1)) {
        if (t_val >= midpoints[j] && t_val <= midpoints[j + 1]) {
          # Linear interpolation
          slope <- (all_K[j + 1] - all_K[j]) / (midpoints[j + 1] - midpoints[j])
          K_per_t[i] <- all_K[j] + slope * (t_val - midpoints[j])
          break
        }
      }
    }
  }
  return(K_per_t)
}

#' Linear \( K(t) \) Function
#'
#' Computes a linear \( K(t) \) function of the form \( K(t) = a + b \cdot t \).
#'
#' @param t Numeric vector of time points.
#' @param k_params Numeric vector with two coefficients: \code{a} (intercept) and \code{b} (slope).
#' @param splits Not used in this function (included for consistency).
#' @return Numeric vector of \( K(t) \) values.
#' @examples
#' K_FN_deg_1(c(1, 2, 3), c(2, 0.5), NULL)  # Returns c(2.5, 3, 3.5)
K_FN_deg_1 <- function(t, k_params, splits) {
  all_K <- unlist(k_params)
  K_per_t <- all_K[1] + all_K[2] * t
  return(K_per_t)
}

#' Quadratic \( K(t) \) Function
#'
#' Computes a quadratic \( K(t) \) function of the form \( K(t) = a + b \cdot t + c \cdot t^2 \).
#'
#' @param t Numeric vector of time points.
#' @param coefs Numeric vector with three coefficients: \code{a}, \code{b}, and \code{c}.
#' @param consts Not used in this function (included for consistency).
#' @return Numeric vector of \( K(t) \) values.
#' @examples
#' K_FN_deg_2(c(1, 2, 3), c(1, -1, 0.5), NULL)  # Returns c(0.5, 3, 7.5)
K_FN_deg_2 <- function(t, coefs, consts) {
  coefs <- unlist(coefs)
  K_per_t <- coefs[1] + coefs[2] * t + coefs[3] * t^2
  return(K_per_t)
}

#' Polynomial \( K(t) \) Function
#'
#' Computes a polynomial \( K(t) \) function of a given degree.
#'
#' @param t Numeric vector of time points.
#' @param coefs Numeric vector of coefficients, starting with the intercept and increasing powers.
#' @param degree Degree of the polynomial.
#' @return Numeric vector of \( K(t) \) values.
#' @examples
#' K_FN_poly(c(1, 2), c(1, 2, 1), 2)  # Returns c(4, 9) for \( K(t) = 1 + 2t + t^2 \)
K_FN_poly <- function(t, coefs, degree){
  if (length(coefs) != degree + 1) {
    stop("The number of coefficients must be equal to the degree of the polynomial plus one.")
  }
  
  # Evaluate the polynomial for each value in t
  sapply(t, function(x) sum(coefs * x^(0:degree)))
}

#' Sinusoidal \( K(t) \) Function
#'
#' Computes a sinusoidal \( K(t) \) function of the form:
#' \[
#' K(t) = A \sin(B (t - C)) + D
#' \]
#'
#' @param t Numeric vector of time points.
#' @param coefs Numeric vector with four coefficients:
#'   - \code{A}: Amplitude,
#'   - \code{B}: Frequency,
#'   - \code{C}: Phase shift,
#'   - \code{D}: Vertical shift.
#' @param consts Not used in this function (included for consistency).
#' @return Numeric vector of \( K(t) \) values.
#' @examples
#' K_FN_sin(c(1, 2, 3), c(1, pi, 0, 0), NULL)  # Computes sinusoidal \( K(t) \)
K_FN_sin <- function(t, coefs, consts) {
  return(coefs[1] * sin(coefs[2] * (t - coefs[3])) + coefs[4])
}


#' Smoothed Logistic \( K(t) \) Function
#'
#' Computes a \( K(t) \) function with transitions between intervals smoothed using a logistic function.
#'
#' @param t Numeric vector of time points.
#' @param k_params Numeric vector with smoothing parameter \code{sigma} and interval-specific \( K \) values.
#' @param splits Numeric vector of time interval boundaries.
#' @return Numeric vector of \( K(t) \) values.
#' @examples
#' K_FN_smoothed_logistic(c(1, 2, 3), c(0.1, 1, 2), c(0, 2, 4))  # Computes smoothed \( K(t) \)
K_FN_smoothed_logistic <- function(t, k_params, splits) {
  sigma = k_params[1]
  all_K <- unlist(k_params[-1])
  
  if (sigma == 0) {
    # Original vertical connection
    t_intervals <- findInterval(t, splits) + 1
    # Handle out-of-bound indices
    t_intervals[t_intervals == 0] <- 1 # Below the first split
    t_intervals[t_intervals > length(all_K)] <- length(all_K) # Above the last split
    K_per_t <- all_K[t_intervals]
  } else {
    # Smooth transition
    for (i in seq_along(t)) {
      # Find the intervals that `t[i]` falls between
      lower_idx <- max(which(splits <= t[i]), 0)
      upper_idx <- min(which(splits > t[i]), length(all_K) + 1)
      
      # Handle edge cases
      if (lower_idx == 0) {
        K_per_t[i] <- all_K[1] # Below the first split
      } else if (upper_idx > length(all_K)) {
        K_per_t[i] <- all_K[length(all_K)] # Above the last split
      } else {
        # Smooth interpolation
        lower_split <- splits[lower_idx]
        upper_split <- splits[upper_idx - 1]
        
        weight <- 1 / (1 + exp(-(t[i] - lower_split) / sigma))
        K_per_t[i] <- all_K[lower_idx] * (1 - weight) + all_K[upper_idx] * weight
      }
    }
  }
  
  return(K_per_t)
}

#' Smoothed Polynomial \( K(t) \) Function
#'
#' Computes a \( K(t) \) function with transitions between intervals smoothed using a polynomial function.
#'
#' @param t Numeric vector of time points.
#' @param k_params Numeric vector with smoothing parameter \code{degree} and interval-specific \( K \) values.
#' @param splits Numeric vector of time interval boundaries.
#' @return Numeric vector of \( K(t) \) values
K_FN_smoothed_poly <- function(t, k_params, splits) {
  all_K <- unlist(k_params)
  degree = splits[1]
  splits = splits[-1]
  
  for (i in seq_along(t)) {
    # Find the lower and upper intervals
    lower_idx <- max(which(splits <= t[i]), 0)
    upper_idx <- min(which(splits > t[i]), length(splits)+1)
    
    if (lower_idx == 0) {
      # Case: t[i] is below the first split
      K_per_t[i] <- all_K[1]
    } else if (upper_idx > length(all_K)) {
      # Case: t[i] is above the last split
      K_per_t[i] <- all_K[length(all_K)]
    } else {
      # Case: t[i] is within an interval
      lower_split <- splits[lower_idx]
      upper_split <- splits[upper_idx]
      
      # Compute the weight w
      if (upper_split != lower_split) {
        w <- (t[i] - lower_split) / (upper_split - lower_split)
      } else {
        w <- 0 # Prevent division by zero in degenerate cases
      }
      
      # Apply polynomial interpolation
      K_per_t[i] <- all_K[lower_idx] * (1 - w^degree) + all_K[upper_idx] * w^degree
    }
  }
  
  return(K_per_t)
}

#' Smoothed Exponential \( K(t) \) Function
#'
#' Computes a \( K(t) \) function with transitions between intervals smoothed using an exponential function.
#'
#' @param t Numeric vector of time points.
#' @param k_params Numeric vector with smoothing parameter \code{alpha} and interval-specific \( K \) values.
#' @param splits Numeric vector of time interval boundaries.
#' @return Numeric vector of \( K(t) \) values
K_FN_smoothed_exp <- function(t, k_params, splits) {
  alpha = k_params[1]
  all_K <- unlist(k_params[-1])
  
  for (i in seq_along(t)) {
    lower_idx <- max(which(splits <= t[i]), 0)
    upper_idx <- min(which(splits > t[i]), length(all_K) + 1)
    
    if (lower_idx == 0) {
      K_per_t[i] <- all_K[1]
    } else if (upper_idx > length(all_K)) {
      K_per_t[i] <- all_K[length(all_K)]
    } else {
      lower_split <- splits[lower_idx]
      upper_split <- splits[upper_idx - 1]
      decay <- exp(-alpha * (t[i] - lower_split))
      K_per_t[i] <- all_K[lower_idx] * decay + all_K[upper_idx] * (1 - decay)
    }
  }
  
  return(K_per_t)
}

#' Smoothed Hyperbolic \( K(t) \) Function
#'
#' Computes a \( K(t) \) function with transitions between intervals smoothed using a hyperbolic function.
#'
#' @param t Numeric vector of time points.
#' @param k_params Numeric vector with smoothing parameter \code{sigma} and interval-specific \( K \) values.
#' @param splits Numeric vector of time interval boundaries.
#' @return Numeric vector of \( K(t) \) values
K_FN_smoothed_tanh <- function(t, k_params, splits) {
  sigma = k_params[1]
  all_K <- unlist(k_params[-1])
  
  for (i in seq_along(t)) {
    lower_idx <- max(which(splits <= t[i]), 0)
    upper_idx <- min(which(splits > t[i]), length(all_K) + 1)
    
    if (lower_idx == 0) {
      K_per_t[i] <- all_K[1]
    } else if (upper_idx > length(all_K)) {
      K_per_t[i] <- all_K[length(all_K)]
    } else {
      lower_split <- splits[lower_idx]
      upper_split <- splits[upper_idx - 1]
      w <- (1 + tanh((t[i] - lower_split) / sigma)) / 2
      K_per_t[i] <- all_K[lower_idx] * (1 - w) + all_K[upper_idx] * w
    }
  }
  
  return(K_per_t)
}

#' Kernel-Weighted \( K(t) \) Function
#'
#' Computes \( K(t) \) as a weighted average of values at specified knots, using Gaussian kernels.
#'
#' @param t Numeric vector of time points.
#' @param K_params Numeric vector of \( K \) values at knots.
#' @param K_fn_consts Numeric vector of knot locations (centers of the kernels).
#' @param bandwidth Bandwidth of the Gaussian kernel. Default is \code{1}.
#' @return Numeric vector of \( K(t) \) values.
#' @examples
#' K_FN_kernel(c(1, 2, 3), c(1, 2, 3), c(0, 2, 4), 1)  # Computes kernel-smoothed \( K(t) \)
K_FN_kernel <- function(t, K_params, K_fn_consts, bandwidth = 1) {
  # K_fn_consts are the time knots (centers of the kernels)
  knots <- K_fn_consts  # Define the time intervals for kernel centers
  K_values <- unlist(K_params)  # Values of K at the knots
  
  # Gaussian kernel function
  gaussian_kernel <- function(x, mu, bw) {
    exp(-0.5 * ((x - mu) / bw)^2) / (sqrt(2 * pi) * bw)
  }
  
  # Compute K(t) as a weighted average of K_params using Gaussian kernels
  K_per_t <- sapply(t, function(ti) {
    weights <- sapply(knots, function(knot) gaussian_kernel(ti, knot, bandwidth))
    sum(weights * K_values) / sum(weights)  # Normalize weights
  })
  
  return(K_per_t)
}

##### ITERATIVE CHANGE POINT DETECTION METHODS #####

#' Objective Function for Parameter Estimation in Hawkes Process
#'
#' Computes the sum of absolute differences between estimated and actual parameter values
#' for a Hawkes process model with change points.
#'
#' @param K Numeric vector of current parameter estimates.
#' @param t_star Numeric vector specifying change points for the process.
#' @param K_fn Function that defines \( K(t) \) based on change points and parameters.
#' @param window List with two elements:
#'   - \code{t}: Numeric vector specifying the observation time window as \code{(start, end)}.
#'   - \code{s}: Numeric vector specifying the spatial window.
#' @param overall_K_hat Numeric, initial overall \( K \) estimate.
#' @param realiz Data frame of observed realizations, with columns \code{t} (time) and other relevant data.
#' @return Numeric value representing the objective function result (sum of absolute differences).
#' @examples
#' obj_fn(c(1, 2), c(5, 10), K_FN_discrete, list(t = c(0, 20), s = NULL), 1.5, data.frame(t = c(1, 2, 3)))
#' @export
obj_fn <- function(K, t_star, K_fn, window, overall_K_hat, realiz) {
  K_per_t <- K_fn(realiz$t, K, t_star)
  return(sum(abs(overall_K_hat - K_per_t)))
}

#' Fit Change Points in a Hawkes Process Model
#'
#' Iteratively estimates parameters for a Hawkes process with given change points.
#'
#' @param t_star Numeric vector of initial change points.
#' @param K_fn Function that defines \( K(t) \) based on change points.
#' @param objective_fn Function to compute the objective, such as \code{obj_fn}.
#' @param mu_hat Numeric, baseline intensity estimate.
#' @param alpha_hat Numeric, branching ratio estimate.
#' @param beta_hat Numeric, decay rate estimate.
#' @param overall_K_hat Numeric, initial estimate of \( K \).
#' @param realiz Data frame of process realizations, with columns \code{t} (time) and other relevant data.
#' @param window List with two elements:
#'   - \code{t}: Numeric vector specifying the time window as \code{(start, end)}.
#'   - \code{s}: Numeric vector specifying the spatial window.
#' @param output Logical, whether to print intermediate results. Default is \code{FALSE}.
#' @return List with:
#'   - \code{obj_fn_val}: Objective function value.
#'   - \code{K}: Estimated parameter values for each segment.
#' @examples
#' fit_t_star_iterative(
#'   t_star = c(5, 10),
#'   K_fn = K_FN_discrete,
#'   objective_fn = obj_fn,
#'   mu_hat = 0.5,
#'   alpha_hat = 0.8,
#'   beta_hat = 1.2,
#'   overall_K_hat = 1.5,
#'   realiz = data.frame(t = c(1, 2, 3)),
#'   window = list(t = c(0, 20), s = NULL)
#' )

fit_t_star_iterative <- function(
  t_star,
  K_fn = K_FN_discrete,
  objective_fn,
  mu_hat,
  alpha_hat,
  beta_hat,
  overall_K_hat,
  realiz,
  window,
  K_max = 2,
  trace = 0,
  output = FALSE) {
  
  if ((any(t_star > window$t[2] - 1e-4)) || (any(t_star < window$t[1] + 1e-4))) {
    return(list(obj_fn_val = -999999))
  }
  
  splits <- c(t_star, window$t[2])
  
  if (!(all(diff(splits) > 0))) return(list(obj_fn_val = -999999))
  
  K_hat <- c()
  for (intrvl in 1:length(splits)) {
    K_hat_varied <- fit_hawkes_varied_K(
      params_init = overall_K_hat, 
      previous_K = K_hat,
      splits = splits,
      mu = mu_hat,
      alpha = alpha_hat,
      beta = beta_hat,
      K_max = K_max,
      trace = trace,
      realiz = realiz[realiz$t < splits[intrvl],],
      window = list(t = c(0, splits[intrvl]), s = window$s),
    )
    
    K_hat <- c(K_hat, K_hat_varied$par)
  }
  
  obj_fn_val <- objective_fn(K_hat, t_star, K_fn, window, overall_K_hat, realiz)
  
  if (output) cat("t_star:\t", t_star, "\nobj_val:", obj_fn_val, "\nK:\t", K_hat, "\n\n")
  
  return(list(obj_fn_val = obj_fn_val, K = K_hat))
}

#' Wrapper Function for Optimizing Change Points
#'
#' A wrapper for optimizing change points using a specified objective function.
#'
#' @param t_star Numeric vector of initial change points.
#' @param ... Additional arguments passed to \code{fit_t_star_iterative}.
#' @return Numeric value of the objective function result.
optim_wrapper_iterative <- function(t_star, ...) {
  result <- fit_t_star_iterative(t_star, ...)
  return(result$obj_fn_val)
}

#' Estimate Change Points in a Hawkes Process
#'
#' Identifies optimized change points and corresponding parameter estimates for a Hawkes process.
#'
#' @param realiz Data frame of observed realizations, with columns \code{t} (time) and other relevant data.
#' @param param_hats Named list of estimated parameters:
#'   - \code{mu}: Baseline intensity.
#'   - \code{alpha}: Branching ratio.
#'   - \code{beta}: Decay rate.
#'   - \code{K}: Initial overall \( K \).
#' @param n_cps Integer, number of change points.
#' @param K_fn Function that defines \( K(t) \) over time.
#' @param objective_fn Function to minimize (e.g., \code{obj_fn}).
#' @param window List with two elements:
#'   - \code{t}: Numeric vector specifying the time window as \code{(start, end)}.
#'   - \code{s}: Numeric vector specifying the spatial window.
#' @param output Logical, whether to print detailed output. Default is \code{FALSE}.
#' @param trace Integer, trace level for optimization. Default is \code{0}.
#' @param maxit Integer, maximum number of iterations for optimization. Default is \code{20}.
#' @return List containing:
#'   - \code{t_star}: Optimized change points.
#'   - \code{K}: Estimated parameters for each segment.
#' @examples
#' get_change_points_iterative(
#'   realiz = data.frame(t = c(1, 2, 3)),
#'   param_hats = list(mu = 0.5, alpha = 0.8, beta = 1.2, K = 1.5),
#'   n_cps = 2,
#'   K_fn = K_FN_discrete,
#'   objective_fn = obj_fn,
#'   window = list(t = c(0, 20), s = NULL)
#' )
get_change_points_iterative <- function(
  realiz,
  param_hats,
  n_cps,
  objective_fn,
  window,
  K_max = 2,
  output = FALSE,
  trace = 0,
  maxit = 20
) {
  t_end = window$t[2]
  
  mu_hat <- param_hats$mu
  alpha_hat <- param_hats$alpha
  beta_hat <- param_hats$beta
  overall_K_hat <- param_hats$K
  
  n_K <- n_cps + 1
  
  if (n_cps == 1) {
    optim_method <- "Brent"
    lower_lim <- window$t[1]
    upper_lim <- window$t[2]
  } else {
    optim_method <- "Nelder-Mead"
    lower_lim <- NULL
    upper_lim <- NULL
  }
  
  t_star_fit <- optim(
    par = seq(t_end/n_K, t_end - t_end/n_K, t_end/n_K),
    fn = optim_wrapper_iterative,
    method = optim_method,
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    objective_fn = obj_fn,
    mu_hat = mu_hat,
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    overall_K_hat = overall_K_hat,
    realiz = realiz,
    window = window,
    K_max = K_max,
    trace = trace,
    output = output,
    lower = lower_lim,
    upper = upper_lim
  )
  
  final_result <- fit_t_star_iterative(
    t_star = t_star_fit$par, 
    objective_fn = obj_fn, 
    mu_hat = mu_hat, 
    alpha_hat = alpha_hat, 
    beta_hat = beta_hat, 
    overall_K_hat = overall_K_hat, 
    realiz = realiz, 
    window = window,
  )
  
  t_star <- t_star_fit$par + window$t[1]
  
  return(list(t_star = t_star[order(t_star)], K = final_result$K))
}

#' Plot change points for a Hawkes process
#'
#' Visualizes change points and estimated K values over time.
#'
#' @param realiz Data frame of observed realizations.
#' @param K_fn Function defining K over time.
#' @param K_true Numeric vector of true K values.
#' @param K_hat Numeric vector of estimated K values.
#' @param t_star_true Numeric vector of true change points.
#' @param t_star_hat Numeric vector of estimated change points.
#' @param bigT Numeric, end time for the process.
#' @return ggplot object visualizing change points and estimated K values.
plot_change_points <- function(
  realiz,
  K_fn,
  K_true,
  K_hat,
  t_star_true,
  t_star_hat,
  bigT
) {
  splits <- c(t_star_true, bigT)
  K_per_t <- K_fn(realiz$t, K_true, t_star_true)
  
  n_K <- length(K_true)
  result_plot <- ggplot(mapping = aes(x = realiz$t, y = K_per_t)) +
    geom_point(alpha = 0.1) +
    ylim(0, 1) +
    geom_segment(aes(x = 0, xend = t_star_hat[1], y = K_hat[1], yend = K_hat[1]), linetype = "dashed", color = "red")  +
    geom_segment(aes(x = t_star_hat[length(t_star_hat)], xend = bigT, y = K_hat[n_K], yend = K_hat[n_K]), linetype = "dashed", color = "red") +
    geom_segment(aes(x = t_star_hat[1], xend = t_star_hat[1], y = K_hat[1], yend = K_hat[2]), linetype = "dashed", color = "blue")
  
  if (n_K > 2) {
    for (i in 2:3) {
      result_plot <- result_plot + geom_segment(aes_string(x = t_star_hat[i-1], xend = t_star_hat[i], y = K_hat[i], yend = K_hat[i]), linetype = "dashed", color = "red")
      result_plot <- result_plot + geom_segment(aes_string(x = t_star_hat[i], xend = t_star_hat[i], y = K_hat[i], yend = K_hat[i+1]), linetype = "dashed", color = "blue")
    }
  }
  
  plot(result_plot)
}


##### GENERAL CHANGE POINT DETECTION #####

#' Fits Change Points for the Hawkes Process
#'
#' Optimizes change points and associated parameter estimates for a Hawkes process
#' based on observed data and given stationary parameter estimates.
#'
#' @param t_star Numeric vector of initial change points.
#' @param n_K Integer, number of intervals for the time-varying \( K(t) \).
#' @param K_fn Function defining \( K(t) \) based on change points.
#' @param objective_fn Objective function to compute goodness of fit, such as \code{obj_fn}.
#' @param stationary_params Named list of stationary parameters:
#'   - \code{mu}: Baseline intensity.
#'   - \code{alpha}: Branching ratio.
#'   - \code{beta}: Decay rate.
#'   - \code{K}: Initial overall \( K \).
#' @param realiz Data frame of observed realizations, with columns \code{t} (time) and other relevant data.
#' @param window List with two elements:
#'   - \code{t}: Numeric vector specifying the observation time window as \code{(start, end)}.
#'   - \code{s}: Numeric vector specifying the spatial window.
#' @param output Logical, if \code{TRUE}, prints intermediate results. Default is \code{FALSE}.
#' @param ... Additional arguments passed to the underlying optimization routines.
#' @return A list containing:
#'   - \code{obj_fn_val}: Value of the objective function.
#'   - \code{K}: Estimated parameters for each interval.
#' @export
fit_t_star <- function(
  t_star,
  n_K,
  K_fn = K_FN_discrete,
  objective_fn = obj_fn,
  stationary_params,
  realiz,
  window,
  output = F,
  ...
) {
  # Penalize if any change point is too close to the window boundary
  if ((any(t_star > window$t[2] - 1e-4)) || (any(t_star < window$t[1] + 1e-4))) {
    return(list(obj_fn_val = -999999))
  }
  
  # Ensure splits are strictly increasing
  if (!(all(diff(t_star) > 0))) return(list(obj_fn_val = -999999))
  
  # Fitting for discrete intervals (more efficient than K_fn)
  K_hat_varied <- fit_hawkes_K_fn(
    K_inits = rep(stationary_params$K, n_K), 
    K_fn = K_fn,
    K_fn_consts = t_star,
    stationary_params = stationary_params,
    realiz = realiz,
    window = window,
    ...
  )
  
  obj_fn_val <- objective_fn(K_hat_varied$par, t_star, K_fn, window, stationary_params$K, realiz)
  
  if (output) cat("t_star:\t", t_star, "\nobj_val:", obj_fn_val, "\nK:\t", K_hat_varied$par, "\n\n")
  
  return(list(obj_fn_val = obj_fn_val, K = K_hat_varied$par))
}

#' Wrapper for Optimization of Change Points in Hawkes Process
#'
#' Calls \code{fit_t_star} with the specified arguments to optimize the objective
#' function for given initial change points.
#'
#' @param t_star Numeric vector of initial change points.
#' @param ... Additional arguments passed to \code{fit_t_star}.
#' @return Value of the objective function from \code{fit_t_star}.
optim_wrapper <- function(t_star, ...) { 
  result <- fit_t_star(t_star, ...)
  return(result$obj_fn_val)
}

#' Estimate Change Points for Hawkes Process with Parameter Variation
#'
#' Identifies the optimal locations of change points and corresponding parameter estimates
#' for a Hawkes process using observed realizations and provided stationary parameter estimates.
#'
#' @param realiz Data frame of observed realizations, with columns \code{t} (time) and other relevant data.
#' @param stationary_params Named list of stationary parameters:
#'   - \code{mu}: Baseline intensity.
#'   - \code{alpha}: Branching ratio.
#'   - \code{beta}: Decay rate.
#'   - \code{K}: Initial overall \( K \).
#' @param n_cps Integer, number of change points to estimate.
#' @param n_K Integer, number of intervals for the time-varying \( K(t) \).
#' @param K_fn Function defining \( K(t) \) over time.
#' @param objective_fn Objective function to minimize, such as \code{obj_fn}.
#' @param window List with two elements:
#'   - \code{t}: Numeric vector specifying the observation time window as \code{(start, end)}.
#'   - \code{s}: Numeric vector specifying the spatial window.
#' @param output Logical, if \code{TRUE}, prints detailed output. Default is \code{FALSE}.
#' @param trace Integer, trace level for optimization diagnostics. Default is \code{0}.
#' @param maxit Integer, maximum iterations for optimization. Default is \code{10}.
#' @param ... Additional arguments passed to the optimization routines.
#' @return A list containing:
#'   - \code{t_star}: Optimized change points (sorted).
#'   - \code{K}: Estimated parameters for each interval.
#' @examples
#' get_change_points(
#'   realiz = data.frame(t = c(1, 2, 3)),
#'   stationary_params = list(mu = 0.5, alpha = 0.8, beta = 1.2, K = 1.5),
#'   n_cps = 2,
#'   n_K = 3,
#'   K_fn = K_FN_discrete,
#'   objective_fn = obj_fn,
#'   window = list(t = c(0, 20), s = NULL)
#' )
#' @export
get_change_points <- function(
  realiz,
  stationary_params,
  n_cps,
  n_K,
  K_fn = K_FN_discrete,
  objective_fn = obj_fn,
  window,
  output = FALSE,
  trace = 0,
  maxit = 10,
  ...
){
  t_end <- window$t[2]
  init_seq <- seq(t_end/(n_cps+1), t_end - t_end/(n_cps+1), t_end/(n_cps+1))
  
  # If only one change point, use "Brent" method in optim as "Nelder-Mead" only
  # for multiple parameters.
  if (n_cps == 1) {
    optim_method <- "Brent"
    lower_lim <- window$t[1]
    upper_lim <- window$t[2]
  } else {
    optim_method <- "Nelder-Mead"
    lower_lim <- -Inf
    upper_lim <- Inf
  }
  
  t_star_fit <- optim(
    par = init_seq,
    fn = optim_wrapper,
    method = optim_method,
    control = list(fnscale = -1, trace = trace, maxit = maxit),
    objective_fn = objective_fn,
    stationary_params = stationary_params,
    realiz = realiz,
    window = window,
    output = output,
    n_K = n_K,
    K_fn = K_fn,
    lower = lower_lim,
    upper = upper_lim,
    ...
  )
  
  final_result <- fit_t_star(
    t_star = t_star_fit$par, 
    objective_fn = objective_fn, 
    stationary_params = stationary_params,
    realiz = onerealiz,
    n_K = n_K,
    K_fn = K_fn,
    window = window,
    ...
  )
  
  t_star <- t_star_fit$par + window$t[1]
  
  return(list(t_star = t_star[order(t_star)], K = final_result$K))
}