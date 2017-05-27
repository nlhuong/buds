#' Fit Bayesian Unidimensional Scaling model to
#'
#' Fitting BUDS model to recover latent coordinates of data points from
#' observed/precomputed pairwise dissimilarities
#'
#' @param D A square matrix of pairwise dissimilarities
#' @param K An integer indicating k-nearest neighbours to compute average
#' distances to kNN then used as weights in the power likelihood.
#' @param method Stan mode to use either "mcmc" or "vb"
#' @param hyperparams A list of BUDS hiperparameters with named components
#' set by default as follows "gamma_tau"= 2.5, "gamma_epsilon" = 2.5, 
#' "gamma_bias" = 2.5, "gamma_rho_sq" = 2.5, "min_sigma" = 0.01
#' @param init_from Initialization option either "random" or from a location
#' along "principal_curve" fitted to 2D representation using PCoA on D
#' @param seed A integer representing a random seed
#' @param max_trials A integer for maximum number of trials to fit data
#' @param ... Other parameters for vb and stan
#'
#' @return The fitted `stanfit` model object
#' @export
fit_buds <- function(D, K = NULL,
                     method = c("vb", "mcmc"),
                     hyperparams = list(
                       "gamma_tau"= 2.5,
                       "gamma_epsilon" = 2.5,
                       "gamma_bias" = 2.5,
                       "gamma_rho" = 2.5,
                       "min_sigma" = 0.03),
                     init_from = c("random", "principal_curve"),
                     seed = 1234, max_trials = 20, ...) {
  # Method can be eithe Variational Bayes or MCMC
  method <- match.arg(method, c("vb", "mcmc"))
  # Parameter initialization
  init <- match.arg(init_from, c("random", "principal_curve"))
  if (init == "principal_curve") {
    offset <- 1e-3
    tau0 <- prin_curve_loc(D)
    tau0 <- (tau0 - min(tau0) + offset) / (max(tau0) - min(tau0) + 2*offset)
    init <- list("tau" = tau0, "bias" = offset, "rho" = 1.0, 
                 "meansd" = hyperparams$min_sigma,
                 "tau_shape1" = 1.0, "tau_shape2" = 1.0)
  }
  # Number of data points
  N <- ncol(D)
  # Convert D matrix to long form and get kNN distances
  distDF <- get_dist_df(D, K = K)
  # Data for stan model
  dist_data <- list(
    "N" = N,
    "Npairs" = nrow(distDF),
    "i_idx" = distDF$i,
    "j_idx" = distDF$j,
    "dvec" = distDF$d,
    "rel_sd" = sqrt(distDF$v/mean(distDF$v))
  )
  stan_data <- utils::modifyList(hyperparams, dist_data)
  # Fit the stan model
  if(!is.numeric(seed)) seed <- sample.int(.Machine$integer.max, 1)
  if(method == "vb") {
    while(max_trials > 0) {
      fit <-try(rstan::vb(stanmodels$buds, data = stan_data, init = init, seed = seed, ...))
      if(class(fit) != "try-error") break
      seed <- sample.int(.Machine$integer.max, 1)
      max_trials <- max_trials - 1
    }
  } else {
    fit <- sampling(stanmodels$buds, data = stan_data, init = init, seed = seed, ...)
  }
  return(list(fit_buds = fit, seed = seed, distDF = distDF, 
              stan_data = stan_data, init = init))
}
