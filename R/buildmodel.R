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
#' set by default as follows "tau_shape1" = 1, "tau_shape2"= 1,
#' "gamma_epsilon" = 2.5, "gamma_bias" = 2.5, "gamma_rho_sq" = 2.5,
#' "min_sigma" = 0.001
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
                       "tau_shape1" = 1,
                       "tau_shape2"= 1,
                       "gamma_epsilon" = 2.5,
                       "gamma_bias" = 2.5,
                       "gamma_rho_sq" = 2.5,
                       "min_sigma" = 0.001),
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
    init <- list("tau" = tau0, "bias" = offset, "rho_sq" = 1.0,
                 "mean_var" = 0.05)
  }
  # Number of data points
  N <- ncol(D)
  # Convert D matrix to long form and get kNN distances
  distDF <- get_dist_df(D, K = K)
  # Compute weights for observations
  weights <- distDF$min_sigma_K
  weights <- weights/max(weights)
  weights <- pmax(weights, rep(0.1, length(weights)))
  # Data for stan model
  dist_data <- list(
    "N" = N,
    "Npairs" = length(distDF$i),
    "i_idx" = distDF$i,
    "j_idx" = distDF$j,
    "dvec" = distDF$d,
    "weight" = weights
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
  return(list(fit_buds = fit, seed = seed, distDF = distDF, weights = weights))
}


#' Reshape distance matrix to long format
#'
#' Internal function for getting distances in a long data.frame format,
#' a format used in fit_buds function
#'
#' @param D A square matrix of pairwise dissimilarities
#' @param K  An integer indicating k-nearest neighbours for computing average 
#' distance to kNN.
#'
#' @return A data frame to be used as input data for BUDS.
get_dist_df <- function(D, K = NULL) {
  N <- nrow(D)
  rownames(D) <- colnames(D) <- 1:N
  distDF <- reshape2::melt(D)
  colnames(distDF) <- c("j", "i", "d")
  # Keep only (n-1)n/2 distances as you assume they are symmetric
  distDF <- distDF[distDF$i < distDF$j, ]
  # Delete all the missing value distances
  distDF <- distDF[!is.na(distDF$d), ]
  # Compute average distance to kNN
  dkNN <- kNN_dist(D, K = K)
  distDF$i_sigma_K <- dkNN$mean[distDF$i]
  distDF$j_sigma_K <- dkNN$mean[distDF$j]
  distDF$max_sigma_K <- pmax(distDF$i_sigma_K, distDF$j_sigma_K)
  distDF$min_sigma_K <- pmin(distDF$i_sigma_K, distDF$j_sigma_K)
  return(distDF)
}

#' k-nearest neighbours distance.
#' 
#' Compute average distance to k-nearest neighbours
#' for each data point from a dissimilarity matrix.
#' This approximates the data densities around 
#' each data point.
#' 
#' @param D A dissimilarity matrix.
#' @param K An integer indicating the number of points 
#' to include in the neighbourhood. If not specified 
#' one tenth of the total number of data points is used.
#' 
#' @return A list with mean distances to k-nearest 
#' neighbours and the associated variance of these
#' k-distances for each data point.
#' 
#' @export
kNN_dist <- function(D, K = NULL) {
  if (is.null(K) || is.na(K)) K <- floor(ncol(D)/10)
  D_kSort <- apply(D, 2, function(x) {sort(x)})
  mean <- apply(D_kSort, 2, function(x){
    mean(x[2:(K+1)]) 
  })
  var <- apply(D_kSort, 2, function(x){
    var(x[2:(K+1)]) 
  })
  return(list(mean = mean, var = var))
}


#' Principal curve locations
#'
#' Estimating locations of data points projected on a principal curve fitted
#' to a PCoA 2D representation
#'
#' @param D A square matrix of pairwise dissimilarities
#' @param ... Other parameters for princurve::principal.curve function.
#'
#' @return A numeric vector of positions along the fitted principal curve.
#' @export
prin_curve_loc <- function(D, ...) {
  pcoa_res <- stats::cmdscale(sqrt(D), k = 2, eig = TRUE)
  X <- as.matrix(pcoa_res$points)
  prin_curve <- princurve::principal.curve(X, plot = FALSE, ...)
  tau <- prin_curve$lambda
  tau <- (tau - min(tau))/(max(tau) - min(tau))
  return(tau)
}
