#' Principal curve locations
#'
#' Estimating locations of data points projected on a principal curve fitted
#' to a PCoA 2D representation
#'
#' @param D A square matrix of pairwise dissimilarities
#' @param \\dots Other parameters for princurve::principal.curve function.
#'
#' @return A numeric vector of positions along the fitted principal curve.
#' @export
prin_curve_loc <- function(D, ...) {
  pcoa_res <- cmdscale(sqrt(D), k = 2, eig = TRUE)
  X <- as.matrix(pcoa_res$points)
  prin_curve <- princurve::principal.curve(X, plot = FALSE, ...)
  tau <- prin_curve$lambda
  tau <- (tau - min(tau))/(max(tau) - min(tau))
  return(tau)
}


#' Data subset along gradient
#'
#' Select a subset of samples evenly spaced along
#' the underlying gradient using the latent coordinates
#' for the data points on a [0, 1] interval.
#'
#' @param tau A numeric vector representing the latent
#' positions of data points.
#' @param nkeep An integer indicating the number of samples
#' to keep.
#'
#' @return A vector with names of selected samples.
#' @export
get_subset_tau <- function(tau, nkeep) {
  if(is.null(names(tau))) names(tau) <- 1:length(tau)
  tau <- sort(tau)
  idx <- sapply(seq(min(tau), max(tau), length.out = nkeep),
                function(x) sum(tau >= x))
  idx <- length(tau) - idx + 1
  return(names(tau)[idx])
}


#' Extract tau
#'
#' Extracts tau parameter from the BUDS fit object. Computes
#' the mode and upper and lower highes posterior density
#' interval.
#'
#' @param budsParams A list of parameters extracted from BUDS stanfit object.
#' @param prob [Optional] A fraction for the probability of HPD interval.
#' Default is 0.95.
#'
#' @importFrom coda HPDinterval
#' @importFrom stats median
#' @return A data frame with columns corresponding
#' to mode tau estimates and the lower and upper HPD.
#'
#' @export
get_tau_df <- function(budsParams, prob = 0.95) {
  tau_mcmc <- coda::mcmc(budsParams$tau)
  tau_mode <- MCMCglmm::posterior.mode(tau_mcmc)
  tau_intv <- HPDinterval(tau_mcmc, prob = prob)
  resdf <- data.frame(tau_mean = colMeans(budsParams$tau),
                      tau_median = apply(budsParams$tau, 2, median),
                      tau_mode = tau_mode,
                      tau_lower = tau_intv[, 1],
                      tau_upper = tau_intv[, 2])
  return(resdf)
}

#' Low dimensional data representation.
#' Compute reduced representations of the data.
#' A wrapper for cmdscale and Rtsne::Rtsne()
#'
#' @param D A pairwise dissimilarity matrix.
#' @param method [Optional] A dimensionality reduction
#' method either "PCoA" or "tSNE". Default is "PCoA".
#' @param dims [Optional] An integer for the number of dimensions
#' for to keep. Default is 2.
#' @param ... other parameters for Rtsne::Rtsne() function
#'
#' @importFrom stats cmdscale
#' @return A data frame with the low-dimensional
#' represenation.
#'
#' @export
low_dim_vis <- function(D, method = "PCoA", dims = 2, ...){
  if(is.null(rownames(D))) {
    if (is.null(colnames(D))) {
      rownames(D) <- colnames(D) <- 1:nrow(D)
    } else {
      rownames(D) <- colnames(D)
    }
  }
  method <- match.arg(method, c("PCoA", "tSNE"))
  if(!(method %in% c("PCoA", "tSNE"))) {
    stop("Only PCoA and tSNE visualizations are now supported.")
  }
  if (dims > ncol(D)) {
    stop("dims cannot be larger than the number of observations.")
  }
  eigs_percent <- NULL
  if (method == "PCoA") {
    pcoa_res <- cmdscale(sqrt(D), k = dims, eig = TRUE)
    eigs_percent <- signif(100 * pcoa_res$eig /sum(pcoa_res$eig), digits = 3)
    Y <- data.frame(pcoa_res$points)
  } else {
    perpl <- 30
    if (nrow(D) - 1 < 3 * perpl) perpl <- floor((nrow(D) - 1)/3)
    tsne_res <- Rtsne::Rtsne(D, dims = dims, is_distance = TRUE, pca = FALSE,
                             perplexity = perpl, ...)
    Y <- data.frame(tsne_res$Y)
  }
  colnames(Y) <- paste0("D", 1:dims)
  return(list(Y=Y, eigs = eigs_percent))
}

#' Posterior samples of noisy dissimilarities
#'
#' Generate copies of noisy dissimilarity matrices
#' according to the estimates of latent coordinates
#' and the BUDS model.
#'
#' @param D A pairwise dissimilarity matrix.
#' @param buds_fit A stanfit object from fit_buds() output.
#' @param B An integer for the number of copies of D to
#' generate.
#' @param min_sigma A numeric parameter for minimum standard deviation
#' for the distance.
#'
#' @importFrom stats dist rgamma
#' @return A list of two list: one with posterior draws of noisy
#' dissimilarities, booD, and one of corresponding posterior
#' tau estimates bootTau.
#'
#' @export
get_D_copies <- function(D, buds_fit, B, min_sigma = 0.03) {
  fitParams <- rstan::extract(buds_fit$fit_buds)
  n <- dim(fitParams$tau)[2]
  nDraws <- dim(fitParams$tau)[1]
  if(B > nDraws) {
    B <- nDraws
  }
  rel_sd <- sqrt(buds_fit$distDF$v/mean(buds_fit$distDF$v))

  D.lst <- lapply(sample(1:nDraws, size = B, replace = FALSE), function(i) {
    delta.mat <- fitParams$bias[i] + fitParams$rho[i] *
      as.matrix(dist(fitParams$tau[i, ]))
    idelta <- delta.mat[lower.tri(delta.mat)]
    sig_sq <- (pmax(fitParams$meansd[i] * rel_sd, min_sigma))^2
    alpha <- idelta^2/sig_sq
    beta <- idelta/sig_sq
    idvec <- sapply(1:length(idelta), function(k)
      rgamma(1, shape = alpha[k], rate = beta[k]))
    iD <- matrix(0, nrow = n, ncol = n)
    iD[lower.tri(iD)] <- idvec
    iD <- iD + t(iD);
    return(iD)
  })
  tau.lst <- as.list(data.frame(t(fitParams$tau[1:B, ])))
  return(list(D.lst = D.lst, tau.lst = tau.lst))
}


#' Get data in a format for DiSTATIS
#'
#' @param D A pairwise dissimilarity matrix.
#' @param D.lst A list of posterior noisy dissimilarity draws.
#' @param tau_mode [Optional] A mode estimate of tau.
#' @param tau.lst  [Optional] A list of tau posterior draws
#' correspodning to each element of D.lst.
#' @param sample_data  [Optional] A data frame with sample data
#'
#' @return A list with a 3D dissimilarity array, bootD,
#' and a list of data frames booData for each
#' slice of bootD and a data frame modeData
#' @export
get_input_for_distatis <- function(D, D.lst,
                                   tau_mode = NA,
                                   tau.lst = NULL,
                                   sample_data = NA) {
  if(nrow(D) != ncol(D)) {
    stop("D must be a square matrix.")
  }
  if(!is.null(tau.lst) & (length(D.lst) != length(tau.lst))) {
    stop("Length of tau.lst and D.lst must match.")
  }
  if(nrow(sample_data) != nrow(D) | length(tau_mode) != ncol(D)) {
    stop("nrow(sample_data), length(tau_mode) and nrow(D) must match.")
  }
  B <- length(D.lst)
  n <- ncol(D)
  if(is.null(tau.lst)) tau.lst <- rep(NA, B)
  # 3D array for DiSTATIS
  D.arr <- array(0, dim=c(n, n, B + 1))
  D.arr[, , 1] <- D
  for(i in 1:B){
    D.arr[, , i+1] <- D.lst[[i]]
  }
  # Corresponding tau data
  data.lst <- lapply(1:B, function(i) {
   data.frame(tau = tau.lst[[i]],
              rank_tau = rank(tau.lst[[i]], ties.method = "first"),
              sample_data)
  })
  mode_data <- data.frame(tau = tau_mode,
                          rank_tau = rank(tau_mode, ties.method = "first"),
                          sample_data)
  data.lst <- c(list(mode_data), data.lst)
  return(list(bootD = D.arr, booData.lst = data.lst, modeData = mode_data))
}


#' Distatis ordination coordiantes
#'
#' Compute ordination coordinates for all copies
#' of D using a three-way MDS method, DiSTATIS.
#'
#' @param bootD A 3D array with posterior dissimilarity samples.
#' @param dims An integer for the number of dimensions for
#' DiSTATIS.
#' @param booData.lst [Optional] A list of length equal dim(bootD)[3]
#' containing data frames corresponding to bootD slices.
#' @param modeData [Optional] A sample data frame for consensus
#' configuration.
#'
#' @return A list with two data frames: partial for DiSTATIS
#' configurations for each slice of bootD separately,
#' and consensus with an "average" configuration.
#'
#' @export
run_distatis <- function(bootD, dims = 2,
                         booData.lst = NULL,
                         modeData = NULL) {
  if(!is.null(booData.lst) & (length(booData.lst) != dim(bootD)[3])) {
    stop("Length of booData.lst and dim(bootD)[3] must match.")
  }
  if(dim(bootD)[1] != dim(bootD)[2]) {
    stop("bootD must be a 3D array [n x n x B] contating B
         square matrices.")
  }
  # Number of data points
  n <- dim(bootD)[1]
  B <- dim(bootD)[3]

  # Compute DiSTATIS representation
  fit_distatis <- DistatisR::distatis(bootD, nfact2keep = dims)
  # Collect configurations for each D in bootD and for
  distatis_coords <- fit_distatis$res4Splus[["PartialF"]]
  distatis_df <- plyr::alply(distatis_coords, 3)
  distatis_df <- lapply(1:length(distatis_df), function(i) {
    data.frame(SampleID = 1:n, booData.lst[[i]], distatis_df[[i]])
  })
  names(distatis_df) <- 1:length(distatis_df) - 1
  distatis_df <- plyr::ldply(distatis_df)

  # Consensus coordinates
  consensus_coords <- fit_distatis$res4Splus[["F"]]
  consensus_df <- data.frame(SampleID = 1:n,
                             modeData, consensus_coords)
  return(list(partial = distatis_df, consensus = consensus_df))
}
