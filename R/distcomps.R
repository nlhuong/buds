#' Compute pairwise correlation-based distance
#'
#' Pearson correlation distance for continuous
#' data: d = (1- cor(x_i, x_j))/2 where cor(x_i, x_j)
#' is the correlation between vector x_i and vector x_j.
#'
#' @param X A data matrix, e.g. gene expression
#' @param method a character string indicating which correlation coefficient
#' is to be computed. One of "pearson" (default), "kendall", or
#' "spearman": can be abbreviated.
#' @param scale A boolean indicating whether to normalize the
#' columns (samples) of the data to the even sum.
#' @param base A numeric value for the shared column sum, if scale is TRUE.
#' @param log_trans A boolean indicating whether to log transform
#' the data prio to distance computation (log(X + 1)). Default is FALSE.
#' @param log_base A number indicating base for log transformation.
#' Default is 10.
#'
#' @importFrom stats cor
#' @return A dissimilarity matrix, D.
#' @export
cor_dist <- function(X, method = "pearson",
                     scale = TRUE, base = 1e6,
                     log_trans = TRUE, log_base = 10){
  if(scale) { # scale samples to the same base level
    X <- apply(X, 2, function(x) x/sum(x)*base)
  }
  if (log_trans) { # log-transform data
    X <- log(X + 1, base = log_base)
  }
  # Compute correlation based distance
  corDist <- (1 - cor(X, method = method))/2
  D <- as.matrix(corDist)
  return(D)
}

#' Compute pairwise generic distances
#'
#' A wrapper for dist function.
#'
#' @param X A data matrix, e.g. gene expression
#' @param method a character string indicating which distance method to use.
#' @param log_trans A boolean indicating whether to log transform
#' the data prio to distance computation (log(X + 1)). Default is FALSE.
#' @param log_base A number indicating base for log transformation.
#' Default is 10.
#' @param min_row_prevalence An integer indicating the minimum
#' prevalence (non-zero occurance) of a feature (species) across
#' samples to be kept for distance computation. Only used for Jaccard distance.
#' @param min_row_sum A scalar indicating the minimum feature (species)
#' summ across columns (samples)  to be kept for distance computation.
#' Only used for Jaccard distance.
#'
#' @importFrom stats dist
#' @return A dissimilarity matrix, D.
#' @export
generic_dist <- function(X,
                         method = c("manhattan", "euclidean", "exp manhattan",
                                    "exp euclidean", "maximum", "canberra",
                                     "binary", "minkowski", "jaccard"),
                         log_trans = TRUE, log_base = 10,
                         min_row_sum = 100, min_row_prevalence = 5){
  method <- match.arg(method, c("manhattan", "euclidean", "exp manhattan",
                                "exp euclidean", "maximum", "canberra",
                                "binary", "minkowski", "jaccard"))
  if(method == "jaccard") {
    X <- X[rowSums(X > 0) >= min_row_prevalence, ]
    X <- X[rowSums(X) >= min_row_sum, ]
    method <- "binary"
  }
  if (log_trans){
    if(method != "binary") {
      X <- log(X + 1, base = log_base)
    }
  }
  # Compute pairwise distances
  D <- as.matrix(dist(t(X), method = gsub("exp ", "", method)))
  if (grepl("exp", method)) D <- D/nrow(X)
  if (grepl("exp", method)) {
    D <- 1 - exp(-D)
  }
  return(D)
}

#' Rank based transform distances to triangular distribution
#'
#' Transforms the dissimularities to follow a triangular
#' distribution on [0, 1] interval to better match distances
#' on scalar coordinates that are approximately uniformly
#' distributed.
#'
#' @param D A dissimilarity matrix.
#' @param threshold A boolean whether threshold ranking should be used.
#' @param increment A numeric scalar for increments for distance ranking.
#'
#' @return A transformed dissimilarity matrix, D.
#' @export
transform_dist <- function(D, threshold = FALSE, increment = NULL) {
  dvec0 <- D[lower.tri(D)]
  if (threshold) {
    if(is.null(increment)) {
      dist_range <- diff(range(dvec0))
      increment <- 10*dist_range/length(dvec0)
    }
    dvec <- floor(dvec0/increment)
    dvec <- factor(dvec, levels = sort(unique(dvec)),
                   labels = 1:length(unique(dvec)))
    dvec <- as.numeric(dvec)
  } else {
    dvec <- rank(dvec0, ties.method = "min")
  }
  dvec <- dvec/max(dvec)
  dvec <- 1 - sqrt(1 - dvec)
  D <- matrix(0, nrow = nrow(D), ncol = ncol(D))
  D[lower.tri(D)] <- dvec
  D <- D + t(D)
  return(D)
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
