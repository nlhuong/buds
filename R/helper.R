#' Indices for k-nearest neighbours.
#' 
#' Function returns a n x K data.frame with every row
#' storing indices of the K nearest neighbour of 
#' a corresponding datapoint.
#' 
#' @param D A dissimilarity matrix.
#' @param K An integer indicating the number of points 
#' to include in the neighbourhood. If not specified 
#' one tenth of the total number of data points is used.
#' 
#' @return A data.frame with.
kNN_idx <- function(D, K) {
  idx <- 1:ncol(D)
  kNN.df <- apply(D, 2, function(idist) {
    idist_order <- order(idist)
    return(idx[idist_order[1:K]])
  })
  kNN.df <- t(kNN.df)
  rownames(kNN.df) <- rownames(D)
  colnames(kNN.df) <- paste0("k.", 0:(K-1))
  return(kNN.df)
}

#' Convert dissimilarity data into input format for BUDS
#'
#' Function returns distances in a long data.frame format,
#' which is an input to fit_buds function, and for obtaining 
#' estimates for distances variances.
#' 
#' For each distance  \eqn{d_{ij} = d(x_i, x_j)} the function computes 
#' the mean and variance defined as the empirical mean and variance 
#' of distances in a set   \eqn{\{d(x_i, x_k): x_k \in kNN of x_j\} 
#' \cup \{d(x_j, x_l): x_l \in kNN of x_i \}}.
#'
#' @param D A square matrix of pairwise dissimilarities
#' @param K  An integer indicating k-nearest neighbours for 
#' computing average distance to kNN.
#'
#' @return A data frame to be used as input data for BUDS.
get_dist_df <- function(D, K){
  kNN.df <- kNN_idx(D, K)
  rownames(D) <- colnames(D) <- 1:ncol(D)
  dist.df <- reshape2::melt(D)
  colnames(dist.df) <- c("i", "j", "d")
  dist.df <- dist.df[dist.df$i < dist.df$j, ]
  mean_var_d <- sapply(1:nrow(dist.df), function(k) {
    i <- dist.df[k, "i"]
    j <- dist.df[k, "j"]
    dist_i_to_jth_kNN <- D[i, kNN.df[j, ]]
    dist_j_to_ith_kNN <- D[j, kNN.df[i, ]]
    dK <- c(dist_i_to_jth_kNN, dist_j_to_ith_kNN)
    # The following is for cases where i is in j's kNN set and vice versa
    dK <- dK[dK > 0] 
    m <- mean(dK)
    v <- var(dK)
    return(c(m = m, v = v))
  })
  dist.df <- cbind(dist.df, t(mean_var_d))
  return(dist.df)
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
