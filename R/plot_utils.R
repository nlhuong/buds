## ----- plot functions -----

#' Plot matrix
#'
#' A simple plot of matrix entries values.
#'
#' @param mat [matrix/ data.frame] (Required) Data matrix to plot.
#' @param rowIdx,colIdx [integer] (Optional) Default \code{NULL}.
#'     Indices of the matrix to plot.
#' @param ... (Optional) Additional plotting options.
#' @return p [ggplot] A ggplot heatmap object with nice defaults.
#' @importFrom ggplot2 ggplot geom_tile aes_string scale_x_discrete
#'   scale_fill_gradient scale_y_discrete geom_raster
#' @importFrom magrittr %>%
#' @export
#' @examples
#' plot_matrix(matrix(rnorm(2000)*100, 100, 200))
plot_matrix <- function(mat, rowIdx = NULL, colIdx = NULL, ...) {
    fill <- NULL
    opts <- merge_heatmap_opts(list(...))
    mat <- as.matrix(mat)
    if (!is.null(rowIdx))
        mat <- mat[rowIdx, ]
    if (!is.null(colIdx))
        mat <- mat[, colIdx]
    if(is.null(rownames(mat))) rownames(mat) <- 1:nrow(mat)
    if(is.null(colnames(mat))) colnames(mat) <- 1:ncol(mat)

    mat_df <- mat %>%
        reshape2::melt(
            varnames = c("row", "col"),
            value.name = "fill"
        )
    if (is.null(opts$x_order)) opts$x_order <- colnames(mat)
    if (is.null(opts$y_order)) opts$y_order <- rownames(mat)
    mat_df[, "col"] <- factor(mat_df[, "col"], levels = opts$x_order)
    mat_df[, "row"] <- factor(mat_df[, "row"], levels = opts$y_order)

    p <- ggplot(mat_df) +
        geom_raster(aes(x = col, y = row, fill = fill)) +
        scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
        scale_fill_gradientn(colors = opts$fill_colors, breaks = opts$fill_breaks) +
        coord_fixed(opts$coord_ratio) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0))

    return(p)
}


#' Plot distance to kNN
#'
#' Computes average distance to K nearest neighbours and plots
#' the values for every data point ordered according to
#' coordinates computed with BUDS.
#'
#' @param D A dissimilarity matrix.
#' @param tau A numeric vector of latent 1D coordinates.
#' @param K [Optional] An integer indicating the number of
#' points to include in the neighbourhood. If not specified
#' one tenth of the total number of data points is used.
#' @param color [Optional] A vector or a character string for
#' point colors.
#' @param color_label [Optional] A character string for the color label.
#' @param error_bars [Optional] Whether to include error bars in the plot.
#' @import ggplot2
#' @importFrom  viridis scale_fill_viridis scale_color_viridis
#' @return A ggplot2 object.
#' @export
plot_kNN_mean_dist <- function(D, tau, K = NULL,  color = "blue",
                               color_label = "Covariate",
                               error_bars = FALSE) {
  kNN_dist_mean <- kNN_dist_sd <- col <- NULL
  if (length(tau) != ncol(D)) {
    stop("Length tau must be equal to nrow(D), i.e. the number of observations.")
  }
  dkNN <- kNN_dist(D, K = K)
  df <- data.frame(tau = tau, col = color,
                   kNN_dist_mean = dkNN$mean,
                   kNN_dist_sd = sqrt(dkNN$var))
  plt <- ggplot(df, aes(x = tau, y = kNN_dist_mean, fill = col, color = col)) +
    geom_point(size = 3, color = "grey77", pch = 21) +
    scale_fill_viridis(name = color_label, discrete = !(is.numeric(color))) +
    scale_color_viridis(name = color_label, discrete = !(is.numeric(color))) +
    ylab("Distance to kNN")
  if(error_bars) {
    plt <- plt + geom_errorbar(aes(ymin = kNN_dist_mean - kNN_dist_sd,
                                   ymax = kNN_dist_mean + kNN_dist_sd))
  }
  return(plt)
}



#' Plot 1D latent coordinates
#'
#' Plot 1D latent coordinates against their ranking
#' or supplied sample covariate data.
#'
#' @param tau_df A data frame with columns "tau", "tau_upper",
#' and "tau_lower".
#' @param covariate [Optional] If supplied taus coordiantes
#' are plotted against thissample  covariate vector.
#' @param covariate_name [Optional] The name of the covariat
#' supplied.
#' @param color [Optional] A vector or a character string for
#' point colors.
#' @param color_label [Optional] A character string for the color label.
#' @param idxBigger [Optional] A vector of integers indicating which
#' data point to highlight (make bigger).
#'
#' @import ggplot2
#' @importFrom viridis scale_color_viridis scale_fill_viridis
#' @return A ggplot2 object.
#' @export
plot_buds_1D <- function(tau_df, covariate = NULL, covariate_name = "covariate",
                         color = "blue", color_label = NULL, idxBigger = NULL) {
tau <- tau_lower <- tau_upper <- x <- NULL
  if(!all(c("tau", "tau_lower", "tau_upper") %in% colnames(tau_df))) {
    stop("'tau_df' must have columns c('tau', 'tau_lower', 'tau_upper')
         specified.")
  }
  if(is.null(covariate)){
    tau_df$x <- rank(tau_df$tau)
  } else {
    tau_df$x <- rank(covariate)
  }
  tau_df$color <- color
  name_x_axis <- ifelse(is.null(covariate), "rank(tau)",
                        paste0("rank(", covariate_name, ")"))
  plt <- ggplot(tau_df, aes(x = x, y = tau, color = color)) +
    geom_errorbar(aes(ymin = tau_lower, ymax = tau_upper),
                  lwd = 0.7) +
    geom_point(aes(fill = color), color = "grey77", pch = 21, size = 2) +
    scale_color_viridis(name = color_label, discrete = (!is.numeric(color))) +
    scale_fill_viridis(name = color_label, discrete = (!is.numeric(color))) +
    xlab(name_x_axis) + ylab("tau")

  if(!is.null(idxBigger)) {
    plt <- plt +
      geom_point(data = tau_df[idxBigger, ], aes(fill = color),
                 color = "grey77", pch = 21, size = 3)
  }
  return(plt)
}


#' Plot data trajectory in 2 or3D
#'
#' Plot 2 or 3D representations of the data with PCoA
#' of t-SNE, together with posterior data trajectories
#' (paths) computed with BUDS.
#'
#' @param budsParams A list of parameters extracted from BUDS stanfit object.
#' @param Y A data frame with 2D representation of the data
#' @param eigs [Optional] A vector of eigenvalues corresponding
#' to columns of Y.
#' @param sample_data [Optional] A data frame with sample data.
#' @param covariate_name [Optional] A charcter string for the column
#' of sample_data to use for coloring points.
#' @param path_col [Optional] A charcter string for the color
#' of the posterior paths. Default is "#2171B5".
#' @param nPaths [Optional] An integer for number of
#' posterior paths to include. Default is 50.
#' @param nCenters [Optional] An integer for number of
#' data points to use for the highlighted "mode-path".
#' Default is 50.
#'
#' @import ggplot2
#' @importFrom plotly plot_ly add_trace add_markers layout
#' @importFrom viridis scale_fill_viridis
#' @return A ggplot2 object.
#' @export
plot_buds_trajectory <- function(budsParams, Y, eigs = NULL,
                                 sample_data = NULL,
                                 covariate_name = "covariate",
                                 path_col = "#2171B5", nPaths = 50,
                                 nCenters = 50){
  D1 <- D2 <- trajectory <- NULL
  nCenters <- min(nCenters, nrow(Y))
  color_data <- 1
  if (!is.null(sample_data) & (covariate_name %in% colnames(sample_data))) {
    color_data <- sample_data[, covariate_name]
  }
  # Extract tau parameter
  tau_mcmc <- coda::mcmc(budsParams$tau)
  tau_mode <- MCMCglmm::posterior.mode(tau_mcmc)
  #names(tau_mode) <- rownames(D)
  tau_samples <- budsParams$tau
  tau_samples <- tau_samples[sample(1:nrow(tau_samples), nPaths), ]
  # Data to plot
  DF <- Y
  colnames(DF) <- paste0("D", 1:ncol(Y))
  if(!is.null(sample_data)) DF <- cbind(DF, sample_data,
                                        color_data = color_data)
  modeDF <- data.frame(tau_mode = tau_mode, DF)
  modeDF <- modeDF[order(tau_mode), ]
  # Pick a subset of samples evenly spaced alon the gradient
  idx <- seq(1, length(tau_mode), length.out = nCenters)
  modeDF <- modeDF[idx, ]
  trajDF <- lapply(1:nrow(tau_samples), function(i) {
    data.frame(trajectory = i,  Y[order(tau_samples[i, ]), ])
  })
  trajDF <- do.call("rbind", trajDF)

  if (ncol(Y) == 2) {
    plt <- ggplot(trajDF, aes(D1, D2)) +
      geom_path(aes(group = trajectory), color = path_col, alpha = 0.1) +
      geom_path(data = modeDF, color = "grey17", lwd = 1.2) +
      geom_point(data = DF, color = "grey77", pch = 21, size = 2,
                 aes(fill = color_data)) +
      geom_point(data = modeDF, color = "grey77", pch = 21, size = 3,
                 aes(fill = color_data)) +
      scale_fill_viridis(name = covariate_name,
                         discrete = !(is.numeric(color_data)))
    if (!is.null(eigs)){
      plt <- plt + xlab(paste0("PC1 [", eigs[1], "%]")) +
        ylab(paste0("PC2 [", eigs[2], "%]")) +
        coord_fixed(ratio = max(0.5, eigs[2]/eigs[1]))
    }
  } else if (ncol(Y) > 2) {
    plt <- plot_ly(type = 'scatter3d', mode = 'lines+markers') %>%
      add_trace(x = modeDF$D1, y = modeDF$D2, z = modeDF$D3,
                color = modeDF$color_data, marker = list(size = 3),
                line = list(width = 2)) %>%
      add_markers(x = DF$D1, y = DF$D2, z = DF$D3, marker = list(size = 1),
                  color = DF$color_data, mode = "markers")
    if (!is.null(eigs)){
      axx <- list(
        title = paste0("PC1 [", eigs[1], "%]")
      )

      axy <- list(
        title = paste0("PC2 [", eigs[2], "%]")
      )

      axz <- list(
        title = paste0("PC3 [", eigs[3], "%]")
      )
      plt <- plt %>%
        layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
    }
  }
  return(plt)
}


#' Plot DiSTATIS representation of the data
#'
#' Uses distatis_df data frame with configurations
#' computed by DiSTATIS for separate dissimilarity matrix
#' copies separetely, and a consensus_df data frame
#' for the "agreement" configuration to plot the data
#' representation with density clouds.
#'
#' @param distatis_df A data frame with partial configurations.
#' @param consensus_df A data frame with a consensus configuration.
#' consensus_df can contain a column "covariate" with sample data
#' that can be used for coloring.
#' @param color_label [Optional] A character string for the color label.
#' @import ggplot2
#' @importFrom  viridis scale_fill_viridis
#' @return A ggplot2 object.
#' @export
plot_distatis <- function(distatis_df, consensus_df, color_label = NULL) {

  Factor.1 <- Factor.2 <- .id <- ..level.. <- NULL
  if (is.na(color_label) | !(color_label %in% colnames(consensus_df))) {
    covariate <- rep(1, nrow(consensus_df))
  } else {
    covariate <- consensus_df[, color_label]
  }
  plt <- ggplot(distatis_df, aes(Factor.1, Factor.2)) +
    geom_density2d(bins = 10, lwd = 0.5) +
    stat_density2d(aes(alpha = ..level..), fill = "#2171B5",
                   n = 20, size = 0.01, geom ="polygon", bins = 20) +
    geom_point(data =  distatis_df %>% filter(.id == 0),
               size = 2, color = "grey67") +
    geom_point(data = consensus_df, size = 3, color = "grey17",
               pch = 21, aes(fill = covariate)) +
    scale_fill_viridis(direction = 1, name = color_label,
                       discrete = !is.numeric(covariate)) +
    scale_alpha(range = c(0, 0.9)) +
    xlim(1.1*min(distatis_df$Factor.1), 1.1*max(distatis_df$Factor.1)) +
    ylim(1.1*min(distatis_df$Factor.2), 1.1*max(distatis_df$Factor.2))
  return(plt)
}


#' Plot data point contours
#'
#' Uses distatis_df data frame with configurations
#' computed by DiSTATIS for separate dissimilarity matrix
#' copies separetely, and a consensus_df data frame
#' for the "agreement" configuration to plot confidence
#' contours for points selected by idx_list.
#'
#' @param distatis_df A data frame with partial configurations.
#' @param consensus_df A data frame with a consensus configuration.
#' consensus_df can contain a column "covariate" with sample data
#' that can be used for coloring.
#' @param idx_list An integer vector of selected data point indices.
#' @param color_label [Optional] A character string for the color label.
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom viridis scale_color_viridis
#' @return A ggplot2 object.
#' @export
plot_point_contours <- function(distatis_df, consensus_df,
                                idx_list = NULL, color_label = NULL) {
  .id <- Factor.1 <- Factor.2 <- SampleID <- covariate <- NULL
  if(!("rank_tau" %in% colnames(consensus_df))) {
    stop("consensus_df must have a column rank_tau indicating the ordering of
         the rows.")
  }
  if (!(color_label %in% colnames(consensus_df))) {
    consensus_df$covariate <- 1
  } else {
    consensus_df$covariate <- consensus_df[, color_label]
  }

  if (!(color_label %in% colnames(distatis_df))) {
    distatis_df$covariate <- 1
  } else {
    distatis_df$covariate <- distatis_df[, color_label]
  }
  n <- nrow(consensus_df)
  if(is.null(idx_list)) {
    idx_list <- floor(seq(5, n-5, length.out = 3))
  }
  chosen_samples <- consensus_df$SampleID[consensus_df$rank_tau %in% idx_list]
  plt <- ggplot(distatis_df, aes(Factor.1, Factor.2, group = SampleID)) +
    stat_density2d(data = distatis_df %>% filter(SampleID %in% chosen_samples),
                   aes(color = covariate),
                   geom ="density2d", bins = 10) +
    geom_point(data =  distatis_df %>% filter(.id == 0),
               size = 2, color = "grey67") +
    geom_point(data = consensus_df, size = 3, aes(color = covariate)) +
    scale_color_viridis(direction = 1, name = color_label,
                        discrete = !is.numeric(distatis_df$covariate)) +
    xlim(1.1*min(distatis_df$Factor.1), 1.1*max(distatis_df$Factor.1)) +
    ylim(1.1*min(distatis_df$Factor.2), 1.1*max(distatis_df$Factor.2))
  return(plt)
}


#' Plot features
#'
#' Plot feature trends along the estimated ordering
#' and include smoothing curves.
#'
#' @param X A data frame of matrix.
#' @param tau A vector with latent 1D coordinates for data points.
#' @param feat_idx An integer vector of selected features.
#' @param log_trans [Optional] A boolean for whether to log
#' transform the data before plotting. Default is FALSE.
#' @importFrom stats median cor
#' @import ggplot2
#' @return A ggplot2 object.
#' @export
plot_features_curves <- function(X, tau, feat_idx = NA,
                                 log_trans = FALSE){
  nFeat <- 5
  Value <- Features <- NULL
  if(!all(names(tau) %in% colnames(X))) {
    names(tau) <- colnames(X)
  }
  X <- median(colSums(X)) * apply(X, 2, function(x) x/sum(x))
  if(any(is.na(feat_idx))) {
    X <- X[rowSums(X) > 0, ]
    tau_cor <- apply(X, 1, function(x) cor(tau, x, method = "spearman"))
    X <- X[order(-abs(tau_cor))[1:nFeat], ]
  } else if (any(feat_idx == "random")){
    X <- X[rowSums(X) > 0, ]
    X <- X[sample(1:nrow(X), nFeat), ]
  } else {
    X <- X[feat_idx, ]
  }
  if (log_trans) {
    X <- log10(X + 1)
  }
  X <- as.matrix(X)
  X.m <- reshape2::melt(X)
  colnames(X.m) <- c("Features", "Samples", "Value")
  X.m$tau <- tau[X.m$Samples]
  plt <- ggplot(X.m, aes(tau, Value, color =factor(Features))) +
    geom_smooth(aes(group = Features)) +
    geom_point(alpha = 0.7) +
    ylab("log10(X + 1)")
  return(plt)
}

#' Plot ordered data matrix
#'
#' Plot a heatmap for a data matrix where columns
#' are reordered according to tau and rows according
#' to either byMean or by maximum value in a window
#' along tau.
#'
#'
#' @param X A data frame of matrix where rows are features and columns
#' are samples.
#' @param tau A vector with latent 1D coordinates for data points.
#' @param log_trans [Optional] A boolean for whether to log
#' transform the data before plotting. Default is FALSE.
#' @param norm_samples A boolean whether the columns (samples) should be
#' normalize to even total sum.
#' @param keep_features [Optional] A list of features (rows) that must be
#' included in the heatmap.
#' @param nfeatures The number of features to keep. Default is 500.
#' @param byMean A boolean whether to order features by their "mean tau"
#' value. If FALSE the features are ordered by maximum in a window
#' along tau.
#' @param window An integer for the window size. Default is 5.
#' @param minsum A for the minimum sum for a feature in a data matrix.
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @return A ggplot2 object.
#' @export
plot_ordered_matrix <- function(X, tau, log_trans = FALSE, norm_samples = FALSE,
                                keep_features = NULL, nfeatures = 500,
                                byMean = TRUE, window = 5, minsum = 0){
  X <- X[rowSums(X) > minsum, ]
  nfeatures <- min(nfeatures, nrow(X))
  if(is.null(keep_features)) keep_features <- sample(1:nrow(X), nfeatures)
  nSamples <- ncol(X)
  X2plot <- X[keep_features, order(tau)]
  normX2plot <- apply(X2plot, 2, function(x) x/sum(x))
  if(byMean) {
    feat_loc <- normX2plot %*% tau
  } else {
    feat_loc <- apply(normX2plot, 1, function(x) {
      x <- c(x)
      y <- cumsum(x)[window:nSamples] - c(0, cumsum(x)[1:(nSamples - window)])
      which.max(y)
    })
  }
  if (norm_samples) X2plot <- normX2plot
  X2plot <- X2plot[order(feat_loc), ]
  rownames(X2plot) <- colnames(X2plot) <- NULL
  if(log_trans) X2plot <- as.matrix(X2plot) + 1
  plt <- plot_matrix(X2plot) + coord_fixed() +
    xlab("Samples") + ylab("Features") +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank())
  if(log_trans) {
    plt <- plt + scale_fill_viridis(name = "Abund.", trans = "log10")
  } else {
    plt <- plt + scale_fill_viridis(name = "Abund.")
  }
  return(plt)
}


#' Plot PCoA
#'
#' Plot a Principal Coordinate Analysis coordinates
#'
#' @param ord_res rordination coordinates
#' @param ord_eig eigenvalues of the ordination
#' @param size size of the segments connecting datapoints to the principal curve
#' @param colData [Optional] covariates to color the points
#' @param colLabel [Optional] label for the covariate
#' @param title [Optional] title of the plot
#' @param prin_curve [Optional] boolean whether to include a principal curve
#' @param edgesCol [Optional] color of the segments. Default "grey".
#' @param pathCol [Optional] color of the path. Default "#2171B5".
#' @param lwd [Optional] path line width
#' @param ... other graphical parameters
#'
#' @import ggplot2
#' @importFrom princurve principal_curve
#' @importFrom viridis scale_color_viridis
#' @return A ggplot2 object.
#' @export
plot_ord <- function(ord_res, ord_eig = NULL, size = 1,
                     colData = NULL, colLabel = "Variable",
                     title = "Ordination plot", prin_curve = FALSE,
                     edgesCol = "grey57", pathCol = "#2171B5",
                     lwd = 1.5, ...) {
    X1 <- X2 <- NULL
    if (!is.null(ord_eig)) {
        ord_eig <- 100 * ord_eig /sum(ord_eig)
        ord_eig <- signif(ord_eig, digits = 3)
    }
    X <- data.frame(ord_res)
    colnames(X) <- paste0("X", 1:ncol(X))
    p <- ggplot(X, aes(X1, X2)) + ggtitle(title) +
        coord_fixed(1)

    if(prin_curve) {
        prin_curve <- principal_curve(as.matrix(X), plot_iterations = FALSE, ...)
        fittedLine <- data.frame(prin_curve$s[prin_curve$ord, ])
        p <- p + geom_path(data = fittedLine, col = pathCol, lwd = lwd) +
            geom_segment(aes(xend = prin_curve$s[, 1], yend = prin_curve$s[, 2]),
                         size = 0.5, col = edgesCol)
    }
    if (!is.null(colData)) {
        p <- p +
            geom_point(
                aes(fill = colData),
                pch = 21, color = "grey50", size = size) +
            scale_fill_viridis(
                name = colLabel, discrete = !is.numeric(colData))
    } else {
        p <- p + geom_point(size = size)
    }
    if (!is.null(ord_eig)){
        eig_ratio =  ord_eig[2]/ord_eig[1]
        p <- p + xlab(paste0("PC1 [", ord_eig[1], "%]")) +
            ylab(paste0("PC2 [", ord_eig[2], "%]"))
    }
    return(list(plot = p, fit.prin_curve = prin_curve))
}


#' Merge default options for a heatmap
#' @param opts [list] (Optional) A partially specified list used to customize
#'   ppearance in ggplot theme(). Options that are already specified will not
#'   be changed, those that are not will be filled in with defaults.
#' @return opts [list]  A version of opts with unspecified options filled in
#'   with defaults.
#' @importFrom viridis viridis
#' @importFrom utils modifyList
#' @export
merge_heatmap_opts <- function(opts = list()) {
  default_opts <- list(
    "x" = "col",
    "y" = "row",
    "fill_colors" = viridis(256),
    "fill_breaks" = NULL,
    "facet_terms" = NULL,
    "facet_scales" = "fixed",
    "facet_space" = "fixed",
    "x_order" = NULL,
    "y_order" = NULL,
    "coord_ratio" = 1,
    "theme_opts" = list()
  )
  modifyList(default_opts, opts)
}

