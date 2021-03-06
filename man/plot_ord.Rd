% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_utils.R
\name{plot_ord}
\alias{plot_ord}
\title{Plot PCoA}
\usage{
plot_ord(ord_res, ord_eig = NULL, size = 1, colData = NULL,
  colLabel = "Variable", title = "Ordination plot",
  prin_curve = FALSE, edgesCol = "grey57", pathCol = "#2171B5",
  lwd = 1.5, ...)
}
\arguments{
\item{ord_res}{rordination coordinates}

\item{ord_eig}{eigenvalues of the ordination}

\item{size}{size of the segments connecting datapoints to the principal curve}

\item{colData}{[Optional] covariates to color the points}

\item{colLabel}{[Optional] label for the covariate}

\item{title}{[Optional] title of the plot}

\item{prin_curve}{[Optional] boolean whether to include a principal curve}

\item{edgesCol}{[Optional] color of the segments. Default "grey".}

\item{pathCol}{[Optional] color of the path. Default "#2171B5".}

\item{lwd}{[Optional] path line width}

\item{...}{other graphical parameters}
}
\value{
A ggplot2 object.
}
\description{
Plot a Principal Coordinate Analysis coordinates
}
