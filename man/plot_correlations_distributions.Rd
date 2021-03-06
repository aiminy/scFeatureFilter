% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_correlations_distributions}
\alias{plot_correlations_distributions}
\title{Produce a density plot of correlation values for each window of feature}
\usage{
plot_correlations_distributions(df, metrics = NULL, vlines = c("mean",
  "median"), facet_ncol = 4)
}
\arguments{
\item{df}{A \code{tibble}, usually the output of \code{\link{correlations_to_densities}}.}

\item{metrics}{Optional. The output of \code{\link{get_mean_median}}. Dashed line will represent
mean or median of the correlation coefficient distributions.}

\item{vlines}{A string, either "mean" or "median". Should the dashed line represent the mean or the median
of the correlation coefficient distributions? Ignored if \code{metrics} is \code{NULL}.}

\item{facet_ncol}{The number of columns to arrange the plots.}
}
\value{
A ggplot2 plot.
}
\description{
Feature by feature correlation values between every windows and the reference to window
of features are visualized as density lines, one facet per comparison.
Two density lines are drown in each facets:
\itemize{
     \item A thin colored line, the correlations between the bin and the reference top bin of features
     \item A thicker blue line with grey error area, the correlations between the bin and the \strong{randomized}
     top bin of features. The lines are not shown if \code{n_random = 0} in \code{\link{correlate_windows}}.
}
}
\examples{
library(magrittr)
myData <- scData_hESC \%>\%
calculate_cvs \%>\%
    define_top_genes(window_size = 100) \%>\%
    bin_scdata(window_size = 1000)

corDistrib <- correlate_windows(myData, n_random = 3)

corDens <- correlations_to_densities(corDistrib)

plot_correlations_distributions(corDens)

metrics <- get_mean_median(corDistrib)

plot_correlations_distributions(corDens, metrics = metrics)

}
\seealso{
\code{\link{correlations_to_densities}}, \code{\link{get_mean_median}}
}
