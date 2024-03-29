% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safe2x2Test.R
\name{computeConfidenceBoundsForDifferenceTwoProportions}
\alias{computeConfidenceBoundsForDifferenceTwoProportions}
\title{Estimate Lower and Upper Bounds on the Confidence Sequence (Interval)
for the Difference Divergence Measure for Two Proportions}
\usage{
computeConfidenceBoundsForDifferenceTwoProportions(
  ya,
  yb,
  precision,
  safeDesign
)
}
\arguments{
\item{ya}{positive observations/ events per data block in group a: a numeric with integer values
between (and including) 0 and \code{na}, the number of observations in group a per block.}

\item{yb}{positive observations/ events per data block in group b: a numeric with integer values
between (and including) 0 and \code{nb}, the number of observations in group b per block.}

\item{precision}{precision of the grid of differences to search over for the lower and upper bounds.}

\item{safeDesign}{a 'safeDesign' object obtained through
\code{\link{designSafeTwoProportions}}}
}
\value{
list with found lower and upper bound.
}
\description{
Estimate Lower and Upper Bounds on the Confidence Sequence (Interval)
for the Difference Divergence Measure for Two Proportions
}
\examples{
balancedSafeDesign <- designSafeTwoProportions(na = 1,
                                               nb = 1,
                                               nBlocksPlan = 10,
                                               alpha = 0.05)
ya <- c(1,1,1,1,1,1,1,1,0,1)
yb <- c(0,0,0,0,1,0,0,0,0,0)
computeConfidenceBoundsForDifferenceTwoProportions(ya = ya,
                                               yb = yb,
                                               precision = 20,
                                               safeDesign = balancedSafeDesign)

}
