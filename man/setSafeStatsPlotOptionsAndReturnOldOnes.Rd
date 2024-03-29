% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{setSafeStatsPlotOptionsAndReturnOldOnes}
\alias{setSafeStatsPlotOptionsAndReturnOldOnes}
\title{Sets 'safestats' Plot Options and Returns the Current Plot Options.}
\usage{
setSafeStatsPlotOptionsAndReturnOldOnes(...)
}
\arguments{
\item{...}{further arguments to be passed to or from methods.}
}
\value{
Returns a list with the user specified plot options.
}
\description{
Sets 'safestats' Plot Options and Returns the Current Plot Options.
}
\examples{
oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
graphics::plot(1:10, 1:10)
setPar <- graphics::par(oldPar)
}
