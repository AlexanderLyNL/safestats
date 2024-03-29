% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logRankTest.R
\name{logrankSingleZ}
\alias{logrankSingleZ}
\title{Helper function computes single component of the logrank statistic}
\usage{
logrankSingleZ(obs0, obs1, y0, y1, ...)
}
\arguments{
\item{obs0}{integer, number of observations in the control group}

\item{obs1}{integer, number of observations in the treatment group}

\item{y0}{integer, total number of participants in the control group}

\item{y1}{integer, total number of participants in the treatment group}

\item{...}{further arguments to be passed to or from methods.}
}
\value{
Returns a list containing at least the following components:
\describe{
  \item{oMinE}{observed minus expected.}
  \item{v}{hypergeometric variance.}
}
}
\description{
Helper function computes single component of the logrank statistic
}
\examples{
y0Vector <- c(6, 4, 4, 1, 0)
y1Vector <- c(6, 6, 5, 2, 2)
obs0Vector <- c(1, 0, 2, 1, 0)
obs1Vector <- c(0, 1, 1, 0, 1)

varVector <- oMinEVector <-y0Vector

for (i in seq_along(y0Vector)) {
  tempResult <- logrankSingleZ(obs0=obs0Vector[i], obs1=obs1Vector[i],
                              y0=y0Vector[i], y1=y1Vector[i])
  oMinEVector[i] <- tempResult[["oMinE"]]
  varVector[i] <- tempResult[["v"]]
}

sum(oMinEVector)/sqrt(sum(varVector))

}
