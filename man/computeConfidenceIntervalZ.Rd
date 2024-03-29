% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zTest.R
\name{computeConfidenceIntervalZ}
\alias{computeConfidenceIntervalZ}
\title{Helper function: Computes the safe confidence sequence for a z-test}
\usage{
computeConfidenceIntervalZ(
  nEff,
  meanObs,
  phiS,
  sigma = 1,
  ciValue = 0.95,
  alternative = "twoSided",
  a = NULL,
  g = NULL
)
}
\arguments{
\item{nEff}{numeric > 0, the effective sample size.}

\item{meanObs}{numeric, the observed mean. For two sample tests this is difference of the means.}

\item{phiS}{numeric > 0, the safe test defining parameter.}

\item{sigma}{numeric > 0 representing the assumed population standard deviation used for the test.}

\item{ciValue}{numeric is the ciValue-level of the confidence sequence. Default ciValue=0.95.}

\item{alternative}{a character string specifying the alternative hypothesis must be one of "twoSided" (default),
"greater" or "less".}

\item{a}{numeric, the centre of the normal prior on population mean (of the normal data). Default
is \code{NULL}, which implies the default choice of setting the centre equal to the null hypothesis.}

\item{g}{numeric > 0, used to define g sigma^2 as the variance of the normal prior on the population
(of the normal data). Default is \code{NULL} in which case g=phiS^2/sigma^2.}
}
\value{
numeric vector that contains the upper and lower bound of the safe confidence sequence
}
\description{
Helper function: Computes the safe confidence sequence for a z-test
}
\examples{
computeConfidenceIntervalZ(nEff=15, meanObs=0.3, phiS=0.2)
}
