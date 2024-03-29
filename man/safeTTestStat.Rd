% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tTest.R
\name{safeTTestStat}
\alias{safeTTestStat}
\title{Computes E-Values Based on the T-Statistic}
\usage{
safeTTestStat(
  t,
  parameter,
  n1,
  n2 = NULL,
  alternative = c("twoSided", "less", "greater"),
  tDensity = FALSE,
  paired = FALSE,
  ...
)
}
\arguments{
\item{t}{numeric that represents the observed t-statistic.}

\item{parameter}{numeric this defines the safe test S, i.e., a likelihood ratio of t distributions with in
the denominator the likelihood with delta = 0 and in the numerator an average likelihood defined by
1/2 time the likelihood at the non-centrality parameter sqrt(nEff)*parameter and 1/2 times the likelihood at
the non-centrality parameter -sqrt(nEff)*parameter.}

\item{n1}{integer that represents the size in a one-sample t-test, (n2=\code{NULL}). When n2 is not \code{NULL},
this specifies the size of the first sample for a two-sample test.}

\item{n2}{an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
\code{NULL} it implies that the t-statistic is based on one-sample.}

\item{alternative}{a character only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE},
then the alternative specified by the design object is used instead.}

\item{tDensity}{Uses the the representation of the safe t-test as the likelihood ratio of t densities.}

\item{paired}{a logical indicating whether you want a paired t-test.}

\item{...}{further arguments to be passed to or from methods.}
}
\value{
Returns a numeric that represent the e10, that is, the e-value in favour of the alternative over the null
}
\description{
A summary stats version of \code{\link{safeTTest}()} with the data replaced by t, n1 and n2, and the
design object by deltaS.
}
\examples{
safeTTestStat(t=1, n1=100, 0.4)
safeTTestStat(t=3, n1=100, parameter=0.3)
}
