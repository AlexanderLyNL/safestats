% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tTest.R
\name{computeEsMinSafeT}
\alias{computeEsMinSafeT}
\title{Helper function: Computes the minimal clinically relevant standardised mean difference for the safe t-test
nPlan and beta.}
\usage{
computeEsMinSafeT(
  nPlan,
  alpha = 0.05,
  beta = 0.2,
  alternative = c("twoSided", "greater", "less"),
  testType = c("oneSample", "paired", "twoSample"),
  lowN = 3,
  highN = 1e+06,
  ratio = 1
)
}
\arguments{
\item{nPlan}{vector of max length 2 representing the planned sample sizes.}

\item{alpha}{numeric in (0, 1) that specifies the tolerable type I error control --independent of n-- that the
designed test has to adhere to. Note that it also defines the rejection rule e10 > 1/alpha.}

\item{beta}{numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both
the sample sizes and deltaS, which defines the test. Note that 1-beta defines the power.}

\item{alternative}{a character string specifying the alternative hypothesis must be one of "twoSided" (default),
"greater" or "less".}

\item{testType}{either one of "oneSample", "paired", "twoSample".}

\item{lowN}{integer minimal sample size of the (first) sample when computing the power due to
optional stopping. Default lowN is set 1.}

\item{highN}{integer minimal sample size of the (first) sample when computing the power due to
optional stopping. Default highN is set 1e6.}

\item{ratio}{numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.}
}
\value{
a list which contains at least nPlan and the phiS the parameter that defines the safe test
}
\description{
Helper function: Computes the minimal clinically relevant standardised mean difference for the safe t-test
nPlan and beta.
}
