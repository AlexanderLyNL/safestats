% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safeS3Methods.R
\name{print.safeTest}
\alias{print.safeTest}
\title{Print Method for Safe Tests}
\usage{
\method{print}{safeTest}(x, digits = getOption("digits"), prefix = "\\t", ...)
}
\arguments{
\item{x}{a safeTest object.}

\item{digits}{number of significant digits to be used.}

\item{prefix}{string, passed to strwrap for displaying the method components.}

\item{...}{further arguments to be passed to or from methods.}
}
\value{
No returned value, called for side effects.
}
\description{
Printing objects of class 'safeTest' modelled after \code{\link[stats]{print.htest}()}.
}
\examples{
safeTTest(rnorm(19), pilot=TRUE)
safeZTest(rnorm(19), pilot=TRUE)
}
