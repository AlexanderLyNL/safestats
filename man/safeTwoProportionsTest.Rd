% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safeProportionsTest.R
\name{safeTwoProportionsTest}
\alias{safeTwoProportionsTest}
\title{Perform a safe test for two proportions}
\usage{
safeTwoProportionsTest(x, y = NULL, testDesign)
}
\arguments{
\item{x}{either a two-dimensional contingency table in matrix form,
or a factor object indicating the groups (equivalent to \code{fisher.test}).}

\item{y}{a factor object indicating the observations; ignored if x is a matrix.
default \code{NULL}.}

\item{testDesign}{a safe test design for two proportions retrieved through \code{\link{designPilotSafeTwoProportions}} or
\code{\link{designSafeTwoProportions}}.}
}
\value{
safe2x2 test result object with the data used for the test, the S-value and the corresponding
safe p-value.
}
\description{
Perform a safe test for two proportions (a 2x2 contingency table test) with a
result object retrieved through one of design functions for two proportions
in this package, \code{\link{designPilotSafeTwoProportions}} or
\code{\link{designSafeTwoProportions}}.
}
\examples{
testDesign <- designPilotSafeTwoProportions(na = 10, nb = 10)
#test with two factor ojects
groups <- factor(c(rep("a", 10), rep("b", 10)))
observed_data <- factor(c(rep("a", 7), rep("b", 3),rep("a", 3), rep("b", 7)))
safeTwoProportionsTest(x = groups, y = observed_data, testDesign = testDesign)

#test with a table
dataTable <- table(groups, observed_data)
safeTwoProportionsTest(x = dataTable, testDesign = testDesign)
}