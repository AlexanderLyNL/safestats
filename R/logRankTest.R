#' Safe Logrank Test
#'
#' Tests the equality of the survival distributions in two independent groups.
#'
#' The functions safeLogrankTest are modelled after \code{\link[coin]{logrank_test}} and consists of two
#' S3 methods, which both make use of the \code{\link{safeLogrankTestCore}}. Details of the arguments are
#' provided in \code{\link[coin]{logrank_test}}.
#'
#' @param object either a formula with as outcome variable an object resulting from \code{\link[survival]{Surv}}, or
#' and object of class "IndependenceProblem" as described in \code{\link[coin]{logrank_test}}.
#' @param designObj a safe logrank design obtained from \code{\link{designSafeLogrank}}.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the number of
#' samples is exactly as planned.
#' @param alpha numeric representing the tolerable type I error rate. This also serves as a decision rule and it was
#' shown that for safe tests S we have P(S > 1/alpha) < alpha under the null.
#' @param ... further arguments to be passed to \code{\link[coin]{logrank_test}}, such as the "ties.method", "type"
#' and distribution.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the z-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{sValue}{the s-value for the safe test.}
#'   \item{confInt}{To be implemented: a safe confidence interval for the hazard ratio.}
#'   \item{estimate}{To be implemented: an estimate of the hazard ratio.}
#'   \item{h0}{the specified hypothesised value of hazard ratio.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user. Currently, only "two.sided".}
#'   \item{testType}{"logrank".}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" obtained from \code{\link{designSafeLogrank}}.}
#'   \item{logrankObj}{an object obtained from \code{\link[coin]{logrank_test}}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Examples taken from coin::logrank_test
#' ## Example data (Callaert, 2003, Tab. 1)
#' callaert <- data.frame(
#' time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
#' group = factor(rep(0:1, c(7, 8)))
#' )
#'
#' designObj <- designSafeLogrank(nPlan=89)
#' ## Data based on exact logrank test using mid-ranks (p = 0.0505)
#' safeLogrankTest(survival::Surv(time) ~ group, data = callaert,
#'                 designObj = designObj)
#' #'
#' ## Data based on exact logrank test using average-scores
#' safeLogrankTest(survival::Surv(time) ~ group, data = callaert,
#'                 ties.method = "average-scores",
#'                 designObj = designObj)
#' }
safeLogrankTest <- function(object, designObj=NULL, pilot=FALSE, alpha=NULL, ...) {
  if (isFALSE(pilot) && is.null(designObj))
    stop("Please provide a safe logrank design object, or run the function with pilot=TRUE. ",
         "A design object results can be acquired from designSafeLogrank().")

  if (isFALSE(pilot) && designObj[["testType"]] != "logrank")
    stop("The design is constructed for logrank tests, please use designSafeLogrank() for this.")

  if (!is.null(designObj) && !is.null(alpha)) {
    warning("Both designObj and alpha given. The alpha in designObj is used, and the given alpha is ignored")
    alpha <- NULL
  }

  logrankObj <- try(
    coin::logrank_test(object, ...)
  )

  safeLogrankTestCore("logrankObj"=logrankObj, "designObj"=designObj, "pilot"=pilot, "alpha"=alpha)
}

#' Core Function of safeLogrankTest
#'
#' Takes as input an object obtained from \code{\link[coin]{logrank_test}} and a design obtained from
#' \code{\link{designSafeLogrank}} to output an object of class "safeTest".
#'
#' @inherit safeLogrankTest
#' @param logrankObj a logrank object obtained from \code{\link[coin]{logrank_test}}.
#' @param ... further arguments to be passed to or from methods.
#'
#' @export
#'
safeLogrankTestCore <- function(logrankObj, designObj=NULL, pilot=FALSE, alpha=NULL, ...) {
  if (!inherits(logrankObj, "ScalarIndependenceTest"))
    stop("The provided logrankObj is not of the right type derived from coin::logrank_test()")

  groupLabel <- names(logrankObj@statistic@x)

  groupLevels <- levels(logrankObj@statistic@x[[groupLabel]])

  alternative <- logrankObj@statistic@alternative

  if (alternative != "two.sided")
    stop("One-sided tests not yet implemented")

  if (length(groupLevels) > 2)
    stop("K-sample log rank test not yet implemented")

  zStat <- unname(logrankObj@statistic@standardizedlinearstatistic)
  names(zStat) <- "z"

  nEvents <- sum(logrankObj@statistic@ytrans < 0)
  names(nEvents) <- "nEvents"

  dataName <- paste0(names(logrankObj@statistic@y), " by ",
                     names(logrankObj@statistic@x), " (",
                     paste(groupLevels, collapse=", "), ")")

  # result <- list("statistic"=NULL, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
  #                "alternative"=alternative, "testType"=NULL, "dataName"=NULL, "mu0"=mu0, "sigma"=sigma)

  h0 <- 1
  names(h0) <- "theta"

  result <- list("statistic"=zStat, "n"=nEvents, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
                 "h0"=h0, "alternative"=alternative, "testType"="logrank", "dataName"=dataName)
  class(result) <- "safeTest"

  if (isTRUE(pilot)) {

    if (is.null(alpha))
      alpha <- 0.05

    designObj <- designSafeLogrank("nPlan"=nEvents, "alpha"=alpha)
    designObj[["pilot"]] <- TRUE
  }

  sValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=nEvents,
                          "n2"=NULL, "alternative"="two.sided", "paired"=FALSE, "sigma"=1)

  result[["sValue"]] <- sValue
  result[["designObj"]] <- designObj
  result[["logrankObj"]] <- logrankObj

  return(result)
}

#' Designs a Safe Logrank Test
#'
#' Designs a safe logrank test experiment for a prespecified tolerable type I error based on planned sample size(s),
#' which are fixed ahead of time. Outputs a list that includes thetaS, i.e., the safe test defining parameter.
#' Computations exploits the asymptotic normality results of the sampling distribution.
#'
#' @inheritParams designSafeZ
#' @param nPlan numeric > 0, targetted number of events.
#' @param h0 a number indicating the hypothesised true value of the hazard ratio under the null, i.e., h0=1.
#' @param thetaMin numeric that defines the minimal relevant hazard ratio, the smallest hazard ratio that we want to
#' detect.
#' @param zApprox logical, default TRUE to use the asymptotic normality results.
#'
#' @return Returns a safeDesign object that includes:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size either (1) specified by the user, or (2) computed based on beta and thetaMin.}
#'   \item{parameter}{the parameter that defines the safe test. Here log(thetaS).}
#'   \item{thetaMin}{the minimally clinically relevant hazard ratio specified by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{"logrank".}
#'   \item{ratio}{default is 1.}
#'   \item{pilot}{\code{FALSE} to indicate that the design is not a pilot study.}
#'   \item{call}{the expression with which this function is called}
#' }
#'
#' @export
#'
#' @examples
#' designSafeLogrank(nPlan=89)
designSafeLogrank <- function(hazardRatioMin=NULL, alpha=0.05, beta=0.2, nEvents=NULL,
                              alternative=c("two.sided", "greater", "less"), h0=1,
                              ratio=1, zApprox=TRUE, tol=1e-5, ...) {
  stopifnot(0 < alpha, alpha < 1)

  alternative <- match.arg(alternative)

  if (zApprox) {
    if (!is.null(hazardRatioMin)) {
      logHazardRatio <-if (alternative=="two.sided") abs(log(hazardRatioMin)) else log(hazardRatioMin)
      meanDiffMin <- sqrt(ratio)/(1+ratio)*logHazardRatio
    } else {
      logHazardRatio <- NULL
      meanDiffMin <- NULL
    }

    safeZObj <- designSafeZ("meanDiffMin"= meanDiffMin, "alpha"=alpha, "beta"=beta,
                            "nPlan"=nEvents, "alternative"=alternative,
                            "h0"=log(h0), "sigma"=1, "testType"="oneSample")
    nEvents <- safeZObj[["nPlan"]]
    safeZObj[["nPlan"]] <- NULL
    safeZObj[["nEvents"]] <- nEvents
    names(safeZObj[["nEvents"]]) <- "nEvents"
    names(safeZObj[["parameter"]]) <- "log(thetaS)"

    safeZObj[["esMin"]] <- logHazardRatio

    if (!is.null(safeZObj[["esMin"]])) {
      names(safeZObj[["esMin"]]) <- switch(alternative,
                                           "two.sided"="log hazard difference at least abs(log(theta))",
                                           "greater"="log hazard ratio at least",
                                           "less"="log hazard ratio less than")
    }

    safeZObj[["testType"]] <- "logrank"
    safeZObj[["paired"]] <- NULL
    names(safeZObj[["h0"]]) <- "log hazard ratio"
    safeZObj[["call"]] <- sys.call()
  }
  result <- safeZObj

  return(result)
}
