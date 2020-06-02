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
#' @param h0 a number indicating the hypothesised true value of the hazard ratio under the null. Default set to 1
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the number of
#' samples is exactly as planned.
#' @param alpha numeric > 0 only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the alpha of
#' the design object is used instead in constructing the decision rule S > 1/alpha.
#' @param alternative a character only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the
#' alternative specified by the design object is used instead.
#' @param ... further arguments to be passed to \code{\link[coin]{logrank_test}}, such as the "ties.method", "type"
#' and distribution.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the z-statistic.}
#'   \item{nEvents}{The number of observed events.}
#'   \item{sValue}{the s-value of the safe test.}
#'   \item{confSeq}{An anytime-valid confidence sequence.}
#'   \item{estimate}{To be implemented: An estimate of the hazard ratio.}
#'   \item{h0}{the specified hypothesised value of hazard ratio.}
#'   \item{testType}{"logrank".}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" obtained from \code{\link{designSafeLogrank}}.}
#'   \item{logrankObj}{an object obtained from \code{\link[coin]{logrank_test}}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' # Examples taken from coin::logrank_test
#' ## Example data (Callaert, 2003, Tab. 1)
#' callaert <- data.frame(
#' time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
#' group = factor(rep(0:1, c(7, 8)))
#' )
#'
#' designObj <- designSafeLogrank(hrMin=1/2, beta=0.2)
#'
#' ## Data based on logrank test using mid-ranks (p = 0.0505)
#' safeLogrankTest(survival::Surv(time) ~ group, data = callaert,
#'                 designObj = designObj)
#'
#' designObj <- designSafeLogrank(hrMin=1/2)
#' ## Data based on logrank test using average-scores
#' safeLogrankTest(survival::Surv(time) ~ group, data = callaert,
#'                 ties.method = "average-scores",
#'                 designObj = designObj)
safeLogrankTest <- function(object, designObj=NULL, h0=1,
                            pilot=FALSE, alpha=NULL, alternative=NULL, ...) {

  if (isFALSE(pilot) && is.null(designObj))
    stop("Please provide a safe logrank design object, or run the function with pilot=TRUE. ",
         "A design object can be obtained by running designSafeLogrank().")

  if (!is.null(designObj)) {
    if (!is.null(alpha))
      warning("Both a design object and an alpha given. The alpha specified by the design object ",
              "is used for the test, and the provided alpha is ignored.")

    if (!is.null(alternative))
      warning("Both a design object and an alternative given. The alternative specified by ",
              "the design object is used for the test, and the provided alternative is ignored.")

    if (names(designObj[["parameter"]]) != "log(thetaS)")
      warning("The provided design is not constructed for the logrank test,",
              "please use designSafeLogrank() instead. The test results might be invalid.")
  }

  logrankObj <- try(
    coin::logrank_test(object, alternative="two.sided", ...)
  )

  safeLogrankTestCore("logrankObj"=logrankObj, "designObj"=designObj, "pilot"=pilot, "alpha"=alpha, "h0"=h0,
                      "alternative"=alternative)
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
#'
safeLogrankTestCore <- function(logrankObj, designObj=NULL, alternative, h0=1,
                                pilot=FALSE, alpha=NULL, ...) {

  if (!inherits(logrankObj, "ScalarIndependenceTest"))
    stop("The provided logrankObj is not of the right type derived from coin::logrank_test()")

  groupLabel <- names(logrankObj@statistic@x)

  groupLevels <- levels(logrankObj@statistic@x[[groupLabel]])

  if (length(groupLevels) > 2)
    stop("K-sample log rank test not yet implemented")

  zStat <- unname(logrankObj@statistic@standardizedlinearstatistic)
  nEvents <- sum(logrankObj@statistic@ytrans < 0)

  dataName <- paste0(names(logrankObj@statistic@y), " by ",
                     names(logrankObj@statistic@x), " (",
                     paste(groupLevels, collapse=", "), ")")

  if (isTRUE(pilot)) {
    if (is.null(alpha))
      alpha <- 0.05

    if (is.null(alternative)) {
      alternative <- "two.sided"
    } else {
      if (!(alternative %in% c("two.sided", "greater", "less")))
        stop('Provided alternative must be one of "two.sided", "greater", or "less".')
    }

    designObj <- designSafeLogrank("hrMin"=NULL, "beta"=NULL, "nEvents"=nEvents, "alpha"=alpha,
                                   "alternative"=alternative)
    designObj[["pilot"]] <- TRUE
  }

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  ratio <- designObj[["ratio"]]

  nEff <- ratio/(1+ratio)^2*nEvents
  meanStat <- zStat/sqrt(nEff)

  # TODO(Alexander): In principle I could replace "estimate"=exp(meanStat)
  result <- list("statistic"=zStat, "n"=nEvents, "estimate"=meanStat, "sValue"=NULL, "confSeq"=NULL,
                 "estimate"=NULL, "h0"=h0, "testType"="logrank", "dataName"=dataName)
  class(result) <- "safeTest"

  names(result[["estimate"]]) <-"hazard ratio"

  zStat <- (zStat - sqrt(nEff)*log(h0))

  sValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=nEff,
                          "n2"=NULL, "alternative"=alternative, "paired"=FALSE, "sigma"=1)

  tempConfSeq <- computeZConfidenceSequence("nEff"=nEff, "meanStat"=meanStat,
                                            "phiS"=abs(designObj[["parameter"]]), "sigma"=1,
                                            "alpha"=alpha, "alternative"=alternative)

  result[["confSeq"]] <- exp(tempConfSeq)

  result[["sValue"]] <- sValue
  result[["designObj"]] <- designObj
  result[["logrankObj"]] <- logrankObj

  names(result[["statistic"]]) <- "z"
  names(result[["n"]]) <- "nEvents"
  names(result[["h0"]]) <- "theta"

  return(result)
}

#' Designs a Safe Logrank Test
#'
#' A designed experiment requires (1) an anticipated number of events nEvents, or even better nPlan, the number of
#' recruited participants of the study, and (2) the parameter of the safe test, i.e., log(thetaS). If nEvents is
#' provided, then only the safe test defining parameter log(thetaS) needs to determined. That resulting log(thetaS)
#' leads to an (approximately) most powerful safe test. Typically, nEvents is unknown and the user has to specify
#' (i) a tolerable type II error beta, and (b) a clinically relevant minimal hazard ratio hrMin. The procedure finds
#' the smallest nEvents for which hrMin is found with power of at least 1 - beta. The computations exploit the
#' asymptotic normal chracterisation of the sampling distribution of the logrank test derived by Schoenfeld (1981).
#'
#' @inheritParams designSafeZ
#' @param nEvents numeric > 0, targetted number of events.
#' @param hrMin numeric that defines the minimal relevant hazard ratio, the smallest hazard ratio that we want to
#' detect.
#' @param zApprox logical, default TRUE to use the asymptotic normality results.
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If the design is
#' based on nPlan, then ratio equals \code{nPlan[2]/nPlan[1]}.
#' @param alternative The null hypothesis of equality of the survival distribution of y in the groups defined by x is
#' tested. When alternative is "two.sided" the null hypothesis theta = h0, where theta = lambda2/lambda1 is the
#' hazard ratio, and compared to theta != h0. If alternative = "less", the null hypothesis is compared to
#' theta <= h0. For h0=1 the alternative specifies that survival in population 1 is higher than in population 2.
#' When alternative = "greater", the null hypothesis is compared to the alternative theta >= h0. For h0=1 the
#' alternative specifies that survival in population 2 is higher than in population 1.
#'
#' @return Returns a safeDesign object that includes:
#'
#' \describe{
#'   \item{nEvents}{the anticipated number of events, either (1) specified by the user, or
#'   (2) computed based on beta and thetaMin.}
#'   \item{parameter}{the parameter that defines the safe test. Here log(thetaS).}
#'   \item{esMin}{the minimally clinically relevant hazard ratio specified by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{"logrank".}
#'   \item{ratio}{default is 1. It defines the ratio between the planned randomisation of
#'   condition 2 over condition 1.}
#'   \item{pilot}{\code{FALSE} to indicate that the design is not a pilot study.}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#' @export
#'
#' @examples
#' designSafeLogrank(nEvents=89)
#' designSafeLogrank(hrMin=0.4, beta=0.05)
#' designSafeLogrank(hrMin=0.4)
designSafeLogrank <- function(hrMin=NULL, beta=NULL, nEvents=NULL,
                              alternative=c("two.sided", "greater", "less"),
                              alpha=0.05, ratio=1, zApprox=TRUE, tol=1e-5, ...) {
  stopifnot(0 < alpha, alpha < 1)

  alternative <- match.arg(alternative)

  if (zApprox) {
    if (!is.null(hrMin)) {
      logHazardRatio <- if (alternative=="two.sided") abs(log(hrMin)) else log(hrMin)
      meanDiffMin <- abs(logHazardRatio)*sqrt(ratio)/(1+ratio)
    } else {
      logHazardRatio <- NULL
      meanDiffMin <- NULL
    }

    # Note(Alexander): I scaled meanDiffMin so I can get nPlan correct. I'll scale back below
    #
    safeZObj <- designSafeZ("meanDiffMin"=meanDiffMin , "beta"=beta, "alpha"=alpha,
                            "nPlan"=nEvents, "alternative"=alternative,
                            "sigma"=1, "testType"="oneSample")
    nEvents <- safeZObj[["nPlan"]]
    safeZObj[["nPlan"]] <- NULL
    safeZObj[["nEvents"]] <- nEvents

    if (!is.null(nEvents))
      names(safeZObj[["nEvents"]]) <- "nEvents"

    safeZObj[["parameter"]] <- safeZObj[["parameter"]]*(1+ratio)/sqrt(ratio)
    names(safeZObj[["parameter"]]) <- "log(thetaS)"

    safeZObj[["esMin"]] <- hrMin

    if (!is.null(safeZObj[["esMin"]]))
      names(safeZObj[["esMin"]]) <- "hazard ratio"

    # if (!is.null(safeZObj[["esMin"]])) {
    #   names(safeZObj[["esMin"]]) <- switch(alternative,
    #                                        "two.sided"="log hazard difference at least abs(log(theta))",
    #                                        "greater"="log hazard ratio at least",
    #                                        "less"="log hazard ratio less than")
    # }

    safeZObj[["testType"]] <- "logrank"
    safeZObj[["paired"]] <- NULL
    safeZObj[["call"]] <- sys.call()
  }
  result <- safeZObj

  return(result)
}
