#' Safe Logrank Test
#'
#' A safe test to test whether there is a difference between two survival curves. This function
#' builds on the Mantel-Cox version of the logrank test computed with \code{\link[survival]{survdiff}}
#' and adds a sign to the statistic based on the output of \code{\link[coin]{logrank_test}}.
#'
#' @param formula a formula expression as for other survival models, of the form Surv(time, status) ~ groupingVariable,
#' see \code{\link[survival]{Surv}} for more details.
#' @param designObj a safe logrank design obtained from \code{\link{designSafeLogrank}}.
#' @param h0 a number indicating the hypothesised true value of the hazard ratio under the null. Default set to 1
#' @param data an optional data frame in which to interpret the variables occurring in survTime and group
#' @param survTime an optional survival time object of class "Surv" created with \code{\link[survival]{Surv}}, or
#' a name of a column in the data set of class "Surv". Does not need specifying if a formula is provided, therefore
#' set to \code{NULL} by default.
#' @param group an optional factor, a grouping variable. Currently, only two levels allowed. Does not need specifying
#' if a formula is provided, therefore set to \code{NULL} by default.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the number of
#' samples is exactly as planned.
#' @param alpha numeric > 0 only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the alpha of
#' the design object is used instead in constructing the decision rule S > 1/alpha.
#' @param alternative a character string only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the
#' alternative specified by the design object is used instead.
#' @param exact a logical indicating whether the exact safe logrank test needs to be performed based on
#' the hypergeometric likelihood
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the z-statistic.}
#'   \item{nEvents}{The number of observed events.}
#'   \item{eValue}{the s-value of the safe test.}
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
#' # Example taken from survival::survdiff
#' designObj <- designSafeLogrank(hrMin=1/2)
#'
#' ovData <- survival::ovarian
#' ovData$survTime <- survival::Surv(ovData$futime, ovData$fustat)
#'
#' safeLogrankTest(formula=survTime~ rx, data=ovData, designObj=designObj)
#'
#' safeLogrankTest(survTime=survTime, group=rx, data=ovData, designObj=designObj)
#'
#' # Examples taken from coin::logrank_test
#' ## Example data (Callaert, 2003, Tab. 1)
#' #'
#' callaert <- data.frame(
#'   time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
#'   group = factor(rep(0:1, c(7, 8)))
#' )
#'
#' designObj <- designSafeLogrank(hrMin=1/2, beta=0.2)
#'
#' safeLogrankTest(survival::Surv(callaert$time)~callaert$group,
#'                 designObj = designObj)
#'
#' safeLogrankTest(survTime=survival::Surv(callaert$time),
#'                 group=callaert$group, designObj = designObj)
safeLogrankTest <- function(formula, designObj=NULL, h0=1, data=NULL, survTime=NULL,
                            group=NULL, pilot=FALSE, alpha=NULL, alternative=NULL, exact=FALSE, ...) {

  # if (missing(formula) && is.null(survTime) && is.null(group))
  #   stop("Please specify a formula. Or provide survTime and group.")

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

  argumentNames <- getArgs()

  if (!missing(formula)) {
    formulaTerms <- terms(formula)

    if (length(attributes(formulaTerms)[["term.labels"]]) > 1)
      stop("Safe log rank test with covariates not yet supported")

    theData <- try(model.frame(formula, data=data, ...))

    if (isTryError(theData))
      stop("Formula could not be converted into a model.frame.")

    if (!is.null(survTime))
      warning("Both a formula and survTime specified. survTime is overwritten")

    if (!is.null(group))
      warning("Both a formula and group specified. group is overwritten")

    survTime <- theData[, 1]
    group <- theData[, 2]

    if (is.null(survTime))
      stop("Could not extract the survival times from the formula")

    if (is.null(group))
      stop("Could not extract the grouping variable from the formula")

    yLabel <- names(theData)[1]
    groupLabel <- names(theData)[2]
  } else {
    if (!missing(data)) {
      survTime <- data[[argumentNames[["survTime"]]]]
      group <- data[[argumentNames[["group"]]]]
    }

    if (is.null(survTime))
      stop("Could not extract the provided survTime variable.")

    if (is.null(group))
      stop("Could not extract the provided group variable.")

    yLabel <- extractNameFromArgs(argumentNames, "survTime")
    groupLabel <- extractNameFromArgs(argumentNames, "group")
  }

  if (is.null(survTime))
    stop("Can't extract survTime from the given input (formula or survTime)")

  if (is.null(group))
    stop("Can't extract group from the given input (formula or survTime)")

  if (!is.factor(group))
    group <- as.factor(group)

  if (!inherits(survTime, "Surv"))
    stop("Provided variable 'survTime' is not of class 'Surv'.",
         "Please preprocess this variable with the Surv function from the survival package.")

  groupLevels <- levels(group)

  if (length(groupLevels) > 2)
    stop("K-sample log rank test not yet implemented")

  dataSetName <- argumentNames[["data"]]

  dataName <- if (is.null(dataSetName)) "" else paste0(dataSetName, ": ")

  dataName <- paste0(dataName, yLabel, " by ", groupLabel, " (",
                     paste(groupLevels, collapse=", "), ")")

  survDiffObj <- survival::survdiff(survTime ~ group)
  nEvents <- sum(survDiffObj[["obs"]])

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

  if (exact) {
    theta0 <- h0

    if (!is.null(designObj[["esMin"]]))
      theta1 <- unname(designObj[["esMin"]])
    else
      theta1 <- unname(exp(designObj[["parameter"]])) # THIS IS NOT CORRECT Here need to account for ratio

    survTimeMatrixTemp <- as.matrix(survTime)
    eventIndex <- which(survTimeMatrixTemp[, 2]==1)

    survTimeMatrix <- survTimeMatrixTemp[eventIndex, ]
    survTimeDf <- as.data.frame(survTimeMatrix)
    survTimeDf[["group"]] <- group[eventIndex]

    nEvents <- dim(survTimeMatrix)[1]
    timeOfEvents <- unique(survTimeDf[["time"]])

    groupLabel0 <- levels(group)[1]
    groupLabel1 <- levels(group)[2]

    #
    if (nEvents != 0) {
      nTotal <- n0 <- n1 <- logEValue <- vector("numeric", length(timeOfEvents))
      nTotal[1] <- length(survTime)
      n0[1] <- sum(group==groupLabel0)
      n1[1] <- sum(group==groupLabel1)
    } else {
      logEValue <- 0
      nTotal <- length(survTime)
      n0 <- sum(group==groupLabel0)
      n1 <- sum(group==groupLabel1)
    }

    # Should go over event times insteas
    for (i in seq_along(timeOfEvents)) {
      currentTime <- timeOfEvents[i]
      currentDf <- survTimeDf[survTimeDf[["time"]]==currentTime, ]

      casesInGroup0 <- sum(currentDf[["group"]]==groupLabel0)
      casesInGroup1 <- sum(currentDf[["group"]]==groupLabel1)

      logP0 <- log(BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                           n=casesInGroup0+casesInGroup1, odds=theta0))

      if (designObj[["alternative"]]=="two.sided") {
        logP1 <- log(1/2*BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                                 n=casesInGroup0+casesInGroup1, odds=theta1) +
                       1/2*BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                                   n=casesInGroup0+casesInGroup1, odds=1/theta1))
      } else if (designObj[["alternative"]]=="less") {
        logP1 <- log(BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                             n=casesInGroup0+casesInGroup1, odds=theta1))
      } else if (designObj[["alternative"]]=="greater") {
        logP1 <- log(BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                             n=casesInGroup0+casesInGroup1, odds=1/theta1))
      }

      logEValue[i] <- logP1-logP0

      n0[i+1] <- n0[i]-casesInGroup0
      n1[i+1] <- n1[i]-casesInGroup1
    }

    eValue <- exp(sum(logEValue))
    names(eValue) <- "e-value"

    result <- list("statistic"=eValue, "n"=nEvents, "estimate"=NULL, "eValue"=eValue,
                   "confSeq"=NULL, "estimate"=NULL, "h0"=h0, "testType"="logrank",
                   "dataName"=dataName, "exact"=TRUE)
    class(result) <- "safeTest"
    result[["designObj"]] <- designObj
  } else {
    nEff <- ratio/(1+ratio)^2*nEvents

    coinObj <- coin::logrank_test(survTime ~ group, alternative="two.sided")
    signZ <- sign(unname(coinObj@statistic@standardizedlinearstatistic))

    zStat <- signZ*sqrt(survDiffObj[["chisq"]])

    meanStat <- zStat/sqrt(nEff)

    result <- list("statistic"=zStat, "n"=nEvents, "estimate"=exp(meanStat), "eValue"=NULL, "confSeq"=NULL,
                   "estimate"=NULL, "h0"=h0, "testType"="logrank", "dataName"=dataName)
    class(result) <- "safeTest"

    names(result[["estimate"]]) <-"hazard ratio"

    zStat <- (zStat - sqrt(nEff)*log(h0))

    eValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=nEff,
                            "n2"=NULL, "alternative"=alternative, "paired"=FALSE, "sigma"=1)

    tempConfSeq <- computeZConfidenceSequence("nEff"=nEff, "meanStat"=meanStat,
                                              "phiS"=abs(designObj[["parameter"]]), "sigma"=1,
                                              "alpha"=alpha, "alternative"="two.sided")

    result[["confSeq"]] <- exp(tempConfSeq)

    result[["eValue"]] <- eValue
    result[["designObj"]] <- designObj
    result[["survDiffObj"]] <- survDiffObj

    names(result[["statistic"]]) <- "z"
  }


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
