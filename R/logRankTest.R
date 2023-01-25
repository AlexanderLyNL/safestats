#' Safe Logrank Test
#'
#' A safe test to test whether there is a difference between two survival curves. This function
#' builds on the Mantel-Cox version of the logrank test.
#'
#' @inheritParams computeLogrankZ
#' @importFrom survival Surv
#'
#' @param formula a formula expression as for other survival models, of the form Surv(time, status) ~ groupingVariable,
#' see \code{\link[survival]{Surv}} for more details.
#' @param designObj a safe logrank design obtained from \code{\link{designSafeLogrank}}.
#' @param data an optional data frame in which to interpret the variables occurring in survTime and group.
#' @param survTime an optional survival time object of class 'Surv' created with \code{\link[survival]{Surv}}, or
#' a name of a column in the data set of class 'Surv'. Does not need specifying if a formula is provided, therefore
#' set to \code{NULL} by default.
#' @param group an optional factor, a grouping variable. Currently, only two levels allowed. Does not need specifying
#' if a formula is provided, therefore set to \code{NULL} by default.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the number of
#' samples is exactly as planned. The default null h0=1 is used, alpha=0.05, and alternative="twoSided" is used.
#' To change these default values, please use \code{\link{designSafeLogrank}}.
#' @param ciValue numeric, represents the ciValue-level of the confidence sequence. Default ciValue=NULL, and
#' ciValue = 1 - alpha, where alpha is taken from the design object.
#' @param exact a logical indicating whether the exact safe logrank test needs to be performed based on
#' the hypergeometric likelihood. Default is \code{TRUE}, if \code{FALSE} then the safe z-test (for Gaussian data)
#' applied to the logrank z-statistic is used instead.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class 'safeTest'. An object of class 'safeTest' is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the summary, i.e., z-statistic or the e-value.}
#'   \item{nEvents}{The number of observed events.}
#'   \item{eValue}{the e-value of the safe test.}
#'   \item{confSeq}{An anytime-valid confidence sequence.}
#'   \item{estimate}{To be implemented: An estimate of the hazard ratio.}
#'   \item{testType}{"logrank".}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" obtained from \code{\link{designSafeLogrank}}.}
#'   \item{sumStats}{a list containing.the time of events, the progression of the risk sets and events.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' # Example taken from survival::survdiff
#'
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
#' designObj <- designSafeLogrank(hrMin=1/2)
#'
#' safeLogrankTest(survival::Surv(callaert$time)~callaert$group,
#'                 designObj = designObj)
#'
#' safeLogrankTest(survTime=survival::Surv(callaert$time),
#'                 group=callaert$group, designObj = designObj)
#'
#' result <- safeLogrankTest(survTime=survival::Surv(callaert$time),
#'                 group=callaert$group, designObj = designObj)
#'
#' result
#'
#' ##  Sequentially
#' # Greater
#' eValueGreater <- exp(cumsum(result$sumStats$logEValueGreater))
#' # Less
#' eValueLess <- exp(cumsum(result$sumStats$logEValueLess))
#'
#' # twoSided
#' eValueTwoSided <- 1/2*eValueGreater+1/2*eValueLess
#'
#' eValueTwoSided
#' result$eValue
#'
#' ###### Example switching between safe exact and safe Gaussian logrank test
#'
#' designObj <- designSafeLogrank(0.8, alternative="less")
#'
#' dat <- safestats::generateSurvData(300, 300, 2, 0.0065, 0.0065*0.8, seed=1)
#' survTime <- survival::Surv(dat$time, dat$status)
#'
#' resultE <- safeLogrankTest(survTime ~ dat$group,
#'                            designObj = designObj)
#'
#' resultG <- safeLogrankTest(survTime ~ dat$group,
#'                            designObj = designObj, exact=FALSE)
#'
#' resultE
#' resultG
#'
#' ###### Example switching between safe exact and safe Gaussian logrank test other side
#'
#' designObj <- designSafeLogrank(1/0.8, alternative="greater")
#'
#' resultE <- safeLogrankTest(survTime ~ dat$group,
#'                            designObj = designObj)
#'
#' resultG <- safeLogrankTest(survTime ~ dat$group,
#'                            designObj = designObj, exact=FALSE)
#'
#' if (log(resultE$eValue) >= 0 && log(resultG$eValue) >= 0 )
#'   stop("one-sided wrong")
#'
safeLogrankTest <- function(formula, designObj=NULL, ciValue=NULL, data=NULL, survTime=NULL,
                            group=NULL, pilot=FALSE, exact=TRUE, computeZ=TRUE, ...) {

  # Check inputs  ----
  #
  if (isFALSE(pilot) && is.null(designObj))
    stop("Please provide a safe logrank design object, or run the function with pilot=TRUE. ",
         "A design object can be obtained by running designSafeLogrank().")

  if (!is.null(designObj)) {
    if (names(designObj[["parameter"]]) != "log(thetaS)" && names(designObj[["parameter"]]) != "thetaS")
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

  # Check data survTime----
  #
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

  survType <- attr(survTime, "type")

  # Check designObj ----
  #
  if (isTRUE(pilot)) {
    alpha <- 0.05
    alternative <- "twoSided"
    h0 <- 1

    survTimeMatrix <- as.matrix(survTime)
    nEvents <- sum(survTimeMatrix[, "status"]==1)

    if (is.null(designObj)) {
      designObj <- designSafeLogrank("hrMin"=NULL, "beta"=NULL, "nEvents"=nEvents, "alpha"=alpha,
                                     "alternative"=alternative, "h0"=h0, "exact"=FALSE)
      designObj[["pilot"]] <- TRUE
    } else {
      warning("The pilot flag is ignored, since a designObj is given",
              "The analysis will be run based on the designObj.")
    }
  }

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  ratio <- designObj[["ratio"]]
  h0 <- designObj[["h0"]]

  thetaS <- designObj[["parameter"]]

  # Note(Alexander): Sign is needed for the safe Gaussian logrank test
  #
  phiS <- log(thetaS)

  # Note(Alexander): Remove sign for safe exact logrank test
  #
  if (thetaS > 1)
    thetaS <- 1/thetaS

  # Compute stats ------
  #
  sumStats <- computeLogrankZ("survObj"=survTime, "group"=group,
                              "computeZ"=computeZ, "computeExactE"=exact,
                              "theta0"=h0, "thetaS"=thetaS)
  nEvents <- sumStats[["nEvents"]]

  # Compute e-value ------
  #
  if (exact) {
    if (designObj[["alternative"]]=="twoSided" || designObj[["alternative"]]=="less")
      eValueLess <- exp(sum(sumStats[["logEValueLess"]]))

    if (designObj[["alternative"]]=="twoSided" || designObj[["alternative"]]=="greater")
      eValueGreater <- exp(sum(sumStats[["logEValueGreater"]]))

    eValue <- switch(designObj[["alternative"]],
                     "greater"=eValueGreater,
                     "less"=eValueLess,
                     "twoSided"=1/2*eValueGreater+1/2*eValueLess)

    names(eValue) <- "e"

    result <- list("n"=nEvents, "estimate"=NULL, "eValue"=eValue,
                   "confSeq"=NULL, "estimate"=NULL, "testType"="eLogrank",
                   "dataName"=dataName, "exact"=TRUE)
    class(result) <- "safeTest"
    result[["designObj"]] <- designObj
  } else {
    nEff <- ratio/(1+ratio)^2*nEvents

    zStat <- sumStats[["z"]]
    sumOMinE <- sumStats[["sumOMinE"]]
    sumVarOMinE <- sumStats[["sumVarOMinE"]]
    
    result <- list("statistic"=zStat, "n"=nEvents, "estimate"=exp(sumOMinE/sumVarOMinE), "eValue"=NULL,
                   "confSeq"=NULL, "testType"="gLogrank", "dataName"=dataName)
    class(result) <- "safeTest"

    names(result[["estimate"]]) <-"hazard ratio"

    # Note(Alexander): This is the same as
    #     zStat <- sqrt(nEff)*(meanObs - meanSlog(h0))
    #
    # but to avoid rounding erros zStat is used instead
    zStat <- zStat - sqrt(nEff)*(log(h0))

    eValue <- safeZTestStat("z"=zStat, "phiS"=phiS, "n1"=nEff,
                            "n2"=NULL, "alternative"=alternative, "paired"=FALSE, "sigma"=1)

    if (is.null(ciValue))
      ciValue <- 1 - designObj[["alpha"]]

    tempConfSeq <- computeConfidenceIntervalZ("nEff"=nEff, "meanObs"=meanObs,
                                              "phiS"=phiS, "sigma"=1,
                                              "ciValue"=ciValue, "alternative"="twoSided")

    result[["ciValue"]] <- ciValue

    result[["confSeq"]] <- exp(tempConfSeq)

    result[["eValue"]] <- eValue

    names(result[["statistic"]]) <- "z"
  }

  sumStats[["ratio"]] <- designObj[["ratio"]]

  names(result[["n"]]) <- "nEvents"
  result[["designObj"]] <- designObj
  result[["sumStats"]] <- sumStats

  return(result)
}


#' @describeIn safeLogrankTest Safe Logrank Test based on Summary Statistic Z
#' All provided data (i.e., z-scores) are assumed to be centred on a hazard ratio = 1, thus, log(hr) = 0 ,
#' and the proper (e.g., hypergeometric) scaling is applied to the data, so sigma = 1. The null hypothesis
#' in the design object pertains to the population and is allowed to differ from log(theta) = 0.
#'
#' @param z numeric representing the observed logrank z statistic.
#' @param nEvents numeric > 0, observed number of events.
#' @param dataNull numeric > 0, the null hypothesis corresponding to the z statistics.
#' By default dataNull = 1 representing equality of the hazard ratio.
#' @param sigma numeric > 0, scaling in the data.
#'
#' @export
safeLogrankTestStat <- function(z, nEvents, designObj, ciValue=NULL,
                                dataNull=1, sigma=1) {

  if (length(z) != length(nEvents))
    stop("The provided number of z-scores and number of events not equal.")

  names(nEvents) <- "nEvents"

  result <- list("statistic"=z, "n"=nEvents, "estimate"=NULL, "eValue"=NULL,
                 "confSeq"=NULL, "testType"="logrank", "dataName"="Logrank z")

  if (is.null(ciValue))
    ciValue <- 1 - designObj[["alpha"]]

  nEff <- designObj[["ratio"]]/(1+designObj[["ratio"]])^2*nEvents

  # Note(Alexander): Assumed the data are centred at log(hr)=0 and standardised,
  # thus, sigma = 1 data scale
  if (length(z)==1) {
    meanObs <- sigma*z/sqrt(nEff) + log(dataNull)
    # TODO(Alexander): For the standard version with h0 = 1 this doesn't matter at all of course
    #                  Do check when dataNull different from population null again
    zStat <- z - sqrt(nEff)/sigma*log(designObj[["h0"]])
  } else {
    meanObs <- sum(sigma*z/sqrt(nEff)+log(dataNull))/sum(nEff)
    zStat <- sqrt(nEff)/sigma*(meanObs - log(designObj[["h0"]]))
  }

  phiS <- log(designObj[["parameter"]])

  eValue <- safeZTestStat("z"=zStat, "phiS"=phiS, "n1"=nEff, "n2"=NULL,
                          "alternative"=designObj[["alternative"]], "paired"=FALSE, "sigma"=1)

  tempConfSeq <- computeConfidenceIntervalZ("nEff"=nEff, "meanObs"=meanObs, "phiS"=phiS,
                                            "sigma"=1, "ciValue"=ciValue, "alternative"="twoSided")

  result[["confSeq"]] <- exp(tempConfSeq)

  result[["eValue"]] <- eValue
  result[["designObj"]] <- designObj

  names(result[["statistic"]]) <- "z"
  class(result) <- "safeTest"
  return(result)
}


#' Designs a Safe Logrank Test Experiment
#'
#' A designed experiment requires (1) an anticipated number of events nEvents, or even better nPlan, the number of
#' participants to be recruited in the study, and (2) the parameter of the safe test, i.e., thetaS. Provided with a
#' clinically relevant minimal hazard ratio hrMin, this function outputs thetaS = hrMin as the safe test defining
#' parameter in accordance to the GROW criterion. If a tolerable type II error beta is provided then nEvents can be
#' sampled. The sampled nEvents is then the smallest nEvents for which hrMin is found with power of at least 1 - beta
#' under optional stopping. If exact equal \code{FALSE}, then the computations exploit the local asymptotic normal
#' approximation to sampling distribution of the logrank test derived by Schoenfeld (1981).
#'
#'
#' @inheritParams designSafeZ
#' @param nEvents numeric > 0, targetted number of events.
#' @param h0 numeric > 0, represents the null hypothesis, default h0=1.
#' @param hrMin numeric that defines the minimal relevant hazard ratio, the smallest hazard ratio that we want to
#' detect.
#' @param exact a logical indicating whether the design should be based on the exact safe logrank test based on the
#' hypergeometric likelihood. Default is \code{TRUE}, if \code{FALSE} then the design is based on a  safe z-test.
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 (Treatment) over condition 1 (Placebo),
#' thus, m1/m0. Note that m1 and m0 are not used to specify ratio. Ratio is only used when \code{zApprox=TRUE}, which
#' ignores m1 and m0.
#' @param parameter numeric > 0 representing the test defining thetaS. Default is NULL, then GROW the choice is used,
#' that is, parameter equals the data generating hazardRatio.
#' @param alternative a character string specifying the alternative hypothesis, which must be one of
#' "twoSided" (default),"greater" or "less". The alternative is pitted against the null hypothesis of equality
#' of the survival distributions. More specifically, let lambda1 be the hazard rate of group 1 (i.e., placebo), and
#' lambda2 the hazard ratio of group 2 (i.e., treatment), then the null hypothesis states that the hazard ratio
#' theta = lambda2/lambda1 = 1. If alternative = "less", the null hypothesis is compared to theta < 1, thus,
#' lambda2 < lambda1, that is, the hazard of group 2 (i.e., treatment) is less than that of group 1 (i.e., placebo),
#' hence, the treatment is beneficial. If alternative = "greater", then the null hypothesis is compared to theta > 1,
#' thus, lambda2 > lambda1, hence, harm.
#' @param m0 Number of subjects in the control group 0/1 at the beginning of the trial, i.e., nPlan[1].
#' @param m1 Number of subjects in the treatment group 1/2 at the beginning of the trial, i.e., nPlan[2].
#' @param parameter Numeric > 0, represents the safe tests defining thetaS. Default NULL so it's decided by the
#' algorithm, typically, this equals hrMin, which corresponds to the GROW choice.
#' @param nSim integer > 0, the number of simulations needed to compute power or the number of events for the exact
#' safe logrank test under continuous monitoring
#' @param groupSizePerTimeFunction A function without parameters and integer output. This function provides the number
#' of events at each time step. For instance, if \code{rpois(1, 7)} leads to a random number of events at each time
#' step.
#' @param nBoot integer > 0 representing the number of bootstrap samples to assess the accuracy of the approximation of
#' power or nEvents for the exact safe logrank test under continuous monitoring
#' @param pb logical, if \code{TRUE}, then show progress bar.
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
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
#'   \item{testType}{"logrank".}
#'   \item{ratio}{default is 1. It defines the ratio between the planned randomisation of
#'   condition 2 over condition 1.}
#'   \item{pilot}{\code{FALSE} to indicate that the design is not a pilot study.}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#' @export
#'
#' @references Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests
#' for comparing survival distributions. Biometrika, 68(1), 316-319.
#'
#' @examples
#' designSafeLogrank(hrMin=0.7)
#' designSafeLogrank(hrMin=0.7, zApprox=TRUE)
#' designSafeLogrank(hrMin=0.7, beta=0.3, nSim=10)
#' designSafeLogrank(hrMin=0.7, nEvents=190, nSim=10)
designSafeLogrank <- function(hrMin=NULL, beta=NULL, nEvents=NULL, h0=1,
                              alternative=c("twoSided", "greater", "less"),
                              alpha=0.05, ratio=1, exact=TRUE, tol=1e-5,
                              m0=50000L, m1=50000L, nSim=1e3L, nBoot=1e4L,
                              parameter=NULL, groupSizePerTimeFunction=returnOne,
                              pb=TRUE, ...) {
  stopifnot(0 < alpha, alpha < 1)

  result <- list()

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  if (!is.null(hrMin))
    hrMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=hrMin, "alternative"=alternative, "esMinName"="hrMin")

  if (!is.null(parameter))
    parameter <- checkAndReturnsEsMinParameterSide("paramToCheck"=parameter, "alternative"=alternative, "esMinName"="thetaS")

  thetaS <- if (is.null(parameter)) hrMin else parameter
  note <- NULL

  nEventsBatch <- nEventsTwoSe <- NULL
  nMean <- nMeanTwoSe <- NULL

  logImpliedTarget <- logImpliedTargetTwoSe <- NULL
  betaTwoSe <- NULL

  bootObjNEvents <- bootObjN1Mean <- bootObjBeta <- bootObjLogImpliedTarget <- NULL

  if (!exact) {
    if (!is.null(hrMin)) {
      logHazardRatio <- if (alternative=="twoSided") abs(log(hrMin)) else log(hrMin)
      meanDiffMin <- logHazardRatio*sqrt(ratio)/(1+ratio)
    } else {
      logHazardRatio <- NULL
      meanDiffMin <- NULL
    }

    # Note(Alexander): I scaled meanDiffMin so I can get nPlan correct.
    # I'll scale back below
    #
    safeZObj <- designSafeZ("meanDiffMin"=meanDiffMin , "beta"=beta,
                            "alpha"=alpha, "nPlan"=nEvents,
                            "alternative"=alternative,
                            "sigma"=1, "testType"="oneSample")

    nEvents <- safeZObj[["nPlan"]]
    safeZObj[["nPlan"]] <- NULL
    safeZObj[["nEvents"]] <- nEvents

    nEventsBatch <- safeZObj[["nPlanBatch"]]
    safeZObj[["nPlanBatch"]] <- NULL
    safeZObj[["nEventsBatch"]] <- nEventsBatch

    nEventsTwoSe <- safeZObj[["nPlanTwoSe"]]
    safeZObj[["nPlanTwoSe"]] <- NULL
    safeZObj[["nEventsTwoSe"]] <- nEventsTwoSe

    if (!is.null(nEventsBatch)) {
      note <- paste0("If it is only possible to look at the data once, ",
                     "then nEvents = ", nEventsBatch, ".")
    }

    safeZObj[["note"]] <- note

    if (!is.null(nEvents))
      names(safeZObj[["nEvents"]]) <- "nEvents"

    safeZObj[["parameter"]] <- if (!is.null(parameter)) {
      parameter
    } else if (!is.null(hrMin)) {
      thetaS
    } else {
      exp(safeZObj[["parameter"]])
    }

    names(safeZObj[["parameter"]]) <- "thetaS"

    if (!is.null(hrMin))
      names(hrMin) <- "hazard ratio"

    safeZObj[["esMin"]] <- hrMin

    safeZObj[["testType"]] <- "gLogrank"
    safeZObj[["paired"]] <- NULL
    safeZObj[["call"]] <- sys.call()

    if (!is.null(h0))
      names(h0) <- "theta"

    safeZObj[["h0"]] <- h0
    result <- safeZObj
    result[["ratio"]] <- ratio
    result[["exact"]] <- exact

    return(result)
  } else {
    designScenario <- NULL

    ratio <- m1/m0

    if (!is.null(hrMin) && !is.null(beta) && is.null(nEvents)) {
      designScenario <- "1a"

      tempResult <- computeLogrankNEvents("hrMin"=hrMin, "beta"=beta, "m0"=m0, "m1"=m1, "alpha"=alpha,
                                          "alternative"=alternative, "nSim"=nSim, "nBoot"=nBoot,
                                          "groupSizePerTimeFunction"=groupSizePerTimeFunction,
                                          "parameter"=parameter, "pb"=pb)
      nEvents <- tempResult[["nEvents"]]

      bootObjNEvents <- tempResult[["bootObjNEvents"]]
      nEventsTwoSe <- 2*bootObjNEvents[["bootSe"]]

      nMean <- tempResult[["n1Mean"]]
      names(nMean) <- "nMean"
      bootObjN1Mean <- tempResult[["bootObjN1Mean"]]
      nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]

      nEventsBatch <- bootObjNEvents[["nEventsBatch"]]

      if (!is.null(nEventsBatch) && is.finite(nEventsBatch)) {
        note <- paste0("If it is only possible to look at the data once, ",
                       "then nEvents = ", nEventsBatch, ".")
      }

    } else if (!is.null(hrMin) && is.null(beta) && is.null(nEvents)) {
      designScenario <- "1b"
    } else if (is.null(hrMin) && is.null(beta) && !is.null(nEvents)) {
      designScenario <- "1c"

      # TODO(Alexander):
      # - For the case with only nEvents,
      # - How to do GROW? Or forget about it.
      #
      #       Do we have a PILOT version for this?

      warning("Designs without minimal clinically relevant hazard ratios not yet implemented")
      # return(designPilotSafeZ("nPlan"=nPlan, "alpha"=alpha, "alternative"=alternative,
      #                         "sigma"=sigma, "kappa"=kappa, "tol"=tol, "paired"=paired))
    } else if (!is.null(hrMin) && is.null(beta) && !is.null(nEvents)) {
      designScenario <- "2"

      tempResult <- computeLogrankBetaFrom("hrMin"=hrMin, "nEvents"=nEvents, "m0"=m0, "m1"=m1, "alpha"=alpha,
                                           "alternative"=alternative, "nSim"=nSim, "nBoot"=nBoot,
                                           "groupSizePerTimeFunction"=groupSizePerTimeFunction,
                                           "parameter"=thetaS, "pb"=pb)

      beta <- tempResult[["beta"]]
      bootObjBeta <- tempResult[["bootObjBeta"]]
      betaTwoSe <- 2*bootObjBeta[["bootSe"]]

      logImpliedTarget <- tempResult[["logImpliedTarget"]]
      bootObjLogImpliedTarget <- tempResult[["bootObjLogImpliedTarget"]]
      logImpliedTargetTwoSe <- 2*bootObjLogImpliedTarget[["bootSe"]]
    } else if (is.null(hrMin) && !is.null(beta) && !is.null(nEvents)) {
      designScenario <- "3"
      designScenario <- NULL

      # TODO(Alexander): Normal approximation. Flag normal approximation
      warning("Designs without minimal clinically relevant hazard ratios not yet implemented")
    }

    if (is.null(designScenario)) {
      stop("Can't design: Please provide this function with either: \n",
           "(1.a) non-null hrMin, non-null beta and NULL nEvents, or \n",
           "(1.b) non-null hrMin, NULL beta, and NULL nEvents, or \n",
           # "(1.c) NULL hrMin, NULL beta, non-null nEvents, or \n",
           "(2) non-null hrMin, NULL beta and non-null nEvents.")
      # "(3) NULL hrMin, non-null beta, and non-null nEvents.")
    }

    if (is.na(hrMin))
      hrMin <- NULL

    if (!is.null(nEvents))
      names(nEvents) <- "nEvents"

    if (!is.null(nEventsBatch))
      names(nEventsBatch) <- "nEventsBatch"

    if (!is.null(hrMin))
      names(hrMin) <- "hazard ratio"

    if (!is.null(thetaS))
      names(thetaS) <- "thetaS"

    if (!is.null(h0))
      names(h0) <- "theta"

    result <- list("parameter"=thetaS, "esMin"=hrMin, "alpha"=alpha, "alternative"=alternative,
                   "h0"=h0, "testType"="eLogrank", "exact"=exact,
                   "ratio"=m1/m0, "pilot"=FALSE,
                   "nPlan"=nEvents, "nPlanTwoSe"=nEventsTwoSe, "nPlanBatch"=nEventsBatch,
                   "nMean"=nMean, "nMeanTwoSe"=nMeanTwoSe,
                   "beta"=beta, "betaTwoSe"=betaTwoSe,
                   "logImpliedTarget"=logImpliedTarget, "logImpliedTargetTwoSe"=logImpliedTargetTwoSe,
                   "bootObjNPlan"=bootObjNEvents, "bootObjBeta"=bootObjBeta,
                   "bootObjLogImpliedTarget"=bootObjLogImpliedTarget, "bootObjN1Mean"=bootObjN1Mean,
                   "call"=sys.call(), "timeStamp"=Sys.time(), "note"=note)

    class(result) <- "safeDesign"

    return(result)
  }
}



#' Helper function computes single component of the logrank statistic
#'
#' @param obs0 integer, number of observations in the control group
#' @param obs1 integer, number of observations in the treatment group
#' @param y0 integer, total number of participants in the control group
#' @param y1 integer, total number of participants in the treatment group
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list containing at least the following components:
#' \describe{
#'   \item{oMinE}{observed minus expected.}
#'   \item{v}{hypergeometric variance.}
#' }
#' @export
#'
#' @examples
#' y0Vector <- c(6, 4, 4, 1, 0)
#' y1Vector <- c(6, 6, 5, 2, 2)
#' obs0Vector <- c(1, 0, 2, 1, 0)
#' obs1Vector <- c(0, 1, 1, 0, 1)
#'
#' varVector <- oMinEVector <-y0Vector
#'
#' for (i in seq_along(y0Vector)) {
#'   tempResult <- logrankSingleZ(obs0=obs0Vector[i], obs1=obs1Vector[i],
#'                               y0=y0Vector[i], y1=y1Vector[i])
#'   oMinEVector[i] <- tempResult[["oMinE"]]
#'   varVector[i] <- tempResult[["v"]]
#' }
#'
#' sum(oMinEVector)/sqrt(sum(varVector))
#'
logrankSingleZ <- function(obs0, obs1, y0, y1, ...) {
  dTotal <- obs0 + obs1
  yTotal <- y0 + y1

  o1 <- obs1
  e1 <- dTotal*y1/yTotal

  if (yTotal==1) {
    variance <- 0
  } else {
    logVar <- log(y0)+log(y1)+log(dTotal)+log(yTotal-dTotal) -
      (2*log(yTotal) + log(yTotal-1))
    variance <- exp(logVar)
  }

  result <- list("oMinE"=o1-e1, "v"=variance)
  return(result)
}

#' Helper function computes single component of the exact logrank e-value
#'
#' @param obs0 integer, number of observations in the control group.
#' @param obs1 integer, number of observations in the treatment group.
#' @param y0 integer, total number of participants in the control group.
#' @param y1 integer, total number of participants in the treatment group.
#' @param y1 integer, total number of participants in the treatment group.
#' @param thetaS numeric > 0 represents the safe test defining (GROW) alternative
#' hypothesis obtained from \code{designSafeLogrank()}.
#' @param theta0 numeric > 0 represents the null hypothesis. Default theta0=1.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list containing at least the following components:
#' \describe{
#'   \item{logP0}{Log likelihood of Fisher's hypergeometric at the null}
#'   \item{logEValueLess}{Log likelihood of Fisher's hypergeometric at the alternative}
#'   \item{logEValueGreater}{Log likelihood of Fisher's hypergeometric at 1/alternative}
#' }
#' @export
#'
#' @examples
#' #'
#' y0Vector <- c(5, 4, 3, 3, 2, 1)
#' y1Vector <- c(5, 5, 4, 2, 2, 0)
#' obs0Vector <- c(1, 1, 0, 1, 0, 1)
#' obs1Vector <- c(0, 0, 1, 0, 1, 0)
#'
#' logEValueGreater <- logEValueLess <- vector("numeric", length(y0Vector))
#'
#' for (i in seq_along(y0Vector)) {
#'   tempResult <- logrankSingleEExact(obs0=obs0Vector[i], obs1=obs1Vector[i],
#'                                     y0=y0Vector[i], y1=y1Vector[i],
#'                                     thetaS=0.7, theta0=1)
#'   logEValueLess[i] <- tempResult[["logEValueLess"]]
#'   logEValueGreater[i] <- tempResult[["logEValueGreater"]]
#' }
#'
#' eValueLess <- exp(sum(logEValueLess))
#' eValueLess #1.116161
#' eValueGreater <- exp(sum(logEValueGreater))
#' eValueGreater # 0.7665818
#' eValue <- 1/2*eValueLess + 1/2*eValueGreater
#' eValue # 0.9413714
#'
logrankSingleEExact <- function(obs0, obs1, y0, y1, thetaS, theta0=1, ...) {
  logP0 <- log(BiasedUrn::dFNCHypergeo(x=obs1, m1=y1, m2=y0, n=obs0+obs1, odds=theta0))
  logPLess <- log(BiasedUrn::dFNCHypergeo(x=obs1, m1=y1, m2=y0, n=obs0+obs1, odds=thetaS))
  logEValueLess <- logPLess-logP0
  logPGreater <- log(BiasedUrn::dFNCHypergeo(x=obs1, m1=y1, m2=y0, n=obs0+obs1, odds=1/thetaS))
  logEValueGreater <- logPGreater-logP0

  result <- list("logP0"=logP0, "logEValueLess"=logEValueLess, "logEValueGreater"=logEValueGreater)
  return(result)
}

#' Computes the sufficient statistics needed to compute 'logrankSingleZ'
#'
#' @param survDataFrame a 'Surv' object converted to a matrix, then to a data.frame
#' @param y0Index vector of integers corresponding to the control group
#' @param y1Index vector of integers corresponding to the treatment group
#' @param timeNow numeric, current time
#' @param timeBefore numeric, previous time
#' @param survType character, either "right" or "counting" (left truncated, right censored)
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list containing at least the following components:
#' \describe{
#'   \item{obs0}{number of observations in the control group.}
#'   \item{obs1}{number of observations in the treatment group.}
#'   \item{y0}{total number of participants in the control group.}
#'   \item{y1}{total number of participants in the treatment group.}#'
#' }
#' @export
#'
#' @examples
#'
#' data <- generateSurvData(nP = 5,
#'                          nT = 5,
#'                          lambdaP = 0.03943723,
#'                          lambdaT = 0.5*0.03943723,
#'                          endTime = 40,
#'                          seed = 2006)
#'
#' survObj <- survival::Surv(data$time, data$status)
#'
#' survDataFrame <- as.data.frame(as.matrix(survObj))
#' y0Index <- which(data$group=="P")
#' y1Index <- which(data$group=="T")
#'
#' timeNow <- 4
#' timeBefore <- 0
#'
#' computeStatsForLogrank(survDataFrame, y0Index, y1Index, timeNow, timeBefore)
#'
#' timeNow <- 13
#' timeBefore <- 4
#'
#' computeStatsForLogrank(survDataFrame, y0Index, y1Index, timeNow, timeBefore)
computeStatsForLogrank <- function(survDataFrame, y0Index, y1Index, timeNow, timeBefore,
                                   survType="right", ...) {
  timeLabel <- switch(survType,
                      "counting"="stop",
                      "right"="time")

  eventIndex <- which(survDataFrame[[timeLabel]]==timeNow &
                        survDataFrame[["status"]]==1)

  obs0 <- length(intersect(eventIndex, y0Index))
  obs1 <- length(intersect(eventIndex, y1Index))

  currentIndex <- which(survDataFrame[[timeLabel]] >= timeNow &
                          survDataFrame[[timeLabel]] > timeBefore)

  y0 <- length(intersect(currentIndex, y0Index))
  y1 <- length(intersect(currentIndex, y1Index))

  result <- list("obs0"=obs0, "obs1"=obs1, "y0"=y0, "y1"=y1)
  return(result)
}


#' Helper function to computes the logrank statistic for 'Surv' objects of type
#' "right" and "counting" with the hypergeometric variance.
#'
#' This function was created to complement \code{\link[survival]{survdiff}} from the
#' 'survival' package, which is restricted to 'Surv' objects of type "right". Most likely
#' \code{\link[survival]{survdiff}} is much faster
#'
#' @param survObj a Surv object that is either of type
#' @param group a grouping factor with 2 levels
#' @param computeZ logical. If \code{TRUE} computes the logrank z-statistic.
#' Default is \code{TRUE}.
#' @param computeExactE logical. If \code{TRUE} computes one-sided exact logrank e-value.
#' Default is \code{FALSE}.
#' @param theta0 numeric > 0 used only for the e-value, i.e., if computeExactE is \code{TRUE}.
#' Default is 1.
#' @param thetaS numeric > 0 used only for the e-value, i.e., if computeExactE is \code{TRUE}.
#' Default is NULL.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list containing at least the following components:
#' \describe{
#'   \item{nEvents}{the number of events.}
#'   \item{z}{the observed logrank statistic.}
#'   \item{oMinEVector}{vector of observed minus expected.}
#'   \item{varVector}{vector of hypergeometric variances.}
#'   \item{stopTimeVector}{vector at which the events occurred.}
#' }
#'
#' @export
#'
#' @examples
#' data <- generateSurvData(nP = 5,
#'                          nT = 5,
#'                          lambdaP = 0.03943723,
#'                          lambdaT = 0.5*0.03943723,
#'                          endTime = 40,
#'                          seed = 2006)
#'
#' survObj <- survival::Surv(data$time, data$status)
#'
#' survObj <- survival::Surv(data$time, data$status)
#'
#' result <- computeLogrankZ(survObj, data$group)
#' result$z
#' sqrt(survival::survdiff(survObj~data$group)$chisq)
computeLogrankZ <- function(survObj, group, computeZ=TRUE, computeExactE=FALSE,
                            theta0=1, thetaS=NULL, ...) {
  result <- list()

  # Check exact requirements -----
  #
  if (computeExactE)
    if (is.null(thetaS))
      stop("Can't compute exact E-value without a designed alternative.",
           "Please check designSafeLogrank.")

  # Get group label info -----
  #
  groupElements <- unique(group)

  if (length(groupElements) > 2)
    stop("Logrank with more than 2 groups not yet implemented")

  if (length(groupElements) < 2)
    stop("Only data of one of the groups")

  groupElementsOrdered <- groupElements[order(groupElements)]

  groupLabel0 <- groupElementsOrdered[1]
  groupLabel1 <- groupElementsOrdered[2]

  survType <- attr(survObj, "type")

  timeLabel <- switch(survType,
                      "right"="time",
                      "counting"="stop")

  survDataFrame <- as.data.frame(as.matrix(survObj))
  survDataFrame[["group"]] <- group

  stopTimeIndeces <- which(survDataFrame[["status"]]==1)

  nEvents <- length(stopTimeIndeces)

  stopTimeVector <- unique(survDataFrame[[timeLabel]][stopTimeIndeces])
  stopTimeVector <- stopTimeVector[order(stopTimeVector)]

  lengthStopTime <- length(stopTimeVector)

  # Init: Resultvectors ----
  #
  varVector <- oMinEVector <- logEValueGreater <- logEValueLess <- logP0 <- NULL

  if (lengthStopTime > 0) {
    y0Vector <- y1Vector <- obs0Vector <- obs1Vector <- rep(NA, length = lengthStopTime)

    if (computeZ)
      varVector <- oMinEVector <- rep(NA, length = lengthStopTime)

    if (computeExactE)
      logEValueGreater <- logEValueLess <- logP0 <- rep(NA, length = lengthStopTime)

    timeBeforeVector <- c(0, stopTimeVector[1:(lengthStopTime-1)])
  } else {
    warning("No events")
    result <- list(n=0, z=0, oMinEVector=NULL, varVector=NULL, stopTimeVector=NULL,
                   y0Vector=NULL, y1Vector=NULL, obs0Vector=NULL, obs1Vector=NULL)
    return(result)
  }

  # Loop data -----
  #
  if (survType=="counting") {
    for (i in seq_along(stopTimeVector)) {
      timeNow <- stopTimeVector[i]
      timeBefore <- timeBeforeVector[i]

      subSurvDataFrame <- survDataFrame[survDataFrame[["start"]] < timeNow, ]

      y0Index <- which(subSurvDataFrame[["group"]]==groupLabel0)
      y1Index <- which(subSurvDataFrame[["group"]]==groupLabel1)

      tempStats <- computeStatsForLogrank("survDataFrame"=subSurvDataFrame,
                                          "y0Index"=y0Index,
                                          "y1Index"=y1Index,
                                          "timeNow"=timeNow,
                                          "timeBefore"=timeBefore,
                                          "survType"="counting")

      y0Vector[i] <- tempStats[["y0"]]
      y1Vector[i] <- tempStats[["y1"]]
      obs0Vector[i] <- tempStats[["obs0"]]
      obs1Vector[i] <- tempStats[["obs1"]]

      if (computeZ) {
        tempResult <- do.call(logrankSingleZ, tempStats)
        oMinEVector[i] <- tempResult[["oMinE"]]
        varVector[i] <- tempResult[["v"]]
      }

      if (computeExactE) {
        tempResult <- logrankSingleEExact(obs0=tempStats[["obs0"]], obs1=tempStats[["obs1"]],
                                          y0=tempStats[["y0"]], y1=tempStats[["y1"]],
                                          theta0=theta0, thetaS=thetaS)
        logP0[i] <- tempResult[["logP0"]]
        logEValueLess[i] <- tempResult[["logEValueLess"]]
        logEValueGreater[i] <- tempResult[["logEValueGreater"]]
      }
    }
  } else if (survType=="right") {
    y0Index <- which(survDataFrame[["group"]]==groupLabel0)
    y1Index <- which(survDataFrame[["group"]]==groupLabel1)

    for (i in seq_along(stopTimeVector)) {
      timeNow <- stopTimeVector[i]
      timeBefore <- timeBeforeVector[i]

      tempStats <- computeStatsForLogrank("survDataFrame"=survDataFrame,
                                          "y0Index"=y0Index,
                                          "y1Index"=y1Index,
                                          "timeNow"=timeNow,
                                          "timeBefore"=timeBefore,
                                          "survType"="right")
      y0Vector[i] <- tempStats[["y0"]]
      y1Vector[i] <- tempStats[["y1"]]
      obs0Vector[i] <- tempStats[["obs0"]]
      obs1Vector[i] <- tempStats[["obs1"]]

      if (computeZ) {
        tempResult <- do.call(logrankSingleZ, tempStats)
        oMinEVector[i] <- tempResult[["oMinE"]]
        varVector[i] <- tempResult[["v"]]
      }

      if (computeExactE) {
        tempResult <- logrankSingleEExact(obs0=tempStats[["obs0"]], obs1=tempStats[["obs1"]],
                                          y0=tempStats[["y0"]], y1=tempStats[["y1"]],
                                          theta0=theta0, thetaS=thetaS)
        logP0[i] <- tempResult[["logP0"]]
        logEValueLess[i] <- tempResult[["logEValueLess"]]
        logEValueGreater[i] <- tempResult[["logEValueGreater"]]
      }
    }
  } else {
    stop("Currently, only Surv type of 'right', and 'counting' (left truncated and right censored) supported")
  }

  # Return results -----
  #
  result <- list("nEvents"=nEvents, "stopTimeVector"=stopTimeVector,
                 "y0Vector"=y0Vector, "y1Vector"=y1Vector,
                 "obs0Vector"=obs0Vector, "obs1Vector"=obs1Vector,
                 "z"=NULL, "sumOMinE"=NULL, "sumVarOMinE"=NULL,
                 "oMinEVector"=oMinEVector, "varVector"=varVector,
                 "logP0"=logP0, "logEValueLess"=logEValueLess,
                 "logEValueGreater"=logEValueGreater)

  if (computeZ) {
    sumOMinE <- sum(oMinEVector)
    sumVarOMinE <- sum(varVector)

    result[["sumOMinE"]] <- sumOMinE
    result[["sumVarOMinE"]] <- sumVarOMinE
    result[["z"]] <- sumOMinE/sqrt(sumVarOMinE)
  }

  return(result)
}







# Sampling functions for design ----

#' Simulate stopping times for the exact safe logrank test
#'
#' @inheritParams designSafeLogrank
#'
#' @param hazardRatio numeric that defines the data generating hazard ratio with which data are sampled.
#' @param nMax An integer. Once nEvents hits nMax the experiment terminates, if it didn't stop due to threshold
#' crossing crossing already. Default set to Inf.
#' @author Muriel Felipe Perez-Ortiz and Alexander Ly
#'
#'
#' @return a list with stoppingTimes and breakVector. Entries of breakVector are 0, 1. A 1 represents stopping
#' due to exceeding nMax, and 0 due to 1/alpha threshold crossing, or running out of participants, which implies
#' that the corresponding stopping time is Inf.
#'
#' @export
#'
#' @examples
#' sampleLogrankStoppingTimes(0.7, nSim=10)
sampleLogrankStoppingTimes <- function(hazardRatio, alpha=0.05, alternative = c("twoSided", "less", "greater"),
                                       m0=5e4L, m1=5e4L, nSim=1e3L, groupSizePerTimeFunction = returnOne,
                                       parameter=NULL, nMax=Inf, pb=TRUE) {

  stopifnot(is.null(parameter) || parameter > 0, alpha > 0, alpha <= 1)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  ## Object that will be returned. A sample of stopping times
  stoppingTimes <- breakVector <- integer(nSim)
  eValuesAtEnd <- numeric(nSim)

  if (is.null(parameter))
    thetaS <- if (hazardRatio > 1) 1/hazardRatio else hazardRatio
  else
    thetaS <- parameter

  if (pb)
    pbSafe <- utils::txtProgressBar(style=3, title="Safe test threshold crossing")

  ## Cycle through simulations
  #
  for (sim in seq_along(stoppingTimes)) {
    ## Reset number of individuals in each group
    y0 <- m0
    y1 <- m1

    nEvents <- 0

    logEValueGreater <- 0
    logEValueLess <- 0

    ## Make events happen in each simulation
    for (group in 1:(y0 + y1)) { ## End point
      groupSize <- min(groupSizePerTimeFunction(),
                       y1 + y0) ## cannot sample more subjects than there are

      obs1 <- rLogrank(n=1, y0=y0, y1=y1, obsTotal=groupSize,
                       theta=hazardRatio)

      obs0 <- groupSize - obs1

      ## If we run out of subjects, we never stopped
      if (y1 - obs1 <= 0 || y0 - obs0 <= 0) {
        stoppingTimes[sim] <- Inf
        break()
      }

      tempResults <- logrankSingleEExact(obs0, obs1, y0, y1, thetaS)

      logEValueGreater <- logEValueGreater + tempResults[["logEValueGreater"]]
      logEValueLess <- logEValueLess + tempResults[["logEValueLess"]]

      y0 <- y0 - obs0
      y1 <- y1 - obs1
      nEvents <- nEvents + groupSize


      evidenceNow <- switch(alternative,
                            "less" = exp(logEValueLess),
                            "greater" = exp(logEValueGreater),
                            "twoSided" = 1/2*exp(logEValueGreater) +
                              1/2*exp(logEValueLess))

      # Note(Alexander): If exceeds 1/alpha threshold then reject normally
      #
      if (evidenceNow >= 1/alpha) {
        eValuesAtEnd[sim] <- evidenceNow
        stoppingTimes[sim] <- nEvents
        break()
      }

      # Note(Alexander): If passed maximum number of events stop.
      #   For power calculations if beyond nEvents, then set to Inf, doesn't matter for the quantile
      #
      if (nEvents >= nMax) {
        eValuesAtEnd[sim] <- evidenceNow
        stoppingTimes[sim] <- nEvents
        breakVector[sim] <- 1
        break()
      }
    }

    if (pb)
      utils::setTxtProgressBar(pbSafe, value=sim/nSim, title="Trials")
  }

  result <- list("stoppingTimes"=stoppingTimes, "breakVector"=breakVector,
                 "eValuesAtEnd"=eValuesAtEnd)
  return(result)
}






#' Helper function: Computes the type II error under optional stopping based on the minimal clinically relevant hazard
#' ratio and the maximum number of nEvents.
#'
#' @inheritParams designSafeLogrank
#' @inheritParams sampleLogrankStoppingTimes
#'
#' @return a list which contains at least beta and an adapted bootObject of class  \code{\link[boot]{boot}}.
#' @author Muriel Felipe Perez-Ortiz and Alexander Ly
#' @export
#'
#' @examples
#' computeLogrankBetaFrom(hrMin=0.7, 300, nSim=10)
computeLogrankBetaFrom <- function(hrMin, nEvents, m0=5e4L, m1=5e4L, alpha=0.05,
                                   alternative = c("twoSided", "greater","less"),
                                   nSim=1e3L, nBoot=1e4L, groupSizePerTimeFunction = returnOne,
                                   parameter=NULL, pb=TRUE) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  tempResult <- sampleLogrankStoppingTimes("hazardRatio"=hrMin, "alternative"=alternative, "alpha"=alpha,
                                           "m0"=m0, "m1"=m1, "nSim"=nSim,
                                           "groupSizePerTimeFunction"=groupSizePerTimeFunction,
                                           "nMax"=nEvents, "parameter"=parameter)

  times <- tempResult[["stoppingTimes"]]

  # Note(Alexander): Break vector is 1 whenever the sample path did not stop
  breakVector <- tempResult[["breakVector"]]

  # Note(Alexander): Setting the stopping time to Inf for these paths doesn't matter for the quantile
  times[as.logical(breakVector)] <- Inf

  bootObjBeta <- computeBootObj("values"=times, "objType"="beta", "nPlan"=nEvents, "nBoot"=nBoot)

  result <- list("beta" = bootObjBeta[["t0"]],
                 "bootObjBeta" = bootObjBeta)

  # TODO(Alexander): Batch version here
  #
  eValuesAtEnd <- tempResult[["eValuesAtEnd"]]

  bootObjLogImpliedTarget <- computeBootObj("values"=eValuesAtEnd, "objType"="logImpliedTarget",
                                            "nBoot"=nBoot)

  result[["logImpliedTarget"]] <- bootObjLogImpliedTarget[["t0"]]
  result[["bootObjLogImpliedTarget"]] <- bootObjLogImpliedTarget

  return(result)
}


#' Helper function: Computes the planned sample size based on the minimal clinical relevant hazard ratio,
#' alpha and beta under optional stopping.
#'
#'
#' @inheritParams designSafeLogrank
#' @inheritParams sampleLogrankStoppingTimes
#' @param digits number of significant digits to be used.
#'
#' @return a list which contains at least nEvents and an adapted bootObject of class  \code{\link[boot]{boot}}.
#' @author Muriel Felipe Perez-Ortiz and Alexander Ly
#'
#' @export
#'
#' @examples
#' computeLogrankNEvents(0.7, 0.2, nSim=10)
computeLogrankNEvents <- function(hrMin, beta, m0=50000, m1=50000, alpha=0.05,
                                  alternative = c("twoSided", "greater","less"),
                                  nSim=1e3L, nBoot=1e3L, groupSizePerTimeFunction = returnOne,
                                  nMax=Inf, parameter=NULL, digits = getOption("digits"), pb=TRUE) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  if (is.infinite(nMax)) {
    if (hrMin >= 0.5 && hrMin <= 2) {
      ratio <- m1/m0

      logHazardRatio <- if (alternative=="twoSided") abs(log(hrMin)) else log(hrMin)
      meanDiffMin <- logHazardRatio*sqrt(ratio)/(1+ratio)

      logThetaS <- if (!is.null(parameter)) log(parameter) else NULL

      tempResult <- computeNPlanBatchSafeZ("meanDiffMin"=meanDiffMin, "beta"=beta,
                                           "alpha"=alpha, "alternative"=alternative,
                                           "testType"="oneSample",
                                           "ratio"=ratio, "parameter"=logThetaS)
      nBatch <- tempResult[["nPlan"]]
    } else {
      nBatch <- nMax
    }
  }

  tempResult <- sampleLogrankStoppingTimes(hazardRatio=hrMin, alternative=alternative, alpha=alpha,
                                           m0=m0, m1=m1, nSim=nSim, groupSizePerTimeFunction=groupSizePerTimeFunction,
                                           nMax=nBatch)

  times <- tempResult[["stoppingTimes"]]

  bootObjNEvents  <- computeBootObj("values"=times, "beta"=beta, "objType"="nPlan", "nBoot"=nBoot)

  nEvents <- ceiling(bootObjNEvents[["t0"]])

  bootObjN1Mean <- computeBootObj("values"=times, "objType"="nMean", "nPlan"=nEvents, "nBoot"=nBoot)

  n1Mean <- ceiling(bootObjN1Mean[["t0"]])

  result <- list("nEvents" = nEvents, "bootObjNEvents" = bootObjNEvents,
                 "n1Mean"=n1Mean, "bootObjN1Mean"=bootObjN1Mean,
                 "nEventsBatch"=nBatch)

  return(result)
}

# Auxilary functions ----

#' Randomly samples from a logrank distribution
#'
#' Draws a number of occurrences in group 1 (treatment) out of
#' obsTotal number of occurrences.
#'
#' @param n integer, number of observations to be sampled.
#' @param y0 Size of the risk set of group 0 (Placebo).
#' @param y1 Size of the risk set of group 1 (Treatment).
#' @param obsTotal Total number of observations.
#' @param theta Odds of group 1 over group 0 (treatment over placebo).
#'
#' @return integer representing the number of occurrences in group 1 out of
#' obsTotal number of occurrences.
#'
#' @export
#'
#' @author Muriel Felipe Perez-Ortiz and Alexander Ly
#'
#' @examples
#' rLogrank(y0=360, y1=89, obsTotal=12, theta=3.14)
#'
rLogrank <- function(n=1, y0, y1, obsTotal, theta) {
  BiasedUrn::rFNCHypergeo(nran = n, # number of rv's to generate
                          m1   = y1,# number of balls in 1st group (treatment)
                          m2   = y0,# number of balls in 2nd group (placebo)
                          n    = obsTotal, # number balls sampled
                          odds = theta) # odds of 1st over 2nd group (treatment over placebo)
}


#' Auxiliary function for sampling of the logrank simulations to return the integer 1 event per time.
#'
#' @return 1
#' @export
#'
#' @examples
#' returnOne()
returnOne <- function() {
  1L
}
