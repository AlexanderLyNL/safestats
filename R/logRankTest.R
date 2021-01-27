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
#' @param ciValue numeric is the ciValue-level of the confidence sequence. Default ciValue=0.95
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
#'
#' # Example with left trunctation due to Judith ter Schure
#'
#' enrollment <- 10     # 5 treatment, 5 placebo
#' lambdaC <- 0.03943723
#' hr1 <- 0.5           # hazard ratio between treatment en placebo group
#' fup <- 40            # folow up of 40 days
#' data <- generateSurvData(nP = 5,
#'                          nT = 5,
#'                          lambdaP = lambdaC,
#'                          lambdaT = hr1*lambdaC,
#'                          endTime = fup,
#'                          seed = 2006)
#'
#' # Add different time of randomisation
#' dateRandStart <- as.Date("2020-05-04")
#' dateRandEnd <- as.Date("2020-05-15")
#'
#' set.seed(2005)
#' data$"dateRand" <- sample(seq.Date(from = dateRandStart, to = dateRandEnd, by = "day"),
#'                           size = enrollment, replace = TRUE)
#' data$"dateEvent/LastFup" <- as.Date(data$dateRand + data$time)
#' data$"dateLastFup" <- as.Date("2020-06-15")
#' data$"participantID" <- 1:nrow(data)
#' data$"participantID"[order(data$"dateRand")] <- 1:nrow(data)
#' data <- data[order(data$"dateRand"), ]
#' data$time <- data$"dateEvent/LastFup" - dateRandStart
#'
#' # Add additional complication with multiple events at the same time
#' data$"dateEvent/LastFup"[data$participantID == 3] <-
#'   data$"dateEvent/LastFup"[data$participantID == 5]
#' data$time[data$participantID == 3] <-
#'   data$"dateEvent/LastFup"[data$participantID == 3] - dateRandStart
#' data$dateRand[data$participantID == 3] <-
#'   data$"dateEvent/LastFup"[data$participantID == 3] -
#'   data$time[data$participantID == 3]
#'
#' # Interim analyses events 1 to 4
#'
#' handResult <- c(-0.8164966, -1.4882057, -0.772088, -0.7502141)
#' #'
#' for (i in 1:4) {
#'   calDate <- sort(data$"dateEvent/LastFup")[i]
#'
#'   dataSoFar <- data[data$dateRand < calDate, ]
#'   dataSoFar$dateLastFup <- calDate
#'   dataSoFar$time <- pmin(dataSoFar$time, calDate - dateRandStart)
#'   dataSoFar$status[dataSoFar$"dateEvent/LastFup" > calDate] <- 1
#'   survObj <- survival::Surv(time = dataSoFar$dateRand - dateRandStart,
#'                             time2 = dataSoFar$time,
#'                             event = dataSoFar$status,
#'                             type = "counting")
#'
#'   interimResult <- safeLogrankTest(survObj ~ dataSoFar$group, designObj = designObj)
#'   interimResult
#'
#' # Compare logrank score to calculations by hand
#'   localTest <- round(interimResult$statistic - handResult[i], 7) == 0
#'
#'   if (!localTest)
#'     stop("Computation of the left-truncated logrank z-score is wrong")
#' }
safeLogrankTest <- function(formula, designObj=NULL, h0=1, ciValue=0.95, data=NULL, survTime=NULL,
                            group=NULL, pilot=FALSE, alpha=NULL, alternative=NULL,
                            exact=FALSE, ...) {

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

  survType <- attr(survTime, "type")

  if (survType=="right") {
    survDiffObj <- survival::survdiff(survTime ~ group)
    nEvents <- sum(survDiffObj[["obs"]])
  } else if (survType=="counting") {
    survDiffObj <- computeLogrankZ("survObj"=survTime, "group"=group)
    nEvents <- survDiffObj[["nEvents"]]
  }


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
      theta1 <- unname(exp(designObj[["parameter"]])) # Think of other point nulls, this only works for null being one

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
      nTotal <- n0 <- n1 <- vector("numeric", length(timeOfEvents))

      if (designObj[["alternative"]]=="two.sided" || designObj[["alternative"]]=="less")
        logEValueLess <- vector("numeric", length(timeOfEvents))

      if (designObj[["alternative"]]=="two.sided" || designObj[["alternative"]]=="greater")
        logEValueGreater <- vector("numeric", length(timeOfEvents))

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

      if (designObj[["alternative"]]=="two.sided" || designObj[["alternative"]]=="less") {
        logPLess <- log(BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                             n=casesInGroup0+casesInGroup1, odds=theta1))
        logEValueLess[i] <- logPLess-logP0
      } else if (designObj[["alternative"]]=="two.sided" || designObj[["alternative"]]=="greater") {
        logPGreater <- log(BiasedUrn::dFNCHypergeo(x=casesInGroup1, m1=n1[i], m2=n0[i],
                                             n=casesInGroup0+casesInGroup1, odds=1/theta1))
        logEValueGreater[i] <- logPGreater-logP0
      }

      n0[i+1] <- n0[i]-casesInGroup0
      n1[i+1] <- n1[i]-casesInGroup1
    }

    if (designObj[["alternative"]]=="two.sided" || designObj[["alternative"]]=="less")
      eValueLess <- exp(sum(logEValueLess))

    if (designObj[["alternative"]]=="two.sided" || designObj[["alternative"]]=="greater")
      eValueGreater <- exp(sum(logEValueGreater))

    eValue <- switch(designObj[["alternative"]],
                     "greater"=eValueGreater,
                     "less"=eValueLess,
                     "two.sided"=1/2*eValueGreater+1/2*eValueLess)

    names(eValue) <- "e"

    result <- list("statistic"=eValue, "n"=nEvents, "estimate"=NULL, "eValue"=eValue,
                   "confSeq"=NULL, "estimate"=NULL, "testType"="logrank",
                   "dataName"=dataName, "exact"=TRUE)
    class(result) <- "safeTest"
    result[["designObj"]] <- designObj
  } else {
    nEff <- ratio/(1+ratio)^2*nEvents

    if (survType=="right") {
      coinObj <- coin::logrank_test(survTime ~ group, alternative="two.sided")
      signZ <- sign(unname(coinObj@statistic@standardizedlinearstatistic))

      zStat <- signZ*sqrt(survDiffObj[["chisq"]])
    } else if (survType=="counting") {
      zStat <- survDiffObj[["z"]]
    }

    meanStat <- zStat/sqrt(nEff)

    result <- list("statistic"=zStat, "n"=nEvents, "estimate"=exp(meanStat), "eValue"=NULL,
                   "confSeq"=NULL, "testType"="logrank", "dataName"=dataName)
    class(result) <- "safeTest"

    names(result[["estimate"]]) <-"hazard ratio"

    # Note(Alexander): This is the same as
    #     zStat <- sqrt(nEff)*(meanStat - meanSlog(h0))
    #
    # but to avoid rounding erros zStat is used instead
    zStat <- zStat - sqrt(nEff)*(log(h0))

    eValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=nEff,
                            "n2"=NULL, "alternative"=alternative, "paired"=FALSE, "sigma"=1)

    tempConfSeq <- computeZConfidenceSequence("nEff"=nEff, "meanStat"=meanStat,
                                              "phiS"=abs(designObj[["parameter"]]), "sigma"=1,
                                              "ciValue"=ciValue, "alternative"="two.sided")

    result[["confSeq"]] <- exp(tempConfSeq)

    result[["eValue"]] <- eValue
    result[["designObj"]] <- designObj
    result[["survDiffObj"]] <- survDiffObj

    names(result[["statistic"]]) <- "z"
  }


  names(result[["n"]]) <- "nEvents"

  return(result)
}


#' @describeIn safeLogrankTest Safe Logrank Test based on Summary Statistic Z
#' All provided data (i.e., z-scores) are assumed to be centred on a hazard ratio = 1, thus, log(hr) = 0 ,
#' and the proper (e.g., hypergeometric) scaling is applied to the data, so sigma = 1. The null hypothesis
#' in the design object pertains to the population and is allowed to differ from log(theta) = 0.
#'
#' @param z numeric representing the observed z statistic.
#' @param nEvents numeric > 0, observed number of events.
#' @param dataNull numeric > 0, the null hypothesis corresponding to the z statistics.
#' By default dataNull = 1 representing
#' @param sigma numeric > 0, scaling in the data
#'
#' @return
#' @export
safeLogrankTestStat <- function(z, nEvents, designObj, ciValue=0.95,
                                alternative=c("two.sided", "greater", "less"),
                                dataNull=1, sigma=1) {
  alternative <- match.arg(alternative)

  if (length(z) != length(nEvents))
    stop("The provided number of z-scores and number of events not equal.")

  names(nEvents) <- "nEvents"

  result <- list("statistic"=z, "n"=nEvents, "estimate"=NULL, "eValue"=NULL,
                 "confSeq"=NULL, "testType"="logrank", "dataName"="Logrank z")

  nEff <- designObj[["ratio"]]/(1+designObj[["ratio"]])^2*nEvents

  # Note(Alexander): Assumed the data are centred at log(hr)=0 and standardised,
  # thus, sigma = 1 data scale
  if (length(z)==1) {
    meanStat <- sigma*z/sqrt(nEff)+log(dataNull)
    zStat <- z - sqrt(nEff)/sigma*log(designObj[["h0"]])
  } else {
    meanStat <- sum(sigma*z/sqrt(nEff)+log(dataNull))/sum(nEff)
    zStat <- sqrt(nEff)/sigma*(meanStat - log(designObj[["h0"]]))
  }

  eValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=nEff,
                          "n2"=NULL, "alternative"=alternative, "paired"=FALSE, "sigma"=1)

  tempConfSeq <- computeZConfidenceSequence("nEff"=nEff, "meanStat"=meanStat,
                                            "phiS"=abs(designObj[["parameter"]]), "sigma"=1,
                                            "ciValue"=ciValue,
                                            "alternative"="two.sided")

  result[["confSeq"]] <- exp(tempConfSeq)

  result[["eValue"]] <- eValue
  result[["designObj"]] <- designObj

  names(result[["statistic"]]) <- "z"
  class(result) <- "safeTest"
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
#' @param h0 numeric > 0, represents the null hypothesis, default h0=1.
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
designSafeLogrank <- function(hrMin=NULL, beta=NULL, nEvents=NULL, h0=1,
                              alternative=c("two.sided", "greater", "less"),
                              alpha=0.05, ratio=1, zApprox=TRUE, tol=1e-5, ...) {
  stopifnot(0 < alpha, alpha < 1)

  result <- list()

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
    safeZObj[["h0"]] <- h0
  }

  result <- safeZObj
  names(result[["h0"]]) <- "theta"

  return(result)
}

#' Helper function computes single component of the logrank statistic
#'
#' @param d0 integer, number of observations in the control group
#' @param d1 integer, number of observations in the treatment group
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
#' d0Vector <- c(1, 0, 2, 1, 0)
#' d1Vector <- c(0, 1, 1, 0, 1)
#'
#' varVector <- oMinEVector <-y0Vector
#'
#' for (i in seq_along(y0Vector)) {
#'   tempResult <- logrankSingle(d0=d0Vector[i], d1=d1Vector[i],
#'                               y0=y0Vector[i], y1=y1Vector[i])
#'   oMinEVector[i] <- tempResult[["oMinE"]]
#'   varVector[i] <- tempResult[["v"]]
#' }
#'
#' sum(oMinEVector)/sqrt(sum(varVector))
#'
logrankSingle <- function(d0, d1, y0, y1, ...) {
  dTotal <- d0 + d1
  yTotal <- y0 + y1

  o1 <- d1
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


#' Computes the sufficient statistics needed to compute logrankSingle
#'
#' @param survDataFrame a Surv object converted to a matrix, then to a data.frame
#' @param y0Index vector of integers corresponding to the control group
#' @param y1Index vector of integers corresponding to the treatment group
#' @param timeNow numeric, current time
#' @param timeBefore numeric, previous time
#' @param survType character, either "right" or "counting" (left truncated, right censored)
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list containing at least the following components:
#' \describe{
#'   \item{d0}{number of observations in the control group.}
#'   \item{d1}{number of observations in the treatment group.}
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

  d0 <- length(intersect(eventIndex, y0Index))
  d1 <- length(intersect(eventIndex, y1Index))

  currentIndex <- which(survDataFrame[[timeLabel]] >= timeNow &
                          survDataFrame[[timeLabel]] > timeBefore)

  y0 <- length(intersect(currentIndex, y0Index))
  y1 <- length(intersect(currentIndex, y1Index))

  result <- list("d0"=d0, "d1"=d1, "y0"=y0, "y1"=y1)
  return(result)
}


#' Helper function to computes the logrank statistic for Surv objects of type
#' "right" and "counting" with the hypergeometric variance.
#'
#' This function was created to complement \code{\link[survival]{survdiff}} from the
#' survival package, which is restricted to Surv objects of type "right". Most likely
#' \code{\link[survival]{survdiff}} is much faster
#'
#' @param survObj a Surv object that is either of type
#' @param group a grouping factor with 2 levels
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
computeLogrankZ <- function(survObj, group, ...) {
  result <- list(nEvents=NULL, z=NULL, oMinEVector=NULL,
                 varVector=NULL, stopTimeVector=NULL)

  # Note(Alexander): Get group label information
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

  # survObj <- survival::aeqSurv(survObj)

  survDataFrame <- as.data.frame(as.matrix(survObj))
  survDataFrame[["group"]] <- group

  # TODO(Alexander): This is probably tricky when status is also right
  #
  stopTimeIndeces <- which(survDataFrame[["status"]]==1)

  # Note(Alexander): Coin provides two numbers for each group one
  nEvents <- length(stopTimeIndeces)

  stopTimeVector <- unique(survDataFrame[[timeLabel]][stopTimeIndeces])
  stopTimeVector <- stopTimeVector[order(stopTimeVector)]

  lengthStopTime <- length(stopTimeVector)

  # TODO(Alexander): Add stash check for old time and subset
  #
  if (lengthStopTime > 0) {
    varVector <- oMinEVector <- rep(NA, length = lengthStopTime)
    timeBeforeVector <- c(0, stopTimeVector[1:(lengthStopTime-1)])
  } else {
    # TODO(Alexander)
    # warning("No observations")
    result <- list(n=0, z=0, oMinEVector=NULL, varVector=NULL, stopTimeVector=NULL)
    return(result)
  }

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
      tempResult <- do.call(logrankSingle, tempStats)

      oMinEVector[i] <- tempResult[["oMinE"]]
      varVector[i] <- tempResult[["v"]]
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
      tempResult <- do.call(logrankSingle, tempStats)

      oMinEVector[i] <- tempResult[["oMinE"]]
      varVector[i] <- tempResult[["v"]]
    }
  } else {
    stop("Currently, only Surv type of 'right', and 'counting' (left truncated and right censored) supported")
  }

  sumOMinE <- sum(oMinEVector)
  sumVarOMinE <- sum(varVector)

  result <- list("nEvents"=nEvents, "z"=sumOMinE/sqrt(sumVarOMinE),
                 "oMinEVector"=oMinEVector, "varVector"=varVector,
                 "sumOMinE"=sumOMinE, "sumVarOMinE"=sumVarOMinE,
                 "stopTimeVector"=stopTimeVector)
}
