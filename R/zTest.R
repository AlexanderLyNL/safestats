# Testing fnts ---------

#' Computes E-Values Based on the Z-Statistic
#'
#' Computes e-values using the z-statistic and the sample sizes  based on the test defining parameter phiS.
#'
#' @inheritParams  designSaviZ
#' @inheritParams  saviZTest
#' @param z numeric that represents the observed z-statistic.
#' @param n1 integer that represents the size in a one-sample Z-test, (n2=\code{NULL}). When n2 is not
#' \code{NULL}, this specifies the size of the first sample for a two-sample test.
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
#' \code{NULL} then the z-statistic is assumed to be one-sample.
#' @param parameter numeric > 0, the savi test defining parameter,
#' see \code{\link{matchEParameterWith}} for details.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an e-value.
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#' @export
#'
#' @examples
#' saviZTestStat(z=3, n1=100, parameter=0.4, eType="grow")
#' saviZTestStat(z=3, n1=100, parameter=0.4^2, eType="eGauss")
#' saviZTestStat(z=3, n1=100, parameter=0.4, eType="eCauchy")
saviZTestStat <- function(
    z, n1, n2=NULL, parameter,
    alternative=c("twoSided", "less", "greater"),
    paired=FALSE, sigma=1,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  nEff <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1 else (1/n1+1/n2)^(-1)

  if (eType=="grow") {
    phiS <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=parameter, "alternative"=alternative,
      "esMinName"="phiS")

    if (alternative=="twoSided") { # two-sided
      logEValue <- -nEff*phiS^2/(2*sigma^2)+lcosh(sqrt(nEff)*phiS/sigma*z)
      result <- list("eValue"=exp(logEValue))
    } else { # one-sided
      result <- list("eValue"=exp(sqrt(nEff)*phiS/sigma*z-nEff*phiS^2/(2*sigma^2)))
    }

    if (result[["eValue"]] < 0) {
      warning("Overflow: e-value smaller than 0")
      result[["eValue"]] <- 2^(-15)
    }
  } else if (eType=="eGauss") {
    g <- parameter

    logResult <- -1/2*log(1+nEff*g)+nEff*g*z^2/(2*(1+nEff*g))

    if (alternative=="twoSided") { # two-sided
      result <- list("eValue"=exp(logResult))
    } else { # one-sided
      if (alternative=="greater") {
        result <- list(
          "eValue"=2*exp(logResult)*
            stats::pnorm(-sqrt(g*nEff/(1+nEff*g))*z, lower.tail = FALSE))
      } else if (alternative=="less") {
        result <- list(
          "eValue"=2*exp(logResult)*
            stats::pnorm(-sqrt(g*nEff/(1+nEff*g))*z, lower.tail = TRUE))
      }
    }
  } else if (eType=="eCauchy") {
    kappaG <- parameter

    if (alternative=="twoSided") { # two-sided
      integrand <- function(g) {
        exp(
          - 1/2*log(1+nEff*g)+nEff*g*z^2/(2*(1+nEff*g))
          - 2*log(g) + stats::dgamma(x=1/g, shape=1/2, rate=kappaG^2/2, log=TRUE)
        )
      }

      tempResult <- stats::integrate(integrand, 0, Inf)
      result <- list("eValue"=tempResult[["value"]],
                     "eValueApproxError"=tempResult[["abs.error"]])
    } else if (alternative %in% c("greater", "less")) {
      wantLowerTail <- if (alternative=="greater") FALSE else TRUE

      integrand <- function(g) {
        2*exp(
          - 1/2*log(1+nEff*g)+nEff*g*z^2/(2*(1+nEff*g))
          - 2*log(g) + stats::dgamma(x=1/g, shape=1/2, rate=kappaG^2/2, log=TRUE)
        )*stats::pnorm(-sqrt(g*nEff/(1+nEff*g))*z, lower.tail = wantLowerTail)
      }

      tempResult <- stats::integrate(integrand, 0, Inf)
      result <- list("eValue"=tempResult[["value"]],
                     "eValueApproxError"=tempResult[["abs.error"]])
    }
  } else if (eType=="mom") {
    g <- parameter

    if (alternative=="twoSided") { # two-sided
      logResult <- - 3/2*log(1+nEff*g) +
        log(1+nEff*g/(1+nEff*g)*z^2) + nEff*g/(1+nEff*g)*z^2/2

      result <- list("eValue"=exp(logResult))
    } else if (alternative %in% c("greater", "less")) {
      momIntegrand <- function(delta) {
        2*g^(-1)*delta^2*
          exp(sqrt(nEff)*z*delta-nEff/2*delta^2+
               stats::dnorm(delta, mean=0, sd=sqrt(g), log=TRUE))
      }

      upperBound <- if (alternative=="greater") Inf else 0
      lowerBound <- if (alternative=="greater") 0 else -Inf

      tempResult <- stats::integrate(momIntegrand,
                              lowerBound, upperBound)
      result <- list("eValue"=tempResult[["value"]],
                     "eValueApproxError"=tempResult[["abs.error"]])
    }
  } else if (eType=="imom") {
    tau <- parameter

    someConstant <- if (alternative=="twoSided") 1 else 2

    iMomIntegrand <- function(delta) {
      someConstant*exp(
        sqrt(nEff)*z*delta-nEff/2*delta^2 +
          1/2*log(tau)-lgamma(1/2)-log(delta^2)-tau/(delta^2)
      )
    }

    upperBound <- if (alternative=="less") 0 else Inf
    lowerBound <- if (alternative=="greater") 0 else -Inf

    tempResult <- stats::integrate(iMomIntegrand, lowerBound, upperBound)

    result <- list("eValue"=tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
  }

  return(result)
}

#' Computes E-Relevance Values Based on the Z-Statistic in favour of the alternative over practical equivalence
#'
#' Evidence for practical equivalence requires the e-value to be small, i.e.
#' smaller than alphaRelevance.
#' If the alternative holds true, i.e. meanDiffTrue >= relevanceSize, then
#' there is no more than alphaRelevance probability of ever seeing
#' eRelevance <= alphaRelevance.
#'
#' @rdname saviZTestStat
#' @inheritParams designSaviZ
#' @export
#'
#' @examples
#' # evidence for the alternative over minimal efficacy
#' saviRelevanceZStat(z=3, n1=100, parameter=0.4)
#' # evidence for minimal efficacy over the alternative
#' saviRelevanceZStat(z=0.35, n1=100, parameter=0.4)
saviRelevanceZStat <- function(
    z, n1, n2=NULL, parameter=NULL,
    alternative=c("twoSided", "less", "greater"),
    paired=FALSE, sigma=1, eType="grow", relevanceSize=NULL, ...) {

  if (is.null(parameter) && is.null(relevanceSize))
    stop("No parameter, nor minimal mean difference given")

  if (is.null(parameter) && !is.null(relevanceSize))
    parameter <- abs(relevanceSize)

  alternative <- match.arg(alternative)

  phiS <- parameter

  nEff <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1 else (1/n1+1/n2)^(-1)

  if (alternative %in% c("twoSided", "greater"))
    sPlus0 <- exp(sqrt(nEff)*phiS/sigma*z-nEff*phiS^2/(2*sigma^2))

  if (alternative %in% c("twoSided", "less"))
    sMin0 <- exp(-sqrt(nEff)*phiS/sigma*z-nEff*phiS^2/(2*sigma^2))

  eValue <- switch(alternative,
                   "twoSided"=max(sPlus0, sMin0),
                   "greater"=sPlus0,
                   "less"=sMin0)

  if (eValue < 0) {
    warning("Overflow: e-value smaller than 0")
    eValue <- 2^(-15)
  }

  result <- list("eValue"=eValue)

  return(result)
}

#' Safe Anytime-Valid Z-Test
#'
#' Savi one- and two-sample Z-tests. Takes as input vector(s) of data and a designObj from
#' \code{\link{designSaviZ}}. The function is modelled after \code{\link[stats]{t.test}()}.
#'
#' @aliases savi.z.test
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param paired a logical indicating whether you want the paired Z-test.
#' @param designObj an object obtained from \code{\link{designSaviZ}()}.
#' @param ciValue numeric representing the confidence level.
#' Default ciValue=NULL yields ciValue = 1 - alpha
#' @param maxRoot Used to bound the candidate set of width of the
#' confidence interval, whenever eType="eCauchy"
#' @param formula a formula of the form lhs ~ rhs where lhs
#' is a numeric variable giving the data values and rhs
#' either 1 for a one-sample or paired test or a factor
#' with two levels giving the corresponding groups. If lhs
#' is of class "Pair" and rhs is 1, a paired test is done
#' @param data an optional matrix or data frame (or similar:
#' see \code{\link[stats]{model.frame}()}) containing the variables in
#' the formula. By default the variables are taken from
#' environment(formula).
#' @param subset an optional vector specifying a subset of
#' observations to be used.
#' @param na.action a function which indicates what should
#' happen when the data contain \code{NA}s. Defaults to
#' getOption("na.action").
#' @param sequential a logical indicating whether a sequential
#' analysis should be performed.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class 'saviTest'. An object of class 'saviTest'
#' is a list containing at least the following components:
#'
#' \describe{
#'   \item{statistic}{the value of the test statistic. Here the z-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{eValue}{the e-value of the savi test.}
#'   \item{confSeq}{a savi confidence interval for the mean (difference).}
#'   \item{estimate}{the estimated means or mean (difference) depending on
#'   whether it was a one-sample test or a two-sample test.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "saviDesign" described in \code{\link{designSaviZ}()}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#'
#' ## Examples with simulated data ----
#' set.seed(1)
#' x <- rnorm(30, mean=0)
#' y <- rnorm(40, mean=0)
#'
#' # Because no designObj is specified, a default
#' # designObj is used with meanDiffMin = 1/2 and sigma=1,
#' # which can be thought of as a medium effect size
#' # according to Cohen.
#'
#' res <- saviZTest(x=x, y=y)
#'
#' # By default sequential=TRUE, because length(x) <= 200.
#' # This allows us to visualise the e-value as a function of
#' # the n1 and associated n2, where the ratio of sample sizes,
#' # ratio=n2/n1 is maintained. Here the e-value at n1=6 uses data
#' # x[1:6] and y[1:ceil(ratio*6)]
#' plot(res)
#'
#' # Plots the confidence sequence
#' plot(res, wantConfSeqPlot=TRUE)
#'
#' # See ?designSaviZ for more info
#' # This designObj also allows for
#' # evidence quantification that the
#' # mean difference is minimal clinically
#' # relevant, here, larger than
#' # relevanceSize=meanDiffMin
#' designObj <- designSaviZ(meanDiffMin=0.7, alpha=0.05,
#'                          alternative="twoSided",
#'                          testType="twoSample",
#'                          relevanceTest=TRUE)
#'
#' res <- saviZTest(x=x, y=y, designObj=designObj)
#'
#' plot(res)
#'
#' # Note that the e-value against relevance falls below alphaRelevance
#' # We can reject the hypothesis that the effect is relevantly larger,
#' # larger than designObj$relevanceTestSim$parameter after a sample size of
#' min(which(res$eRelevanceVec <= designObj$relevanceTestSim$alpha))
#'
#' set.seed(2)
#' x <- rnorm(30, mean=0.6)
#' y <- rnorm(40, mean=0)
#'
#' res <- saviZTest(x, y, designObj=designObj)
#'
#' plot(res)
#' # We could have stopped sampling after
#' min(which(res$eValueVec >= 1/designObj$alpha))
#' # the yellow curve crosses 1/alpha sooner,
#' # but we do **not** compare eRelevance >= 1/alpha, only eRelevance <= alphaRelevance
#'
#' ## Classical example: Student's sleep data -----
#' plot(extra ~ group, data = sleep)
#'
#' designObj <- designSaviZ(meanDiffMin=0.6, sigma=2,
#'                          testType="twoSample")
#'
#' ## Traditional interface
#' with(sleep, saviZTest(extra[group == 1], extra[group == 2],
#'                       designObj=designObj))
#'
#' ## Formula interface
#' saviZTest(extra ~ group, data = sleep, designObj=designObj)
#'
#' ## Formula interface to one-sample test
#' designObj1 <- designSaviZ(meanDiffMin=0.6,
#'                           testType="oneSample",
#'                           sigma=2)
#'
#' saviZTest(extra ~ 1, data = sleep, designObj=designObj1)
#'
#' ## Formula interface to paired test
#' ## The sleep data are actually paired, so could have been in wide format:
#' designObjPaired <- designSaviZ(meanDiffMin=0.6,
#'                                testType="paired",
#'                                sigma=1.4)
#' sleep2 <- reshape(sleep, direction = "wide",
#'                   idvar = "ID", timevar = "group")
#' saviZTest(Pair(extra.1, extra.2) ~ 1, data = sleep2,
#'           designObj=designObjPaired)
saviZTest <- function(x, ...) {
  UseMethod("saviZTest")
}

#' @rdname saviZTest
#' @aliases saviZTest
#' @export
saviZTest.default <- function(
    x, y=NULL, paired=FALSE, designObj=NULL,
    ciValue=NULL, maxRoot=10, sequential=NULL,
    ...) {

  result <- constructSaviTestObj("Z-Test")

  eRelevance <- NULL

  # Vars for sequential analysis
  eValueVec <- NULL
  eRelevanceVec <- NULL
  confSeqMatrix <- NULL
  n1Vec <- NULL
  n2Vec <- NULL
  fpt <- NULL
  fptRelevance <- NULL

  ### Def: test type -------
  if (is.null(y)) {
    testType <- "oneSample"

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    dataName <- deparse1(substitute(x))
  } else {
    testType <- if (paired) "paired" else "twoSample"

    dataName <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  }

  ### Check: designObj ----
  if (is.null(designObj)) {
    designObj <- designSaviZ(0.5, "eType"="mom",
                             "testType"=testType)
    designObj[["pilot"]] <- TRUE

    warningMessage <- paste("No designObj given. Default test computed based",
                            "on a non-local moment prior at +1/2 and - 1/2.")
    warning(warningMessage)
  }

  if (designObj[["testName"]] != "Z-Test")
    warning("The provided design is not constructed for the Z-test,",
            "please use designSaviZ() instead. The test results might be invalid.")

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  sumStats <- computeZTSumStats(
    "x"=x, "y"=y, "sequential"=sequential,
    "varEqual"=TRUE, "paired"=paired,
    "testType"=testType)

  nEff <- meanObs <- n1 <- n2 <- nEffVec <- meanObsVec <- estimate <-
    n <- sdObsVec <- nuVec <- NULL

  list2env(sumStats, envir=environment())

  alpha <- designObj[["alpha"]]
  sigma <- designObj[["sigma"]]
  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]

  if (is.null(ciValue))
    ciValue <- 1-alpha

  if (ciValue < 0 || ciValue > 1)
    stop("Can't make a confidence sequence with ciValue < 0 or ciValue > 1, or alpha < 0 or alpha > 1")

  zStat <- tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sigma)

  if (is.na(zStat))
    stop("Could not compute the z-statistic")

  names(zStat) <- "z"

  ### Compute: eValue ----
  tempResult <- saviZTestStat("z"=zStat, "parameter"=designObj[["parameter"]],
                              "n1"=n1, "n2"=n2, "sigma"=sigma,
                              "alternative"=alternative, "paired"=paired,
                              "eType"=designObj[["eType"]])

  if (designObj[["relevanceTest"]]) {
    tempResultFut <- saviRelevanceZStat(
      "z"=zStat, "parameter"=designObj[["relevanceTestSim"]][["parameter"]],
      "n1"=n1, "n2"=n2, "sigma"=sigma,
      "alternative"=alternative, "paired"=paired,
      "eType"="grow")
    eRelevance <- unname(tempResultFut[["eValue"]])
  }

  ### Compute: confSeq ----
  result[["confSeq"]] <- computeConfidenceIntervalZ(
    "nEff"=nEff, "meanObs"=meanObs, "parameter"=designObj[["parameter"]],
    "sigma"=sigma, "ciValue"=ciValue, "alternative"="twoSided",
    "eType"=designObj[["eType"]], "maxRoot"=maxRoot)

  ### Compute: Sequential ----
  if (sequential) {
    zStatVec <- sqrt(nEffVec)*(meanObsVec-h0)/sigma

    mIter <- length(n1Vec)

    eRelevanceVec <- eValueVec <- numeric(mIter)
    confSeqMatrix <- matrix(nrow=mIter, ncol=2)

    for (i in seq_along(n1Vec)) {
      res <- saviZTestStat("z"=zStatVec[i],
                           "parameter"=designObj[["parameter"]],
                           "n1"=n1Vec[i], "n2"=n2Vec[i], "sigma"=sigma,
                           "alternative"=alternative, "paired"=paired,
                           "eType"=designObj[["eType"]])
      eValueVec[i] <- unname(res[["eValue"]])

      if (designObj[["relevanceTest"]]) {
        relevanceRes <- saviRelevanceZStat("z"=zStatVec[i], "n1"=n1Vec[i], "n2"=n2Vec[i],
                                    "parameter"=designObj[["relevanceTestSim"]][["parameter"]],
                                    "sigma"=sigma, "alternative"=designObj[["alternative"]],
                                    "paired"=paired, "eType"="grow")
        eRelevanceVec[i] <- unname(relevanceRes[["eValue"]])
      }

      kaas <- computeConfidenceIntervalZ(
        "nEff"=nEffVec[i], "meanObs"=meanObsVec[i], "parameter"=designObj[["parameter"]],
        "sigma"=sigma, "ciValue"=ciValue, "alternative"="twoSided",
        "eType"=designObj[["eType"]], "maxRoot"=maxRoot)

      confSeqMatrix[i, ] <- kaas
    }

    tempConfSeq <- c(max(confSeqMatrix[, 1]), min(confSeqMatrix[, 2]))

    if (tempConfSeq[1] >= tempConfSeq[2]) {
      warning("Possible high degree of heterogeneity",
              "leading to an empty running intersection confidence sequence")
    } else if (tempConfSeq[1] < tempConfSeq[2]) {
      result[["confSeq"]] <- tempConfSeq
    }

    fpt <- suppressWarnings(
      min(which(eValueVec >= 1/designObj[["alpha"]]))
    )

    if (designObj[["relevanceTest"]])
      fptRelevance <- suppressWarnings(
        min(which(eRelevanceVec <= designObj[["relevanceTestSim"]][["alpha"]]))
      )
  }

  ### Fill: Result -----
  result[["testType"]] <- testType
  result[["statistic"]] <- zStat
  result[["estimate"]] <- estimate
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["ciValue"]] <- ciValue
  result[["n"]] <- n

  result[["eRelevance"]] <- eRelevance

  result[["eValueVec"]] <- eValueVec
  result[["eRelevanceVec"]] <- eRelevanceVec
  result[["confSeqMatrix"]] <- confSeqMatrix
  result[["n1Vec"]] <- n1Vec
  result[["n2Vec"]] <- n2Vec

  # result[["eType"]] <- eType

  result[["eValue"]] <- tempResult[["eValue"]]
  result[["eValueApproxError"]] <- tempResult[["eValueApproxError"]]

  # Sumstats
  result[["meanObs"]] <- meanObs
  result[["sdObs"]] <- meanObs
  result[["meanObsVec"]] <- meanObsVec
  result[["sdObsVec"]] <- sdObsVec
  result[["nEffVec"]] <- nEffVec
  result[["nuVec"]] <- nuVec
  result[["fpt"]] <- fpt
  result[["fptRelevance"]] <- fptRelevance

  names(result[["statistic"]]) <- "z"

  return(result)
}

#' @rdname saviZTest
#' @aliases saviZTest
#' @export
saviZTest.formula <- function(
    formula, data, subset, na.action, ...) {
  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  wantTwoSample <- TRUE

  if (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L)
    if (formula[[3L]] == 1L)
      wantTwoSample <- FALSE
  else
    stop("'formula' missing or incorrect")

  matchedCall <- match.call(expand.dots = FALSE)

  if (is.matrix(eval(matchedCall[["data"]], parent.frame())))
    matchedCall[["data"]] <- as.data.frame(data)

  # Note: Prepare calling stats::model.frame instead of saviTTest
  #
  matchedCall[[1L]] <- quote(stats::model.frame)
  matchedCall[["..."]] <- NULL

  # Call: stats::model.frame
  #
  modelFrame <- eval(matchedCall,
                     parent.frame())

  # Naming
  dataName <- paste(names(modelFrame),
                    collapse=" by ")

  names(modelFrame) <- NULL
  response <- attr(attr(modelFrame, "terms"),
                   "response")

  if (isTRUE(wantTwoSample)) {
    groupingFactor <- factor(modelFrame[[-response]])

    if (nlevels(groupingFactor) != 2L)
      stop("grouping factor must have exactly 2 levels")

    dataList <- split(modelFrame[[response]], groupingFactor)

    result <- saviZTest("x"=dataList[[1L]], "y"=dataList[[2L]], ...)

    if (length(result[["estimate"]]) == 2L) {
      names(result[["estimate"]]) <- paste("mean in group", levels(groupingFactor))
      names(result[["designObj"]][["h0"]]) <-
        paste("true difference in means between",
              paste("group", levels(groupingFactor), collapse = " and "))
    }
  } else {
    respVar <- modelFrame[[response]]

    if (inherits(respVar, "Pair")) {
      result <- saviZTest("x"=respVar[, 1L], "y"=respVar[, 2L],
                          paired=TRUE, ...)
      firstVar <- substring(dataName,
                            first=6,
                            last=regexpr(",", dataName)-1)
      secondVar <- substring(dataName,
                             first=regexpr(",", dataName)+2,
                             last=regexpr(")", dataName)-1)
      names(result[["estimate"]]) <-
        paste("mean difference between", firstVar, "and", secondVar)
      names(result[["designObj"]][["h0"]]) <-
        paste("true mean difference between",
              paste(c(firstVar, secondVar), collapse = " and "))
    } else {
      result <- saviZTest("x"=respVar, "y"=NULL, ...)
    }
  }

  result[["dataName"]] <- dataName
  return(result)
}

#' @rdname saviZTest
#' @aliases saviZTest
#' @export
savi.z.test <- function(x, y=NULL, paired=FALSE,
                        designObj=NULL, ...) {

  result <- saviZTest("x"=x, "y"=y, "designObj"=designObj,
                      "paired"=paired, ...)
  argumentNames <- getArgs()

  xLabel <- extractNameFromArgs(argumentNames, "x")

  if (is.null(y)) {
    dataName <- xLabel
  } else {
    yLabel <- extractNameFromArgs(argumentNames, "y")
    dataName <- paste(xLabel, "and", yLabel)
  }

  result[["dataName"]] <- dataName
  return(result)
}

#' Helper function: Computes the savi confidence sequence for a Z-test
#'
#' @inheritParams saviZTestStat
#' @inheritParams designSaviZ
#' @inheritParams saviZTest
#'
#' @param nEff numeric > 0, the effective sample size.
#' @param meanObs numeric, the observed mean. For two sample tests this is difference of the means.
#' @param a numeric, only used for eType="credibleInterval". a specifies the centre of the Gaussian
#' prior on population mean.
#' @param g numeric > 0, only used for eType="credibleInterval". g specifies the variance of
#' the Gaussian prior on population mean.
#' @param eType character string, one of "eCauchy", "eGauss", "grow", "freq", "credibleInterval".
#' This is somewhat a misnomer as "freq" and "credibleInterval" do not correspond to e-value tests.
#' "eCauchy" yields an anytime-valid confidence interval based on a Cauchy mixture,
#' whereas "eGauss" yields an anytime-valid confidence interval based on a Gaussian mixture.
#' "grow" yields an anytime-valid confidence interval based on a mixture of point masses
#' at the minimal clinically relevant mean difference. This confidence interval unfortunately
#' does not shrink as the sample size tends to infinity. "freq" yields the standard unsafe frequentist
#' confidence interval, which is not safe. "credibleInterval" yields the Bayesian credible interval
#' based on a conjugate prior as is usual in Bayesian analysis. This interval is also not safe.
#'
#' @return numeric vector that contains the upper and lower bound of the savi confidence sequence
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' computeConfidenceIntervalZ(nEff=15, meanObs=0.3, parameter=0.2)
computeConfidenceIntervalZ <- function(
    nEff, meanObs, parameter, sigma=1, ciValue=0.95,
    alternative="twoSided", a=NULL, g=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow",
            "freq", "credibleInterval"),
    maxRoot=20) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  eType <- match.arg(eType)
  alpha <- 1-ciValue

  if (eType=="freq") {
    width <- sigma/sqrt(nEff)*
      stats::qnorm(alpha/2, lower.tail=FALSE)
  }

  if (eType=="credibleInterval") {
    normalisedCiLower <- stats::qnorm(alpha/2, lower.tail=FALSE)

    if (is.null(a)) {
      warning("Centre of normal prior not given, default to a=0")
      a <- 0
    }

    if (is.null(g)) {
      warning("Variance of normal prior not given, default to g=1")
      g <- 1
    }

    meanMu <- nEff*g/(1+nEff*g)*meanObs + 1/(1+nEff*g)*a
    sdMu <- sqrt(g*sigma^2)/sqrt(1+nEff*g)

    width <- sdMu*normalisedCiLower
  }

  ## Confidence sequences ----

  trivialConfInt <- switch(alternative,
                           "twoSided"=c(-Inf, Inf),
                           "greater"=c(0, Inf),
                           "less"=c(-Inf, 0))

  g <- parameter

  W <- (1+nEff*g)/(nEff*g)*(log(1+nEff*g)-2*log(alpha))

  if (eType=="eGauss") {
    if (!is.null(a) && !is.null(g)) {
      # Note(Alexander): Here normal distribution not centred at null
      if (alternative != "twoSided")
        stop("One-sided confidence sequences for non-zero centred normal priors not implemented.")

      width <- sqrt(sigma^2/nEff*
                      (log(1+nEff*g)-2*log(alpha))+(meanObs-a)^2/(1+nEff*g))
    } else {
      # Note(Alexander): Two-sided normal priors centred at the null
      # One-sided handled numerically

      if (alternative=="twoSided") {
        if (W < 0) return(trivialConfInt)

        width <- sigma/sqrt(nEff)*sqrt(W)
      }
    }
  }

  if (eType=="grow") {
    phiS <- parameter

    if (alternative=="twoSided") {
      acoshTerm <- acosh(exp(nEff*phiS^2/(2*sigma^2))/alpha)

      if (is.infinite(acoshTerm))
        acoshTerm <- log(2)+nEff*phiS^2/(2*sigma^2)-log(alpha)

      width <- sigma^2/(nEff*phiS)*acoshTerm
      # acosh(exp(nEff*phiS^2/(2*sigma^2))/(alpha))
    } else if (alternative %in% c("greater", "less")) {
      width <- sigma^2/(nEff*abs(phiS))*
        log(1/alpha)+abs(phiS)/2
    }
  }

  if (eType=="mom" && alternative=="twoSided") {
    # TODO(Alexander): Perhaps add check for parameter then meanDiffMin
    #
    #   meanDiffMin easier to work with
    #
    g <- parameter

    width <- sigma/sqrt(nEff)*
      sqrt((1+1/(nEff*g))*
             (2*lamW::lambertW0((1+nEff*g)^(3/2)*exp(1/2)/(2*(alpha)))-1))
  }

  if (eType %in% c("eCauchy", "imom") ||
      (eType %in% c("eGauss", "mom") && alternative!="twoSided")) {

    # Note(Alexander): The target function is for the two-sided test,
    # which is why we consider alpha/2
    # Alternatively, we can consider tempAlternative being "twoSided" or "greater"
    if (alternative=="twoSided") {
      ciLogPenaltyFunc <- function(ciValue) 1/(1-ciValue)
    } else if (alternative %in% c("greater", "less")) {
      ciLogPenaltyFunc <- function(ciValue) 1/(2*(1-ciValue))
    }

    targetFunction <- function(z) {
      saviZTestStat(z, "n1"=nEff, "parameter"=parameter, "sigma"=sigma,
                    alternative="twoSided", # this is why we consider 2*(alpha)
                    "eType"=eType)$eValue-ciLogPenaltyFunc(ciValue)
    }

    if (W > 0) {
      lowerB <- max(sqrt(W)-5,0)
      upperB <- max(sqrt(W)+5, 5)
    } else {
      lowerB <- 0
      upperB <- maxRoot
      maxRoot <- 2*maxRoot
    }

    tempResult <- suppressWarnings(
      tryCatch(stats::uniroot(targetFunction, c(lowerB, upperB)),
               error=identity)
    )

    iterationN <- 1

    while ( (inherits(tempResult, "simpleError") || is.null(tempResult)) && iterationN <= 15) {
      iterationN <- iterationN+1

      tempResult <- suppressWarnings(
        tryCatch(stats::uniroot(targetFunction, c(0, maxRoot)),
                 error=identity)
      )

      maxRoot <- maxRoot*2
    }

    if (inherits(tempResult, "simpleError")) {
      return(trivialConfInt)
      # stop("Can't compute the width of the interval")
    }

    width <- sigma/sqrt(nEff)*tempResult$root
  }

  if (alternative=="twoSided") {
    lowerCS <- meanObs - width
    upperCS <- meanObs + width
  } else if (alternative=="greater") {
    lowerCS <- meanObs + width
    upperCS <- Inf
  } else if (alternative=="less") {
    lowerCS <- -Inf
    upperCS <- meanObs - width
  }

  return(unname(c(lowerCS, upperCS)))
}



# Design fnts ---------

#' Design a Frequentist Z-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I
#' and desired power for the frequentist Z-test.
#'
#' @inheritParams designSaviZ
#' @param lowN integer that defines the smallest n of our search space for n.
#' @param highN integer that defines the largest n of our search space for n. This might be the
#' largest n that we are able to fund.
#'
#' @return returns a 'freqZDesign' object.
#' @export
#'
#' @examples
#' freqDesign <- designFreqZ(meanDiffMin = 0.5, highN = 100)
#' freqDesign$nPlan
#' freqDesign2 <- designFreqZ(meanDiffMin = 0.2, lowN = 32, highN = 200)
#' freqDesign2$nPlan
designFreqZ <- function(
    meanDiffMin, alternative=c("twoSided", "greater", "less"),
    alpha=0.05, power=0.8, testType=c("oneSample", "paired", "twoSample"),
    ratio=1, sigma=1, h0=0, kappa=sigma, lowN=3L, highN=100L, ...) {

  beta <- 1-power

  stopifnot(lowN >= 1, highN > lowN, alpha > 0, power >0, power < 1)

  testType <- match.arg(testType)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  result <- list("nPlan"=NA, "esMin"=meanDiffMin, "alpha"=alpha, "power"=power,
                 "lowN"=lowN, "highN"=highN, "testType"=testType, "alternative"=alternative,
                 "h0"=h0)
  class(result) <- "freqZDesign"

  meanDiffMin <- checkAndReturnEsMinParameterSide(
    "paramToCheck"=meanDiffMin, "alternative"=alternative,
    "esMinName"="meanDiffMin")

  n1Plan <- NULL
  n2Plan <- NULL

  if (alternative=="twoSided")
    threshold <- 1-alpha/2
  else if (alternative %in% c("greater", "less"))
    threshold <- 1-alpha

  for (n in seq.int(lowN, highN)) {
    nEff <- if (testType=="twoSample") ratio/(1+ratio)*n else n

    powerZ <- stats::pnorm(stats::qnorm(threshold, mean=0, sd=kappa/sigma),
                           mean=sqrt(nEff)*(meanDiffMin)/sigma, sd=kappa/sigma, lower.tail=FALSE)

    if (powerZ >= (1-beta)) {
      n1Plan <- n

      if (testType=="twoSample")
        n2Plan <- ceil(ratio*n)

      if (testType=="paired")
        n2Plan <- n

      break()
    }
  }


  if (is.null(n1Plan) || is.na(n1Plan))
    return(result)

  if (is.null(n2Plan) || is.na(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  return(result)
}

#' Design a Safe Anytime-Valid Experiment to Test Means with a Z Test
#'
#' A designed experiment requires (1) a sample size nPlan to plan for, and
#' (2) a savi test defining parameter. The design involves alpha and the
#' three quantities: (1) nPlan, (2) power, and (3) a minimal
#' clinically relevant mean difference meanDiffMin.
#' \describe{
#'   \item{Scenario 1.a}{Goal: "nPlan" and optimal E-variable. Given: meanDiffMin and power.}
#'   \item{Scenario 1.b}{Goal: an optimal E-variable. Given: meanDiffMin only.}
#'   \item{Scenario 2}{Goal: "power" and optimal E-variable. Given: meanDiffMin and nPlan.}
#'   \item{Scenario 3.a}{Goal: "meanDiffMin" and optimal E-variable. Given: power and nPlan.}
#'   \item{Scenario 3.b}{Goal: an optimal E-variable. Given: nPlan only.}
#'}
#'
#' Every scenario returns an E-variable adapted to the input. Scenario 1.a,
#' for instance, outputs the parameter of the provided eType (default mom)
#' savi test, see \code{\link{matchEParameterWith}} for details, and nPlan.
#' The nPlan is based on samples paths drawn under meanDiffTrue (if not specified,
#' then meanDiffTrue=meanDiffMin by default). The resulting nPlan corresponds to the
#' power (say 80%) quantile of the first-passage time distribution associated
#' with E crossing threshold 1/alpha.
#'
#' @param meanDiffMin numeric that defines the minimal relevant mean difference,
#' the smallest population mean difference that we would like to detect (with sufficient power).
#' @param power numeric in (0, 1) that specifies the desired power, that is, the targetted
#' chance to stop in favour of the alternative over the null hypothesis, when the alternative
#' holds true. Note that prior to version 0.8.8 power <- 1-beta. The "beta" argument does not
#' need to be specified anymore.
#' @param nPlan optional numeric vector of length at most 2, see scenario 2 and 3 above.
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error and the null
#' rejection rule e >= 1/alpha.
#' @param h0 numeric, representing the null value, default h0=0.
#' @param alternative a character string specifying the alternative hypothesis. Must be one
#' of "twoSided" (default), "greater" or "less".
#' @param sigma numeric > 0 representing the assumed population standard deviation used to
#' scale the data.
#' @param kappa the true population standard deviation. Default kappa=sigma.
#' @param meanDiffTrue numeric, data governing mean difference used for simulations. Default
#' meanDiffTrue=meanDiffMin.
#' @param beta numerical in (0,1). Old parameter now replaced by the power parameter
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1.
#' If testType is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param parameter numeric, an optional savi test defining parameter. Default set to \code{NULL}.
#' and adapts to meanDiffMin and eType, see \code{\link{matchEParameterWith}} for details.
#' @param eType character one of "mom", "grow", "eGauss", and "eCauchy". "mom" is default
#' and uses a non-local moment prior with bump(s) at meanDiffMin, "grow" uses point prior(s) at
#' meanDiffMin, "eGauss" a zero-centred normal prior, "eCauchy" a zero centred Cauchy prior.
#' @param wantSamplePaths logical, if \code{TRUE} then also outputs the sample paths.
#' @param wantSimData logical, if \code{TRUE} then also output the simulated data
#' @param lowEsTrue numeric, lower bound for the candidate set of the
#' targeted minimal clinically relevant effect size for scenario 3.a.
#' @param highEsTrue numeric, upper bound for the candidate set of the
#' targeted minimal clinically relevant effect size for scenario 3.a.
#' @param pb logical, if \code{TRUE}, then show progress bar.
#' @param seed integer, seed number.
#' @param nSim integer > 0, the number of simulations needed to compute power or the number of
#' samples paths for the savi z test under continuous monitoring.
#' @param nBoot integer > 0 representing the number of bootstrap samples to assess the accuracy
#' of the approximations of the power, the number of samples for the savi z test under continuous
#' monitoring,or for the computation of the logarithm of the implied target.
#' @param relevanceTest logical, if \code{TRUE} then impose a rule to stop
#' for minimal efficiency if e <= alphaRelevance. Default \code{FALSE}.
#' @param relevanceSize numeric, the minimal clinical relevant mean
#' difference that we do not want to miss under the alternative.
#' Default relevanceSize=NULL implies relevanceSize=abs(meanDiffMin)
#' @param alphaRelevance numeric, the threshold for relevance test. Taken to be minimum of
#' alpha and 1-power.
#' @param betaDefault numeric, defaulting value for 1-power and alphaRelevance
#' @param highN integer, largest possible sampling horizon. This might be the
#' largest n that we are able to fund, which by default is set to 1e4L.
#' Typically, highN is not used, as the function
#' `computeNPlanBatchSaviZ()` tries to find the sampling horizon.
#' If all fails, then use highN as the sampling horizon.
#' @param wantSampling logical, default TRUE so sampling paths are drawn.
#' For instance, if meanDiffMin and power, are given, then nPlan
#' (scenario 1a) is derived by sampling. Set this to FALSE, whenever we
#' want to run a minimal efficacy test without needing to know nPlan
#' @param ... further arguments to be passed to or from methods.
#'
#'
#' @return Returns a saviDesign object that includes:
#'
#' \describe{
#'   \item{parameter}{the savi test defining parameter, see \code{\link{matchEParameterWith}}.}
# #'   \item{esMin}{the minimally clinically relevant effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
# #'   \item{power}{the desired power specified by the user.}
# #'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
# #'   \item{testType}{any of "oneSample", "paired", "twoSample" effectively provided by the user.}
# #'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
# #'   \item{sigma}{the assumed population standard deviation used for the test provided by the user.}
# #'   \item{kappa}{the true population standard deviation, typically, sigma=kappa.}
# #'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
# #'   ratio between the planned randomisation of condition 2 over condition 1.}
#'   \item{pilot}{logical, specifying whether it's a pilot design, which occurs
#'   when saviZTest is called without a designObj.}
#'   \item{testName}{"Z-Test".}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @examples
#' # Scenario 1.b: Goal: an E-variable
#' designObj <- designSaviZ(meanDiffMin=0.8)
#'
#' # Scenario 1.a: Goal: "nPlan" and optimal E-variable.
#' designObj <- designSaviZ(meanDiffMin=0.8, power=0.6, alpha=0.2,
#'                          alternative="greater", nSim=100)
#'
#' plot(designObj)
#'
#' # Scenario 1a. with relevance testing, also stopping for practically null
#' designObj <- designSaviZ(meanDiffMin=0.8, power=0.6, alpha=0.2,
#'                          alternative="greater", nSim=100,
#'                          relevanceTest=TRUE)
#'
#' plot(designObj)
#'
#' # Scenario 2: Goal: "power" and optimal E-variable
#' designObj <- designSaviZ(meanDiffMin=0.8, nPlan=16, nSim=100)
#'
#' # Scenario 3.a: Goal: "meanDiffMin" and optimal E-variable
#' designObj <- designSaviZ(power=0.7, nPlan=16)
#'
#' # Scenario 3.b: Goal: an optimal E-variable. Given: nPlan only.
#' designObj <- designSaviZ(nPlan=16)
#'
designSaviZ <- function(
    meanDiffMin=NULL, power=NULL, nPlan=NULL,
    alpha=0.05, h0=0, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    meanDiffTrue=NULL, beta=NULL,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    lowEsTrue=0.01, highEsTrue=3,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    relevanceTest=FALSE, relevanceSize=NULL,
    alphaRelevance=NULL, betaDefault=0.2,
    highN=1e4L, wantSampling=TRUE, ...) {

  stopifnot(alpha > 0, alpha < 1, sigma > 0, kappa > 0)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  result <- constructSaviDesignObj("Z-Test")

  # TODO(Alexander): Check this:
  #
  #     checkAndReturnEsMinParameterSide("paramToCheck"=0.4, "esMinName"="kappaG", alternative="greater")
  #
  if (!is.null(parameter)) {
    if (eType=="grow") {
      parameter <- checkAndReturnEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="phiS",
        "alternative"=alternative)
    } else if (eType=="eGauss") {
      parameter <- checkAndReturnEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="g",
        "alternative"=alternative)
    } else if (eType=="eCauchy") {
      parameter <- checkAndReturnEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="kappaG",
        "alternative"=alternative)
    }
  }

  if (!is.null(meanDiffMin)) {
    meanDiffMin <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=meanDiffMin, "esMinName"="meanDiffMin",
      "alternative"=alternative)

    if (is.null(meanDiffTrue))
      meanDiffTrue <- meanDiffMin

    parameter <- matchEParameterWith(
      "parameter"=parameter, "analysisType"="z",
      "esMin"=meanDiffMin, "sigma"=sigma,
      "alternative"=alternative, "eType"=eType)
  }

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (is.null(beta) && !is.null(power))
    beta <- 1-power

  if (relevanceTest) {
    relevanceSize <- matchRelevanceParameterWith(
      "relevanceSize"=relevanceSize, "esMin"=meanDiffMin,
      "esTrue"=meanDiffTrue)

    if (is.null(relevanceSize))
      stop("Can't run a minimal efficacy analysis without relevanceSize or meanDiffMin")

    alphaRelevance <- matchAlphaRelevanceWith(
      "alphaRelevance"=alphaRelevance, "alpha"=alpha,
      "power"=power, "betaDefault"=betaDefault)
  }

  designScenario <- NULL

  tempResult <- list()

  if (!is.null(meanDiffMin) && !is.null(power) && is.null(nPlan) && wantSampling) {
    designScenario <- "1a"

    tempResult <- designSaviZ1aWantNPlan(
      "meanDiffMin"=meanDiffMin, "power"=power,
      "alpha"=alpha, "alternative"=alternative,
      "meanDiffTrue"=meanDiffTrue,
      "sigma"=sigma, "kappa"=kappa, "ratio"=ratio,
      "parameter"=parameter, "testType"=testType,
      "eType"=eType, "wantSamplePaths"=wantSamplePaths,
      "wantSimData"=wantSimData,
      "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot,
      "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
      "highN"=highN, "alphaRelevance"=alphaRelevance, ...)

  } else if (!is.null(meanDiffMin) && !is.null(power) && is.null(nPlan) && isFALSE(wantSampling) ||
             !is.null(meanDiffMin) && is.null(power) && is.null(nPlan)) {
    designScenario <- "1b"

    tempResult <- list("parameter"=parameter,
                       "esMin"=meanDiffMin, "relevanceTest"=relevanceTest)

    if (relevanceTest) {
      relevanceParameter <- matchRelevanceParameterWith(
        "relevanceSize"=relevanceSize,
        "esMin"=meanDiffMin, "esTrue"=meanDiffTrue)

      relevanceTestSim <- list("parameter"=relevanceParameter,
                             "alpha"=alphaRelevance)

      tempResult[["relevanceTestSim"]] <- relevanceTestSim
    }

  } else if (!is.null(meanDiffMin) && is.null(power) && !is.null(nPlan)) {
    designScenario <- "2"

    tempResult <- designSaviZ2WantPower(
      "meanDiffTrue"=meanDiffTrue, "nPlan"=nPlan, "alpha"=alpha,
      "sigma"=sigma, "kappa"=kappa, "alternative"=alternative,
      "testType"=testType, "parameter"=parameter, "meanDiffMin"=meanDiffMin,
      "eType"=eType, "wantSamplePaths"=wantSamplePaths,
      "ratio"=ratio, "wantSimData"=wantSimData,
      "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot,
      "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
      "alphaRelevance"=alphaRelevance, ...)
  } else if (is.null(meanDiffMin) && !is.null(power) && !is.null(nPlan)) {
    designScenario <- "3"

    tempResult <- designSaviZ3WantEsMin(
      "power"=power, "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative,
      "sigma"=sigma, "kappa"=kappa,
      "testType"=testType,
      "parameter"=parameter, "eType"=eType,
      "lowEsTrue"=lowEsTrue, "highEsTrue"=highEsTrue, ...)
  } else if (is.null(meanDiffMin) && is.null(power) && !is.null(nPlan)) {
    #scenario 3b: only nPlan known, find the parameter at which the confidence interval
    # is the most narrow at nPlan

    designScenario <- "3b"

    tempResult <- designSaviZ3bWantParameter(
      "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative, "testType"=testType,
      "parameter"=parameter, "eType"=eType)
  }

  if (is.null(designScenario)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1.a) non-null meanDiffMin, non-null power and NULL nPlan, or \n",
         "(1.b) non-null meanDiffMin, NULL power, and NULL nPlan, or \n",
         "(2) non-null meanDiffMin, NULL power and non-null nPlan, or \n",
         "(3) NULL meanDiffMin, non-null power, and non-null nPlan, or \n",
         "(3.b) NULL meanDiffMin, NULL power, non-null nPlan, or \n")
  }

  ## Fill and name ----

  result <- utils::modifyList(result, tempResult)

  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["sigma"]] <- sigma
  result[["kappa"]] <- kappa
  result[["testType"]] <- testType
  result[["ratio"]] <- ratio
  result[["eType"]] <- eType
  result[["designScenario"]] <- designScenario

  ## Name esMin ----
  esMin <- result[["esMin"]]

  if (!is.null(meanDiffTrue))
    result[["esTrue"]] <- meanDiffTrue

  if (is.na(esMin))
    esMin <- NULL

  if (!is.null(esMin))
    names(esMin) <- "mean difference"

  result[["esMin"]] <- esMin

  ## Name nPlan ----
  nPlan <- result[["nPlan"]]

  if (!is.null(nPlan)) {
    if (designScenario %in% c(2:3, "3b")) {
      n2Plan <- nPlan[2]

      names(nPlan) <- if (is.na(n2Plan)) "n1Plan" else c("n1Plan", "n2Plan")
    }
  }

  result[["nPlan"]] <- nPlan

  ## Name parameter ----
  parameter <- result[["parameter"]]

  if (!is.null(parameter) && is.null(names(parameter))) {
    names(parameter) <- switch(eType,
                               "mom"="gMom",
                               "eGauss"="g",
                               "imom"="tau",
                               "eCauchy"="kappaG",
                               "grow"="phiS")
  }

  result[["parameter"]] <- parameter

  ## Name h0 -----
  names(h0) <- "mu"
  result[["h0"]] <- h0
  result[["relevanceTest"]] <- relevanceTest

  result[["call"]] <- sys.call()

  result <- Filter(Negate(is.null), result)
  class(result) <- "saviDesign"

  return(result)
}

#' Helper function to designing a savi Z-test (output nPlan)
#'
#' Finds the parameter and power when provided with only alpha, esMin, and nPlan
#'
#' @inheritParams designSaviZ
#'
#' @return A list with the parameter and the targeted nPlan amongst other items
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#' @export
#'
#' @examples
#' designSaviZ1aWantNPlan(meanDiffMin=0.9, power=0.7, nSim=10)
designSaviZ1aWantNPlan <- function(
    meanDiffMin, power, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma, beta=NULL,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, meanDiffTrue=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    relevanceTest=FALSE, relevanceSize=NULL, alphaRelevance=NULL,
    highN=1e4L, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  if (is.null(meanDiffTrue)) meanDiffTrue <- meanDiffMin

  power <- matchPowerWith("power"=power, "beta"=beta)

  samplingResult <- computeNPlanSaviZ(
    "meanDiffTrue"=meanDiffTrue, "power"=power, "alpha"=alpha,
    "alternative"=alternative, "sigma"=sigma, "kappa"=kappa, "ratio"=ratio,
    "parameter"=parameter, "testType"=testType, "eType"=eType,
    "wantSamplePaths"=wantSamplePaths,
    "meanDiffMin"=meanDiffMin, "wantSimData"=wantSimData,
    "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "highN"=highN, "alphaRelevance"=alphaRelevance)

  result <- designSavi1aHelper("samplingResult"=samplingResult,
                               "esMin"=meanDiffMin, "power"=power,
                               "beta"=NULL, "ratio"=ratio, "testType"=testType)
  return(result)
}

#' Helper function to designing a Z-test (output power)
#'
#' Finds the parameter and power when provided with only alpha, esMin, and nPlan
#'
#' @inheritParams designSaviZ
#'
#' @return A list with the parameter and power amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviZ2WantPower(meanDiffTrue=0.9, nPlan=7, nSim=10)
designSaviZ2WantPower <- function(
    meanDiffTrue, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    meanDiffMin=meanDiffTrue,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    relevanceTest=FALSE, relevanceSize=NULL, alphaRelevance=NULL,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1
  nPlan <- checkAndReturnNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  samplingResult <- computePowerSaviZ(
    "meanDiffTrue"=meanDiffTrue, "nPlan"=nPlan, "alpha"=alpha,
    "sigma"=sigma, "kappa"=kappa, "alternative"=alternative,
    "testType"=testType, "parameter"=parameter,
    "meanDiffMin"=meanDiffMin, "seed"=seed,
    "eType"=eType, "wantSamplePaths"=wantSamplePaths,
    "wantSimData"=wantSimData,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "alphaRelevance"=alphaRelevance, "nSim"=nSim, "nBoot"=nBoot, "pb"=pb)

  result <- designSavi2Helper("samplingResult"=samplingResult,
                              "esMin"=meanDiffMin, "nPlan"=nPlan, "ratio"=ratio,
                              "testType"=testType)
  return(result)
}

#' Helper function to designing a Z-test (output esMin)
#'
#' Finds the parameter and esMin when provided with only alpha, power, and nPlan
#'
#' @inheritParams designSaviZ
#'
#' @return A list with the parameter and the targeted esMin amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviZ3WantEsMin(power=0.7, nPlan=10)
designSaviZ3WantEsMin <- function(
    power, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL, beta=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1
  nPlan <- checkAndReturnNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  result <- list("parameter"=NULL, "esMin"=NULL,
                 "nPlan"=nPlan, "power"=power, "ratio"=ratio,
                 "note"=NULL)

  meanDiffMin <- tryOrFailWithNA(
    computeMinEsBatchSaviZ(
      "nPlan"=nPlan, "alpha"=alpha, "power"=power, "beta"=NULL, "sigma"=sigma,
      "kappa"=kappa, "alternative"=alternative, "testType"=testType,
      "parameter"=parameter, "eType"=eType)
  )

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="z",
    "esMin"=meanDiffMin, "sigma"=sigma,
    "alternative"=alternative, "eType"=eType
  )

  result[["parameter"]] <- parameter
  result[["esMin"]] <- meanDiffMin

  if (is.na(meanDiffMin))
    result[["note"]] <- "No meanDiffMin found for the provided beta and nPlan"
  else
    result[["note"]] <- "The reported meanDiffMin is based on the batch analysis."

  return(result)
}

#' Helper function to designing a Z-test (output meanDiffMin based on the shortest interval at nPlan)
#'
#' Finds the parameter and meanDiffMin when provided with only alpha, nPlan
#'
#' @inheritParams designSaviZ
#'
#' @return A list with the parameter and the parameter amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviZ3bWantParameter(nPlan=13)
designSaviZ3bWantParameter <- function(
    nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"), ...) {
  # TODO(Alexander): Two-sample and imom don't play well

  defaultErrorText <-
    "Can't compute a design based on alpha and the provided planned sample size(s) alone."

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1
  nPlan <- checkAndReturnNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  n1 <- nPlan[1]
  n2 <- nPlan[2]

  paired <- if (testType=="paired") TRUE else FALSE

  nEff <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1 else (1/n1+1/n2)^(-1)

  if (nEff <= 0)
    stop(defaultErrorText)

  result <- list("parameter"=NULL, "esMin"=NULL,
                 "nPlan"=nPlan, "ratio"=ratio,
                 "note"=NULL)

  gCandidate <- (-1-lamW::lambertWm1(-alpha^2/exp(1)))/nEff

  meanDiffMinCandidate <- sqrt(gCandidate)*sigma

  if (eType=="eGauss") {
    parameter <- gCandidate
    meanDiffMin <- meanDiffMinCandidate
  } else {

    upperMeanDiffMin <- if (eType %in% c("mom", "grow")) 2*meanDiffMinCandidate else max(2*meanDiffMinCandidate, 0.02)

    meanDiffDomain <- seq(meanDiffMinCandidate/4, upperMeanDiffMin, length.out=1e3)
    ciWidths <- rep(Inf, length(meanDiffDomain))

    parameterDomain <- if (eType=="mom") (meanDiffDomain/sigma)^2/2 else meanDiffDomain/sigma

    for (i in seq_along(ciWidths)) {
      # the parameter is already standardised, so sigma=1 should do
      tempRes <- computeConfidenceIntervalZ(nEff=nEff, meanObs=0, sigma=1,
                                            parameter=parameterDomain[i],
                                            eType=eType, ciValue=1-alpha,
                                            alternative="twoSided")
      ciWidths[i] <- tempRes[2]
    }

    if (sum(is.infinite(ciWidths))==length(ciWidths))
      stop(defaultErrorText)

    minIndex <- which.min(ciWidths)

    meanDiffMin <- meanDiffDomain[minIndex]
    parameter <- parameterDomain[minIndex]

    if (minIndex!=1 && is.infinite(ciWidths[minIndex-1])) {
      result[["note"]] <- "Unstable design based on alpha and nPlan alone."
      warning("Unstable: The parameter corresponds to the smallest parameter value",
              "for which the ci width can be calculated. Another eType might yield more stable designs.")
    }
  }

  if (eType=="grow" && alternative=="less") {
    parameter <- -parameter
    meanDiffMin <- -meanDiffMin
  }

  result[["parameter"]] <- parameter
  result[["esMin"]] <- meanDiffMin

  return(result)
}

# Batch design fnts ------

#' Helper function: Computes nPlan based on meanDiffTrue, alpha and power.
#'
#' @inheritParams designSaviZ
#' @inheritParams sampleStoppingTimesSaviZ
#'
#' @return a list which contains at least nPlan and the savi test defining parameter.
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
computeNPlanBatchSaviZ <- function(
    meanDiffTrue, alpha=0.05, power=0.8, sigma=1, kappa=sigma,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    parameter=NULL, meanDiffMin=NULL, beta=NULL,
    highN=1e4L, ratio=1, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  result <- list(nPlan=NULL, parameter=NULL)

  n1Plan <- NULL
  n2Plan <- NULL

  n1OverNEffRatio <- if (testType=="twoSample") (1+ratio)/ratio else 1

  if (is.null(meanDiffMin))
    meanDiffMin <- meanDiffTrue

  meanDiffMin <- suppressWarnings(
    checkAndReturnEsMinParameterSide(
      "paramToCheck"=meanDiffMin, "alternative"=alternative,
      "esMinName"="meanDiffMin")
  )

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="z",
    "esMin"=meanDiffMin, "sigma"=sigma,
    "alternative"=alternative, "eType"=eType
  )

  meanDiffTrue <- abs(meanDiffTrue)

  # Note(Alexander): Computes one-sided grow sample size used
  # to bound the candidate set of nEff
  qB <- stats::qnorm(1-power)

  nTemp <- exp(2*(log(kappa)-log(meanDiffTrue))) *
    (2*qB^2 - 2*qB*sqrt(qB^2+2*sigma^2/kappa^2*log(1/alpha))
     +2*kappa^2/sigma^2*log(1/alpha))

  tempAlternative <- switch(alternative,
                            "twoSided"="twoSided",
                            "greater"="greater",
                            "less"="greater")

  if (eType=="grow" && alternative %in% c("greater", "less") && parameter==abs(meanDiffTrue)) {
    nEff <- nTemp
  } else {
    if (testType=="twoSample") {
      n1Func <- function(nEff) (1+ratio)/ratio*nEff
      n2Func <- function(nEff) (1+ratio)*nEff
    } else if (testType %in% c("oneSample", "paired")) {
      n1Func <- function(nEff) nEff
      n2Func <- function(nEff) NULL
    }

    targetFunction <- function(nEff) {
      saviZTestStat(stats::qnorm("p"=1-power, "mean"=sqrt(nEff)*meanDiffTrue/sigma, "sd"=kappa/sigma),
                    "n1"=n1Func(nEff), "n2"=n2Func(nEff),
                    "parameter"=parameter, "sigma"=sigma,
                    "alternative"=tempAlternative, "eType"=eType)$eValue-1/alpha
    }

    tempResult <- try(stats::uniroot(targetFunction, interval=c(nTemp/2, 3*nTemp)))

    if (isTryError(tempResult))
      tempResult <- try(stats::uniroot(targetFunction, interval=c(1, highN)))

    if (isTryError(tempResult)) {
      warning("Can't compute the batched planned sample size.",
              "Set sampling horizon to highN = ", highN)
      tempResult <- list("root"=highN)
    }

    nEff <- tempResult[["root"]]
  }

  # Two-sample/paired stuff
  if (testType == "twoSample") {
    n1Plan <- ceil(nEff * n1OverNEffRatio)
    n2Plan <- ceil(nEff * n1OverNEffRatio * ratio)
  } else {
    n1Plan <- ceil(nEff)
    n2Plan <- if (testType == "paired") n1Plan else NULL
  }

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  names(parameter) <- switch(eType,
                             "mom"="gMom",
                             "eGauss"="g",
                             "imom"="tau",
                             "eCauchy"="kappaG",
                             "grow"="phiS")

  if (eType=="grow" && alternative=="less" && parameter > 0)
    parameter <- -parameter

  result[["parameter"]] <- parameter

  return(result)
}


#' Helper function: Computes 1-power based on meanDiffTrue and nPlan
#'
#' @inheritParams designSaviZ
#' @inheritParams sampleStoppingTimesSaviZ
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @return numeric that represents the type II error
computeBetaBatchSaviZ <- function(
    meanDiffTrue, nPlan, alpha=0.05, sigma=1, kappa=sigma,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow")) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  meanDiffTrue <- checkAndReturnEsMinParameterSide(
    "paramToCheck"=meanDiffTrue, "alternative"=alternative,
    "esMinName"="meanDiffMin")

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="z",
    "esMin"=meanDiffTrue, "sigma"=sigma,
    "alternative"=alternative, "eType"=eType
  )

  if (eType=="grow") {
    phiS <- parameter
    if (alternative=="twoSided") {
      upperBoundOfIntegral <- exp(4*log(sigma)-log(nEff)-2*log(kappa)-2*log(phiS)) *
        (acosh(exp(nEff*phiS^2/(2*sigma^2))/alpha))^2

      result <- stats::pchisq(upperBoundOfIntegral, df=1, ncp=nEff*meanDiffTrue^2/kappa^2, lower.tail = TRUE)
    } else {
      upperBoundOfIntegral <- sqrt(nEff)*(phiS-2*meanDiffTrue)/(2*kappa) -
        sigma^2*log(alpha)/(kappa*sqrt(nEff))*1/phiS

      result <- stats::pnorm(upperBoundOfIntegral, lower.tail = TRUE)
    }
  } else if (eType=="eGauss") {
    g <- parameter
    if (alternative=="twoSided") {
      upperBoundOfIntegral <- (1+nEff*g)/(nEff*g)*(log(1+nEff*g)-2*log(alpha))

      result <- stats::pchisq(upperBoundOfIntegral, df=1,
                              ncp=nEff*meanDiffTrue^2/kappa^2,
                              lower.tail = TRUE)
    } else {
      stop("Not yet implemented")
    }
  }

  return(result)
}

#' Computes the meanDiffMin that is detectable with probability "power" for given nPlan
#'
#' @inheritParams designSaviZ
#'
#' @return numeric > 0 that represents the minimal detectable mean difference
#' @export
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @examples
#' computeMinEsBatchSaviZ(27)
computeMinEsBatchSaviZ <- function(
    nPlan, alpha=0.05, power=0.8, sigma=1, kappa=sigma,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL, beta=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    lowEsTrue=0.01, highEsTrue=0.002, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  stopifnot(nPlan>0)
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (eType=="eCauchy" && nEff > 1e5 || eType=="imom" && nEff > 1e8 ||
      eType=="mom" && alternative %in% c("greater", "less") && nEff > 1e8) {
    stop('Unable to compute meanDiffMin for eType="', eType, '" with such high sample size(s).')
  }

  if (nEff > 1e9)
    warning("The computed meanDiffMin and parameter might not be reliable")

  if (eType=="mom") {
    paramFunc <- function(deltaTrue) abs(deltaTrue)^2/2
  } else if (eType=="eGauss") {
    paramFunc <- function(deltaTrue) deltaTrue^2
  } else if (eType=="imom") {
    paramFunc <- function(deltaTrue) abs(deltaTrue)
  } else if (eType=="eCauchy") {
    paramFunc <- function(deltaTrue) abs(deltaTrue)
  } else if (eType=="grow") {
    paramFunc <- function(deltaTrue) deltaTrue
  }

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1

  tempAlternative <- switch(alternative,
                            "twoSided"="twoSided",
                            "greater"="greater",
                            "less"="greater")

  if (testType=="twoSample") {
    n1Func <- function(nEff) (1+ratio)/ratio*nEff
    n2Func <- function(nEff) (1+ratio)*nEff
  } else if (testType %in% c("oneSample", "paired")) {
    n1Func <- function(nEff) nEff
    n2Func <- function(nEff) NULL
  }

  targetFunction <- function(deltaTrue) {
    saviZTestStat(stats::qnorm("p"=1-power, "mean"=sqrt(nEff)*deltaTrue, "sd"=kappa/sigma),
                  "n1"=n1Func(nEff), "n2"=n2Func(nEff),
                  "parameter"=paramFunc(deltaTrue), "sigma"=1,
                  "alternative"=tempAlternative, "eType"=eType)$eValue-1/alpha
  }

  unirootBounds <- setLowAndHighEsTrueZ(nEff=nEff, eType=eType, alternative=alternative)

  tempResult <- try(stats::uniroot(targetFunction,
                                   interval=c(unirootBounds[["tempLowEsTrue"]],
                                              unirootBounds[["tempHighEsTrue"]])))

  if (isTryError(tempResult))
    tempResult <- stats::uniroot(targetFunction, interval=c(lowEsTrue, highEsTrue))

  result <- tempResult[["root"]]*sigma

  if (alternative=="less")
    result <- -result

  return(result)
}

#' Helper function to determine the lower and upper bound for the search space for standardised meanDiffMin
#'
#' @inheritParams computeConfidenceIntervalZ
#' @inheritParams designSaviZ
#'
#' @return a list of low and high values for the parameter
#' @export
#'
#' @examples
#' setLowAndHighEsTrueZ(10)
setLowAndHighEsTrueZ <- function(nEff, eType="mom", alternative="twoSided",
                                 lowEsTrue=1e-7, highEsTrue=1e-3) {

  if (nEff >= 1 && nEff <= 1e1) {
    tempLowEsTrue <- 1
    tempHighEsTrue <- 4
  } else if (nEff >= 1e1+1 && nEff <= 1e2) {
    tempLowEsTrue <- 0.3
    tempHighEsTrue <- 1.5
  } else if (nEff >= 1e2+1 && nEff <= 1e3) {
    tempLowEsTrue <- 0.1
    tempHighEsTrue <- 1
  } else if (nEff >= 1e3+1 && nEff <= 1e4) {
    tempLowEsTrue <- 0.03
    tempHighEsTrue <- 0.2
  } else if (nEff >= 1e4+1 && nEff <= 1e5) {
    tempLowEsTrue <- 0.01
    tempHighEsTrue <- 0.05
  } else if (nEff >= 1e5+1 && nEff <= 1e6) {
    tempLowEsTrue <- 0.003
    tempHighEsTrue <- 0.02
  } else if (nEff >= 1e6+1 && nEff <= 1e7) {
    tempLowEsTrue <- 0.001
    tempHighEsTrue <- 0.009
  } else if (nEff >= 1e7+1 && nEff <= 1e8) {
    tempLowEsTrue <- 0.0002
    tempHighEsTrue <- 0.007

    if (eType=="grow")
      tempHighEsTrue <- 0.002

  } else if (nEff >= 1e8+1 && nEff <= 1e9) {
    tempLowEsTrue <- 0.00005
    tempHighEsTrue <- 0.001
  } else {
    tempLowEsTrue <- lowEsTrue
    tempHighEsTrue <- highEsTrue
  }

  result <- list(tempLowEsTrue=tempLowEsTrue, tempHighEsTrue=tempHighEsTrue)
  return(result)
}

# Sampling functions for design ----

#' Simulate stopping times for the savi Z-test
#'
#' @inheritParams designSaviZ
#' @param nMax integer > 0, maximum sample size of the (first) sample in each sample path.
#' @param wantEValuesAtNMax logical. If \code{TRUE}, then compute eValues at nMax. Default \code{FALSE}.
#' @param wantSamplePaths logical. If \code{TRUE}, then output the (stopped) sample paths. Default \code{TRUE}.
#' @param wantSimData logical. If \code{TRUE}, then output the simulated data.
#'
#' @return a list with stoppingTimes and breakVector. Entries of breakVector are 0, 1. A 1 represents stopping
#' due to exceeding nMax, and 0 due to 1/alpha threshold crossing, which implies that in corresponding stopping
#' time is Inf.
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' sampleStoppingTimesSaviZ(0.7, nSim=10, nMax=20)
sampleStoppingTimesSaviZ <- function(
    meanDiffTrue, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    sigma=1, kappa=sigma,
    ratio=1, meanDiffMin=NULL, parameter=NULL, nMax=1e3L,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantEValuesAtNMax=FALSE, power=NULL,
    wantSamplePaths=TRUE, wantSimData=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, relevanceTest=FALSE,
    relevanceSize=NULL, beta=NULL, alphaRelevance=NULL, ...) {

  stopifnot(alpha > 0, alpha <= 1,
            is.finite(nMax),
            is.finite(meanDiffTrue))

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (relevanceTest) {
    relevanceSize <- matchRelevanceParameterWith(
      "relevanceSize"=relevanceSize, "esMin"=meanDiffMin,
      "esTrue"=meanDiffTrue
    )

    alphaRelevance <- matchAlphaRelevanceWith(
      "alphaRelevance"=alphaRelevance, "alpha"=alpha, "beta"=1-power
    )

    stopifnot(alphaRelevance > 0, alphaRelevance < 1)
  }

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  result <- constructSampleStoppingTimesList(
    "nSim"=nSim, "nMax"=nMax,
    "wantEValuesAtNMax"=wantEValuesAtNMax,
    "wantSamplePaths"=wantSamplePaths)

  meanDiffTrue <- abs(meanDiffTrue)

  meanDiffMin <- if (is.null(meanDiffMin)) abs(meanDiffTrue) else abs(meanDiffMin)

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="z",
    "esMin"=meanDiffMin, "sigma"=sigma,
    "alternative"=alternative, "eType"=eType
  )

  relevanceTestSim <- NULL

  if (relevanceTest) {
    relevanceParameter <- matchRelevanceParameterWith(
      "relevanceSize"=relevanceSize,
      "esMin"=meanDiffMin, "esTrue"=meanDiffTrue)

    relevanceTestSim <-
      list("eValuesStopped"=Matrix::sparseVector(x=0, i=1, length=nSim),
           "samplePaths"=result[["samplePaths"]],
           "stoppingTimes"=Matrix::sparseVector(x=0, i=1, length=nSim),
           "parameter"=relevanceParameter, "alpha"=alphaRelevance)
  }

  if (testType=="twoSample" && length(nMax)==1) {
    nMax <- c(nMax, ceil(ratio*nMax))
    n1Max <- nMax[1]
    n2Max <- nMax[2]
    ratio <- nMax[2]/nMax[1]
  } else if (testType %in% c("paired", "oneSample")){
    n1Max <- nMax[1]
    n2Max <- NULL
    nMax <- n1Max
    ratio <- 1
  }

  if (pb)
    pbSavi <- utils::txtProgressBar(style=3, title="Savi test threshold crossing")

  tempN <- defineTTestN("lowN"=1, "highN"=nMax[1], "ratio"=ratio, "testType"=testType)

  n1Vector <- tempN[["n1"]]
  n2Vector <- tempN[["n2"]]
  nEffVector <- tempN[["nEff"]]

  simData <- generateNormalData("nPlan"=nMax, "nSim"=nSim, "deltaTrue"=meanDiffTrue/kappa,
                                "sigma"=kappa, "paired"=FALSE, "seed"=seed)

  for (sim in seq_along(result[["stoppingTimes"]])) {
    if (testType %in% c("oneSample", "paired")) {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/(n1Vector)*cumsum(x1)
      zVector <- sqrt(n1Vector)*x1BarVector/sigma
    } else {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/(n1Vector)*cumsum(x1)
      x1BarVector <- x1BarVector[n1Vector]

      x2 <- simData[["dataGroup2"]][sim, ]
      x2Cumsum <- cumsum(x2)[n2Vector]
      x2BarVector <- 1/n2Vector*x2Cumsum

      zVector <- sqrt(nEffVector)*(x1BarVector - x2BarVector)/sigma
    }

    if (wantEValuesAtNMax) {
      tempResult <- saviZTestStat("z"=zVector[length(zVector)],
                                  "parameter"=parameter,
                                  "n1"=nMax[1], n2=nMax[2],
                                  "alternative"=alternative, "sigma"=sigma,
                                  "eType"=eType)
      result[["eValuesAtNMax"]][sim] <- tempResult[["eValue"]]
    }

    for (j in seq_along(n1Vector)) {
      tempResult <- saviZTestStat("z"=zVector[j], "parameter"=parameter,
                                  "n1"=n1Vector[j], "n2"=n2Vector[j],
                                  "alternative"=alternative, "sigma"=sigma,
                                  "eType"=eType)

      evidenceNow <- tempResult[["eValue"]]

      if (wantSamplePaths)
        result[["samplePaths"]][sim, j] <- evidenceNow

      if (evidenceNow >= 1/alpha) {
        result[["stoppingTimes"]][sim] <- n1Vector[j]
        result[["eValuesStopped"]][sim] <- evidenceNow

        break()
      }

      # TODO(Alexander): Here use reciprocal as evidence for the null
      #
      # if (relevanceTest && !isTRUE(growFutility) && evidenceNow < beta) {
      #   result[["stoppingTimes"]][sim] <- n1Vector[j]
      #   result[["eValuesStopped"]][sim] <- evidenceNow
      #   result[["stoppedVector"]][sim] <- -1
      #
      #   if (wantSamplePaths) {
      #     result[["samplePaths"]][sim, j:nMax[1]] <- evidenceNow
      #   }
      #
      #   break()
      # }

      if (relevanceTest) {
        relevanceRes <- saviRelevanceZStat(
          "z"=zVector[j], "parameter"=relevanceParameter,
          "n1"=n1Vector[j], "n2"=n2Vector[j],
          "alternative"=alternative,
          "sigma"=sigma, "eType"="grow")

        relevanceEValue <- relevanceRes[["eValue"]]

        if (wantSamplePaths)
          relevanceTestSim[["samplePaths"]][sim, j] <- relevanceEValue

        if (relevanceEValue < alphaRelevance) {
          result[["breakVector"]][sim] <- -1
          result[["stoppingTimes"]][sim] <- nMax

          relevanceTestSim[["stoppingTimes"]][sim] <- n1Vector[j]
          relevanceTestSim[["eValuesStopped"]][sim] <- relevanceEValue

          break()
        }
      }

      # Note(Alexander): If passed maximum nPlan[1] stop.
      #   For power calculations if beyond nPlan[1], then set to Inf, doesn't matter for the quantile
      #
      if (n1Vector[j] >= nMax[1]) {
        result[["stoppingTimes"]][sim] <- n1Vector[j]
        result[["breakVector"]][sim] <- 1
        result[["eValuesStopped"]][sim] <- evidenceNow
        break()
      }
    }

    if (pb)
      utils::setTxtProgressBar(pbSavi, "value"=sim/nSim, "title"="Trials")
  }

  if (pb)
    close(pbSavi)

  result[["parameter"]] <- parameter
  result[["n1Vector"]] <- n1Vector
  result[["ratio"]] <- ratio

  if (wantSimData)
    result[["simData"]] <- simData

  if (wantSamplePaths) {
    samplePaths <- result[["samplePaths"]]
    samplePaths[is.na(samplePaths)] <- 0
    result[["samplePaths"]] <- Matrix::Matrix(samplePaths, sparse=TRUE)

    if (relevanceTest) {
      relevanceSamplePaths <- relevanceTestSim[["samplePaths"]]
      relevanceSamplePaths[is.na(relevanceSamplePaths)] <- 0
      relevanceTestSim[["samplePaths"]] <- Matrix::Matrix(relevanceSamplePaths, sparse=TRUE)
    }
  }

  result[["relevanceTestSim"]] <- relevanceTestSim

  return(result)
}

#' Helper function: Computes the power based on meanDiffTrue and nPlan
#'
#' @inheritParams designSaviZ
#' @inheritParams sampleStoppingTimesSaviZ
#'
#' @return a list which contains at least power and an adapted bootObject of class \code{\link[boot]{boot}}.
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' computePowerSaviZ(meanDiffTrue=0.7, 20, nSim=10)
computePowerSaviZ <- function(
    meanDiffTrue, nPlan, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    meanDiffMin=meanDiffTrue,
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    relevanceTest=FALSE, relevanceSize=NULL, alphaRelevance=NULL,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan) == 2) nPlan[2]/nPlan[1] else 1

  if (testType=="twoSample" && length(nPlan)==1) {
    nPlan <- c(nPlan, nPlan)
    warning('testType=="twoSample" specified, but nPlan[2] not provided. nPlan[2] is set to ratio = ', ratio,
            'times nPlan[1] = ', nPlan[2])
  }

  # TODO(Alexander): Revise logic on meanDiffTrue and meanDiffMin
  #
  meanDiffMin <- checkAndReturnEsMinParameterSide(
    "paramToCheck"=meanDiffMin, "alternative"=alternative,
    "esMinName"="meanDiffMin")

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="z",
    "esMin"=meanDiffMin, "sigma"=sigma,
    "alternative"=alternative, "eType"=eType
  )

  samplingResult <- sampleStoppingTimesSaviZ(
    "meanDiffTrue"=meanDiffTrue, "alpha"=alpha,
    "alternative"=alternative, "testType"=testType,
    "sigma"=sigma, "kappa"=kappa,
    "ratio"=ratio, "parameter"=parameter, "nMax"=nPlan,
    "eType"=eType, "wantSimData"=wantSimData,
    "wantEValuesAtNMax"=TRUE, "wantSamplePaths"=wantSamplePaths,
    "pb"=pb, "seed"=seed, "nSim"=nSim,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "alphaRelevance"=alphaRelevance, ...)

  result <- computePowerBootstrapper("samplingResult"=samplingResult,
                                    "parameter"=parameter, "nPlan"=nPlan,
                                    "nBoot"=nBoot)
  return(result)
}


#' Helper function: Computes nPlan based on meanDiffTrue, alpha and power
#'
#'
#' @inheritParams designSaviZ
#' @inheritParams sampleStoppingTimesSaviZ
#'
#' @return a list which contains at least nPlan and an adapted bootObject of class  \code{\link[boot]{boot}}.
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' computeNPlanSaviZ(0.7, 0.2, nSim=10)
computeNPlanSaviZ <- function(
    meanDiffTrue, power=0.8, alpha=0.05,
    alternative=c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    sigma=1, kappa=sigma,
    meanDiffMin=NULL, beta=NULL,
    ratio=1, parameter=NULL, nMax=1e8,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    relevanceTest=FALSE, relevanceSize=NULL, alphaRelevance=NULL,
    highN=1e4L, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  if (is.null(meanDiffMin)) {
    meanDiffTrue <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=meanDiffTrue, "alternative"=alternative,
      "esMinName"="meanDiffTrue")

    meanDiffMin <- meanDiffTrue
  } else {
    meanDiffMin <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=meanDiffMin, "alternative"=alternative,
      "esMinName"="meanDiffMin")
  }

  power <- matchPowerWith("power"=power, "beta"=beta)

  # TODO(Alexander)
  #
  tempObj <- computeNPlanBatchSaviZ(
    "meanDiffTrue"=meanDiffMin, "alpha"=alpha,
    "power"=power, "sigma"=sigma, "kappa"=kappa,
    "alternative"=alternative, "testType"=testType,
    "parameter"=parameter, "ratio"=ratio, "eType"=eType,
    "meanDiffMin"=meanDiffMin, "highN"=highN,
    "wantSimData"=wantSimData, ...)

  nPlanBatch <- tempObj[["nPlan"]]

  if (is.null(parameter))
    parameter <- tempObj[["parameter"]]

  samplingResult <- sampleStoppingTimesSaviZ(
    "meanDiffTrue"=meanDiffTrue, "alpha"=alpha, "power"=power,
    "beta"=NULL, "alternative"=alternative, "testType"=testType,
    "sigma"=sigma, "kappa"=kappa,
    "ratio"=ratio, "parameter"=parameter, "nMax"=nPlanBatch,
    "eType"=eType, "meanDiffMin"=meanDiffMin,
    "wantSamplePaths"=wantSamplePaths, "wantSimData"=wantSimData,
    "pb"=pb, "seed"=seed, "nSim"=nSim,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "alphaRelevance"=alphaRelevance, ...)

  result <- computeNPlanBootstrapper(
    "samplingResult"=samplingResult, "parameter"=parameter,
    "beta"=NULL, "power"=power, "nPlanBatch"=nPlanBatch, "nBoot"=nBoot
  )
  return(result)
}

# Helpers ------


#' Help function to compute the effective sample size based on a length 2 vector of samples
#'
#' @inheritParams designSaviZ
#' @param n vector of length at most 2 representing the sample sizes of the first and second group
#' @param silent logical, if true, then turn off warnings
#'
#' @return a numeric that represents the effective sample size.
computeNEff <- function(n, testType=c("oneSample", "paired", "twoSample"), silent=TRUE) {
  testType <- match.arg(testType)

  stopifnot(all(n > 0))

  if (testType=="oneSample") {
    if (length(n) > 1)
      if (!silent)
        warning("One sample test wanted, but n is longer. Only first element used")

    nEff <- n[1]
  } else if (testType=="paired") {
    nEff <- n[1]

    if (!silent) {
      if (length(n)==1)
        warning("Paired sample design wanted, but n is of length 1. Copied n")

      if (length(n) > 1) {
        n2 <- n[2]
        if (nEff != n2)
          warning("Paired sample design wanted, but n[1] != n[2]. Copied n[2] ignored and n[1] copied")
      }
    }
  } else if (testType=="twoSample") {
    if (length(n) == 1)
      stop("Two sample design specified, but only one sample size provided")

    n1 <- n[1]
    n2 <- n[2]
    nEff <- (1/n1+1/n2)^(-1)
  }
  return(nEff)
}



# Workshop functions -----

#' @rdname pValueZTest
#' @inheritParams saviZTestStat
#' @inheritParams designSaviZ
#'
#' @export
pValueFromZStat <- function(z,
                            alternative=c("twoSided", "less", "greater"),
                            ...) {

  alternative <- match.arg(alternative)

  pValue <- 1-stats::pnorm(abs(z), mean=0, sd=1)

  if (alternative=="twoSided") { # two-sided
    pValue <- 2*pValue
  } else if (alternative=="less") { # one-sided
    if (z >= 0)
      pValue <- 1-pValue
  } else if (alternative=="greater") {
    if (z <= 0)
      pValue <- 1-pValue
  }

  return(unname(pValue))
}

#' Computes the p-value for the Z-test
#'
#'
#' @aliases pValueZTest
#' @inheritParams saviZTest
#'
#' @return pValueTest object
#' @export
#'
#' @examples
#' pValueZTest(rnorm(10))
pValueZTest <- function(x, y=NULL, paired=FALSE, ciValue=NULL,
                        alpha=0.05, sigma=1, h0=0,
                        alternative=c("twoSided", "greater", "less"), ...) {

  result <- list("statistic"=NULL, "n"=NULL, "eValue"=NULL, "confSeq"=NULL, "estimate"=NULL,
                 "testType"=NULL, "dataName"=NULL, "h0"=NULL, "sigma"=NULL, "call"=sys.call())
  class(result) <- "pValueTest"

  alternative <- match.arg(alternative)

  if (is.null(y)) {
    testType <- "oneSample"
    n <- nEff <- n1 <- length(x)
    n2 <- NULL

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    meanObs <- estimate <- mean(x)

    names(estimate) <- "mean of x"
    names(n) <- "n1"
  } else {
    nEff <- n1 <- length(x)
    n2 <- length(y)

    if (paired) {
      if (n1 != n2)
        stop("Data error: Error in complete.cases(x, y): Paired analysis requested, ",
             "but the two samples are not of the same size.")

      testType <- "paired"

      meanObs <- estimate <- mean(x-y)
      names(estimate) <- "mean of the differences"
    } else {
      testType <- "twoSample"

      nEff <- (1/n1+1/n2)^(-1)
      estimate <- c(mean(x), mean(y))
      names(estimate) <- c("mean of x", "mean of y")
      meanObs <- estimate[1]-estimate[2]
    }

    n <- c(n1, n2)
    names(n) <- c("n1", "n2")
  }

  if (is.null(ciValue))
    ciValue <- 1-alpha

  if (ciValue < 0 || ciValue > 1)
    stop("Can't make a confidence sequence with ciValue < 0 or ciValue > 1, or alpha < 0 or alpha > 1")

  zStat <- tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sigma)

  if (is.na(zStat))
    stop("Could not compute the z-statistic")

  names(zStat) <- "z"

  pValue <- pValueFromZStat(zStat, "alternative"=alternative)

  argumentNames <- getArgs()
  xLabel <- extractNameFromArgs(argumentNames, "x")

  if (is.null(y)) {
    dataName <- xLabel
  } else {
    yLabel <- extractNameFromArgs(argumentNames, "y")
    dataName <- paste(xLabel, "and", yLabel)
  }

  result[["testType"]] <- testType
  result[["statistic"]] <- zStat
  result[["estimate"]] <- estimate
  result[["dataName"]] <- dataName
  result[["ciValue"]] <- ciValue
  result[["n"]] <- n

  result[["confInt"]] <- computeConfidenceIntervalZ("nEff"=nEff, "meanObs"=meanObs,
                                                    "sigma"=sigma, "ciValue"=ciValue,
                                                    "alternative"=alternative, "eType"="freq",
                                                    "parameter"=NULL)
  result[["pValue"]] <- pValue

  names(result[["statistic"]]) <- "z"

  return(result)
}

#' Logarithmic marginal likelihood of the normal with conjugate priors
#'
#' @inheritParams subjectiveBfZStat
#' @param n integer, sample size
#' @param x numeric, sample mean
#' @param s2 numeric > 0, sample variance
#' @param a numeric, prior mean of the population mean mu
#' @param b numeric > 0, prior standard deviation of the population mean mu
#'
#'
#' @return numeric, the logarithm of the marginal likelihood
#' @export
#'
#' @examples
#' zMarg(12, 0.3, 2)
zMarg <- function(n, x, s2, a=0, b=1, sigma=1) {
  if (n > 1) {
    result <- -n/2*log(2*pi)+(1-n)*log(sigma)-1/2*log(n*b^2+sigma^2) -
      (n-1)*s2/(2*sigma^2)-n*(x-a)^2/(2*(n*b^2+sigma^2))
  } else {
    result <- -n/2*log(2*pi)+(1-n)*log(sigma)-1/2*log(n*b^2+sigma^2) -
      n*(x-a)^2/(2*(n*b^2+sigma^2))
  }

  return(result)
}

#' A subjective Bayes factor for the two-sample Z-test
#'
#' Based on conjugate priors with a total of 6 hyperparameters.
#'
#' @param x1 numeric, sample mean of group 1
#' @param s21 numeric, sample variance of group 1
#' @param n1 integer sample size of group 1
#' @param x2 numeric, sample mean of group 2
#' @param s22 numeric, sample variance of group 2
#' @param n2 integer sample size of group 2
#' @param sigma numeric > 0, population standard deviation
#' @param a1 numeric, prior mean of the population mean mu1 of group 1
#' @param b1 numeric > 0, prior standard deviation of the population mean mu1 of group 1
#' @param a2 numeric, prior mean of the population mean mu2 of group 2
#' @param b2 numeric > 0, prior standard deviation of the population mean mu2 of group 2
#' @param a0 numeric, prior mean of the overall population mean mu0 of both groups
#' @param b0 numeric > 0, prior standard deviation of the population mean mu0 of both groups
#' @param log logical, default FALSE, if TRUE then return logarithm of the subjective Bayes factor outcome
#'
#' @return numeric > 0 representing the subjective Bayes factor outcome in favour of the alternative over the null
#' @export
#'
#' @examples
#' subjectiveBfZStat(5.2, 2, 3, 3.4, 2, 12)
subjectiveBfZStat <- function(x1, s21, n1, x2, s22, n2, sigma=1,
                              a1=5, b1=2, a2=-5, b2=2, a0=0, b0=1,
                              log=FALSE) {
  nCombined <- n1+n2
  xCombined <- (n1*x1+n2*x2)/nCombined

  if (n1 <= 1 && n2 <= 1) {
    s2Combined <- stats::var(c(x1, x2))
    s21 <- 0
    s22 <- 0
  } else {
    s2Combined <- ((n1-1)*s21+n1*x1^2+(n2-1)*s22+n2*x2^2-(nCombined)*xCombined^2)/(nCombined-1)
  }

  logMarg1 <- zMarg("n"=n1, "x"=x1, "s2"=s21, "a"=a1, "b"=b1, "sigma"=sigma) +
    zMarg("n"=n2, "x"=x2, "s2"=s22, "a"=a2, "b"=b2, "sigma"=sigma)
  logMarg0 <- zMarg("n"=nCombined, "x"=xCombined, "s2"=s2Combined, "a"=a0, "b"=b0, "sigma"=sigma)

  logBf10 <- logMarg1 - logMarg0

  if (isTRUE(log))
    return(logBf10)
  else
    return(exp(logBf10))

}
