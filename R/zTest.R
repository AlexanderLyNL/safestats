#' Computes S-Values Based on the Z-Statistic
#'
#' Computes s-values using the z-statistic and the sample sizes only based on the test defining parameter phiS.
#'
#' @param z numeric that represents the observed z-statistic.
#' @param parameter numeric this defines the safe test S, i.e., a likelihood ratio of z distributions with in the
#' denominator the likelihood with mean difference 0 and in the numerator an average likelihood defined by
#' the likelihood at the parameter value. For the two sided case 1/2 at the parameter value and 1/2 at minus the
#' parameter value.
#' @param n1 integer that represents the size in a one-sample z-test, (n2=\code{NULL}). When n2 is not
#' \code{NULL}, this specifies the size of the first sample for a two-sample test.
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
#' \code{NULL} it implies that the z-statistic is based on one-sample.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less".
#' @param paired a logical, if \code{TRUE} ignores n2, and indicates that a paired z-test is performed.
#' @param sigma numeric, the assumed known standard deviation, default 1.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an s-value.
#'
#' @export
#'
#' @examples
#' safeZTestStat(z=1, n1=100, parameter=0.4)
#' safeZTestStat(z=3, n1=100, parameter=0.3)
safeZTestStat <- function(z, parameter, n1, n2=NULL, alternative=c("two.sided", "less", "greater"),
                          paired=FALSE, sigma=1, ...) {
  alternative <- match.arg(alternative)

  phiS <- parameter

  if (is.null(n2) || is.na(n2) || paired==TRUE)
    nEff <- n1
  else
    nEff <- (1/n1+1/n2)^(-1)

  if (alternative=="two.sided") # two-sided
    result <- exp(-nEff*phiS^2/(2*sigma^2))*cosh(sqrt(nEff)*phiS/sigma*z)
  else # one-sided
    result <- exp(-1/2*(nEff*phiS^2/sigma^2-2*sqrt(nEff)*phiS/sigma*z))

  if (result < 0) {
    warning("Overflow: s-value smaller than 0")
    result <- 2^(-15)
  }

  return(unname(result))
}

#' Computes the Inverse of the Two-Sided Safe Z-Test
#'
#' This helper function is used in \code{\link{designSafeZ}} to find parameter. The function is the (two-sided)
#' inverse of safeZTestStat.
#'
#' @inheritParams safeZTestStat
#' @inheritParams designSafeZ
#'
#' @param nEff numeric > 0, the effective sample size.
#'
#' @return A number that represents a z-value. The function's domain is the positive real line and the range
#' is the real line, i.e., the outcome space of the z-statistic.
#' @export
#'
#' @examples
#'safeZ10Inverse(0.4, n=13)
safeZ10Inverse <- function(parameter, nEff, sigma=1, alpha=0.05) {
  phiS <- parameter
  sigma/(sqrt(nEff)*phiS)*acosh(exp(nEff*phiS^2/(2*sigma^2))/alpha)
}



#' Safe Z-Test
#'
#' Safe one and two sample z-tests on vectors of data. The function is modelled after \code{\link[stats]{t.test}}.
#'
#' @aliases safe.z.test
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param h0 a number indicating the hypothesised true value of the (differences of the) mean(s) under the null.
#' Default 0
#' @param paired a logical indicating whether you want the paired z-test.
#' @param designObj an object obtained from \code{\link{designSafeZ}}, or \code{NULL}, when pilot is set to \code{TRUE}.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the observed
#' number of samples is exactly as planned, and a pilot design object is constructed.
#' @param alpha numeric > 0 only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the alpha of
#' the design object is used instead in constructing the decision rule S > 1/alpha.
#' @param sigma numeric > 0 only used if pilot equals \code{TRUE}, which is then used to construct the design objects.
#' If pilot equals \code{FALSE} the sigma specified by the design object is used instead.
#' @param alternative a character only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the
#' alternative specified by the design object is used instead.
#' @param tol numeric > 0, only used if pilot equals \code{TRUE}, as it then specifies the mesh used to find the test
#' defining parameter to construct a pilot design object.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#'
#' \describe{
#'   \item{statistic}{the value of the test statistic. Here the z-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{sValue}{the s-value of the safe test.}
#'   \item{confInt}{To be implemented: a safe confidence interval for the mean appropriate to the specific alternative
#'   hypothesis.}
#'   \item{estimate}{the estimated mean or difference in means or mean difference depending on whether it was a one-
#'   sample test or a two-sample test.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on whether it was a one-sample
#'   or a two-sample test.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" effectively provided by the user.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" described in \code{\link{designSafeZ}}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeZ(meanDiffMin=0.6, alpha=0.008, beta=0.2,
#' alternative="greater", testType="twoSample", ratio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safeZTest(x, y, designObj=designObj)      #
#'
#' safeZTest(1:10, y = c(7:20), pilot=TRUE, alternative="less")      # s = 7.7543e+20 > 1/alpha
safeZTest <- function(x, y=NULL, h0=0, paired=FALSE, designObj=NULL,
                      pilot=FALSE, alpha=NULL, sigma=NULL,
                      alternative=NULL, tol=1e-05, ...) {

  result <- list("statistic"=NULL, "n"=NULL, "sValue"=NULL, "confSeq"=NULL, "estimate"=NULL,
                 "testType"=NULL, "dataName"=NULL, "h0"=h0, "sigma"=NULL, "call"=sys.call())
  class(result) <- "safeTest"

  if (is.null(designObj) && !pilot)
    stop("Please provide a safe z-test design object, or run the function with pilot=TRUE. ",
         "A design object can be obtained by running designSafeZ().")

  if (!is.null(designObj)) {
    checkDoubleArgumentsDesignObject(designObj, "alternative"=alternative, "alpha"=alpha, "sigma"=sigma)

    if (names(designObj[["parameter"]]) != "phiS")
      warning("The provided design is not constructed for the z-test,",
              "please use designSafeZ() instead. The test results might be invalid.")

  }

  if (is.null(y)) {
    testType <- "oneSample"
    nEff <- n1 <- length(x)
    n2 <- NULL

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    meanStat <- estimate <- mean(x)
    names(estimate) <- "mean of x"
  } else {
    nEff <- n1 <- length(x)
    n2 <- length(y)

    if (paired) {
      if (n1 != n2)
        stop("Data error: Error in complete.cases(x, y): Paired analysis requested, ",
             "but the two samples are not of the same size.")

      testType <- "paired"

      meanStat <- estimate <- mean(x-y)
      names(estimate) <- "mean of the differences"
    } else {
      testType <- "twoSample"

      nEff <- (1/n1+1/n2)^(-1)
      estimate <- c(mean(x), mean(y))
      names(estimate) <- c("mean of x", "mean of y")
      meanStat <- estimate[1]-estimate[2]
    }
  }

  if (pilot) {
    if (is.null(alternative))
      alternative <- "two.sided"
    else {
      if (!(alternative %in% c("two.sided", "greater", "less")))
        stop('Provided alternative must be one of "two.sided", "greater", or "less".')
    }

    if (is.null(alpha))
      alpha <- 0.05

    if (is.null(sigma))
      sigma <- 1

    nPlan <- if (is.null(n2)) n1 else c(n1, n2)
    designObj <- designPilotSafeZ("alpha"=alpha, "nPlan"=nPlan, "alternative"=alternative,
                                  "sigma"=sigma, "paired"=paired, "tol"=tol)
    designObj[["pilot"]] <- TRUE
  }

  alpha <- designObj[["alpha"]]
  sigma <- designObj[["sigma"]]
  alternative <- designObj[["alternative"]]

  zStat <- tryOrFailWithNA(sqrt(nEff)*(meanStat - h0)/sigma)

  if (is.na(zStat))
    stop("Could not compute the z-statistic")

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  sValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=n1, "n2"=n2,
                          "alternative"=alternative, "paired"=paired)

  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  result[["testType"]] <- testType
  result[["statistic"]] <- zStat
  result[["estimate"]] <- estimate
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj

  # result[["confSeq"]] <- computeZConfidenceSequence("nEff"=nEff, "meanStat"=meanStat,
  #                                                   "phiS"=abs(designObj[["parameter"]]),
  #                                                   "sigma"=sigma, "alpha"=alpha,
  #                                                   "alternative"=alternative)

  if (is.null(n2)) {
    result[["n"]] <- n1
    names(result[["n"]]) <- "n1"
  } else {
    result[["n"]] <- c(n1, n2)
    names(result[["n"]]) <- c("n1", "n2")
  }

  result[["sValue"]] <- sValue

  names(result[["h0"]]) <- "mu"
  names(result[["statistic"]]) <- "z"

  return(result)
}

#' Alias for safeZTest
#'
#' @rdname safeZTest
#'
#' @export
safe.z.test <- function(x, y=NULL, h0=0, paired=FALSE, designObj=NULL,
                        pilot=FALSE, alpha=0.05, sigma=1,
                        alternative=c("two.sided", "greater", "less"),
                        tol=1e-05, ...) {
  result <- safeZTest("x"=x, "y"=y, "designObj"=designObj, "alternative"=alternative, "h0"=h0,
                      "paired"=paired, "sigma"=sigma, "pilot"=pilot, "alpha"=alpha, ...)

  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  result[["dataName"]] <- dataName
  return(result)
}

#' Design a Frequentist Z-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I and type II error for the
#' frequentist z-test.
#'
#' @inheritParams designSafeZ
#' @param lowN the smallest number of samples (first group) at which monitoring of the tests begins.
#'
#' @return returns a freqZDesign object.
#' @export
#'
#' @examples
#' freqDesign <- designFreqZ(meanDiffMin = 0.5, highN = 100)
#' freqDesign$nPlan
#' freqDesign2 <- designFreqZ(meanDiffMin = 0.2, lowN = 32, highN = 200)
#' freqDesign2$nPlan
designFreqZ <- function(meanDiffMin, alternative=c("two.sided", "greater", "less"),
                        alpha=0.05, beta=0.2, testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, sigma=1, kappa=sigma, lowN=3L, highN=100L, ...) {

  stopifnot(lowN >= 1, highN > lowN, alpha > 0, beta >0)

  testType <- match.arg(testType)
  alternative <- match.arg(alternative)

  result <- list("nPlan"=NA, "esMin"=meanDiffMin, "alpha"=alpha, "beta"=beta,
                 "lowN"=lowN, "highN"=highN, "testType"=testType, "alternative"=alternative)
  class(result) <- "freqZDesign"

  if (meanDiffMin < 0 && alternative=="greater")
    warning("meanDiffMin < 0, but in the calculations abs(meanDiffMin) is used instead.")

  meanDiffMin <- abs(meanDiffMin)

  n1Plan <- NULL
  n2Plan <- NULL

  if (alternative=="two.sided")
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
        n2Plan <- ceiling(ratio*n)

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

#' Designs a Safe Z-Test Based on Planned Samples nPlan
#'
#' Designs a safe experiment for a prespecified tolerable type I error based on planned sample size(s),
#' which are fixed ahead of time. Outputs a list that includes phiS, i.e., the safe test defining parameter.
#'
#' @inheritParams designSafeZ
#' @param paired logical, if \code{TRUE} then paired z-test.
#'
#' @return Returns a safeDesign object
#' \describe{
#'   \item{nPlan}{the sample size(s) to plan for. Provided by the user.}
#'   \item{parameter}{the safe test defining parameter. Here phiS.}
#'   \item{esMin}{\code{NULL} no minimally clinically relevant effect size provided.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{\code{NULL}, no tolerable type II error specified.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" effectively provided by the user.}
#'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
#'   \item{sigma}{the assumed population standard deviation used for the test provided by the user.}
#'   \item{kappa}{the true population standard deviation, typically, sigma=kappa.}
#'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
#'   ratio between the planned randomisation of condition 2 over condition 1.}
#'   \item{tol}{the step size between parameter values in the candidate space.}
#'   \item{pilot}{logical, specifying whether it's a pilot design.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designPilotSafeZ(nPlan=30, alpha = 0.05)
designPilotSafeZ <- function(nPlan, alternative=c("two.sided", "greater", "less"),
                             alpha=0.05, sigma=1, kappa=sigma, tol=1e-5, paired=FALSE) {

  alternative <- match.arg(alternative)
  stopifnot(all(nPlan > 0))

  if (length(nPlan)==1) {
    n1 <- nPlan
    n2 <- NULL
    names(nPlan) <- "n1Plan"
  } else if (length(nPlan)==2) {
    n1 <- nPlan[1]
    n2 <- nPlan[2]
    names(nPlan) <- c("n1Plan", "n2Plan")
  }

  if (!is.null(n2)) {
    ratio <- n2/n1

    if (paired) {
      if (n1 != n2)
        stop("Paired design specified, but nPlan[1] not equal nPlan[2]")

      testType <- "paired"
      nEff <- n1
    } else {
      nEff <- (1/n1+1/n2)^(-1)
      testType <- "twoSample"
    }
  } else {
    nEff <- n1
    ratio <- 1
    testType <- "oneSample"

    if (isTRUE(paired)) {
      n2 <- n1
      warning("Paired designed specified, but nPlan[2] not provided. nPlan[2] is set to nPlan[1]")
    }
  }

  result <- list("nPlan"=nPlan, "parameter"=NULL, "esMin"=NULL, "alpha"=alpha, "beta"=NULL,
                 "alternative"=alternative, "testType"=testType, "paired"=paired,
                 "sigma"=sigma, "kappa"=kappa, "ratio"=ratio, "tol"=tol, "pilot"=FALSE,
                 "call"=sys.call(), "timeStamp"=Sys.time())

  class(result) <- "safeDesign"

  phiSPlus0 <- sigma*sqrt(2/nEff*log(1/alpha))

  if (alternative == "two.sided") {
    phiS10 <- sigma*sqrt(2/nEff*log(2/alpha))

    candidatePhi <- seq(phiSPlus0, phiS10, by=tol)

    safeZInverseValues <- purrr::map_dbl(candidatePhi, safeZ10Inverse, "nEff"=nEff, "sigma"=sigma, "alpha"=alpha)

    phiIndex <- which.min(safeZInverseValues)

    result[["parameter"]] <- candidatePhi[phiIndex]
  } else {

    if (alternative=="less")
      phiSPlus0 <- -phiSPlus0

    result[["parameter"]] <- phiSPlus0
  }

  names(result[["parameter"]]) <- "phiS"
  return(result)
}

#' Designs a Safe Z Experiment
#'
#' A designed experiment requires (1) a sample size nPlan to plan for, and (2) the parameter of the safe test, i.e.,
#' phiS. If nPlan is provided, then only the safe test defining parameter phiS needs to be determined. That resulting
#' phiS leads to an (approximately) most powerful safe test. Typically, nPlan is unknown and the user has to specify
#' (i) a tolerable type II error beta, and (b) a clinically relevant minimal population mean difference meanDiffMin.
#' The procedure finds the smallest nPlan for which meanDiffMin is found with power of at least 1 - beta.
#'
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control --independent on n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule S10 > 1/alpha.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both "n"
#' and "phiS". Note that 1-beta defines the power.
#' @param meanDiffMin numeric that defines the minimal relevant mean difference, the smallest population mean
#' that we would like to detect.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less".
#' @param nPlan optional numeric vector of length at most 2. When provided, it is used to find the safe test
#' defining parameter phiS. Note that if the purpose is to plan based on nPlan alone, then both meanDiffMin
#' and beta should be set to NULL.
#' @param sigma numeric > 0 representing the assumed population standard deviation used for the test.
#' @param kappa the true population standard deviation. Default kappa=sigma.
#' @param tol a number that defines the stepsizes between the lowParam and highParam.
#' @param highN integer that stops the search process if the lower bound for the candidate nPlan exceeds this
#' number. Default highN is set 8e9 to correspond to the current population of the earth.
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
#' is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param parameter optional test defining parameter. Default set to \code{NULL}.
#' @param grow logical, defaul set to \code{TRUE} so the grow safe test is used in the design.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a safeDesign object that includes:
#'
#' \describe{
#'   \item{nPlan}{the sample size(s) to plan for. Computed based on beta and meanDiffMin, or provided by the user
#'   if known.}
#'   \item{parameter}{the safe test defining parameter. Here phiS.}
#'   \item{esMin}{the minimally clinically relevant effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error specified by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" effectively provided by the user.}
#'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
#'   \item{sigma}{the assumed population standard deviation used for the test provided by the user.}
#'   \item{kappa}{the true population standard deviation, typically, sigma=kappa.}
#'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
#'   ratio between the planned randomisation of condition 2 over condition 1.}
#'   \item{tol}{the step size between parameter values in the candidate space.}
#'   \item{pilot}{logical, specifying whether it's a pilot design.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeZ(meanDiffMin=0.8, alpha=0.08, beta=0.01, alternative="greater")
#'
#' #nPlan known:
#' designObj <- designSafeZ(nPlan = 100, alpha=0.05)
#'
designSafeZ <- function(meanDiffMin=NULL, beta=NULL, nPlan=NULL,
                        alternative=c("two.sided", "greater", "less"),
                        alpha=0.05, sigma=1, kappa=sigma, tol=1e-5, highN=8e9,
                        testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, parameter=NULL, grow=TRUE, ...) {

  stopifnot(alpha > 0, alpha < 1, sigma > 0, kappa > 0)

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  paired <- if (testType=="paired") TRUE else FALSE

  designScenario <- NULL

  tempResult <- list()

  if (!is.null(meanDiffMin) && !is.null(beta) && is.null(nPlan)) {
    # Scenario "1a"
    designScenario <- "1a"

    tempResult <- computeZSafeTestAndNFrom("meanDiffMin"=meanDiffMin, "beta"=beta, "alpha"=alpha,
                                           "sigma"=sigma, "kappa"=kappa, "alternative"=alternative,
                                           "testType"=testType, "tol"=tol, "highN"=highN, "ratio"=ratio,
                                           "parameter"=parameter, "designScenario"=designScenario,
                                           "grow"=grow)
    nPlan <- tempResult[["nPlan"]]
    phiS <- tempResult[["phiS"]]
  } else if (!is.null(meanDiffMin) && is.null(beta) && is.null(nPlan)) {
    # Scenario "1b"
    designScenario <- "1b"

    tempResult <- computeZSafeTestAndNFrom("meanDiffMin"=meanDiffMin, "beta"=beta, "alpha"=alpha,
                                           "sigma"=sigma, "kappa"=kappa, "alternative"=alternative,
                                           "testType"=testType, "tol"=tol, "highN"=highN, "ratio"=ratio,
                                           "parameter"=parameter, "designScenario"=designScenario)
    nPlan <- NULL
    phiS <- tempResult[["phiS"]]
    beta <- NULL
    meanDiffMin <- meanDiffMin
  } else if (is.null(meanDiffMin) && is.null(beta) && !is.null(nPlan)) {
    # Scenario "1c"
    designScenario <- "1c"

    return(designPilotSafeZ("nPlan"=nPlan, "alpha"=alpha, "alternative"=alternative,
                            "sigma"=sigma, "kappa"=kappa, "tol"=tol, "paired"=paired))
  } else if (!is.null(meanDiffMin) && is.null(beta) && !is.null(nPlan)) {
    # Scenario 2
    designScenario <- "2"

    beta <- tryOrFailWithNA(
      computeZBetaFrom("meanDiffMin"=meanDiffMin, "nPlan"=nPlan, "alpha"=alpha, "sigma"=sigma, "kappa"=kappa,
                       "alternative"=alternative, "testType"=testType, "parameter"=parameter)
    )

    phiS <- if (is.null(parameter)) meanDiffMin else parameter

  } else if (is.null(meanDiffMin) && !is.null(beta) && !is.null(nPlan)) {
    # Scenario 3
    designScenario <- "3"

    meanDiffMin <- tryOrFailWithNA(
      computeZMeanDiffMinFrom("nPlan"=nPlan, "alpha"=alpha, "beta"=beta, "sigma"=sigma,
                              "kappa"=kappa, "alternative"=alternative, "testType"=testType,
                              "parameter"=parameter)
    )
    phiS <- if (is.null(parameter)) meanDiffMin else parameter
  }

  if (is.null(designScenario)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1) non-null meanDiffMen, non-null beta and NULL nPlan, or \n",
         "(2) non-null meanDiffMen, NULL beta and non-null meanDiffMen, or \n",
         "(3) NULL meanDiffMen, non-null beta, and non-null nPlan, or \n",
         "(4) non-null meanDiffMen, NULL beta, and NULL nPlan, or \n",
         "(5) NULL meanDiffMen, NULL meanDiffMen, non-null nPlan.")
  }

  if (is.na(meanDiffMin))
    meanDiffMin <- NULL

  if (designScenario %in% 2:3) {
    n2Plan <- nPlan[2]

    names(nPlan) <- if (is.na(n2Plan)) "n1Plan" else c("n1Plan", "n2Plan")
  }

  result <- list("nPlan"=nPlan, "parameter"=phiS, "esMin"=meanDiffMin, "alpha"=alpha, "beta"=beta,
                 "alternative"=alternative, "testType"=testType, "paired"=paired,
                 "sigma"=sigma, "kappa"=kappa,
                 "ratio"=ratio, "pilot"=FALSE, "lowN"=NULL, "highN"=NULL, "call"=sys.call(),
                 "timeStamp"=Sys.time())
  class(result) <- "safeDesign"

  names(result[["esMin"]]) <- "mean difference"
    # switch(alternative, "two.sided"="mean differences at least abs(phi)",
    #        "less"="mean differences smaller than phi",
    #        "greater"="mean differences larger than phi")

  for (neem in c("lowN", "highN", "lowParam", "highParam", "nEffPlan")) {
    value <- tempResult[[neem]]

    if (!is.null(value))
      result[[neem]] <- value
  }

  names(result[["parameter"]]) <- "phiS"
  return(result)
}


#' Helper function: Computes the planned sample size based on the minimal clinical relevant mean
#' difference, alpha and beta
#'
#' @inheritParams  designSafeZ
#' @param designScenario a character string specifying the scenario for which needs designing either
#' "1a" or "1b"
#'
#' @return a list which contains at least nPlan and the phiS the parameter that defines the safe test
#'
#' @examples
#' safestats:::computeZSafeTestAndNFrom(0.4)
#' safestats:::computeZSafeTestAndNFrom(0.4, grow=FALSE)
computeZSafeTestAndNFrom <- function(meanDiffMin, alpha=0.05, beta=0.2, sigma=1, kappa=sigma,
                                     alternative=c("two.sided", "greater", "less"),
                                     testType=c("oneSample", "paired", "twoSample"),
                                     tol=1e-5, highN=8e9, ratio=1, parameter=NULL,
                                     designScenario="1a", grow=TRUE) {
  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  result <- list(nPlan=NULL, phiS=NULL)
  meanDiffMin <- abs(meanDiffMin)

  n1Plan <- NULL
  n2Plan <- NULL

  nEffToN1Ratio <- if (testType=="twoSample") (1+ratio)/ratio else 1

  if (designScenario=="1a") {
    if (grow) {
      phiS <- abs(meanDiffMin)

      if (alternative == "two.sided") {
        criterionFunction <- function(n) {
          lowerTail <- sigma^4/(n*kappa^2*meanDiffMin^2)*acosh(exp((n*meanDiffMin^2)/(2*sigma^2))/alpha)^2
          stats::pchisq(q=lowerTail, df=1, ncp=n*meanDiffMin^2/kappa^2)-beta
        }

        highN <- 2*sigma^2/meanDiffMin^2*log(1e100)
        tempResult <- stats::uniroot(criterionFunction, interval=c(1, highN))
        nEff <- tempResult[["root"]]
      } else {
        qB <- qnorm(beta)

        nEff <- exp(2*(log(kappa)-log(meanDiffMin))) *
          (2*qB^2 - 2*qB*sqrt(qB^2+2*sigma^2/kappa^2*log(1/alpha))+2*kappa^2/sigma^2*log(1/alpha))
      }

      if (testType == "twoSample") {
        n1Plan <- ceiling(nEff * nEffToN1Ratio)
        n2Plan <- ceiling(nEff * nEffToN1Ratio * ratio)
      } else {
        n1Plan <- ceiling(nEff)
        n2Plan <- if (testType == "paired") n1Plan else NULL
      }
    } else {
      # Note(Alexander): Compute one-sided nExact. This provides us with a lower bound on
      # the two-sided test.
      #
      nEffExact <- tryOrFailWithNA(((sigma*sqrt(2*log(1/alpha))-kappa*qnorm(beta))/meanDiffMin)^2)

      if (is.na(nEffExact))
        stop("Something went wrong, couldn't design based on the given input.")

      if (nEffExact > highN)
        stop("More samples needed than highN, which is ", highN)

      if (alternative %in% c("greater", "less")) {
        # Note(Alexander): Here I use nEff exact, not ceiling(nEff), which should be an integer if testType != "twoSample"
        # when testType == "twoSample" I do take nEff <- (1/n1Plan + 1/n2Plan), where n1Plan and n2Plan are integers
        # This means that discriminant D is very close to 0. I tried using the exact nEff, but this performed less well
        # for the actual sample sizes.
        #
        #
        if (testType == "twoSample") {
          n1Plan <- ceiling(nEffExact*nEffToN1Ratio)
          n2Plan <- ceiling(nEffExact*nEffToN1Ratio*ratio)
          nEff <- (1/n1Plan+1/n2Plan)^(-1)
        } else {
          n1Plan <- ceiling(nEffExact)

          if (testType == "paired")
            n2Plan <- n1Plan

          nEff <- ceiling(nEffExact)
        }

        qBeta <- kappa/sigma*qnorm(beta) + sqrt(nEff)*meanDiffMin/sigma
        discriminantD <- max(qBeta^2-2*log(1/alpha), 0)

        phiS <- sigma/sqrt(nEff)*(qBeta + sqrt(discriminantD))
      } else {
        # Two.sided

        nEffExactUpper <- tryOrFailWithNA(
          ((sigma*sqrt(2*log(2/alpha))-kappa*qnorm(beta))/meanDiffMin)^2
        )

        # Note(Alexander): Translate to lower and upper bound in terms of n1
        #
        lowN <- floor(nEffExact*nEffToN1Ratio)
        highN <- ceiling(nEffExactUpper*nEffToN1Ratio)

        # This shouldn't occur
        if (is.na(highN))
          highN <- 2*lowN

        result[["lowN"]] <- lowN
        result[["highN"]] <- highN

        # Note(Alexander): This function is used to
        #   1. Find n and creates is looped downwards, hence the  reflection ! in result
        #   2. Given n find the parameter phiS
        #
        criterionFunctionExact <- function(n, parameter=NULL) {
          if (is.null(parameter))
            parameter <- sqrt(2*sigma^2/n*log(2/alpha))

          zArg <- sigma^4/(kappa^2*n*parameter^2)*(acosh(exp(n*parameter^2/(2*sigma^2))/alpha))^2

          result <- (stats::pchisq(zArg, df=1, ncp=n*meanDiffMin^2/kappa^2) <= beta)

          if (is.null(parameter))
            result <- !result

          return(result)
        }

        candidateN1 <- highN
        candidateNEff <- candidateN1/nEffToN1Ratio

        continueWhile <- TRUE

        # Note(Alexander): Loop backwards to find a smaller n1Plan
        #
        while (continueWhile && candidateN1 > lowN) {
          continueWhile <- criterionFunctionExact(n=candidateNEff, parameter=NULL)

          if (isTRUE(continueWhile)) {
            candidateN1 <- candidateN1 - 1
            candidateNEff <- candidateN1/nEffToN1Ratio
          } else {
            candidateN1 <- candidateN1 + 1
            break()
          }
        }

        nEff <- candidateN1/nEffToN1Ratio

        if (testType=="twoSample") {
          n1Plan <- ceiling(candidateN1)
          n2Plan <- ceiling(n1Plan*ratio)
          nEff <- (1/n1Plan+1/n2Plan)^(-1)
          result[["nEffPlan"]] <- nEff
        } else {
          n1Plan <- ceiling(nEff)

          if (testType=="pairedSample")
            n2Plan <- n1Plan

        }

        # Candidate parameters ---
        phiUmp <- sqrt(2/(sigma^2*nEff)*log(2/alpha))

        chiSqInverseBeta <- stats::qchisq(beta, df=1, ncp=nEff*meanDiffMin^2/kappa^2)
        discriminantD <- max(chiSqInverseBeta-sigma^2/kappa^2*2*log(2/alpha), 0)

        phiSApprox <-kappa/sqrt(nEff)*(sqrt(chiSqInverseBeta)+sqrt(discriminantD))

        # Random lower bound for phi
        # TODO(Alexander): Get a better bounds perhaps
        #
        lowPhi <- min(phiUmp, phiSApprox, meanDiffMin/2)
        highPhi <- if (beta < 1/2) meanDiffMin else 2*meanDiffMin

        candidatePhis <- seq(lowPhi, highPhi, "by"=tol)

        phiIndex <- purrr::detect_index(candidatePhis, criterionFunctionExact, "n"=nEff)

        phiS <- if (phiIndex==0) phiSApprox else candidatePhis[phiIndex]

        result[["lowParam"]] <- lowPhi
        result[["highParam"]] <- highPhi
      } # end two.sided
    }
  } else if (designScenario=="1b") {
    sideConstant <- if (alternative %in% c("greater", "less")) 1 else 2

    phiS <- if (!is.null(parameter)) parameter else meanDiffMin

    nEffExact <- 2*sigma^2*log(sideConstant/alpha)/phiS^2

    if (testType=="twoSample") {
      n1Plan <- ceiling(nEffExact*(1+ratio)/ratio)
      n2Plan <- ceiling(nEffExact*(1+ratio))
      nEff <- (1/n1Plan + 1/n2Plan)^(-1)
    } else {
      nEff <- nEffExact
      n1Plan <- ceiling(nEff)

      if (testType=="paired")
        n2Plan <- n1Plan

    }
    #
    if (alternative=="two.sided") {
      lowerTail <- sigma^4/(nEff*kappa^2*phiS^2)*(acosh(exp(nEff*phiS^2/(2*sigma^2))/alpha))^2
      result[["beta"]] <- stats::pchisq(q=lowerTail, df=1, ncp=nEff*meanDiffMin^2/kappa^2)
    } else {
      lowerTail <- sqrt(nEff)*(phiS-2*meanDiffMin)/(2*kappa) -
        sigma^2*log(alpha)/(kappa*sqrt(nEff))*1/phiS
      result[["beta"]] <- pnorm(lowerTail)
    }
  }

  if (alternative=="less")
    phiS <- - phiS

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  result[["phiS"]] <- phiS

  return(result)
}


#' Help function to compute the effecitve sample size based on a length 2 vector of samples
#'
#' @inheritParams designSafeZ
#' @param n vector of length at most 2 representing the sample sizes of the first and second group
#' @param silent logical, if true, then turn off warnings
#'
#' @return a numeric that represents the effective sample size.
#'
#' @examples
#' safestats:::computeNEff(c(3, 4), testType="twoSample")
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




#' Computes the smallest mean difference that is detectable with chance 1-beta, for the provided
#' sample size
#'
#' @inheritParams designSafeZ
#' @param maxIter maximum number of iterations in the optimisation process for two-sided designs
#'
#' @return numeric > 0 that represents the minimal detectable mean difference
#'
#' @examples
#' safestats:::computeZMeanDiffMinFrom(nPlan=78)
computeZMeanDiffMinFrom <- function(nPlan, alpha=0.05, beta=0.2, sigma=1, kappa=sigma,
                                    alternative=c("two.sided", "greater", "less"),
                                    testType=c("oneSample", "paired", "twoSample"),
                                    parameter=NULL, maxIter=10) {
  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  if (alternative=="two.sided") {
    if (is.null(parameter)) {
      criterionFunction <- function(x) {
        lowerTail <- exp(4*log(sigma)-log(nEff)-2*log(kappa)-2*log(x)) *
          (acosh(exp(nEff*x^2/(2*sigma^2))/alpha))^2

        stats::pchisq(lowerTail, 1, nEff*x^2/kappa^2)-beta
      }
    } else {
      criterionFunction <- function(x) {
        lowerTail <- exp(4*log(sigma)-log(nEff)-2*log(kappa)-2*log(parameter)) *
          (acosh(exp(nEff*parameter^2/(2*sigma^2))/alpha))^2

        stats::pchisq(lowerTail, 1, nEff*x^2/kappa^2)-beta
      }
    }

    mIter <- 1
    tempResult <- 1
    class(tempResult) <- "try-error"

    while (isTryError(tempResult) && mIter <= maxIter) {
      tempResult <- try(uniroot(criterionFunction, c(0, 4*4^(-mIter+1))), silent=TRUE)

      if (isTryError(tempResult))
        mIter <- mIter + 1
      else
        break()
    }

    if (mIter > maxIter || isTryError(tempResult))
      stop("uniroot couldn't find meanDiffMin. Perhaps maxIter not big enough, ",
      "but most likely pchisq underflow")
    else
      result <- tempResult[["root"]]
  } else {
    if (is.null(parameter)) {
      D <- kappa^2*qnorm(beta)^2 - 2 * sigma^2*log(alpha)
      result <- 1/sqrt(nEff)*(-kappa*qnorm(beta) + sqrt(D))
    } else {
      result <- parameter/2-sigma^2*log(alpha)/(nEff*parameter)-kappa/(sqrt(nEff))*qnorm(beta)
    }
  }

  if (alternative=="less")
    result <- -result

  return(result)
}

#' Helper function: Computes the type II error based on the minimal clinically relevant effect size and sample size.
#'
#' @inheritParams designSafeZ
#'
#' @return numeric that represents the type II error
#'
#' @examples
  #' safestats:::computeZBetaFrom(meanDiffMin=0.9, nPlan=12)
computeZBetaFrom <- function(meanDiffMin, nPlan, alpha=0.05, sigma=1, kappa=sigma,
                             alternative=c("two.sided", "greater", "less"),
                             testType=c("oneSample", "paired", "twoSample"),
                             parameter=NULL) {

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  if (is.null(parameter))
    parameter <- meanDiffMin

  if (alternative=="two.sided") {
    lowerTail <- exp(4*log(sigma)-log(nEff)-2*log(kappa)-2*log(parameter)) *
      (acosh(exp(nEff*parameter^2/(2*sigma^2))/alpha))^2

    result <- stats::pchisq(lowerTail, df=1, ncp=nEff*meanDiffMin^2/kappa^2, lower.tail = TRUE)
  } else {
    lowerTail <- sqrt(nEff)*(parameter-2*meanDiffMin)/(2*kappa) -
      sigma^2*log(alpha)/(kappa*sqrt(nEff))*1/parameter

    result <- stats::pnorm(lowerTail, lower.tail = TRUE)
  }

  return(result)
}

#' Helper function: Computes the safe confidence sequence for a z-test
#'
#' @inheritParams safeZTest
#' @param nEff numeric > 0, the effective sample size
#' @param meanStat numeric, the mean statistic, this could be the differences of means as well
#' @param phiS numeric > 0, the safe test defining parameter
#'
#' @return numeric vector that contains the upper and lower bound of the safe confidence sequence
#'
#' @examples
#' safestats:::computeZConfidenceSequence(nEff=15, meanStat=0.3, phiS=0.2)
computeZConfidenceSequence <- function(nEff, meanStat, phiS, sigma=1, alpha=0.05,
                                       alternative="two.sided") {
  # TODO(Alexander): Only for GROW,
  meanDiffMin <- phiS
  g <- meanDiffMin^2/sigma^2

  if (alternative=="two.sided") {
    shift <- sigma*sqrt((1+nEff*g)/(nEff^2*g)*(log(nEff*g)-2*log(alpha)))
    lowerCS <- meanStat - shift
    upperCS <- meanStat + shift
  } else {
    shift <- sigma * sqrt((1+nEff*g)/(nEff^2*g)*(log(nEff*g)-2*log(2*alpha)))

    if (alternative=="greater") {
      lowerCS <- meanStat + shift
      upperCS <- Inf
    } else {
      lowerCS <- -Inf
      upperCS <- meanStat - shift
    }
  }

  # if (alternative=="two.sided") {
  #   shift <- sigma^2/(nEff*phiS)*acosh(exp(nEff*phiS^2/(2*sigma^2))/alpha)
  #   lowerCS <- meanStat - shift
  #   upperCS <- meanStat + shift
  # } else {
  #   shift <- sigma^2/nEff*log(alpha)*1/phiS - phiS/2
  #
  #   if (alternative=="greater") {
  #     lowerCS <- meanStat + shift
  #     upperCS <- Inf
  #   } else {
  #     lowerCS <- -Inf
  #     upperCS <- meanStat - shift
  #   }
  # }
  return(unname(c(lowerCS, upperCS)))
}

