# Testing fnts ---------

#' Computes E-Values Based on the Z-Statistic
#'
#' Computes e-values using the z-statistic and the sample sizes only based on the test defining parameter phiS.
#'
#' @param z numeric that represents the observed z-statistic.
#' @param phiS numeric this defines the safe test S, i.e., a likelihood ratio of z distributions with in the
#' denominator the likelihood with mean difference 0 and in the numerator an average likelihood defined by
#' the likelihood at the parameter value. For the two sided case 1/2 at the parameter value and 1/2 at minus the
#' parameter value.
#' @param n1 integer that represents the size in a one-sample z-test, (n2=\code{NULL}). When n2 is not
#' \code{NULL}, this specifies the size of the first sample for a two-sample test.
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
#' \code{NULL} it implies that the z-statistic is based on one-sample.
#' @param alternative a character string specifying the alternative hypothesis must be one of "twoSided" (default),
#' "greater" or "less".
#' @param paired a logical, if \code{TRUE} ignores n2, and indicates that a paired z-test is performed.
#' @param sigma numeric, the assumed known standard deviation, default 1.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an e-value.
#'
#' @export
#'
#' @examples
#' safeZTestStat(z=1, n1=100, phiS=0.4)
#' safeZTestStat(z=3, n1=100, phiS=0.3)
safeZTestStat <- function(z, phiS, n1, n2=NULL,
                          alternative=c("twoSided", "less", "greater"),
                          paired=FALSE, sigma=1, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  if (is.null(n2) || is.na(n2) || paired==TRUE)
    nEff <- n1
  else
    nEff <- (1/n1+1/n2)^(-1)

  phiS <- checkAndReturnsEsMinParameterSide("paramToCheck"=phiS, "alternative"=alternative,
                                            "esMinName"="phiS")

  if (alternative=="twoSided") { # two-sided
    result <- exp(-nEff*phiS^2/(2*sigma^2))*cosh(sqrt(nEff)*phiS/sigma*z)
  } else { # one-sided
    result <- exp(-1/2*(nEff*phiS^2/sigma^2-2*sqrt(nEff)*phiS/sigma*z))
  }

  if (result < 0) {
    warning("Overflow: e-value smaller than 0")
    result <- 2^(-15)
  }

  return(unname(result))
}

#' Computes the Inverse of the Two-Sided Safe Z-Test
#'
#' This helper function is used in \code{\link{designSafeZ}()} to find parameter. The function is the (two-sided)
#' inverse of 'safeZTestStat'.
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
#' Safe one and two sample z-tests on vectors of data. The function is modelled after \code{\link[stats]{t.test}()}.
#'
#' @aliases safe.z.test
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param paired a logical indicating whether you want the paired z-test.
#' @param designObj an object obtained from \code{\link{designSafeZ}()}, or \code{NULL}, when pilot is set to \code{TRUE}.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the number of
#' samples is exactly as planned. The default null h0=1 is used, alpha=0.05, and alternative="twoSided" is used.
#' To change these default values, please use \code{\link{designSafeZ}()}.
#' @param ciValue numeric is the ciValue-level of the confidence sequence. Default ciValue=NULL, and ciValue = 1 - alpha
#' @param tol numeric > 0, only used if pilot equals \code{TRUE}, as it then specifies the mesh used to find the test
#' defining parameter to construct a pilot design object.
#' @param na.rm a logical value indicating whether \code{NA} values should be stripped before
#' the computation proceeds.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class 'safeTest'. An object of class 'safeTest' is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the test statistic. Here the z-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{eValue}{the e-value of the safe test.}
#'   \item{confInt}{To be implemented: a safe confidence interval for the mean appropriate to the specific alternative
#'   hypothesis.}
#'   \item{estimate}{the estimated mean or difference in means or mean difference depending on whether it was a one-
#'   sample test or a two-sample test.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on whether it was a one-sample
#'   or a two-sample test.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" effectively provided by the user.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" described in \code{\link{designSafeZ}()}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#'
#' designObj <- designSafeZ(meanDiffMin=0.6, alpha=0.008,
#'                          alternative="greater", testType="twoSample",
#'                          ratio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safeZTest(x, y, designObj=designObj)      #
#'
#' safeZTest(1:10, y = c(7:20), pilot=TRUE, alternative="less")      # s = 7.7543e+20 > 1/alpha
safeZTest <- function(x, y=NULL, paired=FALSE, designObj=NULL,
                      pilot=FALSE, ciValue=NULL, tol=1e-05, na.rm=FALSE, ...) {

  result <- list("statistic"=NULL, "n"=NULL, "eValue"=NULL, "confSeq"=NULL, "estimate"=NULL,
                 "testType"=NULL, "dataName"=NULL, "h0"=NULL, "sigma"=NULL, "call"=sys.call())
  class(result) <- "safeTest"

  if (is.null(designObj) && !pilot)
    stop("Please provide a safe z-test design object, or run the function with pilot=TRUE. ",
         "A design object can be obtained by running designSafeZ().")

  if (!is.null(designObj)) {
    if (names(designObj[["parameter"]]) != "phiS")
      warning("The provided design is not constructed for the z-test,",
              "please use designSafeZ() instead. The test results might be invalid.")
  }

  if (is.null(y)) {
    testType <- "oneSample"
    n <- nEff <- n1 <- length(x)
    n2 <- NULL

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    meanObs <- estimate <- mean(x, "na.rm"=na.rm)

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

      meanObs <- estimate <- mean(x-y, "na.rm"=na.rm)
      names(estimate) <- "mean of the differences"
    } else {
      testType <- "twoSample"

      nEff <- (1/n1+1/n2)^(-1)
      estimate <- c(mean(x, "na.rm"=na.rm), mean(y, "na.rm"=na.rm))
      names(estimate) <- c("mean of x", "mean of y")
      meanObs <- estimate[1]-estimate[2]
    }

    n <- c(n1, n2)
    names(n) <- c("n1", "n2")
  }

  if (pilot) {
    if (is.null(designObj)) {
      alternative <- "twoSided"
      alpha <- 0.05
      sigma <- 1

      nPlan <- if (is.null(n2)) n1 else c(n1, n2)
      designObj <- designPilotSafeZ("alpha"=alpha, "nPlan"=nPlan, "alternative"=alternative,
                                    "sigma"=sigma, "paired"=paired, "tol"=tol)
      designObj[["pilot"]] <- TRUE
    } else {
      warning("The pilot flag is ignored, since a designObj is given",
              "The analysis will be run based on the designObj.")
    }
  }

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

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  eValue <- safeZTestStat("z"=zStat, "phiS"=designObj[["parameter"]], "n1"=n1, "n2"=n2,
                          "alternative"=alternative, "paired"=paired)

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
  result[["designObj"]] <- designObj
  result[["ciValue"]] <- ciValue
  result[["n"]] <- n

  result[["confSeq"]] <- computeConfidenceIntervalZ("nEff"=nEff, "meanObs"=meanObs,
                                                    "phiS"=designObj[["parameter"]],
                                                    "sigma"=sigma, "ciValue"=ciValue,
                                                    "alternative"="twoSided")
  result[["eValue"]] <- eValue

  names(result[["statistic"]]) <- "z"

  return(result)
}

#' Alias for safeZTest
#'
#' @rdname safeZTest
#'
#' @export
safe.z.test <- function(x, y=NULL, paired=FALSE, designObj=NULL,
                        pilot=FALSE, tol=1e-05, ...) {

  result <- safeZTest("x"=x, "y"=y, "designObj"=designObj,
                      "paired"=paired, "pilot"=pilot, ...)
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

#' Helper function: Computes the safe confidence sequence for a z-test
#'
#' @inheritParams safeZTest
#' @inheritParams designSafeZ
#' @param nEff numeric > 0, the effective sample size.
#' @param meanObs numeric, the observed mean. For two sample tests this is difference of the means.
#' @param phiS numeric > 0, the safe test defining parameter.
#' @param ciValue numeric is the ciValue-level of the confidence sequence. Default ciValue=0.95.
#' @param a numeric, the centre of the normal prior on population mean (of the normal data). Default
#' is \code{NULL}, which implies the default choice of setting the centre equal to the null hypothesis.
#' @param g numeric > 0, used to define g sigma^2 as the variance of the normal prior on the population
#' (of the normal data). Default is \code{NULL} in which case g=phiS^2/sigma^2.
#' @param pointPrior logical, used to indicate whether a point prior is used that
#' puts half its mass on -phiS and the other half on phiS. Default \code{FALSE}.
#' @param freq logical, used to indicate whether the classical confidence interval
#' needs computing. Default \code{FALSE}.
#'
#' @return numeric vector that contains the upper and lower bound of the safe confidence sequence
#' @export
#'
#' @examples
#' computeConfidenceIntervalZ(nEff=15, meanObs=0.3, phiS=0.2)
computeConfidenceIntervalZ <- function(nEff, meanObs, phiS, sigma=1, ciValue=0.95,
                                       alternative="twoSided", a=NULL, g=NULL,
                                       pointPrior=FALSE, freq=FALSE, credibleInterval=FALSE) {
  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  if (isTRUE(freq)) {
    shift <- sigma/sqrt(nEff)*qnorm((1-ciValue)/2)
    lowerCi <- meanObs + shift
    upperCi <- meanObs - shift
    return(unname(c(lowerCi, upperCi)))
  }

  if (isTRUE(credibleInterval)) {
    normalisedCiLower <- qnorm((1-ciValue)/2)

    if (is.null(a)) {
      warning("Centre of normal prior not given, default to a=0")
      a <- 0
    }

    if (is.null(g)) {
      warning("Variance of normal prior not given, default to g=0")
      g <- 1
    }

    meanMu <- nEff*g/(nEff*g+sigma^2)*meanObs + sigma^2/(nEff*g+sigma^2)*a
    sdMu <- sqrt(g*sigma^2)/sqrt(nEff*g + sigma^2)

    lowerCi <- sdMu*normalisedCiLower+meanMu
    upperCi <- -sdMu*(normalisedCiLower)+meanMu

    return(unname(c(lowerCi, upperCi)))
  }

  if (isTRUE(pointPrior)) {
    shift <- sigma^2/(nEff*phiS)*acosh(exp(nEff*phiS^2/(2*sigma^2))/(1-ciValue))
    lowerCS <- meanObs - shift
    upperCS <- meanObs + shift
  } else {
    if (!is.null(a) && !is.null(g)) {
      # Note(Alexander): Here normal distribution not centred at null
      if (alternative != "twoSided")
        stop("One-sided confidence sequences for non-zero centred normal priors not implemented.")

      shift <- sqrt(sigma^2/nEff*(log(1+nEff*g)-2*log(1-ciValue))+(meanObs-a)^2/(1+nEff*g))
      lowerCS <- meanObs - shift
      upperCS <- meanObs + shift
    } else {
      # Note(Alexander): Here normal distribution centred at the null
      # Here use GROW
      if (is.null(g)) {
        meanDiffMin <- phiS
        g <- meanDiffMin^2/sigma^2
      }

      if (alternative=="twoSided") {
        shift <- sigma/(nEff*sqrt(g))*sqrt((1+nEff*g)*(log(1+nEff*g)-2*log(1-ciValue)))
        lowerCS <- meanObs - shift
        upperCS <- meanObs + shift
      } else {
        shift <- sigma/(nEff*sqrt(g))*sqrt((1+nEff*g)*(log(1+nEff*g)-2*log(2*(1-ciValue))))

        if (alternative=="greater") {
          lowerCS <- meanObs + shift
          upperCS <- Inf
        } else if (alternative=="less") {
          lowerCS <- -Inf
          upperCS <- meanObs - shift
        } else {
          stop('Incorrect specification of the "alternative" argument.')
        }
      }
    }
  }



  return(unname(c(lowerCS, upperCS)))
}


# Design fnts ---------

#' Design a Frequentist Z-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I and type II error for the
#' frequentist z-test.
#'
#' @inheritParams designSafeZ
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
designFreqZ <- function(meanDiffMin, alternative=c("twoSided", "greater", "less"),
                        alpha=0.05, beta=0.2, testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, sigma=1, h0=0, kappa=sigma, lowN=3L, highN=100L, ...) {

  stopifnot(lowN >= 1, highN > lowN, alpha > 0, beta >0)

  testType <- match.arg(testType)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  result <- list("nPlan"=NA, "esMin"=meanDiffMin, "alpha"=alpha, "beta"=beta,
                 "lowN"=lowN, "highN"=highN, "testType"=testType, "alternative"=alternative,
                 "h0"=h0)
  class(result) <- "freqZDesign"

  meanDiffMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=meanDiffMin, "alternative"=alternative,
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
#' @return Returns a 'safeDesign' object
#' \describe{
#'   \item{nPlan}{the sample size(s) to plan for. Provided by the user.}
#'   \item{parameter}{the safe test defining parameter. Here phiS.}
#'   \item{esMin}{\code{NULL} no minimally clinically relevant effect size provided.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{\code{NULL}, no tolerable type II error specified.}
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
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
designPilotSafeZ <- function(nPlan, alternative=c("twoSided", "greater", "less"),
                             alpha=0.05, sigma=1, h0=0, kappa=sigma, tol=1e-5,
                             paired=FALSE, parameter=NULL) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

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

  result <- list("nPlan"=nPlan, "parameter"=parameter, "esMin"=NULL, "alpha"=alpha, "beta"=NULL,
                 "h0"=h0, "sigma"=sigma, "alternative"=alternative, "testType"=testType,
                 "paired"=paired, "kappa"=kappa, "ratio"=ratio, "tol"=tol, "pilot"=FALSE,
                 "call"=sys.call(), "timeStamp"=Sys.time())

  class(result) <- "safeDesign"

  if (is.null(parameter)) {
    phiSPlus0 <- sigma*sqrt(2/nEff*log(1/alpha))

    if (alternative == "twoSided") {
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
  }

  names(result[["parameter"]]) <- "phiS"
  names(result[["h0"]]) <- "mu"
  return(result)
}

#' Designs a Safe Z Experiment
#'
#' A designed experiment requires (1) a sample size nPlan to plan for, and (2) the parameter of the safe test, i.e.,
#' phiS. Provided with a clinically relevant minimal mean difference meanDiffMin, this function outputs
#' phiS = meanDiffMin as the safe test defining parameter in accordance to the GROW criterion.
#' If a tolerable type II error, i.e., beta, is provided then nPlan can be sampled. The sampled nPlan is then
#' the smallest nPlan for which meanDiffMin can be found with power at least 1 - beta under optional stopping.
#'
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control --independent on n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule e10 > 1/alpha.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both "n"
#' and "phiS". Note that 1-beta defines the power.
#' @param meanDiffMin numeric that defines the minimal relevant mean difference, the smallest population mean
#' that we would like to detect.
#' @param alternative a character string specifying the alternative hypothesis must be one of "twoSided" (default),
#' "greater" or "less".
#' @param nPlan optional numeric vector of length at most 2. When provided, it is used to find the safe test
#' defining parameter phiS. Note that if the purpose is to plan based on nPlan alone, then both meanDiffMin
#' and beta should be set to NULL.
#' @param sigma numeric > 0 representing the assumed population standard deviation used for the test.
#' @param h0 numeric, represents the null hypothesis, default h0=0.
#' @param kappa the true population standard deviation. Default kappa=sigma.
#' @param tol a number that defines the stepsizes between the lowParam and highParam.
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
#' is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param parameter optional test defining parameter. Default set to \code{NULL}.
#' @param nSim integer > 0, the number of simulations needed to compute power or the number of samples paths
#' for the safe z test under continuous monitoring.
#' @param nBoot integer > 0 representing the number of bootstrap samples to assess the accuracy of
#' approximation of the power, the number of samples for the safe z test under continuous monitoring,
#' or for the computation of the logarithm of the implied target.
#' @param pb logical, if \code{TRUE}, then show progress bar.
#' @param grow logical, default set to \code{TRUE} so the grow safe test is used in the design.
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
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
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
#' @references Grunwald, de Heide and Koolen (2019) "Safe Testing" <arXiv:1906.07801>
#' @examples
#' designObj <- designSafeZ(meanDiffMin=0.8, alpha=0.08, beta=0.01, alternative="greater")
#'
#' #nPlan known:
#' designObj <- designSafeZ(nPlan = 100, alpha=0.05)
#'
designSafeZ <- function(meanDiffMin=NULL, beta=NULL, nPlan=NULL,
                        alpha=0.05, h0=0, alternative=c("twoSided", "greater", "less"),
                        sigma=1, kappa=sigma, tol=1e-5,
                        testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, parameter=NULL, nSim=1e3L, nBoot=1e3L,
                        pb=TRUE, grow=TRUE, ...) {

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

  if (!is.null(parameter))
    parameter <- checkAndReturnsEsMinParameterSide("paramToCheck"=parameter, "esMinName"="phiS",
                                                   "alternative"=alternative)

  if (!is.null(meanDiffMin))
    meanDiffMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=meanDiffMin, "esMinName"="meanDiffMin",
                                                     "alternative"=alternative)

  paired <- if (testType=="paired") TRUE else FALSE

  designScenario <- NULL
  note <- NULL

  nPlanBatch <- nPlanTwoSe <- NULL
  nMean <- nMeanTwoSe <- NULL

  logImpliedTarget <- logImpliedTargetTwoSe <- NULL
  betaTwoSe <- NULL

  bootObjN1Plan <- bootObjN1Mean <- bootObjBeta <- bootObjLogImpliedTarget <- NULL

  tempResult <- list()

  if (!is.null(meanDiffMin) && !is.null(beta) && is.null(nPlan)) {
    #scenario 1a: delta + power known, calculate nPlan
    designScenario <- "1a"
    phiS <- meanDiffMin

    tempResult <- computeNPlanSafeZ("meanDiffMin"=meanDiffMin, "beta"=beta, "alpha"=alpha, "alternative"=alternative,
                                    "sigma"=sigma, "kappa"=kappa, "ratio"=ratio, "nSim"=nSim,
                                    "nBoot"=nBoot, "pb"=pb, "parameter"=phiS, "testType"=testType)

    nPlanBatch <- tempResult[["nPlanBatch"]]

    bootObjN1Plan <- tempResult[["bootObjN1Plan"]]
    bootObjN1Mean <- tempResult[["bootObjN1Mean"]]

    if (testType=="oneSample") {
      nPlan <- tempResult[["n1Plan"]]
      names(nPlan) <- "nPlan"
      nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]

      nMean <- tempResult[["n1Mean"]]
      names(nMean) <- "nMean"
      nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]

      note <- paste0("If it is only possible to look at the data once, ",
                     "then nPlan = ", nPlanBatch, ".")
    } else if (testType=="paired") {
      nPlan <- c(tempResult[["n1Plan"]], tempResult[["n1Plan"]])
      names(nPlan) <- c("n1Plan", "n2Plan")

      nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
      nPlanTwoSe <- c(nPlanTwoSe, nPlanTwoSe)

      nMean <- c(tempResult[["n1Mean"]], tempResult[["n1Mean"]])
      names(nMean) <- c("n1Mean", "n2Mean")
      nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]
      nMeanTwoSe <- c(nMeanTwoSe, nMeanTwoSe)

      note <- paste0("If it is only possible to look at the data once, ",
                     "then n1Plan = ", nPlanBatch[1], " and n2Plan = ",
                     nPlanBatch[2], ".")
    } else if (testType=="twoSample") {
      nPlan <- c(tempResult[["n1Plan"]], ceiling(ratio*tempResult[["n1Plan"]]))
      names(nPlan) <- c("n1Plan", "n2Plan")
      nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
      nPlanTwoSe <- c(nPlanTwoSe, ratio*nPlanTwoSe)


      nMean <- c(tempResult[["n1Mean"]], ceiling(ratio*tempResult[["n1Mean"]]))
      names(nMean) <- c("n1Mean", "n2Mean")
      nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]
      nMeanTwoSe <- c(nMeanTwoSe, ratio*nMeanTwoSe)

      note <- paste0("If it is only possible to look at the data once, ",
                     "then n1Plan = ", nPlanBatch[1], " and n2Plan = ",
                     nPlanBatch[2], ".")
    }
  } else if (!is.null(meanDiffMin) && is.null(beta) && is.null(nPlan)) {
    designScenario <- "1b"
    phiS <- meanDiffMin

    nPlan <- NULL
    beta <- NULL

    meanDiffMin <- meanDiffMin
  } else if (is.null(meanDiffMin) && is.null(beta) && !is.null(nPlan)) {
    #scenario 1c: only nPlan known, can perform a pilot (no warning though)
    designScenario <- "1c"

    return(designPilotSafeZ("nPlan"=nPlan, "alpha"=alpha, "alternative"=alternative,
                            "sigma"=sigma, "kappa"=kappa, "tol"=tol, "paired"=paired, "parameter"=parameter))
  } else if (!is.null(meanDiffMin) && is.null(beta) && !is.null(nPlan)) {
    # scenario 2: given effect size and nPlan, calculate power and implied target
    designScenario <- "2"

    nPlan <- checkAndReturnsNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

    tempResult <- computeBetaSafeZ("meanDiffMin"=meanDiffMin, "nPlan"=nPlan, "alpha"=alpha, "sigma"=sigma, "kappa"=kappa,
                                   "alternative"=alternative, "testType"=testType, "parameter"=parameter)

    beta <- tempResult[["beta"]]
    bootObjBeta <- tempResult[["bootObjBeta"]]
    betaTwoSe <- 2*bootObjBeta[["bootSe"]]

    logImpliedTarget <- tempResult[["logImpliedTarget"]]
    bootObjLogImpliedTarget <- tempResult[["bootObjLogImpliedTarget"]]
    logImpliedTargetTwoSe <- 2*bootObjLogImpliedTarget[["bootSe"]]

    phiS <- if (is.null(parameter)) meanDiffMin else parameter

  } else if (is.null(meanDiffMin) && !is.null(beta) && !is.null(nPlan)) {
    designScenario <- "3"

    meanDiffMin <- tryOrFailWithNA(
      computeMinEsBatchSafeZ("nPlan"=nPlan, "alpha"=alpha, "beta"=beta, "sigma"=sigma,
                             "kappa"=kappa, "alternative"=alternative, "testType"=testType,
                             "parameter"=parameter)
    )
    phiS <- if (is.null(parameter)) meanDiffMin else parameter

    if (is.na(meanDiffMin))
      note <- "No meanDiffMin found for the provided beta and nPlan"
    else
      note <- "The reported meanDiffMin is based on the batch analysis."
  }

  if (is.null(designScenario)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1.a) non-null meanDiffMin, non-null beta and NULL nPlan, or \n",
         "(1.b) non-null meanDiffMin, NULL beta, and NULL nPlan, or \n",
         "(1.c) NULL meanDiffMin, NULL beta, non-null nPlan, or \n",
         "(2) non-null meanDiffMin, NULL beta and non-null nPlan, or \n",
         "(3) NULL meanDiffMin, non-null beta, and non-null nPlan.")
  }

  if (is.na(meanDiffMin))
    meanDiffMin <- NULL

  if (!is.null(meanDiffMin))
    names(meanDiffMin) <- "mean difference"

  if (designScenario %in% 2:3) {
    n2Plan <- nPlan[2]

    names(nPlan) <- if (is.na(n2Plan)) "n1Plan" else c("n1Plan", "n2Plan")
  }

  result <- list("parameter"=phiS, "esMin"=meanDiffMin, "alpha"=alpha, "alternative"=alternative,
                 "h0"=h0, "testType"=testType, "paired"=paired, "sigma"=sigma, "kappa"=kappa,
                 "ratio"=ratio, "pilot"=FALSE,
                 "nPlan"=nPlan, "nPlanTwoSe"=nPlanTwoSe, "nPlanBatch"=nPlanBatch,
                 "nMean"=nMean, "nMeanTwoSe"=nMeanTwoSe,
                 "beta"=beta, "betaTwoSe"=betaTwoSe,
                 "logImpliedTarget"=logImpliedTarget, "logImpliedTargetTwoSe"=logImpliedTargetTwoSe,
                 "bootObjN1Plan"=bootObjN1Plan, "bootObjBeta"=bootObjBeta,
                 "bootObjLogImpliedTarget"=bootObjLogImpliedTarget, "bootObjN1Mean"=bootObjN1Mean,
                 "call"=sys.call(), "timeStamp"=Sys.time(), "note"=note)
  class(result) <- "safeDesign"

  names(result[["h0"]]) <- "mu"
  names(result[["parameter"]]) <- "phiS"
  return(result)
}



# Batch design fnts ------

#' Helper function: Computes the planned sample size based on the minimal clinical relevant mean
#' difference, alpha and beta.
#'
#' @inheritParams  designSafeZ
#' @param highN integer that defines the largest n of our search space for n. This might be the
#' largest n that we are able to fund.
#'
#' @return a list which contains at least nPlan and the phiS, that is, the parameter that defines
#' the safe test.
computeNPlanBatchSafeZ <- function(meanDiffMin, alpha=0.05, beta=0.2, sigma=1, kappa=sigma,
                                   alternative=c("twoSided", "greater", "less"),
                                   testType=c("oneSample", "paired", "twoSample"),
                                   tol=1e-5, highN=1e6, ratio=1, parameter=NULL,
                                   grow=TRUE) {
  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  result <- list(nPlan=NULL, phiS=NULL)
  meanDiffMin <- abs(meanDiffMin)

  n1Plan <- NULL
  n2Plan <- NULL

  n1OverNEffRatio <- if (testType=="twoSample") (1+ratio)/ratio else 1

  if (grow) {
    phiS <- abs(meanDiffMin)

    if (alternative == "twoSided") {
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
      n1Plan <- ceiling(nEff * n1OverNEffRatio)
      n2Plan <- ceiling(nEff * n1OverNEffRatio * ratio)
    } else {
      n1Plan <- ceiling(nEff)
      n2Plan <- if (testType == "paired") n1Plan else NULL
    }

    phiS <- meanDiffMin
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
      if (testType == "twoSample") {
        n1Plan <- ceiling(nEffExact*n1OverNEffRatio)
        n2Plan <- ceiling(nEffExact*n1OverNEffRatio*ratio)
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
      # twoSided

      nEffExactUpper <- tryOrFailWithNA(
        ((sigma*sqrt(2*log(2/alpha))-kappa*qnorm(beta))/meanDiffMin)^2
      )

      # Note(Alexander): Translate to lower and upper bound in terms of n1
      #
      lowN <- floor(nEffExact*n1OverNEffRatio)
      highN <- ceiling(nEffExactUpper*n1OverNEffRatio)

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
      candidateNEff <- candidateN1/n1OverNEffRatio

      continueWhile <- TRUE

      # Note(Alexander): Loop backwards to find a smaller n1Plan
      #
      while (continueWhile && candidateN1 > lowN) {
        continueWhile <- criterionFunctionExact(n=candidateNEff, parameter=NULL)

        if (isTRUE(continueWhile)) {
          candidateN1 <- candidateN1 - 1
          candidateNEff <- candidateN1/n1OverNEffRatio
        } else {
          candidateN1 <- candidateN1 + 1
          break()
        }
      }

      nEff <- candidateN1/n1OverNEffRatio

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
    } # end twoSided

    if (alternative=="less")
      phiS <- - phiS
  }

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

#' Helper function: Computes the type II error based on the minimal clinically relevant effect size and sample size.
#'
#' @inheritParams designSafeZ
#'
#' @return numeric that represents the type II error
computeBetaBatchSafeZ <- function(meanDiffMin, nPlan, alpha=0.05, sigma=1, kappa=sigma,
                              alternative=c("twoSided", "greater", "less"),
                              testType=c("oneSample", "paired", "twoSample"),
                              parameter=NULL) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  testType <- match.arg(testType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  if (is.null(parameter))
    parameter <- meanDiffMin

  if (alternative=="twoSided") {
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

#' Computes the smallest mean difference that is detectable with chance 1-beta, for the provided
#' sample size
#'
#' @inheritParams designSafeZ
#' @param maxIter maximum number of iterations in the optimisation process for two-sided designs
#'
#' @return numeric > 0 that represents the minimal detectable mean difference
computeMinEsBatchSafeZ <- function(nPlan, alpha=0.05, beta=0.2, sigma=1, kappa=sigma,
                                    alternative=c("twoSided", "greater", "less"),
                                    testType=c("oneSample", "paired", "twoSample"),
                                    parameter=NULL, maxIter=10) {
  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  testType <- match.arg(testType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  if (alternative=="twoSided") {
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

# Sampling functions for design ----

#' Simulate stopping times for the safe z-test
#'
#' @inheritParams designSafeZ
#' @param nMax integer > 0, maximum sample size of the (first) sample in each sample path.
#' @param wantEValuesAtNMax logical. If \code{TRUE} then compute eValues at nMax. Default \code{FALSE}.
#'
#' @return a list with stoppingTimes and breakVector. Entries of breakVector are 0, 1. A 1 represents stopping
#' due to exceeding nMax, and 0 due to 1/alpha threshold crossing, which implies that in corresponding stopping
#' time is Inf.
#'
#' @export
#'
#' @examples
#' sampleStoppingTimesSafeZ(0.7, nSim=10)
sampleStoppingTimesSafeZ <- function(meanDiffMin, alpha=0.05, alternative = c("twoSided", "less", "greater"),
                                     sigma=1, kappa=sigma, nSim=1e3L, nMax=1e3, ratio=1,
                                     testType=c("oneSample", "paired", "twoSample"), parameter=NULL,
                                     wantEValuesAtNMax=FALSE, pb=TRUE) {
  stopifnot(alpha > 0, alpha <= 1, is.finite(nMax))

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  testType <- match.arg(testType)

  ## Object that will be returned. A sample of stopping times
  stoppingTimes <- breakVector <- integer(nSim)
  eValuesStopped <- eValuesAtNMax <- numeric(nSim)

  if (!is.null(parameter)) {
    phiS <- parameter
  } else {
    meanDiffMin <- checkAndReturnsEsMinParameterSide(meanDiffMin, "alternative"=alternative, "esMinName"="meanDiffMin")
    phiS <- meanDiffMin
  }

  if (testType=="twoSample" && length(nMax)==1) {
    nMax <- c(nMax, ceiling(ratio*nMax))
  } else if (testType=="paired" && length(nMax)==2) {
    nMax <- nMax[1]
  }

  if (length(nMax)==2)
    ratio <- nMax[2]/nMax[1]

  if (pb)
    pbSafe <- utils::txtProgressBar(style=3, title="Safe test threshold crossing")

  tempN <- defineTTestN("lowN"=1, "highN"=nMax[1], "ratio"=ratio, "testType"=testType)

  n1Vector <- tempN[["n1"]]
  n2Vector <- tempN[["n2"]]
  nEffVector <- tempN[["nEff"]]

  simData <- generateNormalData("nPlan"=nMax, "nSim"=nSim, "muTrue"=meanDiffMin, "sigmaTrue"=kappa,
                                "paired"=FALSE)

  for (sim in seq_along(stoppingTimes)) {
    if (testType %in% c("oneSample", "paired")) {
      x1 <- simData[["dataGroup1"]][sim, ]
      zVector <- 1/sqrt(n1Vector)*cumsum(x1)

      if (wantEValuesAtNMax) {
        eValuesAtNMax[sim] <- safeZTestStat("z"=zVector[length(zVector)], "phiS"=phiS,
                                            "n1"=nMax[1], n2=NULL,
                                            "alternative"=alternative, "sigma"=sigma)
      }
    } else {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/(n1Vector)*cumsum(x1)
      x1BarVector <- x1BarVector[n1Vector]

      x2 <- simData[["dataGroup2"]][sim, ]
      x2Cumsum <- cumsum(x2)[n2Vector]
      x2BarVector <- 1/n2Vector*x2Cumsum

      zVector <- sqrt(nEffVector)*(x1BarVector - x2BarVector)/sigma

      if (wantEValuesAtNMax) {
        eValuesAtNMax[sim] <- safeZTestStat("z"=zVector[length(zVector)], "phiS"=phiS,
                                            "n1"=nMax[1], n2=nMax[2],
                                            "alternative"=alternative, "sigma"=sigma)
      }
    }

    for (j in seq_along(n1Vector)) {
      evidenceNow <- if (testType %in% c("oneSample", "paired")) {
        safeZTestStat("z"=zVector[j], "phiS"=phiS, "n1"=n1Vector[j], n2=NULL,
                      "alternative"=alternative, "sigma"=sigma)
      } else {
        safeZTestStat("z"=zVector[j], "phiS"=phiS, "n1"=n1Vector[j], "n2"=n2Vector[j],
                      "alternative"=alternative, "sigma"=sigma)
      }

      if (evidenceNow > 1/alpha) {
        stoppingTimes[sim] <- n1Vector[j]
        eValuesStopped[sim] <- evidenceNow
        break()
      }

      # Note(Alexander): If passed maximum nPlan[1] stop.
      #   For power calculations if beyond nPlan[1], then set to Inf, doesn't matter for the quantile
      #
      if (n1Vector[j] >= nMax[1]) {
        stoppingTimes[sim] <- n1Vector[j]
        breakVector[sim] <- 1
        eValuesStopped[sim] <- evidenceNow
        break()
      }
    }

    if (pb)
      utils::setTxtProgressBar(pbSafe, "value"=sim/nSim, "title"="Trials")
  }

  result <- list("stoppingTimes"=stoppingTimes, "breakVector"=breakVector,
                 "eValuesStopped"=eValuesStopped, "eValuesAtNMax"=eValuesAtNMax)
  return(result)
}


#' Helper function: Computes the type II error based on the minimal clinically relevant mean difference and nPlan
#'
#' @inheritParams designSafeZ
#' @inheritParams sampleStoppingTimesSafeZ
#'
#' @return a list which contains at least beta and an adapted bootObject of class  \code{\link[boot]{boot}}.
#' @export
#'
#' @examples
#' computeBetaSafeZ(meanDiffMin=0.7, 20, nSim=10)
computeBetaSafeZ <- function(meanDiffMin, nPlan, alpha=0.05, alternative=c("twoSided", "greater", "less"),
                             sigma=1, kappa=sigma, testType=c("oneSample", "paired", "twoSample"),
                             parameter=NULL, pb=TRUE, nSim=1e3L, nBoot=1e3L) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  testType <- match.arg(testType)

  ratio <- if (length(nPlan) == 2) nPlan[2]/nPlan[1] else 1

  if (testType=="twoSample" && length(nPlan)==1) {
    nPlan <- c(nPlan, nPlan)
    warning('testType=="twoSample" specified, but nPlan[2] not provided. nPlan[2] is set to ratio = ', ratio,
            'times nPlan[1] = ', nPlan[2])
  }

  if (!is.null(parameter)) {
    phiS <- parameter
  } else {
    meanDiffMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=meanDiffMin, "alternative"=alternative,
                                                     "esMinName"="meanDiffMin")
    phiS <- meanDiffMin
  }

  tempResult <- sampleStoppingTimesSafeZ("meanDiffMin"=meanDiffMin, "alpha"=alpha,
                                         "alternative" = alternative, "sigma"=sigma,
                                         "kappa"=kappa, "nSim"=nSim, "nMax"=nPlan,
                                         "ratio"=ratio, "testType"=testType,
                                         "parameter"=phiS, "pb"=pb, "wantEValuesAtNMax"=TRUE)

  times <- tempResult[["stoppingTimes"]]

  # Note(Alexander): Break vector is 1 whenever the sample path did not stop
  breakVector <- tempResult[["breakVector"]]

  # Note(Alexander): Setting the stopping time to Inf for these paths doesn't matter for the quantile
  times[as.logical(breakVector)] <- Inf

  bootObjBeta <- computeBootObj("values"=times, "objType"="beta", "nPlan"=nPlan[1], "nBoot"=nBoot)

  eValuesAtNMax <- tempResult[["eValuesAtNMax"]]

  bootObjLogImpliedTarget <- computeBootObj("values"=eValuesAtNMax, "objType"="logImpliedTarget",
                                            "nBoot"=nBoot)

  result <- list("beta" = bootObjBeta[["t0"]],
                 "bootObjBeta" = bootObjBeta,
                 "logImpliedTarget"=bootObjLogImpliedTarget[["t0"]],
                 "bootObjLogImpliedTarget"=bootObjLogImpliedTarget)

  return(result)
}


#' Helper function: Computes the planned sample size based on the minimal clinical relevant mean
#' difference, alpha and beta
#'
#'
#' @inheritParams designSafeZ
#' @inheritParams sampleStoppingTimesSafeZ
#'
#' @return a list which contains at least nPlan and an adapted bootObject of class  \code{\link[boot]{boot}}.
#'
#' @export
#'
#' @examples
#' computeNPlanSafeZ(0.7, 0.2, nSim=10)
computeNPlanSafeZ <- function(meanDiffMin, beta=0.2, alpha=0.05, alternative = c("twoSided", "less", "greater"),
                              testType=c("oneSample", "paired", "twoSample"), sigma=1, kappa=sigma,
                              ratio=1, nSim=1e3L, nBoot=1e3L, parameter=NULL, pb=TRUE, nMax=1e8) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  testType <- match.arg(testType)

  if (!is.null(parameter)) {
    phiS <- parameter
  } else {
    meanDiffMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=meanDiffMin, "alternative"=alternative,
                                                     "esMinName"="meanDiffMin")
    phiS <- meanDiffMin
  }

  tempObj <- computeNPlanBatchSafeZ("meanDiffMin"=meanDiffMin, "alpha"=alpha,
                                           "beta"=beta, "sigma"=sigma, "kappa"=kappa,
                                           "alternative"=alternative, "testType"=testType,
                                           "parameter"=parameter, "ratio"=ratio)
  nPlanBatch <- tempObj[["nPlan"]]

  samplingResults <- sampleStoppingTimesSafeZ("meanDiffMin"=meanDiffMin, "alpha"=alpha, "alternative" = alternative,
                                              "sigma"=sigma, "kappa"=kappa, "nSim"=nSim, "testType"=testType,
                                              "ratio"=ratio, "nMax"=nPlanBatch, "parameter"=phiS,
                                              "pb"=pb)

  times <- samplingResults[["stoppingTimes"]]

  bootObjN1Plan <- computeBootObj("values"=times, "objType"="nPlan", "beta"=beta, "nBoot"=nBoot)

  n1Plan <- ceiling(bootObjN1Plan[["t0"]])

  bootObjN1Mean <- computeBootObj("values"=times, "objType"="nMean", "nPlan"=n1Plan, "nBoot"=nBoot)

  n1Mean <- ceiling(bootObjN1Mean[["t0"]])

  result <- list("n1Plan" = n1Plan, "bootObjN1Plan" = bootObjN1Plan,
                 "n1Mean"=n1Mean, "bootObjN1Mean"=bootObjN1Mean,
                 "nPlanBatch"=nPlanBatch)

  return(result)
}

# Helpers ------


#' Help function to compute the effective sample size based on a length 2 vector of samples
#'
#' @inheritParams designSafeZ
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
