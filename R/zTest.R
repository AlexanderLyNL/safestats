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
#' The two-sided safe z-test is based on two point priors on the mean differences. It puts half its mass on
#' phiS and the other half on - phiS. To determine the value of phiS, the test is inverted.
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
#' A safe one and two sample z-tests on vectors of data. The function is modelled after \code{\link[stats]{t.test}}.
#'
#' @aliases safe.z.test
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less".
#' @param designObj an object obtained from \code{\link{designSafeZ}}, or \code{NULL}, when pilot is set to \code{TRUE}.
#' @param h0 a number indicating the hypothesised true value of the mean under the null.
#' @param sigma numeric > 0 representing the assumed population standard deviation used for the test.
#' @param paired a logical indicating whether you want the paired z-test.
#' @param confLevel confidence level of the interval. If alpha is given, then set to 1-alpha. Not yet implemented.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the observed
#' number of samples is exactly as planned.
#' @param alpha numeric representing the tolerable type I error rate. This also serves as a decision rule and it was
#' shown that for safe tests S we have P(S > 1/alpha) < alpha under the null.
#' @param tol numeric > 0, the mesh between consecutive parameter values in the candidate space. Only used when pilot
#' is set to \code{TRUE}.
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
#'   \item{sigma}{the assumed population standard deviation provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
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
#' safeZTest(x, y, alternative="greater", designObj=designObj)      #
#'
#' safeZTest(1:10, y = c(7:20), pilot=TRUE)      # s = 7.7543e+20 > 1/alpha
safeZTest <- function(x, y=NULL, designObj=NULL, alternative=c("two.sided", "less", "greater"),
                      h0=0, sigma=1, paired=FALSE, confLevel=1-alpha, pilot=FALSE,
                      alpha=0.05, tol=1e-05, ...) {

  alternative <- match.arg(alternative)

  names(h0) <- "mu"

  result <- list("statistic"=NULL, "n"=NULL, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
                 "alternative"=alternative, "testType"=NULL, "dataName"=NULL, "h0"=h0, "sigma"=sigma,
                 "call"=sys.call())

  class(result) <- "safeTest"

  if (is.null(designObj) && !pilot) {
    stop("No design given and not indicated that this is a pilot study. Run design first and provide ",
         "this to safeZTest/safe.z.test, or run safeZTest with pilot=TRUE")
  }

  if (is.null(y)) {
    testType <- "oneSample"
    n1 <- length(x)
    n2 <- NULL

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    estimate <- mean(x)
    names(estimate) <- "mean of x"
    zStat <- tryOrFailWithNA(sqrt(n1)*(estimate - h0)/sigma)
  } else {
    n1 <- length(x)
    n2 <- length(y)

    if (paired) {
      if (n1 != n2)
        stop("Data error: Error in complete.cases(x, y) : not all arguments have the same length")

      testType <- "paired"

      estimate <- mean(x-y)
      names(estimate) <- "mean of the differences"
      zStat <- tryOrFailWithNA(sqrt(n1)*estimate/sigma)
    } else {
      testType <- "twoSample"

      nEff <- (1/n1+1/n2)^(-1)
      estimate <- c(mean(x), mean(y))
      names(estimate) <- c("mean of x", "mean of y")
      zStat <- tryOrFailWithNA(sqrt(nEff)*(estimate[1]-estimate[2])/sigma)
    }
  }

  result[["testType"]] <- testType

  if (is.na(zStat))
    stop("Could not compute the z-statistic")

  if (pilot) {
    nPlan <- if (is.null(n2)) n1 else c(n1, n2)

    designObj <- designPilotSafeZ("alpha"=alpha, "nPlan"=nPlan, "alternative"=alternative, "h0"=h0,
                                  "sigma"=sigma, "paired"=paired, "tol"=tol)
    designObj[["pilot"]] <- TRUE
  }

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  if (alternative != designObj[["alternative"]])
    warning('The test is designed for a test with direction "', designObj[["alternative"]],
            '", whereas the requested analysis has direction "', alternative, '"')

  sValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=n1, "n2"=n2,
                          "alternative"=alternative, "paired"=paired)

  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  # TODO(Alexander): Add sCIs

  names(zStat) <- "z"

  result[["statistic"]] <- zStat
  result[["estimate"]] <- estimate
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj

  if (is.null(n2)) {
    result[["n"]] <- n1
    names(result[["n"]]) <- "n1"
  } else {
    result[["n"]] <- c(n1, n2)
    names(result[["n"]]) <- c("n1", "n2")
  }
  result[["sValue"]] <- sValue

  return(result)
}

#' Alias for safeZTest
#'
#' @rdname safeZTest
#'
#' @export
safe.z.test <- function(x, y=NULL, designObj=NULL, alternative=c("two.sided", "less", "greater"),
                        h0=0, sigma=1, paired=FALSE, confLevel=1-alpha, pilot=FALSE, alpha=0.05, ...) {
  result <- safeZTest("x"=x, "y"=y, "designObj"=designObj, "alternative"=alternative,
                      "h0"=h0, "paired"=paired, "sigma"=sigma, "confLevel"=confLevel, "pilot"=pilot, "alpha"=alpha, ...)

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
#'
#' @return returns a freqZDesign object.
#' @export
#'
#' @examples
#' designFreqT(0.5)
designFreqZ <- function(alpha=0.05, beta=0.2, meanDiffMin, alternative=c("two.sided", "greater", "less"),
                        lowN=3L, highN=100L, testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, sigma=1, kappa=sigma, ...) {

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
#' @param nPlan vector of max length 2 representing the planned sample sizes.
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
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on whether it was a one-sample
#'   or a two-sample test.}
#'   \item{sigma}{the assumed population standard deviation used for the test provided by the user.}
#'   \item{kappa}{the true population standard deviation, typically, sigma=kappa.}
#'   \item{ratio}{default is 1, only used when testType=="twoSample" and defines n2=ratio*n1.}
#'   \item{tol}{the step size between parameter values in the candidate space.}
#'   \item{pilot}{logical, specifying whether it's a pilot design.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designPilotSafeZ(0.05, nPlan=30)
designPilotSafeZ <- function(nPlan, alpha=0.05, alternative=c("two.sided", "greater", "less"),
                             h0=0, sigma=1, kappa=sigma, tol=1e-5, paired=FALSE) {

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

  names(h0) <- "mu"

  result <- list("nPlan"=nPlan, "parameter"=NULL, "esMin"=NULL, "alpha"=alpha, "beta"=NULL,
                 "alternative"=alternative, "testType"=testType, "paired"=paired,
                 "h0"=h0, "sigma"=sigma, "kappa"=kappa,
                 "ratio"=ratio, "tol"=tol, "pilot"=FALSE, "call"=sys.call())

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

#' Designs a Safe Z Experiment to Test Means with Population Standard Error Assumed to be Known
#'
#' A designed experiment requires (1) a sample size nPlan to plan for, and (2) the parameter of the safe test, i.e.,
#' phiS. If nPlan is provided, then only the safe test defining parameter phiS needs to determined. That resulting phiS
#' leads to an (approximately) most powerful safe test. Typically, nPlan is unknown and the user has to specify
#' (i) a tolerable type II error beta, and (b) a clinically relevant minimal population mean difference meanDiffMin.
#' The procedure finds the smallest nPlan for which meanDiffMin is found with power of at least 1 - beta.
#'
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control --independent on n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule S10 > 1/alpha.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both "n"
#' and "deltaS". Note that 1-beta defines the power.
#' @param meanDiffMin numeric that defines the minimal relevant mean difference, the smallest population mean
#' that we would like to detect.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less".
#' @param nPlan optional numeric vector of length at most 2. When provided, it is used to find the safe test
#' defining parameter phiS.
#' @param h0 a number indicating the hypothesised true value of the mean under the null.
#' @param sigma numeric > 0 representing the assumed population standard deviation used for the test.
#' @param kappa the true population standard deviation. Default kappa=sigma.
#' @param tol a number that defines the stepsizes between the lowParam and highParam.
#' @param highN integer that stops the search process if the lower bound for the candidate nPlan exceeds this
#' number. Default highN is set 8e9 to correspond to the current population of the earth.
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric representing n2/n1. If is.null(n2) then ratio=1.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a safeDesign object that includes:
#'
#' \describe{
#'   \item{nPlan}{the sample size(s) to plan for. Computed based on beta and esMin, or provided by the user if known.}
#'   \item{parameter}{the safe test defining parameter. Here phiS.}
#'   \item{esMin}{the minimally clinically relevant effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error specified by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" effectively provided by the user.}
#'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on whether it was a one-sample
#'   or a two-sample test.}
#'   \item{sigma}{the assumed population standard deviation used for the test provided by the user.}
#'   \item{kappa}{the true population standard deviation, typically, sigma=kappa.}
#'   \item{ratio}{default is 1, only used when testType=="twoSample" and defines n2=ratio*n1.}
#'   \item{tol}{the step size between parameter values in the candidate space.}
#'   \item{pilot}{logical, specifying whether it's a pilot design.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeZ(meanDiffMin=0.8, alpha=0.08, beta=0.01, alternative="greater")
designSafeZ <- function(meanDiffMin=NULL, alpha=0.05, beta=NULL, nPlan=NULL, alternative=c("two.sided", "greater", "less"),
                        h0=0, sigma=1, kappa=sigma, tol=1e-5, highN=8e9,
                        testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, ...) {

  stopifnot(alpha > 0, alpha < 1, sigma > 0, kappa > 0)

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  paired <- if (testType=="paired") TRUE else FALSE

  computeType <- NULL

  if (!is.null(meanDiffMin) && !is.null(beta) && is.null(nPlan)) {
    computeType <- "1a"
  } else if (!is.null(meanDiffMin) && is.null(beta) && !is.null(nPlan)) {
    computeType <- "2"
  } else if (is.null(meanDiffMin) && !is.null(beta) && !is.null(nPlan)) {
    computeType <- "3"
  } else if (!is.null(meanDiffMin) && is.null(beta) && is.null(nPlan)) {
    computeType <- "1b"
  } else if (is.null(meanDiffMin) && is.null(beta) && !is.null(nPlan)) {
    return(designPilotSafeZ("nPlan"=nPlan, "alpha"=alpha, "alternative"=alternative,
                            "h0"=h0, "sigma"=sigma, "kappa"=kappa, "tol"=tol, "paired"=paired))
  }

  if (is.null(computeType)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1) non-null meanDiffMen, non-null beta and NULL nPlan, or \n",
         "(2) non-null meanDiffMen, NULL beta and non-null meanDiffMen, or \n",
         "(3) NULL meanDiffMen, non-null beta, and non-null nPlan, or \n",
         "(4) non-null meanDiffMen, NULL beta, and NULL nPlan, or \n",
         "(5) NULL meanDiffMen, NULL meanDiffMen, non-null nPlan.")
  }

  result <- list("nPlan"=nPlan, "parameter"=NULL, "esMin"=meanDiffMin, "alpha"=alpha, "beta"=beta,
                 "alternative"=alternative, "testType"=testType, "paired"=paired,
                 "h0"=h0, "sigma"=sigma, "kappa"=kappa,
                 "ratio"=ratio, "pilot"=FALSE, "lowN"=NULL, "highN"=NULL, "call"=sys.call())
  class(result) <- "safeDesign"

  if (computeType=="3")
    stop("Not yet implemented")

  # if (computeType==3) {
  #   if (alternative == "two.sided") {
  #     meanDiffMin <- 1
  #   } else {
  #     meanDiffMin <- 1
  #   }
  # }

  meanDiffMin <- abs(meanDiffMin)

  n1Plan <- NULL
  n2Plan <- NULL

  nEffToN1Ratio <- if (testType=="twoSample") (1+ratio)/ratio else 1

  if (computeType=="1a") {
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
        # n1Plan <- ceiling(nTotalExact*nTotalToN1Ratio)
        # n2Plan <- ceiling(nTotalExact*nTotalToN1Ratio*ratio)
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

      if (alternative=="less")
        phiS <- -phiS

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

      criterionFunctionExact <- criterionFunctionFactory("alpha"=alpha, "beta"=beta, "sigma"=sigma,
                                                         "kappa"=kappa, "meanDiffMin"=meanDiffMin,
                                                         "criterionType"="exact")

      candidateN1 <- highN
      candidateNEff <- candidateN1/nEffToN1Ratio

      continueWhile <- TRUE

      # Note(Alexander): Loop backwards to find a smaller n1Plan
      #
      while (continueWhile && candidateN1 > lowN) {
        continueWhile <- criterionFunctionExact(n=candidateNEff, parameter = NULL)
        # candidateNEff, criterionFunctionExact

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

        if (testType=="pairedSampleT")
          n2Plan <- n1Plan

      }

      # Candidate parameters ---
      phiUmp <- sqrt(2/(sigma^2*nEff)*log(2/alpha))

      chiSqInverseBeta <- qchisq(beta, df=1, ncp=nEff*meanDiffMin^2/kappa^2)
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
  } else if (computeType=="1b") {
    sidedConstant <- if (alternative %in% c("greater", "less")) 1 else 2

    phiS <- meanDiffMin

    nEffExact <- 2*sigma^2*log(sidedConstant/alpha)/phiS^2
    # nTotalExact <- 2*sigma^2*log(sidedConstant/alpha)/phiS^2

    if (testType=="twoSample") {
      n1Plan <- ceiling(nEffExact*(1+ratio)/ratio)
      n2Plan <- ceiling(nEffExact*(1+ratio))
      # n1Plan <- ceiling(nTotalExact*nTotalToN1Ratio)
      # n2Plan <- ceiling(nTotalExact*ratio*nTotalToN1Ratio)
      nEff <- (1/n1Plan + 1/n2Plan)^(-1)
      # nTotal <- n1Plan + n2Plan
    } else {
      # nEff <- nTotalExact
      nEff <- nEffExact
      n1Plan <- ceiling(nEff)
      # n1Plan <-ceiling(nTotalExact)
      # nTotal <- n1Plan

      if (testType=="paired")
        n2Plan <- n1Plan

    }

    if (alternative=="two.sided") {
      lowerTail <- sigma^4/(nEff*kappa^2*phiS^2)*(acosh(exp(nEff*phiS^2/(2*sigma^2))/alpha))^2
      result[["beta"]] <- pchisq(q=lowerTail, df=1, ncp=nEff*meanDiffMin^2/kappa^2)

    } else {
      lowerTail <- sqrt(nEff)*(phiS-2*meanDiffMin)/(2*kappa) -
        sigma^2*log(alpha)/(kappa*sqrt(nEff))*1/phiS
      result[["beta"]] <- pnorm(lowerTail)

      if (alternative=="less")
        phiS <- - phiS

    }
  }

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  names(result[["esMin"]]) <- switch(alternative,
                                     "two.sided"="mean differences at least abs(phi)",
                                     "less"="mean differences smaller than phi",
                                     "greater"="mean differences larger than phi")

  names(result[["h0"]]) <- "mu"

  result[["parameter"]] <- phiS
  names(result[["parameter"]]) <- "phiS"
  return(result)
}

#' Builds a Criterion Function Used to Search for the Parameters of a Safe Z-Test
#'
#' Helper function used in the optimisation process.
#'
#' @inheritParams designSafeZ
#' @param criterionType character string any of "lower", "upper", and "exact".
#'
#' @return returns a function to be optimised.
#'
#' @examples
#' safestats:::criterionFunctionFactory(0.05, 0.2, 1, 1, meanDiffMin=0.5, "exact")
criterionFunctionFactory <- function(alpha, beta, sigma, kappa, meanDiffMin,
                                     criterionType=c("lower", "upper", "exact")) {

  result <- switch(criterionType,
                   "upper" = function(n) {
                     pchisq(kappa^2/sigma^2*2*log(2/alpha), df=1, ncp=n*meanDiffMin^2/kappa^2) <= beta
                   },
                   "lower"= function(n) {
                     pchisq(kappa^2/sigma^2*2*log(1/alpha), df=1, ncp=n*meanDiffMin^2/kappa^2) <= beta
                   },
                   "exact"=function(n, parameter=NULL) {
                     if (is.null(parameter))
                       parameter <- sqrt(2*sigma^2/n*log(2/alpha))

                     zArg <- sigma^4/(kappa^2*n*parameter^2)*(acosh(exp(n*parameter^2/(2*sigma^2))/alpha))^2

                     result <- pchisq(zArg, df=1, ncp=n*meanDiffMin^2/kappa^2) <= beta

                     if (is.null(parameter))
                       result <- !result

                     return(result)
                   }
  )
  return(result)
}
