# Testing fnts -------

#' Computes E-Values Based on the T-Statistic
#'
#' A summary stats version of \code{\link{safeTTest}()} with the data replaced by t, n1 and n2, and the
#' design object by deltaS.
#'
#' @param t numeric that represents the observed t-statistic.
#' @param parameter numeric this defines the safe test S, i.e., a likelihood ratio of t distributions with in
#' the denominator the likelihood with delta = 0 and in the numerator an average likelihood defined by
#' 1/2 time the likelihood at the non-centrality parameter sqrt(nEff)*parameter and 1/2 times the likelihood at
#' the non-centrality parameter -sqrt(nEff)*parameter.
#' @param n1 integer that represents the size in a one-sample t-test, (n2=\code{NULL}). When n2 is not \code{NULL},
#' this specifies the size of the first sample for a two-sample test.
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
#' \code{NULL} it implies that the t-statistic is based on one-sample.
#' @param tDensity Uses the the representation of the safe t-test as the likelihood ratio of t densities.
#' @inherit safeTTest
#'
#' @return Returns a numeric that represent the e10, that is, the e-value in favour of the alternative over the null
#'
#' @export
#'
#' @examples
#' safeTTestStat(t=1, n1=100, 0.4)
#' safeTTestStat(t=3, n1=100, parameter=0.3)
safeTTestStat <- function(t, parameter, n1, n2=NULL, alternative=c("two.sided", "less", "greater"), tDensity=FALSE,
                          paired=FALSE, ...) {
  # TODO(Alexander):
  #   One-sided not as stable as two-sided due to hypergeo::genhypergeo for the odd component
  #   1. Use Kummer's transform again (??)
  #   2. Switch to numerical integration. Boundary case
  #
  # safeTTestStat(t=-3.1878, parameter=0.29, n1=315, alternative="greater")
  # safeTTestStat(t=-3.1879, parameter=0.29, n1=315, alternative="greater")
  # safeTTestStat(t=-3.188, parameter=0.29, n1=315, alternative="greater")
  alternative <- match.arg(alternative)
  deltaS <- parameter

  if (is.null(n2) || is.na(n2) || paired==TRUE) {
    nEff <- n1
    nu <- n1-1
  } else {
    nEff <- (1/n1+1/n2)^(-1)
    nu <- n1+n2-2
  }

  # TODO(Alexander): This is not necessarily correct, since the correct one is
  # (2*(y1-y2))^(-1)*(1+sign(y1-y2)*
  #   (pnorm(deltaS/2, sd=1/sqrt(2))-pnorm(-deltaS/2, sd=1/sqrt(2))) #erf(deltaS/2)
  #
  #
  # Problem: t=0, n=1
  #
  # PERHAPS add warning and assume t is then y2-y1, otherwise t must be Inf
  if (nu == 0) {
    if (t==0)
      return(NA)
    else
      return(1)
  }

  a <- t^2/(nu+t^2)
  expTerm <- exp((a-1)*nEff*deltaS^2/2)

  zeroIndex <- abs(expTerm) < .Machine$double.eps
  result <- vector("numeric", length(expTerm))

  zArg <- (-1)*a*nEff*deltaS^2/2
  zArg <- zArg[!zeroIndex]
  # Note(Alexander): This made the vector shorter. Only there where expTerm is non-zero will we evaluate
  # the hypergeometric functions

  aKummerFunction <- Re(hypergeo::genhypergeo(U=-nu/2, L=1/2, zArg))

  if (alternative=="two.sided") {
    result[!zeroIndex] <- expTerm[!zeroIndex] * aKummerFunction
  } else {
    bKummerFunction <- exp(lgamma(nu/2+1)-lgamma((nu+1)/2))*sqrt(2*nEff)*deltaS*t/sqrt(t^2+nu)[!zeroIndex] *
      Re(hypergeo::genhypergeo(U=(1-nu)/2, L=3/2, zArg))
    result[!zeroIndex] <- expTerm[!zeroIndex]*(aKummerFunction + bKummerFunction)
  }

  if (result < 0) {
    warning("Numerical overflow: eValue close to zero. Ratio of t density employed.")
    result <- safeTTestStatTDensity("t"=t, "parameter"=parameter, "nu"=nu,
                                    "nEff"=nEff, "alternative"=alternative)
  }
  return(result)
}

#' safeTTestStat() based on t-densities
#'
#' This is basically just \code{\link{safeTTestStat}()} - 1/alpha. This function is used for root finding for
#' pilot designs.
#'
#' @inheritParams safeTTest
#' @inherit safeTTestStat
#' @param nu numeric > 0 representing the degrees of freedom.
#' @param nEff numeric > 0 representing the effective sample size in a two-sample problem. For one-sample
#' problems this is equal to the sample size.
#'
#' @return Returns a numeric that represent the e10, that is, the e-value in favour of the
#' alternative over the null.
#'
safeTTestStatTDensity <- function(t, parameter, nu, nEff, alternative=c("two.sided", "less", "greater"),
                                  paired=FALSE, ...) {
  deltaS <- parameter

  if (alternative=="two.sided") {
    logTerm1 <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)
    logTerm2 <- stats::dt(t, df=nu, ncp=-sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)

    result <- exp(logTerm1+logTerm2)/2
  } else {
    result <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS)/stats::dt(t, df=nu, ncp=0)
  }

  if (result < 0) {
    warning("Numerical overflow: E-value is essentially zero")
    return(2^(-25))
  }


  return(result)
}

#' safeTTestStat() Subtracted with 1/alpha.
#'
#' This is basically just \code{\link{safeTTestStat}()} - 1/alpha. This function is used for root finding for
#' pilot designs.
#'
#' @inheritParams safeTTest
#' @inherit safeTTestStat
#'
#' @return Returns a numeric that represent the e10 - 1/alpha, that is, the e-value in favour of the
#' alternative over the null - 1/alpha.
#'
safeTTestStatAlpha <- function(t, parameter, n1, n2=NULL, alpha, alternative="two.sided", tDensity=FALSE) {
  safeTTestStat("t"=t, "parameter"=parameter, "n1"=n1, "n2"=n2, "alternative"=alternative, "tDensity"=tDensity) - 1/alpha
}


#' safeTTestStat(), but as a function of the square of t and subtracted with 1/alpha
#'
#' This is basically just \code{\link{safeTTestStat}()} - 1/alpha. This function is used for root finding for
#' pilot designs.
#'
#' @inheritParams safeTTest
#' @inherit safeTTestStat
#'
#' @return Returns a numeric that represent the e10 - 1/alpha, that is, the e-value in favour of the
#' alternative over the null - 1/alpha.
#'
safeT2TestStatAlpha <- function(t, parameter, n1, n2=NULL, alpha, alternative="two.sided", tDensity=FALSE) {
  safeTTestStat("t"=t, "parameter"=parameter, "n1"=n1, "n2"=n2, "alternative"=alternative, "tDensity"=tDensity) - 1/alpha
}


#' Safe Student's T-Test.
#'
#' A safe t-test adapted from  \code{\link[stats]{t.test}()} to perform one and two sample t-tests on vectors of data.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param designObj an object obtained from \code{\link{designSafeT}()}, or \code{NULL}, when pilot
#' equals  \code{TRUE}.
#' @param paired a logical indicating whether you want a paired t-test.
#' @param varEqual a logical variable indicating whether to treat the two variances as being equal. For
#' the moment, this is always \code{TRUE}.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that
#' the number of samples is exactly as planned.
#' @param alpha numeric > 0 only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then
#' the alpha of the design object is used instead in constructing the decision rule S > 1/alpha.
#' @param alternative a character only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE},
#' then the alternative specified by the design object is used instead.
#' @param ciValue numeric is the ciValue-level of the confidence sequence. Default ciValue=NULL,
#' and ciValue = 1 - alpha
#' @param na.rm a logical value indicating whether \code{NA} values should be stripped before
#' the computation proceeds.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the t-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{eValue}{the realised e-value from the safe test.}
#'   \item{confSeq}{A safe confidence interval for the mean appropriate to the specific alternative
#'   hypothesis.}
#'   \item{estimate}{the estimated mean or difference in means or mean difference depending on whether it a one-
#'   sample test or a two-sample test was conducted.}
#'   \item{stderr}{the standard error of the mean (difference), used as denominator in the t-statistic formula.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeTDesign" obtained from \code{\link{designSafeT}()}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.6, alpha=0.008, alternative="greater",
#' testType="twoSample", ratio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safeTTest(x, y, designObj=designObj)      #0.2959334
#'
#' safeTTest(1:10, y = c(7:20), pilot=TRUE)      # s = 658.69 > 1/alpha
safeTTest <- function(x, y=NULL, designObj=NULL, paired=FALSE, varEqual=TRUE,
                      pilot=FALSE, alpha=NULL, alternative=NULL, ciValue=NULL,
                      na.rm=FALSE, ...) {

  result <- list("statistic"=NULL, "n"=NULL, "eValue"=NULL, "confSeq"=NULL, "estimate"=NULL,
                 "alternative"=NULL, "testType"=NULL, "dataName"=NULL, "stderr"=NULL,
                 "call"=sys.call())
  class(result) <- "safeTest"

  if (is.null(designObj) && !pilot) {
    stop("No design given and not indicated that this is a pilot study. Run design first and provide ",
         "this to safeTTest/safe.t.test, or run safeTTest with pilot=TRUE")
  }

  if (!is.null(designObj)) {
    checkDoubleArgumentsDesignObject(designObj, "alternative"=alternative, "alpha"=alpha)

    if (names(designObj[["parameter"]]) != "deltaS")
      warning("The provided design is not constructed for the z-test,",
              "please use designSafeT() instead. The test results might be invalid.")
  }

  n1 <- length(x)

  if (is.null(y)) {
    testType <- "oneSample"
    n <- nEff <- n1
    n2 <- NULL
    nu <- n-1

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    meanObs <- estimate <- mean(x, "na.rm"=na.rm)
    sdObs <- sd(x, "na.rm"=na.rm)

    names(estimate) <- "mean of x"
    names(n) <- "n1"
  } else {
    n2 <- length(y)

    if (paired) {
      if (n1 != n2)
        stop("Data error: Error in complete.cases(x, y): Paired analysis requested, ",
             "but the two samples are not of the same size.")

      testType <- "paired"

      nEff <- n1
      nu <- n1-1

      meanObs <- estimate <- mean(x-y, "na.rm"=na.rm)
      sdObs <- sd(x-y, "na.rm"=na.rm)
      names(estimate) <- "mean of the differences"
    } else {
      testType <- "twoSample"

      nu <- n1+n2-2
      nEff <- (1/n1+1/n2)^(-1)

      sPooledSquared <- ((n1-1)*var(x, "na.rm"=na.rm)+(n2-1)*var(y, "na.rm"=na.rm))/nu

      sdObs <- sqrt(sPooledSquared)

      estimate <- c(mean(x, "na.rm"=na.rm), mean(y, "na.rm"=na.rm))
      names(estimate) <- c("mean of x", "mean of y")
      meanObs <- estimate[1]-estimate[2]
    }

    n <- c(n1, n2)
    names(n) <- c("n1", "n2")
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

    designObj <- designPilotSafeT("nPlan"=n, "alpha"=alpha, "alternative"=alternative, "paired"=paired)
  }

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]
  alpha <- designObj[["alpha"]]

  if (is.null(ciValue))
    ciValue <- 1-alpha

  if (ciValue < 0 || ciValue > 1)
    stop("Can't make a confidence sequence with ciValue < 0 or ciValue > 1, or alpha < 0 or alpha > 1")

  tStat <- tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sdObs)

  if (is.na(tStat))
    stop("Data error: Could not compute the t-statistic")

  names(tStat) <- "t"

  eValue <- safeTTestStat("t"=tStat, "parameter"=designObj[["parameter"]], "n1"=n1,
                          "n2"=n2, "alternative"=alternative, "paired"=paired)

  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  result[["statistic"]] <- tStat
  # result[["parameter"]] <- designObj[["deltaS"]]
  result[["estimate"]] <- estimate
  result[["stderr"]] <- sdObs/nEff
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["testType"]] <- testType
  result[["n"]] <- n
  result[["ciValue"]] <- ciValue

  result[["confSeq"]] <- computeConfidenceIntervalT("meanObs"=meanObs-h0, "sdObs"=sdObs,
                                                    "nEff"=nEff, "nu"=nu,
                                                    "deltaS"=designObj[["parameter"]],
                                                    "ciValue"=ciValue)

  result[["eValue"]] <- eValue

  return(result)
}


#' Alias for safeTTest()
#'
#' @rdname safeTTest
#'
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. For the moment,
#' this is always \code{TRUE}.
#'
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.6, alpha=0.008, alternative="greater",
#' testType="twoSample", ratio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safe.t.test(x, y, alternative="greater", designObj=designObj)      #0.2959334
#'
#' safe.t.test(1:10, y = c(7:20), pilot=TRUE)      # s = 658.69 > 1/alpha
safe.t.test <- function(x, y=NULL, designObj=NULL, paired=FALSE, var.equal=TRUE,
                        pilot=FALSE, alpha=NULL, alternative=NULL, ...) {
  result <- safeTTest("x"=x, "y"=y, "designObj"=designObj, "paired"=paired, "varEqual"=var.equal,
                      "pilot"=pilot, "alpha"=alpha, "alternative"=alternative, ...)
  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  result[["dataName"]] <- dataName
  return(result)
}


#' Helper function: Computes the safe confidence sequence for the mean in a t-test
#'
#' @param nEff numeric > 0, the effective sample size. For one sample test this is just n.
#' @param nu numeric > 0, the degrees of freedom.
#' @param meanObs numeric, the observed mean. For two sample tests this is difference of the means.
#' @param sdObs numeric, the observed standard deviation. For a two-sample test this is the root
#' of the pooled variance.
#' @param deltaS numeric > 0, the safe test defining parameter.
#' @param ciValue numeric is the ciValue-level of the confidence sequence. Default ciValue=0.95.
#' @param g numeric > 0, used as the variance of the normal prior on the population delta
#' Default is \code{NULL} in which case g=delta^2.
#'
#' @return numeric vector that contains the upper and lower bound of the safe confidence sequence
#' @export
#'
#' @examples
#' computeConfidenceIntervalT(meanObs=0.3, sdObs=2, nEff=12, nu=11, deltaS=0.4)
computeConfidenceIntervalT <- function(meanObs, sdObs, nEff, nu, deltaS,
                                       ciValue=0.95, g=NULL) {
  if (is.null(g)) g <- deltaS^2

  trivialConfInt <- c(-Inf, Inf)

  if (nu <= 0) return(trivialConfInt)

  alpha <- 1-ciValue

  numeratorW <- nu*(((1+nEff*g)/alpha^2)^(1/(nu+1))-1)
  denominatorW <- 1-((1+nEff*g)/alpha^2)^(1/(nu+1))/(1+nEff*g)

  W <- numeratorW/denominatorW

  if (W < 0) return(trivialConfInt)

  shift <- sdObs/sqrt(nEff)*sqrt(W)

  lowerCs <- meanObs - shift
  upperCs <- meanObs + shift

  return(unname(c(lowerCs, upperCs)))
}


# Design fnts -------

#' Design a Frequentist T-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I and type II error for the frequentist t-test.
#'
#' @inheritParams designSafeT
#'
#' @return Returns an object of class 'freqTDesign'. An object of class 'freqTDesign' is a list containing at least the
#' following components:
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{esMin}{the minimal clinically relevant standardised effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{lowN}{the smallest n of the search space for n provided by the user.}
#'   \item{highN}{the largest n of the search space for n provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#' }
#' @export
#'
#' @examples
#' designFreqT(0.5)
designFreqT <- function(deltaMin, alpha=0.05, beta=0.2, alternative=c("two.sided", "greater", "less"),
                        h0=0, testType=c("oneSample", "paired", "twoSample"), ...) {
  stopifnot(alpha > 0, beta > 0, alpha < 1, beta < 1)

  testType <- match.arg(testType)
  alternative <- match.arg(alternative)

  alternativeStat <- switch(alternative,
                            "greater"="one.sided",
                            "less"="one.sided",
                            "two.sided"="two.sided")

  testTypeStat <- switch(testType,
                         "twoSample"="two.sample",
                         "oneSample"="one.sample",
                         "paired"="paired")

  tempResult <- stats::power.t.test("delta"=deltaMin, "power"=1-beta, "type"=testTypeStat,
                                    "alternative"=alternativeStat)

  n1Plan <- ceiling(tempResult[["n"]])
  n2Plan <- NULL

  if (testType!="oneSample") n2Plan <- n1Plan

  if (is.null(n2Plan)) {
    nPlan <- n1Plan
    names(nPlan) <- "n1Plan"
  } else {
    nPlan <- c(n1Plan, n2Plan)
    names(nPlan) <- c("n1Plan", "n2Plan")
  }

  result <- list(nPlan=nPlan, "esMin"=deltaMin, "alpha"=alpha, "beta"=beta,
                 "testType"=testType, "alternative"=alternative, "ratio"=1, "h0"=h0)
  class(result) <- "freqTDesign"
  return(result)
}

#' Designs a Safe Experiment to Test Means with a T Test
#'
#' A designed experiment requires (1) a sample size nPlan to plan for, and (2) the parameter of the safe test, i.e.,
#' deltaS. If nPlan is provided, then only the safe test defining parameter deltaS needs to determined. That resulting
#' deltaS leads to an (approximately) most powerful safe test. Typically, nPlan is unknown and the user has to specify
#' (i) a tolerable type II error beta, and (ii) a clinically relevant minimal population standardised effect size
#' deltaMin. The procedure finds the smallest nPlan for which deltaMin is found with power of at least 1 - beta.
#'
#' @param deltaMin numeric that defines the minimal relevant standardised effect size, the smallest effect size that
#' we would the experiment to be able to detect.
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control --independent of n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule e10 > 1/alpha.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both
#' the sample sizes and deltaS, which defines the test. Note that 1-beta defines the power.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less".
#' @param nPlan vector of max length 2 representing the planned sample sizes.
#' @param h0 a number indicating the hypothesised true value of the mean under the null. For the moment h0=0.
#' @param lowN integer minimal sample size of the (first) sample when computing the power due to
#' optional stopping. Default lowN is set 1.
#' @param highN integer minimal sample size of the (first) sample when computing the power due to
#' optional stopping. Default highN is set 1e6.
#' @param lowParam numeric defining the smallest delta of the search space for the test-defining deltaS
#' for scenario 3. Currently not yet in use.
#' @param highParam numeric defining the largest delta of the search space for the test-defining deltaS
#' for scenario 3. Currently not yet in use.
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
#' @param seed integer, seed number.
#' @param pb logical, if \code{TRUE}, then show progress bar.
#' @param ... further arguments to be passed to or from methods, but mainly to perform do.calls.
#'
#' @return Returns an object of class 'safeDesign'. An object of class 'safeDesign' is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{parameter}{the safe test defining parameter. Here deltaS.}
#'   \item{esMin}{the minimal clinically relevant standardised effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on
#'   whether it was a one-sample or a two-sample test.}
#'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
#'   ratio between the planned randomisation of condition 2 over condition 1.}
#'   \item{lowN}{the smallest n of the search space for n provided by the user.}
#'   \item{highN}{the largest n of the search space for n provided by the user.}
#'   \item{lowParam}{the smallest delta of the search space for delta provided by the user.}
#'   \item{highParam}{the largest delta of the search space for delta provided by the user.}
#'   \item{tol}{the step size between lowParam and highParam provided by the user.}
#'   \item{pilot}{\code{FALSE} (default) specified by the user to indicate that the design is not a pilot study.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.8, alpha=0.03, alternative="greater")
#' designObj
#'
#' # "Scenario 1.a": Minimal clinically relevant standarised mean difference and tolerable type
#' # II error also known. Goal: find nPlan.
#' designObj <- designSafeT(deltaMin=0.8, alpha=0.03, beta=0.4, nSim=10, alternative="greater")
#' designObj
#'
#' # "Scenario 2": Minimal clinically relevant standarised mean difference and nPlan known.
#' # Goal: find the power, hence, the type II error of the procedure under optional stopping.
#'
#' designObj <- designSafeT(deltaMin=0.8, alpha=0.03, nPlan=16, nSim=10, alternative="greater")
#' designObj
designSafeT <- function(deltaMin=NULL, beta=NULL, nPlan=NULL, alpha=0.05, h0=0,
                        alternative=c("two.sided", "greater", "less"),
                        lowN=3L, highN=1e6L, lowParam=0.01, highParam=1.5, tol=0.01,
                        testType=c("oneSample", "paired", "twoSample"), ratio=1,
                        nSim=1e3L, nBoot=1e3L, parameter=NULL, pb=TRUE, seed=NULL, ...) {
  stopifnot(alpha > 0, alpha < 1)

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  if (!is.null(parameter))
    parameter <- checkAndReturnsEsMinParameterSide("paramToCheck"=parameter, "esMinName"="deltaS",
                                                   "alternative"=alternative)

  if (!is.null(deltaMin))
    deltaMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=deltaMin, "esMinName"="deltaMin",
                                                  "alternative"=alternative)

  paired <- if (testType=="paired") TRUE else FALSE

  designScenario <- NULL
  note <- NULL

  nPlanBatch <- nPlanTwoSe <- NULL
  logImpliedTarget <- logImpliedTargetTwoSe <- NULL
  betaTwoSe <- NULL

  bootObjN1Plan <- bootObjBeta <- bootObjLogImpliedTarget <- NULL

  tempResult <- list()

  names(h0) <- "mu"

  if (!is.null(deltaMin) && !is.null(beta) && is.null(nPlan)) {
    #scenario 1a: delta + power known, calculate nPlan
    designScenario <- "1a"
    deltaS <- deltaMin

    tempResult <- computeNPlanSafeT("deltaMin"=deltaMin, "beta"=beta, "alpha"=alpha, "alternative"=alternative,
                                    "testType"=testType, "lowN"=lowN, "highN"=highN, "ratio"=ratio,
                                    "seed"=seed, "nSim"=nSim, "nBoot"=nBoot, "parameter"=deltaS, "pb"=pb)

    nPlanBatch <- tempResult[["nPlanBatch"]]

    bootObjN1Plan <- tempResult[["bootObjN1Plan"]]

    if (testType=="oneSample") {
      nPlan <- tempResult[["n1Plan"]]
      names(nPlan) <- "nPlan"
      nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
      note <- paste0("If it is only possible to look at the data once, ",
                     "then nPlan = ", nPlanBatch, ".")
    } else if (testType=="paired") {
      nPlan <- c(tempResult[["n1Plan"]], tempResult[["n1Plan"]])
      names(nPlan) <- c("n1Plan", "n2Plan")

      nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
      nPlanTwoSe <- c(nPlanTwoSe, nPlanTwoSe)
      note <- paste0("If it is only possible to look at the data once, ",
                     "then n1Plan = ", nPlanBatch[1], " and n2Plan = ",
                     nPlanBatch[2], ".")
    } else if (testType=="twoSample") {
      nPlan <- c(tempResult[["n1Plan"]], ceiling(ratio*tempResult[["n1Plan"]]))
      names(nPlan) <- c("n1Plan", "n2Plan")

      nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
      nPlanTwoSe <- c(nPlanTwoSe, ratio*nPlanTwoSe)

      note <- paste0("If it is only possible to look at the data once, ",
                     "then n1Plan = ", nPlanBatch[1], " and n2Plan = ",
                     nPlanBatch[2], ".")
    }
  } else if (!is.null(deltaMin) && is.null(beta) && is.null(nPlan)) {
    designScenario <- "1b"
    deltaS <- deltaMin

    nPlan <- NULL
    beta <- NULL
    deltaMin <- deltaMin
  } else if (is.null(deltaMin) && is.null(beta) && !is.null(nPlan)) {
    #scenario 1c: only nPlan known, can perform a pilot (no warning though)
    designScenario <- "1c"

    return(designPilotSafeT("nPlan"=nPlan, "alpha"=alpha, "h0"=h0, "alternative"=alternative,
                            "lowParam"=lowParam, "tol"=tol, "logging"=TRUE,
                            "paired"=paired))
  } else if (!is.null(deltaMin) && is.null(beta) && !is.null(nPlan)) {
    # scenario 2: given effect size and nPlan, calculate power and implied target
    designScenario <- "2"
    deltaS <- deltaMin

    nPlan <- checkAndReturnsNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

    tempResult <- computeBetaSafeT("deltaMin"=deltaMin, "nPlan"=nPlan, "alpha"=alpha,
                                   "alternative"=alternative, "testType"=testType, "parameter"=deltaS,
                                   "pb"=pb, "nSim"=nSim, "nBoot"=nBoot, "seed"=seed)

    beta <- tempResult[["beta"]]
    bootObjBeta <- tempResult[["bootObjBeta"]]
    betaTwoSe <- 2*bootObjBeta[["bootSe"]]

    logImpliedTarget <- tempResult[["logImpliedTarget"]]
    bootObjLogImpliedTarget <- tempResult[["bootObjLogImpliedTarget"]]
    logImpliedTargetTwoSe <- 2*bootObjLogImpliedTarget[["bootSe"]]
  } else if (is.null(deltaMin) && !is.null(beta) && !is.null(nPlan)) {
    designScenario <- "3"
    stop("Not yet implemented")
  }

  if (is.null(designScenario)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1.a) non-null deltaMin, non-null beta and NULL nPlan, or \n",
         "(1.b) non-null deltaMin, NULL beta, and NULL nPlan, or \n",
         "(1.c) NULL deltaMin, NULL beta, non-null nPlan, or \n",
         "(2) non-null deltaMin, NULL beta and non-null nPlan, or \n",
         "(3) NULL deltaMin, non-null beta, and non-null nPlan.")
  }

  if (is.na(deltaMin))
    deltaMin <- NULL

  if (designScenario %in% 2:3) {
    n2Plan <- nPlan[2]

    names(nPlan) <- if (is.na(n2Plan)) "n1Plan" else c("n1Plan", "n2Plan")
  }

  result <- list("nPlan"=nPlan, "nPlanTwoSe"=nPlanTwoSe, "parameter"=deltaS, "esMin"=deltaMin,
                 "alpha"=alpha, "beta"=beta, "betaTwoSe"=betaTwoSe, "alternative"=alternative,
                 "h0"=h0, "testType"=testType, "paired"=paired,
                 "nPlanBatch"=nPlanBatch, "ratio"=ratio, "pilot"=FALSE,
                 "logImpliedTarget"=logImpliedTarget, "logImpliedTargetTwoSe"=logImpliedTargetTwoSe,
                 "call"=sys.call(), "timeStamp"=Sys.time(),
                 "note"=note, "bootObjN1Plan"=bootObjN1Plan, "bootObjBeta"=bootObjBeta,
                 "bootObjLogImpliedTarget"=bootObjLogImpliedTarget)
  class(result) <- "safeDesign"

  names(result[["esMin"]]) <- "standardised mean difference"
  names(result[["h0"]]) <- "mu"
  names(result[["parameter"]]) <- "deltaS"
  return(result)
}


#' Designs a Safe T-Test Based on Planned Samples nPlan
#'
#' Designs a safe experiment for a prespecified tolerable type I error based on planned sample
#' size(s), which are fixed ahead of time. Outputs a list that includes the deltaS, i.e., the
#' safe test defining parameter.
#'
#' @inheritParams designSafeT
#' @param nPlan the planned sample size(s).
#' @param inverseMethod logical, always \code{TRUE} for the moment.
#' @param paired logical, if \code{TRUE} then paired t-test.
#' @param logging logical, if \code{TRUE}, then add invSToTThresh to output.
#' @param maxIter numeric > 0, the maximum number of iterations of adjustment to the candidate set from
#' lowParam to highParam, if the minimum is not found.
#' @inheritParams replicateTTests
#' @inherit designSafeT
#'
#' @export
#'
#' @examples
#' designPilotSafeT(nPlan=30)
designPilotSafeT <- function(nPlan=50, alpha=0.05, alternative=c("two.sided", "greater", "less"),
                             h0=0, lowParam=0.01, highParam=1.2, tol=0.01, inverseMethod=TRUE,
                             logging=FALSE, paired=FALSE, maxIter=10) {
  # TODO(Alexander): Check relation with forward method, that is, the least conservative test and maximally powered
  # Perhaps trade-off? "inverseMethod" refers to solving minimum of deltaS \mapsto S_{deltaS}^{-1}(1/alpha)
  #
  #
  # TODO(Alexander):
  #       - Add warning when deltaS == lowParam
  #       - Better: Characterise deltaS as a funciton of alpha and n using asymptotic approximation of 1F1
  #
  # TODO(Alexander): Add some bounds on L, or do a presearch
  #
  # stopifnot(n >= 3)
  #
  #     Trick bound the estimation error using Chebyshev: Note do this on log (safeTTestStat)
  #

  alternative <- match.arg(alternative)
  stopifnot(all(nPlan > 0))

  n1 <- nPlan[1]
  n2 <- NULL
  ratio <- 1

  if (length(nPlan)==1) {
    if (isTRUE(paired)) {
      warning("Paired designed specified, but n2 not provided. n2 is set to n1")
      n2 <- n1
      nPlan <- c(n1, n2)
      names(nPlan) <- c("n1Plan", "n2Plan")
      testType <- "paired"
    } else {
      names(nPlan) <- "n1Plan"
      testType <- "oneSample"
    }
  } else if (length(nPlan)==2) {
    n2 <- nPlan[2]
    ratio <- n2/n1

    if (paired) {
      if (n1 != n2) {
        stop("Paired design specified, but nPlan[1] not equal nPlan[2]")
      }
      testType <- "paired"
    } else {
      testType <- "twoSample"
    }
    names(nPlan) <- c("n1Plan", "n2Plan")
  }

  names(h0) <- "mu"

  result <- list("nPlan"=nPlan, "parameter"=NULL, "esMin"=NULL, "alpha"=alpha, "beta"=NULL,
                 "alternative"=alternative, "testType"=testType, "paired"=paired,
                 "h0"=h0, "sigma"=sigma, "kappa"=kappa, "testType"=testType,
                 "ratio"=ratio, "lowParam"=NULL, "highParam"=NULL,
                 "pilot"=FALSE, "call"=sys.call(), "timeStamp"=Sys.time())
  class(result) <- "safeDesign"

  #
  #
  # result <- list("n1Plan"=n1, "n2Plan"=n2, "deltaS"=NA,
  #                "deltaMin"=NULL, "alpha"=alpha, "beta"=NULL,
  #                "lowParam"=lowParam, "highParam"=highParam, "tol"=tol,
  #                "lowN"=NULL, "highN"=NULL, "alternative"=alternative, "testType"=testType,
  #                "ratio"=ratio, "pilot"=TRUE, "call"=sys.call())
  #
  # class(result) <- "safeTDesign"

  if (inverseMethod) {
    invSafeTTestStatAlpha <- function(x) {
      stats::uniroot("f"=safeTTestStatAlpha, "lower"=0, "upper"=3e10, "tol"=1e-9, "parameter"=x,
                     "n1"=n1, "n2"=n2, "alpha"=alpha, "alternative"=alternative)
    }

    # Note(Alexander): tryOrFAilWithNA fixes the problem that there's a possibility
    #   that deltaS too small, the function values are then both negative or both positive.
    #
    invSafeTTestStatAlphaRoot <- function(x){tryOrFailWithNA(invSafeTTestStatAlpha(x)$root)}

    candidateDeltas <- seq(lowParam, highParam, by=tol)

    invSToTThresh <- rep(NA, length(candidateDeltas))

    iter <- 1

    while (all(is.na(invSToTThresh)) && iter <= maxIter) {
      invSToTThresh <- purrr::map_dbl(candidateDeltas, invSafeTTestStatAlphaRoot)

      if (all(is.na(invSToTThresh))) {
        candidateDeltas <- seq(highParam, "length.out"=2*length(candidateDeltas), "by"=tol)
        lowParam <- highParam
        highParam <- candidateDeltas[length(candidateDeltas)]
        iter <- iter + 1
      }
    }

    mPIndex <- which(invSToTThresh==min(invSToTThresh, na.rm=TRUE))

    if (iter > 1) {
      result[["lowParam"]] <- lowParam
      result[["highParam"]] <- highParam
    }

    if (mPIndex==length(candidateDeltas)) {
      # Note(Alexander): Check that mPIndex is not the last one.
      errorMsg <- "The test defining deltaS is equal to highParam. Rerun with do.call on the output object"
      lowParam <- highParam
      highParam <- (length(candidateDeltas)-1)*tol+lowParam
      result[["lowParam"]] <- lowParam
      result[["highParam"]] <- highParam
    } else if (mPIndex==1) {
      errorMsg <- "The test defining deltaS is equal to lowParam. Rerun with do.call on the output object"
      highParam <- lowParam
      lowParam <- highParam-(length(candidateDeltas)-1)*tol
      result[["lowParam"]] <- lowParam
      result[["highParam"]] <- highParam
    }

    if (isTRUE(logging))
      result[["invSToTThresh"]] <- invSToTThresh

    deltaS <- candidateDeltas[mPIndex]

    if (alternative=="less")
      deltaS <- -deltaS

    result[["parameter"]] <- deltaS
    result[["error"]] <- invSafeTTestStatAlpha(candidateDeltas[mPIndex])$estim.prec
  } else {
    # TODO(Alexander): Check relation with forward method, that is, the least conservative test and maximally powered
    # Perhaps trade-off? "inverseMethod" refers to solving minimum of deltaS \mapsto S_{deltaS}^{-1}(1/alpha)
  }
  # TODO(Alexander): By some monotonicity can we only look at the largest or the smallest?
  #
  # designFreqT(deltaMin=deltaMin, alpha=alpha, beta=beta, lowN=lowN, highN=highN)

  names(result[["parameter"]]) <- "deltaS"
  return(result)
}

# Batch design fnts ------

#' Helper function: Computes the planned sample size for the safe t-test based on the minimal clinically
#' relevant standardised effect size, alpha and beta.
#'
#' @inheritParams  designSafeT
#'
#' @return a list which contains at least nPlan and the phiS the parameter that defines the safe test
computeNPlanBatchSafeT <- function(deltaMin, alpha=0.05, beta=0.2,
                                   alternative=c("two.sided", "greater", "less"),
                                   testType=c("oneSample", "paired", "twoSample"),
                                   lowN=3, highN=1e6, ratio=1) {
  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  deltaMin <- abs(deltaMin)
  deltaTrue <- deltaMin
  deltaS <- deltaMin

  result <- list(nPlan=NULL, "deltaS"=deltaS)

  n1Plan <- NULL
  n2Plan <- NULL

  nDef <- defineTTestN("lowN"=lowN, "highN"=highN, "ratio"=ratio, "testType"=testType)

  n1Vector <- nDef[["n1"]]
  n2Vector <- nDef[["n2"]]
  nuVector <- nDef[["nu"]]
  nEffVector <- nDef[["nEff"]]

  if (alternative=="two.sided") {
    for (i in seq_along(n1Vector)) {
      qBeta <- sqrt(stats::qf("p"=beta, "df1"=1, "df2"=nuVector[i], "ncp"=nEffVector[i]*deltaTrue^2))

      eValue <- safeTTestStat("t"=qBeta, "parameter"=deltaMin, "n1"=n1Vector[i], "n2"=n2Vector[i])

      if (eValue > 1/alpha) {
        n1Plan <- n1Vector[i]
        n2Plan <- n2Vector[i]
        break()
      }
    }
  } else {
    for (i in seq_along(n1Vector)) {
      qBeta <- stats::qt("p"=beta, "df"=nuVector[i], "ncp"=sqrt(nEffVector[i])*deltaTrue)

      eValue <- safeTTestStat("t"=qBeta, "parameter"=deltaMin, "n1"=n1Vector[i],
                              "n2"=n2Vector[i], "alternative"="greater")

      if (eValue > 1/alpha) {
        n1Plan <- n1Vector[i]
        n2Plan <- n2Vector[i]
        break()
      }
    }
  }

  if (is.null(n1Plan)) {
    stop("Could not compute a batch sample size. Increase lowN and highN, which are now",
         lowN, " and ", highN, " respectively.")
  }

  if (alternative=="less")
    result[["deltaS"]] <- -deltaS

  if (testType=="paired")
    n2Plan <- n1Plan

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  return(result)
}

#' Helper function: Computes the minimal clinically relevant standardised mean difference for the safe t-test
#' nPlan and beta.
#'
#' @inheritParams  designSafeT
#'
#' @return a list which contains at least nPlan and the phiS the parameter that defines the safe test
computeEsMinSafeT <- function(nPlan, alpha=0.05, beta=0.2,
                              alternative=c("two.sided", "greater", "less"),
                              testType=c("oneSample", "paired", "twoSample"),
                              lowN=3, highN=1e6, ratio=1) {
  stop("Not yet implemented")

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  deltaMin <- abs(deltaMin)
  deltaTrue <- deltaMin
  deltaS <- deltaMin

  result <- list(nPlan=NULL, "deltaS"=deltaS)

  n1Plan <- NULL
  n2Plan <- NULL

  nDef <- defineTTestN("lowN"=lowN, "highN"=highN, "ratio"=ratio, "testType"=testType)

  n1Vector <- nDef[["n1"]]
  n2Vector <- nDef[["n2"]]
  nuVector <- nDef[["nu"]]
  nEffVector <- nDef[["nEff"]]

  if (alternative=="two.sided") {
    for (i in seq_along(n1Vector)) {
      qBeta <- sqrt(stats::qf("p"=beta, "df1"=1, "df2"=nuVector[i], "ncp"=nEffVector[i]*deltaTrue^2))

      eValue <- safeTTestStat("t"=qBeta, "parameter"=deltaMin, "n1"=n1Vector[i], "n2"=n2Vector[i])

      if (eValue > 1/alpha) {
        n1Plan <- n1Vector[i]
        n2Plan <- n2Vector[i]
        break()
      }
    }
  } else {
    for (i in seq_along(n1Vector)) {
      qBeta <- stats::qt("p"=beta, "df"=nuVector[i], "ncp"=sqrt(nEffVector[i])*deltaTrue)

      eValue <- safeTTestStat("t"=qBeta, "parameter"=deltaMin, "n1"=n1Vector[i],
                              "n2"=n2Vector[i], "alternative"="greater")

      if (eValue > 1/alpha) {
        n1Plan <- n1Vector[i]
        n2Plan <- n2Vector[i]
        break()
      }
    }
  }

  if (is.null(n1Plan)) {
    stop("Could not compute a batch sample size. Increase lowN and highN, which are now",
         lowN, " and ", highN, " respectively.")
  }

  if (alternative=="less")
    result[["deltaS"]] <- -deltaS

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  return(result)
}

# Sampling functions for design ----

#' Simulate stopping times for the safe z-test
#'
#' @inheritParams designSafeT
#' @param deltaTrue numeric, the value of the true standardised effect size (test-relevant parameter).
#' This argument is used by `designSafeT()` with `deltaTrue <- deltaMin`
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
#' sampleStoppingTimesSafeT(0.7, nSim=10)
sampleStoppingTimesSafeT <- function(deltaTrue, alpha=0.05, alternative = c("two.sided", "less", "greater"),
                                     testType=c("oneSample", "paired", "twoSample"),
                                     nSim=1e3L, nMax=1e3, ratio=1, #designObj=NULL,
                                     lowN=3L, parameter=NULL, seed=NULL,
                                     wantEValuesAtNMax=FALSE, pb=TRUE) {
  stopifnot(alpha > 0, alpha <= 1, is.finite(nMax))

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  ## Object that will be returned. A sample of stopping times
  stoppingTimes <- breakVector <- integer(nSim)
  eValuesStopped <- eValuesAtNMax <- numeric(nSim)

  # if (!is.null(designObj)) {
  #   parameter <- designObj[["parameter"]]
  #   nMax <- designObj[["nPlanBatch"]]
  #   alternative <- designObj[["alternative"]]
  #   testType <- designObj[["testType"]]
  # }

  if (!is.null(parameter)) {
    deltaS <- parameter
  } else {
    deltaTrue <- checkAndReturnsEsMinParameterSide(deltaTrue, "alternative"=alternative, "esMinName"="deltaTrue")
    deltaS <- deltaTrue
  }

  if (testType=="twoSample" && length(nMax)==1) {
    nMax <- c(nMax, ceiling(ratio*nMax))
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
    pbSafe <- utils::txtProgressBar(style=3, title="Safe test threshold crossing")

  tempN <- defineTTestN("lowN"=1, "highN"=nMax[1], "ratio"=ratio, "testType"=testType)

  n1Vector <- tempN[["n1"]]
  n2Vector <- tempN[["n2"]]
  nEffVector <- tempN[["nEff"]]

  simData <- generateNormalData("nPlan"=nMax, "nSim"=nSim, "deltaTrue"=deltaTrue,
                                "sigmaTrue"=1, "paired"=FALSE, "seed"=seed)

  for (sim in seq_along(stoppingTimes)) {
    if (testType %in% c("oneSample", "paired")) {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/n1Vector*cumsum(x1)
      x1SquareVector <- cumsum(x1^2)
      sX1Vector <- sqrt(1/(n1Vector-1)*(x1SquareVector - n1Vector*x1BarVector^2))

      badIndeces <- which(n1Vector-1 <= 0)
      sX1Vector[badIndeces] <- 1

      tValues <- sqrt(nEffVector)*x1BarVector/sX1Vector

      if (wantEValuesAtNMax) {
        eValuesAtNMax[sim] <- safeTTestStat("t"=tValues[length(tValues)], "parameter"=deltaS,
                                            "n1"=n1Max, "n2"=n2Max,
                                            "alternative"=alternative, "paired"=FALSE)
      }
    } else {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/(n1Vector)*cumsum(x1)
      x1BarVector <- x1BarVector[n1Vector]
      x1SquareVector <- cumsum(x1^2)[n1Vector]

      x2 <- simData[["dataGroup2"]][sim, ]
      x2CumSum <- cumsum(x2)[n2Vector]
      x2BarVector <- 1/(n2Vector)*x2CumSum
      x2SquareVector <- cumsum(x2^2)[n2Vector]

      sPVector <- sqrt(1/(n1Vector+n2Vector-2)*
                         (x1SquareVector-n1Vector*x1BarVector^2 + x2SquareVector - n2Vector*x2BarVector^2))

      badIndeces <- which(n1Vector+n2Vector-2 <= 0)
      sPVector[badIndeces] <- 1

      tValues <- sqrt(nEffVector)*(x1BarVector-x2BarVector)/sPVector

      if (wantEValuesAtNMax) {

        eValuesAtNMax[sim] <- safeTTestStat("t"=tValues[length(tValues)], "parameter"=deltaS,
                                            "n1"=nMax[1], n2=nMax[2],
                                            "alternative"=alternative)
      }
    }

    for (j in seq_along(n1Vector)) {
      evidenceNow <- if (testType %in% c("oneSample", "paired")) {
        safeTTestStat("t"=tValues[j], "parameter"=deltaS,
                      "n1"=n1Vector[j], n2=n2Vector[j],
                      "alternative"=alternative, "paired"=FALSE)
      } else {
        safeTTestStat("t"=tValues[j], "parameter"=deltaS,
                      "n1"=n1Vector[j], n2=n2Vector[j],
                      "alternative"=alternative)
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

  if (pb)
    close(pbSafe)

  result <- list("stoppingTimes"=stoppingTimes, "breakVector"=breakVector,
                 "eValuesStopped"=eValuesStopped, "eValuesAtNMax"=eValuesAtNMax)
  return(result)
}


#' Helper function: Computes the type II error of the safeTTest based on the minimal clinically relevant
#' standardised mean difference and nPlan.
#'
#' @inheritParams designSafeT
#' @inheritParams sampleStoppingTimesSafeT
#'
#' @return a list which contains at least beta and an adapted bootObject of class
#' \code{\link[boot]{boot}}.
#' @export
#'
#' @examples
#' computeBetaSafeT(deltaMin=0.7, 27, nSim=10)
computeBetaSafeT <- function(deltaMin, nPlan, alpha=0.05, alternative=c("two.sided", "greater", "less"),
                             testType=c("oneSample", "paired", "twoSample"), seed=NULL,
                             parameter=NULL, pb=TRUE, nSim=1e3L, nBoot=1e3L) {

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan) == 2) nPlan[2]/nPlan[1] else 1

  if (testType=="twoSample" && length(nPlan)==1) {
    nPlan <- c(nPlan, nPlan)
    warning('testType=="twoSample" specified, but nPlan[2] not provided. nPlan[2] is set to ratio = ', ratio,
            'times nPlan[1] = ', nPlan[2])
  }

  if (!is.null(parameter)) {
    deltaS <- parameter
  } else {
    deltaMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=deltaMin, "alternative"=alternative,
                                                  "esMinName"="deltaMin")
    deltaS <- deltaMin
  }

  tempResult <- sampleStoppingTimesSafeT("deltaTrue"=deltaMin, "alpha"=alpha,
                                         "alternative" = alternative,
                                         "nSim"=nSim, "nMax"=nPlan, "ratio"=ratio,
                                         "testType"=testType, "parameter"=deltaS,
                                         "pb"=pb, "wantEValuesAtNMax"=TRUE, "seed"=seed)

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


#' Helper function: Computes the planned sample size of the safe t-test based on the
#' minimal clinical relevant standardised mean difference.
#'
#'
#' @inheritParams designSafeT
#' @inheritParams sampleStoppingTimesSafeT
#'
#' @return a list which contains at least nPlan and an adapted bootObject of class  \code{\link[boot]{boot}}.
#'
#' @export
#'
#' @examples
#' computeNPlanSafeT(0.7, 0.2, nSim=10)
computeNPlanSafeT <- function(deltaMin, beta=0.2, alpha=0.05, alternative = c("two.sided", "less", "greater"),
                              testType=c("oneSample", "paired", "twoSample"), lowN=3,
                              highN=1e6, ratio=1, nSim=1e3L, nBoot=1e3L, parameter=NULL, pb=TRUE,
                              nMax=1e6, seed=NULL) {

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  if (!is.null(parameter)) {
    deltaS <- parameter
  } else {
    deltaMin <- checkAndReturnsEsMinParameterSide("paramToCheck"=deltaMin, "alternative"=alternative,
                                                  "esMinName"="deltaMin")
    deltaS <- deltaMin
  }

  tempObj <- computeNPlanBatchSafeT("deltaMin"=deltaMin, "alpha"=alpha, "beta"=beta,
                                    "alternative"=alternative, "testType"=testType, "lowN"=lowN,
                                    "highN"=highN, "ratio"=ratio)
  nPlanBatch <- tempObj[["nPlan"]]

  samplingResults <- sampleStoppingTimesSafeT("deltaTrue"=deltaMin, "alpha"=alpha,
                                              "alternative" = alternative, "seed"=seed,
                                              "nSim"=nSim, "nMax"=nPlanBatch, "ratio"=ratio,
                                              "testType"=testType, "parameter"=deltaS, "pb"=pb)

  times <- samplingResults[["stoppingTimes"]]

  bootObjN1Plan <- computeBootObj("values"=times, "objType"="nPlan", "beta"=beta, "nBoot"=nBoot)

  result <- list("n1Plan" = ceiling(bootObjN1Plan[["t0"]]),
                 "bootObjN1Plan" = bootObjN1Plan, "nPlanBatch"=nPlanBatch)

  return(result)
}

# Helper fnts ------

#' Computes a Sequence of (Effective) Sample Sizes
#'
#' Helper function that outputs the sample sizes, effective sample sizes and the degrees of
#' freedom depending on the type of t-test. Also used for z-tests.
#'
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#'
#' @return Returns the sample sizes and degrees of freedom.
defineTTestN <- function(lowN=3, highN=100, ratio=1,
                         testType=c("oneSample", "paired", "twoSample")) {
  testType <- match.arg(testType)

  if (testType %in% c("twoSample")) {
    n1 <- lowN:highN
    n2 <- ceiling(ratio*n1)
    nEff <- ratio/(1+ratio)*n1
    nu <- (1+ratio)*n1-2
  } else if (testType %in% c("oneSample", "paired")) {
    n1 <- lowN:highN
    n2 <- NULL
    nEff <- n1
    nu <- nEff-1
  }
  result <- list("n1"=n1, "n2"=n2, "nEff"=nEff, "nu"=nu)
  return(result)
}

# Data generating fnt ------

#' Generates Normally Distributed Data Depending on the Design
#'
#' The designs supported are "oneSample", "paired", "twoSample".
#'
#' @inheritParams replicateTTests
#' @param muTrue numeric representing the true mean for simulations with a z-test.
#' Default \code{NULL}
#'
#' @return Returns a list of two data matrices contains at least the following components:
#'
#' \describe{
#'   \item{dataGroup1}{a matrix of data dimension nSim by \code{nPlan[1]}.}
#'   \item{dataGroup2}{a matrix of data dimension nSim by \code{nPlan[2]}.}
#' }
#' @export
#'
#' @examples
#' generateNormalData(20, 15, deltaTrue=0.3)
generateNormalData <- function(nPlan, nSim=1000L, deltaTrue=NULL, muGlobal=0, sigmaTrue=1, paired=FALSE,
                               seed=NULL, muTrue=NULL) {
  stopifnot(all(nPlan > 0))

  if ((is.null(deltaTrue) && is.null(muTrue)) || !is.null(deltaTrue) && !is.null(muTrue))
    stop("Please provide either deltaTrue (t-test), or muTrue (z-test).")

  result <- list("dataGroup1"=NULL, "dataGroup2"=NULL)
  set.seed(seed)

  # TODO(Alexander): vector("mode"="list", length=length(nPlan))

  n1Plan <- nPlan[1]

  if (is.null(muTrue))
    muTrue <- deltaTrue*sigmaTrue

  if (length(nPlan)==1) {
    dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muTrue, "sd"=sigmaTrue)
    dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
    dataGroup2 <- NULL
  } else {
    n2Plan <- nPlan[2]

    if (paired) {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muGlobal + muTrue/sqrt(2), "sd"=sigmaTrue)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nSim, "mean"=muGlobal - muTrue/sqrt(2), "sd"=sigmaTrue)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nSim)
    } else {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muGlobal + muTrue/2, "sd"=sigmaTrue)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nSim, "mean"=muGlobal - muTrue/2, "sd"=sigmaTrue)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nSim)
    }
  }

  return(list("dataGroup1"=dataGroup1, "dataGroup2"=dataGroup2))
}

# Vignette helper fnts ------

#' Simulate Early Stopping Experiments for the T Test
#'
#' Applied to a 'safeDesign' object this function empirically shows the performance of safe experiments under
#' optional stopping.
#'
#' @param object A safeDesign obtained obtained from \code{\link{designSafeT}()}.
#' @param nSim integer, number of iterations.
#' @param nsim integer, formally the number of iterations, but by default nsim=nSim
#' @param seed integer, seed number.
#' @param deltaTrue numeric, if NULL, then the minimally clinically relevant standardised effect size is used
#' as the true data generating effect size deltaTrue.
#' @inherit replicateTTests
#'
#' @import stats
#' @export
#'
#' @examples
#' # Design safe test
#' alpha <- 0.05
#' beta <- 0.20
#' deltaMin <- 1
#' designObj <- designSafeT(deltaMin, alpha=alpha, beta=beta, nSim=100)
#'
#' # Design frequentist test
#' freqObj <- designFreqT(deltaMin, alpha=alpha, beta=beta)
#'
#' # Simulate based on deltaTrue=deltaMin
#' simResultsDeltaTrueIsDeltaMin <- simulate(object=designObj, nSim=100)
#'
#' # Simulate based on deltaTrue > deltaMin
#' simResultsDeltaTrueIsLargerThanDeltaMin <- simulate(
#'   object=designObj, nSim=100, deltaTrue=2)
#'
#' # Simulate under the null deltaTrue = 0
#' simResultsDeltaTrueIsNull <- simulate(
#'   object=designObj, nSim=100, deltaTrue=0)
#'
#' simulate(object=designObj, deltraTrue=0, nSim=100, freqOptioStop=TRUE,
#'          nPlanFreq=freqObj$nPlan)
simulate.safeDesign <- function(object, nsim=nSim, seed=NULL, deltaTrue=NULL, muGlobal=0, sigmaTrue=1, lowN=3,
                                safeOptioStop=TRUE, freqOptioStop=FALSE, nPlanFreq=NULL,
                                logging=TRUE, pb=TRUE, nSim=1, ...) {

  if (object[["pilot"]])
    stop("No simulation for unplanned pilot designs")

  if (object[["testType"]] %in% c("oneSample", "paired", "twoSample")) {
    if (is.null(deltaTrue))
      deltaTrue <- object[["parameter"]]

    paired <- if (object[["testType"]]=="paired") TRUE else FALSE

    result <- replicateTTests("nPlan"=object[["nPlan"]], "deltaTrue"=deltaTrue,
                              "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired,
                              "alternative"=object[["alternative"]], "lowN"=lowN, "nSim"=nsim,
                              "alpha"=object[["alpha"]], "beta"=object[["beta"]],
                              "safeOptioStop"=safeOptioStop, "parameter"=object[["parameter"]],
                              "freqOptioStop"=freqOptioStop, "nPlanFreq"=nPlanFreq,
                              "logging"=logging, "seed"=seed, "pb"=pb, ...)

    object <- utils::modifyList(object, result)
    class(object) <- "safeTSim"
    return(object)
  }
}


#' Simulate Early Stopping Experiments
#'
#' Simulate multiple data sets to show the effects of optional testing for safe (and frequentist) tests.
#'
#' @inheritParams designSafeT
#' @param deltaTrue numeric, the value of the true standardised effect size (test-relevant parameter).
#' @param muGlobal numeric, the true global mean of a paired or two-sample t-test. Its value should not
#' matter for the test. This parameter is treated as a nuisance.
#' @param sigmaTrue numeric > 0,the true standard deviation of the data. Its value should not  matter
#' for the test.This parameter treated is treated as a nuisance.
#' @param lowN integer that defines the smallest n of our search space for n.
#' @param nSim the number of replications, that is, experiments with max samples nPlan.
#' @param safeOptioStop logical, \code{TRUE} implies that optional stopping simulation is performed for
#' the safe test.
#' @param parameter numeric, the safe test defining parameter, i.e., deltaS (use designSafeT to find this).
#' @param freqOptioStop logical, \code{TRUE} implies that optional stopping simulation is performed for
#' the frequentist test.
#' @param nPlanFreq the frequentist sample size(s) to plan for. Acquired from \code{\link{designFreqT}()}.
#' @param paired logical, if \code{TRUE} then paired t-test.
#' @param seed To set the seed for the simulated data.
#' @param logging logical, if \code{TRUE}, then return the simulated data.
#' @param pb logical, if \code{TRUE}, then show progress bar.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class "safeTSim". An object of class "safeTSim" is a list containing at
#' least the following components:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{deltaTrue}{the value of the true standardised effect size (test-relevant parameter) provided by
#'   the user.}
#'   \item{muGlobal}{the true global mean of a paired or two-sample t-test (nuisance parameter) provided by
#'   the user.}
#'   \item{paired}{if \code{TRUE} then paired t-test.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{lowN}{the smallest number of samples (first group) at which monitoring of the tests begins.}
#'   \item{nSim}{the number of replications of the experiment.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{parameter}{the parameter (point prior) used in the safe test derived from the design.
#'   Acquired from \code{\link{designSafeT}()}.}
#'   \item{nPlanFreq}{the frequentist planned sample size(s). Acquired from \code{\link{designFreqT}}()}
#'   \item{safeSim}{list with the simulation results of the safe test under optional stopping.}
#'   \item{freqSim}{list with the simulation results of the frequentist test under optional stopping.}
#'}
#'
#' @export
#'
#' @examples
#'
#' # Design safe test
#' alpha <- 0.05
#' beta <- 0.20
#' designObj <- designSafeT(1, alpha=alpha, beta=beta)
#'
#' # Design frequentist test
#' freqObj <- designFreqT(1, alpha=alpha, beta=beta)
#'
#' # Simulate under the alternative with deltaTrue=deltaMin
#' simResults <- replicateTTests(nPlan=designObj$nPlan, deltaTrue=1, parameter=designObj$parameter,
#'                               nPlanFreq=freqObj$nPlan, beta=beta, nSim=250)
#'
#' # Should be about 1-beta
#' simResults$safeSim$powerAtN1Plan
#'
#' # This is higher due to optional stopping
#' simResults$safeSim$powerOptioStop
#'
#' # Optional stopping allows us to do better than n1PlanFreq once in a while
#' simResults$safeSim$probLeqN1PlanFreq
#' graphics::hist(simResults$safeSim$allN, main="Histogram of stopping times", xlab="n1",
#'                breaks=seq.int(designObj$nPlan[1]))
#'
#' # Simulate under the alternative with deltaTrue > deltaMin
#' simResults <- replicateTTests(nPlan=designObj$nPlan, deltaTrue=1.5, parameter=designObj$parameter,
#'                               nPlanFreq=freqObj$nPlan, beta=beta, nSim=250)
#'
#' # Should be larger than 1-beta
#' simResults$safeSim$powerAtN1Plan
#'
#' # This is even higher due to optional stopping
#' simResults$safeSim$powerOptioStop
#'
#' # Optional stopping allows us to do better than n1PlanFreq once in a while
#' simResults$safeSim$probLeqN1PlanFreq
#' graphics::hist(simResults$safeSim$allN, main="Histogram of stopping times", xlab="n1",
#'                breaks=seq.int(designObj$nPlan[1]))
#'
#' # Under the null deltaTrue=0
#' simResults <- replicateTTests(nPlan=designObj$nPlan, deltaTrue=0, parameter=designObj$parameter,
#'                               nPlanFreq=freqObj$nPlan, freqOptioStop=TRUE, beta=beta, nSim=250)
#'
#' # Should be lower than alpha, because if the null is true, P(S > 1/alpha) < alpha for all n
#' simResults$safeSim$powerAtN1Plan
#'
#' # This is a bit higher due to optional stopping, but if the null is true,
#' # then still P(S > 1/alpha) < alpha for all n
#' simResults$safeSim$powerOptioStop
#'
#' # Should be lowr than alpha, as the experiment is performed as was planned
#' simResults$freqSim$powerAtN1Plan
#'
#' # This is larger than alpha, due to optional stopping.
#' simResults$freqSim$powerOptioStop
#' simResults$freqSim$powerOptioStop > alpha
replicateTTests <- function(nPlan, deltaTrue, muGlobal=0, sigmaTrue=1, paired=FALSE,
                            alternative=c("two.sided", "greater", "less"), lowN=3,
                            nSim=1000L, alpha=0.05, beta=0.2,
                            safeOptioStop=TRUE, parameter=NULL,
                            freqOptioStop=FALSE, nPlanFreq=NULL,
                            logging=TRUE, seed=NULL, pb=TRUE, ...) {

  stopifnot(all(nPlan > lowN), lowN > 0, nSim > 0, alpha > 0, alpha < 1,
            any(safeOptioStop, freqOptioStop))

  alternative <- match.arg(alternative)

  n1Plan <- nPlan[1]
  n2Plan <- NULL

  if (length(nPlan)==2) {
    n2Plan <- nPlan[2]

    testType <- if (paired) "paired" else "twoSample"
  } else {
    if (paired) {
      warning("Paired simulations wanted, but length(nPlan)=1. Duplicate n2Plan=n1Plan")
      n2Plan <- n1Plan
      nPlan <- c(n1Plan, n2Plan)
    }

    testType <- if (paired) "paired" else "oneSample"
  }

  result <- list("nPlan"=nPlan, "deltaTrue"=deltaTrue, "muGlobal"=muGlobal, "paired"=paired,
                 "alternative"=alternative, "lowN"=lowN, "nSim"=nSim, "alpha"=alpha, "beta"=beta, "testType"=testType,
                 "parameter"=parameter, "nPlanFreq"=nPlanFreq, safeSim=list(), freqSim=list())
  class(result) <- "safeTSim"

  if (safeOptioStop) {
    if (is.null(parameter)) {
      stop("To simulate safe t-tests results under optional stopping, this function 'replicateTTests' requires ",
           "the specification of the safe test with a parameter. This parameter can be acquired by running ",
           "the 'designSafeT()' function")
    }

    if (paired && n1Plan != n2Plan)
      stop("For a paired t-test n2Plan needs to equal n1Plan")

    safeSim <- list("powerOptioStop"=NA, "powerAtN1Plan"=NA, "nMean"=NA, "probLeqN1PlanFreq"=NA,
                    "probLessNDesign"=NA, "lowN"=NA)

    allSafeN <- rep(n1Plan, "times"=nSim)
    eValues <- safeDecisionAtN <- allSafeDecisions <- vector("mode"="integer", "length"=nSim)
  }

  if (freqOptioStop) {
    if (!safeOptioStop) {
      if (is.null(nPlanFreq)) {
        warning("No nPlanFreq specified, use nPlan instead.")
        n1PlanFreq <- n1Plan
        n2PlanFreq <- n2Plan
      }
    } else {
      n1PlanFreq <- nPlanFreq[1]
      n2PlanFreq <- NULL

      if (length(nPlanFreq)==2) {
        n2PlanFreq <- nPlanFreq[2]
      } else {
        if (paired) {
          warning("Paired simulations wanted, but length(nPlan)=1. Duplicate n2Plan=n1Plan")
          n2PlanFreq <- n1PlanFreq
          nPlanFreq <- c(n1PlanFreq, n2PlanFreq)
          result[["nPlanFreq"]] <- nPlanFreq
        }
      }
    }

    # Note(Alexander): This means that n1Plan and n2Plan refer to the planned samples of the safe tests

    if (is.null(n1PlanFreq)) {
      stop("To simulate frequentist t-tests results under optional stopping, this ",
           "function 'replicateTTests' requires the specification of n1PlanFreq. To figure out how many ",
           "samples one requires in a frequentist test, please run the 'designFreqT' function.")
    }

    if (!is.null(n2Plan) && is.null(n2PlanFreq)) {
      stop("To simulate a two-sample frequentist t-tests results under optional stopping, this",
           "function 'replicateTTests' requires the specification of n1PlanFreq. To figure out how many ",
           "samples one requires in a frequentist test, please run the 'designFreqT' function.")
    }

    if (paired && n1PlanFreq != n2PlanFreq)
      stop("For a paired t-test n2PlanFreq needs to equal n1PlanFreq")

    freqSim <- list("powerOptioStop"=NA, "powerAtN1Plan"=NA, "nMean"=NA, "probLessNDesign"=NA, "lowN"=NA)

    allFreqN <- rep(n1PlanFreq, "times"=nSim)
    pValues <- freqDecisionAtN <- allFreqDecisions <- vector("mode"="integer", "length"=nSim)
  }

  ratio <- if (is.null(n2Plan) || paired) 1 else n2Plan/n1Plan

  someData <- generateNormalData("nPlan"=c(n1Plan, n2Plan), "nSim"=nSim, "deltaTrue"=deltaTrue,
                                 "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed)

  dataGroup1 <- someData[["dataGroup1"]]
  dataGroup2 <- someData[["dataGroup2"]]

  if (safeOptioStop) {
    n1Samples <- seq.int(lowN, n1Plan)
    n2Samples <- if (is.null(n2Plan)) NULL else ceiling(ratio*n1Samples)

    if (pb)
      pbSafe <- utils::txtProgressBar(style=3, title="Safe optional stopping")

    for (iter in seq.int(nSim)) {
      subData1 <- dataGroup1[iter, ]
      subData2 <- dataGroup2[iter, ]

      someT <- unname(stats::t.test("x"=subData1, "y"=subData2, "alternative"=alternative,
                                    "var.equal"=TRUE, "paired"=paired)[["statistic"]])
      someS <- safeTTestStat("t"=someT, "parameter"=parameter, "n1"=n1Plan, "n2"=n2Plan, "alternative"=alternative,
                             "paired"=paired)

      eValues[iter] <- someS

      if (someS >= 1/alpha)
        safeDecisionAtN[iter] <- 1

      for (k in seq_along(n1Samples)) {

        # TODO(Alexander): Perhaps replace by custom t computing to speed things up
        #
        someT <- unname(stats::t.test("x"=subData1[seq.int(n1Samples[k])], "y"=subData2[seq.int(n2Samples[k])],
                                      "alternative"=alternative, "var.equal"=TRUE, "paired"=paired)[["statistic"]])

        someS <- safeTTestStat("n1"=n1Samples[k], "n2"=n2Samples[k], "t"=someT, "parameter"=parameter,
                               "alternative"=alternative, "paired"=paired)

        if (someS >= 1/alpha) {
          allSafeN[iter] <- n1Samples[k]
          allSafeDecisions[iter] <- 1

          eValues[iter] <- someS
          break()
        }
      } # End loop lowN to n1Plan

      if (pb)
        utils::setTxtProgressBar(pbSafe, value=iter/nSim, title="Experiments")

    } # End iterations

    if (pb)
      close(pbSafe)

    safeSim <- list("powerOptioStop"=mean(allSafeDecisions),
                    "powerAtN1Plan"=mean(safeDecisionAtN),
                    "nMean"=mean(allSafeN),
                    "probLessNDesign"=mean(allSafeN < n1Plan),
                    "lowN"=min(allSafeN), "eValues"=eValues
    )

    safeSim[["allN"]] <- allSafeN
    safeSim[["allSafeDecisions"]] <- allSafeDecisions
    safeSim[["allRejectedN"]] <- allSafeN[-which(allSafeN*allSafeDecisions==0)]

    if (!is.null(nPlanFreq))
      safeSim[["probLeqN1PlanFreq"]] <- mean(allSafeN <= nPlanFreq[1])

    if (isTRUE(logging)) {
      safeSim[["dataGroup1"]] <- dataGroup1
      safeSim[["dataGroup2"]] <- dataGroup2
    }

    result[["safeSim"]] <- safeSim
  }

  if (freqOptioStop) {
    # Note(Alexander): Adjust data set
    #
    if (is.null(n2PlanFreq)) {
      ratio <- 1

      if (n1PlanFreq < n1Plan)
        dataGroup1 <- dataGroup1[, seq.int(n1PlanFreq)]

      if (n1PlanFreq > n1Plan) {
        n1Diff <- n1PlanFreq - n1Plan

        someData <- generateNormalData("nPlan"=c(n1Diff, n2Plan), "nSim"=nSim, "deltaTrue"=deltaTrue,
                                       "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed+1)
        dataGroup1 <- cbind(dataGroup1, someData[["dataGroup1"]])
      }
    } else {
      # Note(Alexander): Two-sample case

      ratio <- if (paired) 1 else n2PlanFreq/n1PlanFreq

      if (n1PlanFreq < n1Plan) {
        dataGroup1 <- dataGroup1[, seq.int(n1PlanFreq)]
      } else if (n1PlanFreq > n1Plan) {
        n1Diff <- n1PlanFreq - n1Plan

        someData <- generateNormalData("nPlan"=c(n1Diff, n2PlanFreq), "nSim"=nSim, "deltaTrue"=deltaTrue,
                                       "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed+1)
        dataGroup1 <- cbind(dataGroup1, someData[["dataGroup1"]])
      }

      if (n2PlanFreq < n2Plan) {
        dataGroup2 <- dataGroup2[, seq.int(n2PlanFreq)]
      } else if (n2PlanFreq > n2Plan) {
        n2Diff <- n2PlanFreq - n2Plan

        someData <- generateNormalData("nPlan"=c(n1PlanFreq, n2Diff), "nSim"=nSim, "deltaTrue"=deltaTrue,
                                       "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed+1)
        dataGroup2 <- cbind(dataGroup2, someData[["dataGroup2"]])
      }
    }

    n1Samples <- seq.int(lowN, n1PlanFreq)

    n2Samples <- if (is.null(n2PlanFreq)) NULL else n2Samples <- ceiling(ratio*n1Samples)

    if (pb)
      pbFreq <- utils::txtProgressBar(style=3, title="Frequentist optional stopping")

    for (iter in seq.int(nSim)) {
      subData1 <- dataGroup1[iter, ]
      subData2 <- dataGroup2[iter, ]
      someP <- stats::t.test("x"=subData1, "y"=subData2, "alternative"=alternative,
                             "var.equal"=TRUE, "paired"=paired)[["p.value"]]

      pValues[iter] <- someP

      if (someP < alpha)
        freqDecisionAtN[iter] <- 1

      for (k in seq_along(n1Samples)) {
        someP <- stats::t.test("x"=subData1[seq.int(n1Samples[k])], "y"=subData2[seq.int(n2Samples[k])],
                               "alternative"=alternative, "var.equal"=TRUE, "paired"=paired)[["p.value"]]

        if (someP < alpha) {
          allFreqN[iter] <- n1Samples[k]
          allFreqDecisions[iter] <- 1
          pValues[iter] <- someP
          break()
        }
      } # End loop lowN to n1Plan

      if (pb)
        utils::setTxtProgressBar(pbFreq, value=iter/nSim, title="Experiments")
    } # End iterations

    if (pb)
      close(pbFreq)

    freqSim <- list("powerOptioStop"=mean(allFreqDecisions),
                    "powerAtN1Plan"=mean(freqDecisionAtN),
                    "nMean"=mean(allFreqN),
                    "allFreqDecisions"=allFreqDecisions,
                    "probLessNDesign"=mean(allFreqN < n1PlanFreq),
                    "lowN"=min(allFreqN), "pValues"=pValues
    )

    freqSim[["allN"]] <- allFreqN

    if (safeOptioStop)
      freqSim[["probLeqNSafe"]] <- mean(allFreqN <= n1Plan)

    if (isTRUE(logging)) {
      freqSim[["dataGroup1"]] <- dataGroup1
      freqSim[["dataGroup2"]] <- dataGroup2
    }

    result[["freqSim"]] <- freqSim
  }
  return(result)
}



#' Plots a 'safeTSim' Object
#'
#' @inheritParams plotHistogramDistributionStoppingTimes
#' @param x a 'safeDesign' object acquired from \code{\link{designSafeT}()}.
#' @param y \code{NULL}.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return a histogram object, and called for its side-effect to plot the histogram.
#'
#' @export
#'
#' @examples
#'# Design safe test
#' alpha <- 0.05
#' beta <- 0.20
#' designObj <- designSafeT(1, alpha=alpha, beta=beta)
#'
#' # Design frequentist test
#' freqObj <- designFreqT(1, alpha=alpha, beta=beta)
#'
#' # Simulate under the alternative with deltaTrue=deltaMin
#' simResults <- simulate(designObj, nSim=100)
#'
#' plot(simResults)
#'
#' plot(simResults, showOnlyNRejected=TRUE)
plot.safeTSim <- function(x, y=NULL, showOnlyNRejected=FALSE, nBin=25, ...) {
  plotHistogramDistributionStoppingTimes(x[["safeSim"]],
                                         "nPlan" = x[["nPlan"]][1],
                                         "deltaTrue" = x[["deltaTrue"]],
                                         "showOnlyNRejected" = showOnlyNRejected, "nBin"=nBin)
}



#' Plots the Sample Sizes Necessary for a Tolerable Alpha and Beta as a Function of deltaMin
#'
#' For given tolerable alpha and beta, (1) the planned sample sizes to using a safe test, (2) the
#' frequentist test, and (3) the average sample size necessary due to optional stopping are plotted
#' as a  function of the minimal clinically relevant standardised effect size deltaMin.
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#' @param nMax numeric, the maximum number of samples one has budget for to collect data.
#' @param lowDeltaMin numeric, lowest value for deltaMin of interest
#' @param highDeltaMin numeric, largest value for deltaMin of interest
#' @param stepDeltaMin numeric, step size between lowDeltaMin and highDeltaMin
#' @param freqPlot logical, if \code{TRUE} plot frequentist sample size profiles.
#'
#' @return Returns a list that contains the planned sample size needed for the frequentist and safe tests as a function
#' of the minimal clinically relevant effect sizes. The returned list contains at least the following components:
#'
#' \describe{
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{maxN}{the largest number of samples provided by the user.}
#'   \item{deltaDomain}{vector of the domain of deltaMin.}
#'   \item{allN1PlanFreq}{vector of the planned sample sizes needed for the frequentist test corresponding to
#'   alpha and beta.}
#'   \item{allN1PlanSafe}{vector of the planned sample sizes needed for the safe test corresponding to alpha
#'   and beta.}
#'   \item{allDeltaS}{vector of safe test defining deltaS.}
#' }
#'
#' @export
#'
#' @examples
#' plotSafeTDesignSampleSizeProfile(nSim=1e2L)
plotSafeTDesignSampleSizeProfile <- function(alpha=0.05, beta=0.2, nMax=100, lowDeltaMin=0.1, highDeltaMin=1,
                                             stepDeltaMin=0.1, testType=c("oneSample", "paired", "twoSample"),
                                             alternative=c("two.sided", "greater", "less"), ratio=1, nSim=1e3L,
                                             nBoot=1e3L, seed=NULL, pb=TRUE, freqPlot=FALSE, ...) {
  stopifnot(lowDeltaMin < highDeltaMin, alpha > 0, beta > 0, alpha < 1, beta < 1)

  # Order from high to low
  deltaDomain <- -seq(-highDeltaMin, -lowDeltaMin, by=stepDeltaMin)
  testType <- match.arg(testType)
  alternative <- match.arg(alternative)

  result <- list("alpha"=alpha, "beta"=beta, "nMax"=nMax, "deltaDomain"=deltaDomain)

  lastDeltaIndex <- length(deltaDomain)

  if (lastDeltaIndex < 1)
    stop("Either nMax or deltaDomain is too small. Please decrease lowDeltaMin, increase highDeltaMin, ",
         "or decrease stepDeltaMin.")

  paired <- if (testType=="paired") TRUE else FALSE

  allN1Freq <- allN1PlanSafe <- allN1PlanBatchSafe <- vector("integer", lastDeltaIndex)

  if (pb)
    pbOptioStop <- utils::txtProgressBar("style"=1)

  if (!is.null(seed))
    seed <- 1:lastDeltaIndex+seed

  # 1. Run simulations  ------------

  for (i in seq_along(deltaDomain)) {
    safeDesignObj <- designSafeT("deltaMin"=deltaDomain[i], "alpha"=alpha, "beta"=beta,
                                 "alternative"=alternative, "testType"=testType, "ratio"=ratio,
                                 "nSim"=nSim, "nBoot"=nBoot, "pb"=pb, "seed"=seed[i])

    if (freqPlot) {
      freqObj <- designFreqT("deltaMin"=deltaDomain[i], "alpha"=alpha, "beta"=beta,
                             "alternative"=alternative, "testType"=testType)
      allN1Freq[i] <- freqObj[["nPlan"]][1]
    }

    allN1PlanBatchSafe[i] <- safeDesignObj[["nPlanBatch"]][1]
    allN1PlanSafe[i] <- safeDesignObj[["nPlan"]][1]

    if (pb)
      utils::setTxtProgressBar("pb"=pbOptioStop, "value"=i/lastDeltaIndex)

    if (safeDesignObj[["nPlan"]][1] >= nMax) {
      lastDeltaIndex <- i
      break()
    }
  }

  if (pb)
    close(pbOptioStop)

  deltaDomain <- deltaDomain[1:lastDeltaIndex]
  allN1PlanBatchSafe <- allN1PlanBatchSafe[1:lastDeltaIndex]
  allN1PlanSafe <- allN1PlanSafe[1:lastDeltaIndex]
  allN1Freq <- allN1Freq[1:lastDeltaIndex]

  maxDeltaDomain <- max(deltaDomain)
  minDeltaDomain <- min(deltaDomain)

  # Store in output
  result[["deltaDomain"]] <- deltaDomain
  result[["allN1PlanSafe"]] <- allN1PlanSafe
  result[["allN1PlanBatchSafe"]] <- allN1PlanBatchSafe
  result[["allN1Freq"]] <- allN1Freq


  # 2. Plot Sim  ------
  oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
  on.exit(graphics::par(oldPar))

  graphics::plot(deltaDomain, allN1PlanBatchSafe, type="l", col="blue", lty=2, lwd=2, xlim=c(minDeltaDomain, maxDeltaDomain),
                 ylab="n1", xlab=expression(delta["min"]),
                 main=bquote(~alpha == ~.(alpha) ~ "and" ~beta== ~.(beta)))
  graphics::abline(h=nMax, col="red", lty=2)
  graphics::lines(deltaDomain, allN1PlanSafe, col="black", lwd=2, lty=1)

  if (freqPlot) {
    graphics::lines(deltaDomain, allN1Freq, col="darkgrey", lwd=2, lty=3)
    legendName <- c("nPlan", "nPlanBatch", "nPlanFreq", "nMax")
    legendCol <- c("black", "blue", "darkgrey", "red")
    legendLty <- c(1, 2, 3, 2)
  } else {
    legendName <- c("nPlan", "nPlanBatch", "nMax")
    legendCol <- c("black", "blue", "red")
    legendLty <- c(1, 2, 2)
  }

  graphics::legend("topright", legend = legendName, col = legendCol, lty=legendLty, bty="n")

  return(result)
}
