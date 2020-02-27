#' Safe t-test defined at deltaS based on the t-statistic and the sample sizes
#'
#' @param t numeric that represents the observed t-statistic
#' @param deltaS numeric this defines the safe test S, i.e., a likelihood ratio of t distributions with in the
#' denominator the likelihood with delta = 0 and in the numerator an average likelihood defined by
#' 1/2 time the likelihood at the non-centrality parameter sqrt(nEff)*deltaS and 1/2 times the likelihood at the
#' non-centrality parameter -sqrt(nEff)*deltaS
#' @param n1 integer that represents the size in a one-sample t-test, (n2=NULL). When n2 is not NULL, this specifies
#' the size of the first sample for a two-sample test
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus, NULL it
#' implies that the t-statistic is based on one-sample
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less"
#' @param tDensity Uses the the representation of the safe t-test as the likelihood ratio of t densities
#' @param paired a logical, if TRUE ignores n2, and indicates that a paired t-test is performed
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a safeTest object
#'
#' @export
#'
#' @examples
#' safeTTestStat(t=1, n1=100, 0.4)
#' safeTTestStat(t=3, n1=100, deltaS=0.3)
safeTTestStat <- function(t, deltaS, n1, n2=NULL, alternative=c("two.sided", "less", "greater"), tDensity=FALSE,
                          paired=FALSE, ...) {
  alternative <- match.arg(alternative)

  if (is.null(n2) | paired==TRUE) {
    nEff <- n1
    nu <- n1-1
  } else {
    nEff <- (1/n1+1/n2)^(-1)
    nu <- n1+n2-2
  }

  if (tDensity) {
    if (alternative=="two.sided") {
      logTerm1 <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)
      logTerm2 <- stats::dt(t, df=nu, ncp=-sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)

      result <- exp(logTerm1+logTerm2)/2
    } else {
      result <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS)/stats::dt(t, df=nu, ncp=0)
    }
  } else {
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
  }
  return(result)
}

#' Safe Student's t-test.
#'
#' A safe version of t.test to perform one and two sample t-tests on vectors of data
#'
#' @param x a (non-empty) numeric vector of data values
#' @param y an optional (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less"
#' @param designObj an object from designSafeT, or NULL, when pilot=TRUE
#' @param mu0 a number indicating the hypothesised true value of the mean under the null. For the moment mu0=0
#' @param paired a logical indicating whether you want to paired t-test.
#' @param varEqual a logical variable indicating whether to treat the two variances as being equal. For the moment,
#' this is always TRUE.
#' @param confLevel confidence level of the interval. Not yet implemented
#' @param pilot a logical indicating whether a pilot study is run. If TRUE, it is assumed that the number of samples is
#' exactly as planned.
#' @param alpha numeric representing the tolerable type I error rate. This also serves as a decision rule and it was
#' shown that for safe tests S we have P(S > 1/alpha) < alpha under the null.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a safeTest object
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.6, alpha=0.008, alternative="greater",
#' testType="twoSampleT", sampleSizeRatio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safeTTest(x, y, alternative="greater", designObj=designObj)      #0.2959334
#'
#' safeTTest(1:10, y = c(7:20), pilot=TRUE)      # s = 3121.604 > 1/alpha
safeTTest <- function(x, y=NULL, designObj=NULL, alternative=c("two.sided", "less", "greater"),
                        mu0=0, paired=FALSE, varEqual=TRUE, confLevel=0.95, pilot=FALSE,
                        alpha=0.05, ...) {
  # TODO(Alexander): Generalise mu0 = 0 to other mu0
  alternative <- match.arg(alternative)

  result <- list("statistic"=NULL, "parameter"=NULL, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
                 "mu0"=mu0, "stderr"=NULL, "alternative"=alternative, "testType"=NULL, "dataName"=NULL)
  class(result) <- "safeTResult"

  if (is.null(designObj) && !pilot) {
    stop(paste0("No design given and not indicated that this is a pilot study. Run design first and provide ",
                "this to safeTTest/safe.t.test, or run safeTTest/safe.t.test with pilot=TRUE"))
  }

  freqObject <- try(stats::t.test(x=x, y=y, alternative=alternative, mu=mu0, paired=paired, var.equal=varEqual))
  t <- unname(freqObject[["statistic"]])

  if (isTryError(freqObject))
    stop("Data error: could not compute the t-statistic with t.test: ", freqObject[1])

  if (is.null(y)) {
    n1 <- length(x)
    n2 <- NULL
  } else {
    n1 <- length(x)
    n2 <- length(y)
  }

  if (pilot)
    designObj <- designPilotSafeT("n1"=n1, "n2"=n2, "alpha"=alpha, "alternative"=alternative, "paired"=paired)

  if (designObj[["testType"]]=="oneSampleT") {
    if (!is.null(y)) {
      warning(paste0("The analysis is run on a two-sample or paired sample test, but the design object given is",
                     "made for a one-sample t-test"))
    }
  } else if (designObj[["testType"]]=="pairedSampleT") {
    if (!paired) {
      warning(paste0("The analysis is run on a non-paired two-sample t-test, but the design object given is made for",
                     "a paired sample t-test"))
    }

    if (is.null(y)) {
      warning(paste0("The analysis is run on a one-sample t-test, but the design object given is made for a",
                     "paired sample t-test"))
    }
  } else if (designObj[["testType"]]=="pairedSampleT") {
    if (is.null(y)) {
      warning(paste0("The analysis is run on a one-sample t-test, but the design object given is made for a",
                     "two-sample t-test"))
    }
  }
  #
  # TODO(Alexander): Save result, perhaps save freqObject
  #
  sValue <- safeTTestStat("t"=t, "deltaS"=designObj[["deltaS"]], "n1"=n1, "n2"=n2, "alternative"=alternative, "paired"=paired)

  if (is.null(y)) {
    dataName <- as.character(sys.call())[2]
  } else {
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])
  }

  result[["statistic"]] <- t
  result[["parameter"]] <- freqObject[["parameter"]]
  result[["estimate"]] <- freqObject[["estimate"]]
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["n1"]] <- n1
  result[["n2"]] <- n2
  result[["sValue"]] <- sValue

  return(result)
}

#' Alias for \code{\link{safeTTest}}
#'
#' @inheritParams safeTTest
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. For the moment,
#' this is always TRUE.
#' @param conf.level confidence level of the interval. Not yet implemented
#'
#' @export
safe.t.test <- function(x, y=NULL, designObj=NULL, alternative=c("two.sided", "less", "greater"),
                        mu0=0, paired=FALSE, var.equal=TRUE, conf.level=0.95, pilot=FALSE,
                        alpha=0.05, ...) {
  result <- safeTTest("x"=x, "y"=y, "alternative"=alternative, "designObj"=designObj, "mu0"=mu0, "paired"=paired,
                      "varEqual"=var.equal, "confLevel"=conf.level, "pilot"=pilot, "alpha"=alpha, ...)

  if (is.null(y)) {
    dataName <- as.character(sys.call())[2]
  } else {
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])
  }

  result[["dataName"]] <- dataName
  return(result)
}

#' Computes the number of samples necessary to reach a tolerable type I and type II error for the frequentist t-test
#'
#' @inheritParams designSafeT
#' @return returns a freqDesign object
#' @export
#'
#' @examples
#' designFreqT(0.5)
designFreqT <- function(deltaMin, alpha=0.05, beta=0.2, alternative=c("two.sided", "greater", "less"),
                        lowN=3L, highN=100L, testType=c("oneSampleT", "pairedSampleT", "twoSampleT"),
                        sampleSizeRatio=1, ...) {

  stopifnot(lowN >= 2, highN > lowN, alpha > 0, beta >0)

  testType <- match.arg(testType)
  alternative <- match.arg(alternative)

  result <- list("n1PlanFreq"=NA, "n2PlanFreq"=NULL, "deltaMin"=deltaMin, "alpha"=alpha, "beta"=beta,
                 "lowN"=lowN, "highN"=highN, "testType"=testType, "alternative"=alternative)
  class(result) <- "freqDesign"

  if (deltaMin < 0 && alternative=="greater")
    warning("deltaMin < 0, but in the calculations abs(deltaMin) is used instead.")

  # TODO(Alexander): Also need a warning for deltaMin > 0 and alternative=="less" ?

  deltaMin <- abs(deltaMin)

  if (alternative=="two.sided") {
    threshold <- 1-alpha/2
  } else if (alternative %in% c("greater", "less")) {
    threshold <- 1-alpha
  }


  for (n in seq.int(lowN, highN)) {
    if (testType=="twoSampleT") {
      powerT <- stats::pt(stats::qt(threshold, df=((1+sampleSizeRatio)*n-2), ncp=0), df=(1+sampleSizeRatio)*n-2,
                   ncp=sqrt(sampleSizeRatio/(1+sampleSizeRatio)*n)*deltaMin, lower.tail=FALSE)
    } else {
      powerT <- stats::pt(stats::qt(threshold, df=(n-1), ncp=0), df=(n-1), ncp=sqrt(n)*deltaMin, lower.tail=FALSE)
    }

    if (powerT >= (1-beta)) {
      result[["n1PlanFreq"]] <- n

      if (testType=="twoSampleT")
        result[["n2PlanFreq"]] <- ceiling(sampleSizeRatio*n)

      if (testType=="pairedSampleT")
        result[["n2PlanFreq"]] <- n

      return(result)
      #
      break()
    }
  }
  return(result)
}

#' Basically just safeTTestStat - 1/alpha
#'
#' @inheritParams safeTTestStat
#' @inheritParams safeTTest
#'
safeTTestStatAlpha <- function(t, deltaS, n1, n2=NULL, alpha, alternative="two.sided", tDensity=FALSE) {
  safeTTestStat("t"=t, "deltaS"=deltaS, "n1"=n1, "n2"=n2, "alternative"=alternative, "tDensity"=tDensity) - 1/alpha
}

#' Designs a Safe Experiment to Test Means
#'
#' Designs a safe experiment for a prespecified minimum clinical relevant effect size, tolerable type I and
#' type II error. Outputs a list that includes (1) the deltaS that defines the safe test, and (2) nPlan, the sample
#' size to  plan for.
#'
#' @param deltaMin numeric that defines the minimal relevant effect size, the smallests effect size that we want to
#' detect.
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control --independent on n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule S10 > 1/alpha
#' @param beta numeric in (0, 1) that speficies the tolerable type II error control necessary to calculate both "n"
#' and "deltaS". Note that 1-beta defines the power.
#' @param lowDelta numeric that defines the smallest delta of our search space for the test-defining deltaS
#' @param highDelta numeric that defines the largest delta of our search space for the test-defining deltaS
#' @param tol a number that defines the stepsizes between the lowDelta and highDelta
#' @param lowN integer that defines the smallest n of our search space for n
#' @param highN integer that defines the largest n of our search space for n. This might be the largest n that we
#' are able to fund.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less"
#' @param mu0 a number indicating the hypothesised true value of the mean under the null. For the moment mu0=0
#' @param testType either one of "oneSampleT", "pairedSampleT", "twoSampleT"
#' @param sampleSizeRatio numeric representing n2/n1. If is.null(n2) then sampleSizeRatio=1
#' @param logging logical, if TRUE return altSThreshes
#' @param ... further arguments to be passed to or from methods, but mainly to perform do.calls
#'
#' @return Returns a safeDesign object that includes:
#'
#' \describe{
#'   \item{n2Plan}{the sample size of the second group when testType=="twoSampleT" or "pairedSampleT", otherwise NULL}
#'   \item{nEffPlan}{the resulting effective sample size when testType=="twoSampleT", otherwise non-existing}
#'   \item{deltaS}{the deltaS that defines the safe test}
#'   \item{deltaMin}{the minimal clinical effect size provided by the user}
#'   \item{alpha}{the tolerable type I error provided by the user}
#'   \item{beta}{the tolerable type II error provided by the user}
#'   \item{lowDelta}{the smallest delta of the search space for delta provided by the user}
#'   \item{highDelta}{the largest delta of the search space for delta provided by the user}
#'   \item{tol}{the step size between lowDelta and highDelta provided by the user}
#'   \item{lowN}{the smallest n of the search space for n provided by the user}
#'   \item{highN}{the largest n of the search space for n provided by the user}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user}
#'   \item{testType}{any of "oneSampleT", "pairedSampleT", "twoSampleT" provided by the user}
#'   \item{sampleSizeRatio}{default is 1, only used when testType=="twoSampleT" and defines n2=sampleSizeRatio*n1}
#'   \item{pilot}{FALSE to indicate that the design is not a pilot study}
#'   \item{call}{the expression with which this function is called}
#'   \item{altSThreshes}{if logging=TRUE then shows the s-values at the t-value corresponding to the type II error
#'   under the alternative at deltaMin}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.8, alpha=0.08, beta=0.01, alternative="greater")
#' designObj
designSafeT <- function(deltaMin, alpha=0.05, beta=0.2, alternative=c("two.sided", "greater", "less"),
                        mu0=0, lowDelta=0.01, highDelta=1.5*abs(deltaMin), tol=0.01, lowN=3L, highN=100L,
                        testType=c("oneSampleT", "pairedSampleT", "twoSampleT"), sampleSizeRatio=1,
                        logging=FALSE, ...) {

  stopifnot(alpha > 0, alpha < 1, beta > 0, beta < 1)

  result <- list("n1Plan"=NULL, "n2Plan"=NULL, "mu0"=mu0, "deltaS"=NULL,
                 "deltaMin"=deltaMin, "alpha"=alpha, "beta"=beta,
                 "lowDelta"=lowDelta, "highDelta"=highDelta, "tol"=tol,
                 "lowN"=lowN, "highN"=highN, "alternative"=alternative, "testType"=testType,
                 "sampleSizeRatio"=sampleSizeRatio, "pilot"=FALSE, "call"=sys.call())
  class(result) <- "safeTDesign"

  deltaMin <- abs(deltaMin)

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  sCutOff <- 1/alpha

  nDefinitions <- defineTTestN("lowN"=lowN, "highN"=highN, "sampleSizeRatio"=sampleSizeRatio, "testType"=testType)

  if (testType=="pairedSampleT") {
    paired <- TRUE
  } else {
    paired <- FALSE
  }

  n1 <- nDefinitions[["n1"]]
  n2 <- nDefinitions[["n2"]]
  candidateNEff <- nDefinitions[["candidateNEff"]]
  candidateNu <- nDefinitions[["candidateNu"]]

  if (alternative=="two.sided") {
    candidateFNcp <- candidateNEff*deltaMin^2
  } else {
    candidateTNcp <- sqrt(candidateNEff)*deltaMin
  }

  # TODO(Alexander): Should highDelta just be deltaMin.
  #   Perhaps show that deltaS < deltaMin for alpha, beta. Use monotonicity
  #
  candidateDeltas <- seq(from=lowDelta, to=highDelta, by=tol)

  for (i in seq_along(candidateNEff)) {
    if (alternative=="two.sided") {
      deltaMinThresh <- sqrt(stats::qf("p"=beta, "df1"=1, "df2"=candidateNu[i], "ncp"=candidateFNcp[i])) #*deltaMin^2)
    } else {
      deltaMinThresh <- stats::qt("p"=beta, "df"=candidateNu[i], "ncp"=candidateTNcp[i])
    }
    # TODO(Alexander): Under the assumption that this is unimodal, then stop once the value goes down
    if (testType=="twoSampleT") {
      altSThreshes <- purrr::map_dbl(".x"=candidateDeltas, ".f"=safeTTestStat, "t"=deltaMinThresh,
                                     "n1"=n1[i], "n2"=n2[i], "alternative"=alternative, "paired"=paired)
    } else if (testType %in% c("oneSampleT", "pairedSampleT")) {
      altSThreshes <- purrr::map_dbl(".x"=candidateDeltas, ".f"=safeTTestStat, "t"=deltaMinThresh,
                                     "n1"=candidateNEff[i], "n2"=n2[i], "alternative"=alternative, "paired"=paired)
    }

    if (max(altSThreshes) >= sCutOff) {
      nEff <- candidateNEff[i]
      if (testType=="twoSampleT") {
        result[["n1Plan"]] <- ceiling(n1[i])
        result[["n2Plan"]] <- ceiling(n2[i])
        result[["nEffPlan"]] <- nEff
      } else if (testType %in% c("oneSampleT", "pairedSampleT")) {
        result[["n1Plan"]] <- nEff

        if (testType=="pairedSampleT") {
          result[["n2Plan"]] <- nEff
        }
      }

      deltaIndex <- which(altSThreshes >= sCutOff)[1]
      deltaS <- candidateDeltas[deltaIndex]

      if (alternative=="less")
        deltaS <- -deltaS

      result[["deltaS"]] <- deltaS
      result[["testType"]] <- testType

      if (isTRUE(logging))
        result[["altSThreshes"]] <- altSThreshes

      break()
    }
  }
  return(result)
}


#' Prints a safeTDesign object
#'
#' @param x a safeTDesign object
#' @param ... further arguments to be passed to or from methods.
#'
#' @return prints a safeTDesign object
#' @export
#'
#' @examples
#' safeDesignObj <- designSafeT(0.8)
#' print(safeDesignObj)
print.safeTDesign <- function(x, ...) {
  analysisName <- getNameTestType(testType = x[["testType"]])

  cat("\n")
  cat(paste("       ", analysisName, "\n"))
  cat("\n")

  if (isFALSE(x[["pilot"]])) {
    if (is.null(x[["n2Plan"]])) {
      cat("requires an experiment with a sample size of: ")
      cat("\n")
      cat(paste("    n1Plan =", x[["n1Plan"]]))
      cat("\n")
    } else {
      cat("Requires an experiment with sample sizes: ")
      cat("\n")
      cat(paste("    n1Plan =", x[["n1Plan"]], "and n2Plan =", x[["n2Plan"]]))
      cat("\n")
    }
    cat("to find an effect size of at least: ")
    cat("\n")
    cat("    deltaMin =", round5(x[["deltaMin"]]))
    cat("\n")
    cat("\n")

    cat("with:")
    cat("\n")
    cat("    power = ", 1 - x[["beta"]], " (thus, beta = ", x[["beta"]], ")", sep="")
    cat("\n")

    cat("under the alternative:")
    cat("\n")
    cat("   ", getNameAlternative(x[["alternative"]], x[["testType"]]))
    cat("\n")
    cat("\n")

    cat("Based on the decision rule S > 1/alpha:")
    cat("\n")
    cat("    S > ", round5(1/x[["alpha"]]), sep="")
    cat("\n")

    cat("which occurs with chance less than:")
    cat("\n")
    cat("    alpha =", x[["alpha"]])
    cat("\n")

    cat("under iid normally distributed data and the null hypothesis:")
    cat("\n")
    cat("    mu =", x[["mu0"]])
  } else {
    cat("The experiment is not planned.")
    cat("\n")
    cat("This design object only valid for experiments with:")
    cat("\n")

    if (is.null(x[["n2Plan"]])) {
      cat("    n1 =", x[["n1Plan"]])
      cat("\n")
    } else {
      cat("    n1 =", x[["n1Plan"]], "and n2 =", x[["n2Plan"]])
      cat("\n")
    }
  }
}


#' Prints a safeTDesign object
#'
#' @param x a safeTResult object
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Prints a safeTDesign object
#' @export
#'
#' @examples
#' safeDesignObj <- designSafeT(0.7)
#' safeTTest(rnorm(10), designObj=safeDesignObj)
print.safeTResult <- function(x, ...) {
  designObj <- x[["designObj"]]
  testType <- designObj[["testType"]]

  analysisName <- getNameTestType("testType"=testType)
  alternativeName <- getNameAlternative("alternative"=x[["alternative"]], "testType"=testType)

  cat("\n")
  cat(paste("       ", analysisName, "\n"))
  cat("\n")

  cat("Data:", x[["dataName"]])
  cat("\n")
  cat("sample estimates:")
  cat("\n")
  print(round5(x[["estimate"]]))

  cat("\n")
  cat("Test summary: ")
  cat("t = ", round5(x[["statistic"]]), ", df = ", round5(x[["parameter"]]), ".", sep="")
  cat("\n")

  if (designObj[["pilot"]]) {
    cat("The pilot test is based on an exploratory alpha =", designObj[["alpha"]])
    cat("\n")
    cat("and resulted in:  s-value =", round5(x[["sValue"]]))
    cat("\n")
    cat("Alternative hypothesis:")
  } else {
    cat("The test designed with alpha =", designObj[["alpha"]])
    cat("\n")
    cat("s-value =", round5(x[["sValue"]]), "> 1/alpha =", round5(1/designObj[["alpha"]]), ":",
        x[["sValue"]] > 1/designObj[["alpha"]])
    cat("\n")
    # Iets over n1Plan, n2Plan, etc

    cat("\n")
    if (is.null(designObj[["n2Plan"]])) {
      cat(paste("Experiments required n1Plan =", designObj[["n1Plan"]], "samples."))
    } else {
      cat(paste("Experiments required n1Plan =", designObj[["n1Plan"]], "and n2Plan =",
                designObj[["n2Plan"]], "samples."))
    }
    cat("\n")

    n1Diff <- designObj[["n1Plan"]] - x[["n1"]]

    if (!is.null(designObj[["n2Plan"]])) {
      n2Diff <- designObj[["n2Plan"]] - x[["n2"]]
    } else {
      # Note Dummy
      n2Diff <- 0
    }

    if (n1Diff > 0 || n2Diff > 0) {
      cat("    Note: ")
      if (n1Diff > 0) {
        cat("n1Plan - n1 = ", n1Diff, ", ", sep="")
      }
      if (n2Diff > 0) {
        cat("n2Plan - n2 =", n2Diff)
      }
      cat("\n")
    }

    cat("to guarantee a power = ", round5(1 - designObj[["beta"]]),
        " (beta =", round5(designObj[["beta"]]), ").", sep="")
    cat("\n")
    cat("under the alternative hypothesis:")
  }
  cat("\n")
  cat(alternativeName)
  cat("\n")

  if (isFALSE(designObj[["pilot"]])) {
    cat("and deltaMin =", designObj[["deltaMin"]])
  }

}

#' Helper function that outputs the sample sizes, effect sample sizes and the degrees of freedom depending on
#' the type of t-test.
#'
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#'
#' @return Returns the sample sizes and degrees of freedom
#'
#' @examples
#' \dontrun{
#' defineTTestN()
#' }
defineTTestN <- function(lowN=3, highN=100, sampleSizeRatio=1, testType=c("oneSampleT", "pairedSampleT", "twoSampleT")) {
  testType <- match.arg(testType)

  if (testType=="twoSampleT") {
    n1 <- lowN:highN
    n2 <- sampleSizeRatio*n1
    candidateNEff <- sampleSizeRatio/(1+sampleSizeRatio)*n1
    candidateNu <- (1+sampleSizeRatio)*n1-2
  } else if (testType %in% c("oneSampleT", "pairedSampleT")) {
    n1 <- lowN:highN
    n2 <- NULL
    candidateNEff <- n1
    candidateNu <- candidateNEff-1
  }
  result <- list("n1"=n1, "n2"=n2, "candidateNEff"=candidateNEff, "candidateNu"=candidateNu)
  return(result)
}


#' Pretends that the observed sample sizes are exactly as they were planned for.
#'
#' "Designs" a safe experiment for a prespecified tolerable type I error by pretending that the sample size is "known"
#' and fixed ahead of time. Outputs a list that includes the deltaS that defines the safe test.
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#' @param n1,n2 observed sample sizes
#' @param inverseMethod logical, always TRUE for the moment
#' @param paired logical, if TRUE then paired t-test
#' @param logging ‘logical, if TRUE, then add invSToTThresh to output
#'
#' @return Returns a safeDesign object
#' \describe{
#'   \item{n1Plan}{the sample size to plan for}
#'   \item{n2Plan}{the sample size of the second group when testType=="twoSampleT" or "pairedSampleT", otherwise NULL}
#'   \item{nEffPlan}{the resulting effective sample size when testType=="twoSampleT", otherwise non-existing}
#'   \item{deltaS}{the deltaS that defines the safe test}
#'   \item{deltaMin}{NULL, no deltaMin specified because it's a pilot}
#'   \item{alpha}{the tolerable type I error provided by the user}
#'   \item{beta}{NULL, no tolerable type II error specified}
#'   \item{lowDelta}{the smallest delta of the search space for delta provided by the user}
#'   \item{highDelta}{the largest delta of the search space for delta provided by the user}
#'   \item{tol}{the step size between lowDelta and highDelta provided by the user}
#'   \item{lowN}{NULL}
#'   \item{highN}{NULL}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user}
#'   \item{testType}{any of "oneSampleT", "pairedSampleT", "twoSampleT" provided by the user}
#'   \item{sampleSizeRatio}{default is 1, only used when testType=="twoSampleT" and defines n2=sampleSizeRatio*n1}
#'   \item{pilot}{TRUE to indicate that the design is a pilot study. The assumption is that the sample sizes are as if they were planned for, thus, known in advance.}
#'   \item{call}{the expression with which this function is called}
#'   \item{error}{the error estimated from the inverse function}
#'   \item{invSToTThresh}{if logging=TRUE then shows the inverse of the t-threshold at various test defining deltaS.}
#' }
#' @export
#'
#' @examples
#' designPilotSafeT(n1=30)
designPilotSafeT <- function(n1=50, n2=NULL, alpha=0.05, mu0=0, alternative=c("two.sided", "greater", "less"),
                             lowDelta=0.01, highDelta=1.2, tol=0.01, inverseMethod=TRUE,
                             logging=FALSE, paired=FALSE) {
  # TODO(Alexander): Check relation with forward method, that is, the least conservative test and maximally powered
  # Perhaps trade-off? "inverseMethod" refers to solving minimum of deltaS \mapsto S_{deltaS}^{-1}(1/alpha)
  #
  #
  # TODO(Alexander):
  #       - Add warning when deltaS == lowDelta
  #       - Better: Characterise deltaS as a funciton of alpha and n using asymptotic approximation of 1F1
  #
  # TODO(Alexander): Add some bounds on L, or do a presearch
  #
  # stopifnot(n >= 3)
  #
  #     Trick bound the estimation error using Chebyshev: Note do this on log (safeTTestStat)
  #
  alternative <- match.arg(alternative)
  stopifnot(n1 > 2, n2 > 2)


  if (!is.null(n2)) {
    sampleSizeRatio <- n2/n1

    if (paired) {
      if (n1!=n2) {
        stop("Paired design specified, but n1 not equal n2")
      }
      testType <- "pairedSampleT"
    } else {
      testType <- "twoSampleT"
    }
  } else {
    sampleSizeRatio <- 1
    testType <- "oneSampleT"

    if (isTRUE(paired)) {
      stop("Paired designed specified, but n2 not provided")
    }
  }

  result <- list("n1Plan"=n1, "n2Plan"=n2, "deltaS"=NA,
                 "deltaMin"=NULL, "alpha"=alpha, "beta"=NULL,
                 "lowDelta"=lowDelta, "highDelta"=highDelta, "tol"=tol,
                 "lowN"=NULL, "highN"=NULL, "alternative"=alternative, "testType"=testType,
                 "sampleSizeRatio"=sampleSizeRatio, "pilot"=TRUE, "call"=sys.call())
  class(result) <- "safeTDesign"

  candidateDeltas <- seq(lowDelta, highDelta, by=tol)

  if (inverseMethod) {
    invSafeTTestStatAlpha <- function(x) {
      stats::uniroot("f"=safeTTestStatAlpha, "lower"=0, "upper"=3e10, "tol"=1e-9, "deltaS"=x,
              "n1"=n1, "n2"=n2, "alpha"=alpha, "alternative"=alternative)
    }

    # Note(Alexander): tryOrFAilWithNA fixes the problem that there's a possibility
    #   that deltaS too small, the function values are then both negative or both positive.
    #
    invSafeTTestStatAlphaRoot <- function(x){tryOrFailWithNA(invSafeTTestStatAlpha(x)$root)}
    invSToTThresh <- purrr::map_dbl(candidateDeltas, invSafeTTestStatAlphaRoot)
    mPIndex <- which(invSToTThresh==min(invSToTThresh, na.rm=TRUE))

    if (mPIndex==length(candidateDeltas)) {
      # Note(Alexander): Check that mPIndex is not the last one.
      errorMsg <- "The test defining deltaS is equal to highDelta. Rerun with do.call on the output object"
      lowDelta <- highDelta
      highDelta <- (length(candidateDeltas)-1)*tol+lowDelta
      result[["lowDelta"]] <- lowDelta
      result[["highDelta"]] <- highDelta
    } else if (mPIndex==1) {
      errorMsg <- "The test defining deltaS is equal to lowDelta. Rerun with do.call on the output object"
      highDelta <- lowDelta
      lowDelta <- highDelta-(length(candidateDeltas)-1)*tol
      result[["lowDelta"]] <- lowDelta
      result[["highDelta"]] <- highDelta
    }

    if (isTRUE(logging))
      result[["invSToTThresh"]] <- invSToTThresh

    deltaS <- candidateDeltas[mPIndex]

    if (alternative=="less")
      deltaS <- -deltaS

    result[["deltaS"]] <- deltaS
    result[["error"]] <- invSafeTTestStatAlpha(candidateDeltas[mPIndex])$estim.prec
  } else {
    # TODO(Alexander): Check relation with forward method, that is, the least conservative test and maximally powered
    # Perhaps trade-off? "inverseMethod" refers to solving minimum of deltaS \mapsto S_{deltaS}^{-1}(1/alpha)
  }
  # TODO(Alexander): By some monotonicity can we only look at the largest or the smallest?
  #
  # designFreqT(deltaMin=deltaMin, alpha=alpha, beta=beta, lowN=lowN, highN=highN)
  return(result)
}



#' Plots the sample sizes necessary for a tolerable alpha and beta as a function of deltaMin
#'
#' For given tolerable alpha and beta, as a function of the minimal clinical relevant effect size deltaMin, plots (1) the
#' sample sizes to plan for for a safe test (2) the frequentist test, (3) the average sample size necessary due to
#' optional stopping.
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#' @param maxN numeric, the maximum number of samples one has budget for to collect data
#' @param deltaFactor numeric, a factor to robustify the sequential determination (e.g., from deltaTrue = 0.9, to
#' deltaTrue = 0.8) of lowDelta
#' @param nFactor numeric, a factor to robustify the sequential determination (e.g., from deltaTrue = 0.9, to
#' deltaTrue = 0.8) of highN
#' @param simulateSafeOptioStop logical, if TRUE then provides
#' @param logging logical, if TRUE then output all the safe designs objects including mean n stop if
#' simulateSafeOptioStop==TRUE
#' @param backTest logical, if TRUE it provides the frequentist sample size necessary to attain the power that the
#' safe test attains due to optional stopping
#' @param freqPlot logical, if TRUE plot frequentist sample size profiles

#'
#' @return Plot of the sample size profiles for tolerable type I and type II error, also outputs results object
#' @export
#'
#' @examples
#' plotSafeTDesignSampleSizeProfile()
plotSafeTDesignSampleSizeProfile <- function(alpha=0.05, beta=0.2, maxN=200, lowDelta=0.01, highDelta=1, tol=0.1,
                                             testType=c("oneSampleT", "pairedSampleT", "twoSampleT"),
                                             alternative=c("two.sided", "greater", "less"), sampleSizeRatio=1,
                                             mIter=1000L, deltaFactor=0.5, nFactor=2, simulateSafeOptioStop=FALSE,
                                             logging=FALSE, backTest=FALSE, seedNumber=NULL, freqPlot=FALSE, pb=TRUE,
                                             ...) {

  stopifnot(lowDelta < highDelta, alpha > 0, beta > 0, alpha < 1, beta < 1)

  # Order from high to low
  deltaDomain <- -seq(-highDelta, -lowDelta, by=tol)
  testType <- match.arg(testType)

  result <- list("alpha"=alpha, "beta"=beta, "maxN"=maxN, "deltaDomain"=deltaDomain)

  lastDeltaIndex <- length(deltaDomain)

  if (lastDeltaIndex < 1)
    stop("Either maxN or deltaDomain is too small. Please lower lowDelta or make highDelta larger")


  if (testType=="pairedSampleT") {
    paired <- TRUE
  } else {
    paired <- FALSE
  }

  allN1PlanFreq <- vector("integer", lastDeltaIndex)


  # 1. Freq design  ---------------------------------------------------------------------
  freqDesign <- list("n1PlanFreq"=3)

  for (i in seq.int(lastDeltaIndex)) {
    if (i==1) {
      tempLowN <- 3
    } else {
      # Note(Alexander): Use previous found n1
      # TODO(Alexander): Show that as deltaMin decreases that n1PlanFreq increases
      tempLowN <- freqDesign[["n1PlanFreq"]]
    }

    freqDesign <- designFreqT("deltaMin"=deltaDomain[i], "alpha"=alpha, "beta"=beta, "lowN"=tempLowN,
                              "highN"=maxN, "sampleSizeRatio"=sampleSizeRatio)

    if (is.null(freqDesign[["n1PlanFreq"]]) || is.na(freqDesign[["n1PlanFreq"]])) {
      lastDeltaIndex <- i-1

      # Note(Alexander): Prune
      allN1PlanFreq <- allN1PlanFreq[1:lastDeltaIndex]
      break()
    }

    allN1PlanFreq[i] <- freqDesign[["n1PlanFreq"]]
  }


  # #### 1.a. Plots freq
  # graphics::plot(deltaDomain[seq.int(lastDeltaIndex)], allN1PlanFreq[seq.int(lastDeltaIndex)],
  #                lty=3, lwd=2, type="l", col="darkgrey", ylab="n1", xlab=expression(delta["min"]))
  #
  # abline(h=maxN, col="red", lty=2)
  #
  # legend("topright", legend = c("Freq design", "max n"),
  #        col = c("darkgrey", "red"),
  #        lty = c(3, 2), bty="n")

  # 2. Safe design  ---------------------------------------------------------------------
  #
  allDeltaS <- allN1PlanSafe <- vector("integer", lastDeltaIndex)
  allSafeDesignObj <- vector("list", lastDeltaIndex)

  # TODO(Alexander): Show that for fixed theta that freqN < safeN,
  # TODO(Alexander): Instead, of doing this based on the previous one (deltaMin), try to be faster
  # by going around a guessed quantity based on a factor of safeDesignObj$n/freqDesign$n
  #
  # TODO(Alexander): Show that due to monotonicity that we can take "highDelta"=deltaDomain[i-1]
  for (i in seq.int(lastDeltaIndex)) {
    if (i==1) {
      tempLowDelta <- lowDelta
      tempHighDelta <- highDelta
      tempLowN <- 3
      tempHighN <- maxN
    } else {
      # Note(Alexander): Use previous found deltaS times a correction factor as a lowerbound for the search space of
      # deltaS
      tempLowDelta <- deltaFactor*safeDesignObj[["deltaS"]]/deltaDomain[i-1]*deltaDomain[i]

      # Note(Alexander): Use previous true deltaMin > previous deltaS as an upper bound
      # TODO(Alexander): Show that as deltaMin decreases that deltaS dereases
      tempHighDelta <- deltaDomain[i-1]

      # Note(Alexander): Use previously found design
      # TODO(Alexander): Show that as deltaMin decreases that n1Plan increases
      tempLowN <- safeDesignObj[["n1Plan"]]

      # Note(Alexander): Use previously found n1Plan times a factor as an upper bound
      # TODO(Alexander): Show that as deltaMin decreases that n1Plan increases
      tempHighN <- ceiling(nFactor * safeDesignObj[["n1Plan"]]/allN1PlanFreq[i-1]*allN1PlanFreq[i])
    }

    safeDesignObj <- designSafeT("deltaMin"=deltaDomain[i], "alpha"=alpha, "beta"=beta, "alternative"=alternative,
                                 "lowDelta"=tempLowDelta, "highDelta"=tempHighDelta, "lowN"=tempLowN,
                                 "highN"=tempHighN, "testType"=testType, "sampleSizeRatio"=sampleSizeRatio)

    if (is.null(safeDesignObj[["n1Plan"]]) || is.na(safeDesignObj[["n1Plan"]])) {
      lastDeltaIndex <- i-1
      break()
    }

    # TODO(Alexander): Not necessary anymore, with normal function call
    safeDesignObj[["n1PlanFreq"]] <- allN1PlanFreq[i]
    allSafeDesignObj[[i]] <- safeDesignObj
    allN1PlanSafe[i] <- safeDesignObj[["n1Plan"]]
    allDeltaS[i] <- safeDesignObj[["deltaS"]]
  }

  deltaDomain <- deltaDomain[1:lastDeltaIndex]
  allDeltaS <- allDeltaS[1:lastDeltaIndex]

  # TODO(Alexander) Optional?
  # graphics::plot(deltaDomain, allDeltaS)


  maxDeltaDomain <- max(deltaDomain)
  minDeltaDomain <- min(deltaDomain)


  # Store in output
  result[["deltaDomain"]] <- deltaDomain
  result[["allN1PlanFreq"]] <- allN1PlanFreq
  result[["allN1PlanSafe"]] <- allN1PlanSafe
  result[["allDeltaS"]] <- allDeltaS

  # 2.a. Plot Safe -----
  graphics::par(cex.main=1.5, mar=c(5, 6, 4, 7)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
                font.lab=2, cex.axis=1.3, bty="n", las=1)


  graphics::plot(deltaDomain, allN1PlanSafe, type="l", col="blue", lty=1, lwd=2, xlim=c(minDeltaDomain, maxDeltaDomain),
                 ylab="n1", xlab=expression(delta["min"]),
                 main=bquote(~alpha == ~.(alpha) ~ "and" ~beta== ~.(beta)))

  if (freqPlot) {
    graphics::lines(deltaDomain, allN1PlanFreq, col="darkgrey", lwd=2, lty=3)
    legendName <- c("Safe design", "Freq design", "max n")
    legendCol <- c("blue", "darkgrey", "red")
    legendLty <- c(1, 3, 2)
  } else {
    legendName <- c("Safe design", "max n")
    legendCol <- c("blue", "red")
    legendLty <- c(1, 2)
  }

  graphics::abline(h=maxN, col="red", lty=2)

  graphics::legend("topright", legend = legendName,
                   col = legendCol, lty = legendLty, bty="n")

  # 3. Run simulations  ---------------------------------------------------------------------
  #
  if (simulateSafeOptioStop) {
    allNMean <- allProbLeqNFreq <- vector("integer", lastDeltaIndex)

    if (backTest)
      allNBack <- vector("integer", lastDeltaIndex)

    if (pb)
      pbOptioStop <- utils::txtProgressBar("style"=1)

    for (i in seq.int(lastDeltaIndex)) {
      safeDesignObj <- allSafeDesignObj[[i]]

      if (i==1) {
        tempLowN <- 3
      } else {
        # Note(Alexander): Use for lowN the smallest N found in the previous simulations
        tempLowN <- simObj[["safeSim"]][["lowN"]]
      }

      simObj <- replicateTTests("n1Plan"=safeDesignObj[["n1Plan"]], "n2Plan"=safeDesignObj[["n2Plan"]],
                                "deltaTrue"=deltaDomain[i], "paired"=paired, "alternative"=alternative,
                                "lowN"=tempLowN, "alpha"=alpha, "deltaS"=safeDesignObj[["deltaS"]],
                                "n1PlanFreq"=allN1PlanFreq[i], "pb"=FALSE)
      allNMean[i] <- simObj[["safeSim"]][["nMean"]]
      allProbLeqNFreq[i] <- simObj[["safeSim"]][["probLeqN1PlanFreq"]]

      if (backTest) {
        backFreqDesign <- designFreqT("deltaMin"=deltaDomain[i], "alpha"=alpha,
                                      "beta"=1-safeDesign[["safeSim"]][["powerOptioStop"]],
                                      "alternative"=alternative, "testType"=testType,
                                      "sampleSizeRatio"=sampleSizeRatio)

        safeDesign[["safeSim"]][["nBack"]] <- backFreqDesign[["n1PlanFreq"]]
        allNBack[i] <- backFreqDesign[["n1PlanFreq"]]
      } # End back test

      safeDesignObj <- utils::modifyList(safeDesignObj, simObj)
      allSafeDesignObj[[i]] <- safeDesignObj

      if (pb)
        utils::setTxtProgressBar("pb"=pbOptioStop, "value"=i/lastDeltaIndex)

    } # End looping over deltaDomain as deltaTrue

    if (pb)
      close(pbOptioStop)

    result[["allNMean"]] <- allNMean
    result[["allProbLeqNFreq"]] <- allProbLeqNFreq

    if (backTest)
      result[["allNBack"]] <- allNBack

    # 3.a. Plot Sim  -----
    graphics::par(cex.main=1.5, mar=c(5, 6, 4, 7)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
                  font.lab=2, cex.axis=1.3, bty="n", las=1)

    graphics::plot(deltaDomain, allN1PlanSafe, type="l", col="blue", lty=2, lwd=2, xlim=c(minDeltaDomain, maxDeltaDomain),
                   ylab="n1", xlab=expression(delta["min"]),
                   main=bquote(~alpha == ~.(alpha) ~ "and" ~beta== ~.(beta)))
    graphics::abline(h=maxN, col="red", lty=2)
    graphics::lines(deltaDomain, allNMean, col="black", lwd=2, lty=1)

    if (freqPlot) {
      graphics::lines(deltaDomain, allN1PlanFreq, col="darkgrey", lwd=2, lty=3)
      legendName <- c("Average n", "Safe design", "Freq design", "max n")
      legendCol <- c("black", "blue", "darkgrey", "red")
      legendLty <- c(1, 2, 3, 2)
    } else {
      legendName <- c("Average n", "Safe design", "max n")
      legendCol <- c("black", "blue", "red")
      legendLty <- c(1, 2, 2)
    }

    graphics::legend("topright", legend = legendName, col = legendCol, lty=legendLty, bty="n")
  }

  if (logging)
    result[["allSafeDesignObj"]] <- allSafeDesignObj

  return(result)
}



#' Simulate multiple data sets to show the effects of optional testing for safe (and frequentist) tests.
#'
#' @param n1Plan integer, that defines the maximum number of samples to plan for (according to the safe test,
#' use designSafeT to find this)
#' @param n2Plan optional integer, that defines the maximum number of samples of the second group to plan for
#' @param deltaTrue numeric, the value of the true effect size (test-relevant parameter)
#' @param muGlobal numeric, the true global mean of a paired or two-sample t-test. Its value shouldn't matter for the
#' test. This parameter treated is treated as a nuisance.
#' @param sigmaTrue numeric > 0,the true standard deviation of the data. Its value shouldn't matter for the test.
#' This parameter treated is treated as a nuisance.
#' @param paired logical, true if the simulated data are paired.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less"
#' @param lowN the smallest number of samples (first group) at which monitoring of the tests begins
#' @param mIter the number of replications, that is, experiments with max samples n1Plan and n2Plan
#' @param alpha the tolerable type I error to be conserved. Also defines the decision rule s > 1/alpha, and for
#' frequentist tests the decision rule is p < alpha.
#' @param safeOptioStop logical, TRUE implies that optional stopping simulation is performed for the safe test
#' @param deltaS numeric, the safe test defining deltaS (use designSafeT to find this)
#' @param freqOptioStop logical, TRUE implies that optional stopping simulation is performed for the frequentist test
#' @param n1PlanFreq integer, that defines the maximum number of samples to plan for (according to the frequentist
#' test,use designFreqT to find this)
#' @param n2PlanFreq optional integer, that defines the maximum number of samples of the second group to plan for
#' @param seedNumber To set the seed for the simulated data
#' @param logging logical, if TRUE, then return the sampled sample sizes
#' @param pb logical, if TRUE, then show progress bar
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a safeSim object.
#' @export
#'
#' @examples
#'
#'# Design safe test
#' alpha <- 0.05
#' beta <- 0.20
#' designObj <- designSafeT(1, alpha=alpha, beta=beta)
#'
#' # Design frequentist test
#' freqObj <- designFreqT(1, alpha=alpha, beta=beta)
#'
#' # Simulate under the alternative with deltaTrue=deltaMin
#' simResults <- replicateTTests(n1Plan=designObj$n1Plan, deltaTrue=1, deltaS=designObj$deltaS,
#' n1PlanFreq=freqObj$n1PlanFreq)
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
#' breaks=seq.int(designObj$n1Plan))
#'
#' # Simulate under the alternative with deltaTrue > deltaMin
#' simResults <- replicateTTests(n1Plan=designObj$n1Plan, deltaTrue=1.5, deltaS=designObj$deltaS,
#' n1PlanFreq=freqObj$n1PlanFreq)
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
#' breaks=seq.int(designObj$n1Plan))
#'
#' # Under the null deltaTrue=0
#' simResults <- replicateTTests(n1Plan=designObj$n1Plan, deltaTrue=0, deltaS=designObj$deltaS,
#' n1PlanFreq=freqObj$n1PlanFreq, freqOptioStop=TRUE)
#'
#'# Should be lower than alpha, because if the null is true, P(S > 1/alpha) < alpha for all n
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
#'  simResults$freqSim$powerOptioStop > alpha
replicateTTests <- function(n1Plan, n2Plan=NULL, deltaTrue, muGlobal=0, sigmaTrue=1, paired=FALSE,
                            alternative=c("two.sided", "greater", "less"), lowN=3,
                            mIter=1000L, alpha=0.05,
                            safeOptioStop=TRUE, deltaS=NULL,
                            freqOptioStop=FALSE, n1PlanFreq=NULL, n2PlanFreq=NULL,
                            logging=TRUE, seedNumber=NULL, pb=TRUE, ...) {

  stopifnot(n1Plan > 0, n1Plan > lowN, mIter > 0, alpha > 0, alpha < 1,
            any(safeOptioStop, freqOptioStop))

  alternative <- match.arg(alternative)

  result <- list(n1Plan=n1Plan, n2Plan=n2Plan, deltaTrue=deltaTrue, muGlobal=muGlobal, paired=paired,
                 alternative=alternative, lowN=lowN, mIter=mIter, alpha=alpha,
                 deltaS=deltaS, n1PlanFreq=n1PlanFreq, n2PlanFreq=n2PlanFreq, safeSim=list(), freqSim=list())

  if (safeOptioStop) {
    if (is.null(deltaS)) {
      stop(paste("To simulate safe t-tests results under optional stopping, this function 'replicateTTests' requires",
                 "the specification of the safe test with a deltaS. This deltaS can be found by running",
                 "the 'designSafeT' function")
      )
    }

    if (paired && n1Plan != n2Plan)
      stop("For a paired t-test n2Plan needs to equal n1Plan")

    safeSim <- list(powerOptioStop=NA, powerAtN1Plan=NA, nMean=NA, probLeqN1PlanFreq=NA, probLessNDesign=NA, lowN=NA)

    allSafeN <- rep(n1Plan, times=mIter)
    safeDecisionAtN <- allSafeDecisions <- vector("integer", mIter)
  }

  if (freqOptioStop) {
    if (!safeOptioStop) {
      if (is.null(n1Plan)) {
        warning("No n1PlanFreq specified, use n1Plan instead.")
        n1PlanFreq <- n1Plan
        n2PlanFreq <- n2Plan
      }
    }

    # Note(Alexander): This means that n1Plan and n2Plan refer to the planned samples of the safe tests

    if (is.null(n1PlanFreq)) {
      stop(paste("To simulate frequentist t-tests results under optional stopping, this",
                 "function 'replicateTTests' requires the specification of n1PlanFreq. To figure out how many",
                 "samples one requires in a frequentist test, please run the 'designFreqT' function.")
      )
    }


    if (!is.null(n2Plan) && is.null(n2PlanFreq)) {
      stop(paste("To simulate a two-sample frequentist t-tests results under optional stopping, this",
                 "function 'replicateTTests' requires the specification of n1PlanFreq. To figure out how many ",
                 "samples one requires in a frequentist test, please run the 'designFreqT' function.")
      )
    }

    if (paired && n1PlanFreq != n2PlanFreq)
      stop("For a paired t-test n2PlanFreq needs to equal n1PlanFreq")

    freqSim <- list(powerOptioStop=NA, powerAtN1Plan=NA, nMean=NA, probLessNDesign=NA, lowN=NA)

    allFreqN <- rep(n1PlanFreq, times=mIter)
    freqDecisionAtN <- allFreqDecisions <- vector("integer", mIter)
  }

  set.seed(seedNumber)

  if (is.null(n2Plan)) {
    sampleSizeRatio <- 1

    repData1 <- stats::rnorm("n"=n1Plan*mIter, "mean"=deltaTrue*sigmaTrue, "sd"=sigmaTrue)
    repData1 <- matrix(repData1, "ncol"=n1Plan, "nrow"=mIter)
    repData2 <- NULL
  } else {
    if (paired) {
      sampleSizeRatio <- 1

      repData1 <- stats::rnorm("n"=n1Plan*mIter, "mean"=muGlobal + deltaTrue*sigmaTrue/sqrt(2), "sd"=sigmaTrue)
      repData1 <- matrix(repData1, "ncol"=n1Plan, "nrow"=mIter)
      repData2 <- stats::rnorm("n"=n2Plan*mIter, "mean"=muGlobal - deltaTrue*sigmaTrue/sqrt(2), "sd"=sigmaTrue)
      repData2 <- matrix(repData2, "ncol"=n2Plan, "nrow"=mIter)
    } else {
      sampleSizeRatio <- n2Plan/n1Plan

      repData1 <- stats::rnorm("n"=n1Plan*mIter, "mean"=muGlobal + deltaTrue*sigmaTrue/2, "sd"=sigmaTrue)
      repData1 <- matrix(repData1, "ncol"=n1Plan, "nrow"=mIter)
      repData2 <- stats::rnorm("n"=n2Plan*mIter, "mean"=muGlobal - deltaTrue*sigmaTrue/2, "sd"=sigmaTrue)
      repData2 <- matrix(repData2, "ncol"=n2Plan, "nrow"=mIter)
    }
  }

  if (safeOptioStop) {
    n1Samples <- seq.int(lowN, n1Plan)

    if (is.null(n2Plan)) {
      n2Samples <- NULL
    } else {
      n2Samples <- ceiling(sampleSizeRatio*n1Samples)
    }

    if (pb)
      pbSafe <- utils::txtProgressBar(style=1, title="Safe optional stopping")


    for (iter in seq.int(mIter)) {
      subData1 <- repData1[iter, ]
      subData2 <- repData2[iter, ]

      someT <- unname(stats::t.test("x"=subData1, "y"=subData2, "alternative"=alternative,
                                    "var.equal"=TRUE, "paired"=paired)[["statistic"]])
      someS <- safeTTestStat("t"=someT, "deltaS"=deltaS, "n1"=n1Plan, "n2"=n2Plan, "alternative"=alternative,
                             "paired"=paired)

      if (someS >= 1/alpha)
        safeDecisionAtN[iter] <- 1

      for (k in seq_along(n1Samples)) {

        # TODO(Alexander): Perhaps replace by custom t computing to speed things up
        #
        someT <- unname(stats::t.test("x"=subData1[seq.int(n1Samples[k])], "y"=subData2[seq.int(n2Samples[k])],
                                      "alternative"=alternative, "var.equal"=TRUE, "paired"=paired)[["statistic"]])

        someS <- safeTTestStat("n1"=n1Samples[k], "n2"=n2Samples[k], "t"=someT, "deltaS"=deltaS,
                               "alternative"=alternative, "paired"=paired)

        if (someS >= 1/alpha) {
          allSafeN[iter] <- n1Samples[k]
          allSafeDecisions[iter] <- 1

          break()
        }
      } # End loop lowN to n1Plan

      if (pb)
        utils::setTxtProgressBar(pbSafe, value=iter/mIter, title="Experiments")

    } # End iterations

    if (pb)
      close(pbSafe)

    safeSim <- list(powerOptioStop=mean(allSafeDecisions),
                    powerAtN1Plan=mean(safeDecisionAtN),
                    nMean=mean(allSafeN),
                    probLessNDesign=mean(allSafeN < n1Plan),
                    lowN=min(allSafeN)
    )

    if (logging) {
      safeSim[["allN"]] <- allSafeN
      safeSim[["allSafeDecisions"]] <- allSafeDecisions
      safeSim[["allRejectedN"]] <- allSafeN[-which(allSafeN*allSafeDecisions==0)]
    }

    if (!is.null(n1PlanFreq))
      safeSim[["probLeqN1PlanFreq"]] <- mean(allSafeN <= n1PlanFreq)

    result[["safeSim"]] <- safeSim
  }

  if (freqOptioStop) {
    # Note(Alexander): Adjust data set
    #
    if (is.null(n2Plan)) {
      if (n1PlanFreq < n1Plan) {
        repData1 <- repData1[, seq.int(n1PlanFreq)]
      }

      if (n1PlanFreq > n1Plan) {
        n1Diff <- n1PlanFreq - n1Plan
        repData1Extra <- stats::rnorm("n"=n1Diff*mIter, "mean"=deltaTrue*sigmaTrue, "sd"=sigmaTrue)
        repData1Extra <- matrix(repData1Extra, "ncol"=n1Diff, "nrow"=mIter)
        repData1 <- cbind(repData1, repData1Extra)
      }
    } else {
      # Note(Alexander): Two-sample case
      if (n1PlanFreq < n1Plan) {
        repData1 <- repData1[, seq.int(n1PlanFreq)]
      } else if (n1PlanFreq > n1Plan) {
        n1Diff <- n1PlanFreq - n1Plan

        if (paired) {
          repData1Extra <- stats::rnorm("n"=n1Diff*mIter, "mean"=muGlobal + deltaTrue*sigmaTrue/sqrt(2), "sd"=sigmaTrue)
          repData1Extra <- matrix(repData1Extra, "ncol"=n1Diff, "nrow"=mIter)
          repData1 <- cbind(repData1, repData1Extra)
        } else {
          repData1Extra <- stats::rnorm("n"=n1Diff*mIter, "mean"=muGlobal + deltaTrue*sigmaTrue/2, "sd"=sigmaTrue)
          repData1Extra <- matrix(repData1Extra, "ncol"=n1Diff, "nrow"=mIter)
          repData1 <- cbind(repData1, repData1Extra)
        }
      }

      if (n2PlanFreq < n2Plan) {
        repData2 <- repData2[, seq.int(n2PlanFreq)]
      } else if (n2PlanFreq > n2Plan) {
        n2Diff <- n2PlanFreq - n2Plan

        if (paired) {
          repData2Extra <- stats::rnorm("n"=n2Diff*mIter, "mean"=muGlobal - deltaTrue*sigmaTrue/sqrt(2), "sd"=sigmaTrue)
          repData2Extra <- matrix(repData2Extra, "ncol"=n2Diff, "nrow"=mIter)
          repData2 <- cbind(repData2, repData2Extra)
        } else {
          repData2Extra <- stats::rnorm("n"=n2Diff*mIter, "mean"=muGlobal - deltaTrue*sigmaTrue/2, "sd"=sigmaTrue)
          repData2Extra <- matrix(repData2Extra, "ncol"=n2Diff, "nrow"=mIter)
          repData2 <- cbind(repData2, repData2Extra)
        }
      }
    }

    n1Samples <- seq.int(lowN, n1PlanFreq)

    if (is.null(n2PlanFreq)) {
      n2Samples <- NULL
    } else {
      n2Samples <- ceiling(sampleSizeRatio*n1Samples)
    }

    if (pb)
      pbFreq <- utils::txtProgressBar(style=1, title="Frequentist optional stopping")

    for (iter in seq.int(mIter)) {
      subData1 <- repData1[iter, ]
      subData2 <- repData2[iter, ]
      someP <- stats::t.test("x"=subData1, "y"=subData2, "alternative"=alternative,
                             "var.equal"=TRUE, "paired"=paired)[["p.value"]]

      if (someP < alpha)
        freqDecisionAtN[iter] <- 1

      for (k in seq_along(n1Samples)) {
        someP <- stats::t.test("x"=subData1[seq.int(n1Samples[k])], "y"=subData2[seq.int(n2Samples[k])],
                               "alternative"=alternative, "var.equal"=TRUE, "paired"=paired)[["p.value"]]

        if (someP < alpha) {
          allFreqN[iter] <- n1Samples[k]
          allFreqDecisions[iter] <- 1
          break()
        }
      } # End loop lowN to n1Plan

      if (pb)
        utils::setTxtProgressBar(pbFreq, value=iter/mIter, title="Experiments")
    } # End iterations

    if (pb)
      close(pbFreq)

    freqSim <- list(powerOptioStop=mean(allFreqDecisions),
                    powerAtN1Plan=mean(freqDecisionAtN),
                    nMean=mean(allFreqN),
                    probLessNDesign=mean(allFreqN < n1PlanFreq),
                    lowN=min(allFreqN)
    )

    if (logging)
      freqSim[["allN"]] <- allFreqN

    if (safeOptioStop)
      freqSim[["probLeqNSafe"]] <- mean(allFreqN <= n1Plan)

    result[["freqSim"]] <- freqSim
  }
  return(result)
}
