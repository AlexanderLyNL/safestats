# Testing fnts -------

#' Computes E-Values Based on the T-Statistic
#'
#' A summary stats version of \code{\link{safeTTest}()} with the data replaced by t, n1 and n2, and the
#' design object by deltaS.
#'
#' @inheritParams designSafeT
#' @param t numeric that represents the observed t-statistic.
#' @param parameter numeric > 0, the safe test defining parameter.
#' @param n1 integer that represents the size in a one-sample T-test, (n2=\code{NULL}). When n2 is not \code{NULL},
#' this specifies the size of the first sample for a two-sample test.
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
#' \code{NULL} it implies that the t-statistic is based on one-sample.
#' @param tDensity Uses the the representation of the safe T-test as the likelihood ratio of t densities.
#' @param eType character, type of e-value: "eCauchy" (default), "eGauss", or "grow"
#' @inherit safeTTest
#'
#' @return Returns a numeric that represent the e10, that is, the e-value in favour of the alternative over the null
#'
#' @export
#'
#' @examples
#' safeTTestStat(t=1, n1=100, parameter=0.4)
#' safeTTestStat(t=3, n1=100, parameter=0.3)
safeTTestStat <- function(
    t, n1, n2=NULL, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    ...) {

  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  nEff <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1 else (1/n1+1/n2)^(-1)
  nu <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1-1 else n1+n2-2

  result <- suppressWarnings(
    safeTTestStatNEffNu("t"=t, "nEff"=nEff, "nu"=nu,
                        "parameter"=parameter, "alternative"=alternative,
                        "tDensity"=tDensity, "paired"=paired, "eType"=eType,
                        ...)
  )

  return(result)
}


#' SafeTTestStat based on the t-statistic, nEff and nu
#'
#' @rdname safeTTestStat
#' @inheritParams computeConfidenceIntervalT
#'
#' @param nEff numeric > 0, the effective sample size. For one sample test this is just n.
#' @param nu numeric > 0, the degrees of freedom.
#'
#' @export
safeTTestStatNEffNu <- function(
    t, nEff, nu, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    ...) {

  stopifnot(nEff >= 0, nu >= 0)

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  if (nu == 0) {
    if (t==0)
      return(list("eValue"=1))
    else
      return(list("eValue"=1))
  }

  if (eType=="grow") {
    # TODO(Alexander):
    #   One-sided not as stable as two-sided due to hypergeo::genhypergeo for the odd component
    #   1. Use Kummer's transform again (??)
    #   2. Switch to numerical integration. Boundary case
    #
    # safeTTestStat(t=-3.1878, parameter=0.29, n1=315, alternative="greater")
    # safeTTestStat(t=-3.1879, parameter=0.29, n1=315, alternative="greater")
    # safeTTestStat(t=-3.188, parameter=0.29, n1=315, alternative="greater")

    # TODO(Alexander): Remove in v0.9.0
    #

    deltaS <- parameter
    a <- t^2/(nu+t^2)
    expTerm <- exp((a-1)*nEff*deltaS^2/2)

    zeroIndex <- abs(expTerm) < .Machine$double.eps
    eValues <- vector("numeric", length(expTerm))

    zArg <- (-1)*a*nEff*deltaS^2/2
    zArg <- zArg[!zeroIndex]
    # Note(Alexander): This made the vector shorter. Only there where expTerm is non-zero will we evaluate
    # the hypergeometric functions

    aKummerFunction <- Re(hypergeo::genhypergeo(U=-nu/2, L=1/2, zArg))

    if (alternative=="twoSided") {
      eValues[!zeroIndex] <- expTerm[!zeroIndex] * aKummerFunction
    } else {
      bKummerFunction <- exp(lgamma(nu/2+1)-lgamma((nu+1)/2))*sqrt(2*nEff)*deltaS*t/sqrt(t^2+nu)[!zeroIndex] *
        Re(hypergeo::genhypergeo(U=(1-nu)/2, L=3/2, zArg))
      eValues[!zeroIndex] <- expTerm[!zeroIndex]*(aKummerFunction + bKummerFunction)
    }

    if (eValues < 0) {
      warning("Numerical overflow: eValue close to zero. Ratio of t density employed.")
      eValues <- safeTTestStatTDensity("t"=t, "parameter"=parameter, "nu"=nu,
                                       "nEff"=nEff, "alternative"=alternative)
    }

    result <- list("eValue"=eValues)
    return(result)
  } else if (eType=="eCauchy") {
    kappaG <- parameter

    if (alternative=="twoSided") {
      twoSidedCauchyIntegrand <- function(g) {
        exp(-1/2*log(1+nEff*g)+(nu+1)/2*(log(1+t^2/nu)-log(1+t^2/(nu*(1+nEff*g))))
            - 2*log(g) + stats::dgamma(x=1/g, shape=1/2, rate=kappaG^2/2, log=TRUE))
      }

      tempResult <- stats::integrate(twoSidedCauchyIntegrand, 0, Inf)
    } else {
      oneSidedCauchyIntegrand <- function(delta) {
        2/kappaG*exp(
          stats::dt(t, df=nu, ncp=sqrt(nEff)*delta, log=TRUE)
          -stats::dt(t, df=nu, ncp=0, log=TRUE)
          +stats::dt(delta/kappaG, df=1, ncp=0, log=TRUE)
        )
      }

      if (alternative=="greater") {
        tempResult <- stats::integrate(oneSidedCauchyIntegrand, 0, Inf)
      } else if (alternative=="less") {
        tempResult <- stats::integrate(oneSidedCauchyIntegrand, -Inf, 0)
      }
    }

    result <- list("eValue"=tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
    return(result)
  } else if (eType=="eGauss") {
    g <- parameter

    if (alternative=="twoSided") {
      logResult <- -1/2*log(1+nEff*g)+((nu+1)/2)*(log((1+t^2/nu))-log(1+t^2/(nu*(1+nEff*g))))

      return(list("eValue"=exp(logResult)))
    } else {
      oneSidedGaussIntegrand <- function(delta) {
        2*exp(
          stats::dt(t, df=nu, ncp=sqrt(nEff)*delta, log=TRUE)
          -stats::dt(t, df=nu, ncp=0, log=TRUE)
          +stats::dnorm(delta, mean=0, sd=sqrt(g), log=TRUE)
        )
      }

      if (alternative=="greater") {
        tempResult <- stats::integrate(oneSidedGaussIntegrand, 0, Inf)
      } else if (alternative=="less") {
        tempResult <- stats::integrate(oneSidedGaussIntegrand, -Inf, 0)
      }
    }

    result <- list("eValue"=tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
    return(result)
  } else if (eType=="bayarri") {
    stop("Not yet implemented")

    if (alternative=="twoSided") {
      tempResult <- try(log(1-(1+t^2/(nEff^2-1))^((2-nEff)/2)))

      if (is.na(tempResult))
        return(1)

      logResult <- log((nEff-1)/(nEff-2))+1/2*log(1+nEff)-log(t^2)+nEff/2*log(1+t^2/nu)+tempResult
      result <- exp(logResult)
      return(result)
    }
  } else if (eType=="mom") {
    g <- parameter

    momIntegrand <- function(delta) {
      exp(
        stats::dt(t, df=nu, ncp=sqrt(nEff)*delta, log=TRUE)
        -stats::dt(t, df=nu, ncp=0, log=TRUE)
        +stats::dnorm(delta, mean=0, sd=sqrt(g), log=TRUE)
      )*delta^2/g
    }
    tempResult <- stats::integrate(momIntegrand, -Inf, Inf)

    result <- list("eValue"=tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
    return(result)
  }
}


#' safeTTestStat() based on t-densities
#'
#' This is \code{\link{safeTTestStat}()} based on t-densities instead of
#' hypergeometric functions.
#'
#' @inheritParams safeTTest
#' @inheritParams safeTTestStatNEffNu
#' @rdname safeTTestStat
#'
#' @return Returns a numeric that represent the e10, that is, the e-value in favour of the
#' alternative over the null.
#'
safeTTestStatTDensity <- function(t, parameter, nu, nEff,
                                  alternative=c("twoSided", "less", "greater"),
                                  paired=FALSE, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  deltaS <- parameter

  if (alternative=="twoSided") {
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


#' Safe Student's T-Test.
#'
#' A safe T-test adapted from \code{\link[stats]{t.test}()} to perform one and two sample T-tests on vectors of data.
#'
#' @rdname safeTTest
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param designObj an object obtained from \code{\link{designSafeT}()}, or \code{NULL}, when pilot
#' equals  \code{TRUE}.
#' @param paired a logical indicating whether you want a paired T-test.
#' @param varEqual a logical variable indicating whether to treat the two variances as being equal. For
#' the moment, this is always \code{TRUE}.
#' @param ciValue numeric is the ciValue-level of the confidence sequence. Default ciValue=NULL,
#' and ciValue = 1 - alpha
#' @param maxRoot Used to bound the candidate set of width of the confidence interval.
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
#' # Examples taken from stats::t.test
#'
#' # Test without a designObj is not ideal
#' safeTTest(1:10, y = c(7:20))      # e = 70.454 > 20
#'
#' # See ?designSafeT for more info
#' designObj <- designSafeT(deltaMin=0.6, alpha=0.05,
#'                          alternative="twoSided",
#'                          testType="twoSample")
#'
#' safeTTest(1:10, y = c(7:20), designObj=designObj)
#'
#' # Mimicking the stats::t.test interface.
#' # Standard calls use the camelCased version though
#' safe.t.test(1:10, y = c(7:20), designObj=designObj)
#'
#' ## Classical example: Student's sleep data
#' plot(extra ~ group, data = sleep)
#' ## Traditional interface
#' with(sleep, safeTTest(extra[group == 1], extra[group == 2],
#'                       designObj=designObj))
#'
#' ## Formula interface
#' safeTTest(extra ~ group, data = sleep, designObj=designObj)
#'
#' ## Formula interface to one-sample test
#' designObj1 <- designSafeT(deltaMin=0.6, testType="oneSample")
#'
#' safeTTest(extra ~ 1, data = sleep, designObj=designObj1)
#'
#' ## Formula interface to paired test
#' ## The sleep data are actually paired, so could have been in wide format:
#' designObjPaired <- designSafeT(deltaMin=0.6, testType="paired")
#' sleep2 <- reshape(sleep, direction = "wide",
#'                   idvar = "ID", timevar = "group")
#' safeTTest(Pair(extra.1, extra.2) ~ 1, data = sleep2,
#'           designObj=designObjPaired)
safeTTest <- function(x, ...) {
  UseMethod("safeTTest")
}

#' @describeIn safeTTest Default S3 method
#' @export
safeTTest.default <- function(
    x, y=NULL, designObj=NULL, paired=FALSE,
    varEqual=TRUE, ciValue=NULL,
    maxRoot=10, ...) {

  result <- constructSafeTestObj("T-Test")

  ## Def: test type -------
  if (is.null(y)) {
    testType <- "oneSample"
  } else {
    if (paired) {
      testType <- "paired"
    } else {
      testType <- "twoSample"
    }
  }

  ## Check: designObj ----
  if (is.null(designObj)) {
    designObj <- designSafeT(0.5, "eType"="mom",
                             "testType"=testType)
    designObj[["pilot"]] <- TRUE

    warningMessage <- paste("No designObj given. Default test computed based",
                            "on a non-local prior at +1/2 and -1/2.")
    warning(warningMessage)
  }

  if (designObj[["testName"]] != "T-Test")
    warning("The provided design is not constructed for the t-test,",
            "please use designSafeT() instead. The test results might be invalid.")

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  ## Check: Data -----
  #
  if (is.null(y)) {
    ### One-sample -----
    #
    if (isTRUE(paired))
      stop("Data error: Paired analysis requested without specifying the second variable")

    dataName <- deparse1(substitute(x))
    x <- x[!is.na(x)]

    n <- nEff <- n1 <- length(x)
    n2 <- NULL
    nu <- n-1

    meanObs <- estimate <- mean(x)
    sdObs <- stats::sd(x)

    names(estimate) <- "mean of x"
    names(n) <- "n1"
  } else {
    dataName <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))

    if (isTRUE(paired))
      xGoodIndeces <- yGoodIndeces  <-
        stats::complete.cases(x, y)
    else {
      yGoodIndeces <- !is.na(y)
      xGoodIndeces <- !is.na(x)
    }

    x <- x[xGoodIndeces]
    y <- y[yGoodIndeces]

    n1 <- length(x)
    n2 <- length(y)

    ### Paired ----
    #
    if (isTRUE(paired)) {
      if (n1 != n2)
        stop("Data error: Error in complete.cases(x, y): Paired analysis requested, ",
             "but the two samples are not of the same size.")

      nEff <- n1
      nu <- n1-1
      meanObs <- estimate <- mean(x-y)
      sdObs <- stats::sd(x-y)
      names(estimate) <- "mean of the differences"
    } else {
      ## Two-sample ----
      nEff <- (1/n1+1/n2)^(-1)
      nu <- n1+n2-2

      sPooledSquared <- ((n1-1)*stats::var(x)+(n2-1)*stats::var(y))/nu

      sdObs <- sqrt(sPooledSquared)

      estimate <- c(mean(x), mean(y))
      names(estimate) <- c("mean of x", "mean of y")
      meanObs <- estimate[1]-estimate[2]
    }

    n <- c(n1, n2)
    names(n) <- c("n1", "n2")
  }

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]

  if (is.null(ciValue))
    ciValue <- 1-alpha

  if (ciValue < 0 || ciValue > 1)
    stop("Can't make a confidence sequence with ciValue < 0 or ciValue > 1, or alpha < 0 or alpha > 1")

  tStat <- tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sdObs)

  if (is.na(tStat))
    stop("Data error: Could not compute the t-statistic")

  names(tStat) <- "t"

  ### Compute: eValue ----
  #
  testResult <- suppressWarnings(
    safeTTestStat("t"=tStat, "parameter"=designObj[["parameter"]], "n1"=n1,
                  "n2"=n2, "alternative"=alternative, "paired"=paired,
                  "eType"=designObj[["eType"]])

  )


  ### Compute: confSeq ----
  #
  result[["confSeq"]] <- computeConfidenceIntervalT(
    "meanObs"=meanObs, "sdObs"=sdObs,
    "nEff"=nEff, "nu"=nu,
    "parameter"=designObj[["parameter"]],
    "eType"=designObj[["eType"]], "ciValue"=ciValue, "maxRoot"=maxRoot)

  ### Fill: Result -----
  #
  result[["statistic"]] <- tStat
  result[["estimate"]] <- estimate
  result[["stderr"]] <- sdObs/sqrt(nEff)
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["testType"]] <- testType
  result[["n"]] <- n
  result[["ciValue"]] <- ciValue

  result[["eValue"]] <- testResult[["eValue"]]
  result[["eValueApproxError"]] <- testResult[["eValueApproxError"]]

  names(result[["statistic"]]) <- "t"

  return(result)
}


#' @describeIn safeTTest S3 method for class 'formula'
#' @export
#'
safeTTest.formula <- function(
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

  # Note: Prepare calling stats::model.frame instead of safeTTest
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

    tResult <- safeTTest("x"=dataList[[1L]], "y"=dataList[[2L]], ...)

    if (length(tResult[["estimate"]]) == 2L) {
      names(tResult[["estimate"]]) <- paste("mean in group", levels(groupingFactor))
      names(tResult[["designObj"]][["h0"]]) <-
        paste("true difference in means between",
              paste("group", levels(groupingFactor), collapse = " and "))
    }
  } else {
    respVar <- modelFrame[[response]]

    if (inherits(respVar, "Pair")) {
      tResult <- safeTTest("x"=respVar[, 1L], "y"=respVar[, 2L],
                           paired=TRUE, ...)
      firstVar <- substring(dataName,
                            first=6,
                            last=regexpr(",", dataName)-1)
      secondVar <- substring(dataName,
                             first=regexpr(",", dataName)+2,
                             last=regexpr(")", dataName)-1)
      names(tResult[["estimate"]]) <-
        paste("mean difference between", firstVar, "and", secondVar)
      names(tResult[["designObj"]][["h0"]]) <-
        paste("true mean difference between",
              paste(c(firstVar, secondVar), collapse = " and "))
    } else {
      tResult <- safeTTest("x"=respVar, "y"=NULL, ...)
    }
  }

  tResult[["dataName"]] <- dataName
  return(tResult)
}


#' @describeIn safeTTest Alias for safeTTest
#' @export
safe.t.test <- function(x, y=NULL, paired=FALSE, designObj=NULL, varEqual=TRUE,
                        ciValue=NULL, ...) {
  result <- safeTTest("x"=x, "y"=y, "designObj"=designObj,
                      "paired"=paired, "varEqual"=varEqual,
                      ...)

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


#' Helper function: Computes the safe confidence sequence for the mean in a T-test
#'
#' @inheritParams safeTTestStatNEffNu
#' @inheritParams safeTTest
#' @inheritParams designSafeT
#'
#' @param meanObs numeric, the observed mean. For two sample tests this is difference of the means.
#' @param sdObs numeric, the observed standard deviation. For a two-sample test this is the root
#' of the pooled variance.
#'
#' @return numeric vector that contains the upper and lower bound of the safe confidence sequence
#' @export
#'
#' @examples
#' computeConfidenceIntervalT(meanObs=0.3, sdObs=2, nEff=12, nu=11, parameter=0.4)
computeConfidenceIntervalT <- function(
    meanObs, sdObs, nEff, nu, parameter,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    alternative=c("twoSided", "greater", "less"),
    ciValue=0.95, maxRoot=10) {

  eType <- match.arg(eType)
  alternative <- match.arg(alternative)

  trivialConfInt <- c(-Inf, Inf)

  if (eType=="eGauss" && alternative=="twoSided") {
    g <- parameter

    if (nu <= 0) return(trivialConfInt)

    alpha <- 1-ciValue

    numeratorW <- nu*(((1+nEff*g)/alpha^2)^(1/(nu+1))-1)
    denominatorW <- 1-((1+nEff*g)/alpha^2)^(1/(nu+1))/(1+nEff*g)

    W <- numeratorW/denominatorW

    if (W < 0) return(trivialConfInt)

    width <- sdObs/sqrt(nEff)*sqrt(W)
  } else {

    if (alternative=="twoSided") {
      ciLogPenaltyFunc <- function(ciValue) 1/(1-ciValue)
    } else if (alternative %in% c("greater", "less")) {
      ciLogPenaltyFunc <- function(ciValue) 1/(2*(1-ciValue))
    }

    targetFunction <- function(t) {
      safeTTestStatNEffNu("t"=t, "nEff"=nEff, "nu"=nu, "parameter"=parameter,
                          "eType"=eType)$eValue-ciLogPenaltyFunc(ciValue)
    }

    tempResult <- suppressWarnings(
      try(stats::uniroot(targetFunction, c(0, maxRoot)))
    )

    if (isTryError(tempResult))
      stop("Can't compute the width of the interval")

    width <- sdObs/sqrt(nEff)*tempResult$root
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


# Design fnts -------

#' Design a Frequentist T-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I and type II error for the frequentist T-test.
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
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
#' }
#' @export
#'
#' @examples
#' designFreqT(0.5)
designFreqT <- function(deltaMin, alpha=0.05, beta=0.2,
                        alternative=c("twoSided", "greater", "less"),
                        h0=0, testType=c("oneSample", "paired", "twoSample"), ...) {
  stopifnot(alpha > 0, beta > 0, alpha < 1, beta < 1)

  testType <- match.arg(testType)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  alternativeFreq <- switch(alternative,
                            "greater"="one.sided",
                            "less"="one.sided",
                            "twoSided"="two.sided")

  testTypeFreq <- switch(testType,
                         "twoSample"="two.sample",
                         "oneSample"="one.sample",
                         "paired"="paired")

  tempResult <- stats::power.t.test("delta"=deltaMin, "power"=1-beta, "type"=testTypeFreq,
                                    "alternative"=alternativeFreq)

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
#' @param alternative a character string specifying the alternative hypothesis must be one of "twoSided" (default),
#' "greater" or "less".
#' @param nPlan vector of max length 2 representing the planned sample sizes.
#' @param h0 a number indicating the hypothesised true value of the mean under the null. For the moment h0=0.
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
#' is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param parameter optional numeric test defining parameter. Default set to \code{NULL}.
#' For eType=="eCauchy" the numerator is a mixture with meanDiff/sigma mixed
#' over a Cauchy distribution centred at zero and scale kappaG. For eType=="eGauss"
#' the numerator is a mixture with meanDiff/sigma mixed over a Gaussian centred at
#' zero and variance g. For eType=="grow" the safe test is a likelihood ratio of the
#' non-central t-distributions with in the denominator the likelihood with non-centrality
#' parameter set to 0, and in the numerator an average likelihood defined by the likelihood
#' at the non-centrality parameter value deltaS. For the two sided
#' case 1/2 at -deltaS and 1/2 deltaS.
#' @param lowEsTrue numeric, lower bound for the candidate set of the
#' targeted minimal clinically relevant effect size.
#' Design scenario 3: nPlan and beta given, goal find deltaMin.
#' @param highEsTrue numeric, upper bound for the candidate set of the
#' targeted minimal clinically relevant effect size.
#' Design scenario 3: nPlan and beta given, goal find deltaMin.
#' @param nSim integer > 0, the number of simulations needed to compute power or the number of samples paths
#' for the safe z test under continuous monitoring.
#' @param nBoot integer > 0 representing the number of bootstrap samples to assess the accuracy of
#' approximation of the power, the number of samples for the safe z test under continuous monitoring,
#' or for the computation of the logarithm of the implied target.
#' @param eType character one of "eCauchy", "eGauss", "grow". "eCauchy" yields e-values based on
#' a Cauchy mixture, "eGauss" based on a Gaussian/normal mixture, and "grow" based on a mixture of
#' two point masses at the minimal clinically relevant standardised effect size.
#' @param wantSamplePaths logical, if \code{TRUE} then also outputs the sample paths.
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
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on
#'   whether it was a one-sample or a two-sample test.}
#'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
#'   ratio between the planned randomisation of condition 2 over condition 1.}
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
designSafeT <- function(
    deltaMin=NULL, beta=NULL, nPlan=NULL,
    alpha=0.05, h0=0, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    wantSamplePaths=TRUE,
    lowEsTrue=0.01, highEsTrue=3,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=1e3L, ...) {

  stopifnot(alpha > 0, alpha < 1)

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

  result <- constructSafeDesignObj("T-Test")

  if (!is.null(parameter)) {
    if (eType=="grow") {
      parameter <- checkAndReturnsEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="deltaS",
        "alternative"=alternative)
    } else if (eType %in% "eGauss") {
      parameter <- checkAndReturnsEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="g",
        "alternative"=alternative)
    } else if (eType=="eCauchy") {
      parameter <- checkAndReturnsEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="kappaG",
        "alternative"=alternative)
    }
  }

  if (!is.null(deltaMin)) {
    deltaMin <- checkAndReturnsEsMinParameterSide(
      "paramToCheck"=deltaMin, "esMinName"="deltaMin",
      "alternative"=alternative)
  }

  designScenario <- NULL

  tempResult <- list()

  if (!is.null(deltaMin) && !is.null(beta) && is.null(nPlan)) {
    designScenario <- "1a"

    tempResult <- designSafeT1aWantNPlan(
      "deltaMin"=deltaMin, "beta"=beta,
      "alpha"=alpha, "alternative"=alternative,
      "ratio"=ratio, "parameter"=parameter, testType=testType,
      "eType"=eType, "wantSamplePaths"=wantSamplePaths,
      "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot, ...)
  } else if (!is.null(deltaMin) && is.null(beta) && is.null(nPlan)) {
    designScenario <- "1b"

    if (is.null(parameter)) {
      parameter <- switch(eType,
                          "mom"=deltaMin^2/2,
                          "eGauss"=deltaMin^2,
                          "imom"=abs(deltaMin),
                          "eCauchy"=abs(deltaMin),
                          "grow"=deltaMin)
    }

    tempResult <- list("parameter"=parameter, "esMin"=deltaMin)
  } else if (is.null(deltaMin) && is.null(beta) && !is.null(nPlan)) {
    #scenario 1c: only nPlan known, can perform a pilot (no warning though)
    designScenario <- "1c"
    stop("Still need to do Judith's W function")

    tempResult <- list("note"="TODO")

  } else if (!is.null(deltaMin) && is.null(beta) && !is.null(nPlan)) {
    # scenario 2: given effect size and nPlan, calculate power and implied target
    designScenario <- "2"

    tempResult <- designSafeT2WantBeta(
      "deltaMin"=deltaMin, "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative, "testType"=testType,
      "ratio"=ratio, "parameter"=parameter, "eType"=eType,
      "wantSamplePaths"=wantSamplePaths,
      "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot, ...)
  } else if (is.null(deltaMin) && !is.null(beta) && !is.null(nPlan)) {
    designScenario <- "3"

    tempResult <- designSafeT3WantEsMin(
      "beta"=beta, "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative, "testType"=testType,
      "parameter"=parameter, "eType"=eType,
      "lowEsTrue"=lowEsTrue, "highEsTrue"=highEsTrue)
  }

  if (is.null(designScenario)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1.a) non-null deltaMin, non-null beta and NULL nPlan, or \n",
         "(1.b) non-null deltaMin, NULL beta, and NULL nPlan, or \n",
         "(1.c) NULL deltaMin, NULL beta, non-null nPlan, or \n",
         "(2) non-null deltaMin, NULL beta and non-null nPlan, or \n",
         "(3) NULL deltaMin, non-null beta, and non-null nPlan.")
  }

  # Fill and name ----
  result <- utils::modifyList(result, tempResult)

  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["testType"]] <- testType
  result[["ratio"]] <- ratio
  result[["eType"]] <- eType

  ## Name esMin ----
  esMin <- result[["esMin"]]

  if (is.na(esMin))
    esMin <- NULL

  if (!is.null(esMin))
    names(esMin) <- "standardised mean difference"

  result[["esMin"]] <- esMin

  ## Name nPlan ----
  nPlan <- result[["nPlan"]]

  if (!is.null(nPlan)) {
    if (designScenario %in% 2:3) {
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
                               "grow"="deltaS",
                               "bayarri"="kappaB")
  }

  result[["parameter"]] <- parameter

  ## Name h0 -----
  names(h0) <- "mu"
  result[["h0"]] <- h0

  result[["call"]] <- sys.call()
  return(result)
}


#' Helper function to designing a T-test (output nPlan)
#'
#' Finds the parameter and beta when provided with only alpha, esMin, and nPlan
#'
#' @inheritParams designSafeT
#'
#' @return A list with the parameter and the targeted nPlan amongst other items
#' @export
#'
#' @examples
#' designSafeT1aWantNPlan(deltaMin=0.9, beta=0.7, nSim=10)
designSafeT1aWantNPlan <- function(
    deltaMin, beta, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=1e3L, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  samplingResult <- computeNPlanSafeT(
    "deltaTrue"=deltaMin, "beta"=beta, "alpha"=alpha,
    "alternative"=alternative, "ratio"=ratio,
    "parameter"=parameter, "testType"=testType, "eType"=eType,
    "wantSamplePaths"=wantSamplePaths,
    "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot)


  result <- designSafe1aHelper("samplingResult"=samplingResult,
                               "esMin"=deltaMin, "beta"=beta,
                               "ratio"=ratio, "testType"=testType)
  return(result)
}

#' Helper function to designing a T-test (output beta)
#'
#' Finds the parameter and beta when provided with only alpha, esMin, and nPlan
#'
#' @inheritParams designSafeT
#'
#' @return A list with the parameter and beta amongst other items
#' @export
#'
#' @examples
#' designSafeT2WantBeta(deltaMin=0.9, nPlan=7, nSim=10)
designSafeT2WantBeta <- function(
    deltaMin, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=1e3L, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1
  nPlan <- checkAndReturnsNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  samplingResult <- computeBetaSafeT(
    "deltaTrue"=deltaMin, "nPlan"=nPlan, "alpha"=alpha,
    "alternative"=alternative,
    "testType"=testType, "parameter"=parameter,
    "eType"=eType, "wantSamplePaths"=wantSamplePaths,
    "seed"=seed, "nSim"=nSim, "nBoot"=nBoot, "pb"=pb)

  result <- designSafe2Helper("samplingResult"=samplingResult,
                              "esMin"=deltaMin, "nPlan"=nPlan, "ratio"=ratio,
                              "testType"=c("oneSample", "paired","twoSample"))
  return(result)
}

#' Helper function to designing a T-test (output esMin)
#'
#' Finds the parameter and esMin when provided with only alpha, beta, and nPlan
#'
#' @inheritParams designSafeT
#'
#' @return A list with the parameter and the targeted esMin amongst other items
#' @export
#'
#' @examples
#' designSafeT3WantEsMin(beta=0.7, nPlan=10)
designSafeT3WantEsMin <- function(
    beta, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1

  nPlan <- checkAndReturnsNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  result <- list("parameter"=NULL, "esMin"=NULL,
                 "nPlan"=nPlan, "beta"=beta, "ratio"=ratio,
                 "note"=NULL)

  args(computeMinEsBatchSafeT)

  deltaMin <- tryOrFailWithNA(
    computeMinEsBatchSafeT(
      "nPlan"=nPlan, "alpha"=alpha, "beta"=beta,
      "alternative"=alternative, "testType"=testType,
      "parameter"=parameter, "eType"=eType,
      "lowEsTrue"=lowEsTrue, "highEsTrue"=highEsTrue)
  )

  if (is.null(parameter)) {
    parameter <- switch(eType,
                        "mom"=deltaMin^2/2,
                        "eGauss"=deltaMin^2,
                        "imom"=abs(deltaMin),
                        "eCauchy"=abs(deltaMin),
                        "grow"=deltaMin)
  }

  result[["parameter"]] <- parameter
  result[["esMin"]] <- deltaMin

  if (is.na(deltaMin))
    result[["note"]] <- "No deltaMin found for the provided beta and nPlan"
  else
    result[["note"]] <- "The reported deltaMin is based on the batch analysis."

  return(result)
}


# Batch design fnts ------

#' Helper function: Computes the planned sample size for the safe T-test based on the minimal clinically
#' relevant standardised effect size, alpha and beta.
#'
#' @inheritParams designSafeT
#' @inheritParams sampleStoppingTimesSafeT
#'
#' @return a list which contains at least nPlan and the deltaS the parameter that defines the safe test
computeNPlanBatchSafeT <- function(
    deltaTrue, alpha=0.05, beta=0.2,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    parameter=NULL, ratio=1) {

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

  result <- list(nPlan=NULL, "parameter"=parameter)

  n1Plan <- NULL
  n2Plan <- NULL

  n1OverNEffRatio <- if (testType=="twoSample") (1+ratio)/ratio else 1

  if (is.null(parameter)) {
    deltaTrue <- checkAndReturnsEsMinParameterSide(
      "paramToCheck"=deltaTrue, "alternative"=alternative,
      "esMinName"="deltaTrue")

    parameter <- switch(eType,
                        "mom"=deltaTrue^2/2,
                        "eGauss"=deltaTrue^2,
                        "imom"=abs(deltaTrue),
                        "eCauchy"=abs(deltaTrue),
                        "grow"=abs(deltaTrue))
  }

  deltaTrue <- abs(deltaTrue)

  # Sample size of greater sided Z-test as lower/upper bound for
  # the candidate set of nEff
  qB <- stats::qnorm(beta)

  nTemp <- exp(-2*log(deltaTrue))*
    (2*qB^2 - 2*qB*sqrt(qB^2+2*log(1/alpha))
     +2*log(1/alpha))

  tempAlternative <- switch(alternative,
                            "twoSided"="twoSided",
                            "greater"="greater",
                            "less"="greater")

  if (testType=="twoSample") {
    n1Func <- function(nEff) (1+ratio)/ratio*nEff
    n2Func <- function(nEff) (1+ratio)*nEff
    nuFunc <- function(nEff) (1+ratio)^2/ratio*nEff-2
  } else if (testType %in% c("oneSample", "paired")) {
    n1Func <- function(nEff) nEff
    n2Func <- function(nEff) NULL
    nuFunc <- function(nEff) nEff-1
  }

  targetFunction <- function(nEff) {
    safeTTestStat(
      stats::qt("p"=beta, "df"=nuFunc(nEff), "ncp"=sqrt(nEff)*deltaTrue),
      "n1"=n1Func(nEff), "n2"=n2Func(nEff), "parameter"=parameter, "alternative"=tempAlternative,
      "eType"=eType)$eValue-1/alpha
  }

  tempResult <- try(stats::uniroot(targetFunction, interval=c(nTemp/2, 2*nTemp)))

  if (isTryError(tempResult))
    tempResult <- try(stats::uniroot(targetFunction, interval=c(10*nTemp, 50*nTemp)))

  if (isTryError(tempResult))
    stop("Can't compute the batched planned sample size")

  nEff <- tempResult[["root"]]

  # Two-sample/paired stuff
  if (testType == "twoSample") {
    n1Plan <- ceiling(nEff * n1OverNEffRatio)
    n2Plan <- ceiling(nEff * n1OverNEffRatio * ratio)
  } else {
    n1Plan <- ceiling(nEff)
    n2Plan <- if (testType == "paired") n1Plan else NULL
  }

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  if (eType=="grow" && alternative=="less" && parameter > 0)
    parameter <- -parameter

  names(parameter) <- switch(eType,
                             "mom"="gMom",
                             "eGauss"="g",
                             "imom"="tau",
                             "eCauchy"="kappaG",
                             "grow"="deltaS",
                             "bayarri"="kappaB")

  result[["parameter"]] <- parameter

  return(result)
}

#' Computes the smallest mean difference that is detectable with chance 1-beta, for the provided
#' sample size
#'
#' @inheritParams  designSafeT
#'
#' @return numeric > 0 that represents the minimal detectable mean difference
#' @export
#'
#' @examples
#' computeMinEsBatchSafeT(27)
computeMinEsBatchSafeT <- function(
    nPlan, alpha=0.05, beta=0.2,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

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

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  if (eType=="mom") {
    paramFunc <- function(deltaTrue) deltaTrue^2/2
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

  if (testType %in% c("oneSample", "paired")) {
    n1 <- nPlan[1]
    n2 <- NULL
    nu <- n1-1
  } else if (testType=="twoSample") {
    n1 <- nPlan[1]
    n2 <- nPlan[2]
    nu <- n1+n2-2
  }

  targetFunction <- function(deltaTrue) {
    safeTTestStat(
      stats::qt("p"=beta, "df"=nu, "ncp"=sqrt(nEff)*deltaTrue),
      "n1"=n1, "n2"=n2, "parameter"=paramFunc(deltaTrue),
      "alternative"=tempAlternative, "eType"=eType)$eValue-1/alpha
  }

  if (eType=="grow")  {
    gaussResult <- computeMinEsBatchSafeT(
      "nPlan"=nPlan, "alpha"=alpha, "beta"=beta, "alternative"=tempAlternative,
      testType=testType, eType="eGauss")
  }

  tempResult <- try(stats::uniroot(targetFunction, interval=c(lowEsTrue, highEsTrue)))

  result <- tempResult[["root"]]

  if (alternative=="less")
    result <- -result

  return(result)
}

# Sampling functions for design ----

#' Simulate stopping times for the safe T-test
#'
#' @inheritParams designSafeT
#' @inheritParams generateNormalData
#' @param nMax integer > 0, maximum sample size of the (first) sample in each sample path.
#' @param wantEValuesAtNMax logical. If \code{TRUE} then compute eValues at nMax. Default \code{FALSE}.
#' @param wantSamplePaths logical. If \code{TRUE} then output the (stopped) sample paths. Default \code{TRUE}.
#' @param wantSimData logical. If \code{TRUE}, then output the simulated data.
#' @param lowN integer, smallest sample size (of the first group).
#'
#' @return a list with stoppingTimes and breakVector. Entries of breakVector are 0, 1. A 1 represents stopping
#' due to exceeding nMax, and 0 due to 1/alpha threshold crossing, which implies that in corresponding stopping
#' time is Inf.
#'
#' @export
#'
#' @examples
#' sampleStoppingTimesSafeT(0.7, nSim=10, nMax=20)
sampleStoppingTimesSafeT <- function(
    deltaTrue, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, lowN=3L, nMax=1e8L,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    wantEValuesAtNMax=FALSE,
    wantSamplePaths=TRUE, wantSimData=FALSE,
    pb=TRUE, seed=NULL, nSim=1e3L, ...) {

  stopifnot(alpha > 0, alpha <= 1,
            is.finite(nMax),
            is.finite(deltaTrue))

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

  if (is.null(parameter)) {
    deltaTrue <- checkAndReturnsEsMinParameterSide(
      "paramToCheck"=deltaTrue, "alternative"=alternative,
      "esMinName"="deltaTrue")

    parameter <- switch(eType,
                        "mom"=deltaTrue^2/2,
                        "eGauss"=deltaTrue^2,
                        "imom"=abs(deltaTrue),
                        "eCauchy"=abs(deltaTrue),
                        "grow"=deltaTrue)
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

  for (sim in seq_along(result[["stoppingTimes"]])) {
    if (testType %in% c("oneSample", "paired")) {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/n1Vector*cumsum(x1)
      x1SquareVector <- cumsum(x1^2)
      sX1Vector <- sqrt(1/(n1Vector-1)*(x1SquareVector - n1Vector*x1BarVector^2))

      badIndeces <- which(n1Vector-1 <= 0)
      sX1Vector[badIndeces] <- 1

      tValues <- sqrt(nEffVector)*x1BarVector/sX1Vector
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
    }

    if (wantEValuesAtNMax) {
      tempResult <- safeTTestStat("t"=tValues[length(tValues)],
                                  "parameter"=parameter,
                                  "n1"=nMax[1], n2=nMax[2],
                                  "alternative"=alternative, "eType"=eType)
      result[["eValuesAtNMax"]][sim] <- tempResult[["eValue"]]
    }

    for (j in seq_along(n1Vector)) {
      tempResult <- suppressWarnings(
        safeTTestStat("t"=tValues[j], "parameter"=parameter,
                      "n1"=n1Vector[j], "n2"=n2Vector[j],
                      "alternative"=alternative,
                      "eType"=eType)
      )

      evidenceNow <- tempResult[["eValue"]]

      if (wantSamplePaths)
        result[["samplePaths"]][sim, j] <- evidenceNow

      if (evidenceNow > 1/alpha) {
        result[["stoppingTimes"]][sim] <- n1Vector[j]
        result[["eValuesStopped"]][sim] <- evidenceNow

        if (wantSamplePaths) {
          result[["samplePaths"]][sim, j:nMax[1]] <- evidenceNow
        }
        break()
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
      utils::setTxtProgressBar(pbSafe, "value"=sim/nSim, "title"="Trials")
  }

  if (pb)
    close(pbSafe)


  result[["parameter"]] <- parameter
  result[["n1Vector"]] <- n1Vector
  result[["ratio"]] <- ratio

  if (isTRUE(wantSimData))
    result[["simData"]] <- simData

  return(result)
}


#' Helper function: Computes the type II error of the safeTTest based on the minimal clinically relevant
#' standardised mean difference and nPlan.
#'
#' @inheritParams designSafeT
#' @inheritParams sampleStoppingTimesSafeT
#'
#' @return a list which contains at least beta and an adapted bootObject of class
#' \code{\link[boot]{boot}()}.
#' @export
#'
#' @examples
#' computeBetaSafeT(deltaTrue=0.7, 27, nSim=10)
computeBetaSafeT <- function(
    deltaTrue, nPlan, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=1e3L, ...) {

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

  ratio <- if (length(nPlan) == 2) nPlan[2]/nPlan[1] else 1

  if (testType=="twoSample" && length(nPlan)==1) {
    nPlan <- c(nPlan, nPlan)
    warning('testType=="twoSample" specified, but nPlan[2] not provided. nPlan[2] is set to ratio = ', ratio,
            'times nPlan[1] = ', nPlan[2])
  }

  deltaTrue <- checkAndReturnsEsMinParameterSide(
    "paramToCheck"=deltaTrue, "alternative"=alternative,
    "esMinName"="deltaTrue")

  if (is.null(parameter)) {
    parameter <- switch(eType,
                        "mom"=deltaTrue^2/2,
                        "eGauss"=deltaTrue^2,
                        "imom"=abs(deltaTrue),
                        "eCauchy"=abs(deltaTrue),
                        "grow"=deltaTrue)
  }

  samplingResult <- sampleStoppingTimesSafeT(
    "deltaTrue"=deltaTrue, "alpha"=alpha,
    "alternative" = alternative, "testType"=testType,
    "ratio"=ratio, "parameter"=parameter, "nMax"=nPlan,
    "eType"=eType,
    "wantEValuesAtNMax"=TRUE, "wantSamplePaths"=wantSamplePaths,
    "pb"=pb, "seed"=seed, "nSim"=nSim, ...)

  result <- computeBetaBootstrapper(samplingResult=samplingResult,
                                    parameter=parameter, nPlan=nPlan,
                                    nBoot=nBoot)

  return(result)
}


#' Helper function: Computes the planned sample size of the safe T-test based on the
#' minimal clinical relevant standardised mean difference.
#'
#'
#' @inheritParams designSafeT
#' @inheritParams sampleStoppingTimesSafeT
#'
#' @return a list which contains at least nPlan and an adapted bootObject of class  \code{\link[boot]{boot}()}.
#'
#' @export
#'
#' @examples
#' computeNPlanSafeT(0.7, 0.2, nSim=10)
computeNPlanSafeT <- function(
    deltaTrue, beta=0.2, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, nMax=1e8,
    eType=c("mom", "imom", "eCauchy", "eGauss", "grow", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=1e3L, ...) {

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

  deltaTrue <- checkAndReturnsEsMinParameterSide(
    "paramToCheck"=deltaTrue, "alternative"=alternative,
    "esMinName"="deltaTrue")

  tempObj <- computeNPlanBatchSafeT(
    "deltaTrue"=deltaTrue, "alpha"=alpha, "beta"=beta,
    "alternative"=alternative, "testType"=testType,
    "parameter"=parameter, "ratio"=ratio, "eType"=eType)

  nPlanBatch <- tempObj[["nPlan"]]
  parameter <- tempObj[["parameter"]]

  samplingResult <- sampleStoppingTimesSafeT(
    "deltaTrue"=deltaTrue, "alpha"=alpha,
    "alternative" = alternative, "testType"=testType,
    "ratio"=ratio, "parameter"=parameter, "nMax"=nPlanBatch,
    "eType"=eType,
    "wantSamplePaths"=wantSamplePaths,
    "pb"=pb, "seed"=seed, "nSim"=nSim, ...)

  result <- computeNPlanBootstrapper("samplingResult"=samplingResult,
                                     "parameter"=parameter, "beta"=beta,
                                     "nPlanBatch"=nPlanBatch, "nBoot"=nBoot)
  return(result)
}


# Helper fnts ------

#' Computes a Sequence of (Effective) Sample Sizes
#'
#' Helper function that outputs a sequence of sample sizes, effective sample sizes,
#' and the degrees of freedom depending on the type of T-test. Also used for Z-tests.
#'
#'
#' @inheritParams designSafeT
#'
#' @param lowN integer that defines the smallest n of our search space for n.
#' @param highN integer largest sample size of the (first) sample. Default set to 100.
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
#' @inheritParams designSafeT
#' @inheritParams safeTTest
#'
#' @param meanDiffTrue numeric representing the true mean for simulations with a Z-test.
#' Default \code{NULL}
#' @param muGlobal numeric, population grand mean
#' @param sigmaTrue numeric > 0, population standard deviation
#' @param meanDiffTrue numeric, data governing parameter value
#' @param deltaTrue numeric, the value of the true standardised effect size (test-relevant parameter).
#' This argument is used by `designSafeT()` with `deltaTrue <- deltaMin`
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
generateNormalData <- function(nPlan, nSim=1000L,
                               deltaTrue=NULL, muGlobal=0, sigmaTrue=1,
                               paired=FALSE,
                               seed=NULL, meanDiffTrue=NULL) {
  stopifnot(all(nPlan > 0))

  if ((is.null(deltaTrue) && is.null(meanDiffTrue)) || !is.null(deltaTrue) && !is.null(meanDiffTrue))
    stop("Please provide either deltaTrue (T-test), or meanDiffTrue (Z-test).")

  result <- list("dataGroup1"=NULL, "dataGroup2"=NULL)
  set.seed(seed)

  # TODO(Alexander): vector("mode"="list", length=length(nPlan))

  n1Plan <- nPlan[1]

  if (is.null(meanDiffTrue))
    meanDiffTrue <- deltaTrue*sigmaTrue

  if (length(nPlan)==1) {
    dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=meanDiffTrue, "sd"=sigmaTrue)
    dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
    dataGroup2 <- NULL
  } else {
    n2Plan <- nPlan[2]

    if (paired) {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muGlobal + meanDiffTrue/sqrt(2), "sd"=sigmaTrue)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nSim, "mean"=muGlobal - meanDiffTrue/sqrt(2), "sd"=sigmaTrue)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nSim)
    } else {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muGlobal + meanDiffTrue/2, "sd"=sigmaTrue)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nSim, "mean"=muGlobal - meanDiffTrue/2, "sd"=sigmaTrue)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nSim)
    }
  }

  return(list("dataGroup1"=dataGroup1, "dataGroup2"=dataGroup2))
}

