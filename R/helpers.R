#' Try to evaluate an expression, if not fail with NA (default)
#'
#' @param expr Expression to be evaluated
#' @param value Return value if there is an error, default is NA_real_
#'
#' @return Returns the evaluation of the expression, or value if it doesn't work out
#'
#' @examples
#' \dontrun{
#' tryOrFailWithNA(integrate(exp, -Inf, Inf)[["value"]], NA)
#' tryOrFailWithNA(integrate(exp, 0, 3)[["value"]], NA)
#' }
tryOrFailWithNA <- function(expr, value=NA_real_) {
  tryCatch(
    error=function(cnd) value,
    expr
  )
}


#' Ensure the Truth of R Expressions and returns TRUE if the expressions are not met.
#'
#' This is basically stopifnot{base}, but instead of stopping it returns TRUE. The following descriptions is
#' adapted from stopifnot{base}: If any of the expressions in ... are not all valid, then instead of stopping a TRUE is
#' returned and an error message is printed indicating the first of the elements of ... which were not true.
#'
#' @param ... any number of (logical) R expressions, which should evaluate to TRUE
#'
#' @return Returns TRUE if the provided expressions are not met
#'
#' @examples
#' \dontrun{
#' x <- 0
#' msg <- failIfNot(x > 3)
#' print(msg)
#' }
failIfNot <- function (...) {
  # This is equivalent to
  #
  # tryCatch(error=function(cnd){
  #   return(list("failed"=TRUE, "error"=conditionMessage(cnd)))
  # },
  # stopifnot(...)
  # )
  #
  result <- NULL

  ll <- list(...)
  n <- length(ll)

  if (n == 0L) {
    return(result)
  }

  Dparse <- function(call, cutoff = 60L) {
    ch <- deparse(call, width.cutoff = cutoff)
    if (length(ch) > 1L) {
      paste(ch[1L], "....")
    } else {
      ch
    }
  }

  head <- function(x, n = 6L) {
    x[seq_len(if (n < 0L) max(length(x) + n, 0L) else min(n, length(x)))]
  }

  abbrev <- function(ae, n = 3L) {
    paste(c(head(ae, n), if (length(ae) > n) "...."), collapse = "\n  ")
  }

  mc <- match.call()

  for (i in 1L:n) {
    if (!(is.logical(r <- ll[[i]]) && !anyNA(r) && all(r))) {
      cl.i <- mc[[i + 1L]]
      msg <- if (is.call(cl.i) && identical(cl.i[[1]], quote(all.equal)) &&
                 (is.null(ni <- names(cl.i)) || length(cl.i) == 3L ||
                  length(cl.i <- cl.i[!nzchar(ni)]) == 3L)) {
        sprintf(gettext("%s and %s are not equal:\n  %s"),
                Dparse(cl.i[[2]]), Dparse(cl.i[[3]]), abbrev(r))
      } else {
        sprintf(ngettext(length(r), "%s is not TRUE", "%s are not all TRUE"),
                Dparse(cl.i))
      }

      result <- msg
      return(result)
    }
  }
  return(result)
}

# Some/any checks -----------
#'Checks whether some object evaluates to TRUE for a provided criterion function
#'
#' @param ... objects that need testing
#' @param func function used to evaluate the truth
#'
#' @return Returns TRUE if there's some object that evaluates to TRUE according to func, and FALSE if all func
#' evaluations lead to FALSE
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- "a"
#' z <- NA
#'
#' isSome(x, y, z, func=is.numeric)
#' isSome(x, y, z, func=is.na)
#' isSome(x, y, z, func=is.vector)
#' isSome(x, y, func=is.na)
#' }
isSome <- function(..., func) {
  # TODO: Make these to return at first find
  obj <- list(...)
  return( purrr::some(obj, func) )
}

#'Checks whether some object is NA
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's an NA, FALSE when all objects are not NA.
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- "a"
#' z <- NA
#'
#' isSomeNA(x, y)
#' isSomeNA(x, y, z)
#' }
#'
#'
isSomeNA <- function(...) {
  return(isSome(..., "func"=anyNA, recursive=TRUE))
}

#'Checks whether some object is NULL
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's a NULL, FALSE when all objects are not NULL
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- "a"
#' z <- NULL
#'
#' isSomeNull(x, y)
#' isSomeNull(x, y, z)
#' }
isSomeNull <- function(...) {
  return(isSome(..., func=is.null))
}


#'Checks whether some object is infinite
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's infinite, FALSE when all objects are finite
#'
#'@examples
#'\dontrun{
#' x <- 1
#' y <- "a"
#' z <- 10^(1e10)
#'
#' isSomeInfinite(x, y)
#' isSomeInfinite(x, y, z)
#' isSomeInfinite(x, y, z, w)
#' }
isSomeInfinite <- function(...) {
  isSome(..., func=is.infinite)
}

#'Checks whether some object is TRUE
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's TRUE, FALSE when all objects are not TRUE
#' @examples
#'\dontrun{
#' x <- 1
#' y <- "a"
#' z <- TRUE
#'
#' isSomeTrue(x, y)
#' isSomeTrue(x, y, z)
#' }
isSomeTrue <- function(...) {
  isSome(..., func=isTRUE)
}



#' Checks whether some object is a try error.
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's a try-error, FALSE when all objects are not try-errors
#' @export
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- "a"
#' z <- try(integrate(exp, -Inf, Inf))
#' isTryError(x, y)
#' isTryError(x, y, z)
#' }
#'
isTryError <- function(...) {
  return(isSome(..., func=function(x){inherits(x, "try-error")}))
}


# Every/all checks -----------

#' Checks whether all objects evaluates to TRUE for a provided function criterion
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if any objects are a try error, returns FALSE otherwise
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- "a"
#' z <- NA
#'
#' isEvery(x, y, z, func=is.numeric)
#' isEvery(x, y, z, func=is.na)
#' isEvery(x, y, z, func=is.vector)
#' }
isEvery <- function(..., func) {
  # TODO: Make these to return at first find
  obj <- list(...)
  return(purrr::every(obj, func))
}

#' Checks whether all objects are numeric
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if all objects are numeric, returns FALSE otherwise
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- 2
#' z <- "3"
#'
#' isEveryNumeric(x, y)
#' isEveryNumeric(x, y, z)
#' }
isEveryNumeric <- function(...) {
  # TODO: Make these to return at first find
  return(isEvery(..., func=is.numeric))
}

#' Checks whether all objects are finite
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if all objects are finite, returns FALSE otherwise
#'
#' @examples
#' \dontrun{
#' x <- 1
#' y <- 2
#' z <- 10^(1e10)
#' isEveryFinite(x, y)
#' isEveryFinite(x, y, z)
#' }
isEveryFinite <- function(...) {
  isEvery(..., func=is.finite)
}


# Labelling helpers ----------

#' Gets the name of the analysis
#'
#' @param testType A character string. For the t-tests: "oneSampleT", "pairedSampleT", "twoSampleT".
#'
#' @return Returns a character with the name of the analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' getNameTestType("oneSampleT")
#' }
getNameTestType <- function(testType) {
  nameChar <- switch(testType,
                     "oneSampleT"="Safe One Sample T-Test",
                     "pairedSampleT"="Safe Paired Sample T-Test",
                     "twoSampleT"="Safe Two Sample T-Test")
  return(nameChar)
}

#' Gets the name of alternative
#'
#' @param alternative A character string. "two.sided", "greater", "less".
#' @param testType A character string. For the t-tests: "oneSampleT", "pairedSampleT", "twoSampleT".
#' @return Returns a character with the name of the analysis.
#'
#' @examples
#' \dontrun{
#' getNameTestType("oneSampleT")
#' }
getNameAlternative <- function(alternative=c("two.sided", "greater", "less"), testType) {
  alternative <- match.arg(alternative)

  if (testType=="oneSampleT") {
    trueMeanStatement <- "true mean"
  } else if (testType %in% c("pairedSampleT", "twoSampleT")) {
    trueMeanStatement <- "true difference in means ('x' minus 'y') is"
  } else if (testType == "safe2x2_result") {
	    trueMeanStatement <- "true difference between proportions in group a and b is"
  }

  nameChar <- paste(trueMeanStatement, switch(alternative,
                                              "two.sided"= "not equal to 0",
                                              "greater"= "greater than 0",
                                              "less"= "less than 0")
  )
  return(nameChar)
}

#' Rounds a numeric to 5
#'
#' @param num numeric
#'
#' @return number rounded up to 5 decimal places
#'
#' @examples
#' \dontrun{
#' round5(pi)
#' }
round5 <- function(num) {
  stopifnot(is.numeric(num))
  round(num, 5)
}

# Plot helper -----
#' Sets plot options
#' @param ... further arguments to be passed to or from methods.
#'
#' @examples
#' \dontrun{
#' safestats::setSafeStatsPlotOptions()
#' plot(1:10, 1:10)
#' }
setSafeStatsPlotOptions <- function(...) {
  graphics::par(cex.main=1.5, mar=c(5, 6, 4, 4)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
                font.lab=2, cex.axis=1.3, bty="n", las=1, ...)
}

# Vignette helpers ---------

#' Plots the histogram of stopping times
#'
#' @param safeSim A safeSim object, returned from "replicateTTests" and "simulateSpreadSampleSizeTwoProportions"
#' @param nPlan numeric > 0, the planned sample size (for the first sample)
#' @param deltaTrue numeric, that represents the true underlying effect size delta
#' @param showOnlyNRejected logical, when TRUE discards the cases that are
#' @param nBin numeric > 0, the minimum number of bins in the histogram
#' @param ... further arguments to be passed to or from methods.
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' plotHistogramDistributionStoppingTimes(
#'   simResults$safeSim, nPlan = simResults$n1Plan,
#'   deltaTrue = simResults$deltaTrue)
#' }
plotHistogramDistributionStoppingTimes <- function(safeSim, nPlan, deltaTrue, showOnlyNRejected=FALSE, nBin=25L, ...) {
  if(showOnlyNRejected) {
    dataToPlot <- safeSim[["allRejectedN"]]
  } else {
    dataToPlot <- safeSim[["allN"]]
  }

  nStep <- floor(nPlan/nBin)

  if (nStep == 0)
    nStep <- 1

  maxLength <- ceiling(nPlan/nStep)

  mainTitle <- bquote(~"Spread of stopping times when true difference " == .(round(deltaTrue,2)))
  setSafeStatsPlotOptions()
  graphics::hist(dataToPlot,
                 breaks = nStep*seq.int(maxLength),
                 xlim = c(0, max(safeSim[["allN"]])),
                 xlab = "stopping time (n collected)",
                 main = mainTitle,
                 col = "lightgrey", ...)
}


#' Selectively continue experiments that did not lead to a null rejection for a (safe) t-test
#'
#' @inheritParams replicateTTests
#' @param oldValues vector of "sValues" or "pValues"
#' @param valuesType character either "sValues" or "pValues"
#' @param designObj a safeTDesign object, or NULL if valuesType=="pValues"
#' @param oldData a list of matrices with names "dataGroup1" and "dataGroup2"
#' @param n1Extra integer, that defines the additional number of samples of the first group. If NULL and
#' valuesType=="sValues", then n1Extra <- designObj$n1Plan
#' use designSafeT to find this)
#' @param n2Extra optional integer, that defines the additional number of samples of the second group. If NULL, and
#' valuesType=="sValues", then n2Extra <- designObj$n2Plan
#' @param moreMainText character, additional remarks in the title of the histogram
#' @return a list that includes the continued s or p-values based on the combined data, and a list of the combined
#' data
#' @export
#'
#' @examples
#' alpha <- 0.05
#' mIter <- 1000L
#'
#' designObj <- designSafeT(deltaMin=1, alpha=alpha)
#' oldData <- generateTTestData(n1Plan=designObj$n1Plan, deltaTrue=0, nsim=mIter, seed=1)
#'
#' sValues <- vector("numeric", length=mIter)
#'
#' for (i in seq_along(sValues)) {
#'   sValues[i] <- safeTTest(x=oldData$dataGroup1[i, ], designObj=designObj)$sValue
#' }
#' # First run: 7 false null rejections
#' sum(sValues > 1/alpha)
#'
#' continuedSafe <- selectivelyContinueTTestCombineData(
#'   oldValues=sValues, designObj=designObj, oldData=oldData,
#'   deltaTrue=0, seed=2)
#'
#' # Second run: 8 false null rejections
#' sum(continuedSafe$newValues > 1/alpha)
#'
#' for (i in 1:3) {
#'   sValues <- continuedSafe$newValues
#'   oldData <- continuedSafe$combinedData
#'   continuedSafe <- selectivelyContinueTTestCombineData(
#'    oldValues=sValues, designObj=designObj, oldData=oldData,
#'    deltaTrue=0, seed=i+2)
#'   print(paste("Iteration", i+2))
#'   print("Number of false null rejections")
#'   print(sum(continuedSafe$newValues > 1/alpha))
#' }
selectivelyContinueTTestCombineData <- function(oldValues, valuesType=c("sValues", "pValues"), designObj=NULL,
                                                alternative=c("two.sided", "greater", "less"),
                                                oldData, deltaTrue, alpha=NULL,
                                                n1Extra=NULL, n2Extra=NULL, seed=NULL, paired=FALSE,
                                                muGlobal=0, sigmaTrue=1, moreMainText="") {
  valuesType <- match.arg(valuesType)
  alternative <- match.arg(alternative)

  if (valuesType=="pValues") {
    if (is.null(alpha)) {
      stop('For valuesType="pValues" an alpha is required')
    }

    if (is.null(n1Extra) && is.null(n2Extra)) {
      stop("Can't sample extra data without being specified the number of extra samples")
    }
  } else if (valuesType=="sValues") {
    if (is.null(designObj)) {
      stop("Can't continue safe test analysis without a design object")
    }

    alpha <- designObj[["alpha"]]

    if (is.null(n1Extra) && is.null(n2Extra)) {
      n1Extra <- designObj[["n1Plan"]]
      n2Extra <- designObj[["n2Plan"]]
    }
  }

  if (is.null(oldData[["dataGroup1"]]))
    stop("Can't combine data from old experiment without a matrix of old data referred to as oldData$dataGroup1")

  result <- list("oldValues"=oldValues, "newValues"=NULL, "valuesType"=valuesType, "combinedData"=NULL, "deltaTrue"=deltaTrue,
                 "cal"=sys.call())

  if (valuesType=="sValues") {
    notRejectedIndex <- which(oldValues <= 1/alpha)
  } else if (valuesType=="pValues") {
    notRejectedIndex <- which(oldValues >= alpha)
  }

  if (length(notRejectedIndex)==0)
    stop("All experiments led to a null rejection. Nothing to selectively continue.")

  oldDataGroup1 <- oldData[["dataGroup1"]][notRejectedIndex, ]
  oldDataGroup2 <- oldData[["dataGroup2"]][notRejectedIndex, ]

  newData <- generateTTestData("n1Plan"=n1Extra, "n2Plan"=n2Extra,
                               "deltaTrue"=deltaTrue, "nsim"=length(notRejectedIndex),
                               "paired"=paired, "seed"=seed,
                               "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue)

  dataGroup1 <- cbind(oldDataGroup1, newData[["dataGroup1"]])
  dataGroup2 <- cbind(oldDataGroup2, newData[["dataGroup2"]])

  newValues <- vector("numeric", length(notRejectedIndex))

  if (valuesType=="sValues") {
    for (i in seq_along(newValues)) {
      newValues[i] <- try(safeTTest("x"=dataGroup1[i, ], "y"=dataGroup2[i, ],
                                "designObj"=designObj, "alternative"=alternative,
                                "paired"=paired)$sValue)
    }

    minX <- log(min(oldValues, newValues))
    maxX <- log(max(oldValues, newValues))

    yOld <- log(oldValues[notRejectedIndex])
    yNew <- log(newValues)

    xLabText <- "log(sValues)"
    mainText <- paste("Histogram of s-values", moreMainText)

    threshValue <- log(1/alpha)
  } else if (valuesType=="pValues") {
    for (i in seq_along(newValues)) {
      newValues[i] <- t.test("x"=dataGroup1[i, ], "y"=dataGroup2[i, ],
                             "alternative"=alternative, "paired"=paired)$p.value
    }
    minX <- 0
    maxX <- 1

    yOld <- oldValues[notRejectedIndex]
    yNew <- newValues

    xLabText <- "p-values"
    mainText <- "Histogram of p-values"

    threshValue <- alpha
  }

  oldHist <- graphics::hist("x"=yOld, plot=FALSE)
  newHist <- graphics::hist("x"=yNew, plot=FALSE)

  yMax <- max(oldHist[["counts"]], newHist[["counts"]])

  setSafeStatsPlotOptions()
  graphics::plot(oldHist, "xlim"=c(minX, maxX), ylim=c(0, yMax), "col"="blue",
                 "density"=20, "angle"=45, "xlab"=xLabText, "main"=mainText)
  graphics::plot(newHist, "add"=TRUE, "col"="red", "density"=20, "angle"=-45)
  graphics::abline("v"=threshValue, "col"="grey", "lwd"=2, "lty"=2)

  result[["newValues"]] <- newValues
  result[["combinedData"]] <-list("dataGroup1"=dataGroup1, "dataGroup2"=dataGroup2)

  return(result)
}
