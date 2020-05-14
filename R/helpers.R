# Try helper functions -----

#' Tries to Evaluate an Expression and Fails with \code{NA}
#'
#' The evaluation fails with \code{NA} by default, but it is also able to fail with other values.
#'
#' @param expr Expression to be evaluated.
#' @param value Return value if there is an error, default is \code{NA_real_}.
#'
#' @return Returns the evaluation of the expression, or \code{value} if it doesn't work out.
#'
#' @examples
#' safestats:::tryOrFailWithNA(integrate(exp, -Inf, Inf)[["value"]], NA)
#' safestats:::tryOrFailWithNA(integrate(exp, 0, 3)[["value"]], NA)
tryOrFailWithNA <- function(expr, value=NA_real_) {
  tryCatch(
    error=function(cnd) value,
    expr
  )
}

#' Checks Whether a Vector of Object Inherits from the Class "try-error"
#'
#' Checks whether any of the provided objects contains a try error.
#'
#' @param ... objects that need testing.
#'
#' @return Returns \code{TRUE} if there's some object that's a try-error, \code{FALSE} when all objects are
#' not try-errors.
#'
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- try(integrate(exp, -Inf, Inf))
#' isTryError(x, y)
#' isTryError(x, y, z)
isTryError <- function(...) {
  obj <- list(...)
  tryErrorFunc <- function(x){inherits(x, "try-error")}
  result <- purrr::some(obj, .p=tryErrorFunc)
  return(result)
}



#' Rounds a Numeric to At Most 5 Significant Figures
#'
#' Helper function to round a numeric to 5 significant figures.
#'
#' @param num numeric
#'
#' @return number rounded up to 5 decimal places.
#'
#' @examples
#' safestats:::round5(pi)
round5 <- function(num) {
  stopifnot(is.numeric(num))
  round(num, 5)
}

# Plot helper -----
#' Sets 'safestats' Plot Options and Returns the Current Plot Options.
#'
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list with the user specified plot options.
#'
#' @examples
#' oldPar <- safestats:::setSafeStatsPlotOptionsAndReturnOldOnes()
#' graphics::plot(1:10, 1:10)
#' graphics::par(oldPar)
setSafeStatsPlotOptionsAndReturnOldOnes <- function(...) {
  oldPar <- graphics::par(no.readonly = TRUE)
  graphics::par(cex.main=1.5, mar=c(5, 6, 4, 4)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
                font.lab=2, cex.axis=1.3, bty="n", las=1, ...)
  return(oldPar)
}


# Vignette helpers ---------

#' Plots the Histogram of Stopping Times
#'
#' Helper function to display the histogram of stopping times.
#'
#' @param safeSim A safeSim object, returned from \code{\link{replicateTTests}} and
#' \code{\link{simulateSpreadSampleSizeTwoProportions}}.
#' @param nPlan numeric > 0, the planned sample size(s).
#' @param deltaTrue numeric, that represents the true underlying standardised effect size delta.
#' @param showOnlyNRejected logical, when \code{TRUE} discards the cases that did not reject.
#' @param nBin numeric > 0, the minimum number of bins in the histogram.
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
#' simResults <- replicateTTests(nPlan=designObj$nPlan, deltaTrue=1, parameter=designObj$parameter,
#' nPlanFreq=freqObj$nPlan)
#'
#' plotHistogramDistributionStoppingTimes(
#'   simResults$safeSim, nPlan = simResults$nPlan,
#'   deltaTrue = simResults$deltaTrue)
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
  oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
  on.exit(graphics::par(oldPar))
  graphics::hist(dataToPlot,
                 breaks = nStep*seq.int(maxLength),
                 xlim = c(0, max(safeSim[["allN"]])),
                 xlab = "stopping time (n collected)",
                 main = mainTitle,
                 col = "lightgrey", ...)
}


#' Selectively Continue Experiments that Did Not Lead to a Null Rejection for a (Safe) T-Test
#'
#' Helper function used in the vignette.
#'
#' @inheritParams replicateTTests
#' @param oldValues vector of s-values or p-values.
#' @param valuesType character string either "sValues" or "pValues".
#' @param designObj a safeDesign object obtained from \code{\link{designSafeT}}, or \code{NULL}
#' if valuesType equal "pValues".
#' @param oldData a list of matrices with names "dataGroup1" and "dataGroup2".
#' @param n1Extra integer, that defines the additional number of samples of the first group. If
#' \code{NULL} and valuesType equals "sValues", then n1Extra equals \code{designObj$nPlan[1]}.
#' @param n2Extra optional integer, that defines the additional number of samples of the second group.
#' If \code{NULL}, and valuesType equals "sValues", then n2Extra equals \code{designObj$nPlan[2]}.
#' @param moreMainText character, additional remarks in the title of the histogram.
#'
#' @return a list that includes the continued s or p-values based on the combined data, and a list of the
#' combined data.
#' @export
#'
#' @examples
#' alpha <- 0.05
#' mIter <- 1000L
#'
#' designObj <- designSafeT(deltaMin=1, alpha=alpha)
#' oldData <- generateTTestData(nPlan=designObj$nPlan, deltaTrue=0, nsim=mIter, seed=1)
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
      n1Extra <- designObj[["nPlan"]][1]
      n2Extra <- designObj[["nPlan"]][2]
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

  if (is.null(n2Extra) || is.na(n2Extra))
    tempNPlan <- n1Extra
  else
    tempNPlan <- c(n1Extra, n2Extra)

  newData <- generateTTestData("nPlan"=tempNPlan,
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

  oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
  on.exit(graphics::par(oldPar))

  graphics::plot(oldHist, "xlim"=c(minX, maxX), ylim=c(0, yMax), "col"="blue",
                 "density"=20, "angle"=45, "xlab"=xLabText, "main"=mainText)
  graphics::plot(newHist, "add"=TRUE, "col"="red", "density"=20, "angle"=-45)
  graphics::abline("v"=threshValue, "col"="grey", "lwd"=2, "lty"=2)

  result[["newValues"]] <- newValues
  result[["combinedData"]] <-list("dataGroup1"=dataGroup1, "dataGroup2"=dataGroup2)

  return(result)
}


#' Generate Survival Data which Can Be Analysed With the `survival` Package
#'
#'
#' @param nP integer > 0 representing the number of of patients in the placebo group.
#' @param nT integer > 0 representing the number of of patients in the treatment group.
#' @param alpha numeric > 0, representing the shape parameter of the Weibull distribution.
#' If alpha=1, then data are generated from the exponential, i.e., constant hazard. For alpha > 1
#' the hazard increases, if alpha < 1, the hazard decreases.
#' @param lambdaP The (relative) hazard of the placebo group.
#' @param lambdaT The (relative) hazard of the treatment group.
#' @param seed A seed number.
#' @param nDigits numeric, the number of digits to round of the random time to
#' @param startTime numeric, adds this to the random times. Default 1, so the startTime is not 0, which
#' is the start time of \code{\link[stats]{rweibull}}.
#' @param endTime The endtime of the experiment.
#' @param orderTime logical, if \code{TRUE} then put the data set in increasing order
#' @param competeRatio The ratio of the data that is due to competing risk.
#'
#' @return A data set with time, status and group.
#' @export
#'
#' @examples
#' generateSurvData(800, 800, alpha=1, lambdaP=0.008, lambdaT=0.008/2)
generateSurvData <- function(nP, nT, alpha=1, lambdaP, lambdaT, seed=NULL, nDigits=0,
                             startTime=1, endTime=180, orderTime=TRUE, competeRatio=0) {
  stopifnot(competeRatio >=0, competeRatio < 1, is.numeric(startTime), is.numeric(endTime))
  set.seed(seed)
  data <- list()

  if (competeRatio==0) {
    data[["time"]] <- round(c(stats::rweibull("n" = nP, "shape" = alpha,
                                              "scale" = lambdaP^(-1/alpha)),
                              stats::rweibull("n" = nT, "shape" = alpha,
                                              "scale" = lambdaT^(-1/alpha))) + startTime,
                            "digits" = nDigits)
    data[["status"]] <- 2  # 2 is death
    data[["group"]] <- c(rep("P", times = nP), rep("T", times = nT))
    data <- as.data.frame(data)
    data[["status"]][data[["time"]] > endTime] <- 1  # 1 is censored
    data[["time"]][data[["time"]] > endTime] <- endTime
  } else {
    moreNP <- ceiling((1+competeRatio)*nP)
    moreNT <- ceiling((1+competeRatio)*nT)

    data[["time"]] <- round(c(stats::rweibull("n" = moreNP, "shape" = alpha,
                                              "scale" = lambdaP^(-1/alpha)),
                              stats::rweibull("n" = moreNT, "shape" = alpha,
                                              "scale" = lambdaT^(-1/alpha))) + startTime,
                            "digits" = nDigits)

    data[["status"]] <- 2  # 2 is death
    data[["group"]] <- c(rep("P", times = moreNP), rep("T", times = moreNT))
    data <- as.data.frame(data)
    data[["status"]][data[["time"]] > endTime] <- 0  # 0 is censored
    data[["time"]][data[["time"]] > endTime] <- endTime

    indexOfDeaths <- which(data[["status"]]==2)
    totalNDeaths <- length(indexOfDeaths)
    nCompete <- floor(competeRatio*totalNDeaths)

    if (nCompete < 1)
      nCompete <- 1

    indexOfCompete <- sample(indexOfDeaths, nCompete)

    data[["status"]][indexOfCompete] <- 1
    data[["status"]] <- factor(data[["status"]], 0:2, c("censored", "competing", "death"))
  }

  if (orderTime)
    data <- data[order(data[["time"]]), ]

  return(data)
}
