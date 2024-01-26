# Vignette helper fnts ------

#' Helper function to demonstrate the selective continuation of paired z-tests
#'
#' @inheritParams safeZTest
#' @inheritParams designSafeZ
#' @inheritParams generateNormalData
#'
#' @param n1New Number of samples in the follow-up study
#' @param eValuesOld Matrix of eValues of the previous studies with nSim
#' number of rows
#' @param eOverOld Vector of 0s and 1s, where 1 indicates that
#' @param testName character string either "Z-Test" or "T-Test".
#' @param trackCrossingOld Vector of same number of columns as eValuesOld
#' providing the number of studies that got rejected at the indeced time
#' @param firstPassageTimeOld Vector indicating the first passage time.
#' If Inf, then not yet passed the threshold of 1/alpha
#' @param eStoppedOld The stopped e-value either at the time it crosses
#' the threshold, or at the previous nPlan
#'
#' @return a list with the new e-values amongst other
#' @export
selectivelyContinueZOrTTestData <- function(
    designObj, deltaTrue=NULL, testName=c("Z-Test", "T-Test"),
    n1New, muGlobal, sigma, meanDiffTrue=NULL, nSim=1e3L,
    eValuesOld, eOverOld,
    trackCrossingOld, firstPassageTimeOld,
    eStoppedOld, seed=NULL) {

  n1Old <- length(trackCrossingOld)
  testName <- match.arg(testName)
  alpha <- designObj[["alpha"]]

  if (testName=="Z-Test")
    deltaTrue <- NULL

  if (testName=="T-Test")
    meanDiffTrue <- NULL

  someData <- generateNormalData(
    c(n1New, n1New), muGlobal=muGlobal,
    nSim=nSim, seed=seed,
    meanDiffTrue=meanDiffTrue, deltaTrue=deltaTrue,
    sigmaTrue=sigma)

  statMatrix <- matrix(nrow=nSim, ncol=n1New)

  n1Vector <- 1:n1New
  nuVector <- n1Vector-1

  rejectedIndex <- which(eOverOld==1)
  notRejectedIndex <- which(eOverOld==0)

  for (sim in notRejectedIndex) {
    dataGroup1 <- someData$dataGroup1[sim, ]
    dataGroup2 <- someData$dataGroup2[sim, ]

    differenceScore <- dataGroup1-dataGroup2

    meanDiffVector <- 1/n1Vector*cumsum(differenceScore)

    if (testName=="Z-Test") {
      # The variance of the sum x + (-y) is the sum of the two variances
      # Thus, 2*sigma^2
      sdMeanDiff <- sqrt(2)*sigma
    } else if (testName=="T-Test") {
      sdMeanDiff <- sqrt(
        1/nuVector*(cumsum(differenceScore^2)-n1Vector*meanDiffVector^2)
      )
    }

    statMatrix[sim, ] <-
      sqrt(n1Vector)*meanDiffVector/sdMeanDiff
  }

  if (testName=="T-Test") {
    statMatrix[, 1] <- 0
  }

  # Here we store all the e-values across the
  # number of simulations (nSim) and time (n1)
  eValuesNew <- matrix(nrow=nSim, ncol=n1New)
  eValuesNew[rejectedIndex, ] <- eStoppedOld[rejectedIndex]

  # This indicates whether a simulation yielded e > 1/alpha
  eOverNew <- eOverOld

  # This indicates the first time an experiment yielded e > 1/alpha
  # Default is Inf, which indicates that the e didn't cross 1/alpha
  firstPassageTime <- firstPassageTimeOld

  # This is the e-value at the end time, or whenever
  # it exceeds the threshold of 1/alpha
  eStopped <- eStoppedOld

  # We only continue the not rejected studies
  for (sim in notRejectedIndex) {
    statVector <- statMatrix[sim, ]

    for (i in 1:n1New) {
      # Note the multiplication with the previous e-value
      # stopped at the previous nPlan.
      #
      # If eStopped[sim] = 10, we now only need an
      # additional eValue of 2 to cross the threshold of
      # 1/alpha = 20
      #

      if (testName=="Z-Test") {
        currentEValue <- safeZTestStat(
          statVector[i], parameter=designObj[["parameter"]],
          n1=n1Vector[i], n2=n1Vector[i],
          paired=TRUE,  sigma=sigma, alternative="greater",
          eType=designObj$eType)$eValue * eStoppedOld[sim]
      } else if (testName=="T-Test") {
        currentEValue <- safeTTestStat(
          statVector[i], parameter=designObj[["parameter"]],
          n1=n1Vector[i], n2=n1Vector[i],
          paired=TRUE,  sigma=sigma, alternative="greater",
          eType=designObj$eType)$eValue * eStoppedOld[sim]
      }

      eValuesNew[sim, i] <- currentEValue

      if (currentEValue > 1/alpha && eOverNew[sim]!=1) {
        eOverNew[sim] <- 1
        firstPassageTime[sim] <- i+n1Old
        eStopped[sim] <- currentEValue

        eValuesNew[sim, i:n1New] <- currentEValue
      }

      if (i==n1New && eOverNew[sim]!=1) {
        eStopped[sim] <- currentEValue
      }
    }
  }

  trackCrossing <- integer(n1Old+n1New)

  # Again we carry over the previous 3% type I error, that is,
  # 15 number of false rejections along
  trackCrossing[1:n1Old] <- trackCrossingOld

  for (i in (n1Old+1):(n1Old+n1New)) {
    trackCrossing[i] <-
      sum(firstPassageTime <= i)
  }

  eRejects <- trackCrossing/nSim

  eValues <- cbind(eValuesOld, eValuesNew)
  result <- list("eValues"=eValues, "eOver"=eOverNew,
                 "extraRejections"=sum(eOverNew)-sum(eOverOld),
                 "trackCrossing"=trackCrossing,
                 "firstPassageTime"=firstPassageTime,
                 "eStopped"=eStopped, "eRejects"=eRejects)

  return(result)
}

# Vignette helpers ---------

#' Plots the Histogram of Stopping Times
#'
#' Helper function to display the histogram of stopping times.
#'
#' @param safeSim A safeSim object
#' @param nPlan numeric > 0, the planned sample size(s).
#' @param deltaTrue numeric, that represents the true underlying standardised effect size delta.
#' @param showOnlyNRejected logical, when \code{TRUE} discards the cases that did not reject.
#' @param nBin numeric > 0, the minimum number of bins in the histogram.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return a histogram object, and called for its side-effect to plot the histogram.
#'
#' @export
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

  mainTitle <- bquote(~"Spread of stopping times when true delta " == .(round(deltaTrue,2)))
  oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
  on.exit(graphics::par(oldPar))
  graphics::hist(dataToPlot,
                 breaks = nStep*seq.int(maxLength),
                 xlim = c(0, max(safeSim[["allN"]])),
                 xlab = "stopping time (n collected)",
                 main = mainTitle,
                 col = "lightgrey", ...)
}
