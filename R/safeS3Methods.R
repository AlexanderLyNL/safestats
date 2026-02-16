# S3 helpers ---------
#
#' Gets the Label of the Test
#'
#' Helper function that outputs the name of the analysis.
#'
#' @param testType A character string. For the t-tests: "oneSample", "paired", "twoSample".
#' @param testName The name of the analysis that is performed such as "Z-Test",
#' and "Test of Two Proportions".
#'
#' @return Returns a character string with the name of the analysis.
getNameTestType <- function(testType, testName) {

  testTypeChar <- switch(testType,
                     "oneSample"="One Sample",
                     "paired"="Paired Sample",
                     "twoSample"="Two Sample",
                     "gLogrank"="Gaussian",
                     "eLogrank"="Exact",
                     "logrank"="",
                     "2x2" = "Test of ")
  analysisName <- paste("Savi", testTypeChar, testName)
  return(analysisName)

  # return(paste(nameChar, testName))
}

#' Construct a savi design object to be set in the design function
#'
#' @inheritParams getNameTestType
#'
#' @return a savi design object
#' @export
#'
#' @examples
#' obj <- constructSaviDesignObj("Z-Test")
constructSaviDesignObj <- function(testName) {
  result <- list(
    "parameter"=NULL, "esMin"=NULL, "alpha"=NULL,
    "alternative"=NULL, "h0"=NULL, "pilot"=FALSE,
    "eType"=NULL, "testType"=NULL, "testName"=NULL,
    "nPlanBatch"=NULL, "esTrue"=NULL,
    "nPlan"=NULL, "nPlanTwoSe"=NULL, "bootObjN1Plan"=NULL,
    "nMean"=NULL, "nMeanTwoSe"=NULL, "bootObjN1Mean"=NULL,
    "power"=NULL, "powerTwoSe"=NULL, "bootObjBeta"=NULL,
    "logImpliedTarget"=NULL, "logImpliedTargetTwoSe"=NULL,
    "bootObjLogImpliedTarget"=NULL,
    "beta"=NULL, "betaTwoSe"=NULL,
    "futilityResult"=NULL, "varEqual"=NULL,
    "samplePaths"=NULL, "breakVector"=NULL, "designScenario"=NULL,
    "call"=NULL, "timeStamp"=Sys.time(), "note"=NULL)

  testSpecificList <- list()

  if (testName=="Z-Test") {
    testSpecificList <- list("sigma"=NULL, "kappa"=NULL, "ratio"=NULL, "testName"=testName)
  } else if (testName=="T-Test") {
    testSpecificList <- list("ratio"=NULL, "testName"=testName)
  } else if (testName=="Logrank") {
    testSpecificList <- list("exact"=NULL)
  }

  result <- utils::modifyList(result, testSpecificList)
  class(result) <- "saviDesign"
  return(result)
}


#' Construct a savi test object to be set in the savi testing function
#'
#' @inheritParams getNameTestType
#'
#' @return a savi test object
#' @export
#'
#' @examples
#' obj <- constructSaviTestObj("Z-Test")
constructSaviTestObj <- function(testName) {
  result <- list(
    "statistic"=NULL, "n"=NULL, "eValue"=NULL,
    "confSeq"=NULL, "estimate"=NULL, "ciValue"=FALSE,
    "dataName"=NULL, "alternative"=NULL,
    "testType"=NULL, "designObj"=NULL, "h0"=NULL,
    "eValueVec"=NULL, "confSeqMatrix"=NULL,
    "eValueApproxError"=NULL, "call"=NULL, "note"=NULL)

  testSpecificList <- list()

  if (testName=="Z-Test") {
    testSpecificList <- list("sigma"=NULL, "testName"=testName)
  } else if (testName=="T-Test") {
    testSpecificList <- list("stderr"=NULL, "testName"=testName)
  } else if (testName=="Logrank") {
    testSpecificList <- list("sumStats"=NULL, "testName"=testName)
  }

  result <- utils::modifyList(result, testSpecificList)
  class(result) <- "saviTest"
  return(result)
}
#' Gets the Label of the Alternative Hypothesis
#'
#' Helper function that outputs the alternative hypothesis of the analysis.
#'
#' @param alternative A character string. "twoSided", "greater", "less".
#' @param testType A character string either "oneSample", "paired", "twoSample", "gLogrank", or "eLogrank".
#' @param h0 the value of the null hypothesis
#' @return Returns a character string with the name of the analysis.
getNameAlternative <- function(alternative=c("twoSided", "greater", "less"), testType, h0=0) {
  alternative <- match.arg(alternative)

  if (testType == "oneSample") {
    trueMeanStatement <- "true mean"
  } else if (testType %in% c("twoSample", "paired")
             && names(h0)!="mu") {
    trueMeanStatement <- paste0("\n", names(h0), " is")
  } else if (testType %in% c("paired", "twoSample")) {
    trueMeanStatement <- "true difference in means ('x' minus 'y') is"
  } else if (testType == "2x2") {
    trueMeanStatement <- "true difference between proportions in group a and b is"
  } else if (testType %in% c("gLogrank", "eLogrank", "logrank")) {
    trueMeanStatement <- "true hazard ratio is"
  }

  nameChar <- paste(trueMeanStatement, switch(alternative,
                                              "twoSided"= paste("not equal to", h0),
                                              "greater"= paste("greater than", h0),
                                              "less"= paste("less than", h0))
  )
  return(nameChar)
}

#' Print Method for Savi Test Objects
#'
#' Printing objects of class 'saviTest' modelled after \code{\link[stats]{print.htest}()}.
#'
#' @param x a saviTest object.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to strwrap for displaying the method components.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return No returned value, called for side effects.
#' @export
#'
#' @examples
#' saviTTest(rnorm(19))
print.saviTest <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  designObj <- x[["designObj"]]

  if (is.null(designObj)) {
    print.default(x)
    return()
  }

  if (!is.null(x[["testType"]]) && x[["testType"]] != designObj[["testType"]])
    designObj[["testType"]] <- x[["testType"]]

  testType <- designObj[["testType"]]

  analysisName <- getNameTestType("testType"=testType,
                                  "testName"=designObj[["testName"]])
  alternativeName <- getNameAlternative("alternative"=designObj[["alternative"]],
                                        "testType"=testType, "h0"=designObj[["h0"]])

  cat("\n")
  cat(strwrap(analysisName, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x[["dataName"]], ". ", sep="")

  nObs <- x[["n"]]

  if (!is.null(nObs)) {
    out <- character()
    out <- c(out, paste(names(nObs), "=", format(nObs, digits = max(1L, digits - 2L))))
    cat(paste(out, collapse = ", "), sep="\n")
  }

  estimate <- x[["estimate"]]

  if (!is.null(estimate)) {
    out <- character()
    out <- c(out, paste(names(estimate), "=", format(estimate, digits = max(1L, digits - 2L))))
    cat(paste0("estimates: ", paste(out, collapse = ", "), sep="\n"))
  }

  ciValue <- x[["ciValue"]]
  confSeq <- x[["confSeq"]]

  if (!is.null(confSeq) && !is.null(ciValue)) {
    cat(format(100*(ciValue)), " percent confidence sequence:\n",
        " ", paste(format(x[["confSeq"]][1:2], digits = digits),
                   collapse = " "), "\n", sep = "")
  }
  cat("\n")

  # if (!is.null(ciValue) && !is.null(confSeq)) {
  #   out <- character()
  #   out <- c(out, paste(names(estimate), "=", format(estimate, digits = max(1L, digits - 2L))))
  #   cat(paste0("estimates: ", paste(out, collapse = ", "), sep="\n"))
  # }
  # cat("\n")

  statValue <- x[["statistic"]]
  parameter <- designObj[["parameter"]]
  eValue <- x[["eValue"]]
  eValueApproxError <- x[["eValueApproxError"]]

  alphaString <- format(designObj[["alpha"]], digits = max(1L, digits - 2L))
  eValueString <- format(eValue, digits = max(1L, digits - 2L))
  eThresholdString <- format(1/designObj[["alpha"]], digits = max(1L, digits - 2L))

  out <- character()

  if (!is.null(statValue))
    out <- c(out, paste(names(statValue), "=", format(statValue, digits = max(1L, digits - 2L))))

  out <- c(out, paste(names(parameter), "=", format(parameter, digits = max(1L, digits - 2L))))

  if (!is.null(designObj[["eType"]]))
    out <- c(out, paste("type", "=", designObj[["eType"]]))

  cat(paste0("test: ", paste(out, collapse = ", "), sep="\n"))
  cat("e-value =", eValueString, ">= 1/alpha =", eThresholdString, ":",
      eValue > 1/designObj[["alpha"]])
  cat("\n")
  if (!is.null(eValueApproxError))
    cat("e-value: approx. error = ",
        format(eValueApproxError/eValue*100, digits = max(1L, digits - 2L)),
        "%. ", sep="")
  cat("\n")
  cat("alternative hypothesis:", alternativeName, "\n")

  # if (!is.null(x$conf.int)) {
  #   cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n",
  #       " ", paste(format(x$conf.int[1:2], digits = digits),
  #                  collapse = " "), "\n", sep = "")
  # }

  cat("\n")

  cat("design: ")

  if (designObj[["pilot"]]) {
    cat("the pilot test is based on an exploratory alpha =", alphaString)
    cat("\n")
  } else {
    cat("the test was designed with alpha =", alphaString)
    cat("\n")

    nPlan <- designObj[["nPlan"]]

    if (!is.null(nPlan)) {
      out <- paste(names(nPlan), "=", nPlan)
      cat(paste0("for experiments with ", paste(out, collapse = ", "), sep="\n"))
    }

    betaValue <- designObj[["beta"]]

    if (!is.null(betaValue)) {
      betaString <- format(designObj[["beta"]], digits = max(1L, digits - 2L))
      powerString <- format(1-designObj[["beta"]], digits = max(1L, digits - 2L))

      cat("to guarantee a power = ", powerString,
          " (beta = ", betaString, ")", sep="")
      cat("\n")
    }

    esMin <- designObj[["esMin"]]

    if (!is.null(esMin)) {
      out <- paste0("minimal relevant ", names(esMin), " = ", format(esMin, digits = max(1L, digits - 2L)),
                    " (", designObj[["alternative"]], ")")
      cat("for", out, "\n")
    }

    # TODO(Alexander): Add this note?
    # nDiff <- nPlan - nObs
    # moreNIndex <- which(nDiff > 0)
    #
    # if (length(moreNIndex) > 0) {
    #   nDiffNames <- paste(names(nPlan), "-", names(nObs))
    #   out <- paste(nDiffNames[moreNIndex], "=", nDiff[moreNIndex])
    #   cat(paste0("Note: ", paste(out, collapse = ", "), ".", sep="\n"))
    # }
  }
}


#' Print Method for Savi Test Objects
#'
#' Printing objects of class 'saviTest' modelled after \code{\link[stats]{print.power.htest}()}.
#'
#' @inheritParams print.saviTest
#'
#' @return No returned value, called for side effects.
#' @export
#'
#' @examples
#' designSaviZ(meanDiffMin=0.5)
#' designSaviT(deltaMin=0.5)
#' designSaviLogrank(hrMin=0.7)
print.saviDesign <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  designObj <- x

  if (is.null(designObj[["parameter"]])) {
    print.default(x)
    return()
  }

  testName <- designObj[["testName"]]
  testType <- designObj[["testType"]]

  futility <- if (!is.null(x[["futilityResult"]])) TRUE else FALSE

  note <- designObj[["note"]]

  tempName <- getNameTestType("testType"=testType, "testName"=testName)

  extraChar <- if (futility) " + Futility" else ""

  analysisName <- paste0(tempName, extraChar, " Design")

  cat("\n")
  cat(strwrap(analysisName, prefix = prefix), sep = "\n")
  cat("\n")

  displayList <- list()

  if (futility) {
    designObj[["decision rule 1"]] <- 1/designObj[["alpha"]]
    designObj[["decision rule 2"]] <- designObj[["futilityResult"]][["beta"]]
  } else {
    designObj[["decision rule"]] <- 1/designObj[["alpha"]]
  }

  for (item in c("nPlan", "nEvents", "nMean", "esMin", "alternative",
                 "alternativeRestriction", "power", "beta",
                 "parameter", "alpha",
                 "decision rule", "decision rule 1", "decision rule 2",
                 "logImpliedTarget", "eType")) {
    itemValue <- designObj[[item]]
    itemValueString <- format(itemValue, digits=digits)

    if (!is.null(itemValue)) {
      if (item %in% c("nPlan", "nMean")) {
        itemNeem <- paste0(item, "TwoSe")

        itemTwoSe <- designObj[[itemNeem]]

        if (!is.null(itemTwoSe)) {
          tempNeem <- names(designObj[[item]])

          for (i in seq_along(itemValue)) {
            if (i==1) {
              itemValueString <- paste0(format(itemValue[i], digits=digits), "\U00B1",
                                        format(itemTwoSe[i], digits=digits))
            } else {
              itemValueString <- paste(itemValueString,
                                       paste0(format(itemValue[i], digits=digits), "\U00B1",
                                              format(itemTwoSe[i], digits=digits)),
                                       sep=", ")
            }
          }
          tempNeem <- paste0(names(designObj[[item]]), "\U00B1", "2se")
          displayList[[paste(tempNeem, collapse=", ")]] <- itemValueString
        } else {
          tempNeem <- names(designObj[[item]])
          displayList[[paste(tempNeem, collapse=", ")]] <- itemValue
        }
      } else if (item=="power") {
        powerTwoSe <- designObj[["powerTwoSe"]]
        itemValueString <- format(itemValue, digits=digits)

        if (!is.null(powerTwoSe)) {
          displayList[[paste0("power", "\U00B1", "2se")]] <-
            paste0(itemValueString, "\U00B1",format(powerTwoSe, digits=digits))
        } else {
          displayList[["power"]] <- itemValueString
        }
      } else if (item=="beta") {
        betaTwoSe <- designObj[["betaTwoSe"]]
        itemValueString <- format(1-itemValue, digits=digits)

        if (!is.null(betaTwoSe)) {
          displayList[[paste0("power: (1 - beta)", "\U00B1", "2se")]] <-
            paste0(itemValueString, "\U00B1",format(betaTwoSe, digits=digits))
        } else {
          displayList[["power: 1 - beta"]] <- itemValueString
        }
      } else if (item=="parameter") {
        displayList[[paste("parameter:", names(designObj[["parameter"]]))]] <- itemValueString
      } else if (item=="decision rule") {
        displayList[["decision rule: e-value >= 1/alpha"]] <- itemValueString
      } else if (item=="decision rule 1") {
        displayList[["decision rule: e-value >= 1/alpha"]] <- itemValueString
      } else if (item=="decision rule 2") {
        displayList[["decision rule: e-value <= beta"]] <- itemValueString
      } else if (item=="logImpliedTarget") {
        tempNeem <- "log(implied target)"
        logImpliedTargetTwoSe <- designObj[["logImpliedTargetTwoSe"]]

        if (!is.null(logImpliedTargetTwoSe)) {
          tempNeem <- paste0(tempNeem, "\U00B1", "2se")
          itemValueString <- paste0(itemValueString, "\U00B1",
                                    format(logImpliedTargetTwoSe, digits=digits))
        }

        displayList[[tempNeem]] <- itemValueString
      } else if (item=="esMin") {
        displayList[[paste("minimal", names(itemValue))]] <- itemValueString
      } else if (item == "alternativeRestriction"){
        displayList[["alternative restriction"]] <- itemValueString
      } else if (item == "eType"){
        displayList[["e-variable type"]] <- itemValueString
      } else {
        displayList[[item]] <- itemValueString
      }
    }
  }

  cat(paste(format(names(displayList), width = 20L, justify = "right"),
            format(displayList, digits = digits), sep = " = "), sep = "\n")

  someTime <- designObj[["timeStamp"]]

  if (!is.null(someTime)) {
    cat("\n")
    cat(paste("Timestamp:", format(someTime, usetz = TRUE)))
  }

  if (!is.null(note)) {
    cat("\n")
    nNotes <- length(note)
    if (nNotes == 1) {
      cat("\n", "Note: ", note, "\n", sep = "")
    } else {
      for (i in 1:nNotes) {
        cat("\n", "Note ", i, ": ", note[i], "\n", sep = "")
      }
    }
  }
  # cat("\n")
}


#' Plots the saviDesign object for designs with sample paths
#'
#' @param x designObj
#' @param main character string for the title of plot
#' @param xlab character string for the x-axis
#' @param ylab character string for the y-axis
#' @param xlim vector of length 2 specifying the
#' range of the x-axis
#' @param ylim vector of length 2 specifying the
#' range of the y-axis
#' @param numSamplePaths integer, number of sample paths to plot
#' @param maxNBins integer, maximum number of bins of the histogram
#' @param wantStepLines logical, if TRUE, then plot the sample paths
#' as step functions (realistic).
#' @param wantQuantiles a vector of numerics between zero and one
#' representing the quantile levels
#' @param border the color of the border around the bars.
#' The default is to use blue
#' @param breaks Break points of the histogram see \code{\link[graphics]{hist}()}
#' @param lwd The line width, a positive number, defaulting to 2.
#' @param pch An integer specifying a symbol.
#' @param histInnerColour A colour to be used to fill the bars.
#' @param col colour of the lines
#' @param colQuant colour of the quantiles
#' @param cex size of the labels and the quantile text
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Nothing it only plots
#' @export
#'
plot.saviDesign <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                            xlim=NULL, ylim=NULL, maxNBins=35,
                            numSamplePaths=100, wantStepLines=FALSE,
                            wantQuantiles=NULL, border="#1F78B4E6",
                            breaks=NULL, lwd=2, pch=15,
                            histInnerColour=col,
                            colQuant="#AA0000",
                            col="#A6CEE380",
                            overColourBorder=border,
                            underColour="#FFB90F86",
                            underColourBorder="#FFB90FCC",
                            histInnerColourContinue="#556B2F4D",
                            borderColourContinue="#556B2FCC",
                            cex=1.3, yLabPAdj=-1,
                            wantNotStoppedHist=FALSE,
                            wantNotStoppedSamplePaths=TRUE,
                            wantLegend=TRUE,
                            legendAdjFut=c(-0.1, 0.5, 1),
                            legendAdj=c(0.2, 0.8),
                            legendCexFactor=0.85,
                            pchColourUnder="#FFB90FCC",
                            pchColourOver="#1F78B4E6",...) {

  designScenario <- x[["designScenario"]]

  if (is.null(designScenario) || !(designScenario %in% c("1a", "2")))
    return()

  alpha <- x[["alpha"]]
  nPlan <- x[["nPlan"]][1]

  if (designScenario=="1a") {
    nPlanBatch <- x[["nPlanBatch"]][1]
    breakVector <- x[["breakVector"]]
    stoppingTimes <- x[["bootObjN1Plan"]][["data"]]
  }

  if (designScenario=="2") {

    if (is.null(x[["nPlanBatch"]]))
      x[["nPlanBatch"]] <- nPlan

    nPlanBatch <- x[["nPlanBatch"]]
    samplePaths <- x[["samplePaths"]]
    mIter <- dim(samplePaths)[1]

    firstPassageTimes <- breakVector <- integer(mIter)

    for (j in 1:mIter) {
      firstPassageTimes[j] <- suppressWarnings(
        min(which(samplePaths[j, ] >= 1/alpha))
      )
    }

    stoppingTimes <- firstPassageTimes
    stoppingTimes[is.infinite(firstPassageTimes)] <- nPlanBatch
    breakVector[is.infinite(firstPassageTimes)] <- 1
    breakVector[Matrix::which(x[["breakVector"]]==-1)] <- -1
  }

  if (numSamplePaths==0) {
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()

    maxStoppingTimes <- max(stoppingTimes)
    minStoppingTimes <- min(stoppingTimes)

    if (is.null(breaks)) {
      breaks <- if (maxStoppingTimes-minStoppingTimes+1 > maxNBins)
        maxNBins
      else
        (minStoppingTimes-1):nPlanBatch
    }

    stopHist <- hist(stoppingTimes, freq=FALSE,
                     breaks=breaks,
                     col=histInnerColour,
                     border=border, lwd=lwd, xaxt="n", yaxt="n",
                     xlab="", ylab="", main=main)

    if (is.null(xlab))
      xlab <- "Stopping times"

    axis(side = 2)
    axis(side = 1)
    axis(side = 1, at=c(0, 2*nPlanBatch))

    if (is.null(ylab))
      ylab <- "Density"

    mtext(ylab, side = 2, line = 4, las = 0, cex = cex, adj=0.5)
    mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)

    if (maxStoppingTimes-minStoppingTimes > maxNBins) {
      firstPassageTimes <- stoppingTimes[!breakVector]
      fptHist <- hist(firstPassageTimes, breaks=stopHist[["breaks"]],
                      plot=FALSE)

      lastIndex <- length(stopHist[["counts"]])

      rect(xleft=stopHist[["breaks"]][lastIndex], ybottom=fptHist[["density"]][lastIndex],
           xright=stopHist[["breaks"]][lastIndex+1],
           ytop=stopHist[["density"]][lastIndex],
           col=col, border=border)
    } else {
      lastIndex <- length(stopHist[["counts"]])
      rejectedAtNBatch <- stopHist[["counts"]][lastIndex]-sum(breakVector)

      totalCount <- sum(stopHist[["counts"]])

      rect(xleft=stopHist[["breaks"]][lastIndex], ybottom=rejectedAtNBatch/totalCount,
           xright=stopHist[["breaks"]][lastIndex+1],
           ytop=stopHist[["density"]][lastIndex],
           col=col,
           border=border)
    }

    lines(c(nPlan, nPlan), c(0, max(stopHist[["counts"]])), lwd=lwd, lty=2)
  }

  betaFutility <- NULL
  futResult <- x[["futilityResult"]]
  histFut <- NULL
  futility <- x[["futility"]]

  if (numSamplePaths > 0) {
    fptAll <- stoppingTimes

    # Stop indexes ----
    #
    indexStopAlt <- Matrix::which(breakVector==0)
    indexStopFut <- Matrix::which(breakVector==-1)
    indexStopNot <- Matrix::which(breakVector==1)

    nAlt <- length(indexStopAlt)
    nFut <- length(indexStopFut)
    nNotStopped <- length(indexStopNot)

    nAll <- length(stoppingTimes)

    if (futility) {
      fptFut <- as.numeric(futResult[["stoppingTimes"]][indexStopFut])
      betaFutility <- futResult[["beta"]]

      fptAll[indexStopFut] <- fptFut
    }


    fptAllFinite <- fptAll

    fptStopped <- fptAll
    fptStopped[indexStopNot] <- Inf

    if (isFALSE(wantNotStoppedHist))
      fptAll <- fptStopped


    # Compute hist  ----
    #
    breaksMin <- min(fptAll)

    nMax <- min(nPlanBatch, max(fptAllFinite))

    if (is.null(breaks)) {
      breaks <- if (nMax-breaksMin+1 > maxNBins) {
        maxNBins
      } else {
        (breaksMin-1):nMax
      }
    }

    histAll <- hist(fptAll, plot=FALSE, breaks=breaks)

    histStopped <- hist(fptStopped, plot=FALSE,
                        breaks=breaks)

    y <- histAll[["density"]]
    nB <- length(histAll[["breaks"]])
    yRange <- range(y, 0)

    yMin <- -1*log(3/(2*alpha))

    if (futility)
      yMin <- min(c(yMin, -1*log(3/(2*betaFutility))))

    if (is.null(ylim))
      ylim <- c(yMin, 2.75*log(1/alpha))

    someConstant <- 0.8*(ylim[2]+log(alpha))/yRange[2]
    textHeightQuant <- (ylim[2]+log(alpha))+log(1/alpha)

    # if (!is.null(wantQuantiles))
    # someConstant <- 0.9*0.8*someConstant

    # if (!is.null(wantQuantiles))
    #   someConstant <- someConstant*0.9

    if (is.null(xlim))
      xlim <- c(0, nPlanBatch)

    # Draw canvas ----
    #
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()

    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
         cex.lab = cex, cex.axis = cex, las = 1, main=main,
         xaxt = "n", yaxt = "n", bty = "n", type = "p", pch = pch,
         bg = "grey", ...)

    lines(x=c(0, 1.5*nPlanBatch), y=c(0, 0), col="darkgrey", lwd=lwd, lty=2)

    lines(x=c(0, 1.5*nPlanBatch), y=log(c(1/alpha, 1/alpha)))

    if (futility)
      lines(x=c(0, 1.5*nPlanBatch), y=log(c(betaFutility, betaFutility)))

    if (futility) {
      yLabs <- c(1e-24, alpha, "1", betaFutility, 1/alpha)
      criticalP <- log(c(1e-24, alpha, 1, betaFutility, 1/alpha))
    } else {
      yLabs <- c(1e-24, alpha, "1", 1/alpha)
      criticalP <- log(c(1e-24, alpha, 1, 1/alpha))
    }

    axis(side = 2, at = c(criticalP), tick = TRUE, las = 2, cex.axis = cex,
         labels = yLabs)
    axis(side = 1)
    axis(side = 1, at=c(0, 2*xlim[2]))

    if (is.null(ylab))
      ylab <- "Evidence"

    mtext(ylab, side = 2, line = 2.5, las = 0, cex = cex, adj=0.25, padj=yLabPAdj)

    if (is.null(xlab))
      xlab <- "Sample size"

    mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)

    # Draw hist -----
    #
    if (nFut > nAlt) {
      histInnerColour1 <- underColour
      borderColour1 <- underColourBorder

      histInnerColour2 <- histInnerColour
      borderColour2 <- border
    } else {
      histInnerColour1 <- histInnerColour
      borderColour1 <- border

      histInnerColour2 <- underColour
      borderColour2 <- underColourBorder
    }

    # TODO(Alexander): Here perhaps if for show stopped
    #
    #       BUT THEN ALSO CHECK hist all exists or not
    #
    if (wantNotStoppedHist) {
      rect(histAll[["breaks"]][-nB]+0.5, log(1/alpha),
           histAll[["breaks"]][-1L]+0.5, someConstant*y+log(1/alpha),
           col = histInnerColourContinue,
           border = borderColourContinue,
           lwd=lwd,
           angle = 45, density = NULL, lty = NULL)
    }

    yStopped <- histStopped[["counts"]]/histAll[["counts"]]*y

    rect(histAll[["breaks"]][-nB]+0.5, log(1/alpha),
         histAll[["breaks"]][-1L]+0.5, someConstant*yStopped+log(1/alpha),
         col = histInnerColour1,
         border = borderColour1,
         lwd=lwd,
         angle = 45, density = NULL, lty = NULL)

    if (futility) {
      if (nFut < nAlt) {
        # Futility histogram
        hist2 <- hist(fptFut, plot=FALSE, breaks=histAll[["breaks"]])
      } else {
        # Alt histogram
        fptAlt <- stoppingTimes[indexStopAlt]
        hist2 <- hist(fptAlt, plot=FALSE, breaks=histAll[["breaks"]])
      }

      y2 <- hist2[["counts"]]/histAll[["counts"]]*y

      rect(hist2[["breaks"]][-nB]+0.5, log(1/alpha),
           hist2[["breaks"]][-1L]+0.5, someConstant*y2+log(1/alpha),
           col = histInnerColour2, border = borderColour2, lwd=lwd,
           angle = 45, density = NULL, lty = NULL)
    }

    # Draw sample paths -----
    samplePaths <- x[["samplePaths"]]

    # Alexander: Just to avoid a matrix turning into a vector
    #
    if (nAlt==1)
      indexStopAlt <- c(indexStopAlt, indexStopAlt)

    if (nFut==1)
      indexStopFut <- c(indexStopFut, indexStopFut)

    if (nNotStopped==1)
      indexStopNot <- c(indexStopNot, indexStopNot)

    stoppedPaths <- samplePaths

    if (futility)
      stoppedPaths[indexStopFut, ] <- x[["futilityResult"]][["samplePaths"]][indexStopFut, ]

    nSamplePathsNot <- min(ceil(nNotStopped/nAll*numSamplePaths), nNotStopped)
    nSamplePathsFut <- min(ceil(nFut/nAll*numSamplePaths), nFut)
    nSamplePathsAlt <- min(ceil(nAlt/nAll*numSamplePaths), nAlt)

    if (nFut < nAlt) {
      nSamplePaths <- c(nSamplePathsNot, nSamplePathsAlt, nSamplePathsFut)
      indexes <- list(indexStopNot, indexStopAlt, indexStopFut)
      underColourTemp <- underColourBorder
      overColourTemp <- col
    } else {
      nSamplePaths <- c(nSamplePathsNot, nSamplePathsFut, nSamplePathsAlt)
      indexes <- list(indexStopNot, indexStopFut, indexStopAlt)
      underColourTemp <- underColourBorder
      overColourTemp <- border
    }

    if (isFALSE(wantNotStoppedSamplePaths)) {
      indexes <- indexes[-1]
      nSamplePaths <- nSamplePaths[-1]
    }

    for (i in seq_along(nSamplePaths)) {

      underColourTemp <- adjustcolor(
        underColourTemp, alpha.f=max(1-nSamplePathsFut/nAll, 0.1))
      continueColourTemp <- adjustcolor(
        histInnerColourContinue, alpha.f=max(1-nSamplePathsNot/nAll, 0.1))
      overColourTemp <- adjustcolor(
        overColourTemp, alpha.f=max(1-nSamplePathsAlt/nAll, 0.1))

      if (nSamplePaths[i]==0)
        next()

      drawSamplePaths("fpt"=fptAllFinite[indexes[[i]]], "n"=nSamplePaths[i],
                      "pathsMatrix"=stoppedPaths[indexes[[i]], ], "alpha"=alpha,
                      "beta"=betaFutility,
                      "underColour"= underColourTemp,
                      "continueColour"=continueColourTemp,
                      "col"=overColourTemp,
                      "pch"=pch, "lwd"=lwd, "wantStepLines"=wantStepLines)
    }

    if (is.null(wantQuantiles) && !is.null(x[["power"]]))
      wantQuantiles <- x[["power"]]

    if (!is.null(wantQuantiles) && !isFALSE(wantQuantiles)) {
      mtext("Quantiles [%]", side=2, col=colQuant, cex=cex, adj=0.5, at=textHeightQuant)

      quants <- round(
        stats::quantile(stoppingTimes, wantQuantiles), 2)

      quantileNames <- round(wantQuantiles*100,2)

      for (i in seq_along(quants)) {
        # if (futility) {
        #
        #   if (nFut < nAlt) {
        #     # Futility histogram
        #     hist2 <- hist(fptFut, plot=FALSE, breaks=histAll[["breaks"]])
        #
        #     #
        #     borderColour1
        #     borderColour2
        #
        #     futQuantName <- round(sum(fptFut <= quants[i])/nAll, 0)
        #
        #     text(quantileNames[i], x=quants[i], y=textHeightQuant, col=overColourBorder, cex=legendCexFactor*cex, pos=2)
        #     text(names(quants[i]), x=quants[i], y=textHeightQuant, col=underColourBorder, cex=legendCexFactor*cex, pos=4)
        #
        #   } else {
        #     # Alt histogram
        #     fptAlt <- stoppingTimes[indexStopAlt]
        #     hist2 <- hist(fptAlt, plot=FALSE, breaks=histAll[["breaks"]])
        #   }
        #
        #
        # } else {
        #   text(quantileNames[i], x=quants[i], y=textHeightQuant, col=overColourBorder, cex=cex)
        # }

        text(quantileNames[i], x=quants[i],
             y=textHeightQuant, col=colQuant, cex=cex)

        # TODO(Alexander): perhaps change the lower position as minimum of -log(1/alpha) and log(beta)
        #
        segments(x0=quants[i], y0=-0.9*log(1/alpha), y1=0.95*textHeightQuant, col=colQuant)
        text(quants[i], x=quants[i], y=-log(1/alpha), col=colQuant, cex=cex)
      }
    }

    if (wantLegend) {
      mtext(paste0("True ", names(x[["esMin"]]), " = ", x[["esTrue"]]),
            col=colQuant, side=3, line = 1, las = 1,
            cex = cex*legendCexFactor, adj=0)

      if (futility) {
        mtext(paste0("Alternative: ", round(nAlt/nAll*100, 1), "%"),
              col=border, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdjFut[1])

        mtext(paste0("Futility: ", round(nFut/nAll*100, 1), "%"),
              col=underColourBorder, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdjFut[2])

        mtext(paste0("No decision: ", round(nNotStopped/nAll*100, 1), "%"),
              col=borderColourContinue, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdjFut[3])
      } else {
        mtext(paste0("Alternative: ", round(nAlt/nAll*100, 1), "%"),
              col=border, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdj[1])
        mtext(paste0("No rejection: ", round(nNotStopped/nAll*100, 1), "%"),
              col=borderColourContinue, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdj[2])
      }



      # legend("topleft",
      #        legend=c(paste("Alternative: ", round(nAlt/nAll*100, 1), "%"),
      #                 paste("Futility: ", round(nFut/nAll*100, 1), "%"),
      #                 paste("Continued: ", round(nNotStopped/nAll*100, 1)), "%"),
      #        col=c(col, underColour, histInnerColourContinue),
      #        lty=1, cex=cex, lwd=lwd, box.lty=box.lty, xpd=TRUE)
    }

  }
}





#' Plots the saviTest object for sequential analyses
#'
#' @inheritParams plot.saviDesign
#' @param fillPlot logical, if TRUE then plot the confidence
#' sequence with a background colour
#' @param switchNFill integer, if \code{is.null(fillPlot)}, then
#' fill if the number of samples is smaller than switchNFill
#' @param logScale logical, if TRUE then plot on the logscale
#' @param switchNLog integer, if \code{is.null(logScale)}, then
#' plot on the log scale if the number of samples is larger
#' than switchNLog
#' @param wantConfSeqPlot logical, if \code{TRUE} then plot the
#' confidence sequence instead of the e-value progression
#' @param add logical, default \code{FALSE} so a new plot is made.
#' If \code{TRUE} and \code{wantConfSeqPlot==FALSE} then adds the e-value
#' progression line to an existing plot. If TRUE and
#' \code{wantConfSeqPlot==TRUE} then adds another anytime-valid confidence
#' sequence.
#' @param lineColour The colour to be used for the e-value progression.
#' @param h0Colour Colour to indicate the null hypothesis.
#' @param col The colour for filling the anytime-valid confidence interval.
#' @param border The colour to draw the border.
#' @param density the density of shading lines, in lines per inch.
#' The default value of \code{NULL} means that no shading lines are drawn.
#' A zero value of density means no shading nor filling whereas
#' negative values and \code{NA} suppress shading (and so allow colour filling).
#' @param angle the slope of shading lines, given as an angle in degrees
#' (anti-clockwise).
#' @param fillOddEven logical controlling the polygon shading mode: see
#' \code{\link[graphics]{polygon}()} for details. Default \code{FALSE}.
#' @param runInt logical, if \code{TRUE} (default), then shows the running
#' intersection of the confidence sequence.
#' @param wantFutility logical, if \code{FALSE}, then don't show the
#' e-values for futility. Default \code{wantFutility==NULL}, if
#' \code{designObj[["futility"]]==TRUE} then futility e-values are shown
#' automatically.
#' @return Returns nothing just plots
#' @export
#'
plot.saviTest <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                          xlim=NULL, ylim=NULL, lwd=3, cex=1.3,
                          fillPlot=NULL, switchNFill=1e4,
                          logScale=NULL, switchNLog=30,
                          h0Colour="darkgrey", lineColour="#1F78B4E6",
                          col="#A6CEE380", border="#1F78B4E6",
                          wantConfSeqPlot=FALSE, add=FALSE,
                          density=NULL, angle=45,
                          xaxt=NULL, yaxt=NULL,
                          fillOddEven=FALSE, runInt=TRUE,
                          wantFutility=NULL,
                          underColour="#FFB90F86",
                          underColourBorder="black",
                          # underColourBorder="#FFB90FCC",
                          pchColourUnder="#FFB90FCC",
                          pchColourOver="#1F78B4E6",
                          ...) {
  eValueVec <- x[["eValueVec"]]
  eValueFutVec <- x[["eValueFutVec"]]

  n1Vec <- x[["n1Vec"]]

  designObj <- x[["designObj"]]

  futility <- FALSE

  if (is.null(wantFutility) && designObj[["futility"]] && !is.null(eValueFutVec)) {
    futility <- TRUE
    betaFutility <- designObj[["futilityResult"]][["beta"]]
  }

  if (is.null(n1Vec) || is.null(eValueVec)) {
    warning("Can't plot. No sequential analysis.")
    return()
  }

  if (is.null(xlab)) {
    xlab <- switch(x[["testName"]],
                   "Z-Test"="n1",
                   "T-Test"="n1",
                   "logrank"="Number of events")
  }

  if (isTRUE(wantConfSeqPlot)) {
    confSeqMatrix <- x[["confSeqMatrix"]]
    h0 <- designObj[["h0"]]

    upperLine <- confSeqMatrix[, 2]
    lowerLine <- confSeqMatrix[, 1]

    upperLineFinite <- upperLine[is.finite(upperLine)]
    lowerLineFinite <- lowerLine[is.finite(lowerLine)]

    maxBound <- max(upperLineFinite, na.rm=TRUE)
    minBound <- min(lowerLineFinite, na.rm=TRUE)

    if (maxBound > 0)
      upperLine[is.infinite(upperLine)] <- 2*maxBound
    else if (maxBound <= 0)
      upperLine[is.infinite(upperLine)] <- 1/2*maxBound

    if (minBound > 0)
      lowerLine[is.infinite(lowerLine)] <- 1/2*minBound
    else if (minBound <= 0)
      lowerLine[is.infinite(lowerLine)] <- 2*minBound

    maxX <- max(n1Vec, na.rm=TRUE)

    if (is.null(fillPlot))
      fillPlot <- if (maxX <= switchNFill) TRUE else FALSE

    if (base::isFALSE(add)) {
      oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes();

      if (is.null(logScale))
        logScale <- if (maxX > switchNLog) TRUE else FALSE

      if (is.null(xlim))
        xlim <- c(0.9, maxX)

      if (isTRUE(logScale)) {
        logPlot <- "x"

        if (xlim[1]==0) xlim[1] <- 0.9

      } else {
        logPlot <- ""
      }

      if (is.null(ylim))
        ylim <- c(minBound, maxBound)

      plot(NULL, xlim=xlim, ylim=ylim,
           type="l", xlab = "", ylab = "", cex.lab = cex,
           cex.axis = cex, xaxt="n", yaxt="n", bty="n", log=logPlot)

      lines(c(1, maxX), c(h0, h0), lwd=lwd, lty=2, col=h0Colour)

      if (is.null(xaxt) || xaxt!="n") axis(1)
      if (is.null(yaxt) || yaxt!="n") axis(2)

      if (is.null(ylab))
        ylab <- switch(x[["testName"]],
                       "Z-Test"="mu",
                       "T-Test"="mu",
                       "logrank"="log(hazard ratio)")

      mtext(ylab, side = 2, line = 4, las = 0, cex = cex, adj=0.5)
      mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)
    }

    if (is.null(fillPlot))
      fillPlot <- if (maxX <= switchNFill) TRUE else FALSE

    if (runInt) {
      upperLine <- makeRunningIntersection(upperLine)
      lowerLine <- makeRunningIntersection(lowerLine,
                                           upper=FALSE)
    }

    if (fillPlot) {
      polygon(c(n1Vec, rev(n1Vec)),
              c(upperLine, rev(lowerLine)),
              col=col, border=border, lwd=lwd,
              density = density, angle=angle,
              fillOddEven=fillOddEven)
    } else {
      lines(n1Vec, upperLine, col=border, lwd=lwd)
      lines(n1Vec, lowerLine, col=border, lwd=lwd)
    }
  }

  # e-value plot
  if (!isTRUE(wantConfSeqPlot)) {
    alpha <- designObj[["alpha"]]
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes();

    nInfinite <- sum(is.infinite(eValueVec))

    if (futility) {
      nInfinite <- nInfinite+sum(is.infinite(eValueFutVec))
    }

    if (nInfinite >= 1) {
      warning("Overflow: E-values infinite, removed for plotting")

      finiteIndex <- which(is.finite(eValueVec))

      if (futility) {
        finiteIndex <- intersect(finiteIndex, which(is.finite(eValueFutVec)))
        eValueFutVec <- eValueFutVec[finiteIndex]
      }

      eValueVec <- eValueVec[finiteIndex]
      n1Vec <- n1Vec[finiteIndex]
    }

    maxEValue <- max(eValueVec, na.rm=TRUE)
    minEValue <- min(eValueVec, na.rm=TRUE)

    if (futility) {
      maxEValue <- max(eValueFutVec, maxEValue, na.rm=TRUE)
      minEValue <- min(eValueFutVec, minEValue, na.rm=TRUE)
    }

    if (!isTRUE(add)) {
      rangeEValue <- maxEValue-minEValue

      if (is.null(logScale)) {
        logScale <- FALSE

        logScale <- if (rangeEValue/(1/alpha-1) > 5) TRUE else FALSE

        if (abs(log(minEValue)) > abs(log(maxEValue)))
          logScale <- TRUE
      }

      logPlot <- if (isTRUE(logScale)) "y" else ""

      maxY <- ceil(max(maxEValue, 1/alpha))

      minY <- minEValue

      if (futility)
        minY <- min(minY, betaFutility)

      if (isFALSE(logScale))
        minY <- 0

      lastIndex <- length(n1Vec)
      maxX <- n1Vec[lastIndex]

      if (is.null(xlim))
        xlim <- c(min(n1Vec), maxX)

      if (is.null(ylim))
        ylim <- c(minY, maxY)

      plot(NULL, xlim=xlim, ylim=ylim,
           type="l", xlab = "", ylab = "", cex.lab = cex,
           cex.axis = cex, xaxt="n", yaxt="n", bty="n", log=logPlot)

      threshLine <- c(1/alpha, 1/alpha)
      unitLine <- c(1, 1)

      lines(c(1, maxX), threshLine, lwd=lwd, lty=2, col="grey40")

      if (maxY/threshLine[1] < 10)
        lines(c(1, maxX), unitLine, lwd=lwd, lty=3, col="grey60")

      if (futility)
        lines(c(1, maxX), c(betaFutility, betaFutility),
              lwd=lwd, lty=2, col=underColourBorder)

      if (is.null(xaxt) || xaxt!="n") axis(1)
      if (is.null(yaxt) || yaxt!="n") axis(2)

      if (is.null(ylab))
        ylab <- "e-value"

      mtext(ylab, side = 2, line = 4, las = 0, cex = cex, adj=0.5)
      mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)
    }

    lines(n1Vec, eValueVec, lwd=lwd, col=lineColour)

    if (futility)
      lines(n1Vec, eValueFutVec, lwd=lwd, col=underColour)
  }
}


#' Helper function to draw sample paths
#'
#' @inheritParams designSaviZ
#' @inheritParams plot.saviDesign
#' @param fpt vector, of first passage times
#' @param n integer, number of paths to draw
#' @param pathsMatrix numeric matrix, representing
#' @param continueColour hex colour for the sample paths
#' @param underColour hex colour for the sample paths
#' @param histInnerColourAll
#'
#' @returns
#' @export
#'
#' @examples
drawSamplePaths <- function(fpt, n, pathsMatrix, alpha, betaFutility,
                            col="#1F78B4E6", pch=15, lwd=2,
                            underColour="#FFB90FCC",
                            continueColour="#556B2F4D",
                            pchColourUnder="#FFB90FCC",
                            pchColourOver="#1F78B4E6",
                            wantStepLines=FALSE) {

  for (i in 1:n) {
    pathStopped <- NULL

    stoppedTime <- fpt[i]
    evidenceLine <- pathsMatrix[i, 1:stoppedTime]

    if (evidenceLine[stoppedTime] >= 1/alpha) {
      evidenceLine[stoppedTime] <- 1/alpha
      someColour <- col
      pathStopped <- "alt"
      pchColour <- pchColourOver
    } else if (!is.null(betaFutility) && evidenceLine[stoppedTime] <= betaFutility) {
      evidenceLine[stoppedTime] <- betaFutility
      someColour <- underColour
      pathStopped <- "fut"
      pchColour <- pchColourUnder
    } else {
      someColour <- continueColour
    }


    if (isTRUE(wantStepLines)) {
      xLine <- c(0, rep(1:stoppedTime, each=2))
      yLine <- c(0, 0, rep(log(evidenceLine), each=2))
      yLine <- yLine[-length(yLine)]
    } else {
      xLine <- 0:stoppedTime
      yLine <- c(0, log(evidenceLine))
    }

    lines(xLine, yLine, col=someColour, lwd=lwd, lty=1)

    if (!is.null(pathStopped) && pathStopped %in% c("alt", "fut"))
      points(stoppedTime, log(evidenceLine[stoppedTime]),
             col=pchColour, pch=pch, lwd=lwd)
  }
}
