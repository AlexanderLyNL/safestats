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
  analysisName <- paste("Safe", testTypeChar, testName)
  return(analysisName)

  # return(paste(nameChar, testName))
}

#' Construct a  safe design object to be set in the design function
#'
#' @inheritParams getNameTestType
#'
#' @return a safe design object
#' @export
#'
#' @examples
#' obj <- constructSafeDesignObj("Z-Test")
constructSafeDesignObj <- function(testName) {
  result <- list(
    "parameter"=NULL, "esMin"=NULL, "alpha"=NULL,
    "alternative"=NULL, "h0"=NULL, "pilot"=FALSE,
    "eType"=NULL, "testType"=NULL, "testName"=NULL,
    "nPlanBatch"=NULL,
    "nPlan"=NULL, "nPlanTwoSe"=NULL, "bootObjN1Plan"=NULL,
    "nMean"=NULL, "nMeanTwoSe"=NULL, "bootObjN1Mean"=NULL,
    "beta"=NULL, "betaTwoSe"=NULL, "bootObjBeta"=NULL,
    "logImpliedTarget"=NULL, "logImpliedTargetTwoSe"=NULL,
    "bootObjLogImpliedTarget"=NULL,
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
  class(result) <- "safeDesign"
  return(result)
}


#' Construct a safe test object to be set in the safe testing function
#'
#' @inheritParams getNameTestType
#'
#' @return a safe test object
#' @export
#'
#' @examples
#' obj <- constructSafeTestObj("Z-Test")
constructSafeTestObj <- function(testName) {
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
  class(result) <- "safeTest"
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

#' Print Method for Safe Tests
#'
#' Printing objects of class 'safeTest' modelled after \code{\link[stats]{print.htest}()}.
#'
#' @param x a safeTest object.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to strwrap for displaying the method components.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return No returned value, called for side effects.
#' @export
#'
#' @examples
#' safeTTest(rnorm(19))
print.safeTest <- function (x, digits = getOption("digits"), prefix = "\t", ...) {
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
  cat("e-value =", eValueString, "> 1/alpha =", eThresholdString, ":",
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


#' Print Method for Safe Tests
#'
#' Printing objects of class 'safeTest' modelled after \code{\link[stats]{print.power.htest}()}.
#'
#' @inheritParams print.safeTest
#'
#' @return No returned value, called for side effects.
#' @export
#'
#' @examples
#' designSafeZ(meanDiffMin=0.5)
#' designSafeT(deltaMin=0.5)
#' designSafeLogrank(hrMin=0.7)
print.safeDesign <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  designObj <- x

  if (is.null(designObj[["parameter"]])) {
    print.default(x)
    return()
  }

  testName <- designObj[["testName"]]
  testType <- designObj[["testType"]]

  note <- designObj[["note"]]
  analysisName <- paste(getNameTestType("testType"=testType, "testName"=testName), "Design")

  cat("\n")
  cat(strwrap(analysisName, prefix = prefix), sep = "\n")
  cat("\n")

  designObj[["decision rule"]] <- 1/designObj[["alpha"]]

  displayList <- list()

  for (item in c("nPlan", "nEvents", "nMean", "esMin", "alternative",
                 "alternativeRestriction", "beta", "parameter",
                 "alpha", "decision rule", "logImpliedTarget", "eType")) {
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
        displayList[["decision rule: e-value > 1/alpha"]] <- itemValueString
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


#' Plots the safeDesign object for designs with sample paths
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
#' @param ...
#'
#' @return Nothing it only plots
#' @export
#'
#'
plot.safeDesign <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                            xlim=NULL, ylim=NULL, maxNBins=35,
                            numSamplePaths=100, wantStepLines=FALSE,
                            wantQuantiles=NULL, border="#1F78B4E6",
                            breaks=NULL, lwd=2, pch=15,
                            histInnerColour="#A6CEE380", col="#DAA52066",
                            colQuant="#AA0000", cex=1.3, ...) {

  designScenario <- x[["designScenario"]]

  if (designScenario %in% c("1a", "2")) {
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
          min(which(samplePaths[j, ] > 1/alpha))
        )
      }

      stoppingTimes <- firstPassageTimes
      stoppingTimes[is.infinite(firstPassageTimes)] <- nPlanBatch
      breakVector[is.infinite(firstPassageTimes)] <- 1
    }

    if (numSamplePaths==0) {
      oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()

      maxStoppingTimes <- max(stoppingTimes)
      minStoppingTimes <- min(stoppingTimes)

      if (is.null(breaks)) {
        breaks <- if (maxStoppingTimes-minStoppingTimes > maxNBins)
          maxNBins
        else
          minStoppingTimes:nPlanBatch
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

      lines(c(nPlan, nPlan), c(0, max(stopHist$counts)), lwd=lwd, lty=2)
    }

    if (numSamplePaths > 0) {
      firstPassageTimes <- stoppingTimes

      # These are the sample paths that did not resulted in a rejection
      firstPassageTimes[which(breakVector==1)] <- Inf

      breaksMin <- min(firstPassageTimes)

      if (is.null(breaks)) {
        breaks <- if (nPlanBatch-breaksMin > maxNBins)
          maxNBins
        else
          breaksMin:nPlanBatch
      }

      fptHist <-
        hist(firstPassageTimes, plot=FALSE,
             breaks=breaks)

      oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()

      y <- fptHist[["density"]]
      nB <- length(fptHist$breaks)
      yRange <- range(y, 0)

      if (is.null(ylim))
        ylim <- c(-1*log(3/(2*alpha)), 2.75*log(1/alpha))

      someConstant <- (ylim[2]+log(alpha))/yRange[2]
      textHeightQuant <- (ylim[2]+log(alpha))+log(1/alpha)

      if (!is.null(wantQuantiles))
        someConstant <- 0.8*someConstant

      if (!is.null(wantQuantiles))
        someConstant <- someConstant*0.9

      if (is.null(xlim))
        xlim <- c(0, nPlanBatch)

      plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
           cex.lab = cex, cex.axis = cex, las = 1, main=main,
           xaxt = "n", yaxt = "n", bty = "n", type = "p", pch = pch,
           bg = "grey", ...)

      abline(h = log(1), col = "darkgrey", lwd = lwd, lty = 2)
      abline(h = log(1/alpha))

      criticalP = log(c(alpha, 1, 1/alpha))

      axis(side = 2, at = c(criticalP), tick = TRUE, las = 2, cex.axis = cex,
           labels = c(alpha, "1", 1/alpha))
      axis(side = 1)
      axis(side = 1, at=c(0, 2*xlim[2]))

      if (is.null(ylab))
        ylab <- "Evidence"

      mtext(ylab, side = 2, line = 2.5, las = 0, cex = cex, adj=0.25)

      if (is.null(xlab))
        xlab <- "Sample size"

      mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)


      samplePaths <- x[["samplePaths"]]

      stoppedPaths <- samplePaths[!breakVector, ]

      rect(fptHist$breaks[-nB]+0.5, log(1/alpha),
           fptHist$breaks[-1L]+0.5, someConstant*y+log(1/alpha),
           col = histInnerColour, border = border, lwd=lwd,
           angle = 45, density = NULL, lty = NULL)

      finiteFirstPassageTimes <- firstPassageTimes[!breakVector]

      for (i in 1:numSamplePaths) {
        stoppedTime <- finiteFirstPassageTimes[i]
        evidenceLine <- stoppedPaths[i, 1:stoppedTime]

        if (evidenceLine[stoppedTime] > 1/alpha)
          evidenceLine[stoppedTime] <- 1/alpha

        someScale <- 3*alpha

        evidenceLine <- evidenceLine

        if (isTRUE(wantStepLines)) {
          xLine <- c(0, rep(1:stoppedTime, each=2))
          yLine <- c(0, 0, rep(log(evidenceLine), each=2))
          yLine <- yLine[-length(yLine)]
        } else {
          xLine <- 0:stoppedTime
          yLine <- c(0, log(evidenceLine))
        }

        lines(xLine, yLine, col=col, lwd=lwd, lty=1)

        if (evidenceLine[stoppedTime]==1/alpha)
          points(stoppedTime, log(evidenceLine[stoppedTime]),
                 col=border, pch=pch, lwd=lwd)
      }

      if (!is.null(wantQuantiles)) {
        mtext("quantiles", side=2, col=colQuant, cex=cex, adj=0.5, at=textHeightQuant)

        quants <- stats::quantile(stoppingTimes, wantQuantiles)

        for (i in seq_along(quants)) {
          text(wantQuantiles[i], x=quants[i], y=textHeightQuant, col=colQuant, cex=cex)
          segments(x0=quants[i], y0=-0.9*log(1/alpha), y1=0.95*textHeightQuant, col=colQuant)
          text(quants[i], x=quants[i], y=-log(1/alpha), col=colQuant, cex=cex)
        }
      }

    }
  }
}





#' Plots the safeTest object for sequential analyses
#'
#' @inheritParams plot.safeDesign
#' @param fillPlot logical, if TRUE then plot the confidence
#' sequence with a background colour
#' @param switchNFill integer, if is.null(fillPlot), then
#' fill if the number of samples is smaller than switchNFill
#' @param logScale logical, if TRUE then plot on the logscale
#' @param switchNLog integer, if is.null(logScale), then
#' plot on the log scale if the number of samples is larger
#' than switchNLog
#' @param wantConfSeqPlot logical, if TRUE then plot the
#' confidence sequence instead of the e-value progression
#'
#' @return Returns nothing just plots
#' @export
#'
plot.safeTest <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                          xlim=NULL, ylim=NULL, lwd=3, cex=1.3,
                          fillPlot=NULL, switchNFill=1e4,
                          logScale=NULL, switchNLog=30,
                          lineColour="lemonchiffon4",
                          col="#A6CEE380", border="#1F78B4E6",
                          wantConfSeqPlot=FALSE, ...) {
  eValueVec <- x[["eValueVec"]]
  n1Vec <- x[["n1Vec"]]

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
    h0 <- x[["designObj"]][["h0"]]

    upperLine <- confSeqMatrix[, 2]
    lowerLine <- confSeqMatrix[, 1]

    upperLineFinite <- upperLine[is.finite(upperLine)]
    lowerLineFinite <- lowerLine[is.finite(lowerLine)]

    maxBound <- max(upperLineFinite)
    minBound <- min(lowerLineFinite)

    if (maxBound > 0)
      upperLine[is.infinite(upperLine)] <- 2*maxBound
    else if (maxBound <= 0)
      upperLine[is.infinite(upperLine)] <- 1/2*maxBound

    if (minBound > 0)
      lowerLine[is.infinite(lowerLine)] <- 1/2*minBound
    else if (minBound <= 0)
      lowerLine[is.infinite(lowerLine)] <- 2*minBound

    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes();

    maxX <- max(n1Vec)

    if (is.null(logScale))
      logScale <- if (maxX > switchNLog) TRUE else FALSE

    logPlot <- if (isTRUE(logScale)) "x" else ""

    if (is.null(xlim))
      xlim <- c(0.9, maxX)

    if (is.null(ylim))
      ylim <- c(minBound, maxBound)

    plot(NULL, xlim=xlim, ylim=ylim,
         type="l", xlab = "", ylab = "", cex.lab = cex,
         cex.axis = cex, xaxt="n", yaxt="n", bty="n", log=logPlot)

    if (is.null(fillPlot))
      fillPlot <- if (maxX <= switchNFill) TRUE else FALSE

    if (fillPlot) {
      polygon(c(n1Vec, rev(n1Vec)),
              c(upperLine, rev(lowerLine)),
              col=col, border=border, lwd=lwd,
              density = NULL)
    } else {
      lines(n1Vec, upperLine, col=border, lwd=lwd)
      lines(n1Vec, lowerLine, col=border, lwd=lwd)
    }

    lines(c(1, maxX), c(h0, h0), lwd=lwd, lty=2, col=lineColour)

    axis(1)
    axis(2)

    if (is.null(ylab))
      ylab <- switch(x[["testName"]],
                     "Z-Test"="mu",
                     "T-Test"="mu",
                     "logrank"="log(hazard ratio)")

    mtext(ylab, side = 2, line = 4, las = 0, cex = cex, adj=0.5)
    mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)
  }

  # e-value plot
  if (!isTRUE(wantConfSeqPlot)) {
    alpha <- x[["designObj"]][["alpha"]]
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes();

    maxEValue <- max(eValueVec)
    minEValue <- min(eValueVec)

    if (is.infinite(maxEValue)) {
      warning("Overflow: E-values infinite, removed for plotting")

      finiteIndex <- which(is.finite(eValueVec))
      eValueVec <- eValueVec[finiteIndex]

      n1Vec <- n1Vec[finiteIndex]
      maxEValue <- max(eValueVec)
    }

    rangeEValue <- maxEValue-minEValue

    if (is.null(logScale)) {
      logScale <- FALSE

      logScale <- if (rangeEValue/(1/alpha-1) > 5) TRUE else FALSE

      if (abs(log(minEValue)) > abs(log(maxEValue)))
        logScale <- TRUE
    }

    logPlot <- if (isTRUE(logScale)) "y" else ""

    maxY <- ceiling(max(maxEValue, 1/alpha))
    minY <- if (isTRUE(logScale)) minEValue else 0

    if (is.null(ylim))
      ylim <- c(minY, maxY)

    lastIndex <- length(n1Vec)
    maxX <- n1Vec[lastIndex]

    plot(n1Vec, eValueVec, type="l", lwd=lwd, xlab="", ylim=ylim,
         ylab="", col=border, xaxt="n", yaxt="n", bty="n", log=logPlot)

    threshLine <- c(1/alpha, 1/alpha)
    unitLine <- c(1, 1)

    lines(c(1, maxX), threshLine, lwd=lwd, lty=2, col="grey40")

    if (maxY/threshLine[1] < 10)
      lines(c(1, maxX), unitLine, lwd=lwd, lty=3, col="grey60")

    axis(side = 1)
    axis(side = 2)


    if (is.null(ylab))
      ylab <- "e-value"

    mtext(ylab, side = 2, line = 4, las = 0, cex = cex, adj=0.5)
    mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)

  }
}

