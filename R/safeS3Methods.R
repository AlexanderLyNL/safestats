# S3 helpers ---------
#
#' Gets the Label of the Test
#'
#' Helper function that outputs the name of the analysis.
#'
#' @inheritParams designSaviZ
#' @param testType A character string. For the t-tests: "oneSample", "paired", "twoSample".
#' @param testName The name of the analysis that is performed such as "Z-Test",
#' and "Test of Two Proportions".
#'
#' @return Returns a character string with the name of the analysis.
getNameTestType <- function(testType, testName, relevanceTest=FALSE) {

  if (isFALSE(relevanceTest)) {
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
  } else if (relevanceTest) {
    testTypeChar <- switch(testType,
                           "oneSample"="One Sample",
                           "paired"="Paired Sample",
                           "twoSample"="Two Sample",
                           "gLogrank"="Gaussian",
                           "eLogrank"="Exact",
                           "logrank"="",
                           "2x2" = "Test of ")
    analysisName <- paste("Savi", testTypeChar, "e-Relevance", testName)
    return(analysisName)
  }

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
    "relevanceTestSim"=NULL, "varEqual"=NULL,
    "samplePaths"=NULL, "breakVector"=NULL, "designScenario"=NULL,
    "call"=NULL, "timeStamp"=Sys.time(), "note"=NULL)

  testSpecificList <- list()

  if (testName=="Z-Test") {
    testSpecificList <- list("sigma"=NULL, "kappa"=NULL, "ratio"=NULL, "testName"=testName)
  } else if (testName=="T-Test") {
    testSpecificList <- list("ratio"=NULL, "testName"=testName)
  } else if (testName=="Logrank") {
    testSpecificList <- list("exact"=NULL)
  } else if (testName=="Two Proportions") {
    testSpecificList <- list(
      "betaParameter"=NULL,
      "gaussParameter"=NULL,
      "relevanceTest"=FALSE,
      "runningIntersection"=FALSE,
      "testName"=testName)
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
  } else if (testName=="Two Proportions") {
    testSpecificList <- list("betaPrior"=NULL, "testName"=testName)
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
#' @param runningIntersection logical, if \code{TRUE} then report the minimum
#' of the upper confidence sequence across time, and the maximum of the lower
#' confidence sequence across time
#' @param ... further arguments to be passed to or from methods.
#'
#' @return No returned value, called for side effects.
#' @export
#'
#' @examples
#' saviTTest(rnorm(19))
print.saviTest <- function(x, digits = getOption("digits"), prefix = "\t",
                           runningIntersection=NULL, ...) {
  designObj <- x[["designObj"]]

  if (is.null(designObj)) {
    print.default(x)
    return()
  }

  if (is.null(runningIntersection)) {
    runningIntersection <- designObj[["runningIntersection"]]

    if (is.null(runningIntersection))
      runningIntersection <- FALSE
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
    if (runningIntersection && !is.null(x[["confSeqMatrix"]])) {
      confSeq <- c(max(x[["confSeqMatrix"]][, 1]),
                   min(x[["confSeqMatrix"]][, 2]))
    }

    cat(format(100*(ciValue)), " percent confidence sequence:\n",
        " ", paste(format(confSeq, digits = digits),
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

  # A 2x2 grow design holds no parameter: its esMin is shown instead.
  if (!is.null(parameter))
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
      # The 2x2 nPlan is a list; per-block sizes print as their mean.
      if (is.list(nPlan)) {
        isVector <- lengths(nPlan) > 1L
        nPlan[isVector] <- lapply(nPlan[isVector], mean)
        names(nPlan)[isVector] <- paste("mean", names(nPlan)[isVector])
        nPlan <- unlist(nPlan)
      }
      out <- paste(names(nPlan), "=", sapply(nPlan, format, digits = max(1L, digits - 2L)))
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

  # A 2x2 grow design holds no parameter: its esMin is shown instead.
  if (is.null(designObj[["parameter"]]) && !identical(designObj[["testType"]], "2x2")) {
    print.default(x)
    return()
  }

  testName <- designObj[["testName"]]
  testType <- designObj[["testType"]]

  note <- designObj[["note"]]

  tempName <- getNameTestType("testType"=testType, "testName"=testName)

  # extraChar <- if (relevanceTest) " + e-Relevance Test" else ""

  analysisName <- paste0(tempName, " Design")

  cat("\n")
  cat(strwrap(analysisName, prefix = prefix), sep = "\n")
  cat("\n")

  displayList <- list()

  designObj[["decision rule"]] <- 1/designObj[["alpha"]]

  for (item in c("nPlan", "nEvents", "nMean", "esMin", "alternative",
                 "alternativeRestriction", "power", "beta",
                 "eType", "parameter", "alpha",
                 "decision rule",
                 "logImpliedTarget")) {
    itemValue <- designObj[[item]]
    itemValueString <- format(itemValue, digits=digits)

    if (!is.null(itemValue)) {
      if (item %in% c("nPlan", "nMean")) {
        itemNeem <- paste0(item, "TwoSe")

        itemTwoSe <- designObj[[itemNeem]]

        # The 2x2 nPlan is a list; per-block sizes print as their mean.
        if (is.list(itemValue)) {
          isVector <- lengths(itemValue) > 1L
          itemValue[isVector] <- lapply(itemValue[isVector], mean)
          names(itemValue)[isVector] <- paste("mean", names(itemValue)[isVector])
          itemValue <- unlist(itemValue)
        }

        if (!is.null(itemTwoSe)) {
          # An NA standard error prints that element plainly: the 2x2 nPlan
          # holds na, nb, nBlocksPlan and only the last one is simulated.
          hasSe <- !is.na(itemTwoSe)
          itemValueString <- sapply(itemValue, format, digits=digits)
          itemValueString[hasSe] <- paste0(itemValueString[hasSe], "\U00B1",
                                           sapply(itemTwoSe[hasSe], format, digits=digits))
          tempNeem <- names(itemValue)
          tempNeem[hasSe] <- paste0(tempNeem[hasSe], "\U00B1", "2se")
          displayList[[paste(tempNeem, collapse=", ")]] <- paste(itemValueString, collapse=", ")
        } else {
          tempNeem <- names(itemValue)
          displayList[[paste(tempNeem, collapse=", ")]] <-
            paste(sapply(itemValue, format, digits=digits), collapse=", ")
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
        displayList[["Reject H_0: e-value >= 1/alpha"]] <- itemValueString
      } else if (item == "eType") {
        displayList[["e-variable type"]] <- itemValueString
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
      } else if (item == "alternativeRestriction") {
        displayList[["alternative restriction"]] <- itemValueString
      } else {
        displayList[[item]] <- itemValueString
      }
    }
  }

  cat(paste(format(names(displayList), width = 20L, justify = "right"),
            format(displayList, digits = digits), sep = " = "), sep = "\n")

  # Here done with standard
  relevanceTest <- if (!is.null(x[["relevanceTestSim"]])) TRUE else FALSE

  if (relevanceTest){
    tempName <- getNameTestType("testType"=testType, "testName"=testName, "relevanceTest"=relevanceTest)
    analysisName <- paste0(tempName, " Design")

    cat("\n")
    cat(strwrap(analysisName, prefix = prefix), sep = "\n")
    cat("\n")

    displayList <- list()

    relevanceObj <- designObj[["relevanceTestSim"]]

    relevanceObj[["Reject H_R"]] <- relevanceObj[["alpha"]]

    for (item in c("parameter", "alpha",
                   "Reject H_R")) {
      itemValue <- relevanceObj[[item]]
      itemValueString <- format(itemValue, digits=digits)

      if (!is.null(itemValue)) {
        if (item=="parameter") {
          displayList[[paste("parameter:")]] <- itemValueString
        } else if (item=="alpha") {
          displayList[["alphaR"]] <- itemValueString
        } else if (item=="Reject H_R") {
          displayList[["Reject H_R: e-relevance <= 1/alphaR"]] <- itemValueString
        }
      }
    }

    cat(paste(format(names(displayList), width = 20L, justify = "right"),
              format(displayList, digits = digits), sep = " = "), sep = "\n")
  }

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
#' @param histFewestAtTop logical, if TRUE, then fewest counts
#' (of against the null or against minimal efficacy) at the top
#' of the histogram
#' @param density the density of shading lines, in lines per inch.
#' The default value of NULL means that no shading lines are drawn.
#' Non-positive values of density also inhibit the drawing of shading lines.
#' @param overColour colour of an e-value sequence that cross the
#' upper boundary of 1/alpha
#' @param overColourBorder colour for the border of a rectangle in the histogram
#' corresponding to e-values that cross the upper boundary of 1/alpha
#' @param continueColour colour for lines of e-values that did not (yet) cross
#' a threshold
#' @param continueColourBorder colour for the border of a rectangle in the histogram
#' corresponding to e-values that did not (yet) cross a threshold
#' @param underColour colour of an e-value sequence that cross the
#' lower boundary of alphaRelevance
#' @param underColourBorder colour for the border of a rectangle in the histogram
#' corresponding to e-values that cross the lower boundary of alphaRelevance
#' @param wantNotStoppedHist logical, if \code{TRUE} then show histogram
#' of not stopped paths
#' @param wantNotStoppedSamplePaths logical, if \code{TRUE} then show stopped
#' sample paths of e-values
#' @param wantLegend logical, if \code{TRUE} then show information of the
#' percentage of stopped sample paths
#' @param legendAdjRelevance vector of 3 numerics that allow for the adjustment of
#' the percentages shown at the bottom of the plot, when relevance testing is on
#' @param legendAdj vector of 2 numerics that allow for the adjustment of
#' the percentages shown at the bottom of the plot
#' @param legendCexFactor numeric, additional factor to control the text size of
#' the additional information such as the percentages at the bottom of the plot
#' @param wantTitle logical, if \code{TRUE} then add title
#' @param yLabPAdj numeric to adjust the position of the y-axis label
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Nothing it only plots
#' @export
#'
#'
plot.saviDesign <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                            xlim=NULL, ylim=NULL, maxNBins=40,
                            numSamplePaths=100, wantStepLines=FALSE,
                            wantQuantiles=NULL,
                            breaks=NULL, lwd=2, pch=15,
                            colQuant="#AA0000",
                            overColour="#8FC6E3",
                            overColourBorder="#1F78B4E6",
                            underColour="#E7A72199",
                            underColourBorder="#BD8E17FF",
                            continueColour="#8CA252E6",
                            continueColourBorder="#637939E6",
                            histInnerColour=col,
                            col=overColour,
                            border=overColourBorder,
                            cex=1.3, yLabPAdj=-1,
                            wantNotStoppedHist=FALSE,
                            wantNotStoppedSamplePaths=TRUE,
                            wantLegend=TRUE,
                            legendAdjRelevance=c(-0.1, 0.5, 1),
                            legendAdj=c(0.2, 0.8),
                            legendCexFactor=0.85,
                            histFewestAtTop=FALSE,
                            density=NULL, wantTitle=TRUE, ...) {

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

    stopHist <- graphics::hist(
      stoppingTimes, freq=FALSE,
      breaks=breaks,
      col=histInnerColour,
      border=border, lwd=lwd, xaxt="n", yaxt="n",
      xlab="", ylab="", main=main)

    if (is.null(xlab))
      xlab <- "Stopping times"

    graphics::axis(side = 2)
    graphics::axis(side = 1)
    graphics::axis(side = 1, at=c(0, 2*nPlanBatch))

    if (is.null(ylab))
      ylab <- "Density"

    graphics::mtext(ylab, side = 2, line = 4, las = 0, cex = cex, adj=0.5)
    graphics::mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)

    if (maxStoppingTimes-minStoppingTimes > maxNBins) {
      firstPassageTimes <- stoppingTimes[!breakVector]
      fptHist <- graphics::hist(
        firstPassageTimes, breaks=stopHist[["breaks"]],
        plot=FALSE)

      lastIndex <- length(stopHist[["counts"]])

      graphics::rect(xleft=stopHist[["breaks"]][lastIndex], ybottom=fptHist[["density"]][lastIndex],
           xright=stopHist[["breaks"]][lastIndex+1],
           ytop=stopHist[["density"]][lastIndex],
           col=col, border=border)
    } else {
      lastIndex <- length(stopHist[["counts"]])
      rejectedAtNBatch <- stopHist[["counts"]][lastIndex]-sum(breakVector)

      totalCount <- sum(stopHist[["counts"]])

      graphics::rect(xleft=stopHist[["breaks"]][lastIndex], ybottom=rejectedAtNBatch/totalCount,
           xright=stopHist[["breaks"]][lastIndex+1],
           ytop=stopHist[["density"]][lastIndex],
           col=col,
           border=border)
    }

    graphics::lines(c(nPlan, nPlan), c(0, max(stopHist[["counts"]])), lwd=lwd, lty=2)
  }

  alphaRelevance <- NULL
  relevanceTestSim <- x[["relevanceTestSim"]]
  histRelevance <- NULL
  relevanceTest <- x[["relevanceTest"]]

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

    fptAlt <- stoppingTimes[indexStopAlt]
    fptNot <- stoppingTimes[indexStopNot]

    if (relevanceTest) {
      fptFut <- as.integer(relevanceTestSim[["stoppingTimes"]][indexStopFut])
      alphaRelevance <- relevanceTestSim[["alpha"]]

      fptAll[indexStopFut] <- fptFut
    }

    fptAllFinite <- fptAll

    # fptStopped <- fptAll
    # fptStopped[indexStopNot] <- Inf
    #
    # if (isFALSE(wantNotStoppedHist))
    #   fptAll <- fptStopped

    # Compute hist  ----
    #
    breaksMin <- min(fptAll)

    nMax <- min(nPlanBatch, max(fptAllFinite))

    if (is.null(breaks)) {
      breaks <- if (nMax-breaksMin+1 > maxNBins) {
        seq(breaksMin-1, nMax, length=maxNBins)
      } else {
        (breaksMin-1):nMax
      }
    }

    histAlt <- if (nAlt > 0) graphics::hist(fptAlt, plot=FALSE, breaks=breaks) else
      list(counts=rep(0, length=length(breaks)-1))

    histFut <- if (isTRUE(relevanceTest) && nFut > 0) graphics::hist(fptFut, plot=FALSE, breaks=breaks) else
      list(counts=rep(0, length=length(breaks)-1))

    histNot <- if (isTRUE(wantNotStoppedHist) && nNotStopped > 0) graphics::hist(fptNot, plot=FALSE, breaks=breaks) else
      list(counts=rep(0, length=length(breaks)-1))

    nHistTotal <- sum(histAlt[["counts"]], histFut[["counts"]], histNot[["counts"]])

    breaks <- histAlt[["breaks"]]

    yMaxHist <- max(histAlt[["counts"]], histFut[["counts"]])/nHistTotal

    yMin <- -1*log(3/(2*alpha))

    if (isTRUE(relevanceTest))
      yMin <- min(c(yMin, -1*log(3/(2*alphaRelevance))))

    if (is.null(ylim))
      ylim <- c(yMin, 2.75*log(1/alpha))

    someConstant <- 0.8*(ylim[2]+log(alpha))/yMaxHist
    textHeightQuant <- (ylim[2]+log(alpha))+log(1/alpha)

    if (is.null(xlim))
      xlim <- c(0, nPlanBatch)

    # Draw canvas ----
    #
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()

    plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
         cex.lab = cex, cex.axis = cex, las = 1, main=main,
         xaxt = "n", yaxt = "n", bty = "n", type = "p", pch = pch,
         bg = "grey", ...)

    graphics::lines(x=c(0, 1.5*nPlanBatch), y=c(0, 0), col="darkgrey", lwd=lwd, lty=2)

    graphics::lines(x=c(0, 1.5*nPlanBatch), y=log(c(1/alpha, 1/alpha)))

    if (relevanceTest)
      graphics::lines(x=c(0, 1.5*nPlanBatch), y=log(c(alphaRelevance, alphaRelevance)))

    if (relevanceTest) {
      yLabs <- c(1e-24, alphaRelevance, "1", 1/alpha)
      criticalP <- log(c(1e-24, alphaRelevance, 1, 1/alpha))
    } else {
      yLabs <- c(1e-24, alpha, "1", 1/alpha)
      criticalP <- log(c(1e-24, alpha, 1, 1/alpha))
    }

    graphics::axis(side = 2, at = criticalP, tick = TRUE, las = 2,
                   cex.axis = cex, labels = yLabs)
    graphics::axis(side = 1)
    graphics::axis(side = 1, at=c(0, 2*xlim[2]))

    if (is.null(ylab))
      ylab <- "Evidence"

    graphics::mtext(ylab, side = 2, line = 2.5, las = 0, cex = cex, adj=0.25, padj=yLabPAdj)

    if (is.null(xlab))
      xlab <- "Sample size"

    graphics::mtext(xlab, side = 1, line = 2.5, las = 1, cex = cex)

    # Draw hist -----
    #
    histOrder <- if (nAlt >= nFut) c(2, 1, 3) else 1:3

    if (histFewestAtTop)
      histOrder <- if (nAlt >= nFut) 1:3 else c(2, 1, 3)

    # histOrder <- order(c(nAlt, nFut, nNotStopped), decreasing=decreasingHist)
    lineColours <- c(col, underColour, continueColour)
    borderColours <- c(border, underColourBorder, continueColourBorder)

    countMat <- matrix(c(histAlt[["counts"]],
                         histFut[["counts"]],
                         histNot[["counts"]]), nrow=3, byrow=TRUE)

    countMat <- countMat[histOrder, ]

    cumsumCountMat <- apply(countMat, 2, cumsum)

    for (i in 1:3) {
      someRow <- countMat[i, ]
      tempIndex <- which(countMat[i, ]!=0)

      if (i==1) {
        tempHeightLower <- rep(0, length(breaks)-1)
      } else {
        tempHeightLower <- cumsumCountMat[i-1, ]/nHistTotal
      }

      tempHeightUpper <- cumsumCountMat[i, ]/nHistTotal

      tempColour <- lineColours[histOrder[i]]
      tempColourBorder <- borderColours[histOrder[i]]

      for (j in seq_along(tempIndex)) {
        leftCornerIndex <- tempIndex[j]

        leftCorner <- breaks[leftCornerIndex]+0.5
        rightCorner <- breaks[leftCornerIndex+1]+0.5

        graphics::rect(leftCorner, someConstant*tempHeightLower[leftCornerIndex] + log(1/alpha),
             rightCorner, someConstant*tempHeightUpper[leftCornerIndex] + log(1/alpha),
             col = tempColour, border = tempColourBorder, lwd=2,
             angle = 45, density = density, lty = NULL)
      }
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

    if (relevanceTest)
      stoppedPaths[indexStopFut, ] <- x[["relevanceTestSim"]][["samplePaths"]][indexStopFut, ]

    nLinesNot <- if (isFALSE(wantNotStoppedSamplePaths)) 0 else
      min(ceil(nNotStopped/nAll*numSamplePaths), nNotStopped)

    if (nAlt >= nFut) {
      nLinesFut <- min(ceil(nFut/nAll*numSamplePaths), nFut)
      nLinesAlt <- min(numSamplePaths-nLinesFut-nLinesNot, nAlt)
    } else {
      nLinesAlt <- min(ceil(nAlt/nAll*numSamplePaths), nAlt)
      nLinesFut <- min(numSamplePaths-nLinesAlt-nLinesNot, nFut)
    }

    linesOrder <- order(c(nAlt, nFut, nNotStopped), decreasing=TRUE)

    nLinesToPlot <- c(nLinesAlt, nLinesFut, nLinesNot)[linesOrder]

    lastPoints <- if (!is.null(alphaRelevance)) c(1/alpha, alphaRelevance) else
      c(1/alpha, NULL)

    lastPoints <- c(lastPoints, NULL)[linesOrder]

    lineIndexes <- list(indexStopAlt, indexStopFut, indexStopNot)[linesOrder]

    for (i in 1:3) {
      someLineColour <- lineColours[linesOrder[i]]
      someBorderColour <- borderColours[linesOrder[i]]
      someNLinesToPlot <- nLinesToPlot[i]
      someLastPoint <- lastPoints[i]
      someIndexSet <- lineIndexes[[i]]

      someLineColour <- grDevices::adjustcolor(
        someLineColour,
        alpha.f=max(1-someNLinesToPlot/sum(nLinesToPlot), 0.25))

      if (someNLinesToPlot==0)
        next()

      for (j in 1:someNLinesToPlot) {
        stoppedTime <- fptAllFinite[someIndexSet[j]]

        evidenceLine <- stoppedPaths[someIndexSet[j], 1:stoppedTime]

        if (!is.null(someLastPoint) && !is.na(someLastPoint))
          evidenceLine[stoppedTime] <- someLastPoint

        if (isTRUE(wantStepLines)) {
          xLine <- c(0, rep(1:stoppedTime, each=2))
          yLine <- c(0, 0, rep(log(evidenceLine), each=2))
          yLine <- yLine[-length(yLine)]
        } else {
          xLine <- 0:stoppedTime
          yLine <- c(0, log(evidenceLine))
        }

        graphics::lines(xLine, yLine, col=someLineColour, lwd=lwd, lty=1)

        if (!is.null(someLastPoint) && !is.na(someLastPoint))
          graphics::points(stoppedTime, log(evidenceLine[stoppedTime]),
                 col=someBorderColour, pch=pch, lwd=lwd)
      }
    }

    if (is.null(wantQuantiles) && !is.null(x[["power"]]))
      wantQuantiles <- x[["power"]]

    if (!is.null(wantQuantiles) && !isFALSE(wantQuantiles)) {
      graphics::mtext("Quantiles [%]", side=2, col=colQuant, cex=cex, adj=0.5, at=textHeightQuant)

      badQuantileIndeces <- which(wantQuantiles > nAlt/nAll)

      if (length(badQuantileIndeces) >= 1) {
        warning(paste0("The requested quantile(s) ", wantQuantiles[badQuantileIndeces], " cannot be displayed. Max power in this simulation run: ", nAlt/nAll))
        wantQuantiles <- wantQuantiles[-badQuantileIndeces]
      }


      quants <- round(stats::quantile(stoppingTimes, wantQuantiles), 2)

      quantileNames <- round(wantQuantiles*100, 2)

      for (i in seq_along(quants)) {

        if (relevanceTest) {
          quantFut <- round(sum(fptFut <= quants[i])/nAll*100, 2)

          graphics::text(quantileNames[i], x=quants[i], y=textHeightQuant,
               col=overColourBorder, cex=legendCexFactor*cex, pos=2)
          graphics::text("+", x=quants[i],
               y=textHeightQuant, col=colQuant, cex=cex)
          graphics::text(quantFut, x=quants[i], y=textHeightQuant, col=underColourBorder, cex=legendCexFactor*cex, pos=4)
        } else {
          graphics::text(quantileNames[i], x=quants[i],
               y=textHeightQuant, col=overColourBorder, cex=cex)
        }

        # TODO(Alexander): perhaps change the lower position as minimum of -log(1/alpha) and log(beta)
        #
        graphics::segments(x0=quants[i], y0=-0.9*log(1/alpha), y1=0.95*textHeightQuant, col=colQuant)
        graphics::text(quants[i], x=quants[i], y=-log(1/alpha), col=colQuant, cex=cex)
      }
    }

    if (wantTitle)
      graphics::mtext(paste0("True ", names(x[["esMin"]]), " = ",
                             round(x[["esTrue"]], 3)),
            col=colQuant, side=3, line = 1, las = 1,
            cex = cex*legendCexFactor, adj=0)

    if (wantLegend) {

      if (relevanceTest) {
        graphics::mtext(paste0("Alternative: ",
                               round(nAlt/nAll*100, 1), "%"),
              col=border, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdjRelevance[1])

        graphics::mtext(paste0("Inferior effect: ",
                               round(nFut/nAll*100, 1), "%"),
              col=underColourBorder, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdjRelevance[2])

        graphics::mtext(paste0("No decision: ",
                               round(nNotStopped/nAll*100, 1), "%"),
              col=continueColourBorder, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdjRelevance[3])
      } else {
        graphics::mtext(paste0("Alternative: ",
                               round(nAlt/nAll*100, 1), "%"),
              col=border, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdj[1])
        graphics::mtext(paste0("No rejection: ",
                               round(nNotStopped/nAll*100, 1), "%"),
              col=continueColourBorder, side=1, line = 4, las = 1,
              cex = cex*legendCexFactor, adj=legendAdj[2])
      }



      # legend("topleft",
      #        legend=c(paste("Alternative: ", round(nAlt/nAll*100, 1), "%"),
      #                 paste("Futility: ", round(nFut/nAll*100, 1), "%"),
      #                 paste("Continued: ", round(nNotStopped/nAll*100, 1)), "%"),
      #        col=c(col, underColour, continueColour),
      #        lty=1, cex=cex, lwd=lwd, box.lty=box.lty, xpd=TRUE)
    }
  }
}





#' Plots the saviTest object for sequential analyses
#'
#' @inheritParams plot.saviDesign
#'
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
#' @param runningIntersection logical, if \code{TRUE} then plot the running
#' minimum of the upper confidence sequence across time, and the running
#' maximum of the lower confidence sequence across time
#' @param wantRelevance logical, if \code{FALSE}, then don't show the
#' e-values for relevanceTest. Default \code{wantRelevance==NULL}, if
#' \code{designObj[["relevanceTest"]]==TRUE} then relevance e-values tests are shown
#' automatically.
#' @param xaxt default NULL. If "n" then suppresses plotting of the x-axis.
#' @param yaxt default NULL. If "n" then suppresses plotting of the y-axis.
#'
#' @return Returns nothing just plots
#' @export
#'
plot.saviTest <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                          xlim=NULL, ylim=NULL, lwd=3, cex=1.3,
                          fillPlot=NULL, switchNFill=1e4,
                          logScale=NULL, switchNLog=30,
                          h0Colour="darkgrey",
                          overColour="#8FC6E3",
                          overColourBorder="#1F78B4E6",
                          underColour="#E7A72199",
                          underColourBorder="#BD8E17FF",
                          continueColour="#8CA252E6",
                          continueColourBorder="#637939E6",
                          col=overColour,
                          border=overColourBorder,
                          wantConfSeqPlot=FALSE, add=FALSE,
                          density=NULL, angle=45,
                          xaxt=NULL, yaxt=NULL,
                          fillOddEven=FALSE, runningIntersection=NULL,
                          wantRelevance=NULL,
                          ...) {
  eValueVec <- x[["eValueVec"]]
  eRelevanceVec <- x[["eRelevanceVec"]]

  n1Vec <- x[["n1Vec"]]

  designObj <- x[["designObj"]]

  if (is.null(runningIntersection)) {
    runningIntersection <- designObj[["runningIntersection"]]

    if (is.null(runningIntersection))
      runningIntersection <- FALSE
  }

  relevanceTest <- FALSE

  if ((is.null(wantRelevance) || isTRUE(wantRelevance)) && designObj[["relevanceTest"]] && !is.null(eRelevanceVec)) {
    relevanceTest <- TRUE
    alphaRelevance <- designObj[["relevanceTestSim"]][["alpha"]]
  }

  if (is.null(n1Vec) || is.null(eValueVec)) {
    warning("Can't plot. No sequential analysis.")
    return()
  }

  if (is.null(xlab)) {
    xlab <- switch(x[["testName"]],
                   "Z-Test"="n1",
                   "T-Test"="n1",
                   "logrank"="Number of events",
                   "Two Proportions"="Number of blocks")
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

      if (is.null(xaxt) || xaxt!="n") graphics::axis(1)
      if (is.null(yaxt) || yaxt!="n") graphics::axis(2)

      if (is.null(ylab))
        ylab <- switch(x[["testName"]],
                       "Z-Test"="mu",
                       "T-Test"="mu",
                       "logrank"="log(hazard ratio)",
                       "Two Proportions"=if (x[["designObj"]][["eType"]] == "eBeta") "propDiff" else "logOdds")

      graphics::mtext(ylab, side = 2, line = 4,
                      las = 0, cex = cex, adj=0.5)
      graphics::mtext(xlab, side = 1, line = 2.5,
                      las = 1, cex = cex)
    }

    if (is.null(fillPlot))
      fillPlot <- if (maxX <= switchNFill) TRUE else FALSE

    if (runningIntersection) {
      upperLine <- makeRunningIntersection(upperLine)
      lowerLine <- makeRunningIntersection(lowerLine,
                                           upper=FALSE)
    }

    if (fillPlot) {
      graphics::polygon(c(n1Vec, rev(n1Vec)),
              c(upperLine, rev(lowerLine)),
              col=col, border=border, lwd=lwd,
              density = density, angle=angle,
              fillOddEven=fillOddEven)
    } else {
      graphics::lines(n1Vec, upperLine, col=border, lwd=lwd)
      graphics::lines(n1Vec, lowerLine, col=border, lwd=lwd)
    }

    graphics::lines(c(1, maxX), c(h0, h0), lwd=lwd, lty=2, col=h0Colour)
  }

  # e-value plot
  if (!isTRUE(wantConfSeqPlot)) {
    alpha <- designObj[["alpha"]]
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes();

    nInfinite <- sum(is.infinite(eValueVec))

    if (relevanceTest) {
      nInfinite <- nInfinite+sum(is.infinite(eRelevanceVec))
    }

    if (nInfinite >= 1) {
      warning("Overflow: E-values infinite, removed for plotting")

      finiteIndex <- which(is.finite(eValueVec))

      if (relevanceTest) {
        finiteIndex <- intersect(finiteIndex, which(is.finite(eRelevanceVec)))
        eRelevanceVec <- eRelevanceVec[finiteIndex]
      }

      eValueVec <- eValueVec[finiteIndex]
      n1Vec <- n1Vec[finiteIndex]
    }

    maxEValue <- max(eValueVec, na.rm=TRUE)
    minEValue <- min(eValueVec, na.rm=TRUE)

    if (relevanceTest) {
      maxEValue <- max(eRelevanceVec, maxEValue, na.rm=TRUE)
      minEValue <- min(eRelevanceVec, minEValue, na.rm=TRUE)
    }

    if (minEValue <= .Machine[["double.xmin"]])
      minEValue <- 1e-300

    if (maxEValue >= .Machine[["double.xmax"]])
      maxEValue <- 1e300

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

      if (relevanceTest)
        minY <- min(minY, alphaRelevance)

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

      graphics::lines(c(1, maxX), threshLine, lwd=lwd, lty=2, col="grey40")

      if (maxY/threshLine[1] < 10)
        graphics::lines(c(1, maxX), unitLine, lwd=lwd, lty=3, col="grey60")

      if (relevanceTest)
        graphics::lines(c(1, maxX), c(alphaRelevance, alphaRelevance),
              lwd=lwd, lty=2, col="grey40")

      if (is.null(xaxt) || xaxt!="n") graphics::axis(1)
      if (is.null(yaxt) || yaxt!="n") graphics::axis(2)

      if (is.null(ylab))
        ylab <- "e-value"

      graphics::mtext(ylab, side = 2,
                      line = 4, las = 0, cex = cex, adj=0.5)
      graphics::mtext(xlab, side = 1,
                      line = 2.5, las = 1, cex = cex)
    }

    graphics::lines(n1Vec, eValueVec, lwd=lwd, col=overColour)

    if (relevanceTest)
      graphics::lines(n1Vec, eRelevanceVec, lwd=lwd, col=underColour)
  }
}

# ggplot2 versions -----

#' ggplot2 version of plot.saviTest
#'
#' Plots either the sequence of e-values or the confidence sequence of a
#' saviTest object using ggplot2. Returns a ggplot object, so the plot can be
#' further modified with the usual `+` syntax.
#'
#' @inheritParams plot.saviDesign,plot.saviTest
#'
#' @param x a saviTest object.
#' @param add either FALSE (default) or an existing ggplot object to which the
#'   layers are added (the ggplot analogue of `add=TRUE` in base graphics). If
#'   `TRUE` is passed, the last plot `ggplot2::last_plot()` is used.
#' @param lwd line width in base graphics units; converted to ggplot2 linewidth.
#' @param cex scaling factor for the text size.
#' @param density,angle,fillOddEven ignored; hatching is not supported by
#'   ggplot2 (see the ggpattern package). Kept for compatibility.
#' @param xaxt,yaxt use "n" to suppress the x or y axis.
#' @param xStart where the x-axis starts when `xlim` is not given. Defaults to
#'   1 on a log scale and 0 otherwise.
#'
#' @return a ggplot object.
ggplotSaviTest <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                           xlim=NULL, ylim=NULL, xStart=NULL, lwd=3, cex=1.3,
                           fillPlot=NULL, switchNFill=1e4,
                           logScale=NULL, switchNLog=30,
                           h0Colour="darkgrey",
                           overColour="#8FC6E3",
                           overColourBorder="#1F78B4E6",
                           underColour="#E7A72199",
                           underColourBorder="#BD8E17FF",
                           continueColour="#8CA252E6",
                           continueColourBorder="#637939E6",
                           col=overColour,
                           border=overColourBorder,
                           wantConfSeqPlot=FALSE, add=FALSE,
                           density=NULL, angle=45,
                           xaxt=NULL, yaxt=NULL,
                           fillOddEven=FALSE, runningIntersection=NULL,
                           wantRelevance=NULL,
                           ...) {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("Package 'ggplot2' is needed for this function.")

  eValueVec <- x[["eValueVec"]]
  eRelevanceVec <- x[["eRelevanceVec"]]

  n1Vec <- x[["n1Vec"]]

  designObj <- x[["designObj"]]

  # base lwd is in units of 1/96 inch, ggplot2 linewidth in mm
  lineWidth <- lwd/ggplot2::.pt

  if (isTRUE(add))
    add <- ggplot2::last_plot()

  isNewPlot <- !inherits(add, "ggplot")

  if (is.null(runningIntersection)) {
    runningIntersection <- designObj[["runningIntersection"]]

    if (is.null(runningIntersection))
      runningIntersection <- FALSE
  }

  relevanceTest <- FALSE

  if ((is.null(wantRelevance) || isTRUE(wantRelevance)) && isTRUE(designObj[["relevanceTest"]]) && !is.null(eRelevanceVec)) {
    relevanceTest <- TRUE
    alphaRelevance <- designObj[["relevanceTestSim"]][["alpha"]]
  }

  if (is.null(n1Vec) || is.null(eValueVec)) {
    warning("Can't plot. No sequential analysis.")
    return()
  }

  if (is.null(xlab)) {
    xlab <- switch(x[["testName"]],
                   "Z-Test"="n1",
                   "T-Test"="n1",
                   "logrank"="Number of events",
                   "Two Proportions"="Number of blocks")
  }

  if (isNewPlot) {
    # Mimic bty="n": no frame, and axis lines capped at the outer ticks so
    # they don't meet in a corner, as with graphics::axis()
    p <- ggplot2::ggplot() +
      ggplot2::theme_classic(base_size=11*cex) +
      ggplot2::theme(panel.border=ggplot2::element_blank(),
                     panel.background=ggplot2::element_blank(),
                     plot.background=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=10))) +
      ggplot2::guides(x=ggplot2::guide_axis(cap="both"),
                      y=ggplot2::guide_axis(cap="both"))

    if (!is.null(xaxt) && xaxt=="n")
      p <- p + ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                              axis.ticks.x=ggplot2::element_blank(),
                              axis.line.x=ggplot2::element_blank())

    if (!is.null(yaxt) && yaxt=="n")
      p <- p + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                              axis.ticks.y=ggplot2::element_blank(),
                              axis.line.y=ggplot2::element_blank())

    if (!is.null(main))
      p <- p + ggplot2::ggtitle(main)
  } else {
    p <- add
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

    if (runningIntersection) {
      upperLine <- makeRunningIntersection(upperLine)
      lowerLine <- makeRunningIntersection(lowerLine,
                                           upper=FALSE)
    }

    plotData <- data.frame("n1"=n1Vec, "lower"=lowerLine, "upper"=upperLine)

    if (fillPlot) {
      p <- p + ggplot2::geom_ribbon(data=plotData,
                                    ggplot2::aes(x=.data[["n1"]],
                                                 ymin=.data[["lower"]],
                                                 ymax=.data[["upper"]]),
                                    fill=col, colour=border,
                                    linewidth=lineWidth)
    } else {
      p <- p +
        ggplot2::geom_line(data=plotData,
                           ggplot2::aes(x=.data[["n1"]], y=.data[["upper"]]),
                           colour=border, linewidth=lineWidth) +
        ggplot2::geom_line(data=plotData,
                           ggplot2::aes(x=.data[["n1"]], y=.data[["lower"]]),
                           colour=border, linewidth=lineWidth)
    }

    p <- p + ggplot2::annotate("segment", x=1, xend=maxX, y=h0, yend=h0,
                               linewidth=lineWidth, linetype="dashed",
                               colour=h0Colour)

    if (isNewPlot) {
      if (is.null(logScale))
        logScale <- if (maxX > switchNLog) TRUE else FALSE

      xAxis <- makeAxisLimitsAndBreaks(xlim, defaultLim=c(NA, maxX),
                                       logScale=logScale, start=xStart)
      yAxis <- makeAxisLimitsAndBreaks(ylim, defaultLim=c(minBound, maxBound))

      if (isTRUE(logScale))
        p <- p + ggplot2::scale_x_log10(breaks=xAxis[["breaks"]])
      else
        p <- p + ggplot2::scale_x_continuous(breaks=xAxis[["breaks"]])

      p <- p + ggplot2::scale_y_continuous(breaks=yAxis[["breaks"]])

      if (is.null(ylab))
        ylab <- switch(x[["testName"]],
                       "Z-Test"="mu",
                       "T-Test"="mu",
                       "logrank"="log(hazard ratio)",
                       "Two Proportions"=if (x[["designObj"]][["eType"]] == "eBeta") "propDiff" else "logOdds")

      # coord_cartesian zooms instead of dropping data outside the limits,
      # which mimics base graphics clipping
      p <- p + ggplot2::coord_cartesian(xlim=xAxis[["limits"]],
                                        ylim=yAxis[["limits"]]) +
        ggplot2::labs(x=xlab, y=ylab)
    }

    return(p)
  }

  # e-value plot
  alpha <- designObj[["alpha"]]

  nInfinite <- sum(is.infinite(eValueVec))

  if (relevanceTest) {
    nInfinite <- nInfinite+sum(is.infinite(eRelevanceVec))
  }

  if (nInfinite >= 1) {
    warning("Overflow: E-values infinite, removed for plotting")

    finiteIndex <- which(is.finite(eValueVec))

    if (relevanceTest) {
      finiteIndex <- intersect(finiteIndex, which(is.finite(eRelevanceVec)))
      eRelevanceVec <- eRelevanceVec[finiteIndex]
    }

    eValueVec <- eValueVec[finiteIndex]
    n1Vec <- n1Vec[finiteIndex]
  }

  maxEValue <- max(eValueVec, na.rm=TRUE)
  minEValue <- min(eValueVec, na.rm=TRUE)

  if (relevanceTest) {
    maxEValue <- max(eRelevanceVec, maxEValue, na.rm=TRUE)
    minEValue <- min(eRelevanceVec, minEValue, na.rm=TRUE)
  }

  if (minEValue <= .Machine[["double.xmin"]])
    minEValue <- 1e-300

  if (maxEValue >= .Machine[["double.xmax"]])
    maxEValue <- 1e300

  if (isNewPlot) {
    rangeEValue <- maxEValue-minEValue

    if (is.null(logScale)) {
      logScale <- if (rangeEValue/(1/alpha-1) > 5) TRUE else FALSE

      if (abs(log(minEValue)) > abs(log(maxEValue)))
        logScale <- TRUE
    }

    maxY <- ceiling(max(maxEValue, 1/alpha))

    minY <- minEValue

    if (relevanceTest)
      minY <- min(minY, alphaRelevance)

    if (isFALSE(logScale))
      minY <- 0

    lastIndex <- length(n1Vec)
    maxX <- n1Vec[lastIndex]

    xAxis <- makeAxisLimitsAndBreaks(xlim, defaultLim=c(NA, maxX), start=xStart)
    yAxis <- makeAxisLimitsAndBreaks(ylim, defaultLim=c(minY, maxY),
                                     logScale=logScale)

    threshValue <- 1/alpha

    p <- p + ggplot2::annotate("segment", x=1, xend=maxX,
                               y=threshValue, yend=threshValue,
                               linewidth=lineWidth, linetype="dashed",
                               colour="grey40")

    if (maxY/threshValue < 10)
      p <- p + ggplot2::annotate("segment", x=1, xend=maxX, y=1, yend=1,
                                 linewidth=lineWidth, linetype="dotted",
                                 colour="grey60")

    if (relevanceTest)
      p <- p + ggplot2::annotate("segment", x=1, xend=maxX,
                                 y=alphaRelevance, yend=alphaRelevance,
                                 linewidth=lineWidth, linetype="dashed",
                                 colour="grey40")

    p <- p + ggplot2::scale_x_continuous(breaks=xAxis[["breaks"]])

    if (isTRUE(logScale))
      p <- p + ggplot2::scale_y_log10(breaks=yAxis[["breaks"]])
    else
      p <- p + ggplot2::scale_y_continuous(breaks=yAxis[["breaks"]])

    if (is.null(ylab))
      ylab <- "e-value"

    p <- p + ggplot2::coord_cartesian(xlim=xAxis[["limits"]],
                                      ylim=yAxis[["limits"]]) +
      ggplot2::labs(x=xlab, y=ylab)
  }

  plotData <- data.frame("n1"=n1Vec, "eValue"=eValueVec)

  p <- p + ggplot2::geom_line(data=plotData,
                              ggplot2::aes(x=.data[["n1"]], y=.data[["eValue"]]),
                              linewidth=lineWidth, colour=overColour)

  if (relevanceTest) {
    relevanceData <- data.frame("n1"=n1Vec, "eValue"=eRelevanceVec)

    p <- p + ggplot2::geom_line(data=relevanceData,
                                ggplot2::aes(x=.data[["n1"]], y=.data[["eValue"]]),
                                linewidth=lineWidth, colour=underColour)
  }

  return(p)
}

# Limits and breaks for one axis. A user-supplied lim is used as is. Otherwise
# defaultLim (first entry NA means "use start") is extended to the outermost
# breaks, so that the capped axis line covers the whole range, and the first
# break is put at start.
makeAxisLimitsAndBreaks <- function(lim, defaultLim, logScale=FALSE, start=NULL) {
  extend <- is.null(lim)

  if (extend) {
    if (is.null(start))
      start <- if (isTRUE(logScale)) 1 else 0

    if (isTRUE(logScale) && start <= 0)
      start <- 1

    lim <- defaultLim

    if (is.na(lim[1]))
      lim[1] <- start
    else
      start <- NULL
  } else {
    start <- NULL
  }

  if (isTRUE(logScale)) {
    logLim <- log10(lim)
    lo <- if (is.null(start)) floor(logLim[1]) else logLim[1]
    hi <- ceiling(logLim[2])
    step <- max(1, ceiling((hi-lo)/6))
    breaks <- 10^seq(lo, by=step, length.out=ceiling((hi-lo)/step)+1)

    # Rather than a whole extra decade, end at the first 1, 2 or 5 times a
    # power of ten above the data, e.g. 1, 10, 100, 200 for data up to 180
    if (extend) {
      niceValues <- c(1, 2, 5, 10)*10^floor(logLim[2])
      upper <- niceValues[niceValues >= lim[2]*(1-1e-8)][1]
      breaks <- unique(c(breaks[breaks < upper], upper))
    }
  } else {
    breaks <- pretty(lim)

    if (!is.null(start))
      breaks[1] <- start
  }

  if (extend)
    lim <- range(breaks)
  else
    breaks <- breaks[breaks >= lim[1] & breaks <= lim[2]]

  list("limits"=lim, "breaks"=breaks)
}

#' ggplot2 version of plot.saviDesign
#'
#' Plots the saviDesign object for designs with sample paths using ggplot2.
#' Takes the same arguments as plot.saviDesign, but returns a ggplot object,
#' so the plot can be further modified with the usual `+` syntax.
#'
#' Differences with the base graphics version:
#' * The percentages of stopped sample paths are shown as a ggplot2 legend
#'   below the plot. Hence, `legendAdj` and `legendAdjRelevance` are ignored.
#' * The "Quantiles [%]" label is drawn inside the plot at the top left.
#' * The "True ..." text is the subtitle, and `main` the title.
#' * `density` is ignored, as hatching is not supported by ggplot2 (see the
#'   ggpattern package), and so is `yLabPAdj`.
#'
#'
#' @inheritParams plot.saviDesign
#'
#' @return a ggplot object, or NULL if the design has no sample paths.
ggplotSaviDesign <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                             xlim=NULL, ylim=NULL, maxNBins=40,
                             numSamplePaths=100, wantStepLines=FALSE,
                             wantQuantiles=NULL,
                             breaks=NULL, lwd=2, pch=15,
                             colQuant="#AA0000",
                             overColour="#8FC6E3",
                             overColourBorder="#1F78B4E6",
                             underColour="#E7A72199",
                             underColourBorder="#BD8E17FF",
                             continueColour="#8CA252E6",
                             continueColourBorder="#637939E6",
                             histInnerColour=col,
                             col=overColour,
                             border=overColourBorder,
                             cex=1.3, yLabPAdj=-1,
                             wantNotStoppedHist=FALSE,
                             wantNotStoppedSamplePaths=TRUE,
                             wantLegend=TRUE,
                             legendAdjRelevance=c(-0.1, 0.5, 1),
                             legendAdj=c(0.2, 0.8),
                             legendCexFactor=0.85,
                             histFewestAtTop=FALSE,
                             density=NULL, wantTitle=TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("Package 'ggplot2' is needed for this function.")

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

  # base lwd is in units of 1/96 inch, ggplot2 linewidth in mm, and base text
  # sizes are relative to the 11pt base_size of the theme
  lineWidth <- lwd/ggplot2::.pt
  thinLineWidth <- 1/ggplot2::.pt
  textSize <- 11*cex/ggplot2::.pt

  # Mimic bty="n": no frame, and axis lines that end at the outer ticks
  p <- ggplot2::ggplot() +
    ggplot2::theme_classic(base_size=11*cex) +
    ggplot2::theme(panel.border=ggplot2::element_blank(),
                   panel.background=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=10)))

  if (numSamplePaths==0) {
    maxStoppingTimes <- max(stoppingTimes)
    minStoppingTimes <- min(stoppingTimes)

    if (is.null(breaks)) {
      breaks <- if (maxStoppingTimes-minStoppingTimes+1 > maxNBins)
        maxNBins
      else
        (minStoppingTimes-1):nPlanBatch
    }

    stopHist <- graphics::hist(stoppingTimes, breaks=breaks, plot=FALSE)
    histBreaks <- stopHist[["breaks"]]
    lastIndex <- length(stopHist[["counts"]])

    histData <- data.frame("xmin"=histBreaks[-(lastIndex+1)],
                           "xmax"=histBreaks[-1],
                           "density"=stopHist[["density"]])

    if (maxStoppingTimes-minStoppingTimes > maxNBins) {
      firstPassageTimes <- stoppingTimes[!breakVector]
      fptHist <- graphics::hist(
        firstPassageTimes, breaks=histBreaks,
        plot=FALSE)

      lastBottom <- fptHist[["density"]][lastIndex]
    } else {
      rejectedAtNBatch <- stopHist[["counts"]][lastIndex]-sum(breakVector)
      totalCount <- sum(stopHist[["counts"]])

      lastBottom <- rejectedAtNBatch/totalCount
    }

    if (is.null(xlab))
      xlab <- "Stopping times"

    if (is.null(ylab))
      ylab <- "Density"

    xAxis <- makeAxisLimitsAndBreaks(xlim, defaultLim=c(NA, max(histBreaks)),
                                     start=min(0, histBreaks[1]))
    yAxis <- makeAxisLimitsAndBreaks(ylim, defaultLim=c(0, max(stopHist[["density"]])))

    p <- p +
      ggplot2::geom_rect(data=histData,
                         ggplot2::aes(xmin=.data[["xmin"]], xmax=.data[["xmax"]],
                                      ymin=0, ymax=.data[["density"]]),
                         fill=histInnerColour, colour=border,
                         linewidth=lineWidth) +
      ggplot2::annotate("rect",
                        xmin=histBreaks[lastIndex], xmax=histBreaks[lastIndex+1],
                        ymin=lastBottom, ymax=stopHist[["density"]][lastIndex],
                        fill=col, colour=border, linewidth=thinLineWidth) +
      ggplot2::annotate("segment", x=nPlan, xend=nPlan, y=0, yend=Inf,
                        linewidth=lineWidth, linetype="dashed") +
      ggplot2::scale_x_continuous(breaks=xAxis[["breaks"]]) +
      ggplot2::scale_y_continuous(breaks=yAxis[["breaks"]]) +
      ggplot2::guides(x=ggplot2::guide_axis(cap="both"),
                      y=ggplot2::guide_axis(cap="both")) +
      ggplot2::coord_cartesian(xlim=xAxis[["limits"]], ylim=yAxis[["limits"]]) +
      ggplot2::labs(title=main, x=xlab, y=ylab)

    return(p)
  }

  if (numSamplePaths < 0)
    return()

  alphaRelevance <- NULL
  relevanceTestSim <- x[["relevanceTestSim"]]
  relevanceTest <- isTRUE(x[["relevanceTest"]])

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

  fptAlt <- stoppingTimes[indexStopAlt]
  fptNot <- stoppingTimes[indexStopNot]

  if (relevanceTest) {
    fptFut <- as.integer(relevanceTestSim[["stoppingTimes"]][indexStopFut])
    alphaRelevance <- relevanceTestSim[["alpha"]]

    fptAll[indexStopFut] <- fptFut
  }

  fptAllFinite <- fptAll

  # Compute hist  ----
  #
  breaksMin <- min(fptAll)

  nMax <- min(nPlanBatch, max(fptAllFinite))

  if (is.null(breaks)) {
    breaks <- if (nMax-breaksMin+1 > maxNBins) {
      seq(breaksMin-1, nMax, length=maxNBins)
    } else {
      (breaksMin-1):nMax
    }
  }

  emptyHist <- list(counts=rep(0, length=length(breaks)-1))

  histAlt <- if (nAlt > 0) graphics::hist(fptAlt, plot=FALSE, breaks=breaks) else
    emptyHist

  histFut <- if (isTRUE(relevanceTest) && nFut > 0) graphics::hist(fptFut, plot=FALSE, breaks=breaks) else
    emptyHist

  histNot <- if (isTRUE(wantNotStoppedHist) && nNotStopped > 0) graphics::hist(fptNot, plot=FALSE, breaks=breaks) else
    emptyHist

  nHistTotal <- sum(histAlt[["counts"]], histFut[["counts"]], histNot[["counts"]])

  if (!is.null(histAlt[["breaks"]]))
    breaks <- histAlt[["breaks"]]

  yMaxHist <- max(histAlt[["counts"]], histFut[["counts"]])/nHistTotal

  yMin <- -1*log(3/(2*alpha))

  if (isTRUE(relevanceTest))
    yMin <- min(c(yMin, -1*log(3/(2*alphaRelevance))))

  if (is.null(ylim))
    ylim <- c(yMin, 2.75*log(1/alpha))

  someConstant <- 0.8*(ylim[2]+log(alpha))/yMaxHist
  textHeightQuant <- (ylim[2]+log(alpha))+log(1/alpha)

  # Reference lines ----
  #
  p <- p +
    ggplot2::annotate("segment", x=0, xend=1.5*nPlanBatch, y=0, yend=0,
                      colour="darkgrey", linewidth=lineWidth, linetype="dashed") +
    ggplot2::annotate("segment", x=0, xend=1.5*nPlanBatch,
                      y=log(1/alpha), yend=log(1/alpha), linewidth=thinLineWidth)

  if (relevanceTest)
    p <- p + ggplot2::annotate("segment", x=0, xend=1.5*nPlanBatch,
                               y=log(alphaRelevance), yend=log(alphaRelevance),
                               linewidth=thinLineWidth)

  # Draw hist -----
  #
  histOrder <- if (nAlt >= nFut) c(2, 1, 3) else 1:3

  if (histFewestAtTop)
    histOrder <- if (nAlt >= nFut) 1:3 else c(2, 1, 3)

  lineColours <- c(col, underColour, continueColour)
  borderColours <- c(border, underColourBorder, continueColourBorder)

  countMat <- rbind(histAlt[["counts"]],
                    histFut[["counts"]],
                    histNot[["counts"]])[histOrder, , drop=FALSE]

  cumsumCountMat <- apply(countMat, 2, cumsum)

  heightLowerMat <- rbind(0, cumsumCountMat[-3, , drop=FALSE])/nHistTotal
  heightUpperMat <- cumsumCountMat/nHistTotal

  histRects <- vector("list", 3)

  for (i in 1:3) {
    tempIndex <- which(countMat[i, ]!=0)

    if (length(tempIndex)==0)
      next()

    histRects[[i]] <- data.frame(
      "xmin"=breaks[tempIndex]+0.5,
      "xmax"=breaks[tempIndex+1]+0.5,
      "ymin"=someConstant*heightLowerMat[i, tempIndex] + log(1/alpha),
      "ymax"=someConstant*heightUpperMat[i, tempIndex] + log(1/alpha),
      "fill"=lineColours[histOrder[i]],
      "border"=borderColours[histOrder[i]])
  }

  histRects <- do.call(rbind, histRects)

  if (!is.null(histRects))
    p <- p + ggplot2::geom_rect(data=histRects,
                                ggplot2::aes(xmin=.data[["xmin"]], xmax=.data[["xmax"]],
                                             ymin=.data[["ymin"]], ymax=.data[["ymax"]]),
                                fill=histRects[["fill"]], colour=histRects[["border"]],
                                linewidth=2/ggplot2::.pt)

  # Draw sample paths -----
  #
  stoppedPaths <- x[["samplePaths"]]

  if (relevanceTest)
    stoppedPaths[indexStopFut, ] <- relevanceTestSim[["samplePaths"]][indexStopFut, ]

  nLinesNot <- if (isFALSE(wantNotStoppedSamplePaths)) 0 else
    min(ceiling(nNotStopped/nAll*numSamplePaths), nNotStopped)

  if (nAlt >= nFut) {
    nLinesFut <- min(ceiling(nFut/nAll*numSamplePaths), nFut)
    nLinesAlt <- max(min(numSamplePaths-nLinesFut-nLinesNot, nAlt), 0)
  } else {
    nLinesAlt <- min(ceiling(nAlt/nAll*numSamplePaths), nAlt)
    nLinesFut <- max(min(numSamplePaths-nLinesAlt-nLinesNot, nFut), 0)
  }

  outcomes <- c("alt", "fut", "not")

  linesOrder <- order(c(nAlt, nFut, nNotStopped), decreasing=TRUE)

  nLinesToPlot <- c(nLinesAlt, nLinesFut, nLinesNot)[linesOrder]

  lastPoints <- c(1/alpha, if (relevanceTest) alphaRelevance else NA, NA)[linesOrder]

  lineIndexes <- list(indexStopAlt, indexStopFut, indexStopNot)[linesOrder]

  pathColours <- lineColours

  pathData <- list()
  endPointData <- list()

  for (i in 1:3) {
    someNLinesToPlot <- nLinesToPlot[i]

    if (someNLinesToPlot==0)
      next()

    someLastPoint <- lastPoints[i]
    someIndexSet <- lineIndexes[[i]]

    pathColours[linesOrder[i]] <- grDevices::adjustcolor(
      lineColours[linesOrder[i]],
      alpha.f=max(1-someNLinesToPlot/sum(nLinesToPlot), 0.25))

    for (j in 1:someNLinesToPlot) {
      stoppedTime <- fptAllFinite[someIndexSet[j]]

      evidenceLine <- stoppedPaths[someIndexSet[j], 1:stoppedTime]

      if (!is.na(someLastPoint))
        evidenceLine[stoppedTime] <- someLastPoint

      if (isTRUE(wantStepLines)) {
        xLine <- c(0, rep(1:stoppedTime, each=2))
        yLine <- c(0, 0, rep(log(evidenceLine), each=2))
        yLine <- yLine[-length(yLine)]
      } else {
        xLine <- 0:stoppedTime
        yLine <- c(0, log(evidenceLine))
      }

      pathId <- length(pathData)+1

      pathData[[pathId]] <- data.frame("n"=xLine, "logEvidence"=yLine,
                                       "path"=pathId,
                                       "outcome"=outcomes[linesOrder[i]])

      if (!is.na(someLastPoint))
        endPointData[[pathId]] <- data.frame("n"=stoppedTime,
                                             "logEvidence"=log(someLastPoint),
                                             "border"=borderColours[linesOrder[i]])
    }
  }

  pathData <- do.call(rbind, pathData)
  endPointData <- do.call(rbind, endPointData)

  if (!is.null(pathData)) {
    if (relevanceTest) {
      legendOutcomes <- outcomes
      legendLabels <- c(paste0("Alternative: ", round(nAlt/nAll*100, 1), "%"),
                        paste0("Inferior effect: ", round(nFut/nAll*100, 1), "%"),
                        paste0("No decision: ", round(nNotStopped/nAll*100, 1), "%"))
    } else {
      legendOutcomes <- c("alt", "not")
      legendLabels <- c(paste0("Alternative: ", round(nAlt/nAll*100, 1), "%"),
                        paste0("No rejection: ", round(nNotStopped/nAll*100, 1), "%"))
    }

    # Paths are drawn in the order of their group, that is, in the order above
    p <- p +
      ggplot2::geom_path(data=pathData,
                         ggplot2::aes(x=.data[["n"]], y=.data[["logEvidence"]],
                                      group=.data[["path"]],
                                      colour=.data[["outcome"]]),
                         linewidth=lineWidth) +
      ggplot2::scale_colour_manual(
        values=stats::setNames(pathColours, outcomes),
        limits=outcomes, breaks=legendOutcomes, labels=legendLabels,
        name=NULL,
        guide=ggplot2::guide_legend(
          override.aes=list(colour=borderColours[match(legendOutcomes, outcomes)])))
  }

  if (!is.null(endPointData))
    p <- p + ggplot2::geom_point(data=endPointData,
                                 ggplot2::aes(x=.data[["n"]], y=.data[["logEvidence"]]),
                                 colour=endPointData[["border"]], shape=pch, size=2)

  # Quantiles -----
  #
  if (is.null(wantQuantiles) && !is.null(x[["power"]]))
    wantQuantiles <- x[["power"]]

  if (!is.null(wantQuantiles) && !isFALSE(wantQuantiles)) {
    p <- p + ggplot2::annotate("text", x=-Inf, y=textHeightQuant,
                               label="Quantiles [%]", hjust=0,
                               colour=colQuant, size=textSize)

    badQuantileIndeces <- which(wantQuantiles > nAlt/nAll)

    if (length(badQuantileIndeces) >= 1) {
      warning(paste0("The requested quantile(s) ", wantQuantiles[badQuantileIndeces], " cannot be displayed. Max power in this simulation run: ", nAlt/nAll))
      wantQuantiles <- wantQuantiles[-badQuantileIndeces]
    }

    quants <- unname(round(stats::quantile(stoppingTimes, wantQuantiles), 2))

    quantileNames <- round(wantQuantiles*100, 2)

    if (length(quants) >= 1) {
      if (relevanceTest) {
        quantFut <- vapply(quants, function(q) round(sum(fptFut <= q)/nAll*100, 2),
                           numeric(1))

        # The spaces mimic pos=2 and pos=4 of graphics::text()
        p <- p +
          ggplot2::annotate("text", x=quants, y=textHeightQuant,
                            label=paste0(quantileNames, "  "), hjust=1,
                            colour=overColourBorder, size=legendCexFactor*textSize) +
          ggplot2::annotate("text", x=quants, y=textHeightQuant, label="+",
                            colour=colQuant, size=textSize) +
          ggplot2::annotate("text", x=quants, y=textHeightQuant,
                            label=paste0("  ", quantFut), hjust=0,
                            colour=underColourBorder, size=legendCexFactor*textSize)
      } else {
        p <- p + ggplot2::annotate("text", x=quants, y=textHeightQuant,
                                   label=quantileNames,
                                   colour=overColourBorder, size=textSize)
      }

      # TODO(Alexander): perhaps change the lower position as minimum of -log(1/alpha) and log(beta)
      #
      p <- p +
        ggplot2::annotate("segment", x=quants, xend=quants,
                          y=-0.9*log(1/alpha), yend=0.95*textHeightQuant,
                          colour=colQuant, linewidth=thinLineWidth) +
        ggplot2::annotate("text", x=quants, y=-log(1/alpha), label=quants,
                          colour=colQuant, size=textSize)
    }
  }

  # Axes and labels -----
  #
  if (is.null(xlab))
    xlab <- "Sample size"

  if (is.null(ylab))
    ylab <- "Evidence"

  criticalE <- c(if (relevanceTest) alphaRelevance else alpha, 1, 1/alpha)

  xAxis <- makeAxisLimitsAndBreaks(xlim, defaultLim=c(NA, nPlanBatch), start=0)

  subtitle <- if (wantTitle)
    paste0("True ", names(x[["esMin"]]), " = ", round(x[["esTrue"]], 3))

  # The y-axis line runs from the bottom up to 1/alpha, the x-axis line from 0
  # to the right, as in the base graphics version
  p <- p +
    ggplot2::scale_x_continuous(breaks=xAxis[["breaks"]]) +
    ggplot2::scale_y_continuous(breaks=log(criticalE),
                                labels=as.character(criticalE)) +
    ggplot2::guides(x=ggplot2::guide_axis(cap="lower"),
                    y=ggplot2::guide_axis(cap="upper")) +
    ggplot2::coord_cartesian(xlim=xAxis[["limits"]], ylim=ylim) +
    ggplot2::labs(title=main, subtitle=subtitle, x=xlab, y=ylab) +
    ggplot2::theme(axis.title.y=ggplot2::element_text(hjust=0.25,
                                                      margin=ggplot2::margin(r=10)),
                   plot.subtitle=ggplot2::element_text(colour=colQuant,
                                                       size=ggplot2::rel(legendCexFactor)),
                   legend.position=if (wantLegend) "bottom" else "none",
                   legend.text=ggplot2::element_text(size=ggplot2::rel(legendCexFactor)))

  return(p)
}
