# S3 helpers ---------
#
#' Gets the Label of the Test
#'
#' Helper function that outputs the name of the analysis.
#'
#' @param testType A character string. For the t-tests: "oneSample", "paired", "twoSample".
#' @param parameterName The name of the parameter to identify test performed
#'
#' @return Returns a character string with the name of the analysis.
#'
#' @examples
#' safestats:::getNameTestType("oneSample", "t")
getNameTestType <- function(testType, parameterName) {
  nameChar <- switch(testType,
                     "oneSample"="Safe One Sample",
                     "paired"="Safe Paired Sample",
                     "twoSample"="Safe Two Sample",
                     "gLogrank"="Safe Gaussian",
                     "eLogrank"="Safe Exact",
                     "logrank"="Safe",
                     "2x2" = "Safe Test of Two Proportions")

  testName <- switch(parameterName,
                     "phiS"="Z-Test",
                     "deltaS"="T-Test",
                     "thetaS"="Logrank Test")

  return(paste(nameChar, testName))
}

#' Gets the Label of the Alternative Hypothesis
#'
#' Helper function that outputs the alternative hypothesis of the analysis.
#'
#' @param alternative A character string. "two.sided", "greater", "less".
#' @param testType A character string either "oneSample", "paired", "twoSample", "gLogrank", or "eLogrank".
#' @param h0 the value of the null hypothesis
#' @return Returns a character string with the name of the analysis.
#'
#' @examples
#' safestats:::getNameAlternative("two.sided", testType="oneSample")
getNameAlternative <- function(alternative=c("two.sided", "greater", "less"), testType, h0=0) {
  alternative <- match.arg(alternative)

  if (testType == "oneSample") {
    trueMeanStatement <- "true mean"
  } else if (testType %in% c("paired", "twoSample")) {
    trueMeanStatement <- "true difference in means ('x' minus 'y') is"
  } else if (testType == "2x2") {
    trueMeanStatement <- "true difference between proportions in group a and b is"
  } else if (testType %in% c("gLogrank", "eLogrank", "logrank")) {
    trueMeanStatement <- "true hazard ratio is"
  }

  nameChar <- paste(trueMeanStatement, switch(alternative,
                                              "two.sided"= paste("not equal to", h0),
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
#' @param printConfSeq logical, if \code{TRUE} prints also confidence sequence.
#' Default \code{FALSE}
#' @param ... further arguments to be passed to or from methods.
#'
#' @return No returned value, called for side effects.
#' @export
#'
#' @examples
#' safeTTest(rnorm(19), pilot=TRUE)
#' safeZTest(rnorm(19), pilot=TRUE)
print.safeTest <- function (x, digits = getOption("digits"), prefix = "\t",
                            printConfSeq=FALSE, ...) {
  designObj <- x[["designObj"]]

  if (!is.null(x[["testType"]]) && x[["testType"]] != designObj[["testType"]])
    designObj[["testType"]] <- x[["testType"]]

  testType <- designObj[["testType"]]

  analysisName <- getNameTestType("testType"=testType, "parameterName"=names(designObj[["parameter"]]))
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
    out <- character()
    out <- c(out, paste(names(estimate), "=", format(estimate, digits = max(1L, digits - 2L))))

    cat(format(100 * ciValue), "% confidence sequence: ",
        " ", paste0(format(confSeq[1:2], digits = digits),
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

  alphaString <- format(designObj[["alpha"]], digits = max(1L, digits - 2L))
  eValueString <- format(eValue, digits = max(1L, digits - 2L))
  eThresholdString <- format(1/designObj[["alpha"]], digits = max(1L, digits - 2L))

  out <- character()

  if (!is.null(statValue))
    out <- c(out, paste(names(statValue), "=", format(statValue, digits = max(1L, digits - 2L))))

  out <- c(out, paste(names(parameter), "=", format(parameter, digits = max(1L, digits - 2L))))
  cat(paste0("test: ", paste(out, collapse = ", "), sep="\n"))
  cat("e-value =", eValueString, "> 1/alpha =", eThresholdString, ":",
      eValue > 1/designObj[["alpha"]])
  cat("\n")
  cat("alternative hypothesis:", alternativeName, "\n")

  if (printConfSeq) {
    if (!is.null(x[["confSeq"]])) {
      cat(format(100*(1-designObj[["alpha"]])), " percent confidence sequence:\n",
          " ", paste(format(x[["confSeq"]][1:2], digits = digits),
                     collapse = " "), "\n", sep = "")
    }
  }

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
  testType <- designObj[["testType"]]
  parameterName <- names(designObj[["parameter"]])

  note <- designObj[["note"]]
  analysisName <- paste(getNameTestType("testType"=testType, "parameterName"=parameterName), "Design")

  cat("\n")
  cat(strwrap(analysisName, prefix = prefix), sep = "\n")
  cat("\n")

  designObj[["decision rule"]] <- 1/designObj[["alpha"]]

  displayList <- list()

  for (item in c("nPlan", "nEvents", "esMin", "alternative","alternativeRestriction", "beta", "parameter",
                 "alpha", "decision rule", "logImpliedTarget")) {
    itemValue <- designObj[[item]]
    itemValueString <- format(itemValue, digits=digits)

    if (!is.null(itemValue)) {
      if (item == "nPlan") {
        nPlanTwoSe <- designObj[["nPlanTwoSe"]]

        if (!is.null(nPlanTwoSe)) {
          tempNeem <- names(designObj[["nPlan"]])

          for (i in seq_along(itemValue)) {
            if (i==1) {
              itemValueString <- paste0(format(itemValue[i], digits=digits), "\U00B1",
                                        format(nPlanTwoSe[i]))
            } else {
              itemValueString <- paste(itemValueString,
                                       paste0(format(itemValue[i], digits=digits), "\U00B1",
                                              format(nPlanTwoSe[i])),
                                       sep=", ")
            }
          }
          tempNeem <- paste0(names(designObj[["nPlan"]]), "\U00B1", "2se")
          displayList[[paste(tempNeem, collapse=", ")]] <- itemValueString
        } else {
          tempNeem <- names(designObj[["nPlan"]])
          displayList[[paste(tempNeem, collapse=", ")]] <- itemValue
        }
      } else if (item == "nEvents") {
        nEventsTwoSe <- designObj[["nEventsTwoSe"]]
        tempNeem <- names(designObj[["nEvents"]])

        if (!is.null(nEventsTwoSe)) {
          tempNeem <- paste0(tempNeem, "\U00B1", "2se")

          itemValueString <- paste0(format(itemValue, digits=digits), "\U00B1",
                                        format(nEventsTwoSe))
        } else {
          itemValueString <- paste0(format(itemValue, digits=digits))
        }
        displayList[[paste(tempNeem, collapse=", ")]] <- itemValueString
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

#' Prints a safeTSim Object
#'
#' @param x a 'safeTSim' object.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return No returned value, called for side effects.
#'
#' @export
#'
#' @examples
#' designObj <- designSafeT(1)
#'
#' # Data under deltaTrue=deltaMin
#' simObj <- simulate(designObj, nsim=100)
#' simObj
#'
#' # Data under the null deltaTrue=0
#' simObj <- simulate(designObj, nsim=100, deltaTrue=0, freqOptioStop=TRUE, nPlanFreq=10)
#' simObj
print.safeTSim <- function(x, ...) {
  analysisName <- getNameTestType("testType" = x[["testType"]], "parameterName"=names(x[["parameter"]]))

  if(!is.null(x[["safeSim"]])) {
    cat("\n")
    cat("   Simulations for", analysisName, "\n")
    cat("\n")
  }

  cat("Based on nSim =", x[["nsim"]], "and ")

  cat("if the true effect size is \n")
  cat("    deltaTrue =", x[["deltaTrue"]])
  cat("\n")

  cat("then the safe test optimised to detect an effect size of at least: \n")
  cat("    deltaMin =", x[["esMin"]])
  cat("\n")
  cat("with tolerable type I error rate of ")
  cat("\n")
  cat("    alpha =", x[["alpha"]], "and power: 1-beta =", 1-x[["beta"]])
  cat("\n")
  if (length(x[["nPlan"]])==1) {
    cat("for experiments with planned sample size: \n")
    cat("    n1Plan =", x[["nPlan"]])
  } else {
    cat("For experiments with planned sample sizes: \n")
    cat("    n1Plan =", x[["nPlan"]][1], "and n2Plan =", x[["nPlan"]][2])
  }
  cat("\n")

  cat("\n")
  cat("Is estimated to have a null rejection rate of")
  cat("\n")
  cat("    powerAtNPlan =", x[["safeSim"]][["powerAtN1Plan"]])
  cat("\n")
  cat("at the planned sample sizes.")
  cat("\n")

  freqPowerAtN1Plan <- x[["freqSim"]][["powerAtN1Plan"]]

  if (!is.null(freqPowerAtN1Plan)) {
    cat("For the p-value test:    freqPowerAtNPlan =", freqPowerAtN1Plan)
    cat("\n")
  }
  cat("\n")

  cat("Is estimated to have a null rejection rate of ")
  cat("\n")
  cat("    powerOptioStop =", x[["safeSim"]][["powerOptioStop"]])
  cat("\n")
  cat("under optional stopping, and the average stopping time is:")
  cat("\n")

  if (length(x[["nPlan"]]==1)) {
    cat("    n1Mean =", x[["safeSim"]][["nMean"]])
  } else {
    cat("    n1Mean =", x[["safeSim"]][["nMean"]], "and n2Mean =", x[["ratio"]]*x[["safeSim"]][["nMean"]])
  }
  cat("\n")

  freqPowerOptioStop <- x[["freqSim"]][["powerOptioStop"]]
  if (!is.null(freqPowerAtN1Plan)) {
    cat("For the p-value test:    freqPowerOptioStop =", freqPowerOptioStop)
    cat("\n")
  }
}
