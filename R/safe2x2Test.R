#EXPORT --------------------------------------------------------------------
#' Designs a Safe Experiment to Test Two Proportions in Stream Data
#'
#' The design requires the number of observations one expects to collect in each group in each data block.
#' I.e., when one expects balanced data, one could choose \code{na = nb = 1} and would be allowed to analyze
#' the data stream each time a new observation in both groups has come in. The best results in terms of power
#' are achieved when the data blocks are chosen as small as possible, as this allows for analyzing and updating
#' the safe test as often as possible, to fit the data best.
#' Further, the design requires two out of the following three parameters to be known:
#' \itemize{
#'  \item the power one aims to achieve (\code{1 - beta}),
#'  \item the minimal relevant difference between the groups (\code{delta})
#'  \item the number of blocks planned (\code{nBlocksPlan}),
#' }
#' where the unknown out of the three will be estimated. In the case of an exploratory "pilot" analysis,
#' one can also only provide the number of blocks planned.

#'
#' @param na number of observations in group a per data block
#' @param nb number of observations in group b per data block
#' @param nBlocksPlan planned number of data blocks collected
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both "nBlocksPlan"
#' and "delta". Note that 1-beta defines the power.
#' @param delta a priori minimal relevant divergence between group means b and a, either a numeric between -1 and 1 for
#' no alternative restriction or a restriction on difference, or a real for a restriction on the log odds ratio.
#' @param alternativeRestriction a character string specifying an optional restriction on the alternative hypothesis; must be one of "none" (default),
#' "difference" (difference group mean b minus group b) or "logOddsRatio" (the log odds ratio between group means b and a).
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control --independent on n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule e10 > 1/alpha.
#' @param pilot logical, specifying whether it's a pilot design.
#' @param hyperParameterValues named list containing numeric values for hyperparameters betaA1, betaA2, betaB1 and betaB2, with betaA1 and betaB1 specifying the parameter
#' equivalent to \code{shape1} in \code{stats::dbeta} for groups A and B, respectively, and betaA2 and betaB2 equivalent to \code{shape2}. By default
#' chosen to optimize evidence collected over subsequent experiments (REGRET). Pass in the following format:
#' \code{list(betaA1 = numeric1, betaA2 = numeric2, betaB1 = numeric3, betaB2 = numeric4)}.
#' @param previousSafeTestResult optionally, a previous safe test result can be provided. The posterior
#' of the hyperparameters of this test is then used for the hyperparameter settings. Default NULL.
#' @param M number of simulations used to estimate power or nBlocksPlan. Default \code{1000}.
#'
#' @return Returns a 'safeDesign' object that includes:
#'
#' \describe{
#'   \item{nPlan}{the sample size(s) to plan for. Computed based on beta and meanDiffMin, or provided by the user
#'   if known.}
#'   \item{parameter}{the safe test defining parameter: here the hyperparameters.}
#'   \item{esMin}{the minimally clinically relevant effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error specified by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" based on the \code{alternativeRestriction} provided by the user.}
#'   \item{testType}{here 2x2}
#'   \item{pilot}{logical, specifying whether it's a pilot design.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' #plan for an experiment to detect minimal difference of 0.6 with a balanced design
#' set.seed(3152021)
#' designSafeTwoProportions(na = 1,
#'                          nb = 1,
#'                          alpha = 0.1,
#'                          beta = 0.20,
#'                          delta = 0.6,
#'                          alternativeRestriction = "none",
#'                          M = 1e2)
#'
#' #safe analysis of a pilot: number of samples already known
#' designSafeTwoProportions(na = 1,
#'                           nb = 1,
#'                           nBlocksPlan = 20,
#'                           pilot = TRUE)
#'
#' #specify own hyperparameters
#' hyperParameterValues <- list(betaA1 = 10, betaA2 = 1, betaB1 = 1, betaB2 = 10)
#' designSafeTwoProportions(na = 1,
#'                          nb = 1,
#'                          alpha = 0.1,
#'                          beta = 0.20,
#'                          delta = 0.6,
#'                          hyperParameterValues = hyperParameterValues,
#'                          alternativeRestriction = "none",
#'                          M = 1e2)
#'
#'
designSafeTwoProportions <- function(na, nb,
                                     nBlocksPlan = NULL,
                                     beta = NULL,
                                     delta = NULL,
                                     alternativeRestriction = c("none", "difference", "logOddsRatio"),
                                     alpha = 0.05,
                                     pilot = "FALSE",
                                     hyperParameterValues = NULL,
                                     previousSafeTestResult = NULL,
                                     M = 1e3){
  alternativeRestriction <- match.arg(alternativeRestriction)

  if (alternativeRestriction %in% c("difference", "logOddsRatio") & !is.numeric(delta)) {
    stop("Provide numeric value for divergence measure when testing with restriction: a difference or log Odds ratio")
  }

  if (is.null(hyperParameterValues)) {
    #use the default
    hyperParameterValues <- list(betaA1 = 0.18, betaB1 = (nb/na)*0.18,
                        betaA2 = 0.18, betaB2 = (nb/na)*0.18)
    priorValuesForPrint <- "standard, REGRET optimal"
    note <- "Optimality of hyperparameters only verified for equal group sizes (na = nb = 1)"
  } else {
    #user provided: perform checks
    if (!all(c("betaA1", "betaA2", "betaB1", "betaB2") %in% names(hyperParameterValues))) {
      stop("Provide hyperparameters as a named list for betaA1, betaA2, betaB1 and betaB2, see help file.")
    }

    if (any(hyperParameterValues <= 0)) {
      stop("Provide Beta prior hyperparameter values that yield a proper prior: parameters should be > 0.")
    }
    priorValuesForPrint <- paste(hyperParameterValues, collapse = " ")
    note <- NULL
  }

  if (!is.null(previousSafeTestResult)) {
    #use posterior for hyperparameter settings
    hyperParameterValues <- previousSafeTestResult[["posteriorHyperParameters"]]
    priorValuesForPrint <- paste(hyperParameterValues, collapse = " ")
    note <- "Hyperparameters set according to posterior values from previous test result"
  }

  names(priorValuesForPrint) <- "Beta hyperparameters"
  impliedTarget <- NULL

  #Check each possible design scenario
  if (is.null(nBlocksPlan) && (is.numeric(delta) && delta != 0) && !is.null(beta)) {
    #scenario 1a: delta + power known, calculate nPlan
    nBlocksPlan <- simulateWorstCaseQuantileTwoProportions(delta = delta,
                                            na = na, nb = nb,
                                            priorValues = hyperParameterValues,
                                            alternativeRestriction = alternativeRestriction,
                                            alpha = alpha, beta = beta, M = M)[["worstCaseQuantile"]]
  } else if (!is.null(nBlocksPlan) && !(is.numeric(delta) && delta != 0) && is.null(beta)) {
    #scenario 1c: only nPlan known, can perform a pilot (no warning though)
    pilot <- TRUE
  } else if (!is.null(nBlocksPlan) && (is.numeric(delta) && delta != 0) && is.null(beta)) {
    #scenario 2: given effect size and nPlan, calculate power
    worstCaseSimulationResult <- simulateWorstCaseQuantileTwoProportions(delta = delta,
                                                                         na = na, nb = nb,
                                                                         priorValues = hyperParameterValues,
                                                                         alternativeRestriction = alternativeRestriction,
                                                                         maxSimStoptime = nBlocksPlan,
                                                                         alpha = alpha, beta = 0, M = M)
    beta <- 1 - worstCaseSimulationResult[["worstCasePower"]]
    impliedTarget <- worstCaseSimulationResult[["impliedTarget"]]
  } else if (!is.null(nBlocksPlan) && !(is.numeric(delta) && delta != 0) && !is.null(beta)) {
    #scenario 3: given power and nPlan, calculate minimal effect size to be "detected"
    delta <- simulateWorstCaseDeltaTwoProportions(na = na, nb = nb, priorValues = hyperParameterValues,
                                    alternativeRestriction = alternativeRestriction,
                                    alpha = alpha, beta = beta, M = M, maxSimStoptime = nBlocksPlan)
    if (length(delta) == 0) {
      stop("For this sample size and power, no effect size below deltamax yielded the desired power")
    }
  } else {
    #also includes scenario 1b: only delta known, raise error
    stop("Provide two of nBlocksPlan, delta and power, or only nBlocksPlan for a pilot with default settings.")
  }

  #in the scenario's we accept now, there always is an nBlocksPlan not null
  nPlan <- c(na, nb, nBlocksPlan)
  names(nPlan) <- c("na", "nb", "nBlocksPlan")

  alternative <- switch(alternativeRestriction,
                        "none" = "two.sided",
                        "difference" = ifelse(delta < 0, "less", "greater"),
                        "logOddsRatio" = ifelse(delta < 0, "less", "greater"))

  if (!is.null(delta)) {
    names(delta) <- ifelse(alternativeRestriction == "logOddsRatio", "log odds ratio", "difference")
  }

  testType <- "2x2"

  #might change for the confidence sequences
  h0 <- 0

  result <- list("nPlan"=nPlan,
                 "parameter"= priorValuesForPrint,
                 "betaPriorParameterValues" = hyperParameterValues,
                 "alpha"=alpha,
                 "beta"=beta,
                 "impliedTarget" = impliedTarget,
                 "esMin" = delta,
                 "h0"= h0,
                 "testType"=testType,
                 "alternativeRestriction" = alternativeRestriction,
                 "alternative" = alternative,
                 "pilot" = pilot,
                 "lowN"=NULL,
                 "highN"=NULL,
                 "call"=sys.call(),
                 "timeStamp"=Sys.time(),
                 "note" = note)
  class(result) <- "safeDesign"

  return(result)
}

#' Perform a Safe Test for Two Proportions with Stream Data
#'
#' Perform a safe test for two proportions (a 2x2 contingency table test) with a
#' result object retrieved through the design function for planning an experiment to compare
#' two proportions in this package, \code{\link{designSafeTwoProportions}()}.
#'
#' @param ya positive observations/ events per data block in group a: a numeric with integer values
#' between (and including) 0 and \code{na}, the number of observations in group a per block.
#' @param yb positive observations/ events per data block in group b: a numeric with integer values
#' between (and including) 0 and \code{nb}, the number of observations in group b per block.
#' @param designObj a safe test design for two proportions retrieved through \code{\link{designSafeTwoProportions}()}.
#' @param pilot logical that can be set to true when performing an exploratory analysis
#' without a \code{designObj}; only allows for \code{na = nb = 1}.
#'
#' @return Returns an object of class 'safeTest'. An object of class 'safeTest' is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{n}{The realised sample size(s).}
#'   \item{eValue}{the e-value of the safe test.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" described in \code{\link{designSafeTwoProportions}()}.}
#' }
#'
#' @export
#'
#' @examples
#' #balanced design
#' yb <- c(1,0,1,1,1,0,1)
#' ya <- c(1,0,1,0,0,0,1)
#' safeDesign <- designSafeTwoProportions(na = 1,
#'                                        nb = 1,
#'                                        beta = 0.20,
#'                                        delta = 0.6,
#'                                        alternativeRestriction = "none",
#'                                        M = 1e1)
#' safeTwoProportionsTest(ya = ya, yb = yb, designObj = safeDesign)
#'
#' #pilot
#' safeTwoProportionsTest(ya = ya, yb = yb, pilot = TRUE)
#'
#' #unbalanced design
#' yb <- c(1,0,1,1,1,0,1)
#' ya <- c(2,2,1,2,0,2,2)
#' safeDesign <- designSafeTwoProportions(na = 2,
#'                                        nb = 1,
#'                                        beta = 0.20,
#'                                        delta = 0.6,
#'                                        alternativeRestriction = "none",
#'                                        M = 1e1)
#' safeTwoProportionsTest(ya = ya, yb = yb, designObj = safeDesign)
#'
safeTwoProportionsTest <- function(ya, yb, designObj = NULL, pilot = FALSE) {
  if (is.null(designObj) & !pilot) {
    stop("Please provide a safe 2x2 design object, or run the function with pilot=TRUE.",
         "A design object can be obtained by running designSafeTwoProportions().")
  }

  if (length(ya) != length(yb)) {
    stop("Can only process complete data blocks: provide vectors with numbers of positive observations per timepoint,",
         "see example in helpfile.")
  }

  if (pilot) {
    designObj <- designSafeTwoProportions(na = 1, nb = 1, nBlocksPlan = length(ya),
                                          alternativeRestriction = "none", pilot = TRUE)
  }

  if (any(ya > designObj[["nPlan"]][["na"]] | ya < 0) | any(yb > designObj[["nPlan"]][["nb"]] | yb < 0)) {
    stop("Provided sample sizes within blocks, na and nb, do not match provided ya and yb.")
  }

  eValue <- calculateSequential2x2E(aSample = ya, bSample = yb,
                                    priorValues = designObj[["betaPriorParameterValues"]],
                                    restriction = designObj[["alternativeRestriction"]],
                                    delta = designObj[["esMin"]],
                                    na = designObj[["nPlan"]][["na"]],
                                    nb = designObj[["nPlan"]][["nb"]])

  argumentNames <- getArgs()
  xLabel <- extractNameFromArgs(argumentNames, "ya")
  yLabel <- extractNameFromArgs(argumentNames, "yb")
  dataName <- paste(xLabel, "and", yLabel)
  n <- c(length(ya)*designObj[["nPlan"]][["na"]], length(yb)*designObj[["nPlan"]][["nb"]])
  names(n) <- c("nObsA", "nObsB")

  #calculate the posterior: prior parameters from original design, plus successes and failures
  #seen in this experiment
  posteriorHyperParameters <- list(
    betaA1 = sum(ya) + designObj[["betaPriorParameterValues"]][["betaA1"]],
    betaB1 = length(ya)*designObj[["nPlan"]][["na"]] - sum(ya) + designObj[["betaPriorParameterValues"]][["betaA2"]],
    betaA2 = sum(yb) + designObj[["betaPriorParameterValues"]][["betaB1"]],
    betaB2 = length(yb)*designObj[["nPlan"]][["nb"]] - sum(yb) + designObj[["betaPriorParameterValues"]][["betaB2"]]
  )

  testResult <- list(designObj = designObj,
                     eValue = eValue,
                     dataName = dataName,
                     n = n,
                     posteriorHyperParameters = posteriorHyperParameters)
  class(testResult) <- "safeTest"

  return(testResult)
}

#' Alias for \code{\link{safeTwoProportionsTest}()}
#'
#' @rdname safeTwoProportionsTest
#'
#' @export
safe.prop.test <- function(ya, yb, designObj = NULL, pilot = FALSE) {
  safeTestResult <- tryCatch(safeTwoProportionsTest(ya = ya, yb = yb, designObj = designObj, pilot = pilot),
                  error = function(e){e})

  if (!is.null(safeTestResult[["message"]])) {
    #safeTwoProportionsTest has thrown an error - return neatly, as from this call
    stop(safeTestResult[["message"]])
  } else {
    return(safeTestResult)
  }
}

#' Compare Different Hyperparameter Settings for Safe Tests of Two Proportions.
#'
#' Simulates for a range of divergence parameter values (differences or log odds ratios) the worst-case stopping times
#' (i.e., number of data blocks collected) and expected stopping times needed to achieve the desired power for each hyperparameter setting provided.
#'
#' @inheritParams designSafeTwoProportions
#' @param hyperparameterList list object, its components hyperparameter lists with a format as described in \code{\link{designSafeTwoProportions}()}.
#' @param deltaDesign optional; when using a restricted alternative, the value of the divergence measure used.
#' Either a numeric between -1 and 1 for a restriction on difference, or a real for a restriction on the log odds ratio.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control in the study. Necessary to calculate the
#' worst case stopping time.
#' @param deltamax maximal effect size to calculate power for; between -1 and 1 for designs without restriction or a restriction on difference;
#' real number for a restriction on the log odds ratio. Default \code{0.9}.
#' @param deltamin minimal effect size to calculate power for; between -1 and 1 for designs without restriction or a restriction on difference;
#' real number for a restriction on the log odds ratio. Default \code{0.1}.
#' @param deltaGridSize numeric, positive integer: size of grid of delta values worst case and expected sample sizes are simulated for.
#' @param M number of simulations used to estimate sample sizes. Default \code{100}.
#' @param maxSimStoptime maximal stream length in simulations; when the e value does not reach the rejection threshold before the end of the stream,
#' the maximal stream length is returned as the stopping time. Default \code{1e4}.
#' @param thetaAgridSize numeric, positive integer: size of the grid of probability distributions examined for each delta value to find the
#' worst case sample size over.
#'
#' @return Returns an object of class "safe2x2Sim". An object of class "safe2x2Sim" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{simData}{A data frame containing simulation results with worst case and expected stopping times for each
#'   hyperparameter setting, for the specified or default range of effect sizes.}
#'   \item{alpha}{the significance threshold used in the simulations}
#'   \item{beta}{the type-II error control used in the simulations}
#'   \item{deltaDesign}{the value of restriction on the alternative hypothesis parameter space used for the E variables in the simulations}
#'   \item{restriction}{the type of restriction used for the E variables in the simulation}
#'   \item{hyperparameters}{list of the hyperparameters tested in the simulation}
#' }
#'
#' @export
#'
#' @examples
#' priorList1 <- list(betaA1 = 10, betaA2 = 1, betaB1 = 1, betaB2 = 10)
#' priorList2 <- list(betaA1 = 0.18, betaA2 = 0.18, betaB1 = 0.18, betaB2 = 0.18)
#' priorList3 <- list(betaA1 = 1, betaA2 = 1, betaB1 = 1, betaB2 = 1)
#'
#' simResult <- simulateTwoProportions(
#'   hyperparameterList = list(priorList1, priorList2, priorList3),
#'   alternativeRestriction = "none",
#'   alpha = 0.1, beta = 0.2, na = 1, nb = 1,
#'   deltamax = -0.4, deltamin = -0.9, deltaGridSize = 3,
#'   M = 10
#'   )
#'
#' print(simResult)
#' plot(simResult)
simulateTwoProportions <- function(hyperparameterList,
                                   alternativeRestriction = c("none", "difference", "logOddsRatio"),
                                   deltaDesign = NULL,
                                   alpha,
                                   beta,
                                   na, nb,
                                   deltamax = 0.9,
                                   deltamin = 0.1,
                                   deltaGridSize = 8,
                                   M = 1e2,
                                   maxSimStoptime = 1e4,
                                   thetaAgridSize = 8){

  deltaVec <- seq(deltamax, deltamin, length.out = deltaGridSize)
  resultDataFrame <- data.frame()

  #use list names to index and return the result
  if (is.null(names(hyperparameterList))) {
    names(hyperparameterList) <- paste("setting", 1:length(hyperparameterList), sep = "")
  }

  #for every prior, simulate the worst-case 1 - beta stopping times
  #for a grid of delta values
  for (hyperparameterSet in names(hyperparameterList)) {
    message(paste("Retrieving worst case stopping times for hyperparameter set ", hyperparameterSet))
    for (deltaGenerating in deltaVec) {
      worstCase <- simulateWorstCaseQuantileTwoProportions(na = na, nb = nb,
                                 priorValues = hyperparameterList[[hyperparameterSet]],
                                 alternativeRestriction = alternativeRestriction,
                                 alpha = alpha, beta = beta,
                                 delta = deltaGenerating,
                                 deltaDesign = deltaDesign,
                                 M = M,
                                 maxSimStoptime = maxSimStoptime,
                                 gridSize = thetaAgridSize)[["worstCaseQuantile"]]
      resultDataFrame <- rbind(resultDataFrame,
                               data.frame(hyperparameters = hyperparameterSet,
                                          delta = deltaGenerating,
                                          worstCaseQuantile = worstCase)
                               )
    }
  }

  #now we know the worst case quantiles, retrieve expected stopping times
  message(paste("Retrieving all expected stopping times given the worst case stopping times"))
  for (i in 1:nrow(resultDataFrame)) {
    resultDataFrame[i,"expected"] <- simulateWorstCaseQuantileTwoProportions(na = na, nb = nb,
                                          priorValues = hyperparameterList[[resultDataFrame[i, "hyperparameters"]]],
                                          alternativeRestriction = alternativeRestriction,
                                          alpha = alpha, beta = 0,
                                          delta = resultDataFrame[i, "delta"],
                                          deltaDesign = deltaDesign,
                                          M = M,
                                          #now we stop, each experiment,
                                          #maximally at the worst case stoptime for 80% power
                                          maxSimStoptime = ceiling(resultDataFrame[i,"worstCaseQuantile"]),
                                          gridSize = thetaAgridSize,
                                          #and we calculate the expectation instead of a (1-b) quantile
                                          expectedStopTime = TRUE)[["worstCaseQuantile"]]
  }

  simResult <- list(simdata = resultDataFrame,
                    alpha = alpha,
                    beta = beta,
                    deltaDesign = deltaDesign,
                    restriction = alternativeRestriction,
                    hyperparameters = hyperparameterList
                    )
  class(simResult) <- "safe2x2Sim"
  return(simResult)
}

#' Prints Results of Simulations for Comparing Hyperparameters for Safe Tests of Two Proportions
#'
#' @param x a result object obtained through \code{\link{simulateTwoProportions}()}.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return The data frame with simulation results, called for side effects to pretty print the simulation results.
#'
#' @export
#'
#' @examples
#' priorList1 <- list(betaA1 = 10, betaA2 = 1, betaB1 = 1, betaB2 = 10)
#' priorList2 <- list(betaA1 = 0.18, betaA2 = 0.18, betaB1 = 0.18, betaB2 = 0.18)
#' priorList3 <- list(betaA1 = 1, betaA2 = 1, betaB1 = 1, betaB2 = 1)
#'
#' simResult <- simulateTwoProportions(
#'   hyperparameterList = list(priorList1, priorList2, priorList3),
#'   alternativeRestriction = "none",
#'   alpha = 0.1, beta = 0.2, na = 1, nb = 1,
#'   deltamax = -0.4, deltamin = -0.9, deltaGridSize = 3,
#'   M = 10
#'   )
print.safe2x2Sim <- function(x, ...){
  cat("Simulation results for test of two proportions")
  cat("\n\n")

  cat("Simulations ran with alpha =", x[["alpha"]], "and beta =", x[["beta"]], ".\n")
  if (!is.null(x[["deltaDesign"]])) {
    cat("The alternative hypothesis was restricted based on a",
        x[["restriction"]], "of", x[["deltaDesign"]], ".\n")
  }

  cat("The following hyperparameter settings were evaluated:\n")
  displayList <- list()
  for (i in 1:length(x[["hyperparameters"]])) {
    displaytext <- paste(names(x[["hyperparameters"]][[i]]), unlist(x[["hyperparameters"]][[i]]), sep = " = ",
                         collapse = "; ")
    displayList[[names(x[["hyperparameters"]])[i]]] <- displaytext
  }
  cat(paste(format(names(displayList), width = 20L, justify = "right"),
            format(displayList), sep = ": "), sep = "\n")
  cat("and yielded the following results:\n")
  print(x[["simdata"]], justify = "right")
}

#' Plots Results of Simulations for Comparing Hyperparameters for Safe Tests of Two Proportions
#'
#' @param x a result object obtained through \code{\link{simulateTwoProportions}()}.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Plot data, mainly called for side effects, the plot of simulation results.
#'
#' @export
#'
#' @examples
#' priorList1 <- list(betaA1 = 10, betaA2 = 1, betaB1 = 1, betaB2 = 10)
#' priorList2 <- list(betaA1 = 0.18, betaA2 = 0.18, betaB1 = 0.18, betaB2 = 0.18)
#' priorList3 <- list(betaA1 = 1, betaA2 = 1, betaB1 = 1, betaB2 = 1)
#'
#' simResult <- simulateTwoProportions(
#'   hyperparameterList = list(priorList1, priorList2, priorList3),
#'   alternativeRestriction = "none",
#'   alpha = 0.1, beta = 0.2, na = 1, nb = 1,
#'   deltamax = -0.4, deltamin = -0.9, deltaGridSize = 3,
#'   M = 10
#'   )
#'
#' plot(simResult)
#'
plot.safe2x2Sim <- function(x, ...){
  if (is.null(x[["deltaDesign"]])) {
    mainTitle <- "Worst case and expected stopping times without restriction on H1"
  } else {
    mainTitle <- paste("Worst case and expected stopping times with a restriction on the",
                               x[["restriction"]], "of", round(x[["deltaDesign"]],2))
  }
  subTitle <- bquote(alpha == .(x[["alpha"]]) ~"," ~ beta == .(x[["beta"]]))

  xmin <- min(x[["simdata"]][,"delta"])
  xmax <- max(x[["simdata"]][,"delta"])
  ymin <- 0
  ymax <- ceiling(max(x[["simdata"]][,c("worstCaseQuantile", "expected")]))

  xlab <- paste("divergence value:", ifelse(x[["restriction"]] == "logOddsRatio", "log odds ratio", "difference"))

  graphics::plot(x = 1, type = "n",
                 xlim = c(xmin, xmax),
                 ylim = c(ymin, ymax),
                 xlab = xlab,
                 ylab = "stopping time (m collected)",
                 main = mainTitle,
                 sub = subTitle,
                 col = "lightgrey")

  priorcolors <- grDevices::rainbow(length(unique(x[["simdata"]][,"hyperparameters"])))
  names(priorcolors) <- unique(x[["simdata"]][,"hyperparameters"])

  #first, add the worst case stopping times
  for (hyperparameters in unique(x[["simdata"]][,"hyperparameters"])) {
    plotData <- x[["simdata"]][x[["simdata"]]$hyperparameters == hyperparameters,]
    linecolor <- priorcolors[hyperparameters]
    graphics::lines(x = plotData[,"delta"], y = plotData[,"worstCaseQuantile"], col = linecolor, lty = 2, lwd = 2)
  }

  #then, add the expected stopping times
  for (hyperparameters in unique(x[["simdata"]][,"hyperparameters"])) {
    plotData <- x[["simdata"]][x[["simdata"]]$hyperparameters == hyperparameters,]
    linecolor <- priorcolors[hyperparameters]
    graphics::lines(x = plotData[,"delta"], y = plotData[,"expected"], col = linecolor, lty = 1, lwd = 2)
  }

  graphics::legend(x = "topright", legend = c(names(priorcolors), "worst case", "expected"),
                   col = c(priorcolors, "grey", "grey"),
                   lty = c(rep(2, length(priorcolors)), 2, 1), lwd = 2)
}

#' Estimate Lower and Upper Bounds on the Confidence Sequence (Interval)
#' for the Difference Divergence Measure for Two Proportions
#'
#' @param ya positive observations/ events per data block in group a: a numeric with integer values
#' between (and including) 0 and \code{na}, the number of observations in group a per block.
#' @param yb positive observations/ events per data block in group b: a numeric with integer values
#' between (and including) 0 and \code{nb}, the number of observations in group b per block.
#' @param precision precision of the grid of differences to search over for the lower and upper bounds.
#' @param safeDesign a 'safeDesign' object obtained through
#' \code{\link{designSafeTwoProportions}}
#'
#' @return list with found lower and upper bound.
#' @export
#' @importFrom purrr %>%
#' @importFrom rlang .data
#'
#' @examples
#' balancedSafeDesign <- designSafeTwoProportions(na = 1,
#'                                                nb = 1,
#'                                                nBlocksPlan = 10,
#'                                                alpha = 0.05)
#' ya <- c(1,1,1,1,1,1,1,1,0,1)
#' yb <- c(0,0,0,0,1,0,0,0,0,0)
#' getConfidenceBoundsForDifferenceTwoProportions(ya = ya,
#'                                                yb = yb,
#'                                                precision = 20,
#'                                                safeDesign = balancedSafeDesign)
#'
getConfidenceBoundsForDifferenceTwoProportions <- function(ya,
                                                           yb,
                                                           precision,
                                                           safeDesign){
  na <- safeDesign[["nPlan"]][["na"]]
  nb <- safeDesign[["nPlan"]][["nb"]]
  alpha <- safeDesign[["alpha"]]

  eValuesDeltaGrid <- calculateEValuesForLinearDeltaGrid(ya = ya, yb = yb,
                                                         na = na, nb = nb,
                                                         priorParameters = safeDesign[["betaPriorParameterValues"]],
                                                         precision = precision,
                                                         alpha = alpha,
                                                         runningIntersection = TRUE)
  #include in the CS: not rejected delta values
  ciSummary <- eValuesDeltaGrid %>%
    dplyr::filter(.data[["E"]] < 1/alpha) %>%
    dplyr::summarise(lowerBound = min(.data[["delta"]]), upperBound = max(.data[["delta"]]))

  return(as.list(ciSummary))
}

#' Estimate an upper or lower bound for a safe confidence sequence on the
#' logarithm of the odds ratio for two proportions.
#'
#' @param ya positive observations/ events per data block in group a: a numeric with integer values
#' between (and including) 0 and \code{na}, the number of observations in group a per block.
#' @param yb positive observations/ events per data block in group b: a numeric with integer values
#' between (and including) 0 and \code{nb}, the number of observations in group b per block.
#' @param safeDesign a 'safeDesign' object obtained through
#' \code{\link{designSafeTwoProportions}}
#' @param bound type of bound to calculate; "lower" to get a lower bound on positive delta,
#' "upper" to get an upper bound on negative delta.
#' @param deltaStart starting value of the grid to search over for the bound on the confidence
#' sequence (in practice: the interval). Numeric >0 when searching for a lower bound, numeric < 0
#' when searching for an upper bound.
#' @param deltaStop end value of the grid to search over for the bound on the confidence
#' sequence (in practice: the interval). Numeric >0 when searching for a lower bound, numeric < 0
#' when searching for an upper bound.
#' @param precision precision of the grid between deltaStart and deltaStop.
#'
#' @return numeric: the established lower- or upper bound on the logarithm of the odds
#' ratio between the groups
#'
#' @export
#' @importFrom purrr %>%
#' @importFrom rlang .data
#'
#' @examples
#' balancedSafeDesign <- designSafeTwoProportions(na = 1,
#'                                                nb = 1,
#'                                                nBlocksPlan = 10,
#'                                                alpha = 0.05)
#' #hypothesize OR < 1 (i.e., log OR < 0)
#' ya <- c(1,1,1,1,1,1,1,1,0,1)
#' yb <- c(0,0,0,0,1,0,0,0,0,0)
#' #one-sided CI for OR-, establish upper bound on log odds ratio
#' getConfidenceBoundForLogOddsTwoProportions(ya = ya,
#'                                            yb = yb,
#'                                            safeDesign = balancedSafeDesign,
#'                                            bound = "upper",
#'                                            deltaStart = -0.01,
#'                                            deltaStop = -4,
#'                                            precision = 20)
#'
getConfidenceBoundForLogOddsTwoProportions <- function(ya,
                                                       yb,
                                                       safeDesign,
                                                       bound = c("lower", "upper"),
                                                       deltaStart,
                                                       deltaStop,
                                                       precision){
  na <- safeDesign[["nPlan"]][["na"]]
  nb <- safeDesign[["nPlan"]][["nb"]]
  priorParameters <- safeDesign[["betaPriorParameterValues"]]
  alpha <- safeDesign[["alpha"]]
  bound = match.arg(bound)
  lowerBound = ifelse(bound == "lower", TRUE, FALSE)

  eValuesDeltaGrid <- calculateEValuesForOddsDeltaGrid(ya = ya, yb = yb,
                                                       na = na, nb = nb,
                                                       lowerBound = lowerBound,
                                                       priorParameters = priorParameters,
                                                       precision = precision,
                                                       deltaStart = deltaStart,
                                                       deltaStop = deltaStop,
                                                       alpha = alpha)

  if (all(eValuesDeltaGrid$E < 1/alpha)) {
    warning(paste("No", bound, "bound could be established; try different bound or smaller deltaStart"))
    return(0)
  }

  if (lowerBound) {
    deltaBound <- as.numeric(eValuesDeltaGrid %>%
                               #lower bound: the smallest delta we did not reject
                               dplyr::group_by(.data[["delta"]]) %>%
                               dplyr::filter(.data[["block"]] == max(.data[["block"]])) %>%
                               dplyr::ungroup() %>%
                               dplyr::filter(.data[["E"]] < 1/alpha) %>%
                               dplyr::summarise(bound = min(.data[["delta"]]))
    )
  } else {
    deltaBound <- as.numeric(eValuesDeltaGrid %>%
                               #upper bound: the biggest delta we did not reject
                               dplyr::group_by(.data[["delta"]]) %>%
                               dplyr::filter(.data[["block"]] == max(.data[["block"]])) %>%
                               dplyr::ungroup() %>%
                               dplyr::filter(.data[["E"]] < 1/alpha) %>%
                               dplyr::summarise(bound = max(.data[["delta"]]))
    )
  }

  return(deltaBound)
}

#vignette-----------------------------------------------------------------------
#' Simulate an optional stopping scenario according to a safe design for two proportions
#'
#' @param safeDesign a 'safeDesign' object obtained through \code{\link{designSafeTwoProportions}()}.
#' @param M integer, the number of data streams to sample.
#' @param thetaA Bernoulli distribution parameter in group A
#' @param thetaB Bernoulli distribution parameter in group B
#'
#' @return list with the simulation results of the safe test under optional stopping with the following
#' components:
#'
#' \describe{
#'   \item{powerOptioStop}{Proportion of sequences where H0 was rejected}
#'   \item{nMean}{Mean stopping time}
#'   \item{probLessNDesign}{Proportion of experiments stopped before nBlocksPlan was reached}
#'   \item{lowN}{Minimum stopping time}
#'   \item{eValues}{All achieved E values}
#'   \item{allN}{All stopping times}
#'   \item{allSafeDecisions}{Decisions on rejecting H0 for each M}
#'   \item{allRejectedN}{Stopping times of experiments where H0 was rejected}
#' }
#'
#' @export
#'
#' @examples
#' balancedSafeDesign <- designSafeTwoProportions(na = 1,
#'                                                nb = 1,
#'                                                nBlocksPlan = 30)
#' optionalStoppingSimulationResult <- simulateOptionalStoppingScenarioTwoProportions(
#'   safeDesign = balancedSafeDesign,
#'   M = 1e2,
#'   thetaA = 0.2,
#'   thetaB = 0.5
#' )
simulateOptionalStoppingScenarioTwoProportions <- function(safeDesign,
                                                           M,
                                                           thetaA,
                                                           thetaB){

  stoppingTimes <- stopEs <- numeric(M)

  for (i in 1:M) {
    #For every m, draw a sample of max streamlength and record the time
    #at which we would have stopped
    ya <- rbinom(n = safeDesign[["nPlan"]]["nBlocksPlan"],
                 size = safeDesign[["nPlan"]]["na"],
                 prob = thetaA
    )
    yb <- rbinom(n = safeDesign[["nPlan"]]["nBlocksPlan"],
                 size = safeDesign[["nPlan"]]["nb"],
                 prob = thetaB
    )
    simResult <- calculateSequential2x2E(aSample = ya, bSample = yb,
                                         priorValues = safeDesign[["betaPriorParameterValues"]],
                                         restriction = safeDesign[["alternativeRestriction"]],
                                         #if explicitly passsed deltaDesign (neq delta), use that one for test
                                         #e.g. when studying effect of overestimated/ underestimated effect size
                                         delta = safeDesign[["esMin"]],
                                         na = safeDesign[["nPlan"]]["na"],
                                         nb = safeDesign[["nPlan"]]["nb"],
                                         simSetting = TRUE,
                                         alphaSim = safeDesign[["alpha"]])
    stoppingTimes[i] <- simResult[["stopTime"]]
    stopEs[i] <- simResult[["stopE"]]
  }

  allSafeDecisions <- stopEs >= (1/safeDesign[["alpha"]])
  safeSim <- list("powerOptioStop"= mean(allSafeDecisions),
                  "nMean"= mean(stoppingTimes),
                  "probLessNDesign"= mean(stoppingTimes < safeDesign[["nPlan"]]["nBlocksPlan"]),
                  "lowN"= min(stoppingTimes),
                  "eValues"=stopEs
  )

  safeSim[["allN"]] <- stoppingTimes
  safeSim[["allSafeDecisions"]] <- allSafeDecisions
  safeSim[["allRejectedN"]] <- stoppingTimes[allSafeDecisions]

  return(safeSim)
}

#' Simulate incorrect optional stopping with fisher's exact test's p-value as the
#' stopping rule.
#'
#' @param thetaA Bernoulli distribution parameter in group A
#' @param thetaB Bernoulli distribution parameter in group B
#' @param alpha Significance level
#' @param na number of observations in group a per data block
#' @param nb number of observations in group b per data block
#' @param maxSimStoptime maximal number of blocks to sample in each experiment
#' @param M Number of simulations to carry out, deafult 1e3.
#' @param numberForSeed number for seed to set, default NULL.
#'
#' @return list with stopping times and rejection decisions.
#' @export
#'
#' @examples
#' simulateIncorrectStoppingTimesFisher(thetaA = 0.3,
#'                                      thetaB = 0.3,
#'                                      alpha = 0.05,
#'                                      na = 1,
#'                                      nb = 1,
#'                                      M = 10,
#'                                      maxSimStoptime = 100,
#'                                      numberForSeed = 251)
simulateIncorrectStoppingTimesFisher <- function(thetaA, thetaB, alpha,
                                                 na, nb,
                                                 maxSimStoptime = 1e4,
                                                 M = 1e3, numberForSeed = NULL){

  #setup
  stoppingTimes <- rejections <- numeric(M)
  numberForSeed <- ifelse(is.null(numberForSeed), Sys.time(), numberForSeed)
  set.seed(numberForSeed)

  groupSizeVecA <- (1:maxSimStoptime)*na
  groupSizeVecB <- (1:maxSimStoptime)*nb

  for (m in 1:M) {

    #simulate data
    ya <- rbinom(n = maxSimStoptime, size = na, prob = thetaA)
    yb <- rbinom(n = maxSimStoptime, size = nb, prob = thetaB)

    successAVec <- cumsum(ya)
    successBVec <- cumsum(yb)

    failAVec <- groupSizeVecA - successAVec
    failBVec <- groupSizeVecB - successBVec
    for (i in 5:maxSimStoptime) {
      #new data come in
      successA <- successAVec[i]
      failA <- failAVec[i]
      successB <- successBVec[i]
      failB <- failBVec[i]

      #Fisher's exact test with all data seen so far
      pVal <- tryCatch(fisher.test(matrix(data = c(successA, failA, successB, failB), nrow = 2, byrow = TRUE))$p.value,
                       error = function(e){return(1)})

      #if first significant result, record stopping time
      if (pVal <= alpha) {
        stoppingTimes[m] <- i
        rejections[m] <- 1
        break()
      }
      #if we have reached last iteration, we have not stopped; register max stopping time
      if (i == maxSimStoptime) {
        stoppingTimes[m] <- maxSimStoptime
        rejections[m] <- 0
      }
    }
  }
  return(list(stoppingTimes = stoppingTimes, rejections = rejections))
}

#' Plot bounds of a safe confidence sequence of the difference or log odds ratio for two proportions
#' against the number of data blocks in two data streams ya and yb.
#'
#' @param ya positive observations/ events per data block in group a: a numeric with integer values
#' between (and including) 0 and \code{na}, the number of observations in group a per block.
#' @param yb positive observations/ events per data block in group b: a numeric with integer values
#' between (and including) 0 and \code{nb}, the number of observations in group b per block.
#' @param safeDesign a safe test design for two proportions retrieved through \code{\link{designSafeTwoProportions}()}.
#' @param differenceMeasure the difference measure to construct the confidence interval for:
#' one of "difference" and "odds".
#' @param precision precision of the grid to search over for the confidence sequence bounds.
#' @param deltaStart for the odds difference measure: the (absolute value of the) smallest
#' log odds ratio to assess for in- or exclusion in the confidence sequence. Default 0.001.
#' @param deltaStop for the odds difference measure: the (absolute value of the) highest
#' log odds ratio to assess for in- or exclusion in the confidence sequence. Default 3.
#' @param trueDifference true difference or log odds ratio in groups A and B: added to the plot.
#'
#' @return no return value; called for its side effects, a plot of the confidence sequence.
#' @export
#' @importFrom purrr %>%
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(39413)
#' ya <- rbinom(n = 30, size = 1, prob = 0.1)
#' yb <- rbinom(n = 30, size = 1, prob = 0.8)
#' balancedSafeDesign <- designSafeTwoProportions(na = 1,
#'                                                nb = 1,
#'                                                nBlocksPlan = 30)
#' plotConfidenceSequenceTwoProportions(ya = ya,
#'                                      yb = yb,
#'                                      safeDesign = balancedSafeDesign,
#'                                      differenceMeasure = "difference",
#'                                      precision = 25,
#'                                      trueDifference = 0.7)
#'
#' #log odds ratio difference measure
#' plotConfidenceSequenceTwoProportions(ya = ya,
#'                                      yb = yb,
#'                                      safeDesign = balancedSafeDesign,
#'                                      differenceMeasure = "odds",
#'                                      precision = 25,
#'                                      deltaStop = 5,
#'                                      trueDifference = log(36))
#'
#' #switch ya and yb: observe negative log odds ratio in the data, plot mirrored in x-axis
#' plotConfidenceSequenceTwoProportions(ya = yb,
#'                                      yb = ya,
#'                                      safeDesign = balancedSafeDesign,
#'                                      differenceMeasure = "odds",
#'                                      precision = 25,
#'                                      deltaStop = 5,
#'                                      trueDifference = -log(36))
#'
plotConfidenceSequenceTwoProportions <- function(ya, yb,
                                                 safeDesign,
                                                 differenceMeasure = c("difference", "odds"),
                                                 precision = 100,
                                                 deltaStart = 0.001,
                                                 deltaStop = 3,
                                                 trueDifference = NA){
  differenceMeasure <- match.arg(differenceMeasure)

  if (differenceMeasure == "difference") {
    lowerBounds <- upperBounds <- numeric(length(ya))
    for (m in seq_along(ya)) {
      confidenceInterval <- getConfidenceBoundsForDifferenceTwoProportions(
        ya = ya[1:m],
        yb = yb[1:m],
        precision = precision,
        safeDesign = safeDesign
      )
      lowerBounds[m] <- confidenceInterval[["lowerBound"]]
      upperBounds[m] <- confidenceInterval[["upperBound"]]
    }
    graphics::plot(x = 1:length(ya), y = lowerBounds, ylim = c(-1, 1), col = "blue", type = "l",
                   xlab = "data block number", ylab = "difference",
                   main = "Upper and lower bound confidence sequence for difference")
    graphics::lines(x = 1:length(yb), y = upperBounds, col = "red")
    graphics::abline(h = 0, lty = 2, col = "grey")
    if (!is.na(trueDifference)) {
      graphics::abline(h = trueDifference, lty = 4, col = "black")
    }
  } else {
    positiveLOREstimate <- (sum(yb)/length(yb) - sum(ya)/length(ya)) > 0
    #delta in [0, -infty], we are estimating an upper bound
    if (!positiveLOREstimate) {
      deltaStart <- -1 * deltaStart
      deltaStop <- -1 * deltaStop
    }
    ciValues <- calculateEValuesForOddsDeltaGrid(ya = ya, yb = yb,
                                                 na = safeDesign[["nPlan"]][["na"]],
                                                 nb = safeDesign[["nPlan"]][["nb"]],
                                                 lowerBound = positiveLOREstimate,
                                                 priorParameters = safeDesign[["betaPriorParameterValues"]],
                                                 precision = precision,
                                                 deltaStart = deltaStart,
                                                 deltaStop = deltaStop,
                                                 alpha = safeDesign[["alpha"]])
    if (positiveLOREstimate) {
      plotdfstep <- ciValues %>%
        dplyr::filter(.data[["E"]] < 1/safeDesign[["alpha"]]) %>%
        dplyr::group_by(.data[["block"]]) %>%
        dplyr::summarise(delta = min(.data[["delta"]]))
      #make sure the step function walks until the end of the x axis
      plotdfstepextra <- rbind(plotdfstep,
                               data.frame(
                                 block = length(ya),
                                 delta = max(plotdfstep$delta)
                               )
      )
    } else {
      plotdfstep <- ciValues %>%
        dplyr::filter(.data[["E"]] < 1/safeDesign[["alpha"]]) %>%
        dplyr::group_by(.data[["block"]]) %>%
        dplyr::summarise(delta = max(.data[["delta"]]))
      #make sure the step function walks until the end of the x axis
      plotdfstepextra <- rbind(plotdfstep,
                               data.frame(
                                 block = length(ya),
                                 delta = min(plotdfstep$delta)
                               )
      )
    }

    if (positiveLOREstimate) {
      yLimits <- c(0, max(plotdfstepextra$delta, trueDifference) + 0.1)
    } else {
      yLimits <- c(min(plotdfstepextra$delta, trueDifference) - 0.1, 0)
    }
    graphics::plot(x = plotdfstepextra$block, y = plotdfstepextra$delta,
                   ylim = yLimits, col = "blue", type = "l",
                   xlab = "data block number", ylab = "difference",
                   main = "Confidence sequence for log odds ratio")
    if (!is.na(trueDifference)) {
      graphics::abline(h = trueDifference, lty = 2, col = "grey")
    }
  }
}

#' Simulate the coverage of a safe confidence sequence for differences between proportions
#' for a given distribution and safe design.
#'
#' @param successProbabilityA probability of observing a success in group A.
#' @param trueDelta difference in probability between group A and B.
#' @param safeDesign a safe test design for two proportions retrieved through \code{\link{designSafeTwoProportions}()}.
#' @param precision precision of the grid to search over for the confidence sequence bounds. Default 100.
#' @param M number of simulations to carry out. Default 1000.
#' @param numberForSeed number for seed to set, default NA.
#'
#' @return the proportion of simulations where the trueDelta was included in the confidence sequence.
#' @export
#'
#' @examples
#' balancedSafeDesign <- designSafeTwoProportions(na = 1,
#'                                                nb = 1,
#'                                                nBlocksPlan = 20)
#' simulateCoverageDifferenceTwoProportions(successProbabilityA = 0.2,
#'                                          trueDelta = 0,
#'                                          safeDesign = balancedSafeDesign,
#'                                          M = 100,
#'                                          precision = 25,
#'                                          numberForSeed = 1082021)
simulateCoverageDifferenceTwoProportions <- function(successProbabilityA,
                                                        trueDelta,
                                                        safeDesign,
                                                        precision = 100,
                                                        M = 1000,
                                                        numberForSeed = NA){
  if (!is.na(numberForSeed)) {
    set.seed(numberForSeed)
  }

  successProbabilityB <- successProbabilityA + trueDelta
  trueDeltaIncluded <- logical(M)

  for (simulationNumber in 1:M) {
    yaSim <- rbinom(n = safeDesign[["nPlan"]]["nBlocksPlan"], size = 1, prob = successProbabilityA)
    ybSim <- rbinom(n = safeDesign[["nPlan"]]["nBlocksPlan"], size = 1, prob = successProbabilityB)
    confidenceInterval <- getConfidenceBoundsForDifferenceTwoProportions(
      ya = yaSim,
      yb = ybSim,
      precision = precision,
      safeDesign = safeDesign
    )
    trueDeltaIncluded[simulationNumber] <- (trueDelta >= confidenceInterval[["lowerBound"]]) &
      (trueDelta <= confidenceInterval[["upperBound"]])
  }

  return(mean(trueDeltaIncluded))
}

#NON-EXPORT --------------------------------------------------------------------
#basics ------------------------------------------------------------------------
calculateETwoProportions <- function(na1, na, nb1, nb, thetaA, thetaB, theta0){
  exp(
    na1*log(thetaA) + (na - na1)*log(1-thetaA) + nb1*log(thetaB) + (nb-nb1)*log(1 - thetaB) -
      (na1 + nb1) * log(theta0) - (na + nb - na1 - nb1) * log(1 - theta0)
  )
}

calculateThetaBFromThetaAAndLOR <- function(thetaA, lOR){
  c <- exp(lOR)*thetaA/(1-thetaA)
  return(c/(1+c))
}

logOddsRatio <- function(thetaA, thetaB){
  log(thetaB/(1 - thetaB) * (1 - thetaA)/thetaA)
}

likelihoodTwoProportions<- function(na1, na, nb1, nb, thetaA, thetaB){
  exp(
    na1*log(thetaA) + (na - na1)*log(1-thetaA) + nb1*log(thetaB) + (nb-nb1)*log(1 - thetaB)
  )
}

#No restrictions fncs ---------------------------------------------------------
#NOTE THAT FOR THESE FUNCTIONS TOTALS ARE USED
bernoulliMLTwoProportions <- function(totalSuccess, totalFail, priorSuccess, priorFail){
  (totalSuccess + priorSuccess)/(totalSuccess + totalFail + priorSuccess + priorFail)
}

updateETwoProportions <- function(totalSuccessA, totalFailA, totalSuccessB, totalFailB, na, nb,
                    betaA1, betaA2, betaB1, betaB2){

  thetaA <- bernoulliMLTwoProportions(totalSuccess = totalSuccessA,
                        totalFail = totalFailA,
                        priorSuccess = betaA1,
                        priorFail = betaA2)
  thetaB <- bernoulliMLTwoProportions(totalSuccess = totalSuccessB,
                        totalFail = totalFailB,
                        priorSuccess = betaB1,
                        priorFail = betaB2)

  theta0 <- (na*thetaA + nb*thetaB)/(na+nb)

  return(
    list(
      thetaA = thetaA,
      thetaB = thetaB,
      theta0 = theta0
    )
  )
}

#Restriction on H1 variant fncs ------------------------------------------------
createStartEWithRestrictionTwoProportions <- function(na, nb,
                                        delta,
                                        logOdds,
                                        betaA1,
                                        betaA2,
                                        gridSize = 1e3
                                        ){
  #do not start at 0/ end at 1, because using log/ exp trick on calcualtions
  #later for precision and log(0) raises error
  rhoGrid <- seq(1/gridSize, 1 - 1/gridSize, length.out = gridSize)
  rhoGridDensity <- dbeta(x = rhoGrid, shape1 = betaA1, shape2 = betaA2)
  densityStart <- rhoGridDensity/sum(rhoGridDensity)

  #calculate marginal pred. prob
  if (logOdds) {
    #log odds: theta A in (0,1), no reparameterization needed
    thetaAgrid <- rhoGrid
    thetaBgrid <- sapply(thetaAgrid, calculateThetaBFromThetaAAndLOR, lOR = delta)

    thetaA <- as.numeric(thetaAgrid %*% densityStart)
    thetaB <- calculateThetaBFromThetaAAndLOR(thetaA = thetaA, lOR = delta)
  } else {
    #if delta < 0, add term to reparameterization
    thetaAgrid <- rhoGrid*(1 - abs(delta)) - ifelse(delta < 0, delta, 0)
    thetaBgrid <- thetaAgrid + delta

    thetaA <- as.numeric(thetaAgrid %*% densityStart)
    thetaB <- thetaA + delta
  }

  return(list(posteriorDensity = densityStart,
              thetaAgrid = thetaAgrid,
              thetaBgrid = thetaBgrid,
              thetaA = thetaA,
              thetaB = thetaB,
              theta0 = (na*thetaA + nb*thetaB)/(na+nb)))
}

updateEWithRestrictionTwoProportions <- function(na1, nb1, na, nb, delta, logOdds,
                                   priorDensity, thetaAgrid, thetaBgrid){
  likelihoodTimesPrior <- exp(na1*log(thetaAgrid) + (na - na1)*log(1-thetaAgrid) +
                           nb1*log(thetaBgrid) + (nb - nb1)*log(1-thetaBgrid) +
                           log(priorDensity))

  #normalize
  posteriorDensity <- likelihoodTimesPrior/sum(likelihoodTimesPrior)

  #calculate new marginal pred. probs
  thetaA <- as.numeric(thetaAgrid %*% posteriorDensity)
  if (logOdds) {
    thetaB <- calculateThetaBFromThetaAAndLOR(thetaA = thetaA, lOR = delta)
  } else {
    thetaB <- thetaA + delta
  }

  return(list(posteriorDensity = posteriorDensity,
              thetaAgrid = thetaAgrid,
              thetaBgrid = thetaBgrid,
              thetaA = thetaA,
              thetaB = thetaB,
              theta0 = (na*thetaA + nb*thetaB)/(na+nb)))
}

#Confidence functions ---------------------------------------------------------
calculateKLTwoProportions <- function(candidateThetaA, distanceFunction, delta, na, nb, breveThetaA,breveThetaB){
  candidateThetaB <- distanceFunction(candidateThetaA, delta)

  na1vec <- 0:na
  nb1vec <- 0:nb
  outcomeSpace <- expand.grid(na1vec, nb1vec)

  likelihoodAlternative <- likelihoodTwoProportions(na1 = outcomeSpace[,1], na = na,
                                                    nb1 = outcomeSpace[,2], nb = nb,
                                                    thetaA = breveThetaA, thetaB = breveThetaB
  )
  likelihoodNull <- likelihoodTwoProportions(na1 = outcomeSpace[,1], na = na,
                                             nb1 = outcomeSpace[,2], nb = nb,
                                             thetaA = candidateThetaA, thetaB = candidateThetaB
  )

  sum(likelihoodAlternative * (log(likelihoodAlternative) - log(likelihoodNull)))
}

derivativeKLTwoProportionsLinear <- function(candidateThetaA, delta, na, nb, breveThetaA, breveThetaB, c = 1){
  candidateThetaB <- candidateThetaA + delta

  na*((1 - breveThetaA)/(1 - candidateThetaA) - breveThetaA/candidateThetaA) +
    nb*c*((1 - breveThetaB)/(1 - candidateThetaB) - breveThetaB/candidateThetaB)
}

calculateEValuesForLinearDeltaGrid <- function(ya, yb, na, nb,
                                               priorParameters,
                                               precision = 10,
                                               alpha = 0.05,
                                               runningIntersection = TRUE){
  deltaVec <- seq(-0.99, 0.99, length.out = precision)
  ciEValues <- data.frame()

  for (delta in deltaVec) {
    currentE <- 1
    breveThetaA <- bernoulliMLTwoProportions(totalSuccess = 0,
                                             totalFail = 0,
                                             priorSuccess = priorParameters[["betaA1"]],
                                             priorFail = priorParameters[["betaA2"]])

    breveThetaB <- bernoulliMLTwoProportions(totalSuccess = 0,
                                             totalFail = 0,
                                             priorSuccess = priorParameters[["betaB1"]],
                                             priorFail = priorParameters[["betaB2"]])

    thetaARIPr <- tryOrFailWithNA(stats::uniroot(derivativeKLTwoProportionsLinear,
                                          interval = c(max(c(0, -delta)) + 1e-3, min(c(1,1 - delta))-1e-3),
                                          delta = delta,
                                          na = na, nb = nb,
                                          breveThetaA = breveThetaA,
                                          breveThetaB = breveThetaB)$root)

    #loop over all observed data
    for (i in seq_along(ya)) {
      #if RIPr could not be determined, skip this iteration. E-value stays 1
      if (!is.na(thetaARIPr)) {
        likelihoodAlternative <- likelihoodTwoProportions(na1 = ya[i], na = na,
                                                          nb1 = yb[i], nb = nb,
                                                          thetaA = breveThetaA, thetaB = breveThetaB
        )
        likelihoodRIPr <- likelihoodTwoProportions(na1 = ya[i], na = na,
                                                   nb1 = yb[i], nb = nb,
                                                   thetaA = thetaARIPr, thetaB = thetaARIPr + delta
        )
        currentE <- currentE * likelihoodAlternative/likelihoodRIPr
      }


      #if we reject, we reject this delta FOR EVER in the running intersection
      #do not need to loop over the rest of the data
      if ((currentE >= 1/alpha) & runningIntersection) {
        break()
      }

      #update the E variable for the next data block
      breveThetaA <- bernoulliMLTwoProportions(totalSuccess = sum(ya[1:i]),
                                               totalFail = i*na - sum(ya[1:i]),
                                               priorSuccess = priorParameters[["betaA1"]],
                                               priorFail = priorParameters[["betaA2"]])

      breveThetaB <- bernoulliMLTwoProportions(totalSuccess = sum(yb[1:i]),
                                               totalFail = i*nb - sum(yb[1:i]),
                                               priorSuccess = priorParameters[["betaB1"]],
                                               priorFail = priorParameters[["betaB2"]])

      thetaARIPr <- tryOrFailWithNA(stats::uniroot(derivativeKLTwoProportionsLinear, interval = c(max(c(0, -delta)) + 1e-3, min(c(1,1 - delta))-1e-3),
                                            delta = delta,
                                            na = na, nb = nb,
                                            breveThetaA = breveThetaA,
                                            breveThetaB = breveThetaB)$root)
    }
    ciEValues <- rbind(ciEValues, data.frame(delta = delta, E = currentE))
  }
  return(ciEValues)
}

calculateEValuesForOddsDeltaGrid <- function(ya, yb, na, nb,
                                             lowerBound = TRUE,
                                             priorParameters,
                                             precision = 10,
                                             deltaStart = 0,
                                             deltaStop = 10,
                                             alpha = 0.05,
                                             runningIntersection = TRUE,
                                             stopAfterBoundHasBeenFound = FALSE){
  if (lowerBound & any(c(deltaStart, deltaStop) < 0)) {
    stop("Cannot check for negative bound values when assessing lower bound.")
  }

  if (!lowerBound & any(c(deltaStart, deltaStop) > 0)) {
    stop("Cannot check for positive bound values when assessing upper bound.")
  }

  deltaVector <- seq(deltaStart, deltaStop, length.out = precision)
  ciEValues <- data.frame()

  for (delta in deltaVector) {
    currentE <- 1
    breveThetaA <- bernoulliMLTwoProportions(totalSuccess = 0,
                                             totalFail = 0,
                                             priorSuccess = priorParameters[["betaA1"]],
                                             priorFail = priorParameters[["betaA2"]])

    breveThetaB <- bernoulliMLTwoProportions(totalSuccess = 0,
                                             totalFail = 0,
                                             priorSuccess = priorParameters[["betaB1"]],
                                             priorFail = priorParameters[["betaB2"]])

    #if the point alternative lies within H0(delta), the RIPr and the point alternative should coincide
    if (sign(delta)*logOddsRatio(breveThetaA, breveThetaB) <= sign(delta)*delta) {
      thetaARIPr <- breveThetaA
    } else {
      #otherwise, the RIPr lies on the lOR line
      thetaARIPr <- stats::optim(0.5, fn = calculateKLTwoProportions,
                          method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4,
                          distanceFunction = calculateThetaBFromThetaAAndLOR,
                          delta = delta,
                          na = na, nb = nb,
                          breveThetaA = breveThetaA, breveThetaB = breveThetaB)$par
    }

    for (i in seq_along(ya)) {
      #if point alternative coincides with H0(delta), E value for this block equals 1
      #i.e., the E value remains the same
      if (breveThetaA != thetaARIPr) {
        likelihoodAlternative <- likelihoodTwoProportions(na1 = ya[i], na = na,
                                                          nb1 = yb[i], nb = nb,
                                                          thetaA = breveThetaA, thetaB = breveThetaB
        )

        likelihoodRIPr <- likelihoodTwoProportions(na1 = ya[i], na = na,
                                                   nb1 = yb[i], nb = nb,
                                                   thetaA = thetaARIPr, thetaB = calculateThetaBFromThetaAAndLOR(thetaARIPr, delta)
        )
        currentE <- currentE * likelihoodAlternative/likelihoodRIPr
      }

      ciEValues <- rbind(ciEValues, data.frame(delta = delta, block = i, E = currentE))

      #if we reject, we reject this delta FOR EVER (running intersection)
      if ((currentE >= 1/alpha) & runningIntersection) {
        break()
      }

      #update the E variable
      breveThetaA <- bernoulliMLTwoProportions(totalSuccess = sum(ya[1:i]),
                                               totalFail = i*na - sum(ya[1:i]),
                                               priorSuccess = priorParameters[["betaA1"]],
                                               priorFail = priorParameters[["betaA2"]])

      breveThetaB <- bernoulliMLTwoProportions(totalSuccess = sum(yb[1:i]),
                                               totalFail = i*nb - sum(yb[1:i]),
                                               priorSuccess = priorParameters[["betaA1"]],
                                               priorFail = priorParameters[["betaA2"]])

      #if the point alternative lies within H0(delta), the RIPr and the point alternative should coincide
      if (sign(delta)*logOddsRatio(breveThetaA, breveThetaB) <= sign(delta)*delta) {
        thetaARIPr <- breveThetaA
      } else {
        #otherwise, the RIPr lies on the lOR line
        thetaARIPr <- stats::optim(0.5, fn = calculateKLTwoProportions,
                            method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4,
                            distanceFunction = calculateThetaBFromThetaAAndLOR,
                            delta = delta,
                            na = na, nb = nb,
                            breveThetaA = breveThetaA, breveThetaB = breveThetaB)$par
      }
    }
    #if we want to be fast, we can stop the first time the OR is not significant,
    #if we search in the direction from not extreme -> extreme.
    #Then we have found the bound, where we switch from not include to include
    if (stopAfterBoundHasBeenFound & (currentE <= 1/alpha)) {
        break()
    }
  }
  return(ciEValues)
}

#Main functions ---------------------------------------------------------------
calculateSequential2x2E <- function(aSample, bSample,
                                 restriction = c("none", "difference", "logOddsRatio"),
                                 priorValues,
                                 delta = NULL,
                                 na = 1,
                                 nb = 1,
                                 gridSize = 1e3,
                                 simSetting = FALSE,
                                 alphaSim = 0.05){
  restriction <- match.arg(restriction)

  #these errors should all be caught in the export level function
  #but remain here now for testing purposes
  if (restriction %in% c("difference", "logOddsRatio") & !is.numeric(delta)) {
    stop("Provide numeric value for divergence measure: a difference or log Odds ratio")
  }

  if (length(aSample) != length(bSample)) {
    stop("Can only process complete data blocks: provide vectors with numbers of positive observations per timepoint,",
         "see example in helpfile.")
  }

  if (any(aSample > na | aSample < 0) | any(bSample > nb | bSample < 0)) {
    stop("Provided sample sizes within blocks, na and nb, do not match provided aSample and bSample.")
  }

  #unpack the prior values
  betaA1 <- priorValues[["betaA1"]]
  betaA2 <- priorValues[["betaA2"]]
  betaB1 <- priorValues[["betaB1"]]
  betaB2 <- priorValues[["betaB2"]]

  #set starting E variable
  if (restriction == "difference") {
    eVariable <- createStartEWithRestrictionTwoProportions(na = na, nb = nb,
                                             delta = delta,
                                             logOdds = FALSE,
                                             betaA1 = betaA1, betaA2 = betaA2,
                                             gridSize = gridSize
                                             )

  } else if (restriction == "logOddsRatio") {
    eVariable <- createStartEWithRestrictionTwoProportions(na = na, nb = nb,
                                             delta = delta,
                                             logOdds = TRUE,
                                             betaA1 = betaA1, betaA2 = betaA2,
                                             gridSize = gridSize
                                             )
  } else if (restriction == "none") {
    eVariable <- updateETwoProportions(totalSuccessA = 0, totalFailA = 0,
                         totalSuccessB = 0, totalFailB = 0,
                         na = na, nb = nb,
                         betaA1 = betaA1, betaA2 = betaA2,
                         betaB1 = betaB1, betaB2 = betaB2)

    #for updating without restriction, use totals: store them here
    totalSuccessA <- cumsum(aSample)
    totalSuccessB <- cumsum(bSample)
    groupSizeVecA <- seq_along(totalSuccessA)*na
    groupSizeVecB <- seq_along(totalSuccessB)*nb
    totalFailA <- groupSizeVecA - totalSuccessA
    totalFailB <- groupSizeVecB - totalSuccessB
  }

  currentE <- 1
  for (i in seq_along(aSample)) {
    #use only new data to calculate the new E variable
    newE <- calculateETwoProportions(na1 = aSample[i],
                       na = na,
                       nb1 = bSample[i],
                       nb = nb,
                       thetaA = eVariable[["thetaA"]],
                       thetaB = eVariable[["thetaB"]],
                       theta0 = eVariable[["theta0"]])
    currentE <- newE * currentE

    #in simulation setting, only interested in the stopping time
    if (simSetting & currentE >= (1/alphaSim)) {
      break()
    }

    #after observing the data, also update the E variable
    if (restriction == "none") {
      #updating the E variable without restrictions:
      #using all data seen so far + priorSuccess at the start, new Bernoulli ML
      eVariable <- updateETwoProportions(totalSuccessA = totalSuccessA[i],
                           totalFailA = totalFailA[i],
                           totalSuccessB = totalSuccessB[i],
                           totalFailB = totalFailB[i],
                           na = na, nb = nb,
                           betaA1 = betaA1, betaA2 = betaA2,
                           betaB1 = betaB1, betaB2 = betaB2)
    } else if (restriction == "difference") {
      #updating the E variable with restriction on H1:
      #take product of previous posterior and posterior of NEW data
      eVariable <- updateEWithRestrictionTwoProportions(na1 = aSample[i], nb1 = bSample[i],
                                          na = na, nb = nb,
                                          priorDensity = eVariable[["posteriorDensity"]],
                                          thetaAgrid = eVariable[["thetaAgrid"]],
                                          thetaBgrid = eVariable[["thetaBgrid"]],
                                          delta = delta,
                                          logOdds = FALSE
                                          )
    } else if (restriction == "logOddsRatio") {
      eVariable <- updateEWithRestrictionTwoProportions(na1 = aSample[i], nb1 = bSample[i],
                                          na = na, nb = nb,
                                          priorDensity = eVariable[["posteriorDensity"]],
                                          thetaAgrid = eVariable[["thetaAgrid"]],
                                          thetaBgrid = eVariable[["thetaBgrid"]],
                                          delta = delta,
                                          logOdds = TRUE
                                          )
    }

  }

  if (!simSetting) {
    #we have looped over the entire stream: return the E value
    return(currentE)
  } else {
    return(list(stopTime = i, stopE = currentE))
  }

}

simulateWorstCaseQuantileTwoProportions <- function(na, nb, priorValues,
                                      alternativeRestriction = c("none", "difference", "logOddsRatio"),
                                      alpha,
                                      delta, beta = 0, M = 1e3,
                                      deltaDesign = NULL,
                                      maxSimStoptime = 1e4,
                                      gridSize = 8,
                                      expectedStopTime = FALSE){

  restriction <- match.arg(alternativeRestriction)

  rhoGrid <- seq(1/gridSize, 1 - 1/gridSize, length.out = gridSize)
  if (restriction == "logOddsRatio") {
    #log odds: theta A in (0,1), no reparameterization needed
    thetaAVec <- rhoGrid
  } else {
    #if delta < 0, reparameterize + translate
    thetaAVec <- rhoGrid*(1 - abs(delta)) - ifelse(delta < 0, delta, 0)
  }
  currentWorstCaseQuantile <- 0
  currentWorstCasePower <- currentImpliedTarget <- 1

  message(paste("Simulating E values and stopping times for divergence between groups of ", delta))
  pbSafe <- utils::txtProgressBar(style=1)
  for (t in seq_along(thetaAVec)) {
    stoppingTimes <- stopEs <- numeric(M)

    thetaA <- thetaAVec[t]
    if (restriction == "logOddsRatio") {
      thetaB <- calculateThetaBFromThetaAAndLOR(thetaA, delta)
    } else {
      thetaB <- thetaA + delta
    }

    for (i in 1:M) {
      #For every m, draw a sample of max streamlength and record the time
      #at which we would have stopped
      ya <- rbinom(n = maxSimStoptime, size = na, prob = thetaA)
      yb <- rbinom(n = maxSimStoptime, size = nb, prob = thetaB)
      simResult <- calculateSequential2x2E(aSample = ya, bSample = yb,
                                                  priorValues = priorValues,
                                                  restriction = restriction,
                                                  #if explicitly passsed deltaDesign (neq delta), use that one for test
                                                  #e.g. when studying effect of overestimated/ underestimated effect size
                                                  delta = ifelse(is.null(deltaDesign), delta, deltaDesign),
                                                  na = na,
                                                  nb = nb,
                                                  simSetting = TRUE,
                                                  alphaSim = alpha)
      stoppingTimes[i] <- simResult [["stopTime"]]
      stopEs[i] <- simResult[["stopE"]]
      utils::setTxtProgressBar(pbSafe, value=((t-1)*M+i)/(length(thetaAVec)*M))
    }

    #get the quantile for (1-b) power
    if (expectedStopTime) {
      currentQuantile <- mean(stoppingTimes)
    } else {
      currentQuantile <- quantile(stoppingTimes, probs = 1 - beta)
    }
    currentPower <- mean(stopEs >= 1/alpha)

    #we look for the worst case (1-beta)% stopping time or power: store only that one
    #also store the implied target belonging to the worst case power: exp(Expected [log E-waarde])
    if (currentQuantile >= currentWorstCaseQuantile) {
      currentWorstCaseQuantile <- currentQuantile
    }
    if (currentPower <= currentWorstCasePower) {
      currentWorstCasePower <- currentPower
      currentImpliedTarget <- exp(mean(log(stopEs)))
    }
  }
  close(pbSafe)

  return(list(worstCasePower = currentWorstCasePower,
              worstCaseQuantile = currentWorstCaseQuantile,
              impliedTarget = currentImpliedTarget))
}

simulateWorstCaseDeltaTwoProportions <- function(na, nb, priorValues,
                                   alternativeRestriction = c("none", "difference", "logOddsRatio"),
                                   alpha,
                                   beta, maxSimStoptime,
                                   M = 1e3,
                                   deltaGridSize = 10,
                                   deltamax = 0.99, deltamin = 0.01,
                                   thetaAgridSize = 8){
  deltaVec <- seq(deltamax, deltamin, length.out = deltaGridSize)
  for (deltaIndex in seq_along(deltaVec)) {
    if (simulateWorstCaseQuantileTwoProportions(na = na, nb = nb,
                                 priorValues = priorValues,
                                 alternativeRestriction = alternativeRestriction,
                                 alpha = alpha,
                                 delta = deltaVec[deltaIndex], M = M,
                                 maxSimStoptime = maxSimStoptime,
                                 gridSize = thetaAgridSize)$worstCasePower < (1 - beta)) {
      break()
    }
  }
  return(deltaVec[deltaIndex - 1])
}
