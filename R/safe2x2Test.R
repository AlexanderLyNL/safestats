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
#' chosen to optimize evidence collected over subsequent experiments (REGRET). Pass in the following format: \code{(betaA1 = num1, betaA2 = num2, betaB1 = num3, betaB2 = num4)}.
#' @param M number of simulations used to estimate power or nBlocksPlan. Default \code{1000}.
#'
#' @return Returns a safeDesign object that includes:
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
  names(priorValuesForPrint) <- "Beta hyperparameters"

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
    beta <- 1 - simulateWorstCaseQuantileTwoProportions(delta = delta,
                                         na = na, nb = nb,
                                         priorValues = hyperParameterValues,
                                         alternativeRestriction = alternativeRestriction,
                                         maxSimStoptime = nBlocksPlan,
                                         alpha = alpha, beta = 0, M = M)[["worstCasePower"]]
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
#' two proportions in this package, \code{\link{designSafeTwoProportions}}.
#'
#' @param ya positive observations/ events per data block in group a: a numeric with integer values
#' between (and including) 0 and \code{na}, the number of observations in group a per block.
#' @param yb positive observations/ events per data block in group b: a numeric with integer values
#' between (and including) 0 and \code{nb}, the number of observations in group b per block.
#' @param designObj a safe test design for two proportions retrieved through \code{\link{designSafeTwoProportions}}.
#' @param pilot logical that can be set to true when performing an exploratory analysis
#' without a \code{designObj}; only allows for \code{na = nb = 1}.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{n}{The realised sample size(s).}
#'   \item{eValue}{the e-value of the safe test.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeDesign" described in \code{\link{designSafeTwoProportions}}.}
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

  testResult <- list(designObj = designObj,
                     eValue = eValue,
                     dataName = dataName,
                     n = n)
  class(testResult) <- "safeTest"

  return(testResult)
}

#' Alias for \code{\link{safeTwoProportionsTest}}
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
#' @param hyperparameterList list object, its components hyperparameter lists with a format as described in \code{\link{designSafeTwoProportions}}.
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
#'
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
  for(hyperparameterSet in names(hyperparameterList)){
    message(paste("Retrieving worst case stopping times for hyperparameter set ", hyperparameterSet))
    for(deltaGenerating in deltaVec){
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
  for(i in 1:nrow(resultDataFrame)){
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
#' @param x a result object obtained through \code{\link{simulateTwoProportions}}.
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
#'
#' print(simResult)
#'
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
  for(i in 1:length(x[["hyperparameters"]])){
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
#' @param x a result object obtained through \code{\link{simulateTwoProportions}}.
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
  for(hyperparameters in unique(x[["simdata"]][,"hyperparameters"])){
    plotData <- x[["simdata"]][x[["simdata"]]$hyperparameters == hyperparameters,]
    linecolor <- priorcolors[hyperparameters]
    graphics::lines(x = plotData[,"delta"], y = plotData[,"worstCaseQuantile"], col = linecolor, lty = 2, lwd = 2)
  }

  #then, add the expected stopping times
  for(hyperparameters in unique(x[["simdata"]][,"hyperparameters"])){
    plotData <- x[["simdata"]][x[["simdata"]]$hyperparameters == hyperparameters,]
    linecolor <- priorcolors[hyperparameters]
    graphics::lines(x = plotData[,"delta"], y = plotData[,"expected"], col = linecolor, lty = 1, lwd = 2)
  }

  graphics::legend(x = "topright", legend = c(names(priorcolors), "worst case", "expected"),
                   col = c(priorcolors, "grey", "grey"),
                   lty = c(rep(2, length(priorcolors)), 2, 1), lwd = 2)
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
  for(i in seq_along(aSample)){
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
      break
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
  CurrentWorstCaseQuantile <- 0
  CurrentWorstCasePower <- 1

  message(paste("Simulating E values and stopping times for divergence between groups of ", delta))
  pbSafe <- utils::txtProgressBar(style=1)
  for(t in seq_along(thetaAVec)){
    stoppingTimes <- stopEs <- numeric(M)

    thetaA <- thetaAVec[t]
    if (restriction == "logOddsRatio") {
      thetaB <- calculateThetaBFromThetaAAndLOR(thetaA, delta)
    } else {
      thetaB <- thetaA + delta
    }

    for(i in 1:M){
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
      CurrentQuantile <- mean(stoppingTimes)
    } else {
      CurrentQuantile <- quantile(stoppingTimes, probs = 1 - beta)
    }
    CurrentPower <- mean(stopEs >= 1/alpha)

    #we look for the worst case (1-beta)% stopping time or power: store only that one
    if (CurrentQuantile >= CurrentWorstCaseQuantile) {
      CurrentWorstCaseQuantile <- CurrentQuantile
    }
    if (CurrentPower <= CurrentWorstCasePower) {
      CurrentWorstCasePower <- CurrentPower
    }
  }
  close(pbSafe)

  return(list(worstCasePower = CurrentWorstCasePower,
              worstCaseQuantile = CurrentWorstCaseQuantile))
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
  for(deltaIndex in seq_along(deltaVec)){
    if (simulateWorstCaseQuantileTwoProportions(na = na, nb = nb,
                                 priorValues = priorValues,
                                 alternativeRestriction = alternativeRestriction,
                                 alpha = alpha,
                                 delta = deltaVec[deltaIndex], M = M,
                                 maxSimStoptime = maxSimStoptime,
                                 gridSize = thetaAgridSize)$worstCasePower < (1 - beta)) {
      break
    }
  }
  return(deltaVec[deltaIndex - 1])
}
