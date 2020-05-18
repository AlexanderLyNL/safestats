# 1. Design t-test functions -------

#' Design a Frequentist T-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I and type II error for the frequentist t-test.
#'
#' @inheritParams designSafeT
#'
#' @return Returns an object of class "freqTDesign". An object of class "freqTDesign" is a list containing at least the
#' following components:
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{esMin}{the minimal clinically relevant standardised effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{lowN}{the smallest n of the search space for n provided by the user.}
#'   \item{highN}{the largest n of the search space for n provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#' }
#' @export
#'
#' @examples
#' designFreqT(0.5)
designFreqT <- function(deltaMin, alpha=0.05, beta=0.2, alternative=c("two.sided", "greater", "less"),
                        lowN=3L, highN=100L, testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, ...) {

  stopifnot(lowN >= 2, highN > lowN, alpha > 0, beta >0)

  testType <- match.arg(testType)
  alternative <- match.arg(alternative)

  result <- list(nPlan=NULL, "esMin"=deltaMin, "alpha"=alpha, "beta"=beta,
                 "lowN"=lowN, "highN"=highN, "testType"=testType, "alternative"=alternative)
  class(result) <- "freqTDesign"

  if (deltaMin < 0 && alternative=="greater")
    warning("deltaMin < 0, but in the calculations abs(deltaMin) is used instead.")

  # TODO(Alexander): Also need a warning for deltaMin > 0 and alternative=="less" ?

  deltaMin <- abs(deltaMin)

  n1Plan <- NULL
  n2Plan <- NULL

  if (alternative=="two.sided") {
    threshold <- 1-alpha/2
  } else if (alternative %in% c("greater", "less")) {
    threshold <- 1-alpha
  }

  for (n in seq.int(lowN, highN)) {
    if (testType=="twoSample") {
      someDf <- (1+ratio)*n-2
      someNcp <- sqrt(ratio/(1+ratio)*n)*deltaMin
    } else {
      someDf <- n-1
      someNcp <- sqrt(n)*deltaMin
    }

    powerT <- stats::pt(stats::qt(threshold, df=someDf, ncp=0),
                        df=someDf, ncp=someNcp, lower.tail=FALSE)

    if (powerT >= (1-beta)) {
      n1Plan <- n

      if (testType=="twoSample")
        n2Plan <- ceiling(ratio*n)

      if (testType=="paired")
        n2Plan <- n

      break()
    }
  }

  if (is.null(n1Plan))
    return(result)

  if (is.null(n2Plan)) {
    nPlan <- n1Plan
    names(nPlan) <- "n1Plan"
  } else {
    nPlan <- c(n1Plan, n2Plan)
    names(nPlan) <- c("n1Plan", "n2Plan")
  }

  result[["nPlan"]] <- nPlan
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
#' designed test has to adhere to. Note that it also defines the rejection rule S10 > 1/alpha.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both
#' the sample sizes and deltaS, which defines the test. Note that 1-beta defines the power.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less".
#' @param h0 a number indicating the hypothesised true value of the mean under the null. For the moment h0=0.
#' @param nPlan vector of max length 2 representing the planned sample sizes.
#' @param lowN integer that defines the smallest n of our search space for n.
#' @param highN integer that defines the largest n of our search space for n. This might be the largest n that we
#' are able to fund.
#' @param lowParam numeric that defines the smallest delta of our search space for the test-defining deltaS.
#' @param highParam numeric that defines the largest delta of our search space for the test-defining deltaS.
#' @param tol a number that defines the stepsizes between the lowParam and highParam.
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
#' is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param logging logical, if \code{TRUE} return altSThreshes.
#' @param ... further arguments to be passed to or from methods, but mainly to perform do.calls.
#'
#' @return Returns an object of class "safeDesign". An object of class "safeDesign" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{parameter}{the safe test defining parameter. Here deltaS.}
#'   \item{esMin}{the minimal clinically relevant standardised effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on
#'   whether it was a one-sample or a two-sample test.}
#'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
#'   ratio between the planned randomisation of condition 2 over condition 1.}
#'   \item{lowN}{the smallest n of the search space for n provided by the user.}
#'   \item{highN}{the largest n of the search space for n provided by the user.}
#'   \item{lowParam}{the smallest delta of the search space for delta provided by the user.}
#'   \item{highParam}{the largest delta of the search space for delta provided by the user.}
#'   \item{tol}{the step size between lowParam and highParam provided by the user.}
#'   \item{pilot}{\code{FALSE} (default) specified by the user to indicate that the design is not a pilot study.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.8, alpha=0.08, beta=0.01, alternative="greater")
#' designObj
designSafeT <- function(deltaMin=NULL, alpha=0.05, beta=0.2, alternative=c("two.sided", "greater", "less"),
                        h0=0, nPlan=NULL, lowN=3L, highN=100L, lowParam=0.01, highParam=1.5*abs(deltaMin),
                        tol=0.01, testType=c("oneSample", "paired", "twoSample"),
                        ratio=1, logging=FALSE, ...) {
  stopifnot(alpha > 0, alpha < 1)

  if (is.null(deltaMin) && is.null(nPlan))
    stop("Can't design without (1) beta and deltaMin, or (2) nPlan.")

  if (!is.null(deltaMin)) {
    stopifnot(beta > 0, beta < 1)

    if (!is.null(nPlan)) {
      warning("Both nPlan and deltaMin combined with beta provided. Preference is given to designing with beta; ",
              "nPlan is ignored")
      nPlan <- NULL
    }
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)

  paired <- if (testType=="paired") TRUE else FALSE

  names(h0) <- "mu"

  if (!is.null(nPlan)) {
    return(designPilotSafeT("nPlan"=nPlan, "alpha"=alpha, "h0"=h0, "alternative"=alternative,
                            "lowParam"=lowParam, "tol"=tol, "logging"=logging,
                            "paired"=paired))
  }

  names(deltaMin) <- switch(alternative,
                            "two.sided"="standardised effect sizes at least abs(delta)",
                            "less"="standardised effect sizes smaller than delta",
                            "greater"="standardised effect sizes larger than delta")

  result <- list("nPlan"=NULL, "parameter"=NULL, "esMin"=deltaMin, "alpha"=alpha, "beta"=beta,
                 "alternative"=alternative, "testType"=testType, "paired"=paired, "h0"=h0,
                 "ratio"=ratio, "lowN"=lowN, "highN"=highN, "lowParam"=lowParam,
                 "highParam"=highParam, "tol"=tol, "pilot"=FALSE, "call"=sys.call())
  class(result) <- "safeDesign"

  n1Plan <- NULL
  n2Plan <- NULL

  deltaMin <- abs(deltaMin)

  sCutOff <- 1/alpha

  nDefinitions <- defineTTestN("lowN"=lowN, "highN"=highN, "ratio"=ratio, "testType"=testType)

  n1 <- nDefinitions[["n1"]]
  n2 <- nDefinitions[["n2"]]
  candidateNEff <- nDefinitions[["candidateNEff"]]
  candidateNu <- nDefinitions[["candidateNu"]]

  if (alternative=="two.sided")
    candidateFNcp <- candidateNEff*deltaMin^2
  else
    candidateTNcp <- sqrt(candidateNEff)*deltaMin

  # TODO(Alexander): Should highParam just be deltaMin.
  #   Perhaps show that deltaS < deltaMin for alpha, beta. Use monotonicity
  #
  candidateDeltas <- seq(from=lowParam, to=highParam, by=tol)

  for (i in seq_along(candidateNEff)) {
    if (alternative=="two.sided")
      deltaMinThresh <- sqrt(stats::qf("p"=beta, "df1"=1, "df2"=candidateNu[i], "ncp"=candidateFNcp[i])) #*deltaMin^2)
    else
      deltaMinThresh <- stats::qt("p"=beta, "df"=candidateNu[i], "ncp"=candidateTNcp[i])

    # TODO(Alexander): Under the assumption that this is unimodal, then stop once the value goes down
    if (testType=="twoSample") {
      altSThreshes <- purrr::map_dbl(".x"=candidateDeltas, ".f"=safeTTestStat, "t"=deltaMinThresh,
                                     "n1"=n1[i], "n2"=n2[i], "alternative"=alternative, "paired"=paired)
    } else if (testType %in% c("oneSample", "paired")) {
      altSThreshes <- purrr::map_dbl(".x"=candidateDeltas, ".f"=safeTTestStat, "t"=deltaMinThresh,
                                     "n1"=candidateNEff[i], "n2"=n2[i], "alternative"=alternative, "paired"=paired)
    }

    if (max(altSThreshes) >= sCutOff) {
      nEff <- candidateNEff[i]
      if (testType=="twoSample") {
        n1Plan <- ceiling(n1[i])
        n2Plan <- ceiling(n2[i])
        result[["nEffPlan"]] <- nEff
      } else if (testType %in% c("oneSample", "paired")) {
        n1Plan <- nEff

        if (testType=="paired")
          n2Plan <- nEff

      }

      deltaIndex <- which(altSThreshes >= sCutOff)[1]
      deltaS <- candidateDeltas[deltaIndex]

      if (alternative=="less")
        deltaS <- -deltaS

      result[["parameter"]] <- deltaS
      names(result[["parameter"]]) <- "deltaS"

      if (isTRUE(logging))
        result[["altSThreshes"]] <- altSThreshes

      break()
    }
  }

  if (is.null(n1Plan) || is.null(result[["parameter"]])) {
    warning("Increase deltaMin, or increase highN, currently: ", highN,
            ". Try the function plotSafeTDesignSampleSizeProfile to find  minimal",
            "sample size for deltaMin.")

    result[["lowN"]] <- highN + 1
    result[["highN"]] <- 2*highN
    return(result)
  }

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  return(result)
}


#' Simulate Early Stopping Experiments for the T Test
#'
#' Applied to a safeDesign object this function empirically shows the performance of a safe experiments under
#' optional stopping.
#'
#' @param object A safeDesign obtained obtained from \code{\link{designSafeT}}.
#' @param nsim numeric, number of iterations.
#' @param seed numeric, seed number.
#' @param deltaTrue numeric, if NULL, then the minimally clinically relevant standardised effect size is used
#' as the true data generating effect size deltaTrue.
#' @inherit replicateTTests
#'
#' @import stats
#' @export
#'
#' @examples
#'# Design safe test
#' alpha <- 0.05
#' beta <- 0.20
#' deltaMin <- 1
#' designObj <- designSafeT(deltaMin, alpha=alpha, beta=beta)
#'
#' # Design frequentist test
#' freqObj <- designFreqT(deltaMin, alpha=alpha, beta=beta)
#'
#' # Simulate based on deltaTrue=deltaMin
#' simResultsDeltaTrueIsDeltaMin <- simulate(object=designObj, nsim=100)
#'
#' # Simulate based on deltaTrue > deltaMin
#' simResultsDeltaTrueIsLargerThanDeltaMin <- simulate(
#'   object=designObj, nsim=100, deltaTrue=2)
#'
#' # Simulate under the null deltaTrue = 0
#' simResultsDeltaTrueIsNull <- simulate(
#'   object=designObj, nsim=100, deltaTrue=0)
#'
#' simulate(object=designObj, deltraTrue=0, nsim=100, freqOptioStop=TRUE,
#'          nPlanFreq=freqObj$nPlan)
#'
simulate.safeDesign <- function(object, nsim=1, seed=NULL, deltaTrue=NULL, muGlobal=0, sigmaTrue=1, lowN=3,
                                safeOptioStop=TRUE, freqOptioStop=FALSE, nPlanFreq=NULL,
                                logging=TRUE, pb=TRUE, ...) {

  if (object[["pilot"]])
    stop("No simulation for unplanned pilot designs")

  if (object[["testType"]] %in% c("oneSample", "paired", "twoSample")) {
    if (is.null(deltaTrue))
      deltaTrue <- object[["parameter"]]

    paired <- if (object[["testType"]]=="paired") TRUE else FALSE

    result <- replicateTTests("nPlan"=object[["nPlan"]], "deltaTrue"=deltaTrue,
                              "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired,
                              "alternative"=object[["alternative"]], "lowN"=lowN, "nsim"=nsim,
                              "alpha"=object[["alpha"]], "beta"=object[["beta"]],
                              "safeOptioStop"=safeOptioStop, "parameter"=object[["parameter"]],
                              "freqOptioStop"=freqOptioStop, "nPlanFreq"=nPlanFreq,
                              "logging"=logging, "seed"=seed, "pb"=pb, ...)

    object <- utils::modifyList(object, result)
    class(object) <- "safeTSim"
    return(object)
  }
}


#' Simulate Early Stopping Experiments
#'
#' Simulate multiple data sets to show the effects of optional testing for safe (and frequentist) tests.
#'
#' @inheritParams designSafeT
#' @param deltaTrue numeric, the value of the true standardised effect size (test-relevant parameter).
#' @param muGlobal numeric, the true global mean of a paired or two-sample t-test. Its value should not
#' matter for the test. This parameter is treated as a nuisance.
#' @param sigmaTrue numeric > 0,the true standard deviation of the data. Its value should not  matter
#' for the test.This parameter treated is treated as a nuisance.
#' @param lowN the smallest number of samples (first group) at which monitoring of the tests begins.
#' @param nsim the number of replications, that is, experiments with max samples nPlan.
#' @param safeOptioStop logical, \code{TRUE} implies that optional stopping simulation is performed for
#' the safe test.
#' @param parameter numeric, the safe test defining parameter, i.e., deltaS (use designSafeT to find this).
#' @param freqOptioStop logical, \code{TRUE} implies that optional stopping simulation is performed for
#' the frequentist test.
#' @param nPlanFreq the frequentist sample size(s) to plan for. Acquired from \code{\link{designFreqT}}.
#' @param paired logical, if \code{TRUE} then paired t-test.
#' @param seed To set the seed for the simulated data.
#' @param logging logical, if \code{TRUE}, then return the simulated data.
#' @param pb logical, if \code{TRUE}, then show progress bar.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class "safeTSim". An object of class "safeTSim" is a list containing at
#' least the following components:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{deltaTrue}{the value of the true standardised effect size (test-relevant parameter) provided by
#'   the user.}
#'   \item{muGlobal}{the true global mean of a paired or two-sample t-test (nuisance parameter) provided by
#'   the user.}
#'   \item{paired}{if \code{TRUE} then paired t-test.}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user.}
#'   \item{lowN}{the smallest number of samples (first group) at which monitoring of the tests begins.}
#'   \item{nsim}{the number of replications of the experiment.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{parameter}{the parameter (point prior) used in the safe test derived from the design.
#'   Acquired from \code{\link{designSafeT}}.}
#'   \item{nPlanFreq}{the frequentist planned sample size(s). Acquired from \code{\link{designFreqT}}}
#'   \item{safeSim}{list with the simulation results of the safe test under optional stopping.}
#'   \item{freqSim}{list with the simulation results of the frequentist test under optional stopping.}
#'}
#'
#' @export
#'
#' @examples
#'
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
#' nPlanFreq=freqObj$nPlan, beta=beta, nsim=400)
#'
#' # Should be about 1-beta
#' simResults$safeSim$powerAtN1Plan
#'
#' # This is higher due to optional stopping
#' simResults$safeSim$powerOptioStop
#'
#' # Optional stopping allows us to do better than n1PlanFreq once in a while
#' simResults$safeSim$probLeqN1PlanFreq
#' graphics::hist(simResults$safeSim$allN, main="Histogram of stopping times", xlab="n1",
#' breaks=seq.int(designObj$nPlan[1]))
#'
#' # Simulate under the alternative with deltaTrue > deltaMin
#' simResults <- replicateTTests(nPlan=designObj$nPlan, deltaTrue=1.5, parameter=designObj$parameter,
#' nPlanFreq=freqObj$nPlan, beta=beta, nsim=400)
#'
#' # Should be larger than 1-beta
#' simResults$safeSim$powerAtN1Plan
#'
#' # This is even higher due to optional stopping
#' simResults$safeSim$powerOptioStop
#'
#' # Optional stopping allows us to do better than n1PlanFreq once in a while
#' simResults$safeSim$probLeqN1PlanFreq
#' graphics::hist(simResults$safeSim$allN, main="Histogram of stopping times", xlab="n1",
#' breaks=seq.int(designObj$nPlan[1]))
#'
#' # Under the null deltaTrue=0
#' simResults <- replicateTTests(nPlan=designObj$nPlan, deltaTrue=0, parameter=designObj$parameter,
#' nPlanFreq=freqObj$nPlan, freqOptioStop=TRUE, beta=beta, nsim=400)
#'
#'# Should be lower than alpha, because if the null is true, P(S > 1/alpha) < alpha for all n
#' simResults$safeSim$powerAtN1Plan
#'
#' # This is a bit higher due to optional stopping, but if the null is true,
#' # then still P(S > 1/alpha) < alpha for all n
#' simResults$safeSim$powerOptioStop
#'
#' # Should be lowr than alpha, as the experiment is performed as was planned
#' simResults$freqSim$powerAtN1Plan
#'
#' # This is larger than alpha, due to optional stopping.
#' simResults$freqSim$powerOptioStop
#' simResults$freqSim$powerOptioStop > alpha
replicateTTests <- function(nPlan, deltaTrue, muGlobal=0, sigmaTrue=1, paired=FALSE,
                            alternative=c("two.sided", "greater", "less"), lowN=3,
                            nsim=1000L, alpha=0.05, beta=0.2,
                            safeOptioStop=TRUE, parameter=NULL,
                            freqOptioStop=FALSE, nPlanFreq=NULL,
                            logging=TRUE, seed=NULL, pb=TRUE, ...) {

  stopifnot(all(nPlan > lowN), lowN > 0, nsim > 0, alpha > 0, alpha < 1,
            any(safeOptioStop, freqOptioStop))

  alternative <- match.arg(alternative)

  n1Plan <- nPlan[1]
  n2Plan <- NULL

  if (length(nPlan)==2) {
    n2Plan <- nPlan[2]

    testType <- if (paired) "paired" else "twoSample"
  } else {
    if (paired) {
      warning("Paired simulations wanted, but length(nPlan)=1. Duplicate n2Plan=n1Plan")
      n2Plan <- n1Plan
      nPlan <- c(n1Plan, n2Plan)
    }

    testType <- if (paired) "paired" else "oneSample"
  }

  result <- list("nPlan"=nPlan, "deltaTrue"=deltaTrue, "muGlobal"=muGlobal, "paired"=paired,
                 "alternative"=alternative, "lowN"=lowN, "nsim"=nsim, "alpha"=alpha, "beta"=beta, "testType"=testType,
                 "parameter"=parameter, "nPlanFreq"=nPlanFreq, safeSim=list(), freqSim=list())
  class(result) <- "safeTSim"

  if (safeOptioStop) {
    if (is.null(parameter)) {
      stop("To simulate safe t-tests results under optional stopping, this function 'replicateTTests' requires ",
           "the specification of the safe test with a parameter. This parameter can be acquired by running ",
           "the 'designSafeT()' function")
    }

    if (paired && n1Plan != n2Plan)
      stop("For a paired t-test n2Plan needs to equal n1Plan")

    safeSim <- list("powerOptioStop"=NA, "powerAtN1Plan"=NA, "nMean"=NA, "probLeqN1PlanFreq"=NA,
                    "probLessNDesign"=NA, "lowN"=NA)

    allSafeN <- rep(n1Plan, "times"=nsim)
    sValues <- safeDecisionAtN <- allSafeDecisions <- vector("mode"="integer", "length"=nsim)
  }

  if (freqOptioStop) {
    if (!safeOptioStop) {
      if (is.null(nPlanFreq)) {
        warning("No nPlanFreq specified, use nPlan instead.")
        n1PlanFreq <- n1Plan
        n2PlanFreq <- n2Plan
      }
    } else {
      n1PlanFreq <- nPlanFreq[1]
      n2PlanFreq <- NULL

      if (length(nPlanFreq)==2) {
        n2PlanFreq <- nPlanFreq[2]
      } else {
        if (paired) {
          warning("Paired simulations wanted, but length(nPlan)=1. Duplicate n2Plan=n1Plan")
          n2PlanFreq <- n1PlanFreq
          nPlanFreq <- c(n1PlanFreq, n2PlanFreq)
          result[["nPlanFreq"]] <- nPlanFreq
        }
      }
    }

    # Note(Alexander): This means that n1Plan and n2Plan refer to the planned samples of the safe tests

    if (is.null(n1PlanFreq)) {
      stop("To simulate frequentist t-tests results under optional stopping, this ",
           "function 'replicateTTests' requires the specification of n1PlanFreq. To figure out how many ",
           "samples one requires in a frequentist test, please run the 'designFreqT' function.")
    }

    if (!is.null(n2Plan) && is.null(n2PlanFreq)) {
      stop("To simulate a two-sample frequentist t-tests results under optional stopping, this",
           "function 'replicateTTests' requires the specification of n1PlanFreq. To figure out how many ",
           "samples one requires in a frequentist test, please run the 'designFreqT' function.")
    }

    if (paired && n1PlanFreq != n2PlanFreq)
      stop("For a paired t-test n2PlanFreq needs to equal n1PlanFreq")

    freqSim <- list("powerOptioStop"=NA, "powerAtN1Plan"=NA, "nMean"=NA, "probLessNDesign"=NA, "lowN"=NA)

    allFreqN <- rep(n1PlanFreq, "times"=nsim)
    pValues <- freqDecisionAtN <- allFreqDecisions <- vector("mode"="integer", "length"=nsim)
  }

  ratio <- if (is.null(n2Plan) || paired) 1 else n2Plan/n1Plan

  someData <- generateTTestData("nPlan"=c(n1Plan, n2Plan), "nsim"=nsim, "deltaTrue"=deltaTrue,
                                "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed)

  dataGroup1 <- someData[["dataGroup1"]]
  dataGroup2 <- someData[["dataGroup2"]]

  if (safeOptioStop) {
    n1Samples <- seq.int(lowN, n1Plan)
    n2Samples <- if (is.null(n2Plan)) NULL else ceiling(ratio*n1Samples)

    if (pb)
      pbSafe <- utils::txtProgressBar(style=1, title="Safe optional stopping")

    for (iter in seq.int(nsim)) {
      subData1 <- dataGroup1[iter, ]
      subData2 <- dataGroup2[iter, ]

      someT <- unname(stats::t.test("x"=subData1, "y"=subData2, "alternative"=alternative,
                                    "var.equal"=TRUE, "paired"=paired)[["statistic"]])
      someS <- safeTTestStat("t"=someT, "parameter"=parameter, "n1"=n1Plan, "n2"=n2Plan, "alternative"=alternative,
                             "paired"=paired)

      sValues[iter] <- someS

      if (someS >= 1/alpha)
        safeDecisionAtN[iter] <- 1

      for (k in seq_along(n1Samples)) {

        # TODO(Alexander): Perhaps replace by custom t computing to speed things up
        #
        someT <- unname(stats::t.test("x"=subData1[seq.int(n1Samples[k])], "y"=subData2[seq.int(n2Samples[k])],
                                      "alternative"=alternative, "var.equal"=TRUE, "paired"=paired)[["statistic"]])

        someS <- safeTTestStat("n1"=n1Samples[k], "n2"=n2Samples[k], "t"=someT, "parameter"=parameter,
                               "alternative"=alternative, "paired"=paired)

        if (someS >= 1/alpha) {
          allSafeN[iter] <- n1Samples[k]
          allSafeDecisions[iter] <- 1

          sValues[iter] <- someS
          break()
        }
      } # End loop lowN to n1Plan

      if (pb)
        utils::setTxtProgressBar(pbSafe, value=iter/nsim, title="Experiments")

    } # End iterations

    if (pb)
      close(pbSafe)

    safeSim <- list("powerOptioStop"=mean(allSafeDecisions),
                    "powerAtN1Plan"=mean(safeDecisionAtN),
                    "nMean"=mean(allSafeN),
                    "probLessNDesign"=mean(allSafeN < n1Plan),
                    "lowN"=min(allSafeN), "sValues"=sValues
    )

    safeSim[["allN"]] <- allSafeN
    safeSim[["allSafeDecisions"]] <- allSafeDecisions
    safeSim[["allRejectedN"]] <- allSafeN[-which(allSafeN*allSafeDecisions==0)]

    if (!is.null(nPlanFreq))
      safeSim[["probLeqN1PlanFreq"]] <- mean(allSafeN <= nPlanFreq[1])

    if (isTRUE(logging)) {
      safeSim[["dataGroup1"]] <- dataGroup1
      safeSim[["dataGroup2"]] <- dataGroup2
    }

    result[["safeSim"]] <- safeSim
  }

  if (freqOptioStop) {
    # Note(Alexander): Adjust data set
    #
    if (is.null(n2PlanFreq)) {
      ratio <- 1

      if (n1PlanFreq < n1Plan)
        dataGroup1 <- dataGroup1[, seq.int(n1PlanFreq)]

      if (n1PlanFreq > n1Plan) {
        n1Diff <- n1PlanFreq - n1Plan

        someData <- generateTTestData("nPlan"=c(n1Diff, n2Plan), "nsim"=nsim, "deltaTrue"=deltaTrue,
                                      "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed+1)
        dataGroup1 <- cbind(dataGroup1, someData[["dataGroup1"]])
      }
    } else {
      # Note(Alexander): Two-sample case

      ratio <- if (paired) 1 else n2PlanFreq/n1PlanFreq

      if (n1PlanFreq < n1Plan) {
        dataGroup1 <- dataGroup1[, seq.int(n1PlanFreq)]
      } else if (n1PlanFreq > n1Plan) {
        n1Diff <- n1PlanFreq - n1Plan

        someData <- generateTTestData("nPlan"=c(n1Diff, n2PlanFreq), "nsim"=nsim, "deltaTrue"=deltaTrue,
                                      "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed+1)
        dataGroup1 <- cbind(dataGroup1, someData[["dataGroup1"]])
      }

      if (n2PlanFreq < n2Plan) {
        dataGroup2 <- dataGroup2[, seq.int(n2PlanFreq)]
      } else if (n2PlanFreq > n2Plan) {
        n2Diff <- n2PlanFreq - n2Plan

        someData <- generateTTestData("nPlan"=c(n1PlanFreq, n2Diff), "nsim"=nsim, "deltaTrue"=deltaTrue,
                                      "muGlobal"=muGlobal, "sigmaTrue"=sigmaTrue, "paired"=paired, "seed"=seed+1)
        dataGroup2 <- cbind(dataGroup2, someData[["dataGroup2"]])
      }
    }

    n1Samples <- seq.int(lowN, n1PlanFreq)

    n2Samples <- if (is.null(n2PlanFreq)) NULL else n2Samples <- ceiling(ratio*n1Samples)

    if (pb)
      pbFreq <- utils::txtProgressBar(style=1, title="Frequentist optional stopping")

    for (iter in seq.int(nsim)) {
      subData1 <- dataGroup1[iter, ]
      subData2 <- dataGroup2[iter, ]
      someP <- stats::t.test("x"=subData1, "y"=subData2, "alternative"=alternative,
                             "var.equal"=TRUE, "paired"=paired)[["p.value"]]

      pValues[iter] <- someP

      if (someP < alpha)
        freqDecisionAtN[iter] <- 1

      for (k in seq_along(n1Samples)) {
        someP <- stats::t.test("x"=subData1[seq.int(n1Samples[k])], "y"=subData2[seq.int(n2Samples[k])],
                               "alternative"=alternative, "var.equal"=TRUE, "paired"=paired)[["p.value"]]

        if (someP < alpha) {
          allFreqN[iter] <- n1Samples[k]
          allFreqDecisions[iter] <- 1
          pValues[iter] <- someP
          break()
        }
      } # End loop lowN to n1Plan

      if (pb)
        utils::setTxtProgressBar(pbFreq, value=iter/nsim, title="Experiments")
    } # End iterations

    if (pb)
      close(pbFreq)

    freqSim <- list("powerOptioStop"=mean(allFreqDecisions),
                    "powerAtN1Plan"=mean(freqDecisionAtN),
                    "nMean"=mean(allFreqN),
                    "allFreqDecisions"=allFreqDecisions,
                    "probLessNDesign"=mean(allFreqN < n1PlanFreq),
                    "lowN"=min(allFreqN), "pValues"=pValues
    )

    freqSim[["allN"]] <- allFreqN

    if (safeOptioStop)
      freqSim[["probLeqNSafe"]] <- mean(allFreqN <= n1Plan)

    if (isTRUE(logging)) {
      freqSim[["dataGroup1"]] <- dataGroup1
      freqSim[["dataGroup2"]] <- dataGroup2
    }

    result[["freqSim"]] <- freqSim
  }
  return(result)
}



#' Plots a safeTSim Object
#'
#' @inheritParams plotHistogramDistributionStoppingTimes
#' @param x a "safeDesign" object acquired from \code{\link{designSafeT}}.
#' @param y \code{NULL}.
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
#' simResults <- simulate(designObj, nsim=100)
#'
#' plot(simResults)
#'
#' plot(simResults, showOnlyNRejected=TRUE)
plot.safeTSim <- function(x, y=NULL, showOnlyNRejected=FALSE, nBin=25, ...) {
  plotHistogramDistributionStoppingTimes(x[["safeSim"]],
                                         "nPlan" = x[["nPlan"]][1],
                                         "deltaTrue" = x[["deltaTrue"]],
                                         "showOnlyNRejected" = showOnlyNRejected, "nBin"=nBin)
}

#' Computes a Sequence of (Effective) Sample Sizes
#'
#' Helper function that outputs the sample sizes, effective sample sizes and the degrees of
#' freedom depending on the type of t-test. Also used for z-tests.
#'
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#'
#' @return Returns the sample sizes and degrees of freedom.
#'
#' @examples
#' safestats:::defineTTestN()
defineTTestN <- function(lowN=3, highN=100, ratio=1,
                         testType=c("oneSample", "paired", "twoSample")) {
  testType <- match.arg(testType)

  if (testType %in% c("twoSample")) {
    n1 <- lowN:highN
    n2 <- ceiling(ratio*n1)
    candidateNEff <- ratio/(1+ratio)*n1
    candidateNu <- (1+ratio)*n1-2
  } else if (testType %in% c("oneSample", "paired")) {
    n1 <- lowN:highN
    n2 <- NULL
    candidateNEff <- n1
    candidateNu <- candidateNEff-1
  }
  result <- list("n1"=n1, "n2"=n2, "candidateNEff"=candidateNEff, "candidateNu"=candidateNu)
  return(result)
}



#' Designs a Safe T-Test Based on Planned Samples nPlan
#'
#' Designs a safe experiment for a prespecified tolerable type I error based on planned sample
#' size(s), which are fixed ahead of time. Outputs a list that includes the deltaS, i.e., the
#' safe test defining parameter.
#'
#' @param nPlan the planned sample size(s).
#' @param inverseMethod logical, always \code{TRUE} for the moment.
#' @param paired logical, if \code{TRUE} then paired t-test.
#' @param logging logical, if \code{TRUE}, then add invSToTThresh to output.
#' @param maxIter numeric > 0, the maximum number of iterations of adjustment to the candidate set from
#' lowParam to highParam, if the minimum is not found.
#' @inheritParams replicateTTests
#' @inherit designSafeT
#'
#' @export
#'
#' @examples
#' designPilotSafeT(nPlan=30)
designPilotSafeT <- function(nPlan=50, alpha=0.05, h0=0, alternative=c("two.sided", "greater", "less"),
                             lowParam=0.01, highParam=1.2, tol=0.01, inverseMethod=TRUE,
                             logging=FALSE, paired=FALSE, maxIter=10) {
  # TODO(Alexander): Check relation with forward method, that is, the least conservative test and maximally powered
  # Perhaps trade-off? "inverseMethod" refers to solving minimum of deltaS \mapsto S_{deltaS}^{-1}(1/alpha)
  #
  #
  # TODO(Alexander):
  #       - Add warning when deltaS == lowParam
  #       - Better: Characterise deltaS as a funciton of alpha and n using asymptotic approximation of 1F1
  #
  # TODO(Alexander): Add some bounds on L, or do a presearch
  #
  # stopifnot(n >= 3)
  #
  #     Trick bound the estimation error using Chebyshev: Note do this on log (safeTTestStat)
  #

  alternative <- match.arg(alternative)
  stopifnot(all(nPlan > 0))

  n1 <- nPlan[1]
  n2 <- NULL
  ratio <- 1

  if (length(nPlan)==1) {
    if (isTRUE(paired)) {
      warning("Paired designed specified, but n2 not provided. n2 is set to n1")
      n2 <- n1
      nPlan <- c(n1, n2)
      names(nPlan) <- c("n1Plan", "n2Plan")
      testType <- "paired"
    } else {
      names(nPlan) <- "n1Plan"
      testType <- "oneSample"
    }
  } else if (length(nPlan)==2) {
    n2 <- nPlan[2]
    ratio <- n2/n1

    if (paired) {
      if (n1 != n2) {
        stop("Paired design specified, but nPlan[1] not equal nPlan[2]")
      }
      testType <- "paired"
    } else {
      testType <- "twoSample"
    }
    names(nPlan) <- c("n1Plan", "n2Plan")
  }

  names(h0) <- "mu"

  result <- list("nPlan"=nPlan, "parameter"=NULL, "esMin"=NULL, "alpha"=alpha, "beta"=NULL,
                 "alternative"=alternative, "testType"=testType, "paired"=paired,
                 "h0"=h0, "sigma"=sigma, "kappa"=kappa, "testType"=testType,
                 "ratio"=ratio, "lowParam"=NULL, "highParam"=NULL,
                 "pilot"=FALSE, "call"=sys.call())
  class(result) <- "safeDesign"

  #
  #
  # result <- list("n1Plan"=n1, "n2Plan"=n2, "deltaS"=NA,
  #                "deltaMin"=NULL, "alpha"=alpha, "beta"=NULL,
  #                "lowParam"=lowParam, "highParam"=highParam, "tol"=tol,
  #                "lowN"=NULL, "highN"=NULL, "alternative"=alternative, "testType"=testType,
  #                "ratio"=ratio, "pilot"=TRUE, "call"=sys.call())
  #
  # class(result) <- "safeTDesign"

  if (inverseMethod) {
    invSafeTTestStatAlpha <- function(x) {
      stats::uniroot("f"=safeTTestStatAlpha, "lower"=0, "upper"=3e10, "tol"=1e-9, "parameter"=x,
              "n1"=n1, "n2"=n2, "alpha"=alpha, "alternative"=alternative)
    }

    # Note(Alexander): tryOrFAilWithNA fixes the problem that there's a possibility
    #   that deltaS too small, the function values are then both negative or both positive.
    #
    invSafeTTestStatAlphaRoot <- function(x){tryOrFailWithNA(invSafeTTestStatAlpha(x)$root)}

    candidateDeltas <- seq(lowParam, highParam, by=tol)

    invSToTThresh <- rep(NA, length(candidateDeltas))

    iter <- 1

    while (all(is.na(invSToTThresh)) && iter <= maxIter) {
      invSToTThresh <- purrr::map_dbl(candidateDeltas, invSafeTTestStatAlphaRoot)

      if (all(is.na(invSToTThresh))) {
        candidateDeltas <- seq(highParam, "length.out"=2*length(candidateDeltas), "by"=tol)
        lowParam <- highParam
        highParam <- candidateDeltas[length(candidateDeltas)]
        iter <- iter + 1
      }
    }

    mPIndex <- which(invSToTThresh==min(invSToTThresh, na.rm=TRUE))

    if (iter > 1) {
      result[["lowParam"]] <- lowParam
      result[["highParam"]] <- highParam
    }

    if (mPIndex==length(candidateDeltas)) {
      # Note(Alexander): Check that mPIndex is not the last one.
      errorMsg <- "The test defining deltaS is equal to highParam. Rerun with do.call on the output object"
      lowParam <- highParam
      highParam <- (length(candidateDeltas)-1)*tol+lowParam
      result[["lowParam"]] <- lowParam
      result[["highParam"]] <- highParam
    } else if (mPIndex==1) {
      errorMsg <- "The test defining deltaS is equal to lowParam. Rerun with do.call on the output object"
      highParam <- lowParam
      lowParam <- highParam-(length(candidateDeltas)-1)*tol
      result[["lowParam"]] <- lowParam
      result[["highParam"]] <- highParam
    }

    if (isTRUE(logging))
      result[["invSToTThresh"]] <- invSToTThresh

    deltaS <- candidateDeltas[mPIndex]

    if (alternative=="less")
      deltaS <- -deltaS

    result[["parameter"]] <- deltaS
    result[["error"]] <- invSafeTTestStatAlpha(candidateDeltas[mPIndex])$estim.prec
  } else {
    # TODO(Alexander): Check relation with forward method, that is, the least conservative test and maximally powered
    # Perhaps trade-off? "inverseMethod" refers to solving minimum of deltaS \mapsto S_{deltaS}^{-1}(1/alpha)
  }
  # TODO(Alexander): By some monotonicity can we only look at the largest or the smallest?
  #
  # designFreqT(deltaMin=deltaMin, alpha=alpha, beta=beta, lowN=lowN, highN=highN)

  names(result[["parameter"]]) <- "deltaS"
  return(result)
}



#' Plots the Sample Sizes Necessary for a Tolerable Alpha and Beta as a Function of deltaMin
#'
#' For given tolerable alpha and beta, (1) the planned sample sizes to using a safe test, (2) the
#' frequentist test, and (3) the average sample size necessary due to optional stopping are plotted
#' as a  function of the minimal clinically relevant standardised effect size deltaMin.
#'
#' @inheritParams designSafeT
#' @inheritParams replicateTTests
#' @param maxN numeric, the maximum number of samples one has budget for to collect data.
#' @param deltaFactor numeric, a factor to robustify the sequential determination (e.g., from deltaTrue = 0.9, to
#' deltaTrue = 0.8) of lowParam.
#' @param nFactor numeric, a factor to robustify the sequential determination (e.g., from deltaTrue = 0.9, to
#' deltaTrue = 0.8) of highN.
#' @param simulateSafeOptioStop logical, if \code{TRUE} then provides the simulation for safe tests.
#' @param logging logical, if \code{TRUE} then output all the safe designs objects including mean n stop if
#' simulateSafeOptioStop equal \code{TRUE}.
#' @param backTest logical, if \code{TRUE} it provides the frequentist sample size necessary to attain the power
#' that the safe test attains due to optional stopping.
#' @param freqPlot logical, if \code{TRUE} plot frequentist sample size profiles.
#'
#' @return Returns a list that contains the planned sample size needed for the frequentist and safe tests as a function
#' of the minimal clinically relevant effect sizes. The returned list contains at least the following components:
#'
#' \describe{
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{beta}{the tolerable type II error provided by the user.}
#'   \item{maxN}{the largest number of samples provided by the user.}
#'   \item{deltaDomain}{vector of the domain of deltaMin.}
#'   \item{allN1PlanFreq}{vector of the planned sample sizes needed for the frequentist test corresponding to
#'   alpha and beta.}
#'   \item{allN1PlanSafe}{vector of the planned sample sizes needed for the safe test corresponding to alpha
#'   and beta.}
#'   \item{allDeltaS}{vector of safe test defining deltaS.}
#' }
#'
#' @export
#'
#' @examples
#' plotSafeTDesignSampleSizeProfile(freqPlot=TRUE, backTest=TRUE)
plotSafeTDesignSampleSizeProfile <- function(alpha=0.05, beta=0.2, maxN=200, lowParam=0.01, highParam=1, tol=0.1,
                                             testType=c("oneSample", "paired", "twoSample"), nsim=1000L,
                                             alternative=c("two.sided", "greater", "less"), ratio=1,
                                             deltaFactor=0.5, nFactor=2, simulateSafeOptioStop=FALSE,
                                             logging=FALSE, backTest=FALSE, seed=NULL, freqPlot=FALSE, pb=TRUE,
                                             ...) {

  stopifnot(lowParam < highParam, alpha > 0, beta > 0, alpha < 1, beta < 1)
  # Order from high to low
  deltaDomain <- -seq(-highParam, -lowParam, by=tol)
  testType <- match.arg(testType)

  result <- list("alpha"=alpha, "beta"=beta, "maxN"=maxN, "deltaDomain"=deltaDomain)

  lastDeltaIndex <- length(deltaDomain)

  if (lastDeltaIndex < 1)
    stop("Either maxN or deltaDomain is too small. Please lower lowParam or make highParam larger")

  paired <- if (testType=="paired") TRUE else FALSE

  allN1PlanFreq <- vector("integer", lastDeltaIndex)


  # 1. Freq design  ---------------------------------------------------------------------
  freqDesign <- list("nPlan"=3)

  for (i in seq.int(lastDeltaIndex)) {
    # TODO(Alexander): Show that as deltaMin decreases that n1PlanFreq increases

    tempLowN <- if (i==1) 3 else freqDesign[["nPlan"]]

    freqDesign <- designFreqT("deltaMin"=deltaDomain[i], "alpha"=alpha, "beta"=beta, "lowN"=tempLowN,
                              "highN"=maxN, "ratio"=ratio)

    if (is.null(freqDesign[["nPlan"]][1]) || is.na(freqDesign[["nPlan"]][1])) {
      lastDeltaIndex <- i-1

      # Note(Alexander): Prune
      allN1PlanFreq <- allN1PlanFreq[1:lastDeltaIndex]
      break()
    }

    allN1PlanFreq[i] <- freqDesign[["nPlan"]][1]
  }


  # #### 1.a. Plots freq
  # graphics::plot(deltaDomain[seq.int(lastDeltaIndex)], allN1PlanFreq[seq.int(lastDeltaIndex)],
  #                lty=3, lwd=2, type="l", col="darkgrey", ylab="n1", xlab=expression(delta["min"]))
  #
  # abline(h=maxN, col="red", lty=2)
  #
  # legend("topright", legend = c("Freq design", "max n"),
  #        col = c("darkgrey", "red"),
  #        lty = c(3, 2), bty="n")

  # 2. Safe design  ---------------------------------------------------------------------
  #
  allDeltaS <- allN1PlanSafe <- vector("integer", lastDeltaIndex)
  allSafeDesignObj <- vector("list", lastDeltaIndex)

  # TODO(Alexander): Show that for fixed theta that freqN < safeN,
  # TODO(Alexander): Instead, of doing this based on the previous one (deltaMin), try to be faster
  # by going around a guessed quantity based on a factor of safeDesignObj$n/freqDesign$n
  #
  # TODO(Alexander): Show that due to monotonicity that we can take "highParam"=deltaDomain[i-1]
  for (i in seq.int(lastDeltaIndex)) {
    if (i==1) {
      tempLowParam <- lowParam
      tempHighParam <- highParam
      tempLowN <- 3
      tempHighN <- maxN
    } else {
      # Note(Alexander): Use previous found parameter times a correction factor as a lowerbound for the search space of
      # parameter
      tempLowParam <- deltaFactor*safeDesignObj[["parameter"]]/deltaDomain[i-1]*deltaDomain[i]

      # Note(Alexander): Use previous true deltaMin > previous parameter as an upper bound
      # TODO(Alexander): Show that as deltaMin decreases that parameter dereases
      tempHighParam <- deltaDomain[i-1]

      # Note(Alexander): Use previously found design
      # TODO(Alexander): Show that as deltaMin decreases that n1Plan increases
      tempLowN <- safeDesignObj[["nPlan"]][1]

      # Note(Alexander): Use previously found n1Plan times a factor as an upper bound
      # TODO(Alexander): Show that as deltaMin decreases that n1Plan increases
      tempHighN <- ceiling(nFactor * safeDesignObj[["nPlan"]][1]/allN1PlanFreq[i-1]*allN1PlanFreq[i])
    }


    safeDesignObj <- designSafeT("deltaMin"=deltaDomain[i], "alpha"=alpha, "beta"=beta, "alternative"=alternative,
                                 "lowParam"=tempLowParam, "highParam"=tempHighParam, "lowN"=tempLowN,
                                 "highN"=tempHighN, "testType"=testType, "ratio"=ratio)

    if (is.null(safeDesignObj[["nPlan"]][1]) || is.na(safeDesignObj[["nPlan"]][1])) {
      lastDeltaIndex <- i-1
      break()
    }

    # TODO(Alexander): Not necessary anymore, with normal function call
    safeDesignObj[["n1PlanFreq"]] <- allN1PlanFreq[i]
    allSafeDesignObj[[i]] <- safeDesignObj
    allN1PlanSafe[i] <- safeDesignObj[["nPlan"]][1]
    allDeltaS[i] <- safeDesignObj[["parameter"]]
  }

  deltaDomain <- deltaDomain[1:lastDeltaIndex]
  allDeltaS <- allDeltaS[1:lastDeltaIndex]

  # TODO(Alexander) Optional?
  # graphics::plot(deltaDomain, allDeltaS)


  maxDeltaDomain <- max(deltaDomain)
  minDeltaDomain <- min(deltaDomain)


  # Store in output
  result[["deltaDomain"]] <- deltaDomain
  result[["allN1PlanFreq"]] <- allN1PlanFreq
  result[["allN1PlanSafe"]] <- allN1PlanSafe
  result[["allDeltaS"]] <- allDeltaS

  # 2.a. Plot Safe -----
  oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
  on.exit(graphics::par(oldPar))

  graphics::plot(deltaDomain, allN1PlanSafe, type="l", col="blue", lty=1, lwd=2, xlim=c(minDeltaDomain, maxDeltaDomain),
                 ylab="n1", xlab=expression(delta["min"]),
                 main=bquote(~alpha == ~.(alpha) ~ "and" ~beta== ~.(beta)))

  if (freqPlot) {
    graphics::lines(deltaDomain, allN1PlanFreq, col="darkgrey", lwd=2, lty=3)
    legendName <- c("Safe design", "Freq design", "max n")
    legendCol <- c("blue", "darkgrey", "red")
    legendLty <- c(1, 3, 2)
  } else {
    legendName <- c("Safe design", "max n")
    legendCol <- c("blue", "red")
    legendLty <- c(1, 2)
  }

  graphics::abline(h=maxN, col="red", lty=2)

  graphics::legend("topright", legend = legendName,
                   col = legendCol, lty = legendLty, bty="n")

  # 3. Run simulations  ---------------------------------------------------------------------
  #
  if (simulateSafeOptioStop) {
    allNMean <- allProbLeqNFreq <- vector("integer", lastDeltaIndex)

    if (backTest)
      allNBack <- vector("integer", lastDeltaIndex)

    if (pb)
      pbOptioStop <- utils::txtProgressBar("style"=1)

    for (i in seq.int(lastDeltaIndex)) {
      safeDesignObj <- allSafeDesignObj[[i]]

      tempLowN <- if (i==1) 3 else simObj[["safeSim"]][["lowN"]]

      simObj <- replicateTTests("n1Plan"=safeDesignObj[["n1Plan"]], "n2Plan"=safeDesignObj[["n2Plan"]],
                                "deltaTrue"=deltaDomain[i], "paired"=paired, "alternative"=alternative,
                                "lowN"=tempLowN, "alpha"=alpha, "parameter"=safeDesignObj[["parameter"]],
                                "n1PlanFreq"=allN1PlanFreq[i], "pb"=FALSE, "nsim"=nsim)
      allNMean[i] <- simObj[["safeSim"]][["nMean"]]
      allProbLeqNFreq[i] <- simObj[["safeSim"]][["probLeqN1PlanFreq"]]

      if (backTest) {
        backFreqDesign <- designFreqT("deltaMin"=deltaDomain[i], "alpha"=alpha,
                                      "beta"=1-safeDesign[["safeSim"]][["powerOptioStop"]],
                                      "alternative"=alternative, "testType"=testType,
                                      "ratio"=ratio)

        safeDesign[["safeSim"]][["nBack"]] <- backFreqDesign[["nPlan"]][1]
        allNBack[i] <- backFreqDesign[["nPlan"]][1]
      } # End back test

      safeDesignObj <- utils::modifyList(safeDesignObj, simObj)
      allSafeDesignObj[[i]] <- safeDesignObj

      if (pb)
        utils::setTxtProgressBar("pb"=pbOptioStop, "value"=i/lastDeltaIndex)

    } # End looping over deltaDomain as deltaTrue

    if (pb)
      close(pbOptioStop)

    result[["allNMean"]] <- allNMean
    result[["allProbLeqNFreq"]] <- allProbLeqNFreq

    if (backTest)
      result[["allNBack"]] <- allNBack

    # 3.a. Plot Sim  -----
    oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
    on.exit(graphics::par(oldPar))

    graphics::plot(deltaDomain, allN1PlanSafe, type="l", col="blue", lty=2, lwd=2, xlim=c(minDeltaDomain, maxDeltaDomain),
                   ylab="n1", xlab=expression(delta["min"]),
                   main=bquote(~alpha == ~.(alpha) ~ "and" ~beta== ~.(beta)))
    graphics::abline(h=maxN, col="red", lty=2)
    graphics::lines(deltaDomain, allNMean, col="black", lwd=2, lty=1)

    if (freqPlot) {
      graphics::lines(deltaDomain, allN1PlanFreq, col="darkgrey", lwd=2, lty=3)
      legendName <- c("Average n", "Safe design", "Freq design", "max n")
      legendCol <- c("black", "blue", "darkgrey", "red")
      legendLty <- c(1, 2, 3, 2)
    } else {
      legendName <- c("Average n", "Safe design", "max n")
      legendCol <- c("black", "blue", "red")
      legendLty <- c(1, 2, 2)
    }

    graphics::legend("topright", legend = legendName, col = legendCol, lty=legendLty, bty="n")
  }

  if (logging)
    result[["allSafeDesignObj"]] <- allSafeDesignObj

  return(result)
}



# 2. Data generating helper functions ------

#' Generates Normally Distributed Data Depending on the Design
#'
#' The designs supported are "oneSample", "paired", "twoSample".
#'
#' @inheritParams replicateTTests
#'
#' @return
#' @return Returns a list of two data matrices contains at least the following components:
#'
#' \describe{
#'   \item{dataGroup1}{a matrix of data dimension nsim by \code{nPlan[1]}.}
#'   \item{dataGroup2}{a matrix of data dimension nsim by \code{nPlan[2]}.}
#' }
#' @export
#'
#' @examples
#' generateTTestData(20, 15)
generateTTestData <- function(nPlan, nsim=1000L, deltaTrue=0, muGlobal=0, sigmaTrue=1, paired=FALSE,
                              seed=NULL) {
  stopifnot(all(nPlan > 0))

  result <- list("dataGroup1"=NULL, "dataGroup2"=NULL)
  set.seed(seed)

  # TODO(Alexander): vector("mode"="list", length=length(nPlan))

  n1Plan <- nPlan[1]

  if (length(nPlan)==1) {
    dataGroup1 <- stats::rnorm("n"=n1Plan*nsim, "mean"=deltaTrue*sigmaTrue, "sd"=sigmaTrue)
    dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nsim)
    dataGroup2 <- NULL
  } else {
    n2Plan <- nPlan[2]

    if (paired) {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nsim, "mean"=muGlobal + deltaTrue*sigmaTrue/sqrt(2), "sd"=sigmaTrue)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nsim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nsim, "mean"=muGlobal - deltaTrue*sigmaTrue/sqrt(2), "sd"=sigmaTrue)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nsim)
    } else {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nsim, "mean"=muGlobal + deltaTrue*sigmaTrue/2, "sd"=sigmaTrue)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nsim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nsim, "mean"=muGlobal - deltaTrue*sigmaTrue/2, "sd"=sigmaTrue)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nsim)
    }
  }

  return(list("dataGroup1"=dataGroup1, "dataGroup2"=dataGroup2))
}

# 3. Inference functions -------

#' Safe Student's T-Test.
#'
#' A safe t-test adapted from  \code{\link[stats]{t.test}} to perform one and two sample t-tests on vectors of data.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param designObj an object obtained from \code{\link{designSafeT}}, or \code{NULL}, when pilot equals \code{TRUE}..
#' @param h0 a number indicating the hypothesised true value of the mean under the null. For the moment h0=0
#' @param paired a logical indicating whether you want a paired t-test.
#' @param varEqual a logical variable indicating whether to treat the two variances as being equal. For the moment,
#' this is always \code{TRUE}.
#' @param pilot a logical indicating whether a pilot study is run. If \code{TRUE}, it is assumed that the number of
#' samples is exactly as planned.
#' @param alpha numeric > 0 only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the alpha of
#' the design object is used instead in constructing the decision rule S > 1/alpha.
#' @param alternative a character only used if pilot equals \code{TRUE}. If pilot equals \code{FALSE}, then the
#' alternative specified by the design object is used instead.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class "safeTest". An object of class "safeTest" is a list containing at least the
#' following components:
#'
#' \describe{
#'   \item{statistic}{the value of the t-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{sValue}{the realised s-value from the safe test.}
#'   \item{confSeq}{To be implemented: a safe confidence interval for the mean appropriate to the specific alternative
#'   hypothesis.}
#'   \item{estimate}{the estimated mean or difference in means or mean difference depending on whether it was a one-
#'   sample test or a two-sample test.}
#'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on whether it was a one-sample
#'   or a two-sample test.}
#'   \item{stderr}{the standard error of the mean (difference), used as denominator in the t-statistic formula.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "safeTDesign" obtained from \code{\link{designSafeT}}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.6, alpha=0.008, alternative="greater",
#' testType="twoSample", ratio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safeTTest(x, y, designObj=designObj)      #0.2959334
#'
#' safeTTest(1:10, y = c(7:20), pilot=TRUE)      # s = 658.69 > 1/alpha
safeTTest <- function(x, y=NULL, designObj=NULL, paired=FALSE, varEqual=TRUE,
                      h0=0, pilot=FALSE, alpha=NULL, alternative=NULL, ...) {
  # TODO(Alexander): Generalise h0 = 0 to other h0

  result <- list("statistic"=NULL, "n"=NULL, "sValue"=NULL, "confSeq"=NULL, "estimate"=NULL,
                 "alternative"=NULL, "testType"=NULL, "dataName"=NULL, "h0"=h0, "stderr"=NULL,
                 "call"=sys.call())
  class(result) <- "safeTest"

  if (is.null(designObj) && !pilot) {
    stop("No design given and not indicated that this is a pilot study. Run design first and provide ",
         "this to safeTTest/safe.t.test, or run safeTTest with pilot=TRUE")
  }

  if (!is.null(designObj)) {
    checkDoubleArgumentsDesignObject(designObj, "alternative"=alternative, "alpha"=alpha)

    if (names(designObj[["parameter"]]) != "deltaS")
      warning("The provided design is not constructed for the z-test,",
              "please use designSafeT() instead. The test results might be invalid.")
  }

  n1 <- length(x)

  if (is.null(y)) {
    n <- n1
    names(n) <- "n1"
    testType <- "oneSample"
  } else {
    n2 <- length(y)
    n <- c(n1, n2)
    names(n) <- c("n1", "n2")

    if (paired)
      testType <- "paired"
    else
      testType <- "twoSample"
  }


  if (pilot) {
    if (is.null(alternative))
      alternative <- "two.sided"
    else {
      if (!(alternative %in% c("two.sided", "greater", "less")))
        stop('Provided alternative must be one of "two.sided", "greater", or "less".')
    }

    if (is.null(alpha))
      alpha <- 0.05

    designObj <- designPilotSafeT("nPlan"=n, "alpha"=alpha, "alternative"=alternative, "paired"=paired)
  }

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  alternative <- designObj[["alternative"]]
  alpha <- designObj[["alpha"]]

  freqObject <- try(stats::t.test("x"=x, "y"=y, "alternative"=alternative, "mu"=h0,
                                  "paired"=paired, "var.equal"=varEqual))
  t <- freqObject[["statistic"]]

  if (isTryError(freqObject))
    stop("Data error: could not compute the t-statistic with t.test: ", freqObject[1])

  # TODO(Alexander): Save result, perhaps save freqObject
  #
  sValue <- safeTTestStat("t"=t, "parameter"=designObj[["parameter"]], "n1"=n[1], "n2"=n[2],
                          "alternative"=alternative, "paired"=paired)

  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  result[["statistic"]] <- t
  result[["parameter"]] <- designObj[["deltaS"]]
  result[["estimate"]] <- freqObject[["estimate"]]
  result[["stderr"]] <- freqObject[["stderr"]]
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["freqObject"]] <- freqObject
  result[["testType"]] <- testType
  result[["n"]] <- n
  result[["sValue"]] <- sValue
  result[["h0"]] <- "mu"

  return(result)
}


#' Alias for safeTTest
#'
#' @rdname safeTTest
#'
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. For the moment,
#' this is always \code{TRUE}.
#'
#' @export
#'
#' @examples
#' designObj <- designSafeT(deltaMin=0.6, alpha=0.008, alternative="greater",
#' testType="twoSample", ratio=1.2)
#'
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' safe.t.test(x, y, alternative="greater", designObj=designObj)      #0.2959334
#'
#' safe.t.test(1:10, y = c(7:20), pilot=TRUE)      # s = 658.69 > 1/alpha
safe.t.test <- function(x, y=NULL, designObj=NULL, paired=FALSE, var.equal=TRUE, h0=0,
                        pilot=FALSE, alpha=NULL, alternative=NULL, ...) {
  result <- safeTTest("x"=x, "y"=y, "designObj"=designObj, "h0"=h0, "paired"=paired, "varEqual"=var.equal,
                      "pilot"=pilot, "alpha"=alpha, "alternative"=alternative, ...)
  if (is.null(y))
    dataName <- as.character(sys.call())[2]
  else
    dataName <- paste(as.character(sys.call())[2], "and", as.character(sys.call())[3])

  result[["dataName"]] <- dataName
  return(result)
}

#' Computes S-Values Based on the T-Statistic
#'
#' A summary stats version of \code{\link{safeTTest}} with the data replaced by t, n1 and n2, and the
#' design object by deltaS.
#'
#' @param t numeric that represents the observed t-statistic.
#' @param parameter numeric this defines the safe test S, i.e., a likelihood ratio of t distributions with in
#' the denominator the likelihood with delta = 0 and in the numerator an average likelihood defined by
#' 1/2 time the likelihood at the non-centrality parameter sqrt(nEff)*parameter and 1/2 times the likelihood at
#' the non-centrality parameter -sqrt(nEff)*parameter.
#' @param n1 integer that represents the size in a one-sample t-test, (n2=\code{NULL}). When n2 is not \code{NULL},
#' this specifies the size of the first sample for a two-sample test.
#' @param n2 an optional integer that specifies the size of the second sample. If it's left unspecified, thus,
#' \code{NULL} it implies that the t-statistic is based on one-sample.
#' @param tDensity Uses the the representation of the safe t-test as the likelihood ratio of t densities.
#' @inherit safeTTest
#'
#' @return Returns a numeric that represent the s10, that is, the s-value in favour of the alternative over the null
#'
#' @export
#'
#' @examples
#' safeTTestStat(t=1, n1=100, 0.4)
#' safeTTestStat(t=3, n1=100, parameter=0.3)
safeTTestStat <- function(t, parameter, n1, n2=NULL, alternative=c("two.sided", "less", "greater"), tDensity=FALSE,
                          paired=FALSE, ...) {
  # TODO(Alexander):
  #   One-sided not as stable as two-sided due to hypergeo::genhypergeo for the odd component
  #   1. Use Kummer's transform again (??)
  #   2. Switch to numerical integration. Boundary case
  #
  # safeTTestStat(t=-3.1878, parameter=0.29, n1=315, alternative="greater")
  # safeTTestStat(t=-3.1879, parameter=0.29, n1=315, alternative="greater")
  # safeTTestStat(t=-3.188, parameter=0.29, n1=315, alternative="greater")
  alternative <- match.arg(alternative)
  deltaS <- parameter

  if (is.null(n2) || is.na(n2) || paired==TRUE) {
    nEff <- n1
    nu <- n1-1
  } else {
    nEff <- (1/n1+1/n2)^(-1)
    nu <- n1+n2-2
  }

  # TODO(Alexander): This is not necessarily correct, since the correct one is
  # (2*(y1-y2))^(-1)*(1+sign(y1-y2)*
  #   (pnorm(deltaS/2, sd=1/sqrt(2))-pnorm(-deltaS/2, sd=1/sqrt(2))) #erf(deltaS/2)
  #
  #
  # Problem: t=0, n=1
  #
  # PERHAPS add warning and assume t is then y2-y1, otherwise t must be Inf
  if (nu == 0) {
    if (t==0)
      return(NA)
    else
      return(1)
  }

  if (tDensity) {
    if (alternative=="two.sided") {
      logTerm1 <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)
      logTerm2 <- stats::dt(t, df=nu, ncp=-sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)

      result <- exp(logTerm1+logTerm2)/2
    } else {
      result <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS)/stats::dt(t, df=nu, ncp=0)
    }
  } else {
    a <- t^2/(nu+t^2)
    expTerm <- exp((a-1)*nEff*deltaS^2/2)

    zeroIndex <- abs(expTerm) < .Machine$double.eps
    result <- vector("numeric", length(expTerm))

    zArg <- (-1)*a*nEff*deltaS^2/2
    zArg <- zArg[!zeroIndex]
    # Note(Alexander): This made the vector shorter. Only there where expTerm is non-zero will we evaluate
    # the hypergeometric functions

    aKummerFunction <- Re(hypergeo::genhypergeo(U=-nu/2, L=1/2, zArg))

    if (alternative=="two.sided") {
      result[!zeroIndex] <- expTerm[!zeroIndex] * aKummerFunction
    } else {
      bKummerFunction <- exp(lgamma(nu/2+1)-lgamma((nu+1)/2))*sqrt(2*nEff)*deltaS*t/sqrt(t^2+nu)[!zeroIndex] *
        Re(hypergeo::genhypergeo(U=(1-nu)/2, L=3/2, zArg))
      result[!zeroIndex] <- expTerm[!zeroIndex]*(aKummerFunction + bKummerFunction)
    }
  }

  if (result < 0) {
    warning("Overflow: s-value smaller than 0")
    result <- 2^(-15)
  }
  return(result)
}

#' safeTTestStat Subtracted with 1/alpha.
#'
#' This is basically just \code{\link{safeTTestStat}} - 1/alpha. This function is used for root finding for
#' pilot designs.
#'
#' @inheritParams safeTTest
#' @inherit safeTTestStat
#'
#' @return Returns a numeric that represent the s10 - 1/alpha, that is, the s-value in favour of the
#' alternative over the null - 1/alpha.
#'
safeTTestStatAlpha <- function(t, parameter, n1, n2=NULL, alpha, alternative="two.sided", tDensity=FALSE) {
  safeTTestStat("t"=t, "parameter"=parameter, "n1"=n1, "n2"=n2, "alternative"=alternative, "tDensity"=tDensity) - 1/alpha
}
