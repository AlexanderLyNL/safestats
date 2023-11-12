#' Helper function that extract the results for the design scenario 1a: Target nPlan
#'
#' @param samplingResult output from sampling functions such as computeNPlanSafeZ and computeNPlanSafeT
#' @param esMin numeric that defines the minimal clinically relevant effect size,
#' e.g. meanDiffMin for the z-test, or deltaMin for the t-test.
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control necessary to calculate both "n"
#' and "phiS". Note that 1-beta defines the power.
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
#' is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param testType either one of "oneSample", "paired", "twoSample".
#'
#' @return a list of partial results for the design scenario 1a
#' @export
#'
#' @examples
#'
#' samplingResult <- computeNPlanSafeZ(0.7, nSim=10, nMax=20)
#' result <- designSafe1aHelper(samplingResult, 0.7, 0.2, 1)
designSafe1aHelper <- function(
    samplingResult, esMin, beta, ratio,
    testType=c("oneSample", "paired","twoSample")) {

  testType <- match.arg(testType)

  result <- list("parameter"=NULL, "esMin"=esMin, "beta"=beta,
                 "nPlan"=NULL, "nPlanTwoSe"=NULL, "nPlanBatch"=NULL,
                 "nMean"=NULL, "nMeanTwoSe"=NULL,
                 "bootObjN1Plan"=NULL, "bootObjN1Mean"=NULL,
                 "samplePaths"=NULL, "breakVector"=NULL,
                 "note"=NULL)

  nPlanBatch <- samplingResult[["nPlanBatch"]]
  bootObjN1Plan <- samplingResult[["bootObjN1Plan"]]
  bootObjN1Mean <- samplingResult[["bootObjN1Mean"]]

  if (testType=="oneSample") {
    nPlan <- samplingResult[["n1Plan"]]
    names(nPlan) <- "nPlan"
    nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]

    nMean <- samplingResult[["n1Mean"]]
    names(nMean) <- "nMean"
    nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]

    note <- paste0("If it is only possible to look at the data once, ",
                   "then nPlan = ", nPlanBatch, ".")
  } else if (testType=="paired") {
    nPlan <- c(samplingResult[["n1Plan"]], samplingResult[["n1Plan"]])
    names(nPlan) <- c("n1Plan", "n2Plan")

    nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
    nPlanTwoSe <- c(nPlanTwoSe, nPlanTwoSe)

    nMean <- c(samplingResult[["n1Mean"]], samplingResult[["n1Mean"]])
    names(nMean) <- c("n1Mean", "n2Mean")
    nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]
    nMeanTwoSe <- c(nMeanTwoSe, nMeanTwoSe)

    note <- paste0("If it is only possible to look at the data once, ",
                   "then n1Plan = ", nPlanBatch[1], " and n2Plan = ",
                   nPlanBatch[2], ".")
  } else if (testType=="twoSample") {
    nPlan <- c(samplingResult[["n1Plan"]], ceiling(ratio*samplingResult[["n1Plan"]]))
    names(nPlan) <- c("n1Plan", "n2Plan")
    nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
    nPlanTwoSe <- c(nPlanTwoSe, ratio*nPlanTwoSe)

    nMean <- c(samplingResult[["n1Mean"]], ceiling(ratio*samplingResult[["n1Mean"]]))
    names(nMean) <- c("n1Mean", "n2Mean")
    nMeanTwoSe <- 2*bootObjN1Mean[["bootSe"]]
    nMeanTwoSe <- c(nMeanTwoSe, ratio*nMeanTwoSe)

    note <- paste0("If it is only possible to look at the data once, ",
                   "then n1Plan = ", nPlanBatch[1], " and n2Plan = ",
                   nPlanBatch[2], ".")
  }

  # Fill results

  result[["parameter"]] <- samplingResult[["parameter"]]
  result[["nPlanBatch"]] <- nPlanBatch
  result[["samplePaths"]] <- samplingResult[["samplePaths"]]
  result[["breakVector"]] <- samplingResult[["breakVector"]]

  result[["bootObjN1Plan"]] <- bootObjN1Plan
  result[["bootObjN1Mean"]] <- bootObjN1Mean

  result[["nPlan"]] <- nPlan
  result[["nPlanTwoSe"]] <- nPlanTwoSe
  result[["nMean"]] <- nMean
  result[["nMeanTwoSe"]] <- nMeanTwoSe
  result[["note"]] <- note

  return(result)
}


#' Helper function that extract the results for the design scenario 1a: Target nPlan
#'
#' @param samplingResult output from sampling functions such as computeNPlanSafeZ and computeNPlanSafeT
#' @param esMin numeric that defines the minimal clinically relevant effect size,
#' e.g. meanDiffMin for the z-test, or deltaMin for the t-test.
#' @param nPlan vector of max length 2 representing the planned sample sizes.
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1. If testType
#' is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param testType either one of "oneSample", "paired", "twoSample".
#'
#' @return a list of partial results for the design scenario 1a
#' @export
#'
#' @examples
#'
#' samplingResult <- computeNPlanSafeZ(0.7, nSim=10, nMax=20)
#' result <- designSafe1aHelper(samplingResult, 0.7, 0.2, 1)
designSafe2Helper <- function(
    samplingResult, esMin, nPlan, ratio,
    testType=c("oneSample", "paired","twoSample")) {

  testType <- match.arg(testType)

  result <- list(
    "parameter"=NULL, "esMin"=esMin, "nPlan"=nPlan,
    "beta"=NULL, "betaTwoSe"=NULL, "bootObjBeta"=NULL,
    "logImpliedTarget"=NULL, "logImpliedTargetTwoSe"=NULL,
    "bootObjLogImpliedTarget"=NULL,
    "samplePaths"=NULL, "breakVector"=NULL)

  result[["parameter"]] <- samplingResult[["parameter"]]
  result[["ratio"]] <- ratio

  result[["samplePaths"]] <- samplingResult[["samplePaths"]]
  result[["breakVector"]] <- samplingResult[["breakVector"]]

  bootObjBeta <- samplingResult[["bootObjBeta"]]

  result[["beta"]] <- samplingResult[["beta"]]
  result[["bootObjBeta"]] <- bootObjBeta
  result[["betaTwoSe"]] <- 2*bootObjBeta[["bootSe"]]

  bootObjLogImpliedTarget <- samplingResult[["bootObjLogImpliedTarget"]]

  result[["logImpliedTarget"]] <- samplingResult[["logImpliedTarget"]]
  result[["bootObjLogImpliedTarget"]] <- bootObjLogImpliedTarget
  result[["logImpliedTargetTwoSe"]] <- 2*bootObjLogImpliedTarget[["bootSe"]]

  return(result)
}

# ---------- Boot helpers --------

#' Computes the bootObj for sequential sampling procedures regarding nPlan, beta, the implied target
#'
#' @inheritParams designSafeZ
#' @param values numeric vector. If objType equals "nPlan" or "beta" then values should be stopping times,
#' if objType equals "logImpliedTarget" then values should be eValues.
#' @param nBoot integer > 0 representing the number of bootstrap samples
#' to estimate the uncertainty of various estimates.
#' @param nPlan numeric vector of length at most 2 representing the planned sample size(s).
#' @param objType character string either "nPlan", "nMean", "beta", "betaFromEValues", "expectedStopTime" or "logImpliedTarget".
#'
#' @return bootObj
#' @export
#'
#' @examples
#' computeBootObj(1:100, objType="nPlan", beta=0.3)
computeBootObj <- function(
    values, beta=NULL, nPlan=NULL,
    nBoot=1e3L, alpha=NULL,
    objType=c("nPlan", "nMean", "beta", "betaFromEValues",
              "logImpliedTarget", "expectedStopTime")) {
  objType <- match.arg(objType)

  if (objType=="beta") {
    if (is.null(nPlan) || nPlan <= 0)
      stop("Please provide an nPlan > 0")

    times <- values
    stopifnot(nPlan > 0)

    bootObj <- boot::boot(times,
                          function(x, idx) {
                            1-mean(x[idx] <= nPlan)
                          },  R = nBoot)
  } else if (objType =="betaFromEValues") {
    if (is.null(alpha) || alpha <= 0 || alpha >= 1)
      stop("Please provide an alpha in (0, 1)")

    eValues <- values

    bootObj <- boot::boot(
      data = eValues,
      statistic = function(x, idx) {
        mean(x[idx] >= 1/alpha)
      },
      R = nBoot
    )

  } else if (objType=="nPlan") {
    if (is.null(beta) || beta <= 0 || beta >= 1)
      stop("Please provide a beta in (0, 1)")

    times <- values
    bootObj <- boot::boot(times,
                          function(x, idx) {
                            quantile(x[idx], prob=1-beta, names=FALSE)
                          }, R = nBoot)
  } else if (objType=="nMean") {
    if (is.null(nPlan[1]) || nPlan[1] <= 0)
      stop("Please provide a positive nPlan")

    times <- values

    times[times > nPlan[1]] <- nPlan[1]

    bootObj <- boot::boot(times,
                          function(x, idx) {
                            mean(x[idx])
                          }, R = nBoot)
  } else if (objType=="logImpliedTarget") {
    eValues <- values
    stopifnot(eValues > 0)

    bootObj <- boot::boot(eValues,
                          function(x, idx) {
                            mean(log(x[idx]))
                          }, R = nBoot)
  } else if (objType=="expectedStopTime") {
    times <- values
    bootObj <- boot::boot(times,
                          function(x, idx) {
                            mean(x[idx])
                          }, R = nBoot)
  }

  bootObj[["bootSe"]] <- sd(bootObj[["t"]])
  return(bootObj)
}


#' Helper function to compute uncertainty regarding nPlan estimates
#'
#' @inheritParams designSafe1aHelper
#' @inheritParams computeBootObj
#'
#' @param parameter numeric > 0, the safe test defining parameter.
#' @param nPlanBatch integer, the sample size needed in a batch design
#' to reach the targeted power=1-beta with tolerable type I error alpha
#'
#' @return list with bootstrap objects
#' @export
#'
#' @examples
#' samplingResult <- sampleStoppingTimesSafeT(0.7, nSim=10, nMax=20)
#' result <- computeNPlanBootstrapper(samplingResult, 0.7, 0.2, 20, nBoot=1e2)
computeNPlanBootstrapper <- function(
    samplingResult, parameter,
    beta, nPlanBatch, nBoot) {

  times <- samplingResult[["stoppingTimes"]]

  bootObjN1Plan <- computeBootObj(
    "values"=times, "objType"="nPlan",
    "beta"=beta, "nBoot"=nBoot)

  n1Plan <- ceiling(bootObjN1Plan[["t0"]])

  bootObjN1Mean <- computeBootObj(
    "values"=times, "objType"="nMean",
    "nPlan"=n1Plan, "nBoot"=nBoot)

  n1Mean <- ceiling(bootObjN1Mean[["t0"]])

  result <- list("n1Plan" = n1Plan, "bootObjN1Plan" = bootObjN1Plan,
                 "n1Mean"=n1Mean, "bootObjN1Mean"=bootObjN1Mean,
                 "nPlanBatch"=nPlanBatch, "parameter"=parameter,
                 "samplePaths"=samplingResult[["samplePaths"]],
                 "breakVector"=samplingResult[["breakVector"]])
}



#' Helper function to compute uncertainty regarding nPlan estimates
#'
#' @inheritParams designSafe2Helper
#' @inheritParams computeBootObj
#' @inheritParams computeNPlanBootstrapper
#'
#' @return list with bootstrap objects
#' @export
#'
#' @examples
#' samplingResult <- sampleStoppingTimesSafeT(0.7, nSim=10, nMax=20)
#' result <- computeNPlanBootstrapper(samplingResult, 0.7, 0.2, 20, nBoot=1e2)
computeBetaBootstrapper <- function(
    samplingResult, parameter,
    nPlan, nBoot) {

  times <- samplingResult[["stoppingTimes"]]

  # Note(Alexander): Break vector is 1 whenever the sample path did not stop
  breakVector <- samplingResult[["breakVector"]]

  # Note(Alexander): Setting the stopping time to Inf for these paths doesn't matter for the quantile
  times[as.logical(breakVector)] <- Inf

  bootObjBeta <- computeBootObj(
    "values"=times, "objType"="beta",
    "nPlan"=nPlan[1], "nBoot"=nBoot)

  eValuesAtNMax <- samplingResult[["eValuesAtNMax"]]

  bootObjLogImpliedTarget <- computeBootObj(
    "values"=eValuesAtNMax, "objType"="logImpliedTarget",
    "nBoot"=nBoot)

  result <- list("beta" = bootObjBeta[["t0"]],
                 "bootObjBeta" = bootObjBeta,
                 "logImpliedTarget"=bootObjLogImpliedTarget[["t0"]],
                 "bootObjLogImpliedTarget"=bootObjLogImpliedTarget,
                 "samplePaths"=samplingResult[["samplePaths"]],
                 "breakVector"=samplingResult[["breakVector"]],
                 "parameter"=parameter)

  return(result)
}

#' Construct a list to be set in the sampleStoppingTimes... function
#'
#' @return a list with names
#' @export
#'
#' @examples
#' obj <- constructSampleStoppingTimesObj()
#'
constructSampleStoppingTimesObj <- function(nSim=1e3L, nMax=1e3L,
                                            wantEValuesAtNMax=FALSE,
                                            wantSamplePaths=TRUE) {

  stoppingTimes <- breakVector <- integer(nSim)
  eValuesStopped <- numeric(nSim)

  eValuesAtNMax <- if (wantEValuesAtNMax) numeric(nSim) else NULL
  samplePaths <- if (wantSamplePaths) matrix(nrow=nSim, ncol=nMax[1]) else NULL

  result <- list("parameter"=NULL,
                 "stoppingTimes"=stoppingTimes, "breakVector"=breakVector,
                 "eValuesStopped"=eValuesStopped, "eValuesAtNMax"=eValuesAtNMax,
                 "samplePaths"=samplePaths, "n1Vector"=NULL, "ratio"=NULL,
                 "simData"=NULL)
  return(result)
}
