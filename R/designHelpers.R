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
designSavi1aHelper <- function(
    samplingResult, esMin, power, ratio, beta=NULL,
    testType=c("oneSample", "paired","twoSample")) {

  testType <- match.arg(testType)

  result <- list("parameter"=NULL, "esMin"=esMin, "power"=power,
                 "nPlan"=NULL, "nPlanTwoSe"=NULL, "nPlanBatch"=NULL,
                 "nMean"=NULL, "nMeanTwoSe"=NULL,
                 "bootObjN1Plan"=NULL, "bootObjN1Mean"=NULL,
                 "samplePaths"=NULL, "breakVector"=NULL,
                 "futility"=FALSE, "futilityResult"=NULL,
                 "simData"=NULL, "beta"=beta, "note"=NULL)

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
    nPlan <- c(samplingResult[["n1Plan"]], ceil(ratio*samplingResult[["n1Plan"]]))
    names(nPlan) <- c("n1Plan", "n2Plan")
    nPlanTwoSe <- 2*bootObjN1Plan[["bootSe"]]
    nPlanTwoSe <- c(nPlanTwoSe, ratio*nPlanTwoSe)

    nMean <- c(samplingResult[["n1Mean"]], ceil(ratio*samplingResult[["n1Mean"]]))
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

  result[["futilityResult"]] <- samplingResult[["futilityResult"]]
  result[["simData"]] <- samplingResult[["simData"]]

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
designSavi2Helper <- function(
    samplingResult, esMin, nPlan, ratio,
    testType=c("oneSample", "paired","twoSample")) {

  testType <- match.arg(testType)

  result <- list(
    "parameter"=NULL, "esMin"=esMin, "nPlan"=nPlan,
    "power"=NULL, "powerTwoSe"=NULL, "bootObjPower"=NULL,
    "logImpliedTarget"=NULL, "logImpliedTargetTwoSe"=NULL,
    "bootObjLogImpliedTarget"=NULL, "simData"=NULL,
    "beta"=NULL, "betaTwoSe"=NULL, "bootObjBeta"=NULL,
    "samplePaths"=NULL, "breakVector"=NULL)

  result[["parameter"]] <- samplingResult[["parameter"]]
  result[["ratio"]] <- ratio

  result[["samplePaths"]] <- samplingResult[["samplePaths"]]
  result[["breakVector"]] <- samplingResult[["breakVector"]]

  someBeta <- samplingResult[["beta"]]

  if (!is.null(someBeta)) {
    bootObjBeta <- samplingResult[["bootObjBeta"]]

    result[["beta"]] <- someBeta
    result[["bootObjBeta"]] <- bootObjBeta
    result[["betaTwoSe"]] <- 2*bootObjBeta[["bootSe"]]
  }

  somePower <- samplingResult[["power"]]

  if (!is.null(somePower)) {
    bootObjPower <- samplingResult[["bootObjPower"]]

    result[["power"]] <- somePower
    result[["bootObjPower"]] <- bootObjPower
    result[["powerTwoSe"]] <- 2*bootObjPower[["bootSe"]]
  }

  bootObjLogImpliedTarget <- samplingResult[["bootObjLogImpliedTarget"]]

  result[["logImpliedTarget"]] <- samplingResult[["logImpliedTarget"]]
  result[["bootObjLogImpliedTarget"]] <- bootObjLogImpliedTarget
  result[["logImpliedTargetTwoSe"]] <- 2*bootObjLogImpliedTarget[["bootSe"]]

  result[["simData"]] <- samplingResult[["simData"]]
  result[["futilityResult"]] <- samplingResult[["futilityResult"]]

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
    values, power=NULL, nPlan=NULL,
    nBoot=1e3L, alpha=NULL, beta=NULL,
    objType=c("nPlan", "nMean", "power", "beta", "betaFromEValues",
              "logImpliedTarget", "expectedStopTime")) {
  objType <- match.arg(objType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (objType=="power") {
    if (is.null(nPlan) || nPlan <= 0)
      stop("Please provide an nPlan > 0")

    times <- values
    stopifnot(nPlan > 0)

    bootObj <- try(
      boot::boot(times, function(x, idx) {
        mean(x[idx] <= nPlan)
      },  R = nBoot)
    )

    j <- 1

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(times, function(x, idx) {
          mean(x[idx] <= nPlan)
        },  R = nBoot/2^j)
      )

      j <- j+1
    }
  } else if (objType=="beta") {
    if (is.null(nPlan) || nPlan <= 0)
      stop("Please provide an nPlan > 0")

    times <- values
    stopifnot(nPlan > 0)

    bootObj <- try(
      boot::boot(times, function(x, idx) {
        1-mean(x[idx] <= nPlan)
      },  R = nBoot)
    )

    j <- 1

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(times, function(x, idx) {
          1-mean(x[idx] <= nPlan)
        },  R = nBoot/2^j)
      )

      j <- j+1
    }
  } else if (objType =="betaFromEValues") {
    if (is.null(alpha) || alpha <= 0 || alpha >= 1)
      stop("Please provide an alpha in (0, 1)")

    eValues <- values

    bootObj <- try(
      boot::boot(data = eValues,
                 statistic = function(x, idx) {
                   mean(x[idx] >= 1/alpha)
                 }, R = nBoot)
    )

    j <- 1

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(data = eValues,
                   statistic = function(x, idx) {
                     mean(x[idx] >= 1/alpha)
                   }, R = nBoot/2^j)
      )

      j <- j+1
    }

  } else if (objType=="nPlan") {
    if (is.null(beta)) {
      if (is.null(power) || power <= 0 || power >= 1)
        stop("Please provide a targeted power in (0, 1)")
    }

    times <- values

    bootObj <- try(
      boot::boot(times, function(x, idx) {
        stats::quantile(x[idx], prob=power, names=FALSE)
      } , R = nBoot)
    )

    j <- 1

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(times, function(x, idx) {
          stats::quantile(x[idx], prob=power, names=FALSE)
        } , R = nBoot/2^j)
      )

      j <- j+1
    }

  } else if (objType=="nMean") {
    if (is.null(nPlan[1]) || nPlan[1] <= 0)
      stop("Please provide a positive nPlan")

    times <- values

    times[times > nPlan[1]] <- nPlan[1]

    bootObj <- try(
      boot::boot(times, function(x, idx) {
        mean(x[idx])
      }, R = nBoot)
    )

    j <- 1

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(times, function(x, idx) {
          mean(x[idx])
        }, R = nBoot/2^j)
      )

      j <- j+1
    }
  } else if (objType=="logImpliedTarget") {
    eValues <- values
    stopifnot(eValues > 0)

    bootObj <- try(
      boot::boot(eValues, function(x, idx) {
        mean(log(x[idx]))
      } , R = nBoot)
    )

    j <- 1

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(eValues, function(x, idx) {
          mean(log(x[idx]))
        } , R = nBoot/2^j)
      )
      j <- j+1
    }
  } else if (objType=="expectedStopTime") {
    times <- values
    bootObj <- try(
      boot::boot(times, function(x, idx) {
        mean(x[idx])
      }, R = nBoot)
    )

    while (isTryError(bootObj) && j < 21) {
      bootObj <- try(
        boot::boot(times, function(x, idx) {
          mean(x[idx])
        }, R = nBoot/2^j)
      )
      j <- j+1
    }
  }

  bootObj[["bootSe"]] <- stats::sd(bootObj[["t"]])
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
    power, nPlanBatch, nBoot, beta=NULL) {

  times <- samplingResult[["stoppingTimes"]]

  power <- matchPowerWith("power"=power, "beta"=beta)

  bootObjN1Plan <- computeBootObj(
    "values"=times, "objType"="nPlan",
    "power"=power, "nBoot"=nBoot, "beta"=beta)

  n1Plan <- ceil(bootObjN1Plan[["t0"]])

  futResult <- samplingResult[["futilityResult"]]

  if (!is.null(futResult)) {
    futIndex <- Matrix::which(samplingResult[["breakVector"]]==-1)
    times[futIndex] <- as.numeric(futResult[["stoppingTimes"]])[futIndex]
  }

  bootObjN1Mean <- computeBootObj(
    "values"=times, "objType"="nMean",
    "nPlan"=n1Plan, "nBoot"=nBoot)

  n1Mean <- ceil(bootObjN1Mean[["t0"]])

  result <- list("n1Plan" = n1Plan, "bootObjN1Plan" = bootObjN1Plan,
                 "n1Mean"=n1Mean, "bootObjN1Mean"=bootObjN1Mean,
                 "nPlanBatch"=nPlanBatch, "parameter"=parameter,
                 "samplePaths"=samplingResult[["samplePaths"]],
                 "breakVector"=samplingResult[["breakVector"]],
                 "simData"=samplingResult[["simData"]],
                 "futilityResult"=samplingResult[["futilityResult"]])
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
  times[Matrix::which(breakVector!=0)] <- Inf

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
                 "simData"=samplingResult[["simData"]],
                 "parameter"=parameter, "futilityResult"=samplingResult[["futilityResult"]])

  return(result)
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
computePowerBootstrapper <- function(
    samplingResult, parameter,
    nPlan, nBoot) {

  times <- samplingResult[["stoppingTimes"]]

  # Note(Alexander): Break vector is 1 whenever the sample path did not stop
  breakVector <- samplingResult[["breakVector"]]

  # Note(Alexander): Setting the stopping time to Inf for these paths doesn't matter for the quantile
  times[Matrix::which(breakVector!=0)] <- Inf

  bootObjPower <- computeBootObj(
    "values"=times, "objType"="power",
    "nPlan"=nPlan[1], "nBoot"=nBoot)

  eValuesAtNMax <- samplingResult[["eValuesAtNMax"]]

  bootObjLogImpliedTarget <- computeBootObj(
    "values"=eValuesAtNMax, "objType"="logImpliedTarget",
    "nBoot"=nBoot)

  result <- list("power" = bootObjPower[["t0"]],
                 "bootObjPower" = bootObjPower,
                 "logImpliedTarget"=bootObjLogImpliedTarget[["t0"]],
                 "bootObjLogImpliedTarget"=bootObjLogImpliedTarget,
                 "samplePaths"=samplingResult[["samplePaths"]],
                 "breakVector"=samplingResult[["breakVector"]],
                 "simData"=samplingResult[["simData"]],
                 "parameter"=parameter, "futilityResult"=samplingResult[["futilityResult"]])

  return(result)
}

#' Construct a list to be set in the sampleStoppingTimes... function
#' @param nSim integer > 0, the number of simulations needed to compute power or the number of samples paths
#' for the safe z test under continuous monitoring.
#' @param nMax integer > 0, maximum sample size of the (first) sample in each sample path.
#' @param wantEValuesAtNMax logical. If \code{TRUE}, then compute eValues at nMax. Default \code{FALSE}.
#' @param wantSamplePaths logical. If \code{TRUE}, then output the (stopped) sample paths. Default \code{TRUE}.
#'
#' @return a list with names
#' @export
#'
#' @examples
#' obj <- constructSampleStoppingTimesList()
constructSampleStoppingTimesList <- function(nSim=1e3L, nMax=1e3L,
                                            wantEValuesAtNMax=FALSE,
                                            wantSamplePaths=TRUE) {

  stoppingTimes <- integer(nSim)
  eValuesStopped <- numeric(nSim)
  breakVector <- Matrix::sparseVector(x=0, i=1, length=nSim)

  eValuesAtNMax <- if (wantEValuesAtNMax) numeric(nSim) else NULL
  samplePaths <- if (wantSamplePaths) matrix(nrow=nSim, ncol=nMax[1]) else NULL

  result <- list("parameter"=NULL,
                 "stoppingTimes"=stoppingTimes, "breakVector"=breakVector,
                 "eValuesStopped"=eValuesStopped, "eValuesAtNMax"=eValuesAtNMax,
                 "samplePaths"=samplePaths, "n1Vector"=NULL, "ratio"=NULL,
                 "simData"=NULL)
  return(result)
}


# Match args ----

#' Checks and outputs a threshold for a futility analysis
#'
#' @param betaFutility numeric > 0 and < 1, used to set the threshold of
#' a futility procedure
#' @param beta numeric > 0 and < 1, a tolerable type II error, used to set
#' the threshold if betaFutility is not given
#' @param betaDefault numeric > 0 and < 1 a default value (0.2) to run
#' a futility procedure
#'
#' @returns a betaFutility threshold
#' @export
#'
#' @examples
#' matchBetaFutilityWith(0.3)
matchBetaFutilityWith <- function(betaFutility, power, beta=NULL, betaDefault=0.2) {
  if (!is.null(betaFutility)) {
    stopifnot(betaFutility >0, betaFutility < 1)
    return(betaFutility)
  }

  if (!is.null(power)) {
    stopifnot(power > 0, power < 1)
    return(1-power)
  }

  if (!is.null(beta) && length(beta)!=0) {
    stopifnot(beta >0, beta < 1)
    return(beta)
  }

  warning("To run a futility procedure ",
          "a betaFutility threshold needs to be specified ",
          "by default betaFutility <- ", betaDefault)
  betaFutility <- betaDefault

  return(betaFutility)
}

#' Match the parameter of a savi z or t-test
#'
#' Based on the minimal clinically relevant effect size esMin,
#' sigma (for z-tests), alternative and eType
#'
#' @inherit designSaviZ
#' @param esMin numeric: meanDiffMin for z-tests, or deltaMin for t-tests
#' @param type character representing analysis type "z" or "t"
#'
#' @returns the parameter, a numeric value
#' @export
#'
#' @examples
#' matchParameterZFrom(0.4)
matchEParameterWith <- function(esMin, analysisType=c("z", "t"),
                                sigma=1,
                                alternative=c("twoSided", "greater", "less"),
                                eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
                                parameter=NULL) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  # TODO(Alexander):
  #
  if (!is.null(parameter))
    return(parameter)

  if (analysisType=="z") {
    parameter <- switch(eType,
                        "mom"=1/2*(esMin/sigma)^2,
                        "eGauss"=(esMin/sigma)^2,
                        "imom"=(esMin/sigma)^2,
                        "eCauchy"=abs(esMin/sigma),
                        "grow"=abs(esMin))
  }

  if (analysisType=="t") {
    parameter <- switch(eType,
                        "mom"=1/2*esMin^2,
                        "eGauss"=esMin^2,
                        "imom"=abs(esMin),
                        "eCauchy"=abs(esMin),
                        "grow"=abs(esMin))
  }

  if (eType=="grow") {
    if (alternative=="less")
      parameter <- -parameter
  }
  return(parameter)
}

#' Match the meanDiffMin of a savi z-test
#'
#' Based on the parameter, sigma, alternative and eType
#'
#' @inherit designSaviZ
#'
#' @returns the parameter, a numeric value
#' @export
#'
#' @examples
#' matchMeanDiffMinZFrom(parameter=0.4)
#'
matchEsMinWith <- function(parameter, analysisType=c("z", "t"),
                           sigma=1,
                           alternative=c("twoSided", "greater", "less"),
                           eType=c("mom", "eGauss", "imom", "eCauchy", "grow")) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  parameter <- abs(parameter)

  if (analysisType=="z") {
    esMin <- switch(eType,
                    "mom"=sqrt(2*parameter)*sigma,
                    "eGauss"=sqrt(parameter)*sigma,
                    "imom"=sqrt(parameter)*sigma,
                    "eCauchy"=parameter*sigma,
                    "grow"=abs(parameter))
  }

  if (analysisType=="t") {
    esMin <- switch(eType,
                    "mom"=sqrt(2*parameter),
                    "eGauss"=sqrt(parameter),
                    "imom"=abs(parameter),
                    "eCauchy"=abs(parameter),
                    "grow"=abs(parameter))
  }

  if (alternative=="less" && eType %in% c("grow", "imom", "eCauchy"))
    esMin <- -esMin

  return(esMin)
}

#' Match the parameter of a futility savi z-test
#'
#' Based on the esMinFutility, meanDiffMin, alternative and eType
#'
#' @inherit designSaviZ
#'
#' @returns the parameter, a numeric value
#' @export
#'
#' @examples
#' matchFutilityParameterZFrom(0.4)
matchFutilityParameterWith <- function(esMinFutility, esMin, esTrue) {

  if (!is.null(esMinFutility))
    return(abs(esMinFutility))

  if (!is.null(esMin))
    return(abs(esMin))

  if (!is.null(esTrue))
    return(abs(esTrue))

}

matchPowerWith <- function(power, beta) {
  if (!is.null(power)) {
    if (!is.null(beta))
      warning("Both power and beta specified. Preference given to power")

    return(power)
  }

  if (!is.null(beta)) {
    power <- 1-beta
    return(power)
  }

  return(NULL)
}
