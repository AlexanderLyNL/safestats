#' Title
#'
#' @param object
#' @param designObj
#' @param pilot
#' @param alpha
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
safeLogrankTest <- function(object, designObj=NULL, pilot=FALSE, alpha=NULL, ...) {
  if (isFALSE(pilot) && is.null(designObj))
    stop("Please provide a safe logrank design object, or run the function with pilot=TRUE. ",
         "A design object results can be acquired from designSafeLogrank().")

  if (isFALSE(pilot) && designObj[["testType"]] != "logrank")
    stop("The design is constructed for logrank tests, please use designSafeLogrank() for this.")

  if (!is.null(designObj) && !is.null(alpha)) {
    warning("Both designObj and alpha given. The alpha in designObj is used, and the given alpha is ignored")
    alpha <- NULL
  }

  UseMethod("safeLogrankTest")
}

#' Title
#'
#' @param formula
#' @param data
#' @param subset
#' @param weights
#' @param designObj
#' @param pilot
#' @param alpha
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
safeLogrankTest.formula <- function (formula, data=list(), subset=NULL, weights=NULL, designObj=NULL,
                                     pilot=FALSE, alpha=NULL, ...) {
  logrankObj <- try(
    coin:::logrank_test.formula("formula"=formula, "data"=data, "subset"=subset, "weights"=weights, ...)
  )

  safeLogrankTestCore("logrankObj"=logrankObj, "designObj"=designObj, "pilot"=pilot, "alpha"=alpha)
}


#' Title
#'
#' @param object
#' @param ties.method
#' @param type
#' @param rho
#' @param gamma
#' @param designObj
#' @param pilot
#' @param alpha
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
safeLogrankTest.IndependenceProblem <- function(object, ties.method=c("mid-ranks", "Hothorn-Lausen", "average-scores"),
                                                type=c("logrank", "Gehan-Breslow", "Tarone-Ware", "Prentice",
                                                       "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
                                                       "Fleming-Harrington", "Gaugler-Kim-Liao", "Self"),
                                                rho=NULL, gamma=NULL, designObj=NULL, pilot=FALSE, alpha=NULL, ...) {
  logrankObj <- coin:::logrank_test.IndependenceProblem("object"=object, "ties.method"=ties.method, "type"=type,
                                                        "rho"=rho, "gamma"=gamma, ...)

  safeLogrankTestCore("logrankObj"=logrankObj, "designObj"=designObj, "pilot"=pilot, "alpha"=alpha)
}


#' Title
#'
#' @param logrankObj
#' @param designObj
#' @param pilot
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
safeLogrankTestCore <- function(logrankObj, designObj=NULL, pilot=FALSE, alpha=NULL) {
  if (!inherits(logrankObj, "ScalarIndependenceTest"))
    stop("Unable to compute the relevant test statistic. Debug in coin::logrank_test()")

  groupLevels <- levels(logrankObj@statistic@x[[groupLabel]])

  alternative <- logrankObj@statistic@alternative

  if (alternative != "two.sided")
    stop("One-sided tests not yet implemented")

  if (length(groupLevels) > 2)
    stop("K-sample log rank test not yet implemented")

  zStat <- unname(logrankObj@statistic@standardizedlinearstatistic)
  nEvents <- sum(logrankObj@statistic@ytrans < 0)
  dataName <- names(logrankObj@statistic@y)

  groupLabel <- names(logrankObj@statistic@x)

  result <- list("statistic"=zStat, "n1"=nEvents, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
                 "theta0"=1, "alternative"=alternative, "testType"="logrank", "dataName"=dataName,
                 "groupLabel"=groupLabel, "groupLevels"=groupLevels)
  class(result) <- "safeLogrank"

  if (isTRUE(pilot)) {

    if (is.null(alpha))
      alpha <- 0.05

    designObj <- designSafeLogrank("nPlan"=nEvents, "alpha"=alpha)
  }

  sValue <- safeZTestStat("z"=zStat, "phiS"=designObj[["thetaS"]], "n1"=nEvents,
                          "n2"=NULL, "alternative"="two.sided",
                          "paired"=FALSE, "sigma"=1)

  result[["sValue"]] <- sValue
  result[["designObj"]] <- designObj
  result[["logrankObj"]] <- logrankObj

  return(result)
}

#' Designs a safe logrank rest
#'
#' Designs a safe experiment for a prespecified tolerable type I error combined with (1) a targetted number
#' of events nPlan, or (2) a tolerable type II beta and a minimal clinically relevant hazard ratio thetaMin.
#' Outputs a list that includes the thetaS that defines the safe test.
#'
#' @inheritParams designSafeT
#' @param nPlan numeric > 0, targetted number of events
#' @param thetaMin numeric that defines the minimal relevant hazard ratio, the smallest hazard ratio that we want to
#' detect
#' @param zApprox logical, default TRUE to use the asymptotic normality
#'
#' @return Returns a safeDesign object that includes:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size either (1) specified by the user, or (2) computed based on beta and thetaMin}
#'   \item{thetaS}{the thetaS that defines the safe test}
#'   \item{thetaMin}{the minimally clinically relevant hazard ratio specified by the user}
#'   \item{alpha}{the tolerable type I error provided by the user}
#'   \item{beta}{the tolerable type II error provided by the user}
#'   \item{tol}{the step size between lowTheta and highTheta provided by the user}
#'   \item{alternative}{any of "two.sided", "greater", "less" provided by the user}
#'   \item{testType}{"logrank"}
#'   \item{ratio}{default is 1}
#'   \item{pilot}{FALSE to indicate that the design is not a pilot study}
#'   \item{call}{the expression with which this function is called}
#' }
#'
#' @export
#'
#' @examples
designSafeLogrank <- function(nPlan=50, alpha=0.05, beta=NULL, alternative="two.sided",
                              ratio=1, thetaMin=NULL, zApprox=TRUE, tol=1e-5, ...) {
  stopifnot(0 < alpha, alpha < 1)

  alternative <- match.arg(alternative)

  result <- list("nPlan"=nPlan, "thetaS"=NA, "thetaMin"=thetaMin, "alpha"=alpha, "beta"=beta,
                 # "tol"=tol, "lowN"=NULL, "highN"=NULL,
                 "alternative"=alternative, "testType"="logrank", "ratio"=ratio, "pilot"=FALSE,
                 "call"=sys.call())
  class(result) <- "safeLogrankDesign"

  if (is.null(beta) && is.null(thetaMin)) {
    if (zApprox) {
      safeZObj <- try(designPilotSafeZ("n1"=nPlan, "alpha"=alpha, "alternative"=alternative))
    } else {
      warning("Not yet implemented. Used z approximation instead.")

      safeZObj <- try(designPilotSafeZ("n1"=nPlan, "alpha"=alpha, "alternative"=alternative))
    }

    if (isTryError(safeZObj))
      stop("Design based on the given alpha and nPlan")

    result[["thetaS"]] <- safeZObj[["phiS"]]
  } else {

    if (is.null(nPlan))
      stop("Can't design without targetted nPlan number of events.")

    warning("Not yet implemented. Used z approximation instead.")

    safeZObj <- try(designPilotSafeZ("n1"=nPlan, "alpha"=alpha, "alternative"=alternative))
    result[["thetaS"]] <- safeZObj[["phiS"]]
  }

  return(result)
}

designObj <- designSafeLogrank(nPlan=89, alpha=0.05)
class(designObj)

designObj

bob

zDesign <- designSafeZ(0.5)

print.safeZDesign <- function(x, ...) {
  analysisName <- getNameTestType(testType = x[["testType"]])

  cat("\n")
  cat(paste("       ", analysisName, "\n"))
  cat("\n")

  if (isFALSE(x[["pilot"]])) {
    if (is.null(x[["n2Plan"]])) {
      cat("is powered for an experiment with a sample size of: ")
      cat("\n")
      cat(paste("    n1Plan =", x[["n1Plan"]]))
      cat("\n")
    } else {
      cat("Powered for an experiment with sample sizes: ")
      cat("\n")
      cat(paste("    n1Plan =", x[["n1Plan"]], "and n2Plan =", x[["n2Plan"]]))
      cat("\n")
    }
    cat("to find an effect size of at least: ")
    cat("\n")
    cat("    deltaMin =", round5(x[["phiMin"]]))
    cat("\n")
    cat("\n")

    cat("with:")
    cat("\n")
    cat("    power = ", 1 - x[["beta"]], " (thus, beta = ", x[["beta"]], ")", sep="")
    cat("\n")

    cat("under the alternative:")
    cat("\n")
    cat("   ", getNameAlternative(x[["alternative"]], x[["testType"]]))
    cat("\n")
    cat("\n")

    cat("Based on the decision rule S > 1/alpha:")
    cat("\n")
    cat("    S > ", round5(1/x[["alpha"]]), sep="")
    cat("\n")

    cat("which occurs with chance less than:")
    cat("\n")
    cat("    alpha =", x[["alpha"]])
    cat("\n")

    cat("under iid normally distributed data and the null hypothesis:")
    cat("\n")
    cat("    mu =", x[["mu0"]])
  } else {
    cat("The experiment is not planned.")
    cat("\n")
    cat("This design object only valid for experiments with:")
    cat("\n")

    if (is.null(x[["n2Plan"]])) {
      cat("    n1 =", x[["n1Plan"]])
      cat("\n")
    } else {
      cat("    n1 =", x[["n1Plan"]], "and n2 =", x[["n2Plan"]])
      cat("\n")
    }
  }
}
