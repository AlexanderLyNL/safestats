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
    stop("The provided logrankObj is not of the right type derived from coin::logrank_test()")

  groupLevels <- levels(logrankObj@statistic@x[[groupLabel]])

  alternative <- logrankObj@statistic@alternative

  if (alternative != "two.sided")
    stop("One-sided tests not yet implemented")

  if (length(groupLevels) > 2)
    stop("K-sample log rank test not yet implemented")

  zStat <- unname(logrankObj@statistic@standardizedlinearstatistic)
  names(zStat) <- "z"

  nEvents <- sum(logrankObj@statistic@ytrans < 0)
  names(nEvents) <- "nEvents"

  dataName <- paste0(names(logrankObj@statistic@y), " by ",
                     names(logrankObj@statistic@x), " (",
                     paste(groupLevels, collapse=", "), ")")

  # result <- list("statistic"=NULL, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
  #                "alternative"=alternative, "testType"=NULL, "dataName"=NULL, "mu0"=mu0, "sigma"=sigma)

  result <- list("statistic"=zStat, "n1"=nEvents, "sValue"=NULL, "confInt"=NULL, "estimate"=NULL,
                 "theta0"=1, "alternative"=alternative, "testType"="logrank", "dataName"=dataName)
  class(result) <- "safeTest"

  if (isTRUE(pilot)) {

    if (is.null(alpha))
      alpha <- 0.05

    designObj <- designSafeLogrank("nPlan"=nEvents, "alpha"=alpha)
    designObj[["pilot"]] <- TRUE
  }

  sValue <- safeZTestStat("z"=zStat, "parameter"=designObj[["parameter"]], "n1"=nEvents,
                          "n2"=NULL, "alternative"="two.sided", "paired"=FALSE, "sigma"=1)

  result[["sValue"]] <- sValue
  result[["designObj"]] <- designObj
  result[["logrankObj"]] <- logrankObj

  return(result)
}

#' Designs a safe logrank rest
#'
#' Designs a safe experiment for a prespecified tolerable type I error combined with (1) a targetted number
#' of events nPlan, or (2) a tolerable type II beta and a minimal clinically relevant hazard ratio thetaMin.
#' Outputs a list that includes the parameter that defines the safe test.
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
#'   \item{parameter}{the parameter that defines the safe test}
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
designSafeLogrank <- function(alpha=0.05, beta=NULL, thetaMin=NULL, alternative=c("two.sided", "greater", "less"),
                              nPlan=NULL, ratio=1, zApprox=TRUE, tol=1e-5, ...) {
  stopifnot(0 < alpha, alpha < 1)

  alternative <- match.arg(alternative)

  if (alternative != "two.sided")
    warning("Currently, only the two.sided method is implemented, and therefore used.")

  names(nPlan) <- "nPlan"

  result <- list("nPlan"=nPlan, "parameter"=NULL, "esMin"=thetaMin, "alpha"=alpha, "beta"=beta,
                 "alternative"=alternative, "testType"="logrank", "ratio"=ratio,
                 "pilot"=FALSE, "call"=sys.call())
  class(result) <- "safeDesign"

  if (!is.null(beta))
    warning("Currently, only the method based on nPlan is implemented, and therefore used. beta ignored")

  if (!is.null(thetaMin))
    warning("Currently, only the method based on nPlan is implemented, and therefore used. thetaMin ignored")

  if (is.null(beta) && is.null(thetaMin)) {
    if (zApprox) {
      safeZObj <- try(designPilotSafeZ("alpha"=alpha, "nPlan"=nPlan, "alternative"=alternative))
    } else {
      warning("Currently, only the z-approximation method is implemented, and therefore used.")

      # TODO(Alexander): Replace by hypergeometric method
      #
      safeZObj <- try(designPilotSafeZ("alpha"=alpha, "nPlan"=nPlan, "alternative"=alternative))
    }

    if (isTryError(safeZObj))
      stop("Unable to design the given alpha and nPlan")

  } else {

    if (is.null(nPlan))
      stop("Can't design without targetted nPlan number of events.")

    warning("Currently, only the z-approximation method is implemented, and therefore used.")

    # TODO(Alexander): Replace by most likely a z approximation power analysis based on beta and thetaMin
    #
    safeZObj <- try(designPilotSafeZ("alpha"=alpha, "nPlan"=nPlan, "alternative"=alternative))
  }

  result[["parameter"]] <- safeZObj[["parameter"]]
  names(result[["parameter"]]) <- "thetaS"

  return(result)
}






