# Try helper functions -----

#' Tries to Evaluate an Expression and Fails with \code{NA}
#'
#' The evaluation fails with \code{NA} by default, but it is also able to fail with other values.
#'
#' @param expr Expression to be evaluated.
#' @param value Return value if there is an error, default is \code{NA_real_}.
#'
#' @return Returns the evaluation of the expression, or \code{value} if it doesn't work out.
tryOrFailWithNA <- function(expr, value=NA_real_) {
  tryCatch(
    error=function(cnd) value,
    expr
  )
}

#' Checks Whether a Vector of Object Inherits from the Class 'try-error'
#'
#' Checks whether any of the provided objects contains a try error.
#'
#' @param ... objects that need testing.
#'
#' @return Returns \code{TRUE} if there's some object that's a try-error, \code{FALSE} when all objects are
#' not try-errors.
#'
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- try(integrate(exp, -Inf, Inf))
#' isTryError(x, y)
#' isTryError(x, y, z)
isTryError <- function(...) {
  obj <- list(...)
  tryErrorFunc <- function(x){inherits(x, "try-error")}
  result <- purrr::some(obj, .p=tryErrorFunc)
  return(result)
}

#' Helper function: Get all arguments as entered by the user
#'
#' @return a list of variable names of class "call" that can be changed into names
getArgs <- function() {
  as.list(match.call(definition = sys.function(-1),
                     call = sys.call(-1)))[-1]
}


#' Helper function: Get all names as entered by the user
#'
#' @param list list from which the element needs retrieving
#' @param name character string, name of the item that need retrieving
#'
#' @return returns a character string
extractNameFromArgs <- function(list, name) {
  result <- list[[name]]

  if (inherits(result, "call"))
    result <- as.character(as.expression(result))

  return(result)
}


# Check Consistency function --------

#' Checks consistency between the sided of the hypothesis and the  minimal clinically relevant effect size
#' or safe test defining parameter. Throws an error if the one-sided hypothesis is incongruent with the
#'
#' @inheritParams designSafeZ
#' @param paramToCheck numeric. Either a named safe test defining parameter such as phiS, or thetaS, or a
#' minimal clinically relevant effect size called with a non-null esMinName name
#' @param esMinName provides the name of the effect size. Either "meanDiffMin" for the z-test, "deltaMin" for
#' the t-test, or "hrMin" for the logrank test
#' @param paramDomain Domain of the paramToCheck, typically, positiveNumbers. Default \code{NULL}
#'
#' @return paramToCheck after checking, perhaps with a change in sign
checkAndReturnsEsMinParameterSide <- function(
    paramToCheck, alternative=c("twoSided", "greater", "less"),
    esMinName=c("noName", "meanDiffMin", "phiS",
                "deltaMin", "deltaS",
                "hrMin", "thetaS", "deltaTrue",
                "g", "kappaG"), paramDomain=NULL) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  paramDomain <- match.arg(paramDomain)
  esMinName <- match.arg(esMinName)

  if (alternative == "twoSided") {
    if (esMinName %in% c("meanDiffMin", "deltaMin", "deltaTrue"))
      return(abs(paramToCheck))

    return(paramToCheck)
  }

  if (esMinName=="noName")
    paramName <- NULL
  else
    paramName <- esMinName

  error <- NULL

  if (is.null(paramName)) {
    paramName <- "the safe test defining parameter"
    hypParamName <- "test relevant parameter"
    paramDomain <- "unknown"
  } else if (paramName=="phiS" || esMinName=="meanDiffMin") {
    hypParamName <- "meanDiff"
    paramDomain <- "realNumbers"
  } else if (paramName=="deltaS" || esMinName=="deltaMin"  || esMinName=="deltaTrue") {
    hypParamName <- "delta"
    paramDomain <- "realNumbers"
  } else if (paramName=="thetaS" || esMinName=="hrMin") {
    hypParamName <- "theta"
    paramDomain <- "positiveNumbers"

    error <- if (paramToCheck < 0) "thetaS and hrMin must be positive"
  } else if (paramName=="g") {
    hypParamName <- "g"
    paramDomain <- "positiveNumbers"

    error <- if (paramToCheck < 0) "The parameter g must be positive"
  } else if (paramName=="kappaG") {
    hypParamName <- "kappaG"
    paramDomain <- "positiveNumbers"

    error <- if (paramToCheck < 0) "The parameter kappaG must be positive"
  } else {
    hypParamName <- "testRelevantParameter"
  }

  if (!is.null(error))
    stop(error)

  if (paramDomain=="unknown") {
    nullValue <- "nullValue"

    if (alternative=="greater" && paramToCheck < 0) {
      warning('The safe test defining parameter is incongruent with alternative "greater". ',
              "This safe test parameter is made positive to compare H+: ",
              "test-relevant parameter > 0 against H0 : test-relevant parameter = 0")
      paramToCheck <- -paramToCheck
    }

    if (alternative=="less" && paramToCheck > 0) {
      warning('The safe test defining parameter is incongruent with alternative "less". ',
              "This safe test parameter is made positive to compare H-: ",
              "test-relevant parameter < 0 against H0 : test-relevant parameter = 0")
      paramToCheck <- -paramToCheck
    }

  } else if (paramDomain=="realNumbers") {
    nullValue <- 0

    if (alternative=="greater" && paramToCheck < 0) {
      warning(paramName, ' incongruent with alternative "greater". ',
              paramName, " set to -", paramName, " > 0 in order to compare H+: ",
              hypParamName, " > 0 against H0 : ", hypParamName, " = 0")
      paramToCheck <- -paramToCheck
    }

    if (alternative=="less" && paramToCheck > 0) {
      warning(paramName, ' incongruent with alternative "greater". ',
              paramName, " set to -", paramName, " < 0 in order to compare H-: ",
              hypParamName, " < 0 against H0 : ", hypParamName, " = 0")
      paramToCheck <- -paramToCheck
    }
  } else if (paramDomain=="positiveNumbers") {
    if (alternative=="greater" && paramToCheck < 1) {
      warning(paramName, ' incongruent with alternative "greater". ',
              paramName, " set to 1/", paramName, " > 1 in order to compare H+: ",
              hypParamName, " > 1 against H0 : ", hypParamName, " = 1")

      paramToCheck <- 1/paramToCheck
    }

    if (alternative=="less" && paramToCheck > 1) {
      warning(paramName, ' incongruent with alternative "greater". ',
              paramName, " set to 1/", paramName, " < 1 in order to compare H-: ",
              hypParamName, " < 1 against H0 : ", hypParamName, " = 1")

      paramToCheck <- 1/paramToCheck
    }
  }

  return(paramToCheck)
}

#' Check consistency between nPlan and the testType for one and two-sample z and t-tests
#'
#' @inheritParams designSafeZ
#'
#' @return nPlan a vector of sample sizes of length 1 or 2
#'
checkAndReturnsNPlan <- function(nPlan, ratio=1, testType=c("oneSample", "paired", "twoSample")) {
  if (testType=="twoSample" && length(nPlan)==1) {
    nPlan <- c(nPlan, ratio*nPlan)
    warning('testType=="twoSample" specified, but nPlan[2] not provided. nPlan[2] = ratio*nPlan[1], that is, ',
            nPlan[2], '.')
  } else if (testType=="paired" && length(nPlan)==1) {
    nPlan <- c(nPlan, nPlan)
    warning('testType=="paired" specified, but nPlan[2] not provided. nPlan[2] set to nPlan[1].')
  } else if (testType=="oneSample" && length(nPlan)==2) {
    nPlan <- nPlan[1]
    warning('testType=="oneSample" specified, but two nPlan[2] provided, which is ignored.')
  }
  return(nPlan)
}


# Plot helper -----
#' Sets 'safestats' Plot Options and Returns the Current Plot Options.
#'
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns a list with the user specified plot options.
#'
#' @export
#'
#' @examples
#' oldPar <- setSafeStatsPlotOptionsAndReturnOldOnes()
#' graphics::plot(1:10, 1:10)
#' setPar <- graphics::par(oldPar)
setSafeStatsPlotOptionsAndReturnOldOnes <- function(...) {
  oldPar <- graphics::par(no.readonly = TRUE)
  graphics::par(cex.main=1.5, mar=c(5, 6, 4, 4)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
                font.lab=2, cex.axis=1.3, bty="n", las=1, ...)
  return(oldPar)
}

#' Helper function to check whether arguments are specified in a function at a higher level and already
#' provided in the design object.
#'
#' @param designObj an object of class "safeDesign".
#' @param ... arguments that need checking.
#'
#' @return Returns nothing only used for its side-effects to produces warnings if needed.
#'
#' @export
#'
#' @examples
#' designObj <- designSafeZ(0.4)
#'
#' checkDoubleArgumentsDesignObject(designObj, "alpha"=NULL, alternative=NULL)
#' # Throws a warning
#' checkDoubleArgumentsDesignObject(designObj, "alpha"=0.4, alternative="d")
checkDoubleArgumentsDesignObject <- function(designObj, ...) {

  argsToCheck <- list(...)

  for (neem in names(argsToCheck)) {
    argument <- argsToCheck[[neem]]

    if (!is.null(argument) && argument != designObj[[neem]])
      warning("Both a design object and '", neem, "' provided. The '", neem, "' specified by the design ",
              "object is used for the test, and the provided '", neem, "' is ignored.")

  }
}
