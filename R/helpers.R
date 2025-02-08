#' Add a savi reference
#'
#' @param key character string. firstAuthorYearFirstWord
#'
#' @return a character string
#' @export
#'
#' @examples
#' addCite(grunwald2024safe)
addCite <- function(..., breakLine=TRUE) {
  keys <- as.character(rlang::ensyms(...))

  # TODO(ALEXANDER) CHECK: utils::cite

  refList <- list(
    "grunwald2024safe"="Gr\u00fcnwald, P. D., de Heide, R., & Koolen, W. (2024). Safe testing. \\emph{Journal of the Royal Statistical Society. Series B (Methodological),} \\strong{86}(\\emph{5}), 1091–1128. (With discussions), https://doi.org/10.1093/jrsssb/qkae011.",
    "grunwald2024authors"="Gr\u00fcnwald, P. D., de Heide, R., & Koolen, W. (2024). Authors’ reply to the discussion of ‘Safe testing’. \\emph{Journal of the Royal Statistical Society. Series B (Methodological),} \\strong{86}(\\emph{5}), 1091–1128, https://doi.org/10.1093/jrsssb/qkae069.",
    "ly2024safe"="Ly, A, Boehm, Gr\u00fcnwald, P. D., Ramdas, A., & van Ravenzwaaij, D. (2024). Safe Anytime-Valid Inference: Practical maximally flexible sampling designs for experiments based on e-values. \\emph{PsyArXiv Preprint}, https://doi.org/10.31234/osf.io/h5vae.",
    "schure2024safe"="ter Schure, J., P\u00e9rez-Ortiz, M. F., Ly, A., & Gr\u00fcnwald, P. D. (2024). The Safe Logrank Test: Error control under continuous monitoring with unlimited horizon. \\emph{The New England Journal of Statistics in Data Science}, \\strong{2}(\\emph{2}), 190-214, https://doi.org/10.51387/24-NEJSDS65.",
    "perez2024estatistics"="P\u00e9rez-Ortiz, M. F., Lardy, T., de Heide, R., & Gr\u00fcnwald, P. D. (2024). E-statistics, group invariance and anytime valid testing. \\emph{The Annals of Statistics}, \\strong{52}(\\emph{4}), 1410-1432, http://dx.doi.org/10.1214/24-AOS2394.",
    "turner2024generic"="Turner, R., Ly, A., & Gr\u00fcnwald, P. D. (2024). Generic e-variables for exact sequential k-sample tests that allow for optional stopping. \\emph{Journal of Statistical Planning and inference} \\strong{230}, 106116, https://doi.org/10.1016/j.jspi.2023.106116.",
    "turner2023exact"="Turner, R., & Gr\u00fcnwald, P. D. (2024). Exact anytime-valid confidence intervals for contingency tables and beyond. \\emph{Statistics and Probability letters}, \\strong{198}, 109835, https://doi.org/10.1016/j.spl.2023.109835.",
    "ramdas2023game"="Ramdas, A, Gr\u00fcnwald, P. D., Vovk, V., & Shafer, G. (2023). Game-theoretic statistics and Safe Anytime-Valid Inference. \\emph{Statistical Science}, \\strong{38}(\\emph{4}), 576-597, https://doi.org/10.1214/23-STS894.",
    "wang2025anytime"="Wang, H., & Ramdas, A. (in press). Anytime-valid t-tests and confidence sequences for Gaussian means with unknown variance. \\emph{Sequential Analysis}, https://doi.org/10.48550/arXiv.2310.03722.",
    "schoenfeld1981asymptotic"="Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. \\emph{Biometrika}, \\strong{68}(\\emph{1}), 316-319, https://doi.org/10.2307/2335833."
  )

  refs <- refList[keys]

  refLength <- length(refs)

  res <- character()

  for (i in seq_along(refs)) {
    if (i == refLength && isFALSE(breakLine)) {
      res <- paste(res, refs[[i]])
    } else {
      res <- paste(res, refs[[i]], "</br> </br>")
    }
  }

  return(res)
}

.onAttach <-
  function(libname, pkgname) {
    packageCitation <- "Ly, A., Turner, R. J., Wang, Y., P\u00e9rez-Ortiz, M. F., Boehm, U., ter Schure, J., & Gr\u00fcnwald, P. D. (2024). safestats: Safe anytime-valid inference [Computer software manual]."
    message("Thank you for using safestats! We really appreciate your interest in our work.")
    message("To acknowledge our work, please feel free to cite the package as follows:")
    message(packageCitation)
  }

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
#' or savi test defining parameter. Throws an error if the one-sided hypothesis is incongruent with the
#'
#' @inheritParams designSaviZ
#' @param paramToCheck numeric. Either a named savi test defining parameter such as phiS, or thetaS, or a
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
    paramName <- "the savi test defining parameter"
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
      warning('The savi test defining parameter is incongruent with alternative "greater". ',
              "This savi test parameter is made positive to compare H+: ",
              "test-relevant parameter > 0 against H0 : test-relevant parameter = 0")
      paramToCheck <- -paramToCheck
    }

    if (alternative=="less" && paramToCheck > 0) {
      warning('The savi test defining parameter is incongruent with alternative "less". ',
              "This savi test parameter is made positive to compare H-: ",
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
#' @inheritParams designSaviZ
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
#' designObj <- designSaviZ(0.4)
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

#' A ceiling function with a tolerance at the 13th digit
#'
#' @param x numeric, that needs
#' @param digits integer, position of the digit to round of to
#'
#' @return integer
#' @export
#'
#' @examples
#' ceiling(27/21*21)
#' ceil(27/21*21)
ceil <- function(x, digits=13) {
  ceiling(round(x, digits=digits))
}


#' Helper function to create running intersections
#'
#' @param x vector of numeric, representing a sequence of upper or
#' lower bounds of a confidence sequence
#' @param upper logic, by default \code{TRUE} to construct a
#' running intersection for the upper bound of a sequence with
#' the minimum function. If \code{FALSE}, then use the maximum
#' function for the lower bound.
#'
#' @return a sequence of numerics representing the running minimum
#' or maximum
#
#' @export
#'
#' @examples
#'
#' makeRunningIntersection(c(6, -1, 3, 12))
makeRunningIntersection <- function(x, upper=TRUE) {
  m <- length(x)

  if (m < 2)
    stop("No running intersection for sequence of length less than 2")

  res <- numeric(m)

  if (isTRUE(upper)) {
    currentValue <- Inf
    comparisonFunction <- base::min
  } else {
    currentValue <- -Inf
    comparisonFunction <- base::max
  }


  for (i in seq_along(x)) {
    res[i] <- comparisonFunction(currentValue, x[i])
    currentValue <- res[i]
  }

  return(res)
}
