#' Safe Anytime-Valid One-Way Anova
#'
#' To compute equality of all means in K groups. Still in beta stage.
#'
#' @param formula a formula of the form lhs ~ rhs where lhs is a
#' numeric variable giving the data values and rhs a factor with at
#' least two levels/groups.
#' @param data a data frame (or similar: see model.frame()) containing
#' the variables in the formula. By default the variables are taken
#' from environment(formula).
#' @param subset an optional vector specifying a subset of observations
#' to be used.
#' @param na.action a function which indicates what should happen when
#' the data contain NAs. Defaults to getOption("na.action").
#' @param parameter numeric > 0 representing the scaling the in
#' the multivatiate Gaussian prior
#' @param sigma numeric > 0 if a population standard deviation is assumed
#' to be known. Default \code{NULL} to convey that the it is unknown.
#' @param ... further arguments to be passed to or from methods.
#'
#' @returns a list
#' @export
#'
#' @examples
#' # Data
#' plantData <- data.frame(
#'   treatment = factor(rep(c("control", "fertiliserA", "fertiliserB"), each = 10)),
#'   height = c(
#'     c(12, 14, 11, 13, 12, 15, 13, 11, 14, 12), # Control group
#'     c(18, 20, 19, 22, 21, 20, 17, 19, 21, 23), # Fertiliser A group
#'     c(15, 17, 16, 15, 18, 14, 16, 17, 15, 16)  # Fertiliser B group
#'   )
#' )
#'
#' saviAnovaOneB(height ~ treatment, data = plantData, parameter=1)
saviAnovaOneB <- function(formula, data, subset,
                          na.action, parameter, sigma=NULL, ...) {

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  oneWay <- TRUE

  if (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L)
    if (formula[[3L]] == 1L)
      oneWay <- FALSE
  else
    stop("'formula' missing or incorrect")

  # Rethink this
  #
  # matchedCall <- match.call(expand.dots = FALSE)
  #
  # if (is.matrix(eval(matchedCall[["data"]], parent.frame())))
  #   matchedCall[["data"]] <- as.data.frame(data)
  #
  # # Note: Prepare calling stats::model.frame instead of saviTTest
  # #
  # matchedCall[[1L]] <- quote(stats::model.frame)
  #
  #
  #
  #
  # # Call: stats::model.frame
  # #
  # modelFrame <- eval(matchedCall,
  #                    parent.frame())

  modelFrame <- stats::model.frame(formula=formula, data=data)

  # Naming
  dataName <- paste(names(modelFrame),
                    collapse=" by ")

  names(modelFrame) <- NULL
  response <- attr(attr(modelFrame, "terms"),
                   "response")

  if (isTRUE(oneWay)) {
    groupingFactor <- factor(modelFrame[[-response]])

    K <- nlevels(groupingFactor)

    if (K < 2L)
      stop("grouping factor must have at least 2 levels")

    dataList <- split(modelFrame[[response]], groupingFactor)

    nList <- lapply(dataList, length)
    sumXList <- lapply(dataList, sum)

    # Note(Alexander): Potential for blow up
    #
    sumSquaresList <- lapply(dataList, function(x)sum(x^2))

    # Compute e-value here ----
    #
    anovaResult <- saviAnovaStat("sumX"=sumXList, "n"=nList, "sumSquares"=sumSquaresList,
                                 "sigma"=sigma, "parameter"=parameter)

    # Naming stuff ----
    #
    # if (length(anovaResult[["estimate"]]) == 2L) {
    #   names(anovaResult[["estimate"]]) <- paste("mean in group", levels(groupingFactor))
    #   names(anovaResult[["designObj"]][["h0"]]) <-
    #     paste("true difference in means between",
    #           paste("group", levels(groupingFactor), collapse = " and "))
    # }
  }

  anovaResult[["dataName"]] <- dataName
  return(anovaResult)
}

#' Safe Anytime-Valid One-Way Anova based on summary statistics
#'
#' To compute equality of all means in K groups. Still in beta stage.
#'
#' @inheritParams saviAnovaOneB
#' @param sumX vector of numerics representing the sum of the observations
#' of the K groups
#' @param n vector of integers representing the sample sizes of the K groups
#' @param sumSquares vector of numerics representing the sum of the
#' square of the observations of the K groups
#'
#' @returns a list with the e-value
#' @export
#'
#' @examples
#' saviAnovaStat(1:6, 1:6, (1:6)^2, parameter=1)
saviAnovaStat <- function(sumX, n, sumSquares=NULL,
                          parameter, sigma=NULL) {
  if (is.list(sumX)) sumX <- unlist(sumX)
  if (is.list(n)) n <- unlist(n)

  K <- length(sumX)

  if (length(n)!=K)
    stop("There are K=", K, " cumulative sums, but only ",
         length(n), " number of sample sizes.")

  nPlus <- sum(n)
  qVec <- n/nPlus

  QBig <- diag(qVec) - qVec %*% t(qVec)
  Q <- QBig[(2:K), (2:K)]

  J <- matrix(1/2, nrow=K-1, ncol=K-1)
  diag(J) <- 1
  Id <- diag(K-1)
  N <- nPlus * Q

  xBar <- sumX/n

  xBarDiff <- xBar-xBar[1]
  xBarDiff <- xBarDiff[2:K]

  g <- parameter

  rescaleMatrix <- Id+g*(J %*% N)

  logDet <- -1/2*log(det(rescaleMatrix))

  if (is.null(sigma)) {
    if (is.list(sumSquares)) sumSquares <- unlist(sumSquares)

    ssWithin <- sum(sumSquares - n*xBar^2)

    altScaleMatrix <- N %*% solve(rescaleMatrix)

    logEValue <- logDet+(nPlus-1)/2*(
      log(ssWithin+t(xBarDiff) %*% N %*% xBarDiff)-
        log(ssWithin+t(xBarDiff) %*% altScaleMatrix %*% xBarDiff)
    )
  } else if (!is.null(sigma)) {
    if (sigma <= 0)
      stop("Given sigma non-positive.")

    xBarDiff <- xBarDiff/sigma

    altScaleMatrix <- t(N) %*% solve(rescaleMatrix) %*% J %*% N

    g/2*xBarDiff+t(xBarDiff)

    logEValue <- logDet+g/2*t(xBarDiff) %*% altScaleMatrix %*% xBarDiff
  }

  return(list("eValue"=exp(logEValue)))
}
