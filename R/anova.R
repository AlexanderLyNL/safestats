iets <- function (formula, data, subset, na.action, sigma=NULL, ...) {

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  oneWay <- TRUE

  if (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L)
    if (formula[[3L]] == 1L)
      oneWay <- FALSE
  else
    stop("'formula' missing or incorrect")

  matchedCall <- match.call(expand.dots = FALSE)

  if (is.matrix(eval(matchedCall[["data"]], parent.frame())))
    matchedCall[["data"]] <- as.data.frame(data)

  # Note: Prepare calling stats::model.frame instead of saviTTest
  #
  matchedCall[[1L]] <- quote(stats::model.frame)

  # Call: stats::model.frame
  #
  modelFrame <- eval(matchedCall,
                     parent.frame())

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
    xBarList <- lapply(dataList, mean)

    if (!is.null(sigma))

    varList <- if (is.null(sigma)) NULL else lapply(dataList, stats::var)

    browser()

    # Compute e-value here ----
    #
    anovaResult <- saviTTest("x"=dataList[[1L]], "y"=dataList[[2L]], ...)

    # Naming stuff ----
    #
    if (length(anovaResult[["estimate"]]) == 2L) {
      names(anovaResult[["estimate"]]) <- paste("mean in group", levels(groupingFactor))
      names(anovaResult[["designObj"]][["h0"]]) <-
        paste("true difference in means between",
              paste("group", levels(groupingFactor), collapse = " and "))
    }
  }

  anovaResult[["dataName"]] <- dataName
  return(anovaResult)
}
