# Testing fnts ----





#' Safe anytime-valid 2x2 test for propDiff
#' @noRd
savi2x2TestStatPropDiff <- function(ya, yb, na = NULL, nb = NULL,
                                    designObj = NULL, wantCi = TRUE) {
  # TODO: THESE ARGS CHECKING WILL BE DONE FINAL STEP
  # grow must have a propDiffMin
  # eBeta must not have a propDiffMin, i.e. propDiffMin == NULL and must be twoSided
  result <- constructSaviTestObj("Two Proportions")

  propDiffMin <- designObj[["esMin"]]
  alternative <- designObj[["alternative"]]
  betaParameter <- designObj[["betaParameter"]]
  alpha <- designObj[["alpha"]]

  nBlocks <- length(ya)
  if (is.null(na)) na <- designObj[["nPlan"]][["na"]]
  if (is.null(nb)) nb <- designObj[["nPlan"]][["nb"]]
  if (length(na) == 1L) na <- rep(na, nBlocks)
  if (length(nb) == 1L) nb <- rep(nb, nBlocks)

  # Compute: eValueVec ----
  if (is.null(propDiffMin)) {
    logEValueVec <- computeEValueVecPropDiff(ya, yb, na, nb, betaParameter)
  } else if (alternative == "greater") {
    logEValueVec <- computeEValueVecPropDiff(ya, yb, na, nb, betaParameter,
      propDiff = propDiffMin
    )
  } else {
    # twoSided and grow: equal-weight mixture at +propDiffMin and -propDiffMin
    logEPlus <- computeEValueVecPropDiff(ya, yb, na, nb, betaParameter,
      propDiff = propDiffMin
    )
    logEMinus <- computeEValueVecPropDiff(ya, yb, na, nb, betaParameter,
      propDiff = -propDiffMin
    )
    logEValueVec <- pmax(logEPlus, logEMinus) +
      log1p(exp(-abs(logEPlus - logEMinus))) - log(2)
  }

  eValueVec <- exp(logEValueVec)

  # Compute: confSeqMatrix ----
  if (wantCi) {
    confSetRuns <- computeConfidenceInterval2x2PropDiff(
      ya = ya, yb = yb, na = na, nb = nb,
      betaParameter = betaParameter, alpha = alpha,
      runningIntersection = designObj[["runningIntersection"]]
    )

    # One row per block, as for the other tests: the outermost bounds of that
    # block's union, NA when every candidate is rejected. The hull contains
    # the union, so coverage is kept; only the last block is kept exact.
    block <- factor(rownames(confSetRuns), levels = seq_len(nBlocks))
    confSeqMatrix <- cbind(
      "lowerBound" = as.vector(tapply(confSetRuns[, "lowerBound"], block, min)),
      "upperBound" = as.vector(tapply(confSetRuns[, "upperBound"], block, max))
    )
    lastBlock <- rownames(confSetRuns) == nBlocks

    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSetRuns[lastBlock, c("lowerBound", "upperBound"),
      drop = FALSE
    ]
    result[["ciValue"]] <- 1 - alpha
  }

  # Fill: Result ----
  result[["estimate"]] <- c("thetaA" = sum(ya) / sum(na), "thetaB" = sum(yb) / sum(nb))
  # x-axis of plot.saviTest(): the block index.
  result[["n1Vec"]] <- seq_len(nBlocks)
  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["n"]] <- c("na" = sum(na), "nb" = sum(nb), "nBlocks" = nBlocks)
  result[["betaPrior"]] <- list(
    "betaA1" = betaParameter[["betaA1"]] + sum(ya),
    "betaA2" = betaParameter[["betaA2"]] + sum(na) - sum(ya),
    "betaB1" = betaParameter[["betaB1"]] + sum(yb),
    "betaB2" = betaParameter[["betaB2"]] + sum(nb) - sum(yb)
  )
  result[["designObj"]] <- designObj
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- alternative
  result[["h0"]] <- designObj[["h0"]]
  result[["dataName"]] <- paste(
    deparse1(substitute(ya)), "and",
    deparse1(substitute(yb))
  )
  result[["call"]] <- sys.call()

  return(result)
}

#' Safe anytime-valid 2x2 test for logOdds
#' @noRd
savi2x2TestStatLogOdds <- function(ya, yb, na = NULL, nb = NULL,
                                   designObj = NULL, wantCi = TRUE) {
  # TODO: THESE ARGS CHECKING WILL BE DONE FINAL STEP
  # no restricted alternative for logOdds yet: eGauss or grow, twoSided only
  result <- constructSaviTestObj("Two Proportions")

  betaParameter <- designObj[["betaParameter"]]
  alpha <- designObj[["alpha"]]

  nBlocks <- length(ya)
  if (is.null(na)) na <- designObj[["nPlan"]][["na"]]
  if (is.null(nb)) nb <- designObj[["nPlan"]][["nb"]]
  if (length(na) == 1L) na <- rep(na, nBlocks)
  if (length(nb) == 1L) nb <- rep(nb, nBlocks)

  # Compute: eValueVec ----
  eValueVec <- exp(computeEValueVecLogOdds(ya, yb, na, nb, betaParameter))

  # Compute: confSeqMatrix ----
  if (wantCi) {
    confSetRuns <- computeConfidenceInterval2x2LogOdds(
      ya = ya, yb = yb, na = na, nb = nb,
      betaParameter = betaParameter, alpha = alpha,
      runningIntersection = designObj[["runningIntersection"]]
    )

    block <- factor(rownames(confSetRuns), levels = seq_len(nBlocks))
    confSeqMatrix <- cbind(
      "lowerBound" = as.vector(tapply(confSetRuns[, "lowerBound"], block, min)),
      "upperBound" = as.vector(tapply(confSetRuns[, "upperBound"], block, max))
    )
    lastBlock <- rownames(confSetRuns) == nBlocks

    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSetRuns[lastBlock, c("lowerBound", "upperBound"),
      drop = FALSE
    ]
    result[["ciValue"]] <- 1 - alpha
  }

  # Fill: Result ----
  result[["estimate"]] <- c("thetaA" = sum(ya) / sum(na), "thetaB" = sum(yb) / sum(nb))
  # x-axis of plot.saviTest(): the block index.
  result[["n1Vec"]] <- seq_len(nBlocks)
  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["n"]] <- c("na" = sum(na), "nb" = sum(nb), "nBlocks" = nBlocks)
  result[["betaPrior"]] <- list(
    "betaA1" = betaParameter[["betaA1"]] + sum(ya),
    "betaA2" = betaParameter[["betaA2"]] + sum(na) - sum(ya),
    "betaB1" = betaParameter[["betaB1"]] + sum(yb),
    "betaB2" = betaParameter[["betaB2"]] + sum(nb) - sum(yb)
  )
  result[["designObj"]] <- designObj
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- designObj[["alternative"]]
  result[["h0"]] <- designObj[["h0"]]
  result[["dataName"]] <- paste(
    deparse1(substitute(ya)), "and",
    deparse1(substitute(yb))
  )
  result[["call"]] <- sys.call()

  return(result)
}

# Design fnts ----
# TODO: add "less" when direction is clean in `alternative`
#
#' Design a safe anytime-valid 2x2 test
#' @noRd
designSavi2x2 <- function(
  na, nb,
  propDiffMin = NULL, logOddsMin = NULL,
  alpha = 0.05, power = NULL, h0 = 0,
  alternative = c("twoSided", "greater"),
  eType = c("eBeta", "grow", "eGauss"),
  betaParameter = NULL,
  runningIntersection = NULL
) {
  # TODO: THESE ARGS CHECKING WILL BE DONE FINAL STEP
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  result <- constructSaviDesignObj("Two Proportions")

  if (!is.null(betaParameter)) {
    requiredNames <- names(result[["betaParameter"]])
    result[["betaParameter"]] <- betaParameter[requiredNames]
  }

  # propDiff: eType eBeta or grow. logOdds: eType eGauss or grow. Only one of
  # propDiffMin, logOddsMin is set, matching the effect eType picks out.
  result[["esMin"]] <- if (!is.null(propDiffMin)) propDiffMin else logOddsMin
  result[["parameter"]] <- c(
    "Beta hyperparameters" =
      paste(unlist(result[["betaParameter"]]), collapse = " ")
  )
  result[["eType"]] <- eType
  result[["runningIntersection"]] <- runningIntersection
  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["h0"]] <- c("propDiff" = h0)
  # TODO: add stopping time simulation
  result[["nPlan"]] <- list("na" = na, "nb" = nb)
  result[["testType"]] <- "2x2"
  result[["call"]] <- sys.call()
  result[["timeStamp"]] <- Sys.time()

  return(result)
}


# Confidence Interval ----

#' Anytime-valid confidence sequence for the proportion difference
#' @noRd
computeConfidenceInterval2x2PropDiff <- function(ya, yb, na, nb,
                                                 betaParameter,
                                                 alpha, precision = 100,
                                                 runningIntersection = TRUE) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)

  propDiffGrid <- seq(-1, 1, length.out = precision + 2)[-c(1, precision + 2)]
  logEValues <- numeric(precision)
  inSet <- rep(TRUE, precision)

  # ncol = 2: only lowerBound, upperBound are data columns. Block is the
  # rowname, since a block may contribute zero, one or several rows.
  confSeqMatrix <- matrix(numeric(0),
    ncol = 2,
    dimnames = list(NULL, c("lowerBound", "upperBound"))
  )

  for (i in seq_len(nBlocks)) {
    # Under the running intersection a rejected candidate never returns, so
    # its e-process is not advanced; otherwise every candidate is followed.
    activeCandidates <- if (runningIntersection) which(inSet) else seq_len(precision)

    for (j in activeCandidates) {
      propDiff <- propDiffGrid[j]
      nullThetaA <- solveRIPr2x2PropDiff(
        thetaA = thetas[["thetaA"]][i], thetaB = thetas[["thetaB"]][i],
        na = na[i], nb = nb[i], propDiff = propDiff
      )

      logEValues[j] <- logEValues[j] + logLikelihoodRatioMultiBern(
        ya = ya[i], yb = yb[i], na = na[i], nb = nb[i],
        numeratorThetaA = thetas[["thetaA"]][i],
        numeratorThetaB = thetas[["thetaB"]][i],
        denominatorThetaA = nullThetaA,
        denominatorThetaB = nullThetaA + propDiff,
        log = TRUE
      )
    }

    notRejected <- logEValues < log(1 / alpha)
    inSet <- if (runningIntersection) inSet & notRejected else notRejected

    # Split the non-rejected candidates into runs of neighbours on the grid;
    # each run is one interval of the union.
    runs <- rle(inSet)
    runEnds <- cumsum(runs[["lengths"]])[runs[["values"]]]
    runStarts <- runEnds - runs[["lengths"]][runs[["values"]]] + 1

    newRows <- cbind(
      "lowerBound" = propDiffGrid[runStarts],
      "upperBound" = propDiffGrid[runEnds]
    )
    rownames(newRows) <- rep(i, length(runStarts))

    confSeqMatrix <- rbind(confSeqMatrix, newRows)
  }

  return(confSeqMatrix)
}

#' Anytime-valid confidence sequence for the log odds ratio
#' @noRd
computeConfidenceInterval2x2LogOdds <- function(ya, yb, na, nb,
                                              betaParameter,
                                              alpha, precision = 100,
                                              logOddsBound = 40,
                                              runningIntersection = TRUE) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  # Predictable plug-in alternative: B minus A on the logit scale.
  # TODO: to be decided here
  plugInLogOdds <- stats::qlogis(thetas[["thetaB"]]) -
    stats::qlogis(thetas[["thetaA"]])

  logOddsGrid <- seq(-logOddsBound, logOddsBound,
    length.out = precision + 2
  )[-c(1, precision + 2)]
  logEValues <- numeric(precision)
  inSet <- rep(TRUE, precision)

  # ncol = 2: only lowerBound, upperBound are data columns. Block is the
  # rowname, since a block may contribute zero, one or several rows.
  confSeqMatrix <- matrix(numeric(0),
    ncol = 2,
    dimnames = list(NULL, c("lowerBound", "upperBound"))
  )

  for (i in seq_len(nBlocks)) {
    # Under the running intersection a rejected candidate never returns, so
    # its e-process is not advanced; otherwise every candidate is followed.
    activeCandidates <- if (runningIntersection) which(inSet) else seq_len(precision)

    for (j in activeCandidates) {
      logEValues[j] <- logEValues[j] + logLikelihoodRatioFNCH(
        ya = ya[i], yb = yb[i], na = na[i], nb = nb[i],
        logOdds = plugInLogOdds[i], nullLogOdds = logOddsGrid[j], log = TRUE
      )
    }

    notRejected <- logEValues < log(1 / alpha)
    inSet <- if (runningIntersection) inSet & notRejected else notRejected

    runs <- rle(inSet)
    runEnds <- cumsum(runs[["lengths"]])[runs[["values"]]]
    runStarts <- runEnds - runs[["lengths"]][runs[["values"]]] + 1

    newRows <- cbind(
      "lowerBound" = logOddsGrid[runStarts],
      "upperBound" = logOddsGrid[runEnds]
    )
    rownames(newRows) <- rep(i, length(runStarts))

    confSeqMatrix <- rbind(confSeqMatrix, newRows)
  }

  return(confSeqMatrix)
}

# Helpers ----

## propDiff ----
# TODO: check if this works for ya yb na nb of vector at the same length!
logLikelihoodRatioMultiBern <- function(ya, yb, na, nb,
                                        numeratorThetaA, numeratorThetaB,
                                        denominatorThetaA, denominatorThetaB,
                                        log = FALSE, ...) {
  successesA <- ya * (log(numeratorThetaA) - log(denominatorThetaA))
  successesB <- yb * (log(numeratorThetaB) - log(denominatorThetaB))

  failuresA <- (na - ya) * (log1p(-numeratorThetaA) - log1p(-denominatorThetaA))
  failuresB <- (nb - yb) * (log1p(-numeratorThetaB) - log1p(-denominatorThetaB))

  logEValueVec <- cumsum(successesA + failuresA + successesB + failuresB)

  if (log) logEValueVec else exp(logEValueVec)
}

# Cumulative log e-value vector against thetaA = thetaB; length nBlocks.
computeEValueVecPropDiff <- function(ya, yb, na, nb, betaParameter,
                                     propDiff = NULL) {
  thetas <- if (is.null(propDiff)) {
    predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  } else {
    predictiveThetas2x2PropDiff(ya, yb, na, nb, betaParameter,
                                propDiff = propDiff
    )
  }
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]
  thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

  logLikelihoodRatioMultiBern(
    ya = ya, yb = yb, na = na, nb = nb,
    numeratorThetaA = thetaA, numeratorThetaB = thetaB,
    denominatorThetaA = thetaNull, denominatorThetaB = thetaNull,
    log = TRUE
  )
}
# Predictable plug-in for the numerator, for block i given the counts of
# blocks 1 to i - 1 only. predictiveThetas2x2(): the independent Beta
# posterior means of thetaA and thetaB. predictiveThetas2x2PropDiff(): the
# numerator lives on the curve thetaB - thetaA = propDiff, and thetaA is the
# posterior mean under a grid posterior on that curve.
predictiveThetas2x2 <- function(ya, yb, na, nb, betaParameter) {
  nBlocks <- length(ya)
  betaA1 <- betaParameter[["betaA1"]]
  betaA2 <- betaParameter[["betaA2"]]
  betaB1 <- betaParameter[["betaB1"]]
  betaB2 <- betaParameter[["betaB2"]]

  # Successes and sizes accumulated over blocks 1 to i - 1.
  previousYa <- c(0, cumsum(ya))[seq_len(nBlocks)]
  previousYb <- c(0, cumsum(yb))[seq_len(nBlocks)]
  previousNa <- c(0, cumsum(na))[seq_len(nBlocks)]
  previousNb <- c(0, cumsum(nb))[seq_len(nBlocks)]

  return(list(
    "thetaA" = (betaA1 + previousYa) / (betaA1 + betaA2 + previousNa),
    "thetaB" = (betaB1 + previousYb) / (betaB1 + betaB2 + previousNb)
  ))
}

predictiveThetas2x2PropDiff <- function(ya, yb, na, nb, betaParameter,
                                        propDiff, nWeight = 1e3L) {
  nBlocks <- length(ya)

  # output placeholder
  thetaA <- numeric(nBlocks)

  # Only the thetaA prior is used
  betaA1 <- betaParameter[["betaA1"]]
  betaA2 <- betaParameter[["betaA2"]]

  # thetaA are restricted by propDiff
  rho <- seq(1 / nWeight, 1 - 1 / nWeight, length.out = nWeight)
  thetaAGrid <- max(0, -propDiff) + rho * (1 - abs(propDiff))
  thetaBGrid <- thetaAGrid + propDiff

  logThetaA <- log(thetaAGrid)
  logOneMinusThetaA <- log1p(-thetaAGrid)
  logThetaB <- log(thetaBGrid)
  logOneMinusThetaB <- log1p(-thetaBGrid)

  # Un-normalised log posterior weights, shifted so their maximum is 0: the
  # largest weight is then exactly 1 and the sum can neither underflow nor
  # overflow, however many blocks have been seen.
  logWeights <- (betaA1 - 1) * log(rho) + (betaA2 - 1) * log1p(-rho)
  logWeights <- logWeights - max(logWeights)

  for (i in seq_len(nBlocks)) {
    # posterior mean for thetaA
    weights <- exp(logWeights)
    thetaA[i] <- sum(thetaAGrid * weights) / sum(weights)

    logWeights <- logWeights +
      ya[i] * logThetaA + (na[i] - ya[i]) * logOneMinusThetaA +
      yb[i] * logThetaB + (nb[i] - yb[i]) * logOneMinusThetaB
    logWeights <- logWeights - max(logWeights)
  }

  list("thetaA" = thetaA, "thetaB" = thetaA + propDiff)
}

# Find the means that minimize the KL between the alternative and null
# The null is H0: thetaB - thetaA = propDiff
# TODO: this can be a cubic function solver
solveRIPr2x2PropDiff <- function(thetaA, thetaB, na, nb, propDiff) {
  derivativeKL <- function(nullThetaA) {
    nullThetaB <- nullThetaA + propDiff
    na * ((1 - thetaA) / (1 - nullThetaA) - thetaA / nullThetaA) +
      nb * ((1 - thetaB) / (1 - nullThetaB) - thetaB / nullThetaB)
  }

  # The derivative is infinite at the edges, so search just inside them.
  stats::uniroot(derivativeKL,
    lower = max(0, -propDiff) + 1e-12, upper = min(1, 1 - propDiff) - 1e-12,
    tol = 1e-12
  )[["root"]]
}

## logOdds ----

# Cumulative log e-value against thetaA = thetaB (logOdds scale); length
# nBlocks. Predictable plug-in alternative for block i: logit(thetaB) -
# logit(thetaA) from the Beta posterior means given blocks 1..i-1.
computeEValueVecLogOdds <- function(ya, yb, na, nb, betaParameter) {
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  plugInLogOdds <- stats::qlogis(thetas[["thetaB"]]) -
    stats::qlogis(thetas[["thetaA"]])

  cumsum(mapply(
    logLikelihoodRatioFNCH,
    ya = ya, yb = yb, na = na, nb = nb, logOdds = plugInLogOdds,
    MoreArgs = list(log = TRUE)
  ))
}


# Conditional likelihood ratio at a fixed logOdds (B minus A) against the
# null nullLogOdds, given the block's total successes ya + yb: under either
# value yb is Fisher's noncentral hypergeometric with odds exp(logOdds) on group
# B (central at 0).
logLikelihoodRatioFNCH <- function(ya, yb, na, nb, logOdds,
                                   nullLogOdds = 0, log = FALSE) {
  totalSuccesses <- ya + yb
  logLR <- log(BiasedUrn::dFNCHypergeo(yb, nb, na, totalSuccesses, exp(logOdds))) -
    log(BiasedUrn::dFNCHypergeo(yb, nb, na, totalSuccesses, exp(nullLogOdds)))

  if (log) logLR else exp(logLR)
}

# Log partition function of Fisher's noncentral hypergeometric distribution
fnchLogPartition <- function(na, nb, totalSuccesses, logOdds) {
  if (logOdds == 0) {
    return(lchoose(na + nb, totalSuccesses))
  }

  feasibleSuccesses <-
    max(0, totalSuccesses - nb):min(na, totalSuccesses)
  logTerms <- lchoose(na, feasibleSuccesses) +
    lchoose(nb, totalSuccesses - feasibleSuccesses) +
    logOdds * feasibleSuccesses

  # log-sum-exp, shifted by the largest term against overflow
  maxLogTerm <- max(logTerms)
  maxLogTerm + log(sum(exp(logTerms - maxLogTerm)))
}

# UMP plug-in for one table: the logOdds on the side of the alternative at which
#   KL(FNCH(logOdds) || FNCH(nullLogOdds)) = log(1 / alpha),
# with KL = (logOdds - nullLogOdds) * E_logOdds[yb]
#           - fnchLogPartition(logOdds) + fnchLogPartition(nullLogOdds),
# the mean E_logOdds[yb] taken from BiasedUrn (odds = exp(logOdds) on B). At
# nullLogOdds = 0 the last term is lchoose(na + nb, ya + yb).
# The KL is 0 at nullLogOdds and increases away from it, but is bounded by
# -log P0(ya at its feasible extreme), so the equation may have no root:
# then NULL is returned and the caller uses the trivial e-factor 1.
solveUmpLogOdds <- function(na, nb, totalSuccesses, alpha,
                          alternative = c("greater", "less"),
                          nullLogOdds = 0,
                          searchBound = 100) {
  alternative <- match.arg(alternative)

  # logOdds is B minus A, so the weighted count is yb: group B goes first.
  klMinusTarget <- function(logOdds) {
    (logOdds - nullLogOdds) *
      BiasedUrn::meanFNCHypergeo(nb, na, totalSuccesses, exp(logOdds)) -
      fnchLogPartition(nb, na, totalSuccesses, logOdds) +
      fnchLogPartition(nb, na, totalSuccesses, nullLogOdds) + log(alpha)
  }

  # The KL is bounded, so the target may be unreachable within the search
  # interval; uniroot() would error on equal signs, hence the explicit check.
  bounds <- if (alternative == "greater") {
    c(nullLogOdds, nullLogOdds + searchBound)
  } else {
    c(nullLogOdds - searchBound, nullLogOdds)
  }
  if (klMinusTarget(bounds[1]) * klMinusTarget(bounds[2]) > 0) {
    return(NULL)
  }

  stats::uniroot(klMinusTarget,
    lower = bounds[1], upper = bounds[2],
    tol = 1e-10
  )[["root"]]
}
