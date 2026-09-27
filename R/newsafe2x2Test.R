# Testing fnts ----

#' Safe anytime-valid 2x2 test for propDiff
#'
#' `eType = "eBeta"`: the unrestricted numerator of Decision 3 (`twoSided`
#' only). `eType = "grow"`: the numerator restricted to
#' `thetaB - thetaA = propDiffMin` (Decision 6), `"greater"` only. Only
#' eBeta gets a confidence interval (Decision 23), on all blocks.
#' @noRd
savi2x2TestStatPropDiff <- function(ya, yb,
                                    designObj = NULL, wantCi = TRUE) {
  # TODO: THESE ARGS CHECKING WILL BE DONE FINAL STEP
  # grow must have a propDiffMin; eBeta must have none and be twoSided
  propDiffMin <- designObj[["esMin"]]
  alternative <- designObj[["alternative"]]
  betaParameter <- designObj[["betaParameter"]]
  alpha <- designObj[["alpha"]]
  eType <- designObj[["eType"]]
  na <- designObj[["nPlan"]][["na"]]
  nb <- designObj[["nPlan"]][["nb"]]

  nBlocks <- length(ya)
  if (length(na) == 1L) na <- rep(na, nBlocks)
  if (length(nb) == 1L) nb <- rep(nb, nBlocks)

  # Data checks: one count and size per block, counts within their sizes.
  if (length(yb) != nBlocks || length(na) != nBlocks ||
      length(nb) != nBlocks) {
    stop("ya, yb, na and nb must have one value per block")
  }
  counts <- c(ya, yb, na, nb)
  if (!all(is.finite(counts)) || any(counts %% 1 != 0) ||
      any(c(ya, yb) < 0) || any(c(na, nb) < 1) ||
      any(ya > na) || any(yb > nb)) {
    stop("ya, yb must be integers in 0..na, 0..nb; na, nb positive integers")
  }

  # Compute: eValueVec ----
  # Numerator: predictable thetas from blocks 1..i-1 only.
  # eBeta (twoSided): independent Beta posterior means.
  # grow (greater only for now): learnt on the curve thetaB - thetaA =
  # propDiffMin. grow + twoSided: TODO later.
  thetas <- switch(eType,
    "eBeta" = predictiveThetas2x2(ya, yb, na, nb, betaParameter),
    "grow" = predictiveThetas2x2PropDiff(ya, yb, na, nb, betaParameter,
      propDiff = propDiffMin
    ),
    stop("eType ", eType, " is not implemented for propDiff")
  )
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]

  # Null: projection onto thetaA = thetaB, the size-weighted pooled mean.
  thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

  # Cumulative log likelihood ratio; dbinom() handles theta at 0 or 1.
  logEValueVec <- cumsum(
    stats::dbinom(ya, na, thetaA, log = TRUE) +
      stats::dbinom(yb, nb, thetaB, log = TRUE) -
      stats::dbinom(ya, na, thetaNull, log = TRUE) -
      stats::dbinom(yb, nb, thetaNull, log = TRUE)
  )

  # Compute: confSeq ----
  result <- constructSaviTestObj("Two Proportions")

  # eBeta only, on all blocks (Decision 23).
  if (wantCi && eType == "eBeta") {
    result[["confSeq"]] <- computeConfidenceInterval2x2PropDiff(
      ya, yb, na, nb, thetaA, thetaB, alpha
    )
    result[["ciValue"]] <- 1 - alpha
  }

  # Fill: Result ----
  eValueVec <- exp(logEValueVec)
  result[["estimate"]] <- c(
    "thetaA" = sum(ya) / sum(na), "thetaB" = sum(yb) / sum(nb)
  )
  # x-axis of plot.saviTest(): the block index.
  result[["n1Vec"]] <- seq_len(nBlocks)
  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["n"]] <- c("na" = sum(na), "nb" = sum(nb), "nBlocks" = nBlocks)
  # Beta posterior of the observed blocks, the prior a further block would use.
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
#'
#' Every block is conditioned on its total successes (Decision 7).
#' `eType = "grow"`: the fixed alternative `logOddsMin` (Decision 15),
#' `"greater"` only. `eType = "eGauss"`: N(0, 1) mixture on a fixed grid
#' (Decision 18), `"twoSided"` only. Only eGauss gets a confidence interval
#' (Decision 17), inverting its e-process on point nulls; grow gets none.
#' @noRd
savi2x2TestStatLogOdds <- function(ya, yb,
                                   designObj = NULL, wantCi = TRUE) {
  # TODO: THESE ARGS CHECKING WILL BE DONE FINAL STEP
  # grow must have a logOddsMin > 0
  nBlocks <- length(ya)
  logOddsMin <- designObj[["esMin"]]
  alternative <- designObj[["alternative"]]
  alpha <- designObj[["alpha"]]
  eType <- designObj[["eType"]]
  na <- designObj[["nPlan"]][["na"]]
  nb <- designObj[["nPlan"]][["nb"]]
  if (length(na) == 1L) na <- rep(na, nBlocks)
  if (length(nb) == 1L) nb <- rep(nb, nBlocks)

  # Data checks: one count and size per block, counts within their sizes.
  if (length(yb) != nBlocks || length(na) != nBlocks ||
      length(nb) != nBlocks) {
    stop("ya, yb, na and nb must have one value per block")
  }
  counts <- c(ya, yb, na, nb)
  if (!all(is.finite(counts)) || any(counts %% 1 != 0) ||
      any(c(ya, yb) < 0) || any(c(na, nb) < 1) ||
      any(ya > na) || any(yb > nb)) {
    stop("ya, yb must be integers in 0..na, 0..nb; na, nb positive integers")
  }

  result <- constructSaviTestObj("Two Proportions")

  # Numerator: cumulative log likelihood of blocks 1..i under the alternative.
  logP0 <- stats::dhyper(yb, nb, na, ya + yb, log = TRUE)
  logPCum <- switch(eType,
    # grow: the fixed alternative logOddsMin (Decision 15).
    grow = cumsum(logLikelihoodFNCHVec(ya, yb, na, nb, logOddsMin)),
    # eGauss: N(0, 1) prior on a fixed logOdds grid, twoSided (Decision 18).
    eGauss = {
      logOddsGrid <- seq(-20, 20, length.out = 2000)
      logPrior <- stats::dnorm(logOddsGrid, log = TRUE)
      logPrior <- logPrior - max(logPrior) -
        log(sum(exp(logPrior - max(logPrior))))
      # nBlocks x grid: FNCH log density of yb at every grid logOdds.
      logPGrid <- t(mapply(function(ya, yb, na, nb) {
        k <- max(0, ya + yb - na):min(nb, ya + yb)
        logTerms <- outer(logOddsGrid, k) +
          rep(lchoose(nb, k) + lchoose(na, ya + yb - k),
              each = length(logOddsGrid))
        shift <- apply(logTerms, 1, max)
        logTerms[, yb - k[1] + 1] - shift - log(rowSums(exp(logTerms - shift)))
      }, ya = ya, yb = yb, na = na, nb = nb))
      # Cumulate over blocks, then mix over the grid (log-sum-exp per block).
      logMix <- sweep(
        matrix(apply(logPGrid, 2, cumsum), nrow = nBlocks), 2, logPrior, "+"
      )
      shift <- apply(logMix, 1, max)
      shift + log(rowSums(exp(logMix - shift)))
    },
    stop("eType ", eType, " is not implemented for logOdds")
  )
  logEvalueVec <- logPCum - cumsum(logP0)

  # confSeq ----
  # The alternative is fixed before the data, so the log e-value on all
  # blocks against the point null logOdds = delta is the numerator's total
  # log likelihood minus the FNCH log likelihood at delta (Decision 17). It
  # is convex in delta, so the kept set {delta : f < log(1 / alpha)} is one
  # interval or empty.
  fLogOdds <- function(logPTotal, logOdds, alpha) {
    logPTotal -
      sum(logLikelihoodFNCHVec(ya, yb, na, nb, logOdds)) -
      log(1 / alpha)
  }

  # grow: no confidence interval until its construction is agreed.
  if (wantCi && eType == "eGauss") {
    logOddsBound <- 40
    # optimize() and uniroot() vary the one argument left unnamed, logOdds.
    minimiser <- stats::optimize(fLogOdds,
      interval = c(-logOddsBound, logOddsBound),
      logPTotal = logPCum[nBlocks], alpha = alpha
    )[["minimum"]]

    if (fLogOdds(logPCum[nBlocks], minimiser, alpha) >= 0) {
      # Even the conditional MLE is rejected: the set is empty.
      confSeq <- c("lowerBound" = NA_real_, "upperBound" = NA_real_)
    } else {
      # Still inside at the search edge: report the edge itself, the
      # interval is only claimed on (-logOddsBound, logOddsBound).
      lowerBound <- if (fLogOdds(logPCum[nBlocks], -logOddsBound, alpha) < 0) {
        -logOddsBound
      } else {
        stats::uniroot(fLogOdds,
          lower = -logOddsBound, upper = minimiser,
          logPTotal = logPCum[nBlocks], alpha = alpha
        )[["root"]]
      }
      upperBound <- if (fLogOdds(logPCum[nBlocks], logOddsBound, alpha) < 0) {
        logOddsBound
      } else {
        stats::uniroot(fLogOdds,
          lower = minimiser, upper = logOddsBound,
          logPTotal = logPCum[nBlocks], alpha = alpha
        )[["root"]]
      }
      confSeq <- c("lowerBound" = lowerBound, "upperBound" = upperBound)
    }

    result[["confSeq"]] <- confSeq
    result[["ciValue"]] <- 1 - alpha
  }

  # Fill: Result ----
  eValueVec <- exp(logEvalueVec)
  result[["estimate"]] <- c(
    "thetaA" = sum(ya) / sum(na), "thetaB" = sum(yb) / sum(nb)
  )
  # x-axis of plot.saviTest(): the block index.
  result[["n1Vec"]] <- seq_len(nBlocks)
  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["n"]] <- c("na" = sum(na), "nb" = sum(nb), "nBlocks" = nBlocks)
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

# Design fnts ----
# TODO: add "less" when direction is clean in `alternative`

#' Design a safe anytime-valid 2x2 test
#'
#' `eType` picks the effect: `"eBeta"` (propDiff) is unrestricted and
#' `"twoSided"`; `"grow"` plugs in whichever of `propDiffMin`, `logOddsMin`
#' is set and is `"greater"` only. Inputs are assumed valid.
#' @noRd
designSavi2x2 <- function(
  na, nb, propDiffMin = NULL, logOddsMin = NULL,
  alpha = 0.05, power = NULL, h0 = 0,
  alternative = c("twoSided", "greater"),
  eType = c("eBeta", "grow", "eGauss"),
  betaParameter = NULL,
  runningIntersection = NULL
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  result <- constructSaviDesignObj("Two Proportions")

  # NULL keeps the constructor's defaults.
  if (!is.null(betaParameter)) {
    result[["betaParameter"]] <- betaParameter
  }
  if (!is.null(runningIntersection)) {
    result[["runningIntersection"]] <- runningIntersection
  }

  # TODO: what if both are NULL
  result[["esMin"]] <- if (!is.null(propDiffMin)) propDiffMin else logOddsMin
  result[["parameter"]] <- c(
    "Beta hyperparameters" =
      paste(unlist(result[["betaParameter"]]), collapse = " ")
  )
  result[["eType"]] <- eType
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

#' Anytime-valid confidence interval for propDiff on all blocks
#'
#' Inverts the eBeta test on point nulls `thetaB - thetaA = propDiff`: the
#' numerator is the test's own predictable `thetaA`, `thetaB` (one per
#' block, from blocks `1..i-1`), and each block's null `thetaA` is its
#' projection `solveRIPr2x2PropDiff()` onto that line. The log e-value is
#' convex in `propDiff`, so the kept set `{propDiff : f < log(1 / alpha)}`
#' is one interval or empty.
#'
#' @return Named numeric `c(lowerBound, upperBound)`; `-1` or `1` when that
#'   edge is still inside, both `NA` when the set is empty.
#' @noRd
computeConfidenceInterval2x2PropDiff <- function(ya, yb, na, nb,
                                                 thetaA, thetaB, alpha) {
  # f(propDiff): log e-value on all blocks against propDiff, minus the
  # threshold log(1 / alpha).
  fPropDiff <- function(propDiff) {
    nullThetaA <- mapply(solveRIPr2x2PropDiff,
      thetaA = thetaA, thetaB = thetaB, na = na, nb = nb,
      MoreArgs = list(propDiff = propDiff)
    )
    sum(
      stats::dbinom(ya, na, thetaA, log = TRUE) +
        stats::dbinom(yb, nb, thetaB, log = TRUE) -
        stats::dbinom(ya, na, nullThetaA, log = TRUE) -
        stats::dbinom(yb, nb, nullThetaA + propDiff, log = TRUE)
    ) - log(1 / alpha)
  }

  # The projection needs thetaA strictly inside its range, so search just
  # inside (-1, 1).
  eps <- 1e-9
  minimiser <- stats::optimize(fPropDiff,
    interval = c(-1 + eps, 1 - eps)
  )[["minimum"]]

  # Even the minimiser is rejected: the set is empty.
  if (fPropDiff(minimiser) >= 0) {
    return(c("lowerBound" = NA_real_, "upperBound" = NA_real_))
  }

  # Still inside at the edge: the bound is the parameter limit.
  lowerBound <- if (fPropDiff(-1 + eps) < 0) {
    -1
  } else {
    stats::uniroot(fPropDiff, lower = -1 + eps, upper = minimiser)[["root"]]
  }
  upperBound <- if (fPropDiff(1 - eps) < 0) {
    1
  } else {
    stats::uniroot(fPropDiff, lower = minimiser, upper = 1 - eps)[["root"]]
  }

  return(c("lowerBound" = lowerBound, "upperBound" = upperBound))
}

# Helpers ----

## propDiff ----
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

# Per-block conditional log likelihood at one logOdds (B minus A): given the
# block's total ya + yb, yb is Fisher's noncentral hypergeometric with odds
# exp(logOdds) on group B. A vector of length nBlocks; dFNCHypergeo() takes
# scalar sizes, hence the loop over blocks.
logLikelihoodFNCHVec <- function(ya, yb, na, nb, logOdds) {
  mapply(function(ya, yb, na, nb) {
    log(BiasedUrn::dFNCHypergeo(yb, nb, na, ya + yb, exp(logOdds)))
  }, ya = ya, yb = yb, na = na, nb = nb)
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
