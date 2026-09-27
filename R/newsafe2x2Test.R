# Testing fnts ----

#' Safe anytime-valid 2x2 test for propDiff
#'
#' `eType = "eBeta"`: the unrestricted numerator of Decision 3 (`twoSided`
#' only). `eType = "grow"`: the numerator restricted to
#' `thetaB - thetaA = propDiffMin` (Decision 6), `"greater"` only. Only
#' eBeta gets a confidence interval (Decision 23), on all blocks;
#' `wantConfidenceSequence = TRUE` adds the blockwise `confSeqMatrix`
#' (Decisions 24, 26).
#' @noRd
savi2x2TestStatPropDiff <- function(ya, yb,
                                    designObj = NULL, wantCi = TRUE,
                                    wantConfidenceSequence = FALSE, ciValue = NULL) {
  # grow must have a propDiffMin; eBeta must have none and be twoSided
  propDiffMin <- designObj[["esMin"]]
  alternative <- designObj[["alternative"]]
  betaParameter <- designObj[["betaParameter"]]
  runningIntersection <- designObj[["runningIntersection"]]
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

  # confidence interval or sequences only in eBeta
  # use 1 - alpha unless user specified
  ciValue <- ifelse(is.null(ciValue), 1 - designObj[["alpha"]], ciValue)
  result[["ciValue"]] <- ciValue

  if (wantConfidenceSequence && eType == "eBeta") {
    confSeqMatrix <- computeConfidenceSequence2x2PropDiff(
      ya, yb, na, nb, betaParameter, 1 - ciValue, runningIntersection
    )
    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSeqMatrix[nBlocks, ]
  } else if (wantCi && eType == "eBeta") {
    # One confidence interval on last block
    result[["confSeq"]] <- computeConfidenceInterval2x2PropDiff(
      ya, yb, na, nb, betaParameter, 1 - ciValue
    )
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
#' (Decision 17), inverting its e-process on point nulls, with the first
#' block's factor the UMP plug-in (Decision 28); grow gets none.
#' @noRd
savi2x2TestStatLogOdds <- function(ya, yb,
                                   designObj = NULL, wantCi = TRUE, ciValue = NULL) {
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
    # The first block's factor is the UMP plug-in instead of the prior
    # mixture; the posterior still absorbs block 1 (Decision 28).
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
      logMixCum <- shift + log(rowSums(exp(logMix - shift)))
      # Block 1: UMP plug-in solved from its total at the design's alpha;
      # NULL (KL target out of reach) is the trivial e-factor 1.
      logOddsUmp <- solveUmpLogOdds(na[1], nb[1], ya[1] + yb[1], alpha,
                                    "greater")
      logPUmp <- if (is.null(logOddsUmp)) {
        logP0[1]
      } else {
        logLikelihoodFNCHVec(ya[1], yb[1], na[1], nb[1], logOddsUmp)
      }
      logPUmp + logMixCum - logMixCum[1]
    },
    stop("eType ", eType, " is not implemented for logOdds")
  )
  logEvalueVec <- logPCum - cumsum(logP0)

  # Compute: confSeq ----
  # grow: no confidence interval until its construction is agreed.
  # eGauss inverts its own e-process on point nulls (Decisions 17, 28);
  # the level is ciValue, 1 - alpha unless the user specified one.
  ciValue <- ifelse(is.null(ciValue), 1 - designObj[["alpha"]], ciValue)
  result[["ciValue"]] <- ciValue
  if (wantCi && eType == "eGauss") {
    result[["confSeq"]] <- computeConfidenceInterval2x2LogOdds(
      ya, yb, na, nb, logPCum[nBlocks], 1 - ciValue
    )
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
#' numerator is the predictable Beta posterior means `predictiveThetas2x2()`
#' (block `i` from blocks `1..i-1`), and each block's null `thetaA` is its
#' projection `solveRIPr2x2PropDiff()` onto that line. The log e-value is
#' convex in `propDiff`, so the kept set `{propDiff : f < log(1 / alpha)}`
#' is one interval or empty, and searching only inside `domain` returns its
#' intersection with `domain`.
#'
#' @param domain `c(lower, upper)` inside `[-1, 1]`, the candidates searched;
#'   a previous interval gives the running intersection.
#' @return Named numeric `c(lowerBound, upperBound)`; the `domain` edge when
#'   that edge is still inside, both `NA` when the set is empty.
#' @noRd
computeConfidenceInterval2x2PropDiff <- function(ya, yb, na, nb,
                                                 betaParameter, alpha,
                                                 domain = c(-1, 1)) {
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]

  # product of the seq-RIPr e-variable against the H0: thetaB = thetaA + propDiff
  # f: find zero points against 1/alpha
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

  # The projection needs thetaA strictly inside its range, so stay just
  # inside (-1, 1).
  eps <- 1e-9
  lower <- max(domain[1], -1 + eps)
  upper <- min(domain[2], 1 - eps)
  minimiser <- stats::optimize(fPropDiff,
    interval = c(lower, upper)
  )[["minimum"]]

  # min > 1 / alpha, no confidence interval found
  if (fPropDiff(minimiser) >= 0) {
    warning("No confidence interval is found!")
    return(c("lowerBound" = -1, "upperBound" = 1))
  }

  # Still inside at the edge: the bound is the domain edge.
  lowerBound <- if (fPropDiff(lower) < 0) {
    domain[1]
  } else {
    stats::uniroot(fPropDiff, lower = lower, upper = minimiser)[["root"]]
  }
  upperBound <- if (fPropDiff(upper) < 0) {
    domain[2]
  } else {
    stats::uniroot(fPropDiff, lower = minimiser, upper = upper)[["root"]]
  }

  return(c("lowerBound" = lowerBound, "upperBound" = upperBound))
}

#' Anytime-valid confidence sequence for propDiff, one row per block
#'
#' Row `i` inverts the eBeta test on blocks `1..i` against point nulls
#' `thetaB - thetaA = propDiff`, on a fixed grid of candidates whose
#' cumulative log e-process is advanced once per block (Decision 26). The
#' grid has 2000 candidates; a warning says when its step is coarser than
#' `sdMax`, the worst-case Wald standard deviation of the difference on the
#' observed totals, the scale of the narrowest (last) interval. The null
#' `thetaA` of every
#' candidate is the root of a cubic, found for all candidates at once by
#' bisection. The kept candidates form one run (convexity, Decision 13) and
#' each bound is refined outward by a secant, so the reported interval
#' contains the exact one.
#'
#' @return `nBlocks x 2` matrix, `lowerBound` and `upperBound`; `-1` / `1`
#'   when the outermost candidate is kept, both `NA` when the set is empty.
#'   With `runningIntersection` the rows are nested and stay `NA` after the
#'   first empty row.
#' @noRd
computeConfidenceSequence2x2PropDiff <- function(ya, yb, na, nb,
                                                 betaParameter, alpha,
                                                 runningIntersection) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]
  logThreshold <- log(1 / alpha)

  # 2000 candidates strictly inside (-1, 1). The final interval has width
  # on the scale of sdMax; warn when the grid step is coarser than that.
  nGrid <- 2000L
  sdMax <- sqrt(1 / (4 * sum(na)) + 1 / (4 * sum(nb)))
  if (ceiling(2 / sdMax) > nGrid) {
    warning("The confidence sequence grid step ", 2 / nGrid,
      " is coarser than the standard deviation scale ", signif(sdMax, 3),
      " of propDiff on these totals; the bounds are conservative")
  }
  grid <- seq(-1, 1, length.out = nGrid + 2L)[-c(1L, nGrid + 2L)]

  # Feasible null thetaA per candidate: both thetaA and thetaA + delta in
  # [0, 1].
  gridLower <- pmax(0, -grid)
  gridUpper <- pmin(1, 1 - grid)

  # Cumulative log e-process per candidate; which candidates may still be
  # kept, and which take part in the update (the same, plus two guard
  # nodes on each side for the secant under the running intersection).
  logEValues <- numeric(nGrid)
  candidate <- rep(TRUE, nGrid)
  active <- rep(TRUE, nGrid)
  previous <- c(-1, 1)

  confSeqMatrix <- matrix(NA_real_, nBlocks, 2,
    dimnames = list(NULL, c("lowerBound", "upperBound"))
  )

  for (i in seq_len(nBlocks)) {
    delta <- grid[active]

    # Null thetaA for every active candidate: the KL derivative of
    # solveRIPr2x2PropDiff() times x (1 - x)(x + delta)(1 - x - delta) is a
    # cubic in x, negative then positive across its single root on the
    # feasible interval. 52 bisection steps reach machine precision.
    lower <- gridLower[active]
    upper <- gridUpper[active]
    for (step in seq_len(52L)) {
      x <- (lower + upper) / 2
      cubic <- na[i] * (x - thetaA[i]) * (x + delta) * (1 - x - delta) +
        nb[i] * (x + delta - thetaB[i]) * x * (1 - x)
      positive <- cubic > 0
      upper[positive] <- x[positive]
      lower[!positive] <- x[!positive]
    }
    nullThetaA <- (lower + upper) / 2

    # Block i's log likelihood ratio term against each active candidate.
    logEValues[active] <- logEValues[active] +
      stats::dbinom(ya[i], na[i], thetaA[i], log = TRUE) +
      stats::dbinom(yb[i], nb[i], thetaB[i], log = TRUE) -
      stats::dbinom(ya[i], na[i], nullThetaA, log = TRUE) -
      stats::dbinom(yb[i], nb[i], nullThetaA + delta, log = TRUE)

    # Kept candidates; one run by convexity.
    f <- logEValues - logThreshold
    kept <- candidate & f < 0
    if (!any(kept)) {
      if (runningIntersection) break
      next
    }
    keptRange <- range(which(kept))

    # Refine each bound outward: the secant through the first rejected node
    # and its outer neighbour lies below the convex f outside them, so its
    # root is outside the true boundary. At the grid edge the bound is the
    # rejected node itself, or -1 / 1 when even the outermost node is kept.
    lowerBound <- if (keptRange[1] == 1L) {
      -1
    } else {
      j <- keptRange[1] - 1L
      if (j == 1L) {
        grid[j]
      } else {
        grid[j] - f[j] * (grid[j] - grid[j - 1L]) / (f[j] - f[j - 1L])
      }
    }
    upperBound <- if (keptRange[2] == nGrid) {
      1
    } else {
      j <- keptRange[2] + 1L
      if (j == nGrid) {
        grid[j]
      } else {
        grid[j] - f[j] * (grid[j + 1L] - grid[j]) / (f[j + 1L] - f[j])
      }
    }

    # Running intersection: a candidate that leaves never returns, so drop
    # it from the update, keeping the two nodes beyond the run on each side
    # up to date for the secant; intersect with the previous row, and once
    # empty the remaining rows stay NA.
    if (runningIntersection) {
      candidate <- kept
      active <- seq_len(nGrid) >= keptRange[1] - 2L &
        seq_len(nGrid) <= keptRange[2] + 2L
      lowerBound <- max(lowerBound, previous[1])
      upperBound <- min(upperBound, previous[2])
      if (lowerBound > upperBound) break
      previous <- c(lowerBound, upperBound)
    }
    confSeqMatrix[i, ] <- c(lowerBound, upperBound)
  }

  return(confSeqMatrix)
}

#' Anytime-valid confidence interval for logOdds on all blocks
#'
#' Inverts the conditional test on point nulls `logOdds = delta`. The
#' alternative is fixed before the data, so the log e-value on all blocks
#' against `delta` is `logPTotal`, the numerator's total log likelihood,
#' minus the FNCH log likelihood of the data at `delta` (Decision 17). That
#' log likelihood is concave in `delta`, so `f` is convex with its minimum
#' at the conditional MLE and the kept set `{delta : f < log(1 / alpha)}`
#' is one interval or empty. The interval is only claimed on
#' `(-logOddsBound, logOddsBound)`.
#'
#' @param logPTotal The numerator's cumulative log likelihood after the
#'   last block, e.g. `logPCum[nBlocks]` of `savi2x2TestStatLogOdds()`.
#' @return Named numeric `c(lowerBound, upperBound)`; `-logOddsBound` /
#'   `logOddsBound` when that edge is still inside, and the whole range
#'   with a warning when the set is empty.
#' @noRd
computeConfidenceInterval2x2LogOdds <- function(ya, yb, na, nb, logPTotal,
                                                alpha, logOddsBound = 40) {
  # product of the conditional e-variable against H0: logOdds = delta
  # f: find zero points against 1 / alpha
  fLogOdds <- function(logOdds) {
    logPTotal -
      sum(logLikelihoodFNCHVec(ya, yb, na, nb, logOdds)) -
      log(1 / alpha)
  }

  minimiser <- stats::optimize(fLogOdds,
    interval = c(-logOddsBound, logOddsBound)
  )[["minimum"]]

  # min > 1 / alpha, no confidence interval found
  if (fLogOdds(minimiser) >= 0) {
    warning("No confidence interval is found!")
    return(c("lowerBound" = -logOddsBound, "upperBound" = logOddsBound))
  }

  # Still inside at the search edge: the bound is the edge itself.
  lowerBound <- if (fLogOdds(-logOddsBound) < 0) {
    -logOddsBound
  } else {
    stats::uniroot(fLogOdds, lower = -logOddsBound, upper = minimiser)[["root"]]
  }
  upperBound <- if (fLogOdds(logOddsBound) < 0) {
    logOddsBound
  } else {
    stats::uniroot(fLogOdds, lower = minimiser, upper = logOddsBound)[["root"]]
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
