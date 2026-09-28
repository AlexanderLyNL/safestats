# Test functions ----

# Conditional e-factor of one block at a UMP plug-in logOdds: the log
# likelihood of yb given the block's total under the alternative (FNCH with
# odds exp(logOdds) on group B) minus the same under the null
# (hypergeometric). twoSided averages the two one-sided e-factors.
savi2x2TestStatUmp <- function(ya, yb, na, nb, alpha,
                               alternative = c("twoSided", "greater", "less")) {
  alternative <- match.arg(alternative)
  logLikelihoodNull <- stats::dhyper(yb, nb, na, ya + yb, log = TRUE)

  # UMP plug-in on each side (Decision 8): the logOdds at which the
  # conditional KL against the null reaches log(1 / alpha), solved from the
  # block's total only. NULL means the target is out of reach; the plug-in
  # is then the null itself, logOdds = 0, i.e. the trivial e-factor 1.
  logOddsPositive <- solveUmpLogOdds(na, nb, ya + yb, alpha, "greater")
  if (is.null(logOddsPositive)) logOddsPositive <- 0
  logOddsNegative <- solveUmpLogOdds(na, nb, ya + yb, alpha, "less")
  if (is.null(logOddsNegative)) logOddsNegative <- 0

  logLikelihoodPositive <- logLikelihoodFNCH(ya, yb, na, nb, logOddsPositive)
  logLikelihoodNegative <- logLikelihoodFNCH(ya, yb, na, nb, logOddsNegative)

  switch(alternative,
    "twoSided" = 0.5 * exp(logLikelihoodPositive - logLikelihoodNull) +
      0.5 * exp(logLikelihoodNegative - logLikelihoodNull),
    "greater" = exp(logLikelihoodPositive - logLikelihoodNull),
    "less" = exp(logLikelihoodNegative - logLikelihoodNull)
  )
}

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

  if (nBlocks == 1L) {
    warnings("There is only 1 table, switched to conditional e-variable")
  }

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

  # Compute: eValueVec ----
  if (eType == "grow") {
    # greater only for now; twoSided later. The blockwise e-process carries
    # the UMP factor of block 1 itself.
    logEValueVec <- logEProcess2x2PropDiffGrow(
      ya, yb, na, nb, betaParameter, propDiffMin, alpha
    )
  } else if (eType == "eBeta") {
    # Block 1's UMP conditional e-factor multiplies the whole e-process.
    eValueUmp <- savi2x2TestStatUmp(
      ya[1], yb[1], na[1], nb[1], alpha, alternative
    )

    # Numerator (twoSided): independent Beta posterior means of thetaA and
    # thetaB from blocks 1..i-1 only.
    thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
    thetaA <- thetas[["thetaA"]]
    thetaB <- thetas[["thetaB"]]

    # Null: projection onto thetaA = thetaB, the size-weighted pooled mean.
    thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

    # Cumulative log likelihood of blocks 1..i under the null (denominator)
    # and under the alternative (numerator); dbinom() handles theta at 0 or 1.
    logLikelihoodNull <- cumsum(
      stats::dbinom(ya, na, thetaNull, log = TRUE) +
        stats::dbinom(yb, nb, thetaNull, log = TRUE)
    )
    logLikelihoodAlternative <- cumsum(
      stats::dbinom(ya, na, thetaA, log = TRUE) +
        stats::dbinom(yb, nb, thetaB, log = TRUE)
    )
    logEValueVec <- logLikelihoodAlternative - logLikelihoodNull + log(eValueUmp)
  } else {
    stop("eType ", eType, " is not implemented for propDiff")
  }

  # Compute: confSeq ----
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
#' `wantConfidenceSequence = TRUE` adds the blockwise `confSeqMatrix`
#' (Decision 30), row `i` the same inversion on blocks `1..i`.
#' @noRd
savi2x2TestStatLogOdds <- function(ya, yb,
                                   designObj = NULL, wantCi = TRUE,
                                   wantConfidenceSequence = FALSE,
                                   ciValue = NULL) {
  # grow must have a logOddsMin > 0
  nBlocks <- length(ya)
  logOddsMin <- designObj[["esMin"]]
  alternative <- designObj[["alternative"]]
  alpha <- designObj[["alpha"]]
  eType <- designObj[["eType"]]
  runningIntersection <- designObj[["runningIntersection"]]
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

  # Compute: eValueVec ----
  # Cumulative log likelihood of blocks 1..i under the null (denominator,
  # hypergeometric given each block's total) and under the alternative
  # (numerator).
  logLikelihoodNull <- cumsum(stats::dhyper(yb, nb, na, ya + yb, log = TRUE))
  logLikelihoodAlternative <- switch(eType,
    # grow: the fixed alternative logOddsMin (Decision 15).
    grow = cumsum(logLikelihoodFNCH(ya, yb, na, nb, logOddsMin)),
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
        logLikelihoodNull[1]
      } else {
        logLikelihoodFNCH(ya[1], yb[1], na[1], nb[1], logOddsUmp)
      }
      logPUmp + logMixCum - logMixCum[1]
    },
    stop("eType ", eType, " is not implemented for logOdds")
  )
  logEValueVec <- logLikelihoodAlternative - logLikelihoodNull

  # Compute: confSeq ----
  # grow: no confidence interval until its construction is agreed.
  # eGauss inverts its own e-process on point nulls (Decisions 17, 28);
  # the level is ciValue, 1 - alpha unless the user specified one.
  ciValue <- ifelse(is.null(ciValue), 1 - designObj[["alpha"]], ciValue)
  result[["ciValue"]] <- ciValue
  if (wantConfidenceSequence && eType == "eGauss") {
    # Row i inverts the e-process on blocks 1..i: every numerator factor is
    # fixed given its block's total or predictable, so
    # logLikelihoodAlternative[i] is the numerator at block i just as its
    # last element is at the end (Decision 30). Recomputed from scratch per
    # block, quadratic in nBlocks.
    confSeqMatrix <- matrix(NA_real_, nBlocks, 2,
      dimnames = list(NULL, c("lowerBound", "upperBound"))
    )
    domain <- c(-40, 40)
    for (i in seq_len(nBlocks)) {
      # The empty set is reported with a warning; here it is an NA row.
      row <- tryCatch(
        computeConfidenceInterval2x2LogOdds(
          ya[1:i], yb[1:i], na[1:i], nb[1:i], logLikelihoodAlternative[i],
          1 - ciValue,
          domain = domain
        ),
        warning = function(w) c("lowerBound" = NA_real_, "upperBound" = NA_real_)
      )
      if (runningIntersection) {
        # A value that leaves never returns: search inside the previous
        # row, and once empty the remaining rows stay NA.
        if (anyNA(row)) break
        domain <- row
      }
      confSeqMatrix[i, ] <- row
    }
    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSeqMatrix[nBlocks, ]
  } else if (wantCi && eType == "eGauss") {
    result[["confSeq"]] <- computeConfidenceInterval2x2LogOdds(
      ya, yb, na, nb, logLikelihoodAlternative[nBlocks], 1 - ciValue
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

# Design functions ----
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
#' is one interval or empty. The interval is only claimed on `domain`,
#' and searching only inside `domain` returns its intersection with
#' `domain` (convexity).
#'
#' @param logPTotal The numerator's cumulative log likelihood after the
#'   last block, `logLikelihoodAlternative[nBlocks]` of
#'   `savi2x2TestStatLogOdds()`.
#' @param domain `c(lower, upper)`, the candidates searched; a previous
#'   interval gives the running intersection.
#' @return Named numeric `c(lowerBound, upperBound)`; the `domain` edge
#'   when that edge is still inside, and the whole `domain` with a warning
#'   when the set is empty.
#' @noRd
computeConfidenceInterval2x2LogOdds <- function(ya, yb, na, nb, logPTotal,
                                                alpha, domain = c(-40, 40)) {
  # product of the conditional e-variable against H0: logOdds = delta
  # f: find zero points against 1 / alpha
  fLogOdds <- function(logOdds) {
    logPTotal -
      sum(logLikelihoodFNCH(ya, yb, na, nb, logOdds)) -
      log(1 / alpha)
  }

  minimiser <- stats::optimize(fLogOdds, interval = domain)[["minimum"]]

  # min > 1 / alpha, no confidence interval found
  if (fLogOdds(minimiser) >= 0) {
    warning("No confidence interval is found!")
    return(c("lowerBound" = domain[1], "upperBound" = domain[2]))
  }

  # Still inside at the search edge: the bound is the domain edge.
  lowerBound <- if (fLogOdds(domain[1]) < 0) {
    domain[1]
  } else {
    stats::uniroot(fLogOdds, lower = domain[1], upper = minimiser)[["root"]]
  }
  upperBound <- if (fLogOdds(domain[2]) < 0) {
    domain[2]
  } else {
    stats::uniroot(fLogOdds, lower = minimiser, upper = domain[2])[["root"]]
  }

  return(c("lowerBound" = lowerBound, "upperBound" = upperBound))
}

# Helpers: propDiff ----

# Predictable plug-in for the numerator of the eBeta test, for block i given
# the counts of blocks 1 to i - 1 only: the independent Beta posterior means
# of thetaA and thetaB.
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

# Cumulative log e-process of the grow test on propDiff ("greater"), block
# by block. The numerator lives on the curve thetaB - thetaA = propDiffMin:
# block i uses the posterior mean of thetaA under a grid posterior on that
# curve given blocks 1..i-1, the null is the pooled projection onto
# thetaA = thetaB, and block 1's UMP conditional e-factor multiplies the
# whole process. Returns logEValueVec of length nBlocks; with earlyStopping
# the loop stops at the first block whose value reaches log(1 / alpha), so
# the vector is shorter and its length is the stopping time.
logEProcess2x2PropDiffGrow <- function(ya, yb, na, nb, betaParameter,
                                       propDiffMin, alpha,
                                       earlyStopping = FALSE,
                                       nWeight = 1e3L) {
  nBlocks <- length(ya)
  logThreshold <- if (earlyStopping) log(1 / alpha) else Inf

  # Only the thetaA prior is used: rho = thetaA rescaled to (0, 1).
  betaA1 <- betaParameter[["betaA1"]]
  betaA2 <- betaParameter[["betaA2"]]

  # Grid on the curve; thetaA is restricted to (0, 1 - propDiffMin).
  rho <- seq(1 / nWeight, 1 - 1 / nWeight, length.out = nWeight)
  thetaAGrid <- rho * (1 - propDiffMin)
  thetaBGrid <- thetaAGrid + propDiffMin
  logThetaA <- log(thetaAGrid)
  logOneMinusThetaA <- log1p(-thetaAGrid)
  logThetaB <- log(thetaBGrid)
  logOneMinusThetaB <- log1p(-thetaBGrid)

  # Un-normalised log posterior weights, shifted so their maximum is 0: the
  # largest weight is then exactly 1 and the sum can neither underflow nor
  # overflow, however many blocks have been seen.
  logWeights <- (betaA1 - 1) * log(rho) + (betaA2 - 1) * log1p(-rho)
  logWeights <- logWeights - max(logWeights)

  # Block 1's UMP conditional e-factor, in the direction of the alternative.
  eValueUmp <- savi2x2TestStatUmp(
    ya[1], yb[1], na[1], nb[1], alpha, "greater"
  )

  # Cumulative log likelihood of blocks 1..i under the null (denominator)
  # and under the alternative (numerator).
  logLikelihoodNull <- 0
  logLikelihoodAlternative <- log(eValueUmp)
  logEValueVec <- numeric(0)

  for (i in seq_len(nBlocks)) {
    # Numerator: posterior mean of thetaA given blocks 1..i-1, on the curve.
    weights <- exp(logWeights)
    thetaA <- sum(thetaAGrid * weights) / sum(weights)
    thetaB <- thetaA + propDiffMin
    # Null: projection onto thetaA = thetaB, the size-weighted pooled mean.
    thetaNull <- (na[i] * thetaA + nb[i] * thetaB) / (na[i] + nb[i])

    logLikelihoodNull <- logLikelihoodNull +
      stats::dbinom(ya[i], na[i], thetaNull, log = TRUE) +
      stats::dbinom(yb[i], nb[i], thetaNull, log = TRUE)
    logLikelihoodAlternative <- logLikelihoodAlternative +
      stats::dbinom(ya[i], na[i], thetaA, log = TRUE) +
      stats::dbinom(yb[i], nb[i], thetaB, log = TRUE)
    logEValueVec[i] <- logLikelihoodAlternative - logLikelihoodNull

    if (logEValueVec[i] >= logThreshold) break

    # Only now add block i to the posterior, so block i + 1 is predicted from
    # the past alone.
    logWeights <- logWeights +
      ya[i] * logThetaA + (na[i] - ya[i]) * logOneMinusThetaA +
      yb[i] * logThetaB + (nb[i] - yb[i]) * logOneMinusThetaB
    logWeights <- logWeights - max(logWeights)
  }

  logEValueVec
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

# Helpers: logOdds ----

# Per-block conditional log likelihood at one logOdds (B minus A): given the
# block's total ya + yb, yb is Fisher's noncentral hypergeometric with odds
# exp(logOdds) on group B. A vector of length nBlocks; dFNCHypergeo() takes
# scalar sizes, hence the loop over blocks.
logLikelihoodFNCH <- function(ya, yb, na, nb, logOdds) {
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

# Sampling functions for design ----

#' Simulate stopping times of the propDiff grow test
#'
#' Decision 29. `propDiffMin` (`> 0`, `"greater"`) is both the grow plug-in
#' and the data-generating effect: data lie on `thetaB = thetaA +
#' propDiffMin` at `nTheta` baselines `thetaA`, and `nPlan` is the worst
#' `power` quantile of the stopping time over those baselines.
#'
#' @return A list: `thetaA`, `thetaB`, `stoppingTimes` (`nTheta x nSim`,
#'   `Inf` when a path never crosses `1 / alpha`), `nPlan`, `worstCaseIndex`.
#' @noRd
sampleStoppingTimesSavi2x2 <- function(
  propDiffMin, logOddsMin = NULL, na, nb, power, alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  eType = c("grow"),
  betaParameter = NULL, nTheta = 8L, nSim = 1e3L, nMax = 1e4L, nBoot = 1e4L,
  seed = NULL, wantEValuesAtNMax = FALSE,
  wantSamplePaths = TRUE, wantSimData = TRUE, pb = TRUE
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  # Only the grow test on propDiff, "greater", is planned for now.
  if (is.null(propDiffMin) || !is.null(logOddsMin)) {
    stop("sampleStoppingTimesSavi2x2 plans propDiffMin only for now")
  }
  if (alternative != "greater") {
    stop("sampleStoppingTimesSavi2x2 is designed for alternative = 'greater' only")
  }
  stopifnot(
    propDiffMin > 0, propDiffMin < 1, alpha > 0, alpha < 1,
    power > 0, power < 1, na >= 1, nb >= 1, is.finite(nMax)
  )

  if (is.null(betaParameter)) {
    betaParameter <- constructSaviDesignObj("Two Proportions")[["betaParameter"]]
  }
  # Reproducible by default: 2026 unless a seed is given.
  set.seed(if (is.null(seed)) 2026 else seed)

  # Baselines on the curve thetaB = thetaA + propDiffMin: thetaA runs over
  # its feasible range (0, 1 - propDiffMin) at nTheta equally spaced
  # interior points.
  rhoTheta <- seq(1 / (nTheta + 1), nTheta / (nTheta + 1), length.out = nTheta)
  thetaATrue <- rhoTheta * (1 - propDiffMin)
  thetaBTrue <- thetaATrue + propDiffMin

  logThreshold <- log(1 / alpha)
  stoppingTimes <- matrix(Inf, nrow = nTheta, ncol = nSim)

  if (pb) {
    pbSavi <- utils::txtProgressBar(style = 3, title = "Sampling worst-case stopping time")
  }

  for (k in seq_len(nTheta)) {
    for (sim in seq_len(nSim)) {
      if (pb) {
        utils::setTxtProgressBar(
          pbSavi, "value" = ((k - 1) * nSim + sim) / (nTheta * nSim), "title" = "Trials"
        )
      }

      ya <- stats::rbinom(nMax, na, thetaATrue[k])
      yb <- stats::rbinom(nMax, nb, thetaBTrue[k])

      # The test's own grow e-process, stopped at the first crossing of
      # 1 / alpha: the length of the returned vector is the stopping time,
      # unless the path ran through all nMax blocks without crossing.
      logEValueVec <- logEProcess2x2PropDiffGrow(
        ya, yb, rep(na, nMax), rep(nb, nMax), betaParameter, propDiffMin,
        alpha, earlyStopping = TRUE
      )
      if (logEValueVec[length(logEValueVec)] >= logThreshold) {
        stoppingTimes[k, sim] <- length(logEValueVec)
      }
    }
  }

  if (pb) close(pbSavi)

  # Planned block count: the power quantile of the stopping time at the
  # hardest baseline. type = 1 is an order statistic, so it is a realised
  # stopping time, finite exactly when at least a fraction power of that
  # baseline's paths crossed 1 / alpha within nMax (never-crossing paths are
  # Inf).
  quantiles <- apply(stoppingTimes, 1, stats::quantile, probs = power,
    names = FALSE, type = 1
  )
  worstCaseIndex <- which.max(quantiles)
  nPlan <- ceiling(quantiles[worstCaseIndex])

  if (!is.finite(nPlan)) {
    fractionNeverCrossed <- mean(!is.finite(stoppingTimes[worstCaseIndex, ]))
    warning(sprintf(paste(
      "the %g quantile of the stopping time is Inf: %.1f%% of the paths at",
      "thetaA = %.3f never cross 1/alpha at nMax = %g, try increasing nMax or propDiffMin"
    ), power, 100 * fractionNeverCrossed, thetaATrue[worstCaseIndex], nMax))
  }

  list(
    "thetaA" = thetaATrue, "thetaB" = thetaBTrue,
    "stoppingTimes" = stoppingTimes,
    "nPlan" = nPlan,
    "worstCaseIndex" = worstCaseIndex
  )
}
