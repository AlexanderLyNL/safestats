# Test functions ----

#' Computes Conditional E-Value with UMP Log Odds Ratio
savi2x2TestStatUmp <- function(ya, yb, na, nb, alpha,
                               alternative = c("twoSided", "greater", "less")) {
  alternative <- match.arg(alternative)
  logLikelihoodNull <- stats::dhyper(ya, na, nb, ya + yb, log = TRUE)

  # solve logOdds on each side with oneSided test
  sides <- switch(alternative,
    "twoSided" = c("greater", "less"),
    "greater" = "greater",
    "less" = "less"
  )
  eValue <- 0
  for (side in sides) {
    logOdds <- solveUmpLogOdds(na, nb, ya + yb, alpha, side)
    if (is.null(logOdds)) logOdds <- 0
    logLikelihoodAlternative <- logLikelihoodFNCH(ya, yb, na, nb, logOdds)
    eValue <- eValue + exp(logLikelihoodAlternative - logLikelihoodNull)
  }

  # twoSided: 1/2 less + 1/2 greater
  eValue / length(sides)
}

#' Safe anytime-valid 2x2 test for propDiff
#'
#' Tests `propDiff = 0` on a stream of 2x2 tables, one table per block, where
#' `propDiff = thetaB - thetaA` (B minus A). The numerator for block `i`
#' uses blocks `1..i-1` only. `designObj[["eType"]]` picks the e-process:
#'
#' - `"eBeta"`: two independent Beta posterior means for `thetaA` and
#'   `thetaB`, tested against the pooled null mean. It is `"twoSided"` only.
#'   The first table's e-value is replaced by the UMP conditional e-value.
#' - `"grow"`: `thetaA` and `thetaB` are restricted to
#'   `thetaB - thetaA = propDiffMin`. The alternative is `"greater"`, or
#'   `"twoSided"`, which averages the cumulative e-values at `+propDiffMin`
#'   and `-propDiffMin`. The first table's e-value is replaced by the UMP
#'   conditional e-value.
#'
#' Only `"eBeta"` gives a confidence interval or sequence.
#'
#' @param ya,yb Successes in groups A and B, one nonnegative integer per
#'   block, in observation order.
#' @param designObj A `saviDesign` for `propDiff` from `designSavi2x2()`.
#'   It supplies `alpha`, `alternative`, `eType`, `esMin`, `betaParameter`,
#'   `runningIntersection` and the group sizes per block `nPlan$na` and
#'   `nPlan$nb`. A scalar size is repeated for every block.
#' @param wantCi `TRUE` for one confidence interval on all blocks
#'   (`"eBeta"` only).
#' @param wantConfidenceSequence `TRUE` for a confidence sequence with one
#'   row per block (`"eBeta"` only). It takes precedence over `wantCi`.
#' @param ciValue Confidence level; the default `NULL` gives `1 - alpha`,
#'   with `alpha` from `designObj`.
#'
#' @return A `saviTest` list with `testName = "Two Proportions"`,
#'   `testType = "2x2"` and:
#'
#' - `eValueVec`: the realised cumulative e-value after each block, on
#'   blocks `1..i` for element `i` (not the blockwise e-factors).
#' - `eValue`: the last element of `eValueVec`, the e-value on all blocks.
#'   Reject when it is at least `1 / alpha`.
#' - `n`: `c(na = sum(na), nb = sum(nb), nBlocks = length(ya))`.
#' - `n1Vec`: `seq_len(nBlocks)`, the block index, used for plotting only.
#' - `estimate`: `c(thetaA, thetaB)`, the pooled observed proportions
#'   `sum(ya) / sum(na)` and `sum(yb) / sum(nb)`.
#' - `confSeq`: named `c(lowerBound, upperBound)` for `propDiff` at level
#'   `ciValue` on all blocks. With `wantConfidenceSequence`, this is the last
#'   row of `confSeqMatrix`. Both bounds are `NA` when the set is empty. It
#'   is `NULL` for `"grow"`, or when neither `wantCi` nor
#'   `wantConfidenceSequence` is `TRUE`.
#' - `confSeqMatrix`: an `nBlocks x 2` matrix of `lowerBound` and
#'   `upperBound`. Row `i` is the interval on blocks `1..i`. With
#'   `runningIntersection`, the rows are nested and stay `NA` after the first
#'   empty row. It is present only with `wantConfidenceSequence` and
#'   `"eBeta"`.
#' - `ciValue`: the confidence level used.
#' - `betaPrior`: `list(betaA1, betaA2, betaB1, betaB2)`, the design's Beta
#'   prior updated with all blocks. This is the prior that a next block would
#'   use.
#' - `alternative`, `h0`, `designObj`: copied from the design.
#' - `dataName`: the deparsed `ya` and `yb` arguments.
#' - `call`: the matched call.
#'
#' The constructor's other fields (`statistic`, `eValueApproxError`,
#' `note`) stay `NULL`.
#' @noRd
savi2x2TestStatPropDiff <- function(ya, yb,
                                    designObj = NULL, wantCi = TRUE,
                                    wantConfidenceSequence = FALSE, ciValue = NULL) {
  # grow must have a propDiffMin; eBeta must have none and be twoSided
  propDiffMin <- unname(designObj[["esMin"]])
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
    warnings("There is only 1 table, switched to UMP conditional e-variable")
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
    logEValueVec <- logEValueVec2x2PropDiffGrow(
      ya, yb, na, nb, betaParameter, propDiffMin, alpha, alternative
    )
  } else if (eType == "eBeta") {
    # Block 1's UMP conditional e-factor replaces its plug-in factor.
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
    # Replace block 1's plug-in factor by the UMP factor: blocks 2..i keep
    # their predictable factors, so the product stays an e-process.
    logEValueVec <- logLikelihoodAlternative - logLikelihoodNull
    logEValueVec <- logEValueVec - logEValueVec[1] + log(eValueUmp)
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
#' `eType = "grow"`: the fixed alternative `logOddsMin` (Decision 15), or
#' for `"twoSided"` the average of the processes at `+/-logOddsMin`
#' (Decision 38),
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
  logOddsMin <- unname(designObj[["esMin"]])
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
  logLikelihoodNull <- cumsum(stats::dhyper(ya, na, nb, ya + yb, log = TRUE))
  logLikelihoodAlternative <- switch(eType,
    # grow: the fixed alternative logOddsMin (Decision 15); twoSided is
    # the log of the plain average of the cumulative likelihoods at
    # +logOddsMin and -logOddsMin, shifted by their max (Decision 38).
    grow = {
      logPlus <- cumsum(logLikelihoodFNCH(ya, yb, na, nb, logOddsMin))
      if (alternative == "twoSided") {
        logMinus <- cumsum(logLikelihoodFNCH(ya, yb, na, nb, -logOddsMin))
        shift <- pmax(logPlus, logMinus)
        shift + log(0.5 * (exp(logPlus - shift) + exp(logMinus - shift)))
      } else {
        logPlus
      }
    },
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
      # Cumulate over blocks (cumsum down each grid column; matrix() keeps
      # a single block as a 1-row matrix), then add the log prior weight to
      # every row of its column: sweep(x, 2, v, "+") adds v[j] to column j.
      # Row i then holds log(prior * likelihood of blocks 1..i) on the grid,
      # mixed over the grid by the log-sum-exp below.
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
    domain <- c(-40, 40) # TODO: maybe add bounds print or else
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
# TODO: add h0 != 0 situation

#' Design a safe anytime-valid 2x2 test
#'
#' `eType` picks the e-variable and, with it, the effect measure. Both
#' effects are signed B minus A and anchored on `thetaA`: `propDiff =
#' thetaB - thetaA` and `logOdds = logit(thetaB) - logit(thetaA)`.
#'
#' - `"eBeta"` (propDiff) and `"eGauss"` (logOdds) are unrestricted and
#'   twoSided only; they take no `propDiffMin` or `logOddsMin` and no
#'   planning, and a one-sided `alternative` is ignored with a warning.
#' - `"grow"` plugs in exactly one of `propDiffMin`, `logOddsMin` as the
#'   fixed alternative (Decisions 15, 34, 38), `"greater"` or `"twoSided"`.
#'   Planning exists for `propDiff` only (Decisions 37, 43): `power` alone
#'   plans the block count at the hardest baseline, `nBlocksPlan` alone
#'   evaluates the worst-case power there, and both without `propDiffMin`
#'   find the minimal detectable `propDiff` (Decision 42). Both with
#'   `propDiffMin` errors; neither gives the design without simulation.
#'   `logOddsMin` with `power` or `nBlocksPlan` errors: its worst case is
#'   set by the baseline grid, not by the effect.
#'
#' @param na,nb Planned group sizes per block, one positive integer each.
#' @param nBlocksPlan Planned block count at which the worst-case power, or
#'   with `power` the minimal `propDiff`, is evaluated (`"grow"` on
#'   `propDiff` only).
#' @param propDiffMin,logOddsMin Minimal effect for `"grow"`, at most one of
#'   them: `propDiffMin` strictly inside `(0, 1)`, `logOddsMin` finite and
#'   `> 0`. Stored as `esMin`; `propDiffMin` is found from `power` and
#'   `nBlocksPlan` when both `*Min` are `NULL`.
#' @param alpha Significance level; the test rejects at `1 / alpha`.
#' @param power Target power (`"grow"` on `propDiff` only). Plans
#'   `nBlocksPlan` when that is `NULL`, the minimal `propDiff` when it is
#'   given.
#' @param h0 The null value of `propDiff`; only `0` is designed.
#' @param alternative `"twoSided"` or `"greater"`; `"less"` is not designed
#'   yet.
#' @param eType `"eBeta"`, `"grow"` or `"eGauss"`, see Details.
#' @param betaParameter `list(betaA1, betaA2, betaB1, betaB2)`, the Beta
#'   prior shapes on `thetaA` and `thetaB`; `NULL` keeps the constructor's
#'   default of `0.18` each.
#' @param runningIntersection `TRUE` to intersect each row of the blockwise
#'   confidence sequence with the previous one; `NULL` keeps the
#'   constructor's `FALSE`.
#' @param nTheta,nSim,nBoot,nMax,seed,wantSamplePaths,pb Simulation settings
#'   passed to [sampleStoppingTimesSavi2x2()] via [computeNPlanSavi2x2()] or
#'   [computePowerSavi2x2()]: baselines per curve, paths per baseline,
#'   bootstrap resamples, block cap per path, seed (`NULL` is `2026`),
#'   whether to keep the e-value paths, and the progress bar.
#'
#' @return A `saviDesign` with `testName = "Two Proportions"`, `testType =
#'   "2x2"`, `h0 = c(propDiff = h0)`, `esMin`, `eType`, `alpha`,
#'   `alternative`, `betaParameter`, `parameter` (the prior summarised for
#'   printing), `runningIntersection` and `nPlan = list(na, nb)`, with a
#'   third element `nBlocksPlan` when planned or given. With `power`:
#'   `designScenario = "1a"`, `power` as the target, `nPlanTwoSe = c(NA,
#'   NA, 2 * bootSe)`, `bootObjNBlocksPlan`, `nMean`, `nMeanTwoSe`,
#'   `bootObjNMean`. With `nBlocksPlan`: `designScenario = "2"`, `power`
#'   (the worst case), `powerTwoSe`, `bootObjPower`. Both also carry
#'   `worstCaseIndex`, `worstCaseThetaA`, `worstCaseThetaB`, `breakVector`
#'   and `samplePaths`. With `power` and `nBlocksPlan` but no `*Min`:
#'   `designScenario = "3"`, `esMin` the minimal detectable `propDiff`,
#'   `power` as the target, and no simulation summaries.
#' @noRd
designSavi2x2 <- function(
  na, nb, nBlocksPlan = NULL,
  propDiffMin = NULL, logOddsMin = NULL,
  alpha = 0.05, power = NULL, h0 = 0,
  alternative = c("twoSided", "greater"),
  eType = c("eBeta", "grow", "eGauss"),
  betaParameter = NULL,
  runningIntersection = NULL,
  nTheta = 8L, nSim = 1e3L, nBoot = nSim, nMax = 1e4L, seed = NULL,
  wantSamplePaths = FALSE, pb = TRUE
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

  # Named so print() and plot() show which effect the minimal value is on.
  result[["esMin"]] <- if (!is.null(propDiffMin)) {
    c("propDiff" = propDiffMin)
  } else if (!is.null(logOddsMin)) {
    c("logOdds" = logOddsMin)
  }

  # Dispatch (Decisions 37, 42, 43). grow on propDiff is the only case with
  # sampling. With propDiffMin: power alone plans the block count at the
  # hardest baseline, nBlocksPlan alone evaluates the worst-case power
  # there, both is contradictory. Without a *Min: power and nBlocksPlan
  # together find the minimal detectable propDiff. logOddsMin takes no
  # planning: its worst case over the baselines is set by how close the
  # outermost baseline sits to 0 or 1, where the pair approaches the null,
  # not by logOddsMin. eBeta and eGauss take no *Min and no planning, and
  # are twoSided only.
  if (!is.null(propDiffMin) && !is.null(logOddsMin)) {
    stop("supply propDiffMin or logOddsMin, not both")
  }
  wantEsMin <- FALSE
  if (eType == "grow") {
    if (!is.null(logOddsMin) && (!is.null(power) || !is.null(nBlocksPlan))) {
      stop("no planning on logOdds, plan with propDiffMin")
    }
    if (is.null(propDiffMin) && is.null(logOddsMin)) {
      if (is.null(power) || is.null(nBlocksPlan)) {
        stop("eType = 'grow' needs propDiffMin or logOddsMin, or both power and nBlocksPlan to find the minimal propDiff")
      }
      wantEsMin <- TRUE
    } else if (!is.null(power) && !is.null(nBlocksPlan)) {
      stop("supply power (to find nBlocksPlan) or nBlocksPlan (to find power), not both")
    }
  } else {
    if (!is.null(propDiffMin) || !is.null(logOddsMin)) {
      stop("eType = '", eType, "' takes no propDiffMin or logOddsMin")
    }
    if (!is.null(power) || !is.null(nBlocksPlan)) {
      stop("no sampling for eType = '", eType, "'; power and nBlocksPlan need eType = 'grow'")
    }
    if (alternative != "twoSided") {
      warning("eType = '", eType, "' is twoSided; alternative = '", alternative, "' is ignored")
    }
  }
  if (wantEsMin) {
    esMin <- computeEsMinSavi2x2(
      na = na, nb = nb, nBlocksPlan = nBlocksPlan, power = power,
      alpha = alpha, alternative = alternative,
      betaParameter = result[["betaParameter"]], nTheta = nTheta,
      nSim = nSim, seed = seed, pb = pb
    )
    # NA: the worst-case power never reaches the target on the search
    # bounds, so no minimal effect can be reported.
    if (is.na(esMin)) {
      stop(sprintf(paste(
        "no minimal propDiff found: at nBlocksPlan = %g the worst-case power",
        "does not reach %g for any value in (0.01, 0.9); try a larger",
        "nBlocksPlan or a smaller power"
      ), nBlocksPlan, power))
    }
    result[["designScenario"]] <- "3"
    result[["esMin"]] <- c("propDiff" = esMin)
    result[["power"]] <- power
  } else if (!is.null(power)) {
    planning <- computeNPlanSavi2x2(
      propDiffMin = propDiffMin, na = na, nb = nb,
      power = power, alpha = alpha, alternative = alternative,
      betaParameter = result[["betaParameter"]], nTheta = nTheta,
      nSim = nSim, nBoot = nBoot, nMax = nMax, seed = seed,
      wantSamplePaths = wantSamplePaths, pb = pb
    )
    result[["designScenario"]] <- "1a"
    result[["power"]] <- power
    nBlocksPlan <- planning[["nPlan"]]
    # na, nb are planned, not simulated: no standard error for them.
    result[["nPlanTwoSe"]] <- c(NA, NA, 2 * planning[["bootObjNPlan"]][["bootSe"]])
    result[["bootObjNBlocksPlan"]] <- planning[["bootObjNPlan"]]
    result[["nMean"]] <- c("nMean" = planning[["nMean"]])
    result[["nMeanTwoSe"]] <- 2 * planning[["bootObjNMean"]][["bootSe"]]
    result[["bootObjNMean"]] <- planning[["bootObjNMean"]]
  } else if (!is.null(nBlocksPlan)) {
    planning <- computePowerSavi2x2(
      propDiffMin = propDiffMin, na = na, nb = nb,
      nBlocks = nBlocksPlan, alpha = alpha, alternative = alternative,
      betaParameter = result[["betaParameter"]], nTheta = nTheta,
      nSim = nSim, nBoot = nBoot, seed = seed,
      wantSamplePaths = wantSamplePaths, pb = pb
    )
    result[["designScenario"]] <- "2"
    result[["power"]] <- planning[["power"]]
    result[["powerTwoSe"]] <- 2 * planning[["bootObjPower"]][["bootSe"]]
    result[["bootObjPower"]] <- planning[["bootObjPower"]]
  }
  if (!wantEsMin && (!is.null(power) || !is.null(nBlocksPlan))) {
    worstCaseIndex <- planning[["worstCaseIndex"]]
    result[["worstCaseIndex"]] <- worstCaseIndex
    result[["worstCaseThetaA"]] <- planning[["thetaA"]][worstCaseIndex]
    result[["worstCaseThetaB"]] <- planning[["thetaB"]][worstCaseIndex]
    result[["breakVector"]] <- planning[["breakVector"]]
    result[["samplePaths"]] <- planning[["samplePaths"]]
  }
  result[["parameter"]] <- c(
    "Beta hyperparameters" =
      paste(unlist(result[["betaParameter"]]), collapse = " ")
  )
  result[["eType"]] <- eType
  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["h0"]] <- c("propDiff" = h0)
  result[["nPlan"]] <- list("na" = na, "nb" = nb)
  if (!is.null(nBlocksPlan)) {
    result[["nPlan"]][["nBlocksPlan"]] <- unname(nBlocksPlan)
  }
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

    # Null thetaA for every active candidate, one uniroot call each. Could
    # be replaced by a vectorised bisection over all candidates at once (the
    # KL derivative times x (1 - x)(x + delta)(1 - x - delta) is a cubic in
    # x with a single root on the feasible interval), which is 6 to 15
    # times faster; kept as is for readability.
    nullThetaA <- mapply(solveRIPr2x2PropDiff, delta,
      MoreArgs = list(thetaA = thetaA[i], thetaB = thetaB[i],
        na = na[i], nb = nb[i])
    )

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

# Cumulative log e-process of the grow test on propDiff, block by block.
# The numerator lives on the curve thetaB = thetaA + propDiffMin: block i
# uses the posterior mean of thetaA under a grid posterior on that curve
# given blocks 1..i-1, the null is the pooled projection onto
# thetaA = thetaB, and block 1's plug-in ratio is replaced by its UMP
# conditional e-factor (the posterior still absorbs block 1). "twoSided" runs the same process on the curve
# thetaB = thetaA - propDiffMin as well and averages the two cumulative
# e-values as processes (Decision 34). Returns logEValueVec of length
# nBlocks; with earlyStopping the loop stops at the first block whose value
# reaches log(1 / alpha), so the vector is shorter and its length is the
# stopping time.
logEValueVec2x2PropDiffGrow <- function(ya, yb, na, nb, betaParameter,
                                       propDiffMin, alpha,
                                       alternative = c("twoSided", "greater", "less"),
                                       earlyStopping = FALSE,
                                       nWeight = 1e3L) {
  alternative <- match.arg(alternative)
  if (alternative == "less") {
    stop("alternative = 'less' is not designed yet for the grow test on propDiff")
  }
  nBlocks <- length(ya)
  logThreshold <- log(1 / alpha)
  # The design requires propDiffMin > 0; abs() is for a later signed value.
  propDiffMin <- abs(propDiffMin)

  # One-sided processes to run: the plus curve, and for "twoSided" the minus
  # curve as well. Each side keeps its own grid posterior and likelihoods.
  signs <- if (alternative == "twoSided") c(1, -1) else 1
  nSides <- length(signs)

  # Only the thetaA prior is used: rho = thetaA rescaled to (0, 1).
  betaA1 <- betaParameter[["betaA1"]]
  betaA2 <- betaParameter[["betaA2"]]

  # Grid on each curve, one column per side: thetaA is restricted to
  # (0, 1 - propDiffMin) on the plus curve and (propDiffMin, 1) on the
  # minus curve.
  rho <- seq(1 / nWeight, 1 - 1 / nWeight, length.out = nWeight)
  thetaAGrid <- matrix(rho * (1 - propDiffMin), nWeight, nSides)
  thetaAGrid[, signs < 0] <- propDiffMin + thetaAGrid[, signs < 0]
  thetaBGrid <- thetaAGrid + rep(signs * propDiffMin, each = nWeight)
  logThetaA <- log(thetaAGrid)
  logOneMinusThetaA <- log1p(-thetaAGrid)
  logThetaB <- log(thetaBGrid)
  logOneMinusThetaB <- log1p(-thetaBGrid)

  # Un-normalised log posterior weights, shifted so their maximum is 0: the
  # largest weight is then exactly 1 and the sum can neither underflow nor
  # overflow, however many blocks have been seen. The prior on rho is the
  # same on both curves.
  logWeights <- (betaA1 - 1) * log(rho) + (betaA2 - 1) * log1p(-rho)
  logWeights <- matrix(logWeights - max(logWeights), nWeight, nSides)

  # Block 1's UMP conditional e-factor per side, in that side's direction.
  # It depends on block 1 alone and replaces block 1's plug-in ratio after
  # the loop; inside the loop it is used only by the early-stopping check.
  logEValueUmp <- numeric(nSides)
  for (s in seq_len(nSides)) {
    logEValueUmp[s] <- log(savi2x2TestStatUmp(
      ya[1], yb[1], na[1], nb[1], alpha,
      if (signs[s] > 0) "greater" else "less"
    ))
  }

  # The plain plug-in process, Bayesian updating from start to finish:
  # every block, block 1 included, takes the posterior-mean plug-in given
  # blocks 1..i-1 and is added to the posterior afterwards. Row i holds
  # the cumulative log e-value of blocks 1..i per side.
  logLikelihoodNull <- numeric(nSides)
  logLikelihoodAlternative <- numeric(nSides)
  logEValueSides <- matrix(NA_real_, nBlocks, nSides)

  for (i in seq_len(nBlocks)) {
    for (s in seq_len(nSides)) {
      # Numerator: posterior mean of thetaA given blocks 1..i-1, on the curve.
      weights <- exp(logWeights[, s])
      thetaA <- sum(thetaAGrid[, s] * weights) / sum(weights)
      thetaB <- thetaA + signs[s] * propDiffMin
      # Null: projection onto thetaA = thetaB, the size-weighted pooled mean.
      thetaNull <- (na[i] * thetaA + nb[i] * thetaB) / (na[i] + nb[i])

      logLikelihoodNull[s] <- logLikelihoodNull[s] +
        stats::dbinom(ya[i], na[i], thetaNull, log = TRUE) +
        stats::dbinom(yb[i], nb[i], thetaNull, log = TRUE)
      logLikelihoodAlternative[s] <- logLikelihoodAlternative[s] +
        stats::dbinom(ya[i], na[i], thetaA, log = TRUE) +
        stats::dbinom(yb[i], nb[i], thetaB, log = TRUE)
    }
    logEValueSides[i, ] <- logLikelihoodAlternative - logLikelihoodNull

    # Early stopping looks at the replaced process (block 1 swapped for the
    # UMP factor, see below), averaged over the sides on the log scale.
    if (earlyStopping) {
      logEValueReplaced <- logEValueSides[i, ] - logEValueSides[1, ] + logEValueUmp
      logEValueMax <- max(logEValueReplaced)
      if (logEValueMax + log(mean(exp(logEValueReplaced - logEValueMax))) >= logThreshold) {
        logEValueSides <- logEValueSides[seq_len(i), , drop = FALSE]
        break
      }
    }

    # Only now add block i to each posterior, so block i + 1 is predicted
    # from the past alone.
    logWeights <- logWeights +
      ya[i] * logThetaA + (na[i] - ya[i]) * logOneMinusThetaA +
      yb[i] * logThetaB + (nb[i] - yb[i]) * logOneMinusThetaB
    # Re-centre each side's column at its maximum (sweep(x, 2, v) subtracts
    # v[j] from column j): the largest weight is again exactly 1, so the
    # posterior mean's exp() and sum() cannot underflow however many
    # blocks have been seen. A per-column constant leaves the mean unchanged.
    logWeights <- sweep(logWeights, 2, apply(logWeights, 2, max))
  }

  # Replace block 1's plug-in ratio by its UMP factor: the cumulative
  # process of each side shifts by the same constant from block 1 on, so
  # replacing logEValueSides[1, ] alone would leave the later rows wrong.
  # sweep(x, 2, v, "+") adds v[j], side j's log UMP factor minus its block-1
  # plug-in log ratio, to every row of column j; row 1 becomes the log UMP
  # factor itself.
  logEValueSides <- sweep(logEValueSides, 2, logEValueUmp - logEValueSides[1, ], "+")

  # Average of the one-sided cumulative e-values per block, on the log
  # scale shifted by the larger one so exp() cannot overflow. One side:
  # the value itself.
  logEValueMax <- apply(logEValueSides, 1, max)
  logEValueVec <- logEValueMax +
    log(rowMeans(exp(logEValueSides - logEValueMax)))

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

# Per-block conditional log likelihood of ya at one logOdds (A minus B):
# given the block's total ya + yb, ya is Fisher's noncentral hypergeometric
# with odds exp(logOdds) on group A. A vector of length nBlocks;
# dFNCHypergeo() takes scalar sizes, hence the loop over blocks.
logLikelihoodFNCH <- function(ya, yb, na, nb, logOdds) {
  mapply(function(ya, yb, na, nb) {
    log(BiasedUrn::dFNCHypergeo(ya, na, nb, ya + yb, exp(logOdds)))
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
#
# Returns the root, or NULL when no logOdds on that side reaches the target.
# The KL is 0 at nullLogOdds and increases away from it towards its
# supremum -log P0(yb at its feasible extreme): yb = min(nb, totalSuccesses)
# for "greater", yb = max(0, totalSuccesses - na) for "less". NULL is
# returned exactly when that supremum is at most log(1 / alpha), i.e. when
# even the most extreme table under this total has null probability at
# least alpha, so no one-block test at level alpha can reject. Typical cases
# are totalSuccesses = 0 or na + nb (one feasible table, KL = 0), and small
# blocks: na = nb = 1 gives supremum log(2). The caller then uses the
# trivial e-factor 1; this is by design, not a numerical failure.
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

#' Simulate stopping times of the 2x2 grow test
#'
#' Decisions 29, 34, 43. `propDiffMin` (`> 0`) is the grow plug-in and the
#' data-generating effect: data lie on `thetaB = thetaA + propDiffMin` at
#' `nTheta` baselines `thetaA`, and for `"twoSided"` on `thetaB = thetaA -
#' propDiffMin` at `nTheta` more. There is no planning on `logOdds`
#' (Decision 43).
#' `nPlan` is the worst `power` quantile of the stopping time over all
#' baselines; with `power = NULL` the quantile step is skipped and `nPlan`,
#' `worstCaseIndex` are `NULL` (Decision 36).
#'
#' @return A list: `thetaA`, `thetaB`, and one row per baseline in
#'   `stoppingTimes` (`Inf` when a path never crosses `1 / alpha`),
#'   `breakVector` (`0` crossed, `1` reached `nMax`), `eValuesStopped`;
#'   `samplePaths` (a list of `nSim x nMax` sparse matrices, or `NULL`),
#'   `n1Vector` (the block index), `nPlan`,
#'   `worstCaseIndex`.
#' @noRd
sampleStoppingTimesSavi2x2 <- function(
  propDiffMin, na, nb, power = NULL, alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  eType = c("grow"),
  betaParameter = NULL, nTheta = 8L, nSim = 1e3L, nMax = 1e4L, nBoot = 1e4L,
  seed = NULL, wantEValuesAtNMax = FALSE,
  wantSamplePaths = FALSE, wantSimData = TRUE, pb = TRUE
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  # "less" is not designed yet.
  if (alternative == "less") {
    stop("alternative = 'less' is not designed yet!")
  }
  stopifnot(
    length(propDiffMin) == 1, propDiffMin > 0, propDiffMin < 1,
    alpha > 0, alpha < 1, is.null(power) || (power > 0 && power < 1),
    na >= 1, nb >= 1, is.finite(nMax)
  )

  if (is.null(betaParameter)) {
    betaParameter <- constructSaviDesignObj("Two Proportions")[["betaParameter"]]
  }
  # Reproducible by default: 2026 unless a seed is given.
  set.seed(if (is.null(seed)) 2026 else seed)

  # TODO: a lot of time wasted near the boundary
  # Baselines: thetaA at nTheta equally spaced interior points of its
  # feasible range on the curve thetaB = thetaA + propDiffMin, thetaA over
  # (0, 1 - propDiffMin); twoSided adds the curve thetaB = thetaA -
  # propDiffMin, thetaA over (propDiffMin, 1), since the test is not
  # symmetric under a group swap when na != nb or the Beta priors differ.
  rhoTheta <- seq(1 / (nTheta + 1), nTheta / (nTheta + 1), length.out = nTheta)
  thetaATrue <- rhoTheta * (1 - propDiffMin)
  thetaBTrue <- thetaATrue + propDiffMin
  if (alternative == "twoSided") {
    thetaATrue <- c(thetaATrue, propDiffMin + rhoTheta * (1 - propDiffMin))
    thetaBTrue <- c(thetaBTrue, thetaATrue[-seq_len(nTheta)] - propDiffMin)
  }
  nBaselines <- length(thetaATrue)

  logThreshold <- log(1 / alpha)
  naVec <- rep(na, nMax)
  nbVec <- rep(nb, nMax)
  stoppingTimes <- matrix(Inf, nrow = nBaselines, ncol = nSim)
  breakVector <- matrix(1L, nrow = nBaselines, ncol = nSim)
  eValuesStopped <- matrix(NA_real_, nrow = nBaselines, ncol = nSim)
  samplePaths <- if (wantSamplePaths) vector("list", nBaselines) else NULL

  if (pb) {
    pbSavi <- utils::txtProgressBar(style = 3, title = "Sampling worst-case stopping time")
  }

  for (k in seq_len(nBaselines)) {
    if (wantSamplePaths) {
      samplePaths[[k]] <- matrix(0, nrow = nSim, ncol = nMax)
    }
    for (sim in seq_len(nSim)) {
      if (pb) {
        utils::setTxtProgressBar(
          pbSavi, "value" = ((k - 1) * nSim + sim) / (nBaselines * nSim), "title" = "Trials"
        )
      }

      # The test's own grow e-process, cut at the first crossing of
      # 1 / alpha: the length of the vector is the stopping time, unless the
      # path ran through all nMax blocks without crossing.
      ya <- stats::rbinom(nMax, na, thetaATrue[k])
      yb <- stats::rbinom(nMax, nb, thetaBTrue[k])
      logEValueVec <- logEValueVec2x2PropDiffGrow(
        ya, yb, naVec, nbVec, betaParameter, propDiffMin,
        alpha, alternative, earlyStopping = TRUE
      )

      nStopped <- length(logEValueVec)
      eValuesStopped[k, sim] <- exp(logEValueVec[nStopped])
      if (logEValueVec[nStopped] >= logThreshold) {
        stoppingTimes[k, sim] <- nStopped
        breakVector[k, sim] <- 0L
      }
      if (wantSamplePaths) {
        # The path up to the crossing, then the crossed value held to nMax.
        samplePaths[[k]][sim, ] <-
          exp(logEValueVec)[pmin(seq_len(nMax), nStopped)]
      }
    }
    if (wantSamplePaths) {
      samplePaths[[k]] <- Matrix::Matrix(samplePaths[[k]], sparse = TRUE)
    }
  }

  if (pb) close(pbSavi)

  # Planned block count: the power quantile of the stopping time at the
  # hardest baseline. type = 1 is an order statistic, so it is a realised
  # stopping time, finite exactly when at least a fraction power of that
  # baseline's paths crossed 1 / alpha within nMax (never-crossing paths are
  # Inf).
  nPlan <- NULL
  worstCaseIndex <- NULL
  if (!is.null(power)) {
    quantiles <- apply(stoppingTimes, 1, stats::quantile, probs = power,
      names = FALSE, type = 1
    )
    worstCaseIndex <- which.max(quantiles)
    nPlan <- ceiling(quantiles[worstCaseIndex])
  }

  if (!is.null(nPlan) && !is.finite(nPlan)) {
    fractionNeverCrossed <- mean(!is.finite(stoppingTimes[worstCaseIndex, ]))
    warning(sprintf(paste(
      "the %g quantile of the stopping time is Inf: %.1f%% of the paths at",
      "thetaA = %.3f never cross 1/alpha at nMax = %g, try increasing nMax",
      "or propDiffMin"
    ), power, 100 * fractionNeverCrossed, thetaATrue[worstCaseIndex], nMax))
  }

  list(
    "thetaA" = thetaATrue, "thetaB" = thetaBTrue,
    "stoppingTimes" = stoppingTimes,
    "breakVector" = breakVector,
    "eValuesStopped" = eValuesStopped,
    "samplePaths" = samplePaths,
    "n1Vector" = seq_len(nMax),
    "nPlan" = nPlan,
    "worstCaseIndex" = worstCaseIndex
  )
}


#' Worst-case power of the 2x2 grow test at a planned block count
#'
#' Decision 36. Runs [sampleStoppingTimesSavi2x2()] with `nMax = nBlocks`
#' and reports, per baseline, the fraction of paths that cross `1 / alpha`
#' within `nBlocks`; the worst case is the smallest of these, since the
#' test must meet its target whatever `thetaA` is.
#'
#' @param nBlocks Planned block count at which the test is evaluated.
#' @inheritParams sampleStoppingTimesSavi2x2
#'
#' @return A list: `power` (the worst-case power), `powerVec` (one per
#'   baseline), `worstCaseIndex`, `bootObjPower` (a [boot::boot()] object
#'   on the worst baseline, with `bootSe`), `nBlocks`, and the sampler's
#'   `thetaA`, `thetaB`, `stoppingTimes`, `breakVector`, `eValuesStopped`,
#'   `samplePaths`, `n1Vector`.
#' @noRd
computePowerSavi2x2 <- function(
  propDiffMin, na, nb, nBlocks, alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  betaParameter = NULL, nTheta = 8L, nSim = 1e3L, nBoot = nSim,
  seed = NULL, wantSamplePaths = FALSE, pb = TRUE
) {
  alternative <- match.arg(alternative)
  stopifnot(length(nBlocks) == 1, is.finite(nBlocks), nBlocks >= 1)

  samplingResult <- sampleStoppingTimesSavi2x2(
    propDiffMin = propDiffMin, na = na, nb = nb,
    power = NULL, alpha = alpha, alternative = alternative,
    betaParameter = betaParameter, nTheta = nTheta, nSim = nSim,
    nMax = nBlocks, seed = seed, wantSamplePaths = wantSamplePaths, pb = pb
  )

  # Power per baseline: the fraction of paths that crossed 1 / alpha within
  # nBlocks (a never-crossing path has stopping time Inf). The worst case is
  # the smallest.
  stoppingTimes <- samplingResult[["stoppingTimes"]]
  powerVec <- rowMeans(stoppingTimes <= nBlocks)
  worstCaseIndex <- which.min(powerVec)

  bootObjPower <- computeBootObj(
    values = stoppingTimes[worstCaseIndex, ], objType = "power",
    nPlan = nBlocks, nBoot = nBoot
  )

  list(
    "power" = powerVec[worstCaseIndex],
    "powerVec" = powerVec,
    "worstCaseIndex" = worstCaseIndex,
    "bootObjPower" = bootObjPower,
    "nBlocks" = nBlocks,
    "thetaA" = samplingResult[["thetaA"]],
    "thetaB" = samplingResult[["thetaB"]],
    "stoppingTimes" = stoppingTimes,
    "breakVector" = samplingResult[["breakVector"]],
    "eValuesStopped" = samplingResult[["eValuesStopped"]],
    "samplePaths" = samplingResult[["samplePaths"]],
    "n1Vector" = samplingResult[["n1Vector"]]
  )
}


#' Worst-case planned block count of the 2x2 grow test
#'
#' Decision 36. Runs [sampleStoppingTimesSavi2x2()] and reports its `nPlan`,
#' the `power` quantile of the stopping time at the hardest baseline
#' (Decision 33), with a bootstrap SE and the mean stopping time of paths
#' capped at `nPlan`, both at that baseline.
#'
#' @inheritParams sampleStoppingTimesSavi2x2
#'
#' @return A list: `nPlan` (the worst-case block count, `Inf` with a warning
#'   when the worst baseline crossed too rarely), `nPlanVec` (the quantile
#'   per baseline), `worstCaseIndex`, `bootObjNPlan`, `nMean`,
#'   `bootObjNMean` ([boot::boot()] objects on the worst baseline, `NULL`
#'   when `nPlan` is `Inf`), and the sampler's `thetaA`, `thetaB`,
#'   `stoppingTimes`, `breakVector`, `eValuesStopped`, `samplePaths`,
#'   `n1Vector`.
#' @noRd
computeNPlanSavi2x2 <- function(
  propDiffMin, na, nb, power = 0.8, alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  betaParameter = NULL, nTheta = 8L, nSim = 1e3L, nBoot = nSim,
  nMax = 1e4L, seed = NULL, wantSamplePaths = FALSE, pb = TRUE
) {
  alternative <- match.arg(alternative)
  stopifnot(!is.null(power), power > 0, power < 1)

  samplingResult <- sampleStoppingTimesSavi2x2(
    propDiffMin = propDiffMin, na = na, nb = nb,
    power = power, alpha = alpha, alternative = alternative,
    betaParameter = betaParameter, nTheta = nTheta, nSim = nSim,
    nMax = nMax, seed = seed, wantSamplePaths = wantSamplePaths, pb = pb
  )

  stoppingTimes <- samplingResult[["stoppingTimes"]]
  nPlan <- samplingResult[["nPlan"]]
  worstCaseIndex <- samplingResult[["worstCaseIndex"]]
  # The same order-statistic quantile per baseline as the sampler's nPlan.
  nPlanVec <- apply(stoppingTimes, 1, stats::quantile, probs = power,
    names = FALSE, type = 1
  )

  # Simulation uncertainty at the worst baseline only: the bootstrap
  # quantile, and the mean stopping time with paths capped at nPlan. A
  # never-crossing path (Inf) has no finite bootstrap, so both stay NULL.
  bootObjNPlan <- NULL
  bootObjNMean <- NULL
  nMean <- NA_real_
  if (is.finite(nPlan)) {
    bootObjNPlan <- computeBootObj(
      values = stoppingTimes[worstCaseIndex, ], objType = "nPlan",
      power = power, nBoot = nBoot
    )
    bootObjNMean <- computeBootObj(
      values = stoppingTimes[worstCaseIndex, ], objType = "nMean",
      nPlan = nPlan, nBoot = nBoot
    )
    nMean <- ceiling(bootObjNMean[["t0"]])
  }

  list(
    "nPlan" = nPlan,
    "nPlanVec" = nPlanVec,
    "worstCaseIndex" = worstCaseIndex,
    "bootObjNPlan" = bootObjNPlan,
    "nMean" = nMean,
    "bootObjNMean" = bootObjNMean,
    "thetaA" = samplingResult[["thetaA"]],
    "thetaB" = samplingResult[["thetaB"]],
    "stoppingTimes" = stoppingTimes,
    "breakVector" = samplingResult[["breakVector"]],
    "eValuesStopped" = samplingResult[["eValuesStopped"]],
    "samplePaths" = samplingResult[["samplePaths"]],
    "n1Vector" = samplingResult[["n1Vector"]]
  )
}


#' Minimal detectable propDiff of the 2x2 grow test
#'
#' Decisions 40, 43. The smallest `propDiffMin` at which the worst-case
#' power of [computePowerSavi2x2()] at `nBlocks = nBlocksPlan` reaches
#' `power`, found by [stats::uniroot()] on `propDiffBounds`. Every candidate
#' is simulated with the same `seed`, so the target is a deterministic step
#' function of the candidate. There is no planning on `logOdds`.
#'
#' @param nBlocksPlan Planned block count at which the test is evaluated.
#' @param propDiffBounds Search interval for `propDiffMin`, strictly inside
#'   `(0, 1)`.
#' @param tol Tolerance of the root on the `propDiff` scale.
#' @inheritParams sampleStoppingTimesSavi2x2
#'
#' @return A single numeric: the minimal `propDiffMin`, or `NA` when the
#'   worst-case power minus `power` has no sign change on `propDiffBounds`
#'   (still below the target at the upper bound, or already above it at the
#'   lower bound). No bootstrap object, as for [computeMinEsBatchSaviT()].
#' @noRd
computeEsMinSavi2x2 <- function(
  na, nb, nBlocksPlan, power = 0.8, alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  betaParameter = NULL, nTheta = 8L, nSim = 1e3L, seed = NULL, pb = TRUE,
  propDiffBounds = c(0.01, 0.9), tol = 1e-5
) {
  alternative <- match.arg(alternative)
  bounds <- propDiffBounds
  stopifnot(
    length(nBlocksPlan) == 1, is.finite(nBlocksPlan), nBlocksPlan >= 1,
    power > 0, power < 1, length(bounds) == 2, bounds[1] > 0,
    bounds[1] < bounds[2], bounds[2] < 1
  )

  # Worst-case power minus the target, at a candidate propDiffMin. The same
  # seed for every candidate makes this deterministic in the candidate; the
  # baselines are rescaled to (0, 1 - propDiffMin) inside the sampler, so
  # the worst case is taken afresh each time.
  targetFunction <- function(propDiffMin) {
    computePowerSavi2x2(
      propDiffMin = propDiffMin, na = na, nb = nb, nBlocks = nBlocksPlan,
      alpha = alpha, alternative = alternative,
      betaParameter = betaParameter, nTheta = nTheta, nSim = nSim,
      seed = seed, pb = pb
    )[["power"]] - power
  }

  # No sign change means the target is out of reach on the bracket (or
  # already met at its lower end); report NA rather than a spurious edge.
  targetAtBounds <- c(targetFunction(bounds[1]), targetFunction(bounds[2]))
  if (targetAtBounds[1] >= 0 || targetAtBounds[2] < 0) {
    return(NA_real_)
  }

  stats::uniroot(targetFunction, interval = bounds,
    f.lower = targetAtBounds[1], f.upper = targetAtBounds[2], tol = tol
  )[["root"]]
}
