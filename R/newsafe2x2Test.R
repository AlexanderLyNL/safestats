# Test functions ----

#' Computes Conditional E-Value with UMP Log Odds Ratio
savi2x2TestStatUmp <- function(
  ya, yb, na, nb, alpha,
  alternative = c("twoSided", "greater", "less")
) {
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
    if (is.null(logOdds)) {
      # No one-block test at level alpha can reject on this side, so there
      # is no UMP e-value; the caller keeps block 1 of the plain process.
      return(NULL)
    }
    logLikelihoodAlternative <- logLikelihoodFNCH(ya, yb, na, nb, logOdds)
    eValue <- eValue + exp(logLikelihoodAlternative - logLikelihoodNull)
  }

  # twoSided: 1/2 less + 1/2 greater
  eValue / length(sides)
}

#' Safe anytime-valid 2x2 test for propDiff or logOdds
#'
#' Tests a null of no effect on a sequence of 2x2 tables, one table per
#' block. `designObj[["eType"]]` picks the e-process and, with it, the effect.
#' Both effects are signed A minus B: `propDiff = thetaA - thetaB` and
#' `logOdds = logit(thetaA) - logit(thetaB)`, so `"greater"` means group A
#' has the larger proportion. The numerator for block `i` uses blocks
#' `1..i-1` only, except that eGauss mixes over the grid posterior including
#' block `i`, which is valid because it conditions on the block's total.
#'
#' - `"eBeta"` (propDiff): two independent Beta posterior means for `thetaA`
#'   and `thetaB`, tested against the pooled null mean. `"twoSided"` only.
#' - `"grow"` (propDiff or logOdds, whichever minimal effect the design
#'   holds): the alternative is restricted to the signed minimal effect,
#'   positive with `"greater"`, negative with `"less"`. `"twoSided"` averages
#'   the cumulative e-values at both signs of its magnitude.
#' - `"eGauss"` (logOdds): the conditional FNCH likelihood mixed under the
#'   design's `gaussParameter` Normal prior on a logOdds grid, restricted to
#'   the side of a one-sided `alternative`, against the hypergeometric null.
#'
#' In every case the first table's e-value is replaced by the UMP
#' conditional e-value when one exists at level `alpha`; otherwise block 1
#' keeps its plain factor. Only `"eBeta"` and `"eGauss"` give a confidence
#' interval or sequence, on their own effect.
#'
#' @param ya Successes in groups A, a numeric vector of nonnegative integers.
#' @param yb Successes in groups B, a numeric vector of nonnegative integers.
#' @param designObj an object obtained from \code{\link{designSavi2x2}}.
#' @param wantCi default `TRUE` a confidence interval on the last block
#' @param wantConfidenceSequence `TRUE` for a confidence sequence with one
#'   row per block (`"eBeta"` and `"eGauss"` only). It takes precedence over
#'   `wantCi`.
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
#' - `confSeq`: named `c(lowerBound, upperBound)` for the design's effect
#'   at level `ciValue` on all blocks. With `wantConfidenceSequence`, this is
#'   the last row of `confSeqMatrix`. Both bounds are `NA` when the set is
#'   empty. It is `NULL` for `"grow"`, or when neither `wantCi` nor
#'   `wantConfidenceSequence` is `TRUE`.
#' - `confSeqMatrix`: an `nBlocks x 2` matrix of `lowerBound` and
#'   `upperBound`. Row `i` is the interval on blocks `1..i`. With
#'   `runningIntersection`, the rows are nested and stay `NA` after the first
#'   empty row. It is present only with `wantConfidenceSequence`.
#' - `ciValue`: the confidence level used.
#' - `betaParameter`: `list(betaA1, betaA2, betaB1, betaB2)`, the design's
#'   Beta prior updated with all blocks. This is the prior that a next block
#'   would use.
#' - `alternative`, `h0`, `designObj`: copied from the design.
#' - `dataName`: the deparsed `ya` and `yb` arguments.
#' - `call`: the matched call.
#'
#' The constructor's other fields (`statistic`, `eValueApproxError`,
#' `note`) stay `NULL`.
#' @noRd
savi2x2TestStat <- function(
  ya,
  yb,
  designObj = NULL,
  wantCi = TRUE,
  wantConfidenceSequence = FALSE,
  ciValue = NULL
) {
  esMin <- unname(designObj[["esMin"]])
  alternative <- designObj[["alternative"]]
  betaParameter <- designObj[["betaParameter"]]
  gaussParameter <- designObj[["gaussParameter"]]
  runningIntersection <- designObj[["runningIntersection"]]
  alpha <- designObj[["alpha"]]
  eType <- designObj[["eType"]]
  na <- designObj[["nPlan"]][["na"]]
  nb <- designObj[["nPlan"]][["nb"]]
  # The effect follows from eType, or for grow from which minimal effect
  # the design holds; the design names esMin after it.
  effect <- switch(eType,
    "eBeta" = "propDiff",
    "eGauss" = "logOdds",
    "grow" = names(designObj[["esMin"]])
  )

  nBlocks <- length(ya)
  if (length(na) == 1L) {
    na <- rep(na, nBlocks)
  }
  if (length(nb) == 1L) {
    nb <- rep(nb, nBlocks)
  }

  # Checking: data ----
  # design handle all the other argument
  if (nBlocks == 1L) {
    warnings("There is only 1 table, switched to UMP conditional e-variable")
  }

  if (length(yb) != nBlocks || length(na) != nBlocks || length(nb) != nBlocks) {
    stop("ya, yb, na and nb must have one value per block")
  }
  counts <- c(ya, yb, na, nb)
  if (
    !all(is.finite(counts)) ||
      any(counts %% 1 != 0) ||
      any(c(ya, yb) < 0) ||
      any(c(na, nb) < 1) ||
      any(ya > na) ||
      any(yb > nb)
  ) {
    stop("ya, yb must be integers in 0..na, 0..nb; na, nb positive integers")
  }

  result <- constructSaviTestObj("Two Proportions")

  # Compute: eValueVec ----
  # UMP conditional e-value for block 1
  eValueUmp <- savi2x2TestStatUmp(
    ya[1],
    yb[1],
    na[1],
    nb[1],
    alpha,
    alternative
  )

  # The plain cumulative log e-process of the chosen e-variable, block 1
  # included; each helper is a self-contained construction.
  logEValueVec <- switch(paste(eType, effect),
    "eBeta propDiff" = logEValueVec2x2PropDiffEBeta(
      ya,
      yb,
      na,
      nb,
      betaParameter
    ),
    "grow propDiff" = logEValueVec2x2PropDiffGrow(
      ya,
      yb,
      na,
      nb,
      betaParameter,
      esMin,
      alpha,
      alternative
    ),
    "grow logOdds" = logEValueVec2x2LogOddsGrow(
      ya,
      yb,
      na,
      nb,
      esMin,
      alternative
    ),
    "eGauss logOdds" = logEValueVec2x2LogOddsEGauss(
      ya,
      yb,
      na,
      nb,
      gaussParameter,
      alternative
    )
  )
  # Replace block 1 by the UMP e-value
  # logEValueVec[1] should be log(1) = 0 but write it out for clarity
  # No UMP e-value at this alpha (NULL): block 1 keeps its plain factor.
  if (!is.null(eValueUmp)) {
    logEValueVec <- logEValueVec - logEValueVec[1] + log(eValueUmp)
  } else {
    warning(
      "no UMP e-value exists for the first block at alpha = ",
      alpha,
      "; block 1 keeps the plain e-value"
    )
  }

  # Compute: confSeq ----
  # Only for eBeta and eGauss, use 1 - alpha unless specified
  ciValue <- ifelse(is.null(ciValue), 1 - designObj[["alpha"]], ciValue)
  result[["ciValue"]] <- ciValue

  if (eType == "grow" && (wantCi || wantConfidenceSequence)) {
    warning(
      "Confidence interval/sequences only available for eType = ",
      "'eBeta' or 'eGauss'"
    )
  } else if (wantConfidenceSequence) {
    confSeqMatrix <- if (effect == "propDiff") {
      computeConfidenceSequence2x2PropDiff(
        ya,
        yb,
        na,
        nb,
        betaParameter,
        1 - ciValue,
        runningIntersection
      )
    } else {
      computeConfidenceSequence2x2LogOdds(
        ya,
        yb,
        na,
        nb,
        1 - ciValue,
        runningIntersection,
        gaussParameter,
        alternative
      )
    }
    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSeqMatrix[nBlocks, ]
  } else if (wantCi) {
    # One confidence interval on all blocks
    result[["confSeq"]] <- if (effect == "propDiff") {
      computeConfidenceInterval2x2PropDiff(
        ya,
        yb,
        na,
        nb,
        betaParameter,
        1 - ciValue
      )
    } else {
      # The eGauss numerator on all blocks, without the UMP replacement
      logNumerator <- logEValueVec2x2LogOddsEGauss(
        ya,
        yb,
        na,
        nb,
        gaussParameter,
        alternative
      ) +
        cumsum(stats::dhyper(ya, na, nb, ya + yb, log = TRUE))
      computeConfidenceInterval2x2LogOdds(
        ya,
        yb,
        na,
        nb,
        logNumerator[nBlocks],
        1 - ciValue
      )
    }
  }

  # Fill: Result ----
  # TODO: what happens with overflow
  eValueVec <- exp(logEValueVec)
  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["estimate"]] <- c(
    "thetaA" = sum(ya) / sum(na),
    "thetaB" = sum(yb) / sum(nb)
  )
  # x-axis of plot.saviTest(): the block index.
  result[["n"]] <- c("na" = sum(na), "nb" = sum(nb), "nBlocks" = nBlocks)
  # Beta posterior of the observed blocks, the prior a further block would use.
  # TODO: should I update the prior when eType = "grow"?
  result[["betaParameter"]] <- list(
    "betaA1" = betaParameter[["betaA1"]] + sum(ya),
    "betaA2" = betaParameter[["betaA2"]] + sum(na) - sum(ya),
    "betaB1" = betaParameter[["betaB1"]] + sum(yb),
    "betaB2" = betaParameter[["betaB2"]] + sum(nb) - sum(yb)
  )
  result[["designObj"]] <- designObj
  result[["n1Vec"]] <- seq_len(nBlocks)
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- alternative
  result[["h0"]] <- designObj[["h0"]]
  result[["dataName"]] <- paste(
    deparse1(substitute(ya)),
    "and",
    deparse1(substitute(yb))
  )
  result[["call"]] <- sys.call()

  return(result)
}

# Design functions ----
# TODO: add h0 != 0 situation

#' Design a safe anytime-valid 2x2 test
#'
#' `eType` picks the e-variable and, with it, the effect measure. Both
#' effects are signed A minus B and anchored on `thetaB`: `propDiff =
#' thetaA - thetaB` and `logOdds = logit(thetaA) - logit(thetaB)`, so
#' `"greater"` means group A has the larger proportion, as `x - y > 0`
#' does in [stats::t.test()].
#'
#' - `"eBeta"` (propDiff) and `"eGauss"` (logOdds) are unrestricted; a
#'   `propDiffMin`, `logOddsMin` or `power` is an error for them, while
#'   `nBlocksPlan` is kept as the planned block count. `"eBeta"` is twoSided
#'   only and a one-sided `alternative` is ignored with a warning; `"eGauss"`
#'   restricts its `gaussParameter` prior to the side of a one-sided
#'   `alternative`.
#' - `"grow"` plugs in exactly one of `propDiffMin`, `logOddsMin` as the
#'   fixed alternative.
#'   `"greater"` needs a positive value, `"less"` a negative one, or the
#'   alternative would contain the null and no GROW e-variable exists.
#'   `"twoSided"` uses the magnitude.
#'   Planning exists for `propDiff` only (Decisions 37, 43): `power` alone
#'   plans the block count at the hardest baseline, `nBlocksPlan` alone
#'   evaluates the worst-case power there, and both without `propDiffMin`
#'   find the minimal detectable `propDiff` (Decision 42). Both with
#'   `propDiffMin` errors. Without a minimal effect and without both, grow
#'   warns and continues as `"eBeta"`, dropping a lone `power`.
#'   `logOddsMin` with `power` errors: its worst case is set by the baseline
#'   grid, not by the effect; `nBlocksPlan` is kept as the planned count.
#'
#' @param na number of observations in group a per data block
#' @param nb number of observations in group b per data block
#' @param nBlocksPlan planned number of data blocks collected
#' @param propDiffMin minimal difference in proportions: `thetaA - thetaB`
#' @param logOddsMin minimal log odds ratio: `logOddsA - logOddsB`, that is
#'   `logit(thetaA) - logit(thetaB)`
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error
#' @param power numeric in (0, 1) that specifies the desired power
#' @param h0 0
#' @param alternative `"twoSided"`, `"greater"` or `"less"`, the direction
#'   of the effect A minus B under the alternative.
#' @param eType `"eBeta"`, `"grow"` or `"eGauss"`, see Details.
#' @param betaParameter `list(betaA1, betaA2, betaB1, betaB2)`, the Beta
#'   prior shapes on `thetaA` and `thetaB`; `NULL` keeps the constructor's
#'   default of `0.18` each.
#' @param gaussParameter `list(mean, sd)`, the Normal prior on `logOdds` for
#'   `"eGauss"`, restricted to the grid `(-20, 20)` and to the side of a
#'   one-sided `alternative`; `NULL` means `list(mean = 0, sd = 1)`.
#' @param runningIntersection `TRUE` to intersect each row of the blockwise
#'   confidence sequence with the previous one; `NULL` keeps the
#'   constructor's `FALSE`.
#' @param nSim,nBoot,nMax,seed,wantSamplePaths,pb Simulation settings
#'   passed to [sampleStoppingTimesSavi2x2()] via [computeNPlanSavi2x2()] or
#'   [computePowerSavi2x2()]: paths per baseline,
#'   bootstrap resamples, block cap per path, seed (`NULL` is `2026`),
#'   whether to keep the e-value paths, and the progress bar.
#'
#' @return A `saviDesign` with `testName = "Two Proportions"`, `testType =
#'   "2x2"`, `h0 = c(propDiff = h0)`, `esMin`, `eType`, `alpha`,
#'   `alternative`, `betaParameter`, `gaussParameter`, `parameter` (the prior summarised for
#'   printing), `runningIntersection` and `nPlan = list(na, nb)`, with a
#'   third element `nBlocksPlan` when planned or given. With `power`:
#'   `designScenario = "1a"`, `power` as the target, `nPlanTwoSe = c(NA,
#'   NA, 2 * bootSe)`, `bootObjNBlocksPlan`, `nMean`, `nMeanTwoSe`,
#'   `bootObjNMean`. With `nBlocksPlan`: `designScenario = "2"`, `power`
#'   (the worst case), `powerTwoSe`, `bootObjPower`. Both also carry
#'   `worstCaseThetaA`, `worstCaseThetaB`, `breakVector`
#'   and `samplePaths`. With `power` and `nBlocksPlan` but no `*Min`:
#'   `designScenario = "3"`, `esMin` the minimal detectable `propDiff`,
#'   `power` as the target, and no simulation summaries.
#' @noRd
designSavi2x2 <- function(
  na,
  nb,
  nBlocksPlan = NULL,
  propDiffMin = NULL,
  logOddsMin = NULL,
  alpha = 0.05,
  power = NULL,
  h0 = 0,
  alternative = c("twoSided", "greater", "less"),
  eType = c("eBeta", "grow", "eGauss"),
  betaParameter = NULL,
  gaussParameter = NULL,
  runningIntersection = NULL,
  nSim = 1e3L,
  nBoot = 1e3L,
  nMax = 1e4L,
  seed = NULL,
  wantSamplePaths = FALSE,
  pb = TRUE
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  result <- constructSaviDesignObj("Two Proportions")

  # Fill: result ----
  if (!is.null(runningIntersection)) {
    result[["runningIntersection"]] <- runningIntersection
  }
  result[["eType"]] <- eType
  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["h0"]] <- 0

  # Checking: arg ----
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1)")
  }
  if (
    !is.null(power) &&
      (length(power) != 1L || !is.finite(power) || power <= 0 || power >= 1)
  ) {
    stop("power must be a single number in (0, 1)")
  }

  # na nb are finite positive integer
  # must be same length, if length > 1, then that is the nBlocksPlan
  if (
    !is.numeric(na) ||
      !is.numeric(nb) ||
      length(na) < 1L ||
      length(nb) < 1L ||
      !all(is.finite(c(na, nb))) ||
      any(c(na, nb) %% 1 != 0) ||
      any(c(na, nb) < 1)
  ) {
    stop("na and nb must be finite positive integers")
  }
  if (length(na) != length(nb)) {
    stop(
      "na and nb must have the same length: one value each, or one per block"
    )
  }
  if (
    !is.null(nBlocksPlan) &&
      (length(nBlocksPlan) != 1L ||
        !is.finite(nBlocksPlan) ||
        nBlocksPlan %% 1 != 0 ||
        nBlocksPlan < 1)
  ) {
    stop("nBlocksPlan must be a single positive integer")
  }
  if (length(na) > 1L) {
    if (!is.null(nBlocksPlan) && nBlocksPlan != length(na)) {
      stop(
        "nBlocksPlan = ",
        nBlocksPlan,
        " does not match the ",
        length(na),
        " blocks given by na and nb"
      )
    }
    nBlocksPlan <- length(na)
  }

  # TODO: add parameter for eGauss
  if (!is.null(betaParameter)) {
    result[["betaParameter"]] <- betaParameter
  } else if (length(na) == 1L) {
    result[["betaParameter"]] <- list(
      "betaA1" = 1 / (2 * na),
      "betaA2" = 1 / (2 * na),
      "betaB1" = 1 / (2 * nb),
      "betaB2" = 1 / (2 * nb)
    )
  } else {
    stop("na nb are vectors")
  }
  # Gaussian prior on logOdds for eGauss: N(mean, sd), restricted to the
  # helper's grid (-20, 20) and to the side of a one-sided alternative.
  if (is.null(gaussParameter)) {
    gaussParameter <- list("mean" = 0, "sd" = 1)
  }
  if (
    !is.list(gaussParameter) ||
      !all(c("mean", "sd") %in% names(gaussParameter)) ||
      length(gaussParameter[["mean"]]) != 1L ||
      length(gaussParameter[["sd"]]) != 1L ||
      !is.finite(gaussParameter[["mean"]]) ||
      !is.finite(gaussParameter[["sd"]]) ||
      abs(gaussParameter[["mean"]]) >= 20 ||
      gaussParameter[["sd"]] <= 0
  ) {
    stop(
      "gaussParameter must be list(mean, sd) with a finite mean in (-20, 20) ",
      "and a finite sd > 0"
    )
  }
  result[["gaussParameter"]] <- gaussParameter
  # Simulation settings: positive integers, and the cap nMax at least the
  # planned block count. Checked here so a bad value fails before sampling.
  for (setting in c("nSim", "nBoot", "nMax")) {
    value <- get(setting)
    if (length(value) != 1L || !is.finite(value) || value %% 1 != 0 || value < 1) {
      stop(setting, " must be a single positive integer")
    }
  }
  if (!is.null(nBlocksPlan) && nMax < nBlocksPlan) {
    stop("nMax = ", nMax, " is below the planned block count ", nBlocksPlan)
  }

  # At most one minimal effect
  # propDiffMin in (-1, 1) = thetaA - thetaB
  # logOddsMin finite = log(oddsA / oddsB) = logit(thetaA) - logit(thetaB)
  if (!is.null(propDiffMin) && !is.null(logOddsMin)) {
    stop("supply propDiffMin or logOddsMin, not both")
  }
  if (
    !is.null(propDiffMin) &&
      (length(propDiffMin) != 1L ||
        !is.finite(propDiffMin) ||
        abs(propDiffMin) >= 1)
  ) {
    stop("propDiffMin must be a single number in (-1, 1)")
  }
  if (
    !is.null(logOddsMin) &&
      (length(logOddsMin) != 1L || !is.finite(logOddsMin))
  ) {
    stop("logOddsMin must be a single finite number")
  }
  # Named so print() and plot() show which effect the minimal value is on.
  esMin <- if (!is.null(propDiffMin)) {
    c("propDiff" = propDiffMin)
  } else if (!is.null(logOddsMin)) {
    c("logOdds" = logOddsMin)
  }
  effect <- names(esMin)
  # 0 is no restriction, pls use eBeta or eGauss
  if (!is.null(esMin) && esMin == 0) {
    stop(
      effect,
      "Min = 0 is no restriction; use eType = '",
      if (effect == "propDiff") "eBeta" else "eGauss",
      "'"
    )
  }

  # esMin must < 0 for alternative == "less"
  if (!is.null(esMin) && alternative == "less" && esMin > 0) {
    stop(
      effect,
      "Min = ",
      esMin,
      " with alternative = 'less' puts the null ",
      effect,
      " = 0 inside the alternative ",
      effect,
      " <= ",
      esMin,
      "; for 'less' supply ",
      effect,
      "Min < 0, or use alternative = ",
      "'greater'. A null at ",
      effect,
      " = ",
      esMin,
      " is h0 != 0, which is not designed"
    )
  }
  # esMin must > 0 for alternative == "greater"
  if (!is.null(esMin) && alternative == "greater" && esMin < 0) {
    stop(
      effect,
      "Min = ",
      esMin,
      " with alternative = 'greater' puts the null ",
      effect,
      " = 0 inside the alternative ",
      effect,
      " >= ",
      esMin,
      "; for 'greater' supply ",
      effect,
      "Min > 0, or use alternative = ",
      "'less'. A null at ",
      effect,
      " = ",
      esMin,
      " is h0 != 0, which is not designed"
    )
  }
  # twoSided tests both signs of the value, so require the positive one:
  # the stored esMin then always carries the sign of the direction tested
  if (!is.null(esMin) && alternative == "twoSided" && esMin < 0) {
    stop(
      effect,
      "Min = ",
      esMin,
      " with alternative = 'twoSided'; supply ",
      effect,
      "Min > 0, both signs of it are tested"
    )
  }
  result[["esMin"]] <- esMin

  # Effect and planning must fit eType (design: checks). eBeta and eGauss
  # take no minimal effect and no planning; grow needs a minimal effect, or
  # both power and nBlocksPlan to find one, else it continues as eBeta.
  if (eType != "grow" && !is.null(esMin)) {
    stop(
      effect,
      "Min needs eType = 'grow'; eType = '",
      eType,
      "' has no minimal effect"
    )
  }
  if (eType != "grow" && !is.null(power)) {
    stop(
      "power needs eType = 'grow' with propDiffMin; eType = '",
      eType,
      "' has no planning"
    )
  }
  if (eType == "grow" && is.null(esMin) && (is.null(power) || is.null(nBlocksPlan))) {
    warning(
      "eType = 'grow' needs propDiffMin or logOddsMin, or both power and ",
      "nBlocksPlan to find the minimal propDiff; using eType = 'eBeta'",
      if (!is.null(power)) " and dropping power"
    )
    eType <- "eBeta"
    result[["eType"]] <- eType
    power <- NULL
  }
  if (eType == "grow" && !is.null(propDiffMin) && !is.null(power) && !is.null(nBlocksPlan)) {
    stop("with propDiffMin supply power or nBlocksPlan, not both")
  }
  if (eType == "grow" && !is.null(logOddsMin) && !is.null(power)) {
    stop("no planning on logOdds: power needs propDiffMin")
  }

  # Planning only for propDiffMin
  planning <- NULL
  if (eType == "grow" && !is.null(propDiffMin) && !is.null(power)) {
    # Scenario 1: propDiffMin + power -> worst-case stopping time
    planning <- computeNPlanSavi2x2(
      propDiffMin = propDiffMin,
      na = na,
      nb = nb,
      power = power,
      alpha = alpha,
      alternative = alternative,
      betaParameter = result[["betaParameter"]],
      nSim = nSim,
      nBoot = nBoot,
      nMax = nMax,
      seed = seed,
      wantSamplePaths = wantSamplePaths,
      pb = pb
    )
    result[["designScenario"]] <- "1a"
    result[["power"]] <- power
    nBlocksPlan <- planning[["nPlan"]]
    # na, nb are planned, not simulated: no standard error for them.
    result[["nPlanTwoSe"]] <- c(
      NA,
      NA,
      2 * planning[["bootObjNPlan"]][["bootSe"]]
    )
    result[["bootObjNBlocksPlan"]] <- planning[["bootObjNPlan"]]
    result[["nMean"]] <- c("nMean" = planning[["nMean"]])
    result[["nMeanTwoSe"]] <- 2 * planning[["bootObjNMean"]][["bootSe"]]
    result[["bootObjNMean"]] <- planning[["bootObjNMean"]]
  } else if (eType == "grow" && !is.null(propDiffMin) && !is.null(nBlocksPlan)) {
    # Scenario 2: propDiffMin + nBlocksPlan -> worst-case power
    planning <- computePowerSavi2x2(
      propDiffMin = propDiffMin,
      na = na,
      nb = nb,
      nBlocks = nBlocksPlan,
      alpha = alpha,
      alternative = alternative,
      betaParameter = result[["betaParameter"]],
      nSim = nSim,
      nBoot = nBoot,
      seed = seed,
      wantSamplePaths = wantSamplePaths,
      pb = pb
    )
    result[["designScenario"]] <- "2"
    result[["power"]] <- planning[["power"]]
    result[["powerTwoSe"]] <- 2 * planning[["bootObjPower"]][["bootSe"]]
    result[["bootObjPower"]] <- planning[["bootObjPower"]]
  } else if (eType == "grow" && is.null(esMin)) {
    # Scenario 3: power + nBlocksPlan -> the minimal detectable propDiff,
    esMin <- computeEsMinSavi2x2(
      na = na,
      nb = nb,
      nBlocksPlan = nBlocksPlan,
      power = power,
      alpha = alpha,
      alternative = alternative,
      betaParameter = result[["betaParameter"]],
      nSim = nSim,
      seed = seed,
      pb = pb
    )
    # NA: the worst-case power never reaches the target on the search
    # bounds, so no minimal effect can be reported.
    if (is.na(esMin)) {
      stop(sprintf(
        paste(
          "no minimal propDiff found: at nBlocksPlan = %g the worst-case power",
          "does not reach %g for any magnitude in (0.01, 0.9); try a larger",
          "nBlocksPlan or a smaller power"
        ),
        nBlocksPlan,
        power
      ))
    }
    result[["designScenario"]] <- "3"
    result[["esMin"]] <- c("propDiff" = esMin)
    result[["power"]] <- power
  }
  # Scenarios 1 and 2 keep the worst baseline and its simulated paths.
  if (!is.null(planning)) {
    result[["worstCaseThetaA"]] <- planning[["worstCaseThetaA"]]
    result[["worstCaseThetaB"]] <- planning[["worstCaseThetaB"]]
    result[["breakVector"]] <- planning[["breakVector"]]
    result[["samplePaths"]] <- planning[["samplePaths"]]
  }

  result[["parameter"]] <- if (eType == "eGauss") {
    c("Gaussian prior (mean, sd)" = paste(unlist(gaussParameter), collapse = " "))
  } else {
    c(
      "Beta hyperparameters" = paste(
        unlist(result[["betaParameter"]]),
        collapse = " "
      )
    )
  }

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

#' Computes the savi confidence interval for propDiff
#' @noRd
computeConfidenceInterval2x2PropDiff <- function(
  ya,
  yb,
  na,
  nb,
  betaParameter,
  alpha,
  domain = c(-1, 1)
) {
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]

  # product of the seq-RIPr e-variable against the H0: thetaA = thetaB + propDiff
  # f: find zero points against 1/alpha
  # this is a convex function
  fPropDiff <- function(propDiff) {
    nullThetaA <- mapply(
      solveRIPr2x2PropDiff,
      thetaA = thetaA,
      thetaB = thetaB,
      na = na,
      nb = nb,
      MoreArgs = list(propDiff = propDiff)
    )
    sum(
      stats::dbinom(ya, na, thetaA, log = TRUE) +
        stats::dbinom(yb, nb, thetaB, log = TRUE) -
        stats::dbinom(ya, na, nullThetaA, log = TRUE) -
        stats::dbinom(yb, nb, nullThetaA - propDiff, log = TRUE)
    ) -
      log(1 / alpha)
  }

  # The projection needs thetaA strictly inside its range, so stay just
  # inside (-1, 1).
  eps <- 1e-9
  lower <- max(domain[1], -1 + eps)
  upper <- min(domain[2], 1 - eps)
  minimiser <- stats::optimize(fPropDiff, interval = c(lower, upper))[["minimum"]]

  # min > 1 / alpha, no confidence interval found
  if (fPropDiff(minimiser) >= 0) {
    warning("No confidence interval is found! return a non-informative CI: (-1,1)")
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

#' Computes the savi confidence sequences for propDiff using a fixed grid
#' @noRd
computeConfidenceSequence2x2PropDiff <- function(
  ya,
  yb,
  na,
  nb,
  betaParameter,
  alpha,
  runningIntersection,
  nGrid = 2000L
) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]
  logThreshold <- log(1 / alpha)

  # 2000 candidates strictly inside (-1, 1). The final interval has width
  # on the scale of sdMax; warn when the grid step is coarser than that.
  nGrid <- nGrid
  sdMax <- sqrt(1 / (4 * sum(na)) + 1 / (4 * sum(nb)))
  if (ceiling(2 / sdMax) > nGrid) {
    warning(
      "The confidence sequence grid step ",
      2 / nGrid,
      " is coarser than the standard deviation scale ",
      signif(sdMax, 3),
      " of propDiff on these totals; the bounds are conservative"
    )
  }
  grid <- seq(-1, 1, length.out = nGrid + 2L)[-c(1L, nGrid + 2L)]

  # Cumulative log e-process per candidate; which candidates may still be
  # kept, and which take part in the update (the same, plus two guard
  # nodes on each side for the secant under the running intersection).
  logEValues <- numeric(nGrid)
  candidate <- rep(TRUE, nGrid)
  active <- rep(TRUE, nGrid)
  previous <- c(-1, 1)

  confSeqMatrix <- matrix(
    NA_real_,
    nBlocks,
    2,
    dimnames = list(NULL, c("lowerBound", "upperBound"))
  )

  for (i in seq_len(nBlocks)) {
    delta <- grid[active]

    # find RIPr thetaA for all active delta
    nullThetaA <- mapply(
      solveRIPr2x2PropDiff,
      delta,
      MoreArgs = list(
        thetaA = thetaA[i],
        thetaB = thetaB[i],
        na = na[i],
        nb = nb[i]
      )
    )

    # Block i's log likelihood ratio term against each active candidate
    logEValues[active] <- logEValues[active] +
      stats::dbinom(ya[i], na[i], thetaA[i], log = TRUE) +
      stats::dbinom(yb[i], nb[i], thetaB[i], log = TRUE) -
      stats::dbinom(ya[i], na[i], nullThetaA, log = TRUE) -
      stats::dbinom(yb[i], nb[i], nullThetaA - delta, log = TRUE)

    # Kept candidates; one run by convexity
    f <- logEValues - logThreshold
    kept <- candidate & f < 0
    if (!any(kept)) {
      if (runningIntersection) {
        break
      }
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
      if (lowerBound > upperBound) {
        break
      }
      previous <- c(lowerBound, upperBound)
    }
    confSeqMatrix[i, ] <- c(lowerBound, upperBound)
  }

  return(confSeqMatrix)
}

#' Anytime-valid confidence interval for logOdds on all blocks
#'
#' Find the confidence interval by conditional likelihood ratio between the
#' alternative and the null (logOdds in domain)
#' The confidence interval is the set of logOdds values for which the
#' conditional likelihood ratio is below 1 / alpha.
#'
#' @param logPTotal The numerator's cumulative log likelihood after the
#'   last block: the eGauss log numerator, the log e-value plus the
#'   cumulative hypergeometric log likelihood.
#' @param domain `c(lower, upper)`, the candidates searched; a previous
#'   interval gives the running intersection.
#' @return Named numeric `c(lowerBound, upperBound)`; the `domain` edge
#'   when that edge is still inside, and the whole `domain` with a warning
#'   when the set is empty.
#' @noRd
computeConfidenceInterval2x2LogOdds <- function(
  ya,
  yb,
  na,
  nb,
  logPTotal,
  alpha,
  domain = c(-40, 40)
) {
  # product of the conditional e-variable against H0: logOdds = delta
  # f: find zero points against 1 / alpha
  # f is convex in logOdds
  fLogOdds <- function(logOdds) {
    logPTotal -
      sum(logLikelihoodFNCH(ya, yb, na, nb, logOdds)) -
      log(1 / alpha)
  }

  minimiser <- stats::optimize(fLogOdds, interval = domain)[["minimum"]]

  # min > 1 / alpha, no confidence interval found
  if (fLogOdds(minimiser) >= 0) {
    warning(
      "No confidence interval is found on (",
      round(domain[1], 4),
      ",",
      round(domain[2], 4),
      ")"
    )
    return(c("lowerBound" = domain[1], "upperBound" = domain[2]))
  }

  # Still inside at the search edge: the bound is the domain edge.
  lowerBound <- if (fLogOdds(domain[1]) < 0) {
    warning("Cannot find lowerBound for logOdds, return NA")
    NA_real_
  } else {
    stats::uniroot(fLogOdds, lower = domain[1], upper = minimiser)[["root"]]
  }
  upperBound <- if (fLogOdds(domain[2]) < 0) {
    warning("Cannot find upperBound for logOdds, return NA")
    NA_real_
  } else {
    stats::uniroot(fLogOdds, lower = minimiser, upper = domain[2])[["root"]]
  }

  return(c("lowerBound" = lowerBound, "upperBound" = upperBound))
}

#' Compute the eGauss confidence sequence for logOdds
#'
#' Row i inverts the e-process on blocks 1..i: every numerator factor is
#' fixed given its block's total or predictable, so the eGauss numerator at
#' block i is the numerator of the test on blocks 1..i just as its last
#' element is at the end (Decision 30). Recomputed from scratch per block,
#' quadratic in nBlocks.
#'
#' @param runningIntersection `TRUE` searches each row inside the previous
#'   one; an empty row stays empty thereafter.
#' @param gaussParameter,alternative The design's eGauss prior and side, passed
#'   to [logEValueVec2x2LogOddsEGauss()].
#' @return An `nBlocks x 2` matrix of `lowerBound` and `upperBound`; an
#'   empty set is an `NA` row.
#' @noRd
computeConfidenceSequence2x2LogOdds <- function(
  ya,
  yb,
  na,
  nb,
  alpha,
  runningIntersection,
  gaussParameter = NULL,
  alternative = c("twoSided", "greater", "less")
) {
  # THIS IS SLOW, WE MIGHT WANT TO STOP EARLY!
  alternative <- match.arg(alternative)
  nBlocks <- length(ya)
  # The eGauss numerator on blocks 1..i, without the UMP replacement
  logNumerator <- logEValueVec2x2LogOddsEGauss(
    ya,
    yb,
    na,
    nb,
    gaussParameter,
    alternative
  ) +
    cumsum(stats::dhyper(ya, na, nb, ya + yb, log = TRUE))

  confSeqMatrix <- matrix(
    NA_real_,
    nBlocks,
    2,
    dimnames = list(NULL, c("lowerBound", "upperBound"))
  )
  domain <- c(-40, 40) # TODO: maybe add bounds print or else
  for (i in seq_len(nBlocks)) {
    # The empty set is reported with a warning; here it is an NA row.
    row <- tryCatch(
      computeConfidenceInterval2x2LogOdds(
        ya[1:i],
        yb[1:i],
        na[1:i],
        nb[1:i],
        logNumerator[i],
        alpha,
        domain = domain
      ),
      warning = function(w) c("lowerBound" = NA_real_, "upperBound" = NA_real_)
    )
    if (runningIntersection) {
      # if NA is found in previous confSeq, break
      if (anyNA(row)) {
        break
      }
      # update the domain for logOdds for next table
      domain <- row
    }
    confSeqMatrix[i, ] <- row
  }

  confSeqMatrix
}

# Helpers: propDiff ----

# propDiff = thetaA - thetaB always

# Learn the thetaA and thetaB with two independent Beta
# TODO: is there a betaA1 * na term? check with Peter
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

# a vector of logEValue that is the cumsum of likelihood ratio process
logEValueVec2x2PropDiffEBeta <- function(ya, yb, na, nb, betaParameter) {
  thetas <- predictiveThetas2x2(ya, yb, na, nb, betaParameter)
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]

  # null theta is the pooled mean, of length nBlocks
  thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

  # Cumulative log likelihood ratio
  # we use dbinom instead na * log(theta) + (na - ya) * log(1 - theta)
  # to avoid manually handle NaNs
  # binom terms cancels out in logLikelihoodNull and logLikelihoodAlternative
  # TODO: future might switch to y * log(theta) + (n - y) * log1p(-theta)
  # once theta is guarded
  logLikelihoodNull <- cumsum(
    stats::dbinom(ya, na, thetaNull, log = TRUE) +
      stats::dbinom(yb, nb, thetaNull, log = TRUE)
  )
  logLikelihoodAlternative <- cumsum(
    stats::dbinom(ya, na, thetaA, log = TRUE) +
      stats::dbinom(yb, nb, thetaB, log = TRUE)
  )

  logLikelihoodAlternative - logLikelihoodNull
}

# Cumulative log e-process of the grow propDiff test. The alternative is
# restricted to the curve thetaA = thetaB + propDiffMin (sign from
# alternative), which leaves one free parameter; a grid posterior on it gives
# the predictable plug-in for each block. Element i covers blocks 1..i;
# earlyStopping cuts the vector at the first crossing of 1 / alpha.
logEValueVec2x2PropDiffGrow <- function(ya, yb, na, nb, betaParameter,
  propDiffMin, alpha, alternative = c("twoSided", "greater", "less"),
  earlyStopping = FALSE,
  nWeight = 1e3L
) {
  alternative <- match.arg(alternative)
  nBlocks <- length(ya)
  logThreshold <- log(1 / alpha)

  # Only the magnitude is used
  propDiffMin <- abs(propDiffMin)

  # +1 (greater): thetaA = thetaB + propDiffMin
  # -1 (less):    thetaA = thetaB - propDiffMin
  # twoSided runs both and averages them
  signs <- switch(alternative,
    "twoSided" = c(1, -1),
    "greater" = 1,
    "less" = -1
  )
  nSides <- length(signs)

  # smaller theta in (0, 1 - propDiffMin)
  # larger theta in (propDiffMin, 1)
  # free parameter rho in (0, 1) is the smaller one rescaled:
  #   thetaSmall = rho * (1 - propDiffMin),  thetaLarge = thetaSmall + propDiffMin
  rho <- seq(1 / nWeight, 1 - 1 / nWeight, length.out = nWeight)
  thetaSmall <- rho * (1 - propDiffMin)
  thetaLarge <- thetaSmall + propDiffMin

  # nWeight x nSides grids, one column per side.
  thetaAGrid <- vapply(
    signs,
    function(sign) if (sign > 0) thetaLarge else thetaSmall,
    numeric(nWeight)
  )
  thetaBGrid <- vapply(
    signs,
    function(sign) if (sign > 0) thetaSmall else thetaLarge,
    numeric(nWeight)
  )
  logThetaA <- log(thetaAGrid)
  logOneMinusThetaA <- log1p(-thetaAGrid)
  logThetaB <- log(thetaBGrid)
  logOneMinusThetaB <- log1p(-thetaBGrid)

  # Prior: Beta(1, 1) on rho, uniform prior on
  logWeights <- stats::dbeta(rho, 1, 1,log = TRUE)
  logWeights <- matrix(logWeights - max(logWeights), nWeight, nSides)

  # Cumulative log likelihoods per side. Block i plugs in the posterior mean
  # given blocks 1..i-1 (block 1 uses the prior) and joins the posterior
  # afterwards. Row i holds the cumulative log e-value of blocks 1..i.
  logLikelihoodNull <- numeric(nSides)
  logLikelihoodAlternative <- numeric(nSides)
  logEValueSides <- matrix(NA_real_, nBlocks, nSides)

  for (i in seq_len(nBlocks)) {
    for (s in seq_len(nSides)) {
      # calculate thetaA mean
      thetaA <- stats::weighted.mean(thetaAGrid[, s], exp(logWeights[, s]))
      # thetaB follows by restriction
      thetaB <- thetaA - signs[s] * propDiffMin

      # GRO thetaA = thetaB
      thetaNull <- (na[i] * thetaA + nb[i] * thetaB) / (na[i] + nb[i])

      logLikelihoodNull[s] <- logLikelihoodNull[s] +
        stats::dbinom(ya[i], na[i], thetaNull, log = TRUE) +
        stats::dbinom(yb[i], nb[i], thetaNull, log = TRUE)
      logLikelihoodAlternative[s] <- logLikelihoodAlternative[s] +
        stats::dbinom(ya[i], na[i], thetaA, log = TRUE) +
        stats::dbinom(yb[i], nb[i], thetaB, log = TRUE)
    }

    logEValueSides[i, ] <- logLikelihoodAlternative - logLikelihoodNull

    # Early stopping tests the side-averaged cumulative e-value
    if (earlyStopping) {
      logEValueMax <- max(logEValueSides[i, ])
      logEValue <- logEValueMax +
        log(mean(exp(logEValueSides[i, ] - logEValueMax)))
      if (logEValue >= logThreshold) {
        logEValueSides <- logEValueSides[seq_len(i), , drop = FALSE]
        break
      }
    }

    # Only now add block i to each posterior, so block i + 1 is predicted
    # from the past alone. Binomial coefficients are constant over the grid
    # and are left out.
    logWeights <- logWeights +
      ya[i] * logThetaA +
      (na[i] - ya[i]) * logOneMinusThetaA +
      yb[i] * logThetaB +
      (nb[i] - yb[i]) * logOneMinusThetaB
    # Re-centre each column at its maximum (a per-column constant leaves the
    # posterior mean unchanged) so exp() stays safe as blocks accumulate.
    logWeights <- sweep(logWeights, 2, apply(logWeights, 2, max))
  }

  # twoSided: average the two cumulative e-values
  # One side: the value itself
  logEValueMax <- apply(logEValueSides, 1, max)
  logEValueMax +
    log(rowMeans(exp(logEValueSides - logEValueMax)))
}

# Find the means that minimize the KL between the alternative and null
# The null is H0: thetaA - thetaB = propDiff
# TODO: this can be a cubic function solver
solveRIPr2x2PropDiff <- function(
  thetaA,
  thetaB,
  na,
  nb,
  propDiff,
  tol = 1e-12
) {
  derivativeKL <- function(nullThetaA) {
    nullThetaB <- nullThetaA - propDiff
    na * ((1 - thetaA) / (1 - nullThetaA) - thetaA / nullThetaA) +
      nb * ((1 - thetaB) / (1 - nullThetaB) - thetaB / nullThetaB)
  }

  # nullThetaA must keep both nullThetaA and nullThetaA - propDiff in
  # (0, 1); the derivative is infinite at the edges, so search just inside.
  stats::uniroot(
    derivativeKL,
    lower = max(0, propDiff) + tol,
    upper = min(1, 1 + propDiff) - tol,
    tol = tol
  )[["root"]]
}

# Worst-case baseline of the propDiff grow test: the pair on the curve of
# each alternative at which the e-process grows slowest, so a block count
# planned there holds for every baseline. The expected log e-increment per
# block is the KL divergence of the truth from its pooled null projection,
#   R(theta) = nLow KL(theta || theta0) + nHigh KL(theta + d || theta0),
#   theta0   = theta + nHigh d / (nLow + nHigh),
# with theta the lower proportion, nLow its group size and nHigh the size
# of the group at theta + d. Its minimiser is the root of
#   nLow logit(theta) + nHigh logit(theta + d) = (nLow + nHigh) logit(theta0),
# summed over blocks when the sizes vary per block. Equal sizes give exactly
# the midpoint (1 - d) / 2; in general the root is
# (1 - d) / 2 + d (nHigh - nLow) / (6 (nLow + nHigh)) + O(d^3). The rate
# ignores the first-block UMP factor, which moves the simulated worst case
# when one group is very small (2x2-review.md, R9).
#
# "greater": thetaA = thetaB + d, so thetaB is the lower proportion and
# nLow = nb. "less": thetaA = thetaB - d, so thetaA is the lower one and
# nLow = na. "twoSided" gives both, greater first. One row per curve with
# columns thetaA, thetaB.
solveWorstCaseTheta2x2PropDiff <- function(
  propDiffMin,
  na,
  nb,
  alternative = c("twoSided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  d <- abs(propDiffMin)
  stopifnot(length(d) == 1, d > 0, d < 1, all(na >= 1), all(nb >= 1))

  lowerProportion <- function(nLow, nHigh) {
    if (all(nLow == nHigh)) {
      return((1 - d) / 2)
    }
    stationarity <- function(theta) {
      theta0 <- theta + nHigh * d / (nLow + nHigh)
      sum(
        nLow * stats::qlogis(theta) +
          nHigh * stats::qlogis(theta + d) -
          (nLow + nHigh) * stats::qlogis(theta0)
      )
    }
    # logit is infinite at the edges, so search just inside them.
    stats::uniroot(stationarity, c(1e-9, 1 - d - 1e-9), tol = 1e-10)[["root"]]
  }

  greater <- NULL
  less <- NULL
  if (alternative != "less") {
    thetaB <- lowerProportion(nb, na)
    greater <- data.frame("thetaA" = thetaB + d, "thetaB" = thetaB)
  }
  if (alternative != "greater") {
    thetaA <- lowerProportion(na, nb)
    less <- data.frame("thetaA" = thetaA, "thetaB" = thetaA + d)
  }
  rbind(greater, less)
}

# Helpers: logOdds ----

# logOdds is always log(oddsA / oddsB) = logit(thetaA) - logit(thetaB)
# it tells the difference between odds in group A against group B.
# conditioning on ya + yb, ya is Fisher's noncentral hypergeometric with odds
# exp(logOdds) on group A.

# A vector of conditional log likelihood at one logOdds
# TODO: what if ya, yb, na, nb are vector
logLikelihoodFNCH <- function(ya, yb, na, nb, logOdds) {
  mapply(
    function(ya, yb, na, nb) {
      log(BiasedUrn::dFNCHypergeo(ya, na, nb, ya + yb, exp(logOdds)))
    },
    ya = ya,
    yb = yb,
    na = na,
    nb = nb
  )
}

# a vector of logEValue that is the cumsum of likelihood ratio process
logEValueVec2x2LogOddsGrow <- function(
  ya,
  yb,
  na,
  nb,
  logOddsMin,
  alternative = c("twoSided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  nBlocks <- length(ya)

  logOddsMin <- abs(logOddsMin)
  signs <- switch(alternative,
    "twoSided" = c(1, -1),
    "greater" = 1,
    "less" = -1
  )

  logLikelihoodNull <- cumsum(stats::dhyper(ya, na, nb, ya + yb, log = TRUE))

  # nBlocks x nSides: cumulative conditional log likelihood at sign * logOddsMin
  logLikelihoodSides <- matrix(
    vapply(
      signs,
      function(sign) cumsum(logLikelihoodFNCH(ya, yb, na, nb, sign * logOddsMin)),
      numeric(nBlocks)
    ),
    nrow = nBlocks
  )
  logEValueSides <- logLikelihoodSides - logLikelihoodNull

  # row maximum over the sides, for the log-mean-exp below
  logEValueMax <- do.call(
    pmax,
    lapply(seq_len(ncol(logEValueSides)), function(j) logEValueSides[, j])
  )
  logEValueMax + log(rowMeans(exp(logEValueSides - logEValueMax)))
}

# a vector of logEValue that is the cumsum of conditional likelihood ratio process
logEValueVec2x2LogOddsEGauss <- function(
  ya,
  yb,
  na,
  nb,
  gaussParameter = NULL,
  alternative = c("twoSided", "greater", "less"),
  logOddsGrid = seq(-20, 20, length.out = 2000)
) {
  alternative <- match.arg(alternative)
  if (is.null(gaussParameter)) {
    gaussParameter <- list("mean" = 0, "sd" = 1)
  }
  nBlocks <- length(ya)
  nGrid <- length(logOddsGrid)

  # logPrior: length nGrid. N(mean, sd) on the grid, restricted to the side of
  # the alternative (greater: logOdds > 0, less: logOdds < 0, twoSided: all),
  # normalised on the log scale; an excluded point has weight exp(-Inf) = 0.
  logPrior <- stats::dnorm(
    logOddsGrid,
    gaussParameter[["mean"]],
    gaussParameter[["sd"]],
    log = TRUE
  )
  onSide <- switch(alternative,
    "twoSided" = rep(TRUE, nGrid),
    "greater" = logOddsGrid > 0,
    "less" = logOddsGrid < 0
  )
  logPrior[!onSide] <- -Inf
  logPrior <- logPrior - max(logPrior) - log(sum(exp(logPrior - max(logPrior))))

  # nBlocks x grid: FNCH log density of ya at every grid logOdds, the
  # odds exp(logOdds) on group A; k runs over the feasible ya.
  # compute all the likelihood for all nBlocks and all
  # logPGrid: nBlocks x nGrid, block i's conditional log likelihood at grid point k
  logPGrid <- t(mapply(
    function(ya, yb, na, nb) {
      k <- max(0, ya + yb - nb):min(na, ya + yb)
      logTerms <- outer(logOddsGrid, k) +
        rep(
          lchoose(na, k) + lchoose(nb, ya + yb - k),
          each = length(logOddsGrid)
        )
      shift <- apply(logTerms, 1, max)
      logTerms[, ya - k[1] + 1] - shift - log(rowSums(exp(logTerms - shift)))
    },
    ya = ya,
    yb = yb,
    na = na,
    nb = nb
  ))

  # Cumulate over blocks (cumsum down each grid column; matrix() keeps
  # a single block as a 1-row matrix), then add the log prior weight to
  # every row of its column: sweep(x, 2, v, "+") adds v[j] to column j.
  # Row i then holds log(prior * likelihood of blocks 1..i) on the grid,
  # mixed over the grid by the log-sum-exp below.
  # logMix: nBlocks x nGrid, log(prior_k * likelihood of blocks 1..i at k)
  logMix <- sweep(
    matrix(apply(logPGrid, 2, cumsum), nrow = nBlocks),
    2,
    logPrior,
    "+"
  )
  # shift: length nBlocks, the row maxima; logNumerator: length nBlocks, the
  # log of each row's sum over the grid
  shift <- apply(logMix, 1, max)
  logNumerator <- shift + log(rowSums(exp(logMix - shift)))

  logLikelihoodNull <- cumsum(stats::dhyper(ya, na, nb, ya + yb, log = TRUE))

  logNumerator - logLikelihoodNull
}

#' Log partition function of Fisher's noncentral hypergeometric distribution
#'
#' @param na nonnegative integer, group size of A
#' @param nb nonnegative integer, group size of B
#' @param totalSuccesses nonnegative integer, the total successes ya + yb in the table, at most na + nb
#'
#' @return a numeric
#' @noRd
fnchLogPartition <- function(na, nb, totalSuccesses, logOdds) {
  if (logOdds == 0) {
    return(lchoose(na + nb, totalSuccesses))
  }

  # k: number of success is group A
  feasibleSuccesses <- max(0, totalSuccesses - nb):min(na, totalSuccesses)

  # sum: choose(na, k) * choose(nb, totalSuccesses - k) * exp(logOdds * k)
  logTerms <- lchoose(na, feasibleSuccesses) +
    lchoose(nb, totalSuccesses - feasibleSuccesses) +
    logOdds * feasibleSuccesses

  # log-sum-exp, shifted by the largest term against overflow
  maxLogTerm <- max(logTerms)
  maxLogTerm + log(sum(exp(logTerms - maxLogTerm)))
}

#' Compute the Log Odds Ratio that is Uniformly Most Powerful for One Sided Test
solveUmpLogOdds <- function(
  na,
  nb,
  totalSuccesses,
  alpha,
  alternative = c("greater", "less"),
  logOddsNull = 0,
  searchBound = 100
) {
  alternative <- match.arg(alternative)

  # f: KL(logOdds || logOddsNull) - log(1/alpha) is uniroot on each side
  # logOdds is logOddsA - logOddsB
  klMinusTarget <- function(logOdds) {
    (logOdds - logOddsNull) *
      BiasedUrn::meanFNCHypergeo(na, nb, totalSuccesses, exp(logOdds)) -
      fnchLogPartition(na, nb, totalSuccesses, logOdds) +
      fnchLogPartition(na, nb, totalSuccesses, logOddsNull) +
      log(alpha)
  }

  # greater, root should be on + side, logOddsA > logOddsB
  bounds <- if (alternative == "greater") {
    c(logOddsNull, logOddsNull + searchBound)
  } else {
  # less, root should be on - side, logOddsA < logOddsB
    c(logOddsNull - searchBound, logOddsNull)
  }

  # uniroot() would error on equal signs
  if (klMinusTarget(bounds[1]) * klMinusTarget(bounds[2]) > 0) {
    msg <- sprintf(
      "No root for UMP logOdds at alpha = %s on bounds (%.3f, %.3f).
      Try decreasing alpha or increasing searchBound!",
      alpha, bounds[1], bounds[2]
    )
    warning(msg, call. = FALSE)
    return(NULL)
  }

  stats::uniroot(
    klMinusTarget,
    lower = bounds[1],
    upper = bounds[2],
    tol = 1e-10
  )[["root"]]
}

# Sampling functions for design ----

#' Simulate stopping times of the 2x2 grow test
#'
#' Decisions 29, 34, 43. `propDiffMin` is the grow plug-in and the
#' data-generating effect; only its magnitude is used, the curve(s) follow
#' `alternative`: data lie on `thetaA = thetaB + |propDiffMin|` for
#' `"greater"`, on `thetaA = thetaB - |propDiffMin|` for `"less"`, and on
#' both curves for `"twoSided"`, at the worst-case baseline of each curve
#' from [solveWorstCaseTheta2x2PropDiff()]. There is no planning on
#' `logOdds` (Decision 43).
#' `nPlan` is the worst `power` quantile of the stopping time over the
#' curves; with `power = NULL` the quantile step is skipped and `nPlan`,
#' `worstCaseThetaA`, `worstCaseThetaB` are `NULL` (Decision 36).
#'
#' @return A list: `thetaA`, `thetaB`, and one row per baseline in
#'   `stoppingTimes` (`Inf` when a path never crosses `1 / alpha`),
#'   `breakVector` (`0` crossed, `1` reached `nMax`), `eValuesStopped`;
#'   `samplePaths` (a list of `nSim x nMax` sparse matrices, or `NULL`),
#'   `n1Vector` (the block index), `nPlan`, and the baseline it was taken
#'   at, `worstCaseThetaA`, `worstCaseThetaB`.
#' @noRd
sampleStoppingTimesSavi2x2 <- function(
  propDiffMin,
  na,
  nb,
  power = NULL,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  eType = c("grow"),
  betaParameter = NULL,
  nSim = 1e3L,
  nMax = 1e4L,
  nBoot = 1e4L,
  seed = NULL,
  wantEValuesAtNMax = FALSE,
  wantSamplePaths = FALSE,
  wantSimData = TRUE,
  pb = TRUE
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  propDiffMin <- abs(propDiffMin)
  stopifnot(
    length(propDiffMin) == 1,
    propDiffMin > 0,
    propDiffMin < 1,
    alpha > 0,
    alpha < 1,
    is.null(power) || (power > 0 && power < 1),
    na >= 1,
    nb >= 1,
    is.finite(nMax)
  )
  # One size per group is repeated over the nMax simulated blocks; one size
  # per block needs nMax equal to the block count.
  if (length(na) > 1L || length(nb) > 1L) {
    stopifnot(length(na) == length(nb), nMax == length(na))
  }

  if (is.null(betaParameter)) {
    betaParameter <- constructSaviDesignObj("Two Proportions")[[
      "betaParameter"
    ]]
  }

  set.seed(if (is.null(seed)) 2026 else seed)

  # given propDiffMin, worstCaseTheta make stopping time largest
  worstCaseTheta <- solveWorstCaseTheta2x2PropDiff(propDiffMin, na, nb, alternative)
  thetaATrue <- worstCaseTheta[["thetaA"]]
  thetaBTrue <- worstCaseTheta[["thetaB"]]
  nBaselines <- length(thetaATrue)

  logThreshold <- log(1 / alpha)
  naVec <- rep_len(na, nMax)
  nbVec <- rep_len(nb, nMax)
  stoppingTimes <- matrix(Inf, nrow = nBaselines, ncol = nSim)
  breakVector <- matrix(1L, nrow = nBaselines, ncol = nSim)
  eValuesStopped <- matrix(NA_real_, nrow = nBaselines, ncol = nSim)
  samplePaths <- if (wantSamplePaths) vector("list", nBaselines) else NULL

  message("Simulating ", nSim, " stopping times for propDiffMin = ", propDiffMin)

  if (pb) {
    pbSavi <- utils::txtProgressBar(
      style = 1,
      title = "Sampling worst-case stopping time"
    )
  }

  # 1 for oneSided, 2 for twoSided
  for (k in seq_len(nBaselines)) {
    if (wantSamplePaths) {
      samplePaths[[k]] <- matrix(0, nrow = nSim, ncol = nMax)
    }
    for (sim in seq_len(nSim)) {
      if (pb) {
        utils::setTxtProgressBar(
          pbSavi,
          "value" = ((k - 1) * nSim + sim) / (nBaselines * nSim)
        )
      }

      # stop at the first time crossing 1/alpha
      ya <- stats::rbinom(nMax, naVec, thetaATrue[k])
      yb <- stats::rbinom(nMax, nbVec, thetaBTrue[k])
      logEValueVec <- logEValueVec2x2PropDiffGrow(
        ya,
        yb,
        naVec,
        nbVec,
        betaParameter,
        propDiffMin,
        alpha,
        alternative,
        earlyStopping = TRUE
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

  if (pb) {
    close(pbSavi)
  }

  # Finding percentage of stopping times according to power
  nPlan <- NULL
  worstCaseThetaA <- NULL
  worstCaseThetaB <- NULL
  if (!is.null(power)) {
    quantiles <- apply(
      stoppingTimes,
      1,
      stats::quantile,
      probs = power,
      names = FALSE,
      type = 1
    )

    # For twoSided, find the worst of "less" and "greater"
    worstSign <- which.max(quantiles)
    nPlan <- ceiling(quantiles[worstSign])
    worstCaseThetaA <- thetaATrue[worstSign]
    worstCaseThetaB <- thetaBTrue[worstSign]
  }

  if (!is.null(nPlan) && !is.finite(nPlan)) {
    fractionNeverCrossed <- mean(!is.finite(stoppingTimes[worstSign, ]))
    warning(sprintf(
      paste(
        "the %g quantile of the stopping time is Inf: %.1f%% of the paths at",
        "thetaA = %.3f never cross 1/alpha at nMax = %g, try increasing nMax",
        "or propDiffMin"
      ),
      power,
      100 * fractionNeverCrossed,
      worstCaseThetaA,
      nMax
    ))
  }

  list(
    "thetaA" = thetaATrue,
    "thetaB" = thetaBTrue,
    "stoppingTimes" = stoppingTimes,
    "breakVector" = breakVector,
    "eValuesStopped" = eValuesStopped,
    "samplePaths" = samplePaths,
    "n1Vector" = seq_len(nMax),
    "nPlan" = nPlan,
    "worstCaseThetaA" = worstCaseThetaA,
    "worstCaseThetaB" = worstCaseThetaB
  )
}


#' Worst-case power of the 2x2 grow test at a planned block count
#'
#' Runs [sampleStoppingTimesSavi2x2()] with `nMax = nBlocks`
#' and reports, per baseline, the fraction of paths that cross `1 / alpha`
#' within `nBlocks`; the worst case is the smallest of these, since the
#' test must meet its target whatever `thetaA` is.
#'
#' @param nBlocks Planned block count at which the test is evaluated.
#' @inheritParams sampleStoppingTimesSavi2x2
#'
#' @return A list: `power` (the worst-case power), `powerVec` (one per
#'   curve), `worstCaseThetaA`, `worstCaseThetaB` (its baseline),
#'   `bootObjPower` (a [boot::boot()] object on that baseline, with
#'   `bootSe`), `nBlocks`, and the sampler's
#'   `thetaA`, `thetaB`, `stoppingTimes`, `breakVector`, `eValuesStopped`,
#'   `samplePaths`, `n1Vector`.
#' @noRd
computePowerSavi2x2 <- function(
  propDiffMin,
  na,
  nb,
  nBlocks,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  betaParameter = NULL,
  nSim = 1e3L,
  nBoot = nSim,
  seed = NULL,
  wantSamplePaths = FALSE,
  pb = TRUE
) {
  alternative <- match.arg(alternative)
  stopifnot(length(nBlocks) == 1, is.finite(nBlocks), nBlocks >= 1)

  samplingResult <- sampleStoppingTimesSavi2x2(
    propDiffMin = propDiffMin,
    na = na,
    nb = nb,
    power = NULL,
    alpha = alpha,
    alternative = alternative,
    betaParameter = betaParameter,
    nSim = nSim,
    nMax = nBlocks,
    seed = seed,
    wantSamplePaths = wantSamplePaths,
    pb = pb
  )

  # Power per curve: the fraction of paths that crossed 1 / alpha within
  # nBlocks (a never-crossing path has stopping time Inf).
  # The worst case is the smallest
  stoppingTimes <- samplingResult[["stoppingTimes"]]
  powerVec <- rowMeans(stoppingTimes <= nBlocks)
  # find the worst case for twoSided: "less" or "greater"
  worstSign <- which.min(powerVec)

  # the row index only selects its paths for the bootstrap.
  bootObjPower <- computeBootObj(
    values = stoppingTimes[worstSign, ],
    objType = "power",
    nPlan = nBlocks,
    nBoot = nBoot
  )

  list(
    "power" = powerVec[worstSign],
    "powerVec" = powerVec,
    "worstCaseThetaA" = samplingResult[["thetaA"]][worstSign],
    "worstCaseThetaB" = samplingResult[["thetaB"]][worstSign],
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
#' Runs [sampleStoppingTimesSavi2x2()] and reports its `nPlan`,
#' the `power` quantile of the stopping time at the hardest baseline
#' with a bootstrap SE and the mean stopping time of paths
#' capped at `nPlan`, both at that baseline.
#'
#' @inheritParams sampleStoppingTimesSavi2x2
#'
#' @return A list: `nPlan` (the worst-case block count, `Inf` with a warning
#'   when the worst baseline crossed too rarely), `nPlanVec` (the quantile
#'   per curve), `worstCaseThetaA`, `worstCaseThetaB` (the worst baseline),
#'   `bootObjNPlan`, `nMean`,
#'   `bootObjNMean` ([boot::boot()] objects on the worst baseline, `NULL`
#'   when `nPlan` is `Inf`), and the sampler's `thetaA`, `thetaB`,
#'   `stoppingTimes`, `breakVector`, `eValuesStopped`, `samplePaths`,
#'   `n1Vector`.
#' @noRd
computeNPlanSavi2x2 <- function(
  propDiffMin,
  na,
  nb,
  power = 0.8,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  betaParameter = NULL,
  nSim = 1e3L,
  nBoot = nSim,
  nMax = 1e4L,
  seed = NULL,
  wantSamplePaths = FALSE,
  pb = TRUE
) {
  alternative <- match.arg(alternative)
  stopifnot(!is.null(power), power > 0, power < 1)

  samplingResult <- sampleStoppingTimesSavi2x2(
    propDiffMin = propDiffMin,
    na = na,
    nb = nb,
    power = power,
    alpha = alpha,
    alternative = alternative,
    betaParameter = betaParameter,
    nSim = nSim,
    nMax = nMax,
    seed = seed,
    wantSamplePaths = wantSamplePaths,
    pb = pb
  )

  stoppingTimes <- samplingResult[["stoppingTimes"]]
  nPlan <- samplingResult[["nPlan"]]
  # The same order-statistic quantile per curve as the sampler's nPlan; its
  # largest row holds the worst baseline's paths for the bootstraps.
  nPlanVec <- apply(
    stoppingTimes,
    1,
    stats::quantile,
    probs = power,
    names = FALSE,
    type = 1
  )
  worstSign <- which.max(nPlanVec)

  # Simulation uncertainty at the worst baseline only: the bootstrap
  # quantile, and the mean stopping time with paths capped at nPlan. A
  # never-crossing path (Inf) has no finite bootstrap, so both stay NULL.
  bootObjNPlan <- NULL
  bootObjNMean <- NULL
  nMean <- NA_real_
  if (is.finite(nPlan)) {
    bootObjNPlan <- computeBootObj(
      values = stoppingTimes[worstSign, ],
      objType = "nPlan",
      power = power,
      nBoot = nBoot
    )
    bootObjNMean <- computeBootObj(
      values = stoppingTimes[worstSign, ],
      objType = "nMean",
      nPlan = nPlan,
      nBoot = nBoot
    )
    nMean <- ceiling(bootObjNMean[["t0"]])
  }

  list(
    "nPlan" = nPlan,
    "nPlanVec" = nPlanVec,
    "worstCaseThetaA" = samplingResult[["worstCaseThetaA"]],
    "worstCaseThetaB" = samplingResult[["worstCaseThetaB"]],
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
#' The smallest `propDiffMin` at which the worst-case
#' power of [computePowerSavi2x2()] at `nBlocks = nBlocksPlan` reaches
#' `power`, found by [stats::uniroot()] on `propDiffBounds`. Every candidate
#' is simulated with the same `seed`, so the target is a deterministic step
#' function of the candidate. There is no planning on `logOdds`.
#'
#' @param nBlocksPlan Planned block count at which the test is evaluated.
#' @param propDiffBounds Search interval for the magnitude of `propDiffMin`,
#'   strictly inside `(0, 1)`.
#' @param tol Tolerance of the root on the `propDiff` scale.
#' @inheritParams sampleStoppingTimesSavi2x2
#'
#' @return A single numeric: the minimal `propDiffMin`, negative for
#'   `"less"`, or `NA` when the worst-case power minus `power` has no sign
#'   change on `propDiffBounds` (still below the target at the upper bound,
#'   or already above it at the lower bound). No bootstrap object, as for
#'   [computeMinEsBatchSaviT()].
#' @noRd
computeEsMinSavi2x2 <- function(
  na,
  nb,
  nBlocksPlan,
  power = 0.8,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  betaParameter = NULL,
  nSim = 1e3L,
  seed = NULL,
  pb = TRUE,
  propDiffBounds = c(0.01, 0.9),
  tol = 1e-5
) {
  alternative <- match.arg(alternative)
  bounds <- propDiffBounds
  stopifnot(
    length(nBlocksPlan) == 1,
    is.finite(nBlocksPlan),
    nBlocksPlan >= 1,
    power > 0,
    power < 1,
    length(bounds) == 2,
    bounds[1] > 0,
    bounds[1] < bounds[2],
    bounds[2] < 1
  )

  # Worst-case power minus the target, at a candidate propDiffMin. The same
  # seed for every candidate makes this deterministic in the candidate; the
  # worst-case baseline of each curve is solved afresh inside the sampler
  # for each candidate.
  targetFunction <- function(propDiffMin) {
    computePowerSavi2x2(
      propDiffMin = propDiffMin,
      na = na,
      nb = nb,
      nBlocks = nBlocksPlan,
      alpha = alpha,
      alternative = alternative,
      betaParameter = betaParameter,
      nSim = nSim,
      seed = seed,
      pb = pb
    )[["power"]] -
      power
  }

  # No sign change means the target is out of reach on the bracket (or
  # already met at its lower end); report NA rather than a spurious edge.
  targetAtBounds <- c(targetFunction(bounds[1]), targetFunction(bounds[2]))
  if (targetAtBounds[1] >= 0 || targetAtBounds[2] < 0) {
    return(NA_real_)
  }

  esMin <- stats::uniroot(
    targetFunction,
    interval = bounds,
    f.lower = targetAtBounds[1],
    f.upper = targetAtBounds[2],
    tol = tol
  )[["root"]]
  # The search is on the magnitude; the sign follows the alternative.
  if (alternative == "less") -esMin else esMin
}
