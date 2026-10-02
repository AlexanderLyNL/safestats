# Test functions ----

#' UMP conditional e-value of one 2x2 table
#'
#' The e-value of the conditional test on one table: the likelihood of
#' [logLikelihoodFNCH()] at the uniformly most powerful log odds ratio of
#' [solveUmpLogOdds()] over the hypergeometric null. `"twoSided"` averages
#' the two one-sided e-values. The test uses it in place of block 1's
#' plain factor.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#'
#' @return A numeric, or `NULL` when no one-sided UMP log odds ratio exists
#'   at level `alpha`.
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

#' Safe Anytime-Valid Test of Two Proportions
#'
#' Tests `thetaA = thetaB` on a stream of 2x2 tables, one per block, with
#' the e-process chosen by `designObj[["eType"]]`. Both effects are A
#' minus B, `propDiff = thetaA - thetaB` and `logOdds = logit(thetaA) -
#' logit(thetaB)`, so `"greater"` means group A has the larger proportion,
#' as in [stats::t.test()]. Block 1 of every e-process is replaced by the
#' UMP conditional e-value of [savi2x2TestStatUmp()] when one exists at
#' level `alpha`; otherwise it keeps its plain factor, with a warning.
#'
#' - `"eBeta"` (`propDiff`): independent Beta posterior means of `thetaA`
#'   and `thetaB` from the previous blocks, against the pooled null mean.
#' - `"grow"` (`propDiff` or `logOdds`, whichever minimal effect the design
#'   holds): the alternative restricted to the signed minimal effect;
#'   `"twoSided"` averages the cumulative e-processes of both signs.
#' - `"eGauss"` (`logOdds`): the conditional likelihood of Fisher's
#'   noncentral hypergeometric distribution mixed under the design's Normal
#'   prior on a `logOdds` grid, against the hypergeometric null. The
#'   mixture may include the current block because it conditions on the
#'   block's total.
#'
#' Only `"eBeta"` and `"eGauss"` give a confidence interval or sequence,
#' each on its own effect.
#'
#' @param ya positive observations/ events per data block in group A: a
#'   numeric with integer values between (and including) 0 and `na`, the
#'   number of observations in group A per block.
#' @param yb positive observations/ events per data block in group B: a
#'   numeric with integer values between (and including) 0 and `nb`, the
#'   number of observations in group B per block.
#' @param designObj an object obtained from [designSavi2x2()], which also
#'   supplies `na` and `nb`.
#' @param wantCi default `FALSE`, compute a confidence interval.
#' @param wantConfidenceSequence logical that can be set to true when the
#'   user wants a savi confidence sequence to be estimated, one row per
#'   block; takes precedence over `wantCi`.
#' @param ciValue numeric representing the confidence level.
#'   Default ciValue=NULL yields ciValue = 1 - alpha
#'
#' @return Returns an object of class 'saviTest'. An object of class 'saviTest'
#' is a list containing at least the following components:
#'
#' \describe{
#'   \item{n}{The realised sample size(s):
#'   `c(na = sum(na), nb = sum(nb), nBlocks = length(ya))`.}
#'   \item{eValue}{the e-value of the savi test on all blocks; reject when
#'   it is at least `1 / alpha`.}
#'   \item{eValueVec}{the realised e-values after each block, element `i`
#'   on blocks `1..i`.}
#'   \item{n1Vec}{the block index, used for plotting.}
#'   \item{estimate}{the estimated proportions `c(thetaA, thetaB)`, pooled
#'   over all blocks.}
#'   \item{confSeq}{a savi confidence interval for the design's effect at
#'   level `ciValue` on all blocks; `NA` when the set is empty, `NULL` for
#'   "grow" or when no interval was requested.}
#'   \item{confSeqMatrix}{with `wantConfidenceSequence`, the savi confidence
#'   sequence as an `nBlocks x 2` matrix, row `i` the interval on blocks
#'   `1..i`; `confSeq` is its last row. With `runningIntersection` the
#'   rows are nested and stay `NA` after the first empty row.}
#'   \item{ciValue}{the confidence level used.}
#'   \item{betaParameter}{for "eBeta" only, `list(betaA1, betaA2, betaB1,
#'   betaB2)`, the Beta prior updated with all blocks.}
#'   \item{alternative}{any of "twoSided", "greater", "less" copied from
#'   the design.}
#'   \item{h0}{the null value copied from the design.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "saviDesign" described in
#'   [designSavi2x2()].}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(turner2024generic)`
#'   `r addCite(turner2023exact)`
#'
#' @export
#'
#' @examples
#' designObj <- designSavi2x2(na = 10, nb = 10, eType = "eBeta")
#' ya <- c(8, 7, 9, 6)
#' yb <- c(4, 5, 3, 5)
#' savi2x2TestStat(ya, yb, designObj = designObj)
savi2x2TestStat <- function(
  ya,
  yb,
  designObj = NULL,
  wantCi = FALSE,
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
    warning("There is only 1 table, switched to UMP conditional e-variable")
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

  # Prior: a user-supplied Beta prior keeps its own block 1; the default
  # prior has block 1 replaced by the UMP conditional e-value.
  replaceUmpAtFirstBlock <- TRUE
  if (eType == "eBeta") {
    if (!is.null(betaParameter)) {
      replaceUmpAtFirstBlock <- FALSE
      message("betaParameter given: block 1 keeps the eBeta e-value")
    } else {
      # Default shapes 1 / (2 n): a vague prior on the scale of one block. The
      # prior only acts before block 1's data, so varying sizes take block 1's.
      if (any(na != na[1]) || any(nb != nb[1])) {
        warning(
          "betaParameter defaults to 1 / (2 * na[1]) and 1 / (2 * nb[1]) on ",
          "the first block's sizes na = ",
          na[1],
          ", nb = ",
          nb[1]
        )
      }
      betaParameter <- list(
        "betaA1" = 1 / (2 * na[1]),
        "betaA2" = 1 / (2 * na[1]),
        "betaB1" = 1 / (2 * nb[1]),
        "betaB2" = 1 / (2 * nb[1])
      )
    }
  }
  if (eType == "eGauss") {
    # Gaussian prior on logOdds: N(mean, sd), restricted to the helper's grid
    # (-20, 20) and to the side of a one-sided alternative.
    if (is.null(gaussParameter)) {
      gaussParameter <- list("mean" = 0, "sd" = 1)
    }
    result[["gaussParameter"]] <- gaussParameter
  }

  # Compute: eValueVec ----
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
      esMin,
      alpha,
      alternative,
      betaParameter
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
  # eGauss only: the numerator of the plain process on blocks 1..i, the
  # log e-process plus the hypergeometric log likelihood it was divided by.
  # The logOdds interval and sequence invert it, so it is kept before the
  # UMP replacement of block 1.
  logNumerator <- if (eType == "eGauss") {
    logEValueVec + cumsum(stats::dhyper(ya, na, nb, ya + yb, log = TRUE))
  }
  # Replace block 1 by the UMP e-value, unless a custom Beta prior opted out.
  # logEValueVec[1] should be log(1) = 0 but write it out for clarity
  # No UMP e-value at this alpha (NULL): block 1 keeps its plain factor.
  if (replaceUmpAtFirstBlock) {
    eValueUmp <- savi2x2TestStatUmp(
      ya[1],
      yb[1],
      na[1],
      nb[1],
      alpha,
      alternative
    )
    if (is.null(eValueUmp)) {
      warning(
        "no UMP e-value exists for the first block at alpha = ",
        alpha,
        "; block 1 keeps the plain e-value"
      )
    } else {
      logEValueVec <- logEValueVec - logEValueVec[1] + log(eValueUmp)
    }
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
        logNumerator,
        1 - ciValue,
        runningIntersection
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
  # eBeta only: the Beta posterior of the observed blocks, the prior a further
  # block would use. grow's prior lives on its curve, with no conjugate
  # update to store; eGauss holds no Beta prior.
  if (eType == "eBeta") {
    result[["betaParameter"]] <- list(
      "betaA1" = betaParameter[["betaA1"]] + sum(ya),
      "betaA2" = betaParameter[["betaA2"]] + sum(na) - sum(ya),
      "betaB1" = betaParameter[["betaB1"]] + sum(yb),
      "betaB2" = betaParameter[["betaB2"]] + sum(nb) - sum(yb)
    )
  }
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

#' Design a Safe Anytime-Valid Experiment to Test Two Proportions
#'
#' A designed experiment requires (1) a number of data blocks nBlocksPlan
#' to plan for, and (2) a savi test defining parameter: `eType` and, with
#' it, the effect measure. Both effects are A minus B and anchored on
#' `thetaB`: `propDiff = thetaA - thetaB` and `logOdds = logit(thetaA) -
#' logit(thetaB)`, so `"greater"` means group A has the larger proportion,
#' as `x - y > 0` does in [stats::t.test()].
#'
#' - `"eBeta"` (`propDiff`) and `"eGauss"` (`logOdds`) are unrestricted: a
#'   minimal effect is an error, a `power` is dropped with a warning, and
#'   `nBlocksPlan` is kept as the planned block count. Each reads one
#'   prior, `betaParameter` or `gaussParameter`; `NULL` keeps the default,
#'   which [savi2x2TestStat()] fills in. A one-sided `alternative` acts in
#'   the UMP first block, and for `"eGauss"` also restricts the prior to
#'   that side.
#' - `"grow"` plugs in exactly one of `propDiffMin`, `logOddsMin` as the
#'   fixed alternative and reads no prior. As in [designSaviZ()], a value
#'   whose sign contradicts a one-sided `alternative` is flipped with a
#'   warning; `"twoSided"` uses the magnitude. Zero is an error.
#'
#' For `"grow"` on `propDiff` the design involves alpha and the three
#' quantities: (1) nBlocksPlan, (2) power, and (3) a minimal relevant
#' difference in proportions propDiffMin, simulated at the worst-case
#' baseline of [sampleStoppingTimesSavi2x2()].
#' \describe{
#'   \item{Scenario 1a}{Goal: "nBlocksPlan" and optimal E-variable. Given: propDiffMin and power.}
#'   \item{Scenario 2}{Goal: "power" and optimal E-variable. Given: propDiffMin and nBlocksPlan.}
#'   \item{Scenario 3}{Goal: "propDiffMin" and optimal E-variable. Given: power and nBlocksPlan.}
#' }
#' A minimal effect alone designs without simulation. Any other
#' combination is an error, except that `logOddsMin` drops a `power` with
#' a warning. Per-block size vectors supply `nBlocksPlan` when it is `NULL`.
#' Other arguments are taken as given; only `alpha` and `power` in
#' `(0, 1)`, equal lengths of `na` and `nb`, at most one minimal effect,
#' and the minimal effects' ranges are checked.
#'
#' @param na number of observations in group A per data block, one value
#'   or one per block.
#' @param nb number of observations in group B per data block, one value
#'   or one per block.
#' @param nBlocksPlan planned number of data blocks collected, see
#'   scenario 2 and 3 above.
#' @param propDiffMin numeric in (-1, 1) that defines the minimal relevant
#'   difference in proportions `thetaA - thetaB`, the smallest difference
#'   that we would like to detect (with sufficient power).
#' @param logOddsMin numeric that defines the minimal relevant log odds
#'   ratio `logit(thetaA) - logit(thetaB)`.
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error
#'   and the null rejection rule e >= 1/alpha.
#' @param power numeric in (0, 1) that specifies the desired power, that
#'   is, the targetted chance to stop in favour of the alternative over the
#'   null hypothesis, when the alternative holds true.
#' @param h0 numeric, representing the null value, default h0=0. Only
#'   `h0 = 0` is currently supported.
#' @param alternative a character string specifying the alternative
#'   hypothesis. Must be one of "twoSided" (default), "greater" or "less",
#'   where "greater" means that group A has the larger proportion.
#' @param eType character one of "eBeta", "grow", and "eGauss". "eBeta" is
#'   default and uses a Beta prior on each proportion, "grow" uses a point
#'   prior at the minimal effect, "eGauss" a normal prior on the log odds
#'   ratio.
#' @param betaParameter `list(betaA1, betaA2, betaB1, betaB2)`, the Beta
#'   prior shapes on `thetaA` and `thetaB` for "eBeta", each a single
#'   finite positive number; `NULL` means `1 / (2 * na)` and `1 / (2 * nb)`,
#'   taking block 1's sizes with a warning when they vary by block. A given
#'   prior keeps its own block 1: the test does not replace it by the UMP
#'   e-value. For "grow" on propDiff, the test puts `Beta(betaA1, betaA2)` on
#'   `thetaA` rescaled to its feasible interval on the curve, B's shapes
#'   unused, and `NULL` means uniform; planning ignores it.
#' @param gaussParameter `list(mean, sd)`, the Normal prior on `logOdds`
#'   for "eGauss", restricted to the grid `(-20, 20)` and to the side of a
#'   one-sided `alternative`; `NULL` means `list(mean = 0, sd = 1)`.
#' @param runningIntersection logical, if `TRUE` then intersect each row
#'   of the blockwise confidence sequence with the previous one; `NULL`
#'   keeps the constructor's `FALSE`.
#' @param nSim integer > 0, the number of simulations needed to compute
#'   power or the number of samples paths for the savi 2x2 test under
#'   continuous monitoring.
#' @param nBoot integer > 0 representing the number of bootstrap samples
#'   to assess the accuracy of the approximations of the power, or the
#'   number of blocks for the savi 2x2 test under continuous monitoring.
#' @param nMax integer > 0, maximum number of data blocks in each sample path.
#' @param seed integer, seed number. Default seed=NULL yields seed=2026.
#' @param wantSamplePaths logical, if `TRUE` then also outputs the sample paths.
#' @param pb logical, if `TRUE`, then show progress bar.
#'
#' @return Returns a saviDesign object that includes:
#'
#' \describe{
#'   \item{nPlan}{the planned sample size(s): `list(na, nb)`, with
#'   `nBlocksPlan` when given or planned.}
#'   \item{parameter}{for "eBeta" and "eGauss", the prior as one named
#'   string, or `"default"` when `NULL`; absent for "grow", whose minimal
#'   effect is `esMin`.}
#'   \item{esMin}{the minimal relevant effect size provided by the user, or
#'   found in scenario 3, signed and named `propDiff` or `logOdds`; `NULL`
#'   for "eBeta" and "eGauss".}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{power}{the desired power provided by the user, or the worst-case
#'   power found in scenario 2.}
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
#'   \item{eType}{any of "eBeta", "grow", "eGauss" provided by the user.}
#'   \item{h0}{the null value, 0.}
#'   \item{betaParameter}{for "eBeta" and "grow" on propDiff, the Beta prior
#'   given, `NULL` for the default.}
#'   \item{gaussParameter}{for "eGauss", the Normal prior given, `NULL` for
#'   the default.}
#'   \item{runningIntersection}{logical, as provided by the user.}
#'   \item{designScenario}{"1a", "2" or "3" when planned.}
#'   \item{nPlanTwoSe, bootObjNBlocksPlan, nMean, nMeanTwoSe, bootObjNMean}{
#'   scenario 1a: the bootstrap summaries of the planned block count and of
#'   the mean stopping time.}
#'   \item{powerTwoSe, bootObjPower}{scenario 2: the bootstrap summaries of
#'   the power.}
#'   \item{worstCaseThetaA, worstCaseThetaB, breakVector, samplePaths}{
#'   scenarios 1a and 2: the worst baseline and its simulated paths.}
#'   \item{testType}{here 2x2}
#'   \item{testName}{"Two Proportions".}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(turner2024generic)`
#'
#' @export
#'
#' @examples
#' # Unrestricted test on propDiff with the default Beta prior
#' designSavi2x2(na = 10, nb = 10, eType = "eBeta")
#'
#' # Grow on a minimal difference, no planning
#' designSavi2x2(na = 10, nb = 10, propDiffMin = 0.2, eType = "grow")
#'
#' # Scenario 2: worst-case power at a planned block count
#' designSavi2x2(na = 10, nb = 10, propDiffMin = 0.3, nBlocksPlan = 8,
#'               eType = "grow", nSim = 50, pb = FALSE)
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

  # Arguments are taken as given; only cheap inconsistencies are caught.
  stopifnot(
    "alpha must lie in (0, 1)" = alpha > 0 && alpha < 1,
    "power must lie in (0, 1)" = is.null(power) || (power > 0 && power < 1),
    "na and nb must have the same length" = length(na) == length(nb),
    "supply propDiffMin or logOddsMin, not both" =
      is.null(propDiffMin) || is.null(logOddsMin),
    "propDiffMin must be nonzero in (-1, 1)" =
      is.null(propDiffMin) || (propDiffMin != 0 && abs(propDiffMin) < 1),
    "logOddsMin must be a nonzero finite number" =
      is.null(logOddsMin) || (is.finite(logOddsMin) && logOddsMin != 0)
  )
  # Fix nBlocksPlan if not
  if (length(na) > 1L && is.null(nBlocksPlan)) {
    nBlocksPlan <- length(na)
  }

  result <- constructSaviDesignObj("Two Proportions")

  # Fill: result ----
  if (!is.null(runningIntersection)) {
    result[["runningIntersection"]] <- runningIntersection
  }
  result[["eType"]] <- eType
  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["h0"]] <- 0

  # Minimal effect: signed by the alternative. A sign contradicting a
  # one-sided alternative is flipped with a warning, twoSided takes the
  # magnitude. Named so print() and plot() show which effect it is on.
  esMin <- NULL
  if (!is.null(propDiffMin)) {
    propDiffMin <- checkAndReturnEsMinParameterSide(
      propDiffMin,
      alternative,
      "propDiffMin"
    )
    esMin <- c("propDiff" = propDiffMin)
  }
  if (!is.null(logOddsMin)) {
    logOddsMin <- checkAndReturnEsMinParameterSide(
      logOddsMin,
      alternative,
      "logOddsMin"
    )
    esMin <- c("logOdds" = logOddsMin)
  }
  result[["esMin"]] <- esMin

  # eBeta and eGauss: unrestricted, no planning, one prior each ----
  if (eType != "grow" && !is.null(esMin)) {
    stop(
      "a minimal effect needs eType = 'grow'; eType = '",
      eType,
      "' is unrestricted"
    )
  }
  if (eType != "grow" && !is.null(power)) {
    warning("eType = '", eType, "' has no planning: power is dropped")
    power <- NULL
  }

  # The prior as given; NULL means the default, which savi2x2TestStat fills
  # in on the data's block sizes.
  if (eType == "eBeta") {
    result[["betaParameter"]] <- betaParameter
  }
  if (eType == "eGauss") {
    result[["gaussParameter"]] <- gaussParameter
  }
  # grow on propDiff: the test's prior on the curve, uniform when NULL; the
  # planners ignore it.
  if (eType == "grow" && is.null(logOddsMin)) {
    result[["betaParameter"]] <- betaParameter
  }

  # grow: a minimal effect, or power and nBlocksPlan to find one ----
  # Planning exists on propDiff only (design: planning); each branch below
  # is a row of its table.
  planning <- NULL
  if (eType == "grow") {
    if (!is.null(logOddsMin) && !is.null(power)) {
      warning("no planning on logOdds: power is dropped")
      power <- NULL
    }
    if (!is.null(propDiffMin) && !is.null(power) && is.null(nBlocksPlan)) {
      # Scenario 1a: propDiffMin + power -> worst-case stopping time
      planning <- computeNPlanSavi2x2(
        propDiffMin = propDiffMin,
        na = na,
        nb = nb,
        power = power,
        alpha = alpha,
        alternative = alternative,
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
    } else if (
      !is.null(propDiffMin) && is.null(power) && !is.null(nBlocksPlan)
    ) {
      # Scenario 2: propDiffMin + nBlocksPlan -> worst-case power
      planning <- computePowerSavi2x2(
        propDiffMin = propDiffMin,
        na = na,
        nb = nb,
        nBlocks = nBlocksPlan,
        alpha = alpha,
        alternative = alternative,
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
    } else if (is.null(esMin) && !is.null(power) && !is.null(nBlocksPlan)) {
      # Scenario 3: power + nBlocksPlan -> the minimal detectable propDiff
      esMin <- computeEsMinSavi2x2(
        na = na,
        nb = nb,
        nBlocksPlan = nBlocksPlan,
        power = power,
        alpha = alpha,
        alternative = alternative,
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
    } else if (!is.null(esMin) && is.null(power)) {
      # A minimal effect alone, or logOddsMin with nBlocksPlan: no simulation
    } else {
      stop(
        "can't design with eType = 'grow': give a minimal effect alone, ",
        "propDiffMin with power or with nBlocksPlan, or power and ",
        "nBlocksPlan without a minimal effect; per-block sizes fix nBlocksPlan"
      )
    }
    # Scenarios 1a and 2 keep the worst baseline and its simulated paths.
    if (!is.null(planning)) {
      result[["worstCaseThetaA"]] <- planning[["worstCaseThetaA"]]
      result[["worstCaseThetaB"]] <- planning[["worstCaseThetaB"]]
      result[["breakVector"]] <- planning[["breakVector"]]
      result[["samplePaths"]] <- planning[["samplePaths"]]
    }
  }

  # eBeta and eGauss only: the prior as one named string, names and values
  # comma-separated, so the shared print methods show it side by side as
  # they do nPlan; a NULL prior shows as "default". grow has none: its
  # defining quantity is esMin, which the print methods already show.
  if (eType != "grow") {
    prior <- if (eType == "eBeta") betaParameter else gaussParameter
    result[["parameter"]] <- if (is.null(prior)) {
      stats::setNames(
        "default",
        if (eType == "eBeta") "betaParameter" else "gaussParameter"
      )
    } else {
      stats::setNames(
        paste(vapply(prior, format, character(1)), collapse = ", "),
        paste(names(prior), collapse = ", ")
      )
    }
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

#' Helper function: Computes the savi confidence interval for propDiff
#'
#' Inverts the eBeta e-process: the predictive Beta numerator of
#' [predictiveThetas2x2()] against, per block, the reverse information
#' projection of [solveRIPr2x2PropDiff()] onto each candidate `propDiff`.
#' The interval is the set of candidates whose e-value stays below
#' `1 / alpha`, found by one minimum and two roots on `domain`.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#' @param alpha numeric in (0, 1), one minus the confidence level.
#' @param domain `c(lower, upper)`, the candidates searched.
#'
#' @return numeric vector that contains the lower and upper bound of the
#'   savi confidence interval, named `lowerBound`, `upperBound`; the
#'   `domain` edge when that edge is still inside, and `c(-1, 1)` with a
#'   warning when the set is empty.
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

#' Helper function: Computes the savi confidence sequence for propDiff
#'
#' Row `i` is the interval of [computeConfidenceInterval2x2PropDiff()] on
#' blocks `1..i`, found by advancing one e-process per grid candidate and
#' reading off the run of candidates below `1 / alpha`, each bound refined
#' outward by a secant.
#'
#' @inheritParams computeConfidenceInterval2x2PropDiff
#' @inheritParams designSavi2x2
#' @param nGrid integer > 0, the number of candidates strictly inside `(-1, 1)`.
#'
#' @return an `nBlocks x 2` matrix that contains the lower and upper bound
#'   of the savi confidence sequence per block; an empty set is an `NA`
#'   row, and with `runningIntersection` the rows are nested and stay `NA`
#'   after the first empty row.
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

#' Helper function: Computes the savi confidence interval for logOdds
#'
#' Inverts the eGauss numerator: the interval is the set of `logOdds`
#' values at which the numerator over the conditional likelihood of
#' [logLikelihoodFNCH()] stays below `1 / alpha`, found by one minimum and
#' two roots on `domain`.
#'
#' @inheritParams computeConfidenceInterval2x2PropDiff
#' @param logPTotal numeric, the eGauss log numerator on all blocks: the
#'   plain cumulative log e-value plus the cumulative hypergeometric log
#'   likelihood.
#'
#' @return numeric vector that contains the lower and upper bound of the
#'   savi confidence interval, named `lowerBound`, `upperBound`; `NA` with
#'   a warning for a bound that lies beyond `domain`, and the whole
#'   `domain` with a warning when the set is empty.
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

#' Helper function: Computes the savi confidence sequence for logOdds
#'
#' Row `i` is the interval of [computeConfidenceInterval2x2LogOdds()] on
#' blocks `1..i`. The eGauss numerator at block `i` is that of the test on
#' blocks `1..i`, as every factor is fixed given its block's total or
#' predictable, so the test's numerator vector serves all rows; the root
#' finding is repeated per prefix.
#'
#' @inheritParams computeConfidenceInterval2x2LogOdds
#' @inheritParams designSavi2x2
#' @param logNumerator numeric vector, the eGauss log numerator on blocks
#'   `1..i` for element `i`: the plain cumulative log e-process, before the
#'   UMP replacement of block 1, plus the cumulative hypergeometric log
#'   likelihood.
#'
#' @return an `nBlocks x 2` matrix that contains the lower and upper bound
#'   of the savi confidence sequence per block; an empty set is an `NA`
#'   row, and with `runningIntersection` the rows stay `NA` after the first
#'   empty row.
computeConfidenceSequence2x2LogOdds <- function(
  ya,
  yb,
  na,
  nb,
  logNumerator,
  alpha,
  runningIntersection
) {
  # THIS IS SLOW, WE MIGHT WANT TO STOP EARLY!
  nBlocks <- length(ya)
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

# TODO: is there a betaA1 * na term? check with Peter

#' Predictive Beta posterior means of thetaA and thetaB
#'
#' Independent Beta posteriors on `thetaA` and `thetaB` with prior shapes
#' `betaParameter`, updated with the blocks before each block: block `i`
#' uses blocks `1..i-1` only, and block 1 the prior.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#'
#' @return A list of numeric vectors `thetaA` and `thetaB`, one value per block.
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

#' Plain eBeta e-process for propDiff
#'
#' The cumulative log likelihood ratio of the predictive Beta means of
#' [predictiveThetas2x2()] over the pooled null mean
#' `(na * thetaA + nb * thetaB) / (na + nb)`, block 1 included.
#'
#' @inheritParams predictiveThetas2x2
#'
#' @return A numeric vector, element `i` the log e-value on blocks `1..i`.
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

#' Plain grow e-process for propDiff
#'
#' The alternative is restricted to the curve `thetaA = thetaB +
#' propDiffMin`, signed by `alternative`, which leaves one free parameter;
#' its posterior on a grid of `nWeight` points under a Beta prior, uniform
#' by default, gives the predictable plug-in for each block, against the
#' pooled null mean.
#' `"twoSided"` runs both signs of the magnitude and averages their
#' cumulative e-processes.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#' @param betaParameter `NULL` or a list with `betaA1`, `betaA2`: the prior
#'   `Beta(betaA1, betaA2)` on `thetaA` rescaled to its feasible interval on
#'   the curve; B's shapes are unused. `NULL` means `Beta(1, 1)`, uniform.
#' @param earlyStopping logical, if `TRUE` cut the vector at the first
#'   crossing of `1 / alpha`.
#' @param nWeight integer > 0, the number of grid points.
#'
#' @return A numeric vector, element `i` the log e-value on blocks `1..i`.
logEValueVec2x2PropDiffGrow <- function(ya, yb, na, nb,
  propDiffMin, alpha, alternative = c("twoSided", "greater", "less"),
  betaParameter = NULL,
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

  # Prior: Beta(1, 1) on rho, uniform prior on the line
  # A given betaParameter puts Beta(betaA1, betaA2) on rho instead. thetaA
  # increases with rho on both curves, so this is A's prior rescaled.
  if (is.null(betaParameter)) {
    betaParameter <- list("betaA1" = 1, "betaA2" = 1)
  }
  logWeights <- stats::dbeta(
    rho,
    betaParameter[["betaA1"]],
    betaParameter[["betaA2"]],
    log = TRUE
  )
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

# TODO: this can be a cubic function solver

#' Reverse information projection onto a fixed propDiff
#'
#' The null pair `(nullThetaA, nullThetaA - propDiff)` closest in KL
#' divergence to `(thetaA, thetaB)` for one block: the root of the
#' derivative of `na KL(thetaA || nullThetaA) + nb KL(thetaB || nullThetaB)`.
#'
#' @inheritParams designSavi2x2
#' @param thetaA,thetaB numeric in (0, 1), the alternative's proportions.
#' @param propDiff numeric in (-1, 1), the null difference `thetaA - thetaB`.
#' @param tol numeric > 0, the tolerance of [stats::uniroot()], also kept
#'   from the edges of the search interval.
#'
#' @return A numeric, `nullThetaA`.
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

#' Worst-case baseline of the propDiff grow test
#'
#' The pair on the curve of each alternative at which the e-process grows
#' slowest, so a block count planned there holds for every baseline. The
#' expected log e-increment per block is the KL divergence of the truth
#' from its pooled null projection,
#' ```
#' R(theta) = nLow KL(theta || theta0) + nHigh KL(theta + d || theta0),
#' theta0   = theta + nHigh d / (nLow + nHigh),
#' ```
#' with `theta` the lower proportion, `nLow` its group size and `nHigh` the
#' size of the group at `theta + d`. Its minimiser is the root of
#' ```
#' nLow logit(theta) + nHigh logit(theta + d) = (nLow + nHigh) logit(theta0),
#' ```
#' summed over blocks when the sizes vary per block. Equal sizes give
#' exactly the midpoint `(1 - d) / 2`; in general the root is
#' `(1 - d) / 2 + d (nHigh - nLow) / (6 (nLow + nHigh)) + O(d^3)`. The rate
#' ignores the first-block UMP factor, which moves the simulated worst case
#' when one group is very small.
#'
#' `"greater"`: `thetaA = thetaB + d`, so `thetaB` is the lower proportion
#' and `nLow = nb`. `"less"`: `thetaA = thetaB - d`, so `thetaA` is the
#' lower one and `nLow = na`. `"twoSided"` gives both, greater first.
#'
#' @inheritParams designSavi2x2
#'
#' @return A data frame with one row per curve and columns `thetaA`, `thetaB`.
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

# TODO: what if ya, yb, na, nb are vector

#' Conditional log likelihood of 2x2 tables at one logOdds
#'
#' Given its total `ya + yb`, `ya` is Fisher's noncentral hypergeometric
#' with odds `exp(logOdds)` on group A; `logOdds = 0` is the hypergeometric
#' null. Evaluated on the log scale, binomial coefficients plus `logOdds`
#' times the feasible count normalised by log-sum-exp, so extreme `logOdds`
#' give finite values where the plain density underflows to 0.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#' @param logOdds numeric, one log odds ratio `logit(thetaA) - logit(thetaB)`.
#'
#' @return A numeric vector, one log density per block.
logLikelihoodFNCH <- function(ya, yb, na, nb, logOdds) {
  mapply(
    function(ya, yb, na, nb) {
      k <- max(0, ya + yb - nb):min(na, ya + yb)
      logTerms <- lchoose(na, k) + lchoose(nb, ya + yb - k) + logOdds * k
      shift <- max(logTerms)
      logTerms[ya - k[1] + 1] - shift - log(sum(exp(logTerms - shift)))
    },
    ya = ya,
    yb = yb,
    na = na,
    nb = nb
  )
}

#' Plain grow e-process for logOdds
#'
#' The cumulative conditional log likelihood ratio of [logLikelihoodFNCH()]
#' at `logOddsMin`, signed by `alternative`, over the hypergeometric null.
#' `"twoSided"` runs both signs of the magnitude and averages their
#' cumulative e-processes.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#'
#' @return A numeric vector, element `i` the log e-value on blocks `1..i`.
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

#' Plain eGauss e-process for logOdds
#'
#' The conditional likelihood of Fisher's noncentral hypergeometric
#' distribution on blocks `1..i`, mixed under the Normal prior
#' `gaussParameter` on `logOddsGrid` restricted to the side of a one-sided
#' `alternative`, over the hypergeometric null. The mixture for block `i`
#' includes block `i`, which is valid because each factor conditions on
#' its block's total.
#'
#' @inheritParams savi2x2TestStat
#' @inheritParams designSavi2x2
#' @param logOddsGrid numeric vector, the grid the prior is normalised on.
#'
#' @return A numeric vector, element `i` the log e-value on blocks `1..i`.
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
#' @inheritParams designSavi2x2
#' @inheritParams logLikelihoodFNCH
#' @param totalSuccesses nonnegative integer, the total successes `ya + yb`
#'   in the table, at most `na + nb`.
#'
#' @return A numeric, the log of the sum over the feasible `k` of
#'   `choose(na, k) * choose(nb, totalSuccesses - k) * exp(logOdds * k)`.
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

#' Uniformly most powerful log odds ratio of a one-sided conditional test
#'
#' Solves the conditional KL divergence of the alternative at `logOdds`
#' from `logOddsNull` equal to `log(1 / alpha)`, searching up to
#' `searchBound` from `logOddsNull` in the direction of `alternative`.
#'
#' @inheritParams designSavi2x2
#' @inheritParams fnchLogPartition
#' @param alternative `"greater"` or `"less"`, the side searched.
#' @param logOddsNull numeric, the null log odds ratio.
#' @param searchBound numeric > 0, the width of the search interval.
#'
#' @return A numeric, or `NULL` with a warning when the target is not
#'   reached on the search interval.
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

#' Simulate stopping times for the savi 2x2 test
#'
#' Draws `nSim` streams of tables on the curve `thetaA = thetaB +
#' |propDiffMin|` for `"greater"`, `thetaA = thetaB - |propDiffMin|` for
#' `"less"`, and on both curves for `"twoSided"`, at the worst-case
#' baseline of each curve from [solveWorstCaseTheta2x2PropDiff()], and
#' records when the plain grow e-process first crosses `1 / alpha`.
#' `nPlan` is the worst `power` quantile of the stopping time over the
#' curves; with `power = NULL` it is skipped. There is no planning on
#' `logOdds`.
#'
#' @inheritParams designSavi2x2
#' @param eType "grow", the only e-variable with planning.
#' @param nBoot integer > 0, the number of bootstrap samples; not used by
#'   the sampler itself, see [computePowerSavi2x2()] and
#'   [computeNPlanSavi2x2()].
#' @param wantEValuesAtNMax logical. If `TRUE`, then compute eValues at
#'   nMax. Default `FALSE`. Not yet implemented.
#' @param wantSimData logical. If `TRUE`, then output the simulated data.
#'   Not yet implemented.
#'
#' @return a list with `thetaA`, `thetaB` (one per curve) and, one row per
#'   curve, `stoppingTimes` and `breakVector`. Entries of `breakVector` are
#'   0, 1. A 1 represents stopping due to exceeding `nMax`, which implies
#'   that the corresponding stopping time is `Inf`, and 0 due to `1/alpha`
#'   threshold crossing. Further `eValuesStopped`, `samplePaths` (a list of
#'   `nSim x nMax` sparse matrices, or `NULL`), `n1Vector` (the block
#'   index), `nPlan`, and the baseline it was taken at, `worstCaseThetaA`,
#'   `worstCaseThetaB`.
sampleStoppingTimesSavi2x2 <- function(
  propDiffMin,
  na,
  nb,
  power = NULL,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
  eType = c("grow"),
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


#' Helper function: Computes the power of the savi 2x2 test based on propDiffMin and nBlocks
#'
#' Runs [sampleStoppingTimesSavi2x2()] with `nMax = nBlocks` and reports,
#' per curve, the fraction of paths that cross `1 / alpha` within
#' `nBlocks`; the worst case is the smallest, since the test must meet its
#' target whatever the baseline is.
#'
#' @inheritParams designSavi2x2
#' @param nBlocks integer > 0, the planned number of data blocks.
#'
#' @return a list which contains at least `power` (the worst case over the
#'   curves) and an adapted bootObject `bootObjPower` of class
#'   [boot::boot()] on that baseline; further `powerVec` (one per curve),
#'   `worstCaseThetaA`, `worstCaseThetaB`, `nBlocks`, and the sampler's
#'   `thetaA`, `thetaB`, `stoppingTimes`, `breakVector`, `eValuesStopped`,
#'   `samplePaths`, `n1Vector`.
computePowerSavi2x2 <- function(
  propDiffMin,
  na,
  nb,
  nBlocks,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
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


#' Helper function: Computes nBlocksPlan based on propDiffMin, alpha and power
#'
#' Runs [sampleStoppingTimesSavi2x2()] and reports its `nPlan`, the
#' `power` quantile of the stopping time at the hardest baseline, with a
#' bootstrap SE and the mean stopping time of paths capped at `nPlan`,
#' both at that baseline.
#'
#' @inheritParams designSavi2x2
#'
#' @return a list which contains at least `nPlan` (the worst-case block
#'   count, `Inf` with a warning when the worst baseline crossed too
#'   rarely) and adapted bootObjects `bootObjNPlan` and `bootObjNMean` of
#'   class [boot::boot()] on the worst baseline (`NULL` when `nPlan` is
#'   `Inf`); further `nPlanVec` (the quantile per curve), `nMean`,
#'   `worstCaseThetaA`, `worstCaseThetaB`, and the sampler's `thetaA`,
#'   `thetaB`, `stoppingTimes`, `breakVector`, `eValuesStopped`,
#'   `samplePaths`, `n1Vector`.
computeNPlanSavi2x2 <- function(
  propDiffMin,
  na,
  nb,
  power = 0.8,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
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


#' Computes the smallest detectable propDiffMin with power probability, for the provided number of blocks
#'
#' The smallest `propDiffMin` at which the worst-case power of
#' [computePowerSavi2x2()] at `nBlocks = nBlocksPlan` reaches `power`,
#' found by [stats::uniroot()] on `propDiffBounds`. Every candidate is
#' simulated with the same `seed`, so the target is a deterministic step
#' function of the candidate. There is no planning on `logOdds`.
# FIXME: should i use different seed?
#'
#' @inheritParams designSavi2x2
#' @param propDiffBounds search interval for the magnitude of `propDiffMin`,
#'   strictly inside `(0, 1)`.
#' @param tol tolerance of the root on the `propDiff` scale.
#'
#' @return numeric that represents the minimal detectable difference in
#'   proportions, negative for "less", or `NA` when not found.
computeEsMinSavi2x2 <- function(
  na,
  nb,
  nBlocksPlan,
  power = 0.8,
  alpha = 0.05,
  alternative = c("twoSided", "less", "greater"),
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
