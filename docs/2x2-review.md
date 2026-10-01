# 2×2 review: deferred findings

Reviewed **2026-09-28**, against commit **`e4464b587ea6`**; documented
**2026-09-29**. All findings are **unresolved and deferred by the user**.
This record does not authorize fixes or tests. Agree any changed statistical
contract in the [current design](2x2-design.md) first.

| ID | Priority | Finding | Status |
| --- | --- | --- | --- |
| [R1](#r1) | P1 | Assumed convexity causes confidence-interval undercoverage | Deferred |
| [R3](#r3) | P1 | FNCH density underflow corrupts intervals and e-values | Deferred |
| [R4](#r4) | P1 | A grid can miss the entire accepted confidence set | Deferred |
| [R5](#r5) | P2 | Minimal-effect search can return an effect below target power | Deferred |
| [R6](#r6) | P2 | Planning and bootstrap use different quantiles | Deferred |
| [R7](#r7) | P2 | Design plotting cannot consume the new sample-path structure | Deferred |
| [R9](#r9) | P2 | The block-1 UMP factor moves the worst case off the analytic baseline | Deferred |

P1 affects statistical validity or numerical correctness; P2 affects planning
or API behavior. Expand a finding for reproduction details. Run snippets from
the package root after `devtools::load_all(".")`; each section is independent.
Results below are observations from the review, not a claim that the rest of
the module has been validated. Line numbers refer to the reviewed commit;
function names are the durable reference.

<a id="r1"></a>
<details>
<summary>R1 — The single-interval assumption breaks coverage</summary>

**Impact.** The projected log e-value need not be convex in `propDiff`.
Searching for one minimum and one pair of roots can discard a separate
accepted component. Numerical spot checks in old Decision 13 were not a proof.

```r
prior <- list(betaA1 = 1, betaA2 = 999, betaB1 = 17, betaB2 = 83)
interval <- function(a, b) {
  suppressWarnings(computeConfidenceInterval2x2PropDiff(
    a, b, 181, 18, prior, alpha = .05
  ))
}
interval(27, 14)
delta <- .5
q <- solveRIPr2x2PropDiff(.001, .17, 181, 18, delta)
dbinom(27, 181, .001, log = TRUE) + dbinom(14, 18, .17, log = TRUE) -
  dbinom(27, 181, q, log = TRUE) - dbinom(14, 18, q + delta, log = TRUE)

# Partial exact enumeration: probability of excluded true delta = .5.
tables <- expand.grid(a = 15:40, b = 8:17)
miss <- mapply(function(a, b) {
  ci <- interval(a, b)
  ci[1] > .5 || ci[2] < .5
}, tables$a, tables$b)
sum(dbinom(tables$a, 181, .15) * dbinom(tables$b, 18, .65) * miss)
```

**Observed.** The interval is `[-0.564836, 0.023235]`, but the log e-value
at `.5` is `-9.775783`, below `log(20)`. The partial enumeration gives
`0.6186155`: at least **61.86% noncoverage** for a nominal 95% interval.
This is a lower bound from the listed tables, not an estimate of full
noncoverage from simulation.

**Location.** `R/newsafe2x2Test.R:584–586`,
`computeConfidenceInterval2x2PropDiff()`; related convexity assumptions also
appear in `computeConfidenceSequence2x2PropDiff()`.
**Next design question.** Return every accepted component, or construct and
justify a conservative hull without assuming convexity?

</details>

<a id="r3"></a>
<details>
<summary>R3 — Taking log after the density has underflowed loses information</summary>

**Impact.** `log(BiasedUrn::dFNCHypergeo(...))` can become `-Inf` even when
the mathematical log density is finite. Two-sided averaging can then encounter
`-Inf - -Inf`; interval optimization can return the full search domain.

```r
d <- designSavi2x2(300, 300, eType = "eGauss")
savi2x2TestStat(150, 150, designObj = d)$confSeq
d <- designSavi2x2(1000, 1000, logOddsMin = 5, eType = "grow")
savi2x2TestStat(500, 500, designObj = d, wantCi = FALSE)$eValue

# Stable one-block reference, evaluated directly on the log scale.
logP <- function(eta) {
  k <- 0:300
  w <- lchoose(300, k) + lchoose(300, 300 - k) + eta * k
  2 * lchoose(300, 150) + eta * 150 - max(w) -
    log(sum(exp(w - max(w))))
}
eta <- solveUmpLogOdds(300, 300, 300, .05, "greater")
f <- function(delta) logP(eta) - logP(delta) - log(20)
c(uniroot(f, c(-40, 0))$root, uniroot(f, c(0, 40))$root)
```

**Observed.** eGauss reports `[-40, 40]` with numerical warnings; the stable
reference gives approximately `[-0.566224, 0.566224]`. The grow call gives
`NaN`. The reference is diagnostic evidence, not an implemented replacement.
**Location.** `R/newsafe2x2Test.R:965–968`, `logLikelihoodFNCH()`;
`savi2x2TestStat()` also needs explicit infinite-term handling.
**Next design question.** Which stable log-density method and boundary
conventions should the public calculation use?

</details>

<a id="r4"></a>
<details>
<summary>R4 — No accepted grid points does not establish an empty set</summary>

**Impact.** An accepted interval narrower than the grid spacing can lie between
nodes. The sequence reports an empty set even though its warning describes
the bounds as conservative.

```r
d <- designSavi2x2(1e8, 1e8, eType = "eBeta")
savi2x2TestStat(0, 0, designObj = d,
  wantConfidenceSequence = TRUE)$confSeq
savi2x2TestStat(0, 0, designObj = d)$confSeq
```

**Observed.** The sequence gives `NA` bounds, while root finding gives
approximately `[-0.000183, 0.000183]`. At `propDiff = 0`, the e-value is
exactly one, so that candidate must remain accepted.
**Location.** `R/newsafe2x2Test.R:682–685`,
`computeConfidenceSequence2x2PropDiff()`.
**Next design question.** How should a grid miss be distinguished from genuine
emptiness, with a conservative fallback that also respects R1?

</details>

<a id="r5"></a>
<details>
<summary>R5 — The reported minimal effect may not meet target power</summary>

**Impact.** `uniroot()` is applied to a simulated step function. Its returned
point need not pass the power target; a zero plateau also does not identify
the smallest passing effect.

```r
d <- designSavi2x2(5, 5, nBlocksPlan = 5, power = .51,
  eType = "grow", nTheta = 1, nSim = 10, seed = 1, pb = FALSE)
d$esMin
computePowerSavi2x2(unname(d$esMin), 5, 5, nBlocks = 5,
  nTheta = 1, nSim = 10, seed = 1, pb = FALSE)$power
```

**Observed.** Reported effect `0.4581119`; power using the same simulation
settings is `.50`, below `.51`.
**Location.** `R/newsafe2x2Test.R:1371–1373`, `computeEsMinSavi2x2()`.
**Next design question.** Define the target as a smallest verified passing
effect to a stated resolution, and decide what uncertainty to report.

</details>

<a id="r6"></a>
<details>
<summary>R6 — Bootstrap uncertainty targets a different quantile</summary>

**Impact.** The sampler uses the type-1 quantile, while the general bootstrap
uses R's default type 7. With censored stopping times, the original statistic
and its bootstrap counterpart can disagree even about finiteness.

```r
d <- designSavi2x2(5, 5, propDiffMin = .7, power = .8, eType = "grow",
  nTheta = 1, nSim = 10, nBoot = 20, nMax = 5, seed = 4, pb = FALSE)
d$nPlan$nBlocksPlan
d$bootObjNBlocksPlan$t0
d$bootObjNBlocksPlan$bootSe
```

**Observed.** Planned blocks `4`, bootstrap original statistic `Inf`,
bootstrap standard error `NaN`.
**Location.** `R/newsafe2x2Test.R:1289–1292`, `computeNPlanSavi2x2()`;
`R/designHelpers.R:283`, `computeBootObj()`.
**Next design question.** Align the quantile definition and define meaningful
uncertainty when bootstrap resamples remain censored, within the 2×2 scope.

</details>

<a id="r7"></a>
<details>
<summary>R7 — Design plotting expects a matrix but receives a list</summary>

**Impact.** The design now retains one sample-path matrix per baseline;
`plot.saviDesign()` still assumes one matrix.

```r
d <- designSavi2x2(5, 5, propDiffMin = .5, nBlocksPlan = 3,
  eType = "grow", nTheta = 1, nSim = 10, nBoot = 10,
  pb = FALSE, wantSamplePaths = TRUE)
plot(d)
```

**Observed.** `invalid 'length' argument` at `integer(mIter)` because
`dim(samplePaths)[1]` is `NULL` for the list.
**Location.** `R/safeS3Methods.R:613–616`, `plot.saviDesign()`.
**Next design question.** Show the worst baseline, a selected baseline, or
several? Reconcile planned block-count fields when agreeing plotting behavior.

</details>


<a id="r9"></a>
<details>
<summary>R9 — The block-1 UMP factor moves the worst case off the analytic baseline</summary>

**Impact.** Planning simulates one baseline per curve, the minimiser of the
asymptotic growth rate (`solveWorstCaseTheta2x2PropDiff`). The test replaces
block 1 by the UMP conditional e-factor, whose expected log value under the
alternative is negative and depends on the baseline: about −1.8 nats for
10:2 and −2.2 nats for 5:5 at `d = 0.1` against a threshold of `log 20 ≈ 3`.
For very small unequal groups that handicap peaks away from the analytic
baseline, so the single-baseline plan is short.

**Evidence (2026-09-29/10-01, 300–600 paths, `alpha = 0.05`, B-minus-A
convention of that date, mirrored here).** With the UMP block, the 80%
stopping-time quantile at the analytic baseline versus the worst of an
11-point grid: 10:2 and 2:10 at `d = 0.1` are 30% short (158 vs 226; 149 vs
214), 20:5 is 16% short (109 vs 129); 3:7, all balanced designs and every
design at `d = 0.3` agree within noise. The fine-grid worst case sits at
`thetaLow ≈ 0.31` for 10:2 (analytic 0.46). Without the UMP replacement the
same designs agree within noise everywhere, with the same stopping-time
quantile in both group orders. The prior (0.18, Jeffreys, uniform, skewed)
has no visible effect. The exact adjustment that reproduced the shifted peak
was `argmax (log(1/alpha) − E[log UMP_1]) / R(theta)`, with the expectation
an exact sum over block-1 tables.

**Location.** `sampleStoppingTimesSavi2x2()` and
`solveWorstCaseTheta2x2PropDiff()`; the UMP replacement in
`savi2x2TestStat()`.

**Next design question.** Add the block-1 term to the baseline objective, fall
back to a grid when `min(na, nb)` is small, or revisit the UMP replacement,
whose expected cost is a large fraction of the threshold at small effects.

</details>

## Validation boundary

The prior review used direct R reproductions. Its prescribed command,
`Rscript -e 'testthat::test_local(".", filter = "2x2")'`, reported
**“No test files found.”** No fixes or tests were added for this documentation
change. Broader validation and any regression tests remain separate work to
agree when the user resumes the findings.
