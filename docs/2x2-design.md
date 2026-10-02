# 2×2: current design

This is the current agreement for [the new module](../R/newsafe2x2Test.R),
condensed from Decisions 1–43 at `e4464b587ea6:AGENTS.md`.
Read [working rules](../AGENTS.md) and
[deferred review findings](2x2-review.md) alongside it. Update the relevant
section when a change is agreed; Git preserves superseded decisions.

**Status:** several agreed constructions have known defects. The user has
deferred fixes. Describing them here does not establish validity or approve
a replacement; the issue register separates evidence from proposed remedies.
No fix is currently authorized; the next step is the user's choice of a
finding or design question, whose replacement contract is agreed here first.

## Supported modes

Both effects are **A minus B**, anchored on `thetaB`:
`thetaA = thetaB + propDiff`, or
`thetaA = plogis(qlogis(thetaB) + logOdds)`.
This follows `t.test(x, y)`, where `"greater"` means `x - y > 0`: here
`"greater"` means group A has the larger proportion or odds.
The testing null is equality of proportions (`h0 = 0`). `"greater"`
describes the alternative's direction, not a composite null.

| Effect | `eType` | Alternative | Minimal effect | Interval / sequence | Planning |
|---|---|---|---|---|---|
| `propDiff` | `eBeta` | `twoSided`, `greater`, `less` | None | Yes, defects R1/R4 | No |
| `propDiff` | `grow` | `twoSided`, `greater`, `less` | signed `propDiffMin`, `0 < |propDiffMin| < 1` | No | Yes |
| `logOdds` | `eGauss` | `twoSided`, `greater`, `less` | None | Yes | No |
| `logOdds` | `grow` | `twoSided`, `greater`, `less` | signed finite `logOddsMin != 0` | No | No |

**Sign rule for grow.** The minimal effect is a signed value on the
A-minus-B scale of its effect. `greater` tests a positive value and learns
on `thetaA` above `thetaB`; `less` tests a negative value and learns on
`thetaA` below `thetaB`. The constructor passes `propDiffMin` or
`logOddsMin` through `checkAndReturnEsMinParameterSide`, as the z and t
designs do: a value whose sign contradicts a one-sided alternative is
flipped with a warning naming the effect, and `twoSided` takes the
magnitude and runs both signs of it, so the stored `esMin` always shows the
sign of the direction tested. Zero is an error (no restriction), as is a
value outside `(-1, 1)` for propDiff or a non-finite one for logOdds. The
legacy `safe2x2Test.R` inferred the alternative from the sign of `delta` and
measured B minus A; here the alternative sets the sign and the effect is A
minus B. `less` on `(ya, yb, na, nb)` with effect `d` equals `greater` on
`(yb, ya, nb, na)` with `-d`: exactly on logOdds, and on propDiff because
both curves rescale the same `rho`. The direction is validated once, in the
constructor; the test, helper and sampler code apply it through group A
only, so `greater`, a positive `signs` entry, a positive `logOdds`, and
odds `exp(logOdds)` on A all mean the same thing and no formula negates the
effect except the `twoSided` mirror.
The constructor takes arguments as given and checks only cheap
inconsistencies: `alpha` and `power` lie in `(0, 1)`, `na` and `nb` have
equal length, at most one minimal
effect is given, `propDiffMin` is nonzero in `(-1, 1)` and `logOddsMin` is
nonzero; `nSim`, `nBoot`, `nMax`, the priors and `runningIntersection` are
not checked. The flow is linear: fix `nBlocksPlan` from per-block sizes,
sign the minimal effect, then apply the eType's rules. eBeta and eGauss stop
on a minimal effect, drop a `power` with a warning, keep `nBlocksPlan` as
the planned count (stored and printed, read only by scenarios 2 and 3), and
read one prior each, `betaParameter` for eBeta and `gaussParameter` for
eGauss, stored as given (`NULL` for the default, which `savi2x2TestStat`
fills in on the data's block sizes); the prior the eType does not
read is ignored and not stored. grow reads no prior and ignores a supplied
one; it plans as in the planning table below, and any other combination of
`power` and `nBlocksPlan` is an error, except that `logOddsMin` drops a
`power` with a warning. `betaParameter` is
`list(betaA1, betaA2, betaB1, betaB2)`, exactly these names, each a single
finite positive number; `NULL` means `1 / (2 * na)` and `1 / (2 * nb)`,
taking block 1's sizes with a warning when sizes vary by block; a given
`betaParameter` opts out of the UMP first block (below). All three
eTypes accept every `alternative`; eBeta and eGauss use a one-sided one in
the UMP first block, and eGauss also restricts its prior to that side.
`gaussParameter` is `list(mean, sd)` with finite `mean` inside the grid
`(-20, 20)` and finite `sd > 0`; `NULL` means `list(mean = 0, sd = 1)`.

## Inputs and result objects

`designSavi2x2(na, nb, ...)` takes planned finite positive integer block
sizes as numeric vectors of equal length: one value each for a constant
size, or one per block, which fixes `nBlocksPlan` at their length (a
supplied `nBlocksPlan` is taken as given; lists are rejected). It also takes
optional `propDiffMin` or `logOddsMin` (never both), `alpha = .05`,
`power = NULL`, `nBlocksPlan = NULL`, `h0 = 0`, `alternative`, `eType`,
`betaParameter = NULL`, `gaussParameter = NULL`, and `runningIntersection = NULL`.
Simulation settings are `nSim = 1000`, `nBoot = nSim`,
`nMax = 10000`, `seed = NULL`, `wantSamplePaths = FALSE`, `pb = TRUE`.

The `saviDesign` has `testName = "Two Proportions"`, `testType = "2x2"`,
`h0 = c(propDiff = h0)`, and:

- `nPlan = list(na, nb)`, adding named `nBlocksPlan` when given or planned;
  the print methods show per-block size vectors as `mean na`, `mean nb`.
- `esMin`, named by the supplied effect; no `effectMeasure` field.
- `betaParameter = list(betaA1, betaA2, betaB1, betaB2)`, eBeta only:
  positive success and failure shapes, exactly these names, as given;
  `NULL` (the default) means `1 / (2 na)`, `1 / (2 nb)` on block 1's sizes
  (a warning when sizes vary).
- `gaussParameter = list(mean, sd)`, eGauss only: the Normal prior on
  `logOdds`, exactly these names, as given; `NULL` (the default) means
  `list(mean = 0, sd = 1)`.
  The prior the eType does not read is absent (`NULL`).
- `parameter`, eBeta and eGauss only, for printing: the prior as one named
  string formatted in the constructor, names and values comma-separated
  (`betaA1, betaA2, betaB1, betaB2 = 0.05, 0.05, 0.05, 0.05`), or `"default"`
  named `betaParameter` or `gaussParameter` when the prior is `NULL`. grow
  stores none, since `esMin` already shows its minimal effect; the shared
  print methods skip a missing 2x2 `parameter`.
  Constructor defaults
  `runningIntersection = FALSE`, `relevanceTest = FALSE`; `NULL` arguments
  preserve these defaults. `nMax` is a simulation cap, not a design field.

The single entry point is `savi2x2TestStat(ya, yb, designObj = NULL,
wantCi = FALSE, wantConfidenceSequence = FALSE, ciValue = NULL)`. It reads
the effect from the design: `propDiff` for eBeta, `logOdds` for eGauss, and
for grow the name of `esMin`. The shared checks, the UMP replacement of
block 1, and the result fill are written once; each of the four e-processes
is a helper returning the plain cumulative log e-process, and the
confidence code dispatches on the effect.
Counts are finite nonnegative integers, one per block in observation order,
and no larger than that block's positive integer group size.
`designSavi2x2` and `savi2x2TestStat` are exported; every helper has a
roxygen page but stays internal.

**Unresolved API mismatch:** the approved observed-size contract allows
`na` and `nb` independently as `NULL` (planned size), a broadcast scalar,
or a vector of block length. Current signatures omit these arguments and
read sizes only from the design. Observed sizes were allowed to differ
without warning. Likewise, `designObj = NULL` has no current pilot fallback.
These implementation gaps do not cancel the approved contracts.

Results are `saviTest` objects with cumulative `eValueVec`, its last value
`eValue`, `n = c(na = sum(na), nb = sum(nb), nBlocks = length(ya))`,
`estimate = c(thetaA, thetaB)` of pooled observed proportions, and
`n1Vec = seq_len(nBlocks)` for plotting. eBeta also stores `betaParameter`,
the posterior shapes after all blocks and hence the prior for a next block.
Do not store per-block sizes, `sumStats`, or `eFactorVec`. A `NULL` constructor
placeholder may be absent because `modifyList()` drops it.

## E-process constructions

All prediction for block `i` uses blocks `1..i-1` only. Shared names are
`logLikelihoodAlternative`, `logLikelihoodNull`, and `logEValueVec`;
the last is the difference of cumulative log likelihoods.

**PropDiff eBeta.** `logEValueVec2x2PropDiffEBeta` returns the plain
process. `predictiveThetas2x2` computes independent Beta posterior
means using prior shapes, previous successes, and previous cumulative sizes.
The denominator uses pooled probability
`(na * thetaA + nb * thetaB) / (na + nb)`; sum the two binomial log-density
ratios over blocks. With the default prior, block 1's plug-in factor is
**replaced** by the UMP conditional factor, as for grow; the posterior still
absorbs block 1. A given `betaParameter` keeps block 1's plug-in factor,
with a message: a user prior already defines block 1's numerator. This
opt-out is eBeta only; grow and eGauss always replace.

**PropDiff grow.** `logEValueVec2x2PropDiffGrow` learns on
`thetaA = thetaB + propDiffMin` with the signed value: 1000 interior grid
points in `thetaA`'s feasible interval, `(d, 1)` for `d > 0` and
`(0, 1 + d)` for `d < 0`, with a uniform prior, `Beta(1, 1)`, on its rescaling
to `(0,1)`; grow reads no prior, and neither the helper nor the sampler and
planners take a prior argument. Each
block uses the grid posterior mean, the pooled denominator, then updates
weights on the log scale.
Two-sided grow runs separate positive/negative curves and averages their
**cumulative processes**, not their blockwise factors. The helper returns
the plain plug-in process; `savi2x2TestStat` then **replaces**
block 1 of that (averaged) process by `savi2x2TestStatUmp` at the test's
`alternative`, as for eBeta; the posterior still absorbs block 1.
Multiplying the two factors was invalid. `earlyStopping = TRUE` cuts the
helper's returned log vector at the plain averaged process's first crossing
of `log(1 / alpha)`.

**Conditional logOdds.** Given `ya + yb`, the weighted count is `ya`:
FNCH has odds `exp(logOdds)` on A; the null is hypergeometric.
`logLikelihoodFNCH` returns per-block conditional log densities, evaluated
on the log scale (binomial coefficients plus `logOdds` times the feasible
count, normalized by log-sum-exp), so they are finite at every finite
`logOdds`.
`logEValueVec2x2LogOddsGrow` uses the magnitude of `logOddsMin` with the sign
of `alternative`, or for `twoSided` averages the cumulative e-processes at
both signs.
`logEValueVec2x2LogOddsEGauss` returns the plain eGauss process: cumulative
likelihoods mixed under the design's `gaussParameter` prior Normal(mean, sd),
restricted to the side of `alternative` (`greater`: log-odds above 0, `less`:
below 0, `twoSided`: both) and normalized on 2000 equally spaced log-odds
values in `[-20,20]`; later factors use the grid posterior including
block 1. For both eTypes `savi2x2TestStat`
**replaces** block 1 of the cumulative process by `savi2x2TestStatUmp` at
the design's `alternative`, as for grow on propDiff.

`savi2x2TestStatUmp(ya, yb, na, nb, alpha, alternative)` returns one plain
conditional e-factor, averaging the two one-sided factors for `twoSided`.
`solveUmpLogOdds` solves conditional KL equal to `log(1/alpha)`, searching
up to 100 from `logOddsNull = 0` in the chosen direction. It returns
`NULL` exactly when the target is unreachable: the KL's supremum,
`-log P0(yb at its extreme)` for that side, is at most `log(1/alpha)`, e.g.
totals `0` or `na + nb`, or `na = nb = 1`. Dependence on the current total is
allowed because the conditional test conditions on it. On `NULL` for any
side, `savi2x2TestStatUmp` returns `NULL` and `savi2x2TestStat` keeps block 1
of the plain process instead of replacing it, with a warning.

## Confidence intervals and sequences

Default confidence level is `1 - alpha`; `ciValue` can override it.
Only eBeta/eGauss supply intervals. `wantConfidenceSequence` takes precedence
over `wantCi`, produces `confSeqMatrix` with `nBlocks × 2` columns
`lowerBound`, `upperBound`, and stores its last row as `confSeq`.
Otherwise `wantCi` stores one named two-element `confSeq` on all data.
Running intersection is optional; an empty row stays empty thereafter.

- **PropDiff interval:** `computeConfidenceInterval2x2PropDiff` takes the
  observed vectors, `betaParameter`, `alpha`, `domain = c(-1,1)`.
  It inverts the predictable Beta numerator against per-block reverse
  information projections onto each candidate difference, excluding the
  test's UMP first-block factor. Current root finding uses one minimum and two roots;
  failure to find an interval warns and returns `[-1,1]`.
  The asserted convexity/single-interval justification is false
  ([R1](2x2-review.md#r1)); its replacement needs agreement.
- **PropDiff sequence:** `computeConfidenceSequence2x2PropDiff` advances
  candidate log e-processes once per block on 2000 interior points.
  It uses outward secant bounds and optional intersection, with `NA` empty
  rows. The claimed conservative guarantee is unresolved: the convexity
  premise fails, and a grid can miss an entire accepted set
  ([R4](2x2-review.md#r4)).
- **LogOdds:** `computeConfidenceInterval2x2LogOdds(..., logPTotal, alpha,
  domain = c(-40,40))` inverts the test's own numerator, the eGauss process
  at the design's `gaussParameter` and `alternative`, minus candidate
  conditional log likelihood. It uses `optimize` and `uniroot`, reports
  domain edges when accepted, and warns/returns the domain when empty.
  `computeConfidenceSequence2x2LogOdds(..., logNumerator, alpha,
  runningIntersection)` repeats this on prefixes, optionally restricting
  the next domain to the previous row; empty rows are `NA`. The test
  computes `logNumerator` once, the plain cumulative log e-process kept
  before the UMP replacement of block 1 plus the cumulative hypergeometric
  log likelihood, and passes the vector to the sequence and its last
  element to the interval. Runtime is quadratic in the root finding.

Grow's interval construction is deferred; code nevertheless assigns its
`ciValue` metadata, contrary to the earlier no-`ciValue` agreement.

## Planning on propDiff grow

| Given to `designSavi2x2` | Scenario | Output |
|---|---|---|
| Minimal effect only | None | Design without simulation |
| `propDiffMin`, `power` | `1a` | Worst-baseline planned block count |
| `propDiffMin`, `nBlocksPlan` | `2` | Worst-baseline power |
| `power`, `nBlocksPlan`, no minimal effect | `3` | Minimal detectable `propDiff` |

All three together error. Per-block size vectors supply `nBlocksPlan`, so
they run scenarios 2 and 3 on exactly those sizes and make `1a` an error. `sampleStoppingTimesSavi2x2` simulates grow at
one baseline per curve the test runs: the worst case from
`solveWorstCaseTheta2x2PropDiff`, the root of
`nLow logit(theta) + nHigh logit(theta + d) = n logit(theta0)` with `theta`
the lower proportion, `theta0` its pooled null projection and `d = |propDiffMin|`,
summed over blocks for per-block sizes; equal sizes give exactly
`(1 - d) / 2`, and in general `(1 - d) / 2 + d (nHigh - nLow) / (6 n) + O(d^3)`.
This minimises the asymptotic growth rate `R(theta) = nLow KL(theta || theta0)
+ nHigh KL(theta + d || theta0)`; simulation of the plain plug-in process
(no UMP block) on fine `thetaA` grids found its stopping-time mean and
80/90% quantiles at this root within noise of the grid maximum for size
ratios up to 100:1 and `d` up to 0.5, with a flat plateau of width about
0.2 to 0.3 around it. Two-sided runs both curves, which differ with unequal
sizes. The rate ignores the block-1 UMP factor; with it, the
simulated worst case moves for very small unequal groups at small `d`
([R9](2x2-review.md#r9)). Agreed next step: double-check against an analytic
stopping-time approximation, `E[tau] ~ log(1 / alpha) n / (2 d^2 na nb)` plus
the plug-in learning cost.
Paths stop at `E >= 1/alpha`; noncrossing times are `Inf`, not `nMax`.
`seed = NULL` means 2026. Outputs include baseline probabilities and
baseline-by-path matrices `stoppingTimes`, `breakVector` (0 crossed,
1 capped), `eValuesStopped`; optional `samplePaths` is a list of sparse
matrices, with crossed values repeated through the cap.

`computePowerSavi2x2` takes the minimum crossing fraction.
`computeNPlanSavi2x2` takes the maximum type-1 stopping-time quantile,
warning when infinite. Bootstrap summaries use the selected worst baseline;
the planned-count and capped-mean bootstraps are omitted when `nPlan = Inf`.
Scenarios 1a/2 store their bootstrap objects, relevant `*TwoSe` fields,
the worst baseline (`worstCaseThetaA`, `worstCaseThetaB`), `breakVector`, and optional paths. Finite-sample
quantile/bootstrap disagreement is deferred ([R6](2x2-review.md#r6)).

`computeEsMinSavi2x2` searches magnitudes in `(0.01,0.9)` with fixed
simulation seed and `uniroot(tol = 1e-5)`, returning the negative value for
`less`; no suitable bracket returns `NA`, which the design turns into an error. Scenario 3 stores the effect and targets, no simulation
summaries. The discrete objective does not guarantee a passing/minimal
answer ([R5](2x2-review.md#r5)).

LogOdds planning was dropped: as baseline probabilities approach 0 or 1,
a fixed log odds ratio approaches equality on the probability scale.
Unrestricted worst-case planning is therefore controlled by the baseline
grid rather than the effect. `logOddsMin` with planning arguments errors.

## Remaining boundaries

General type/restriction checks, `h0 != 0`, public S3 wrappers,
and grow effect dispatch are unfinished. Do not infer support from an
argument being accepted. Sampler-level `nBoot` and data/e-value-at-cap flags
remain unimplemented. `plot.saviTest` uses block indices and `confSeqMatrix`;
design plotting cannot yet consume per-baseline path lists
([R7](2x2-review.md#r7)). Future work must first agree the affected contract.
