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

Both effects are **B minus A**, anchored on `thetaA`:
`thetaB = thetaA + propDiff`, or
`thetaB = plogis(qlogis(thetaA) + logOdds)`.
The testing null is equality of proportions (`h0 = 0`). `"greater"`
describes the alternative's direction, not a composite null.

| Effect | `eType` | Alternative | Minimal effect | Interval / sequence | Planning |
|---|---|---|---|---|---|
| `propDiff` | `eBeta` | `twoSided` | None | Yes, defects R1/R4 | No |
| `propDiff` | `grow` | `greater`, `twoSided` | `0 < propDiffMin < 1` | No | Yes |
| `logOdds` | `eGauss` | `twoSided` | None | Yes, defect R3 | No |
| `logOdds` | `grow` | `greater`, `twoSided` | Finite `logOddsMin > 0` | No | No |

`"less"` is supported by the one-block UMP helper only. A one-sided
alternative supplied to eBeta/eGauss is supposed to warn and be ignored;
eBeta currently retains its effect on the first factor ([R8](2x2-review.md#r8)).

## Inputs and result objects

`designSavi2x2(na, nb, ...)` takes planned positive integer block sizes,
optional `propDiffMin` or `logOddsMin` (never both), `alpha = .05`,
`power = NULL`, `nBlocksPlan = NULL`, `h0 = 0`, `alternative`, `eType`,
`betaParameter = NULL`, and `runningIntersection = NULL`.
Simulation settings are `nTheta = 8`, `nSim = 1000`, `nBoot = nSim`,
`nMax = 10000`, `seed = NULL`, `wantSamplePaths = FALSE`, `pb = TRUE`.

The `saviDesign` has `testName = "Two Proportions"`, `testType = "2x2"`,
`h0 = c(propDiff = h0)`, and:

- `nPlan = list(na, nb)`, adding named `nBlocksPlan` when given or planned.
- `esMin`, named by the supplied effect; no `effectMeasure` field.
- `betaParameter = list(betaA1, betaA2, betaB1, betaB2)`: positive success
  and failure shapes, exactly these names; constructor defaults `.18` each.
- `parameter`: named prior summary for printing. Constructor defaults
  `runningIntersection = FALSE`, `relevanceTest = FALSE`; `NULL` arguments
  preserve these defaults. `nMax` is a simulation cap, not a design field.

The two effect-specific entry points are `savi2x2TestStatPropDiff` and
`savi2x2TestStatLogOdds`, currently both taking
`(ya, yb, designObj = NULL, wantCi = TRUE,
wantConfidenceSequence = FALSE, ciValue = NULL)`.
Counts are finite nonnegative integers, one per block in observation order,
and no larger than that block's positive integer group size.

**Unresolved API mismatch:** the approved observed-size contract allows
`na` and `nb` independently as `NULL` (planned size), a broadcast scalar,
or a vector of block length. Current signatures omit these arguments and
read sizes only from the design. Observed sizes were allowed to differ
without warning. Likewise, `designObj = NULL` has no current pilot fallback.
These implementation gaps do not cancel the approved contracts.

Results are `saviTest` objects with cumulative `eValueVec`, its last value
`eValue`, `n = c(na = sum(na), nb = sum(nb), nBlocks = length(ya))`,
`estimate = c(thetaA, thetaB)` of pooled observed proportions, and
`n1Vec = seq_len(nBlocks)` for plotting. PropDiff also stores `betaPrior`,
the posterior shapes after all blocks and hence the prior for a next block.
Do not store per-block sizes, `sumStats`, or `eFactorVec`. A `NULL` constructor
placeholder may be absent because `modifyList()` drops it.

## E-process constructions

All prediction for block `i` uses blocks `1..i-1` only. Shared names are
`logLikelihoodAlternative`, `logLikelihoodNull`, and `logEValueVec`;
the last is the difference of cumulative log likelihoods.

**PropDiff eBeta.** `predictiveThetas2x2` computes independent Beta posterior
means using prior shapes, previous successes, and previous cumulative sizes.
The denominator uses pooled probability
`(na * thetaA + nb * thetaB) / (na + nb)`; sum the two binomial log-density
ratios over blocks. Block 1's plug-in factor is **replaced** by the UMP
conditional factor, as for grow; the posterior still absorbs block 1.

**PropDiff grow.** `logEValueVec2x2PropDiffGrow` learns on
`thetaB = thetaA + propDiffMin`: 1000 interior grid points in `thetaA`'s
feasible interval, with `Beta(betaA1, betaA2)` on its rescaling to `(0,1)`.
Only A's prior shapes are used. Each block uses the grid posterior mean,
the pooled denominator, then updates weights on the log scale.
Block 1's plug-in factor is **replaced** by the corresponding UMP factor;
the posterior still absorbs block 1. Multiplying the two factors was invalid.
Two-sided grow runs separate positive/negative curves and averages their
**cumulative processes**, not their blockwise factors. `earlyStopping = TRUE`
cuts the returned log vector at the replaced, averaged process's first
crossing of `log(1 / alpha)`.

**Conditional logOdds.** Given `ya + yb`, the weighted count is `yb`:
FNCH has odds `exp(logOdds)` on B; the null is hypergeometric.
`logLikelihoodFNCH` returns per-block conditional log densities.
Grow uses a fixed `logOddsMin`, or averages cumulative likelihoods at both
signs; it has no first-block UMP modification.
eGauss mixes cumulative likelihoods under Normal(0,1), normalized on 2000
equally spaced log-odds values in `[-20,20]`. Its first predictive factor
is replaced by UMP `"greater"` even though the mode is two-sided; later
factors use the grid posterior including block 1.

`savi2x2TestStatUmp(ya, yb, na, nb, alpha, alternative)` returns one plain
conditional e-factor, averaging the two one-sided factors for `twoSided`.
`solveUmpLogOdds` solves conditional KL equal to `log(1/alpha)`, searching
up to 100 from `nullLogOdds = 0` in the chosen direction. It returns
`NULL` (factor 1) exactly when the target is unreachable: the KL's supremum,
`-log P0(yb at its extreme)` for that side, is at most `log(1/alpha)`, e.g.
totals `0` or `na + nb`, or `na = nb = 1`. Dependence on the current total is
allowed because the conditional test conditions on it. Numerical FNCH
underflow remains unresolved ([R3](2x2-review.md#r3)).

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
  domain = c(-40,40))` inverts the test's own numerator minus candidate
  conditional log likelihood. It uses `optimize` and `uniroot`, reports
  domain edges when accepted, and warns/returns the domain when empty.
  The sequence repeats this on prefixes, optionally restricting the next
  domain to the previous row; empty rows are `NA`. Runtime is quadratic.

Grow's interval construction is deferred; code nevertheless assigns its
`ciValue` metadata, contrary to the earlier no-`ciValue` agreement.

## Planning on propDiff grow

| Given to `designSavi2x2` | Scenario | Output |
|---|---|---|
| Minimal effect only | None | Design without simulation |
| `propDiffMin`, `power` | `1a` | Worst-baseline planned block count |
| `propDiffMin`, `nBlocksPlan` | `2` | Worst-baseline power |
| `power`, `nBlocksPlan`, no minimal effect | `3` | Minimal detectable `propDiff` |

All three together error. `sampleStoppingTimesSavi2x2` simulates grow on
`nTheta` equally spaced interior baselines of the feasible positive curve;
two-sided adds the negative curve. Both curves matter with unequal group
sizes or priors. “Worst case” means worst **on this finite grid**.
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
worst-baseline identifiers, `breakVector`, and optional paths. Finite-sample
quantile/bootstrap disagreement is deferred ([R6](2x2-review.md#r6)).

`computeEsMinSavi2x2` searches `(0.01,0.9)` with fixed simulation seed and
`uniroot(tol = 1e-5)`; no suitable bracket returns `NA`, which the design
turns into an error. Scenario 3 stores the effect and targets, no simulation
summaries. The discrete objective does not guarantee a passing/minimal
answer ([R5](2x2-review.md#r5)).

LogOdds planning was dropped: as baseline probabilities approach 0 or 1,
a fixed log odds ratio approaches equality on the probability scale.
Unrestricted worst-case planning is therefore controlled by the baseline
grid rather than the effect. `logOddsMin` with planning arguments errors.

## Remaining boundaries

General type/restriction checks, `h0 != 0`, public S3 wrappers/exports,
and grow effect dispatch are unfinished. Do not infer support from an
argument being accepted. Sampler-level `nBoot` and data/e-value-at-cap flags
remain unimplemented. `plot.saviTest` uses block indices and `confSeqMatrix`;
design plotting cannot yet consume per-baseline path lists
([R7](2x2-review.md#r7)). Future work must first agree the affected contract.
