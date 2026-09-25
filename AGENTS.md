# Agent instructions

## Objective

Rewrite the safe anytime-valid 2x2 (two-proportion) test as a small, direct
implementation. Priorities, in order: statistical validity, explicit
conventions, numerically stable code, focused tests, then everything else.

Behaviour of the older code is not a requirement. Keep it only when it is
known to be intentional and correct. The `cond` branch holds the previous
attempt and may be consulted with `git show cond:<path>`, never merged.

## Scope

In scope: `R/newsafe2x2Test.R` (all new code goes here; the legacy
`R/safe2x2Test.R` is left untouched until the rewrite replaces it), its
tests under `tests/testthat/`, its roxygen and generated `man/` pages, and
the smallest 2x2-related edits to `NAMESPACE`, `DESCRIPTION`,
`R/safeS3Methods.R` and `R/deprecate.R`.

Out of scope: the t-test, z-test, log-rank test, the general design
framework, unrelated S3 methods, vignettes unless asked, package-wide style.

## Workflow

Each feature is agreed before it is coded:

1. the request is described in plain words;
2. it is reviewed and questions are settled;
3. the agreed contract is recorded under **Decisions** below;
4. only then is the R code written (tests only when asked).

Do not implement anything not yet recorded under Decisions. If code and a
decision disagree, the decision is the reference; change the code or reopen
the decision, never silently diverge.

## Conventions

- Groups are `A` and `B`. `ya`, `yb` are success counts, `na`, `nb` group
  sizes. Counts are integer-valued, non-negative and at most the group size.
- Two effect measures, named exactly `propDiff` (difference of
  proportions) and `logOdds` (log odds ratio). Both are signed **B minus A**
  and both are anchored on `thetaA`: group B is always derived from group A
  and the effect, never the other way round.
  - `thetaB = thetaA + propDiff`, so `propDiff = thetaB - thetaA`.
  - `thetaB = plogis(qlogis(thetaA) + logOdds)`, so
    `logOdds = logit(thetaB) - logit(thetaA)`.
  The legacy spellings `difference`, `linearDifference` and `logOddsRatio`
  are not used anywhere in new code.
- `alternative` takes `c("twoSided", "greater", "less")`, in that order.
- Distinguish a blockwise e-factor from the cumulative e-process in names,
  documentation and tests. Step `i` uses only data through step `i`.
- Work on the log scale wherever underflow or overflow is plausible. Handle
  zero-probability cases explicitly; never let `NaN` or `Inf - Inf` decide.
- Mirror the layout and naming of `R/tTest.R`: a pure `...Stat` function,
  then the S3 generic with `.default` and `.formula`, then the `savi.*.test`
  alias, then design and sampling functions.

## Testing

Write tests only when asked. Do not test basic scaffolding such as object
fields, default values or name assignment. When asked, test against the
public API, not the internals. Prefer small exact tests: a hand-computable
table, one block versus the first step of a multi-block call, `log = TRUE`
versus the log of the plain output, a group swap with the mirrored
alternative, boundary and impossible inputs, and seeded simulation output. No Cartesian test grids.
Never loosen a tolerance before the discrepancy is understood.

Run 2x2 tests with `Rscript -e 'testthat::test_local(".", filter = "2x2")'`.
Regenerate docs with `LANG=en_US.UTF-8 Rscript -e 'devtools::document()'`;
the C locale corrupts unrelated `.Rd` files.

## Usage script

`~/Downloads/local-run.R` (outside the repository) mimics how users will
call the 2x2 design and test. Keep it working when the public API changes,
and run it from the package root after `devtools::load_all()` as a quick
end-to-end check.

## R style

2-space indent, `<-`, spaces around operators, `camelCase`, `match.arg()`
for enumerated arguments, guard clauses over nesting, roxygen2 Markdown on
exported functions with effect direction and return semantics stated.
Comments explain statistical intent, not syntax.

## Decisions

Agreed contracts, one entry per feature, newest last. Each entry states the
function, its inputs, its output, and the convention it relies on.

### 1. Sample-size vocabulary

`nSim` is the number of simulated paths, `nBoot` the number of bootstrap
resamples. `nPlan` is `list(na, nb)`, the planned per-group block size
(reopened by Decision 12: block count is observed data, not planned, so
`nPlan` no longer carries `nBlocks`; `plot.saviDesign` reading `nPlan[1]`
as `nBlocks` is stale and needs revisiting). The 2x2 design has no `nMax`,
`highN` or `nPlanBatch` (left `NULL`), and no `thetaA` scenario.

### 2. Design and test object fields

Two effect measures, `propDiff` and `logOdds`, both first-class; there is
still no `effectMeasure` field (reopened by Decision 12: `eType` alone
picks the effect and the e-variable construction: `"eBeta"` or `"grow"`
for `propDiff`, `"eGauss"` or `"grow"` for `logOdds`; `"grow"`'s dispatch
between the two is not yet designed). `designSavi2x2(na, nb, propDiffMin,
logOddsMin, alpha, power, h0, alternative, eType, betaParameter,
runningIntersection)` returns a
`saviDesign` with `testName = "Two Proportions"`, `testType = "2x2"`,
`h0 = c(propDiff = h0)`, and:

- `betaParameter`: `list(betaA1, betaA2, betaB1, betaB2)`, the
  success and failure shapes of the Beta priors on `thetaA` and `thetaB`.
  The default, all four `0.18`, lives in `constructSaviDesignObj("Two
  Proportions")`; the design function replaces it only when the argument is
  not `NULL`, and errors on other names.
- `parameter`: a length-one named string summarising the prior, for
  printing.
- `runningIntersection`: `TRUE` or `FALSE`, default `TRUE` in
  `constructSaviDesignObj("Two Proportions")`. The test reads it for the
  confidence sequence (Decision 4); `plot()` and `print()` read it as for
  the other tests.

Reserved for when the conditional e-variable returns: its prior on `logOdds`
defaults to mean `0` and sd `1`.

`constructSaviTestObj("Two Proportions")` gets no `sumStats` or
`eFactorVec`. Note that `modifyList()` drops `NULL` placeholders, so a
field declared `NULL` in a constructor is absent until a function sets it.

### 3. Test function, unrestricted two-sided case

Superseded by the `propDiff`/`logOdds` split of Decision 12
(`savi2x2TestStatPropDiff`, `savi2x2TestStatLogOdds`); the numerator/
denominator construction below is unchanged for `propDiff`.

`savi2x2Test(ya, yb, designObj = NULL)`: `ya`, `yb` are per-block success
counts in observation order; `na`, `nb`, `alternative`, `h0` and the prior
all come from `designObj` (`NULL` gives a pilot `designSavi2x2()` with a
warning). No `ciValue` or confidence sequence yet.

- Numerator for block `i`: Beta posterior means of `thetaA`, `thetaB` given
  blocks `1..i-1` only. Denominator: the common
  `(na * thetaA + nb * thetaB) / (na + nb)`, the projection onto
  `thetaA = thetaB`.
- `eValueVec` is the cumulative e-process, `eValue` its last element.
- `n = c(na, nb, nBlocks)` in totals; `estimate` holds both observed
  proportions and `propDiff` (B minus A); `betaParameter` on the result is
  the Beta posterior after the last block (same field name as the design's
  prior, but holding the posterior; flagged as ambiguous, not yet resolved).
- Errors, until each is designed: non-`NULL` `esMin`, `alternative` other
  than `"twoSided"`, `h0 != 0`. The design rejects non-positive Beta shapes.

### 4. Confidence sequence for propDiff

`computeConfidenceInterval2x2PropDiff(ya, yb, na, nb, betaParameter,
alpha, precision = 100)` inverts the test on `precision` equally spaced
candidates strictly inside `(-1, 1)`. Each candidate `propDiff` is a point
null with its own e-process: numerator the predictable Beta posterior mean
(shared helper `predictiveThetas2x2`), denominator its reverse information
projection onto `thetaB - thetaA = propDiff`, found by `uniroot` on the KL
derivative (`solveRIPr2x2PropDiff`). `runningIntersection = TRUE` (the
default): a candidate leaves for good at `1/alpha`, and its e-process is
not advanced further. `FALSE`: every candidate's e-process is advanced each
block and block `i` keeps the candidates whose current e-value is below
`1/alpha`, so the sets need not be nested. Each run of consecutive non-rejected
candidates is one interval, and the confidence set is the union of the
intervals; min and max over the whole set are not taken, since that would
fill holes. Returns a two-column matrix, `lowerBound` and `upperBound`
(reopened by Decision 12: block is the rowname, not a data column, since
`ncol` stays 2): one row per block without holes, as for the z-test,
several rows for a block with holes, none for a block with everything
rejected.

`savi2x2TestStatPropDiff(..., wantCi = TRUE)` passes the design's
`runningIntersection` through and stores, as the other tests do, a two-column
`nBlocks x 2` `confSeqMatrix` (`lowerBound`, `upperBound`; other code reads it
by position): row `i` is the outermost bounds of block `i`'s union, `NA` for a
fully rejected block. This hull contains the union, so coverage holds, and
the rows stay nested under the running intersection. The exact union is kept only for the last block, as
`confSeq` (a `k x 2` matrix of `lowerBound`, `upperBound`; `k = 1` without
holes), with `ciValue = 1 - alpha` (no separate `ciValue` argument).
`logLikelihoodRatioMultiBern` (renamed from `savi2x2TestStat`, Decision 12)
returns the cumulative e-process, on the log scale when `log = TRUE`.

### 5. Plotting with plot.saviTest

`plot.saviTest` works on 2x2 results unchanged in its general logic:

- `savi2x2Test` sets `n1Vec = seq_len(nBlocks)`, the block index, as the
  x-axis (legacy name the plot reads).
- `constructSaviDesignObj("Two Proportions")` defaults `relevanceTest =
  FALSE`; `NULL` makes the plot's `&&` error.
- The label switches map `"Two Proportions"` to `"Number of blocks"` (x) and
  `"propDiff"` (confidence-sequence y).
- For `wantConfSeqPlot = TRUE` the plot reads `confSeqMatrix` like any
  other test's; no 2x2 branch. The hull per block (Decision 4) fills holes
  in the picture; `confSeq` stays the exact union.

### 6. Restricted alternative on propDiff

Only `propDiff` restrictions; `logOdds` is out of scope. `designSavi2x2`
accepts `propDiffMin` as `NULL` or one number strictly inside `(0, 1)`,
stored as `esMin`. Allowed combinations; everything else errors:

- `propDiffMin = NULL`, `"twoSided"`: the unrestricted test of Decision 3.
- `propDiffMin > 0`, `"greater"`: the numerator is restricted to the
  curve `thetaB - thetaA = propDiffMin`.
- `propDiffMin > 0`, `"twoSided"`: the e-process is the average of the
  two cumulative e-processes restricted at `+propDiffMin` and
  `-propDiffMin` (averaged as processes, not per block).
- `"less"`, or `"greater"` without `propDiffMin`, is not designed yet.

The restricted numerator builds on `learnPredictiveThetas` from `cond`:
`predictiveThetas2x2PropDiff(..., propDiff, nWeight = 1000)`, next to
`predictiveThetas2x2(...)` for the Beta posterior means of Decision 3;
`computeEValueVecPropDiff(..., propDiff = NULL)` (renamed from
`logEProcess2x2PlugIn`, Decision 12) picks between them and returns
`eValueVec` (length `nBlocks`) directly, not a list. The free
coordinate `rho` is `thetaA` rescaled to its feasible interval
`(max(0, -propDiff), min(1, 1 - propDiff))`, on `nWeight` equally spaced
grid points strictly inside `(0, 1)`, with the prior `Beta(betaA1, betaA2)`;
`betaB*` are not used. Block `i` uses the posterior mean of `thetaA` given blocks `1..i-1`,
and `thetaB = thetaA + propDiff`. The weights are updated on the log scale.
The denominator stays the pooled projection onto `thetaA = thetaB`, so the
null is always the point `thetaA = thetaB`; `"greater"` names the direction
of the alternative, not a composite null. The confidence sequence keeps the
unrestricted numerator. `nWeight` is not a design field.

### 7. Plug-in conditional e-factor on logOdds

`logLikelihoodRatioFNCH(ya, yb, na, nb, logOdds, log = FALSE)`
returns the conditional e-factor of **one** block (scalar counts), not a
cumulative e-process; `savi2x2CondStat(ya, yb, na, nb, logOdds, ...)` applies
it blockwise to vectors and returns the per-block e-factors (its `eType` and
`alternative` dispatch is not designed yet). It conditions on the block's total `ya + yb`: under
the null `yb` is hypergeometric, under `logOdds` (B minus A, anchored on
`thetaA`) it is Fisher's noncentral hypergeometric with odds `exp(logOdds)` on
group B, and the log e-factor is `log(dFNCHypergeo(yb, nb, na, ya + yb,
exp(logOdds)))` minus the same at `logOdds = 0`, on the log scale when `log =
TRUE`. Weighting `yb` is what gives the B-minus-A sign; an earlier version
weighted `ya` and so measured A minus B.
`logOdds` is one finite number
supplied by the caller (plug-in, e.g. GROW or UMP); a prior on `logOdds` is
reserved for later.

`fnchLogPartition(na, nb, totalSuccesses, logOdds)` is the log of
`sum_k choose(na, k) choose(nb, totalSuccesses - k) exp(logOdds * k)` over the
feasible `k`, computed by a max-shifted log-sum-exp; at `logOdds = 0` it
returns `lchoose(na + nb, totalSuccesses)` exactly.

### 8. UMP plug-in conditional e-factor on logOdds

`solveUmpLogOdds(na, nb, totalSuccesses, alpha, alternative = c("greater",
"less"), nullLogOdds = 0, searchBound = 100)` finds the UMP plug-in `logOdds` for
**one** block by `uniroot` on `(nullLogOdds, nullLogOdds + 100)` for `"greater"`
and `(nullLogOdds - 100, nullLogOdds)` for `"less"`: the `logOdds` where
`KL(FNCH(logOdds) || FNCH(nullLogOdds)) = log(1/alpha)`, with
`KL = (logOdds - nullLogOdds) * mean - fnchLogPartition(nb, na, logOdds) +
fnchLogPartition(nb, na, nullLogOdds)`, the FNCH mean of `yb` coming from
`BiasedUrn::meanFNCHypergeo(nb, na, totalSuccesses, exp(logOdds))`. The KL is
bounded, so when the target is out of reach (e.g. a degenerate conditional
distribution, or tiny tables) the solver returns `NULL` and the caller uses
the trivial e-factor `1`. The e-factor itself is
`logLikelihoodRatioFNCH` at the solved `logOdds`; the former wrapper
`savi2x2UmpStat` is removed and a `"twoSided"` rule (previously the plain
average of the two sides) is to be decided with the `eType` dispatch of
`savi2x2CondStat`. `alpha` is the target level of the one-shot test, not a
design field yet.

### 9. Per-block group sizes

Block sizes are observed data, like `ya` and `yb`. `designSavi2x2(na, nb)`
is unchanged: one positive integer each, the planned sizes, kept in
`nPlan` (`list(na, nb)`, Decision 12). `savi2x2TestStatPropDiff(ya, yb, na
= NULL, nb = NULL, designObj = NULL, ...)` and `savi2x2TestStatLogOdds`
(same signature): each of `na`, `nb` is `NULL` (the design's planned
size), one number (broadcast to `nBlocks`) or a vector of length `nBlocks`;
any other length errors. Sizes are positive integers with `ya[i] <= na[i]`
and `yb[i] <= nb[i]` per block. Observed sizes may differ from the planned
ones without a warning.

Internals always receive full-length vectors and block `i` uses `na[i]`,
`nb[i]`: the plug-in learners, the pooled projection, the RIPr solver and
the confidence sequence. The Beta posterior means of `predictiveThetas2x2`
divide by the cumulative size of blocks `1..i-1`, not `na * (i - 1)`.
Output: `n = c(na = sum(na), nb = sum(nb), nBlocks)`, the result's
`betaParameter` uses `sum(na)`, `sum(nb)`; the per-block vectors
are not stored on the result. The logOdds conditional e-factor is untouched.

### 10. Confidence sequence for logOdds

`logLikelihoodRatioFNCH(ya, yb, na, nb, logOdds, nullLogOdds = 0,
log = FALSE)` gains the null: the log e-factor is `log(dFNCHypergeo(yb, nb,
na, ya + yb, exp(logOdds)))` minus the same at `logOdds = nullLogOdds`, unchanged
at `nullLogOdds = 0`.

`computeConfidenceInterval2x2LogOdds(ya, yb, na, nb, betaParameter,
alpha, precision = 100, logOddsBound = 40)` inverts the conditional test on
`precision` equally spaced candidates strictly inside `(-logOddsBound,
logOddsBound)`. Each candidate is the `nullLogOdds` of its own e-process, the
product over blocks of the conditional e-factors. The plug-in alternative
for block `i` is predictable: `logit(thetaB) - logit(thetaA)` from the Beta
posterior means of `predictiveThetas2x2` given blocks `1..i-1` (the prior
means for block 1, so a symmetric prior gives the trivial factor 1 there).
The same plug-in serves every candidate. `runningIntersection`, runs and the
returned two-column, block-rowname matrix are as in Decision 4.
Wired into `savi2x2TestStatLogOdds` (Decision 12) via
`computeEValueVecLogOdds`, the same plug-in cumulated blockwise.

### 11. FNCH conditional e-factor via BiasedUrn

`conditionalEValueFixedAlternative` (Decisions 7, 10) is renamed
`logLikelihoodRatioFNCH` and computed as a log likelihood ratio of
`BiasedUrn::dFNCHypergeo` at `logOdds` over `nullLogOdds`, not via
`fnchLogPartition`; same arguments, same value (checked against the old
implementation to ~1e-9, floating-point noise). `dFNCHypergeo` is
near-constant time in block size where `fnchLogPartition`'s log-sum-exp is
`O(n)`, so this is the faster route for large blocks; `fnchLogPartition`
itself stays, since `solveUmpLogOdds` still needs it for the KL solve.

### 12. propDiff/logOdds split, and naming cleanup

`priorHyperParameters` is renamed `betaParameter` on the design
(`constructSaviDesignObj`); the test result's posterior (`constructSaviTestObj`)
is its own field, `betaPrior` (it is the posterior of the observed blocks
but also the prior a further block would use), resolving the earlier
name collision. `logOR`/`LogOR` is renamed `logOdds`/
`LogOdds` throughout, spelled like `propDiff` rather than abbreviated
(`solveUmpLogOR` → `solveUmpLogOdds`,
`computeConfidenceInterval2x2LogOR` → `computeConfidenceInterval2x2LogOdds`,
`logOR`/`nullLogOR` params → `logOdds`/`nullLogOdds`).

`savi2x2Test` is replaced by two effect-specific functions, both
`(ya, yb, na = NULL, nb = NULL, designObj = NULL, wantCi = TRUE)` per
Decision 9: `savi2x2TestStatPropDiff` (Decision 3's numerator/denominator,
Decision 6's restriction) and `savi2x2TestStatLogOdds` (Decision 7's
conditional e-factor, cumulated via `computeEValueVecLogOdds`, and
Decision 10's confidence sequence). Neither does its own restriction or
type checking yet (`# TODO: THESE ARGS CHECKING WILL BE DONE FINAL STEP`);
`savi2x2TestStatLogOdds` has no restricted alternative yet (`logOdds` stays
out of scope per Decision 6).

`savi2x2TestStat` is renamed `logLikelihoodRatioMultiBern`.
`logEProcess2x2PlugIn` is renamed `computeEValueVecPropDiff` and confirmed
to return the plain `eValueVec` (length `nBlocks`), not a list.

`designSavi2x2(na, nb, propDiffMin = NULL, logOddsMin = NULL, alpha, power,
h0, alternative, eType, betaParameter, runningIntersection)`: `eType` is
`"eBeta"`, `"grow"` or `"eGauss"` (Decision 2); at most one of
`propDiffMin`, `logOddsMin` is supplied and whichever is set is stored as
`esMin` (still no `effectMeasure` field — `eType` and which `*Min` is set
together imply the effect, `"grow"`'s dispatch between the two not yet
designed). `nPlan` is stored as `list(na, nb)`, not `c(nBlocks, na, nb)`
(reopens Decision 1): block count is observed per call, not planned, so it
no longer belongs in `nPlan`; `plot.saviDesign`, which reads `nPlan[1]` as
`nBlocks`, is not yet updated for this.

`computeConfidenceInterval2x2PropDiff` and `...LogOdds` (Decision 4, 10)
always return a two-column matrix, `lowerBound`/`upperBound` (`ncol = 2`):
block is the rowname, not a `"block"` data column, since a block may
contribute zero, one, or several rows (holes) and `rbind` needs a fixed
column count throughout. Callers read the block index with
`rownames(confSetRuns)`, e.g. `rownames(confSetRuns) == nBlocks`.
