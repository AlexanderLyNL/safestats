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
- `runningIntersection`: `TRUE` or `FALSE`, default `FALSE` in
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

Shelved in part: the two-sided `"grow"` average at `+/-``propDiffMin` is
removed for now, together with `logMeanExpPair`; `"grow"` is `"greater"`
only. The contract below is kept for when it returns.

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

### 13. Root-finding confidence sequence for propDiff

`computeConfidenceInterval2x2PropDiff(ya, yb, na, nb, betaParameter, alpha,
runningIntersection = TRUE)` replaces Decision 4's grid for `propDiff`
(`computeConfidenceInterval2x2LogOdds`, Decision 10, is untouched and stays
grid-based): no `precision` argument, no candidate grid. The eBeta
numerator (`predictiveThetas2x2`) does not depend on the candidate
`propDiff`, so for data through block `i` the cumulative log e-process, as
a function of the candidate `delta`,

    f_i(delta) := cumsum(computeEValueVecPropDiff(ya[1:i], yb[1:i], na[1:i],
      nb[1:i], betaParameter, propDiff = delta))[i]

is a sum over `j = 1..i` of the KL-projection terms of
`solveRIPr2x2PropDiff()`, each convex in `delta`; `f_i` is therefore convex
with one minimum. Checked numerically across balanced, extreme and
near-degenerate tables: `f_i` always has exactly one sign change in its
derivative, confirming this. Consequently `{delta : f_i(delta) <
log(1 / alpha)}` is always a single interval or empty inside `(-1, 1)` —
the multi-run ("holes") case of Decisions 4 and 12 cannot arise for this
numerator and is dropped for `propDiff` only.

Per block `i`: try a cheap interior guess first, the posterior mean of
`thetaB - thetaA` given blocks `1..i` (unlike the numerator's `thetas`,
which holds out block `i`); only if `f_i` there is `>= log(1 / alpha)` does
it fall back to locating the minimiser with `stats::optimize()` on
`(-1 + eps, 1 - eps)`. If `f_i` at the minimiser is still `>= log(1 /
alpha)` the block's raw interval is empty; otherwise find the lower and
upper bound with `stats::uniroot()` bracketed on either side of it (or
`-1`/`1` directly when `f_i` at that edge is already below the threshold).
No warm-starting or bracket reuse across blocks — simplicity over speed;
this recomputes the cumulative e-process from scratch at every block, so
runtime is quadratic in `nBlocks` (a few minutes at `nBlocks = 1000` in
`~/Downloads/local-run.R`), left as-is for now.

`runningIntersection = TRUE` (default): block `i`'s interval is the
intersection of its raw interval with block `i - 1`'s returned interval, so
a value that leaves the set never returns and the sequence stays nested,
matching Decision 4's semantics without per-candidate state. `FALSE`: each
block's interval is its raw root-finding result on the data seen so far,
independent of other blocks. Output shape is unchanged from Decision 12: a
two-column matrix, `lowerBound`/`upperBound`, block as the rowname, at most
one row per block (never several, per the argument above), absent for a
fully rejected block. `savi2x2TestStatPropDiff`'s wiring (hull into
`confSeqMatrix`, exact union as `confSeq` for the last block) needs no
change, since the hull of a single row is that row.

### 14. One confidence interval on all data for propDiff

`computeConfidenceInterval2x2PropDiff(ya, yb, na, nb, betaParameter, alpha)`
takes the full observed vectors and returns a single interval for the data
as a whole: a named numeric `c(lowerBound, upperBound)`, both `NA` when the
set is empty. No block loop, no `runningIntersection` argument, no matrix,
no block index. The candidate `delta` is kept when

    f(delta) = sum_{i=1}^{t} log [ p(ya_i, yb_i | thetaA_i, thetaB_i)
                                 / p(ya_i, yb_i | thetaA_i*(delta), thetaA_i*(delta) + delta) ]

is below `log(1 / alpha)`, with `thetaA_i, thetaB_i` the predictive Beta
posterior means from blocks `1..i-1` (`predictiveThetas2x2`) and
`thetaA_i*(delta)` the projection `solveRIPr2x2PropDiff(thetaA_i, thetaB_i,
na_i, nb_i, delta)`. `f` is convex (Decision 13), so the root finding is:
minimiser by `stats::optimize()` on `(-1, 1)` (no interior guess, for
simplicity), `stats::uniroot()` on each side, `-1` or
`1` when the edge is already inside.

`savi2x2TestStatPropDiff(wantCi = TRUE)` stores that vector as `confSeq`
and sets `ciValue = 1 - alpha`; it no longer sets `confSeqMatrix`, so
`plot.saviTest(wantConfSeqPlot = TRUE)` has nothing to draw for `propDiff`
for now. The blockwise sequence and its running intersection (Decision 13)
are shelved, not contradicted. `computeConfidenceInterval2x2LogOdds` and
`savi2x2TestStatLogOdds` are untouched.

### 15. GROW plug-in on logOdds

Shelved in part: the two-sided `"grow"` average at `+/-``logOddsMin` is
removed for now, together with `logMeanExpPair`; `"grow"` is `"greater"`
only. The contract below is kept for when it returns.

`designSavi2x2(logOddsMin, eType = "grow", alternative)`: `logOddsMin` is
`NULL` or one finite number strictly greater than `0`, stored as `esMin`
(Decision 12). `savi2x2TestStatLogOdds` reads `designObj[["eType"]]`;
only `"grow"` uses `esMin`, every other `eType` keeps the predictable
posterior-mean plug-in of Decision 10. Allowed combinations for `"grow"`;
everything else errors once argument checking is done:

- `logOddsMin > 0`, `"greater"`: block `i`'s e-factor is
  `logLikelihoodRatioFNCH(ya[i], yb[i], na[i], nb[i], logOdds =
  logOddsMin)`, the conditional e-factor of Decision 7 at the fixed
  alternative, cumulated over blocks.
- `logOddsMin > 0`, `"twoSided"`: the e-process is the average of the two
  cumulative e-processes at `+logOddsMin` and `-logOddsMin` (averaged as
  processes, not per block, as in Decision 6).
- `"less"` is not designed yet, as for `propDiff`.

`computeEValueVecLogOdds(ya, yb, na, nb, logOdds)` cumulates the
conditional e-factors at the fixed `logOdds` (the same value in every
block); the predictable posterior-mean plug-in it used to cumulate
(Decision 10) now lives only inside `computeConfidenceInterval2x2LogOdds`. The null is the point `thetaA = thetaB` (`nullLogOdds =
0`) throughout; the confidence sequence (Decision 10) keeps the
unrestricted plug-in, as Decision 6 does for `propDiff`. The two-sided
average is computed on the log scale by the shared helper
`logMeanExpPair(logX, logY)`, also used by `savi2x2TestStatPropDiff`.

### 16. Gaussian-mixture conditional e-process on logOdds (eGauss)

Shelved: all eGauss code is removed for now (`eType` is `"eBeta"` or
`"grow"`, and `savi2x2TestStatLogOdds` takes `"grow"` only). The contract
below is kept for when it returns.

`designSavi2x2(eType = "eGauss")` with `logOddsMin = NULL` and
`alternative = "twoSided"`; any `*Min` or one-sided alternative with
`"eGauss"` errors once argument checking is done. The prior on `logOdds`
is Normal(`priorMean = 0`, `priorSd = 1`) (reserved in Decision 2), held
on `nWeight = 2000` equally spaced grid points on `[-logOddsBound,
logOddsBound]`, `logOddsBound = 20`, and normalised over the grid. These
three are arguments of the helper with the stated defaults, not design
fields (as `nWeight` in Decision 6).

`computeEValueVecLogOddsGauss(ya, yb, na, nb, priorMean = 0, priorSd = 1,
nWeight = 2000L, logOddsBound = 20)` returns the cumulative log e-process
(length `nBlocks`): for block `i` it is
`logSumExp(priorLogWeights + sum_{j <= i} logLR_j(grid))`, where
`logLR_j(logOdds)` is the conditional log likelihood ratio of Decision 7 at
`logOdds` against `nullLogOdds = 0`, evaluated on the whole grid at once via
`logLikelihoodRatioFNCHGrid(ya, yb, na, nb, logOddsGrid)` (the
exponential-family form `logOdds * yb - psi(logOdds) + psi(0)`, `psi` the
log partition over the feasible `yb`, shifted column-wise against
overflow). This equals the reference `computeEGaussGrid` in
`~/Downloads/safe2x2TestCond.R`, which averages each block's likelihood
ratio under the running grid posterior, since the product of those
blockwise factors telescopes to the prior mixture of the cumulative
likelihood ratio; the reference weights `ya` (A minus B), so it is
reproduced by our code with the groups swapped. Blockwise e-factors are
`diff(c(0, logEValueVec))`, on the log scale.

`savi2x2TestStatLogOdds` dispatches on `eType`: `"eGauss"` uses this
mixture, `"grow"` Decision 15, any other `eType` errors.
`savi2x2TestStatPropDiff` likewise: `"eBeta"` or `"grow"`, else an error.
The confidence sequence (Decision 10) keeps its predictable plug-in for
every `eType`. The two test functions keep their size handling and
result filling inline, without shared helpers.

### 17. One confidence interval on all data for logOdds (grow)

`computeConfidenceInterval2x2LogOdds(ya, yb, na, nb, logEValue, alpha,
logOddsBound = 40)` replaces Decision 10's grid: no `betaParameter`, no
`precision`, no `runningIntersection`, no plug-in of its own, no block
loop. It inverts the **test's own** grow e-process on point nulls
`nullLogOdds = delta`, and returns a single named numeric
`c(lowerBound, upperBound)` for the data as a whole, as Decision 14 does
for `propDiff`; both `NA` when the set is empty.

The key identity: the grow alternative is fixed in advance, so the
conditional e-factor of block `i` against `delta` factors as
`LR_i(logOddsMin | 0) / LR_i(delta | 0)` (Decision 7's
`logLikelihoodRatioFNCH` with `nullLogOdds = delta`). Cumulated, the log
e-process against `delta` is

    f(delta) = logEValue - S(delta),
    S(delta) = sum_i logLikelihoodRatioFNCH(ya_i, yb_i, na_i, nb_i,
                                            logOdds = delta, nullLogOdds = 0, log = TRUE),

with `logEValue` the test's final cumulative log e-value against `0`
(the last element of `logEValueVec`, which `savi2x2TestStatLogOdds` already
holds). `S` is shared by the `+logOddsMin` and `-logOddsMin` processes, so
it also factors out of the two-sided average (Decision 15): one formula
for `"greater"` and `"twoSided"`. `S(delta)` is the conditional
log-likelihood of a common log odds ratio `delta`, concave in `delta`, so
`f` is convex with its minimiser at the conditional MLE, and the kept set
`{delta : f(delta) < log(1 / alpha)}` is one interval or empty. Root
finding as in Decision 14: `stats::optimize()` on
`(-logOddsBound, logOddsBound)`, `stats::uniroot()` on each side of the
minimiser. When `f` is still below the threshold at `-logOddsBound` or
`logOddsBound` the bound is reported as that edge, `-logOddsBound` or
`logOddsBound` (the range the interval is claimed on), not `-Inf`/`Inf`;
the empty set stays `NA`. At odds
`exp(40)` every block sits at its extreme feasible `yb`, so `f` has
integer slope there and can only stay below `log(1 / alpha)` if that
slope is `0`, in which case the set is genuinely unbounded. Blocks whose
conditional distribution is degenerate (a single feasible `yb`) have
`S_i = 0` for every `delta` and are skipped, an exact speed-up.

`savi2x2TestStatLogOdds(wantCi = TRUE)` computes this only for
`eType = "grow"`, stores the vector as `confSeq` and sets
`ciValue = 1 - alpha`; it no longer sets `confSeqMatrix`, so, as for
`propDiff`, `plot.saviTest(wantConfSeqPlot = TRUE)` has nothing to draw.
`"eGauss"` gets no confidence interval for now (the same identity holds
for any alternative fixed before the data, mixtures included, so it can
be wired in later with the same call). The design's
`runningIntersection` is now read by neither test. Decision 10's
blockwise sequence is shelved, not contradicted.

### 18. eGauss on logOdds, fixed-grid marginal

Reopens Decision 16, `"twoSided"` only. `designSavi2x2(eType = "eGauss")`
with `logOddsMin = NULL`. Hardcoded inside `savi2x2TestStatLogOdds`: prior
Normal(0, 1) on 2000 equally spaced `logOdds` in `[-20, 20]`, normalised
over the grid. The numerator after block `i` is the grid mixture of the
cumulative conditional likelihood,
`logPCum[i] = logSumExp(logPrior + sum_{j <= i} logP_j(grid))`, with
`logP_j` the FNCH log density of `yb_j` (B minus A), computed per block for
the whole grid at once (`outer` over grid and feasible `yb`, shifted by the
row max). No posterior loop: the product of the predictive factors
telescopes to this mixture. `logEValueVec = logPCum - cumsum(logP0)`. The
confidence interval is Decision 17 unchanged, with the numerator total
`logPCum[nBlocks]` in place of grow's. Checked against `cond`'s
`computeEGaussGrid` with the groups swapped. No new helpers.

### 19. No confidence interval for grow on logOdds

Reopens Decision 17 for grow: `savi2x2TestStatLogOdds(eType = "grow")`
sets no `confSeq` or `ciValue`, whatever `wantCi`; its construction is
still under discussion. eGauss keeps Decision 17's interval (Decision 18).

### 20. propDiff e-process inline in the test function

`savi2x2TestStatPropDiff` computes its e-process in three inline steps:
the predictable numerator thetas (`predictiveThetas2x2` for `"eBeta"`,
twoSided; `predictiveThetas2x2PropDiff` at `propDiffMin` for `"grow"`,
greater only, twoSided later), the pooled projection
`(na * thetaA + nb * thetaB) / (na + nb)` as the null, and
`cumsum` of the binomial log likelihood ratio via `stats::dbinom(log =
TRUE)`. `logLikelihoodRatioMultiBern` and `computeEValueVecPropDiff`
(Decisions 4, 6, 12) are removed; `computeConfidenceInterval2x2PropDiff`
sums the same `dbinom` terms inline. Output unchanged to ~1e-15.

### 21. propDiff confidence interval inline, eBeta only

Replaces Decision 14's function: `computeConfidenceInterval2x2PropDiff` is
removed and `savi2x2TestStatPropDiff` builds the interval inline, as
`savi2x2TestStatLogOdds` does, only for `eType = "eBeta"` (grow gets no
`confSeq` or `ciValue`). `fPropDiff(thetaA, thetaB, propDiff, alpha)` is
the log e-value on all blocks against `thetaB - thetaA = propDiff`, minus
`log(1 / alpha)`: the test's own eBeta numerator thetas, and per block the
null `thetaA` from `solveRIPr2x2PropDiff`. One interval for the last block
(all data): `optimize()` on `(-1, 1)`, `uniroot()` on each side, `-1`/`1`
when the edge is still inside, `NA` when empty. Output unchanged to ~1e-15.

### 22. Data checks and estimate in both test functions

`savi2x2TestStatPropDiff` and `savi2x2TestStatLogOdds`, after broadcasting
`na`, `nb` from the design, stop unless `ya`, `yb`, `na`, `nb` all have
length `nBlocks`, are finite integers, `0 <= ya <= na`, `0 <= yb <= nb`
and `na, nb >= 1`. The other argument checks stay for the final step.
Both return `estimate = c(thetaA, thetaB)`, the observed pooled
proportions (`propDiff` dropped from the estimate).

### 23. propDiff confidence interval isolated again

Reverses the inline placement of Decision 21, same statistics:
`computeConfidenceInterval2x2PropDiff(ya, yb, na, nb, thetaA, thetaB,
alpha)`, under `# Confidence Interval ----`, takes the test's own eBeta
numerator thetas (one per block) instead of `betaParameter`, so they are
computed once. Inside, `fPropDiff(propDiff)` is the log e-value on all
blocks against `thetaB - thetaA = propDiff` (per-block
`solveRIPr2x2PropDiff` nulls) minus `log(1 / alpha)`. Returns
`c(lowerBound, upperBound)`: `-1`/`1` when the edge is still inside, both
`NA` when empty. `savi2x2TestStatPropDiff` calls it only for eBeta.

### 24. Blockwise confidence sequence for propDiff

`savi2x2TestStatPropDiff(..., wantConfidenceSequence = FALSE)`: when
`TRUE` and `eType = "eBeta"`, it stores `confSeqMatrix`, `nBlocks x 2`
(`lowerBound`, `upperBound`), row `i` from
`computeConfidenceInterval2x2PropDiff` on blocks `1..i` (the test's thetas
restricted to the first `i` blocks, which are predictable, so unchanged).
The design's `runningIntersection` (default `FALSE`) intersects row `i`
with row `i - 1`, so the rows are nested and a fully rejected row is `NA`
from then on; `FALSE` keeps each row as computed. `confSeq` is then the
last row and `ciValue = 1 - alpha`; otherwise `wantCi` gives Decision 23's
interval on all data. Runtime is quadratic in `nBlocks` (about 0.4 ms x
nBlocks^2: ~7 min at 1000 blocks), hence off by default. Reopens the
shelved blockwise sequence of Decision 13 on top of Decision 23.

### 25. propDiff interval takes betaParameter and a domain

`computeConfidenceInterval2x2PropDiff(ya, yb, na, nb, betaParameter,
alpha, domain = c(-1, 1))`: computes the eBeta numerator thetas itself
with `predictiveThetas2x2` (cheap; replaces Decision 23's `thetaA`,
`thetaB` arguments) and searches only inside `domain` (kept `1e-9` inside
`(-1, 1)` for the projection). By convexity the result is the kept set
intersected with `domain`; a bound still inside at the edge is reported
as that `domain` edge, the empty set as `NA`. In Decision 24's sequence,
`runningIntersection = TRUE` passes the previous row as `domain` instead
of intersecting afterwards, and stops at the first empty row (the rest
stay `NA`); `FALSE` always uses `c(-1, 1)`. Running-intersection rows
match the previous code to ~3e-5, `uniroot()`'s default tolerance.

### 26. Grid-cumulated confidence sequence for propDiff

Reopens Decision 4's grid, with the projection vectorised, for the
blockwise sequence only; Decision 23's single interval on all data and
`solveRIPr2x2PropDiff` are unchanged.
`computeConfidenceSequence2x2PropDiff(ya, yb, na, nb, betaParameter,
alpha, runningIntersection)` returns the `nBlocks x 2` matrix
(`lowerBound`, `upperBound`) that `savi2x2TestStatPropDiff(
wantConfidenceSequence = TRUE)` stores as `confSeqMatrix` (Decision 24),
now linear in `nBlocks`:

- Grid: `nGrid = 2000` candidates `delta` equally spaced strictly inside
  `(-1, 1)`, always. `sdMax = sqrt(1 / (4 sum(na)) + 1 / (4 sum(nb)))`,
  the worst-case Wald standard deviation of the difference on the observed
  totals, is the scale of the last (narrowest) interval; when it would call
  for more than 2000 points (`ceiling(2 / sdMax) > 2000`, i.e. the step
  `1e-3` is coarser than `sdMax`) the grid stays at 2000 and a warning
  says the bounds are conservative. `nGrid` is not a design field.
- Cumulation: the log e-process against every candidate is held on the
  grid and block `i` adds its term once. The null `thetaA` for all
  candidates at once is the single root on the feasible interval of the
  cubic `na (x - thetaA)(x + delta)(1 - x - delta) + nb (x + delta -
  thetaB) x (1 - x)` (the KL derivative of `solveRIPr2x2PropDiff` times
  `x (1 - x)(x + delta)(1 - x - delta)`), found by vectorised bisection.
- Bounds: the kept candidates form one run (convexity, Decision 13).
  Each bound is refined outward by the secant through the first rejected
  node and its outer neighbour: it lies below the convex `f` outside them,
  so its root brackets the true boundary from outside and the reported
  interval contains the exact one. `-1`/`1` when the outermost candidate
  is kept, `NA` when none is.
- `runningIntersection = TRUE`: row `i` is the intersection with row
  `i - 1`, rejected candidates are dropped from the update (they never
  return), and the rows after the first empty one stay `NA`. `FALSE`:
  every candidate stays active and each row is its own raw interval.

Matches the root-finding rows of Decision 25 to ~3e-5 (`uniroot`'s
tolerance) on ordinary tables, conservative by at most one grid step on
extreme ones; about 3 s at 1000 blocks against ~7 min.

### 27. logOdds confidence interval isolated

The inline eGauss interval of Decisions 17 and 18 moves out of
`savi2x2TestStatLogOdds` into `computeConfidenceInterval2x2LogOdds(ya,
yb, na, nb, logPTotal, alpha, logOddsBound = 40)`, under `# Confidence
Interval ----`, as Decision 23 did for `propDiff`. `logPTotal` is the
numerator's total log likelihood after the last block (`logPCum[nBlocks]`),
and `fLogOdds(delta)` is `logPTotal` minus the FNCH log likelihood at
`delta` minus `log(1 / alpha)`; root finding is unchanged. Returns
`c(lowerBound, upperBound)`, `-logOddsBound` / `logOddsBound` when the
edge is still inside. The empty set is reported as the whole range with a
warning, as `computeConfidenceInterval2x2PropDiff` now reports `c(-1, 1)`
(reopens the `NA` of Decisions 17 and 23). `savi2x2TestStatLogOdds` calls
it for eGauss only; grow still gets no interval (Decision 19).

### 28. UMP plug-in at the first block of eGauss on logOdds

No new `eType`. `savi2x2TestStatLogOdds(eType = "eGauss")` (Decision 18)
replaces only the **first block's** e-factor: instead of the prior mixture
(which for the symmetric N(0, 1) prior is a wasted block) it uses the
conditional e-factor at the UMP plug-in `logOddsUmp = solveUmpLogOdds(na[1],
nb[1], ya[1] + yb[1], alpha, "greater")` (Decision 8), the `logOdds > 0` at
which the conditional KL against the null reaches `log(1 / alpha)`. It
depends on block 1 only through its total, which the FNCH e-factor
conditions on, so it is a valid conditional e-variable. The Bayesian
updating is unchanged: the grid posterior still absorbs block 1, and blocks
`2..i` keep their predictive mixture factors, so

    logPCum[1] = log dFNCH(yb_1; nb_1, na_1, ya_1 + yb_1, exp(logOddsUmp)),
    logPCum[i] = logPCum[1] + logMixCum[i] - logMixCum[1]   (i >= 2),

with `logMixCum` Decision 18's grid mixture of the cumulative conditional
likelihood. When `solveUmpLogOdds` returns `NULL` (the KL target is out of
reach, e.g. a degenerate total) block 1's factor is `1`, i.e. `logPCum[1] =
logPNull[1]`. The plug-in is solved at the design's `alpha`. `"greater"` is
used for the first block even though eGauss is `"twoSided"`; the two-sided
UMP rule (Decision 8) is still open, as is the same replacement for
`propDiff` and for grow (grow is untouched: fixed `logOddsMin` in every
block).

The confidence interval is Decision 17's identity unchanged: every factor
is fixed given its block's total or predictable, so
`computeConfidenceInterval2x2LogOdds(ya, yb, na, nb, logPCum[nBlocks],
1 - ciValue)` inverts the new numerator directly, and the kept set is still
one interval. The inversion level is `1 - ciValue`, the test's `ciValue`
argument (default `1 - alpha`), for eGauss.

### 29. Stopping-time simulation for propDiff grow (stub)

`sampleStoppingTimesSavi2x2(propDiffMin, power, na, nb, alpha = 0.05,
betaParameter = NULL, nTheta = 8L, nSim = 1e3L, maxBlocks = 1e4L, seed =
NULL)`, under `# Sampling functions for design ----`. Plans the grow test on
`propDiff`, `"greater"` only: `propDiffMin` is both the numerator plug-in
and the data-generating effect, as in the t-test design and `cond`. Data
lie on the curve `thetaB = thetaA + propDiffMin` at `nTheta` baselines
`thetaA = rho * (1 - propDiffMin)`, `rho` equally spaced strictly inside
`(0, 1)`. Per baseline, `nSim` paths of `maxBlocks` blocks, `ya ~
Binom(na, thetaA)`, `yb ~ Binom(nb, thetaB)`; each path runs the same
e-process as `savi2x2TestStatPropDiff(eType = "grow")` and stops at the
first block with `E >= 1 / alpha`, `Inf` when none. One standalone
per-path block loop, no chunking. Returns a list: `thetaA`, `thetaB`
(length `nTheta`), `stoppingTimes` (`nTheta x nSim`), `nPlan` the
`ceiling` of the largest `power` quantile over the baselines (the worst
case), `worstCaseIndex`. No bootstrap SE and no eBeta planning yet; the
design does not call it yet. Body is a stub until the contract is
confirmed.

### 30. Blockwise confidence sequence for logOdds (naive)

`savi2x2TestStatLogOdds(..., wantConfidenceSequence = FALSE)`: when
`TRUE` and `eType = "eGauss"`, it stores `confSeqMatrix`, `nBlocks x 2`
(`lowerBound`, `upperBound`), row `i` from
`computeConfidenceInterval2x2LogOdds(ya[1:i], yb[1:i], na[1:i], nb[1:i],
logPCum[i], 1 - ciValue)`, Decision 27's inversion on blocks `1..i`.
Decision 17's identity holds at every block, not only the last: each
numerator factor is fixed given its block's total (the UMP plug-in of
block 1, Decision 28) or predictable (the mixture factors of blocks
`2..i`), so the log e-process against `delta` after block `i` is
`logPCum[i] - S_i(delta)`, and the test already holds `logPCum` for every
block. Each row's kept set is one interval or empty (concavity of `S_i`).
`computeConfidenceInterval2x2LogOdds` gains `domain = c(-40, 40)` in place
of `logOddsBound`, as Decision 25's `domain` for `propDiff`: the search
range, whose edge is reported when still inside, and the whole `domain`
with a warning when empty. In the sequence an empty row is `NA` (the
warning is caught, not shown). The design's `runningIntersection`
(default `FALSE`) passes the previous row as `domain`, so the rows are
nested, and stops at the first empty row (the rest stay `NA`); `FALSE`
searches `(-40, 40)` every block. `confSeq` is then the last row,
`ciValue` unchanged; without the flag `wantCi` gives Decision 27's
interval. Grow still gets nothing (Decision 19). One `optimize` and up to
two `uniroot` per block, each summing `i` FNCH log densities: quadratic
in `nBlocks`, hence off by default. `plot.saviTest` labels the
confidence-sequence axis `logOdds` when the design's `eType` is not
`"eBeta"`. The linear-time grid cumulation (Decision 26's route, reusing
the eGauss `logPGrid`) is the follow-up.
