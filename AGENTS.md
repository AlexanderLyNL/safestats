# Working agreements

## Start here

For the 2x2 rewrite, read [the current design](docs/2x2-design.md) before
working, and the relevant entries in [the review findings](docs/2x2-review.md)
before changing affected code. The design records approved behavior and the
next agreed action; the review records unresolved defects. Neither a TODO nor
an issue is permission to implement a new feature or a deferred fix.

## Objective and scope

Rewrite the anytime-valid two-proportion test as a small, direct implementation.
Priorities: statistical validity, explicit conventions, numerical stability,
focused tests, then everything else. Preserve legacy behavior only when correct.
The `cond` branch is reference material: use `git show cond:<path>`, never merge it.

New 2x2 code belongs in `R/newsafe2x2Test.R`; leave `R/safe2x2Test.R` untouched.
Scope includes its tests, roxygen/man pages, these working documents, and minimal
2x2 edits to `NAMESPACE`, `DESCRIPTION`, `R/safeS3Methods.R`, and `R/deprecate.R`.
Other tests, general design infrastructure, unrelated S3 methods, vignettes, and
package-wide cleanup are out of scope unless requested.

## Agree, record, implement

Describe a feature in plain language, settle questions with the user, update its
current contract in `docs/2x2-design.md`, then implement. Do not append another
chronological decision for routine edits. User authorization persists across turns.
If code and contract disagree, fix the code or reopen the contract; do not silently
change the intended behavior. A documented mathematical claim can be wrong:
record counterexamples and seek an agreed replacement, rather than enforcing it.
Distinguish proofs, exact checks, and simulation evidence.

## Statistical conventions

- Groups are A and B; `ya`, `yb` are successes and `na`, `nb` are group sizes.
  Counts are finite nonnegative integers no larger than their positive sizes.
- Effects are exactly `propDiff` and `logOdds`, both B minus A, anchored on A:
  `thetaB = thetaA + propDiff`; `thetaB = plogis(qlogis(thetaA) + logOdds)`.
  Do not introduce legacy effect spellings in new code.
- Order alternatives as `twoSided`, `greater`, `less` wherever supported.
  The design specifies which combinations are currently available.
- Distinguish blockwise e-factors from cumulative e-processes. Predictive
  numerators use previous blocks; a conditional factor may use its conditioned total.
- Average cumulative processes, not their blockwise factors. Multiplying two
  factors based on the same block needs a validity argument.
- Work on the log scale where needed; handle zero probabilities and infinities
  explicitly. Numerical spot checks do not prove convexity or coverage.

## Implementation and verification

Use 2-space indentation, `<-`, spaces around operators, camelCase, `match.arg()`,
and guard clauses. Comments explain statistical intent. Exported roxygen uses
Markdown and states effect direction and return semantics. Follow `R/tTest.R`'s
layout: statistic, S3 generic/default/formula, alias, design, then sampling.

Write tests only when asked. Prefer small exact public-API checks over grids:
hand-computable tables, first-block consistency, log/plain agreement, mirrored
group swaps, boundaries, and seeded simulation. Skip scaffolding tests.
Never loosen tolerances before understanding a discrepancy.
Run existing 2x2 checks with `Rscript -e 'testthat::test_local(".", filter = "2x2")'`;
report a missing suite as missing, not passing. Regenerate roxygen with
`LANG=en_US.UTF-8 Rscript -e 'devtools::document()'` to avoid locale corruption.
When the API changes, keep `~/Downloads/local-run.R` working and run it from the
package root after `devtools::load_all()`; report inaccessible files honestly.

## Documentation

Keep this file about one page and the design about 2–4 pages. Update the design
in place; keep detailed issue evidence in the review. Use Git for history, not a
session log; the old decisions remain available through
`git show e4464b587ea6:AGENTS.md`.
