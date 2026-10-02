test_that("FNCH weights are 1, 8, 4 at odds 2", {
  # probability of FNCH
  p <- c(1, 8, 4) / 13
  expect_equal(BiasedUrn::dFNCHypergeo(0:2, 2, 2, 2, 2), p)

  partitionConst <- 0
  totalSuccesses <- 2
  na <- 2
  nb <- 2
  kMin <- max(0, totalSuccesses - nb)
  kMax <- min(na, totalSuccesses)
  odds <- 2
  for (k in seq(kMin, kMax)) {
    # choose(na, k) choose(na, totalSuccess-k) 2^k
    partitionConst <- partitionConst + choose(2, k) * choose(2, k) * odds**k
  }
  expect_equal(
    fnchLogPartition(
      na = na,
      nb = nb,
      totalSuccesses = totalSuccesses,
      logOdds = log(odds)
    ),
    log(partitionConst)
  )
  expect_equal(
    logLikelihoodFNCH(
      ya = 0:na,
      yb = nb:0,
      na = na,
      nb = totalSuccesses,
      logOdds = log(odds)
    ),
    log(p)
  )
})

test_that("UMP recovers odds 3 for a single pair", {
  # Null P(A succeeds | total = 1) = 1/2; at odds 3 it is 3/4.
  # Choose alpha from KL so that the root is exactly log(3).
  # (logOdds - logOddsNull) * mean - fnchLogPartition(logOdds) + fnchLogPartition(logOddsNull) = log(1/alpha)
  #
  # na = nb = 1, total 1: K in {0, 1} with P_delta(K = 1) = e^delta / (1 + e^delta),
  # so Z(delta) = 1 + e^delta, Z(0) = 2 and E_delta[K] = P_delta(K = 1).
  # At delta = log(3): E[K] = 3/4, and the KL from the null is
  #   log(3) * 3/4 - log(4) + log(2) = ...
  # Setting log(1/alpha) equal to it makes log(3) the root; the "less" root
  # is -log(3) by the symmetry na = nb. The e-value of the table (1, 0) is
  #   P_log3(K = 1) / P_0(K = 1) = (3/4) / (1/2) = 3/2.
  mean <- BiasedUrn::dFNCHypergeo(1, 1, 1, 1, 3)
  alpha <- exp(-(3 / 4 * log(3 / 2) + 1 / 4 * log(1 / 2)))
  expect_equal(
    solveUmpLogOdds(
      na = 1,
      nb = 1,
      totalSuccesses = 1,
      alpha = alpha,
      alternative = "greater"
    ),
    log(3)
  )
  expect_equal(solveUmpLogOdds(1, 1, 1, alpha, "less"), -log(3))
  expect_equal(savi2x2TestStatUmp(1, 0, 1, 1, alpha, "greater"), 3 / 2)
})

test_that("log-odds grow multiplies by 3/2 for each A-only success", {
  # Each block is the table (1, 0) with na = nb = 1, total 1. At logOdds
  # log(3) the conditional probability of ya = 1 is 3 / (1 + 3) = 3/4; under
  # the hypergeometric null it is 1/2. Every block multiplies by
  # (3/4) / (1/2) = 3/2, so the cumulative e-value after t blocks is (3/2)^t.
  ya <- na <- nb <- rep(1, 5)
  yb <- rep(0, 5)
  actual <- logEValueVec2x2LogOddsGrow(ya, yb, na, nb, log(3), "greater")
  expect_equal(exp(actual), (3 / 2)^(1:5))
})

test_that("eGauss learns from five A-only successes", {
  ya <- na <- nb <- rep(1, 5)
  yb <- rep(0, 5)
  actual <- logEValueVec2x2LogOddsEGauss(
    ya,
    yb,
    na,
    nb,
    gaussParameter = list(mean = 0, sd = 1),
    alternative = "twoSided",
    logOddsGrid = seq(-8, 8, length.out = 161)
  )
  # same five tables ya = 1, yb = 0, na = nb = 1:
  # 1. One table at log odds z. Given the total 1, ya is 0 or 1 with FNCH
  #    weights choose(1, 0) choose(1, 1) e^0 = 1 and choose(1, 1) choose(1, 0) e^z = e^z,
  #    so P_z(ya = 1) = e^z / (1 + e^z) = plogis(z). The null z = 0 gives 1/2.
  # 2. The likelihood ratio of one table is plogis(z) / (1/2) = 2 plogis(z);
  #    t identical tables give (2 plogis(z))^t.
  # 3. eGauss averages this over the prior
  # 4. t = 1 is exactly 1: z and w are symmetric about 0 and
  #    plogis(z) + plogis(-z) = 1, so sum(w * 2 * plogis(z)) = sum(w) = 1.
  # 5. t >= 2 has no closed form. The grid sum is a Riemann sum for
  #    E[(2 plogis(Z))^t] with Z ~ N(0, 1), which agrees to 1e-14:
  #      integrate(function(z) dnorm(z) * (2 * plogis(z))^t, -Inf, Inf)
  #
  # E[(2 * plogis(Z))^t], Z ~ N(0,1): first is 1 by symmetry.
  # The remaining values are numerical Gaussian integrals, not hand formulas.
  expected <- c(
    1,
    1.173516143432372,
    1.520548430297116,
    2.105509009062085,
    3.057222176662974
  )

  # How `expected` is calculated (step 3), in the style of the first block.
  logOddsGrid <- seq(-8, 8, length.out = 161)
  # prior N(0, 1) evaluated on the grid and normalised to sum to 1
  priorWeights <- dnorm(logOddsGrid) / sum(dnorm(logOddsGrid))
  expectedByHand <- numeric(5)
  for (t in 1:5) {
    # one table: plogis(z) / (1/2) = 2 plogis(z); t tables: (2 plogis(z))^t;
    # then average over the prior
    expectedByHand[t] <- sum(priorWeights * (2 * plogis(logOddsGrid))^t)
  }
  expect_equal(expectedByHand, expected, tolerance = 1e-12)
  expect_equal(exp(actual), expected, tolerance = 1e-10)
})

test_that("uniform Beta priors update after one A-only success", {
  prior <- list(betaA1 = 1, betaA2 = 1, betaB1 = 1, betaB2 = 1)
  ya <- na <- nb <- c(1, 1)
  yb <- c(0, 0)
  # Beta(1, 1) posterior mean (1 + y) / (2 + n): block 1 uses the prior,
  # 1/2 for both groups; block 2 has seen ya = 1 of na = 1 and yb = 0 of
  # nb = 1, so thetaA = 2/3 and thetaB = 1/3.
  theta <- predictiveThetas2x2(ya, yb, na, nb, prior)
  expect_equal(theta$thetaA, c(1 / 2, 2 / 3))
  expect_equal(theta$thetaB, c(1 / 2, 1 / 3))
  # Block 1 contributes 1; block 2 contributes (2/3)^2 / (1/2)^2.
  # Block 1: thetaA = thetaB = 1/2, so the pooled thetaNull is 1/2 as well
  # and the ratio is 1.
  # Block 2: numerator thetaA^ya (1 - thetaB)^(nb - yb) = (2/3)(2/3) = 4/9;
  # pooled thetaNull = (1 * 2/3 + 1 * 1/3) / 2 = 1/2 gives (1/2)(1/2) = 1/4;
  # the ratio is (4/9) / (1/4) = 16/9.
  actual <- logEValueVec2x2PropDiffEBeta(ya, yb, na, nb, prior)
  expect_equal(exp(actual), c(1, 16 / 9))
})

# TODO: not zero differences cases
test_that("RIPr with zero difference is the pooled mean", {
  expect_equal(
    solveRIPr2x2PropDiff(
      thetaA = 3 / 4,
      thetaB = 1 / 4,
      na = 1,
      nb = 1,
      propDiff = 0
    ),
    1 / 2
  )
})
