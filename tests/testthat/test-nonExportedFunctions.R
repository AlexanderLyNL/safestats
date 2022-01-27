test_that("checkAndReturnsEsMinParameterSide throws a warning", {
  ### Use in the design stage
  meanDiffMin <- 0.4
  paramChecked <- safestats:::checkAndReturnsEsMinParameterSide(meanDiffMin,
                                                                esMin="meanDiffMin",
                                                                alternative="two.sided")
  # Returns absolute value of meanDiffMin
  paramChecked

  expect_equal("object"=paramChecked, "expected"=0.4)

  # Invokes warnings
  paramChecked <- safestats:::checkAndReturnsEsMinParameterSide(meanDiffMin,
                                                                esMin="meanDiffMin",
                                                                alternative="greater")
  paramChecked == meanDiffMin
  #
  # ### Use in the execution stage
  phiS <- -0.3
  expect_warning(safestats:::checkAndReturnsEsMinParameterSide(phiS,
                                                               alternative="greater"))
})

test_that("checkAndReturnsNPlan throws a warning", {
  expect_warning(safestats:::checkAndReturnsNPlan(nPlan=5, testType="twoSample"))
})

test_that("tryOrFailWithNA returns correct value", {
  result <- safestats:::tryOrFailWithNA(integrate(exp, -Inf, Inf)[["value"]], NA)
  expect_equal("object"=result, "expected"=NA)

  result <- safestats:::tryOrFailWithNA(integrate(exp, 0, 3)[["value"]], NA)
  expect_equal("object"=result, "expected"=19.085537)
})

test_that("getArgs returns the arguments correctly", {
  foo <- function(x, y) {
    safestats:::getArgs()
  }

  result <- foo(x="3", y=df)

  expect_equal(object=result$x, expected="3")
  expect_equal(object=class(result$y), expected="name")
})

test_that("extractNameFromArgs returns the arguments correctly", {
  foo <- function(x, y) {
    safestats:::getArgs()
  }

  result <- foo(x="3", y=df)

  expect_equal(object=result$x, expected="3")
  expect_equal(object=class(result$y), expected="name")
})

test_that("extractNameFromArgs returns the arguments correctly", {
  foo <- function(x, y) {
    safestats:::getArgs()
  }

  bar <- foo(x="3", y=rnorm(10))
  result <- safestats:::extractNameFromArgs(bar, "y")
  expect_equal(object=result, expected="rnorm(10)")
})

test_that("getNameTestType returns the correct name", {
  result <- safestats:::getNameTestType("oneSample", "deltaS")
  expect_equal(object=result, expected="Safe One Sample T-Test")
})

test_that("getNameAlternative returns the correct alternative", {
  result <- safestats:::getNameAlternative("two.sided", testType="oneSample")
  expect_equal(object=result, expected="true mean not equal to 0")
})

test_that("computeNPlanBatchSafeT returns the correct batch sample size", {
  result <- safestats:::computeNPlanBatchSafeT(0.4)
  nPlan <- 88
  names(nPlan) <- "n1Plan"
  expectedResult <- list(nPlan=nPlan, deltaS=0.4)
  expect_equal(object=result, expected=expectedResult)
})

test_that("computeEsMinSafeT throws an error", {
  expect_error(safestats:::computeEsMinSafeT(3))
})

test_that("defineTTestN returns correct list", {
  result <- safestats:::defineTTestN()
  expectedResult <- list(n1=3:100, n2=NULL, nEff=3:100, nu=2:99)
  expect_equal(object=result, expected=expectedResult)
})

test_that("computeNPlanBatchSafeZ returns correct list", {
  result <- safestats:::computeNPlanBatchSafeZ(0.4)

  nPlan <- 85
  names(nPlan) <- "n1Plan"
  expectedResult <- list(nPlan=nPlan, phiS=0.4)

  expect_equal(object=result, expected=expectedResult)

  result <- safestats:::computeNPlanBatchSafeZ(0.4, grow=FALSE)

  nPlan <- 80
  names(nPlan) <- "n1Plan"

  expectedResult <- list(nPlan=nPlan, phiS=0.2691, lowN=67, highN=80, lowParam=0.2, highParam=0.4)

  expect_equal(object=result, expected=expectedResult)
})


test_that("computeBetaBatchSafeZ returns correct batch beta", {
  result <- safestats:::computeBetaBatchSafeZ(meanDiffMin=0.9, nPlan=12)
  expectedResult <- 0.35359335
  expect_equal(object=result, expected=expectedResult)
})


test_that("computeMinEsBatchSafeZ returns correct batch minEs", {
  result <- safestats:::computeMinEsBatchSafeZ(nPlan=78)
  expectedResult <- 0.41726942
  expect_equal(object=result, expected=expectedResult)
})


test_that("computeNEff returns correct effective sample size", {
  result <- safestats:::computeNEff(c(3, 4), testType="twoSample")
  expectedResult <- 1.7142857
  expect_equal(object=result, expected=expectedResult)
})
