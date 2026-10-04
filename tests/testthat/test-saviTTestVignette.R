test_that("The paired t-test design objects and saviTTest yield the same numbers as in the vignette", {
  alpha <- 0.05
  power <- 0.8
  deltaMin <- 12/(sqrt(2)*15)
  sigma <- 15

  designObj <- designSaviT(deltaMin=deltaMin, alpha=alpha,
                           power=power, sigma=sigma,
                           alternative="greater",
                           testType="paired", seed=1, pb=FALSE)

  expect_equal("object"=unname(designObj$parameter), "expected"=0.16)

  expect_equal("object"=unname(designObj$nPlan), "expected"=c(36, 36))
  expect_equal("object"=unname(designObj$nPlanBatch[1]-designObj$nPlan[1]),
               "expected"=8)

  designObj2 <- designSaviT(deltaMin=deltaMin, alpha=alpha,
                            nPlan=c(30, 30), sigma=sigma,
                            alternative="greater",
                            testType="paired", seed=1, pb=FALSE)

  expect_equal("object"=designObj2$power, "expected"=0.719)

  designObj3 <- designSaviT(nPlan=c(50, 50),
                            alpha=alpha, power=power,
                            sigma=sigma,
                            alternative="greater",
                            testType="paired")

  expect_equal("object"=round(unname(designObj3$esMin), 7),
               "expected"=0.5245391)

})

test_that("The paired saviTTest yield the same numbers as in the vignette", {
  alpha <- 0.05
  power <- 0.8
  deltaMin <- 12/(sqrt(2)*15)
  sigma <- 15

  designObj <- designSaviT(deltaMin=deltaMin, alpha=alpha,
                           sigma=sigma,
                           alternative="greater",
                           testType="paired", seed=1, pb=FALSE)

  # saviTTest ------
  set.seed(1)
  preData <- rnorm(n=36, mean=120, sd=15)
  postData <- rnorm(n=36, mean=120, sd=15)
  res <- saviTTest(x=preData, y=postData,
                   designObj=designObj, paired=TRUE, sequential=FALSE)

  expect_equal("object"=round(unname(res$eValue), 8),
               "expected"=0.02019529)
  expect_equal("object"=res$confSeq[1], "expected"=-12.79129930)
  expect_equal("object"=res$confSeq[2], "expected"=7.63151051)
  expect_equal("object"=unname(res$statistic), "expected"=-0.783633940)
  expect_equal("object"=res$meanObs, "expected"=-2.57989440)

  nSim <- 1000

  set.seed(1)
  eValues <- replicate(n=nSim, expr={
    preData <- rnorm(n=36, mean=120,
                     sd=15)
    postData <- rnorm(n=36, mean=120,
                      sd=15)
    saviTTest(x=preData, y=postData,
              designObj=designObj,
              paired=TRUE, sequential=FALSE)$eValue}
  )

  expect_equal("object"=sum(eValues > 20), "expected"=3)
})

test_that("The sampleStoppingTimesSaviT yields the same numbers as in the vignette", {
  alpha <- 0.05
  power <- 0.8
  deltaMin <- 12/(sqrt(2)*15)
  sigma <- 15

  alpha <- 0.05
  power <- 0.8
  deltaMin <- 12/(sqrt(2)*15)
  sigma <- 15
  nSim <- 1000

  designObj <- designSaviT(deltaMin=deltaMin, alpha=alpha,
                           sigma=sigma,
                           alternative="greater",
                           testType="paired", seed=1, pb=FALSE)

  # Null ----
  simDeltaTrueIsDeltaMin <-
    sampleStoppingTimesSaviT(
      deltaTrue=deltaMin, alternative="greater",
      testType="paired", sigma=sigma,
      nMax=36, seed=1,
      parameter=designObj$parameter,nSim=nSim)


  expect_equal("object"=sum(simDeltaTrueIsDeltaMin$stoppingTimes < 36/2), "expected"=412)
  expect_equal("object"=mean(simDeltaTrueIsDeltaMin$stoppingTimes), "expected"=21.672)

  # Alternative ----
  simDeltaTrueLargerDeltaMin <-
    sampleStoppingTimesSaviT(
      deltaTrue=1.2*deltaMin, alternative="greater",
      testType="paired", sigma=sigma,
      nMax=36, seed=1,
      parameter=designObj$parameter,nSim=nSim)

  expect_equal("object"=sum(simDeltaTrueLargerDeltaMin$eValuesStopped >=20),
               "expected"=935)
  expect_equal("object"=mean(simDeltaTrueLargerDeltaMin$stoppingTimes),
               "expected"=17.517)
})

