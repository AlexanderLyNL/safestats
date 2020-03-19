## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 4,
  fig.width = 8
)

## ----install, eval=FALSE------------------------------------------------------
#  install.packages("safestats")

## ----devtools, eval=FALSE-----------------------------------------------------
#  devtools::install_github("AlexanderLyNL/safestats", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
library(safestats)

## -----------------------------------------------------------------------------
alpha <- 0.05
beta <- 0.2
deltaMin <- 9/(sqrt(2)*15)

designObj <- designSafeT(deltaMin=deltaMin, alpha=alpha, beta=beta,
                         alternative="greater", testType="pairedSampleT")
designObj

## ---- fig.height=4, fig.width=8-----------------------------------------------
# Recall:
# alpha <- 0.05
# beta <- 0.2

result <- plotSafeTDesignSampleSizeProfile(alpha=alpha, beta=beta, 
                                           maxN=100, 
                                           testType="pairedSampleT")

## -----------------------------------------------------------------------------
set.seed(1)
preData <- rnorm(n=designObj$n1Plan, mean=120, sd=15)
postData <- rnorm(n=designObj$n2Plan, mean=120, sd=15)
# Thus, the true delta is 0:
# deltaTrue <- (120-120)/(sqrt(2)*15)     

## -----------------------------------------------------------------------------
safeTTest(x=preData, y=postData, alternative = "greater",
          designObj=designObj, paired=TRUE)

## -----------------------------------------------------------------------------
safe.t.test(x=preData, y=postData, alternative = "greater",
            designObj=designObj, paired=TRUE)

## -----------------------------------------------------------------------------
# alpha <- 0.05

set.seed(1)
sValues <- replicate(n=1000, expr={
  preData <- rnorm(n=designObj[["n1Plan"]], mean=120, sd=15)
  postData <- rnorm(n=designObj[["n2Plan"]], mean=120, sd=15)
  safeTTest(x=preData, y=postData, alternative = "greater",
            designObj=designObj,paired=TRUE)$sValue}
)

mean(sValues > 20)
mean(sValues > 20) < alpha

## -----------------------------------------------------------------------------
set.seed(1)
preData <- rnorm(n=designObj[["n1Plan"]], mean=120, sd=15)
postData <- rnorm(n=designObj[["n2Plan"]], mean=111, sd=15)

safeTTest(x=preData, y=postData, alternative = "greater",
          designObj=designObj, paired=TRUE)

## -----------------------------------------------------------------------------
# Recall:
# alpha <- 0.05
# beta <- 0.2
power <- 1-beta

set.seed(1)
sValues <- replicate(n=1000, expr={
  preData <- rnorm(n=designObj[["n1Plan"]], mean=120, sd=15)
  postData <- rnorm(n=designObj[["n2Plan"]], mean=111, sd=15)
  safeTTest(x=preData, y=postData, alternative = "greater", designObj=designObj,
            paired=TRUE)$sValue})
mean(sValues > 1/alpha)
mean(sValues > 1/alpha) >= power

## -----------------------------------------------------------------------------
designObj

## -----------------------------------------------------------------------------
# Recall:
# alpha <- 0.05
# beta <- 0.2
# deltaMin <- 9/(sqrt(2)*15)      # = 0.42
simResultDeltaTrueIsDeltaMin <- simulate(object=designObj, nsim=1000L,
                                         seed=1, deltaTrue=deltaMin,
                                         muGlobal=120, sigmaTrue=15)
simResultDeltaTrueIsDeltaMin

## -----------------------------------------------------------------------------
plot(simResultDeltaTrueIsDeltaMin)

## -----------------------------------------------------------------------------
plot(simResultDeltaTrueIsDeltaMin, showOnlyNRejected=TRUE)

## -----------------------------------------------------------------------------
# Recall:
# alpha <- 0.05
# beta <- 0.2
# deltaMin <- 9/(sqrt(2)*15)      # = 0.42
deltaTrueLarger <- 0.6

simResultDeltaTrueLargerThanDeltaMin <- simulate(object=designObj,
                                                 nsim=1000L, seed=1,
                                                 deltaTrue=deltaTrueLarger,
                                                 muGlobal=120, sigmaTrue=15)
simResultDeltaTrueLargerThanDeltaMin

## -----------------------------------------------------------------------------
plot(simResultDeltaTrueLargerThanDeltaMin)

## -----------------------------------------------------------------------------
# Recall:
# alpha <- 0.05
# beta <- 0.2
# deltaMin <- 9/(sqrt(2)*15)      # = 0.42

freqDesignObj <- designFreqT(deltaMin=deltaMin, alpha=alpha, beta=beta,
                             alternative="greater", testType="pairedSampleT")

simResultDeltaTrueIsZero <- simulate(object=designObj, nsim=1000L, seed=1,
                                     deltaTrue=0, freqOptioStop=TRUE,
                                     n1PlanFreq=freqDesignObj$n1PlanFreq,
                                     n2PlanFreq=freqDesignObj$n2PlanFreq,
                                     muGlobal=120, sigmaTrue=15)
simResultDeltaTrueIsZero

## -----------------------------------------------------------------------------
dataBatch1 <- generateTTestData(n1Plan=freqDesignObj$n1PlanFreq,
                               n2Plan=freqDesignObj$n2PlanFreq, 
                               deltaTrue=0, nsim=1000, paired=TRUE, seed=1,
                               muGlobal=120, sigmaTrue=15)

pValuesBatch1 <- vector("numeric", length=1000)

for (i in seq_along(pValuesBatch1)) {
  pValuesBatch1[i] <- t.test(x=dataBatch1$dataGroup1[i, ], 
                             y=dataBatch1$dataGroup2[i, ], 
                             alternative="greater", paired=TRUE)$p.value
}
mean(pValuesBatch1 > alpha)
sum(pValuesBatch1 < alpha)

## -----------------------------------------------------------------------------
selectivelyContinueDeltaTrueIsZeroWithP <-
  selectivelyContinueTTestCombineData(oldValues=pValuesBatch1,
                                      valuesType="pValues", 
                                      alternative="greater", 
                                      oldData=dataBatch1,
                                      deltaTrue=0,
                                      n1Extra=freqDesignObj$n1PlanFreq,
                                      n2Extra=freqDesignObj$n2PlanFreq,
                                      alpha=alpha,
                                      seed=2, paired=TRUE,
                                      muGlobal=120, sigmaTrue=15)

## -----------------------------------------------------------------------------
pValuesBatch1To2 <- selectivelyContinueDeltaTrueIsZeroWithP$newValues
sum(pValuesBatch1To2 < alpha)

## -----------------------------------------------------------------------------
dataBatch1 <- list(dataGroup1=simResultDeltaTrueIsZero$safeSim$dataGroup1,
                   dataGroup2=simResultDeltaTrueIsZero$safeSim$dataGroup2)

sValuesBatch1 <- simResultDeltaTrueIsZero$safeSim$sValues
sum(sValuesBatch1 > 1/alpha)

## -----------------------------------------------------------------------------
selectivelyContinueDeltaTrueIsZero <- 
  selectivelyContinueTTestCombineData(oldValues=sValuesBatch1,
                                      designObj=designObj,
                                      alternative="greater", 
                                      oldData=dataBatch1,
                                      deltaTrue=0,
                                      seed=2, paired=TRUE,
                                      muGlobal=120, sigmaTrue=15,
                                      moreMainText="Batch 1-2")

## -----------------------------------------------------------------------------
sValuesBatch1To2 <- selectivelyContinueDeltaTrueIsZero$newValues
sum(sValuesBatch1To2 > 1/alpha)
length(sValuesBatch1To2)

## -----------------------------------------------------------------------------
for (j in 1:3) {
  oldSValues <- selectivelyContinueDeltaTrueIsZero$newValues
  oldData <- selectivelyContinueDeltaTrueIsZero$combinedData
  
  selectivelyContinueDeltaTrueIsZero <- 
  selectivelyContinueTTestCombineData(oldValues=oldSValues,
                                      designObj=designObj,
                                      alternative="greater", 
                                      oldData=oldData,
                                      deltaTrue=0,
                                      seed=2+j, paired=TRUE, 
                                      muGlobal=120, sigmaTrue=15,
                                      moreMainText=paste("Batch: 1 to", j+2))
  print(paste("Batch: 1 to", j+2))
  print(paste("Number of rejections:",
              sum(selectivelyContinueDeltaTrueIsZero$newValues > 1/alpha)))
}

## -----------------------------------------------------------------------------
dataBatch1 <- list(
  dataGroup1=simResultDeltaTrueIsDeltaMin$safeSim$dataGroup1,
  dataGroup2=simResultDeltaTrueIsDeltaMin$safeSim$dataGroup2
)

sValuesBatch1 <- simResultDeltaTrueIsDeltaMin$safeSim$sValues
sum(sValuesBatch1 > 1/alpha)

## -----------------------------------------------------------------------------
selectivelyContinueDeltaTrueIsDeltaMin <- 
  selectivelyContinueTTestCombineData(oldValues=sValuesBatch1,
                                      designObj=designObj,
                                      alternative="greater", 
                                      oldData=dataBatch1,
                                      deltaTrue=deltaMin,
                                      seed=2, paired=TRUE, muGlobal=120,
                                      sigmaTrue=15)

## -----------------------------------------------------------------------------
sValuesBatch1To2 <- selectivelyContinueDeltaTrueIsDeltaMin$newValues
sum(sValuesBatch1To2 > 1/alpha)

## -----------------------------------------------------------------------------
sValuesOri <- simResultDeltaTrueIsZero$safeSim$sValues

## -----------------------------------------------------------------------------
# Needs to be larger than 1/designObj$n1Plan to have at least two samples 
# in the replication attempt
someConstant <- 1.2

repData <- generateTTestData(n1Plan=ceiling(someConstant*designObj$n1Plan),
                             n2Plan=ceiling(someConstant*designObj$n2Plan), 
                             deltaTrue=0, nsim=1000, 
                             muGlobal=90, sigmaTrue=6,
                             paired=TRUE, seed=2)

sValuesRep <- vector("numeric", length=1000)

for (i in seq_along(sValuesRep)) {
  sValuesRep[i] <- safeTTest(x=repData$dataGroup1[i, ], 
                          y=repData$dataGroup2[i, ], 
                          designObj=designObj,
                          alternative="greater", paired=TRUE)$sValue
}
sValuesMultiply <- sValuesOri*sValuesRep
mean(sValuesMultiply > 1/alpha)

## -----------------------------------------------------------------------------
sValuesOri <- simResultDeltaTrueIsDeltaMin$safeSim$sValues

## -----------------------------------------------------------------------------
# Needs to be larger than 1/designObj$n1Plan to have at least two samples 
# in the replication attempt
someConstant <- 1.2

repData <- generateTTestData(n1Plan=ceiling(someConstant*designObj$n1Plan),
                             n2Plan=ceiling(someConstant*designObj$n2Plan), 
                             deltaTrue=deltaMin, nsim=1000, 
                             muGlobal=110, sigmaTrue=50,
                             paired=TRUE, seed=2)

sValuesRep <- vector("numeric", length=1000)

for (i in seq_along(sValuesRep)) {
  sValuesRep[i] <- safeTTest(x=repData$dataGroup1[i, ], 
                          y=repData$dataGroup2[i, ], 
                          designObj=designObj,
                          alternative="greater", paired=TRUE)$sValue
}
sValuesMulti <- sValuesOri*sValuesRep
mean(sValuesMulti > 1/alpha)

## -----------------------------------------------------------------------------
safeDesignProportions <- designSafeTwoProportions(deltaMin=0.3, alpha=0.05,
                                                  beta=0.20, lowN=100,
                                                  numberForSeed = 5227)

## -----------------------------------------------------------------------------
safeDesignProportions$n.star

## -----------------------------------------------------------------------------
sampleExample <- as.table(matrix(c(10, safeDesignProportions[["na"]]-10, 40,
                                   safeDesignProportions[["nb"]]-40), 
                                 byrow=TRUE, nrow=2))
colnames(sampleExample) <- c(0, 1)
sampleExample

## -----------------------------------------------------------------------------
safeTwoProportionsTest(x = sampleExample, testDesign = safeDesignProportions)

## -----------------------------------------------------------------------------
plotResult <- plotSafeTwoProportionsSampleSizeProfile(alpha=0.05,
                                                      beta=0.20,
                                                      highN=200, 
                                                      maxN=100,
                                                      numberForSeed=5222)

## -----------------------------------------------------------------------------
safePilotDesign <- designPilotSafeTwoProportions(na = 50, nb = 50)

## -----------------------------------------------------------------------------
set.seed(5224)

optionalStoppingTrueMeanIsDesign <- 
  simulateSpreadSampleSizeTwoProportions(
    safeDesign=safeDesignProportions, M=1000,
    parametersDataGeneratingDistribution=c(0.3, 0.6))

plotHistogramDistributionStoppingTimes(
  optionalStoppingTrueMeanIsDesign, 
  nPlan=safeDesignProportions[["n.star"]], 
  deltaTrue = 0.3)

## -----------------------------------------------------------------------------
#power achieved:
mean(optionalStoppingTrueMeanIsDesign$rejected == 1)

## -----------------------------------------------------------------------------
set.seed(5224)

optionalStoppingTrueDifferenceBig <- 
  simulateSpreadSampleSizeTwoProportions(
    safeDesign=safeDesignProportions, M=1000, 
    parametersDataGeneratingDistribution = c(0.2, 0.9))

plotHistogramDistributionStoppingTimes(
  optionalStoppingTrueDifferenceBig, nPlan=safeDesignProportions[["n.star"]],
  deltaTrue = 0.7)

## -----------------------------------------------------------------------------
#power achieved:
mean(optionalStoppingTrueDifferenceBig$rejected == 1)

## -----------------------------------------------------------------------------
set.seed(5224)

optionalStoppingTrueMeanNull <- 
  simulateSpreadSampleSizeTwoProportions(
    safeDesign=safeDesignProportions, M=1000, 
    parametersDataGeneratingDistribution = c(0.5, 0.5))

plotHistogramDistributionStoppingTimes(
  optionalStoppingTrueMeanNull, 
  nPlan=safeDesignProportions[["n.star"]], 
  deltaTrue = 0)

## -----------------------------------------------------------------------------
# The rate of false null rejections remained under alpha=0.05
mean(optionalStoppingTrueMeanNull$rejected == 1)

## -----------------------------------------------------------------------------
set.seed(5224)

fisher_result <- simulateFisherSpreadSampleSizeOptionalStopping(
  deltaDesign=0.5, alpha=0.05, nDesign=safeDesignProportions$n.star, 
  power=0.8, M=100, parametersDataGeneratingDistribution=c(0.5, 0.5))

mean(fisher_result$rejected == 1)

## -----------------------------------------------------------------------------
notRejectedIndex <- which(optionalStoppingTrueMeanIsDesign$rejected==FALSE)
sValuesNotRejected <- optionalStoppingTrueMeanIsDesign$s_values[notRejectedIndex]
nullNotRejectedIndex <- which(optionalStoppingTrueMeanNull$rejected == FALSE)
sValuesNotRejectedNull <- optionalStoppingTrueMeanNull$s_values[nullNotRejectedIndex]

## ---- echo = FALSE------------------------------------------------------------
trueHist <- graphics::hist(x = sValuesNotRejected, plot = FALSE)
nullHist <- graphics::hist(x = sValuesNotRejectedNull, plot = FALSE)
yMax <- max(trueHist[["counts"]], nullHist[["counts"]])
graphics::par(cex.main=1.5, mar=c(5, 6, 4, 4)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
              font.lab=2, cex.axis=1.3, bty="n", las=1)
graphics::plot(nullHist, xlim = c(0, max(sValuesNotRejected, sValuesNotRejectedNull)), 
               freq = FALSE, col = "blue", density = 20, angle = 45, xlab = "s-values", 
               main = "Histogram of s-values where null not rejected")
graphics::plot(trueHist, add = TRUE, freq = FALSE, col = "red", density = 20, 
               angle = -45)
graphics::legend(x = "topright", legend = c("True delta: null", "True delta: design"), fill = c("blue", "red"))

## ----optionalContinuation2x2--------------------------------------------------
continueIndex <- which(optionalStoppingTrueMeanIsDesign$s_values < 20 & 
                         optionalStoppingTrueMeanIsDesign$s_values > 10)

interestingSValues <-
  optionalStoppingTrueMeanIsDesign$s_values[continueIndex]

newSValues <- 
  simulateOptionalContinuationTwoProportions(
    interestingSValues, nFollowUp=40, 
    parametersDataGeneratingDistribution=c(0.3, 0.6))

mean(newSValues>=20)

## ----optionalContinuation2x2Null----------------------------------------------
continueIndex <- optionalStoppingTrueMeanNull$s_values < 20 & 
  optionalStoppingTrueMeanNull$s_values > 1

interestingSValues <-optionalStoppingTrueMeanNull$s_values[continueIndex]

newSValues <- 
  simulateOptionalContinuationTwoProportions(
    interestingSValues, nFollowUp=40, 
    parametersDataGeneratingDistribution=c(0.5, 0.5))

mean(newSValues>=20)

## -----------------------------------------------------------------------------
safeDesignProportionsOneSided <- 
  designSafeTwoProportions(deltaMin=0.5, alternative="greater",
                           numberForSeed = 291202)

## -----------------------------------------------------------------------------
sampleExampleGreater <- 
  as.table(matrix(c(5, safeDesignProportionsOneSided[["na"]]-5, 19,
                    safeDesignProportionsOneSided[["nb"]]-19), 
                  byrow=TRUE, nrow=2))

colnames(sampleExampleGreater) <- c(0,1)
sampleExampleGreater

## -----------------------------------------------------------------------------
safeTwoProportionsTest(x=sampleExampleGreater, 
                       testDesign=safeDesignProportionsOneSided)

## -----------------------------------------------------------------------------
sampleExampleLesser <- 
  as.table(matrix(c(safeDesignProportionsOneSided[["na"]]-5, 5,
                    safeDesignProportionsOneSided[["nb"]]-19, 19), 
                  byrow=TRUE, nrow=2))

colnames(sampleExampleGreater) <- colnames(sampleExampleLesser) <- c(0,1)
sampleExampleLesser

## -----------------------------------------------------------------------------
safeTwoProportionsTest(x=sampleExampleLesser,
                       testDesign=safeDesignProportionsOneSided)

## -----------------------------------------------------------------------------
safeDesignProportionsImbalanced <- 
  designSafeTwoProportions(deltaMin=0.3, alpha=0.05, beta=0.20, lowN=120,
                           sampleSizeRatio=2)
safeDesignProportionsImbalanced

