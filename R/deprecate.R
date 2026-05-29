# Deprecate -------

#' @describeIn designSavi1aHelper Deprecated version of designSavi1aHelper
#' @usage NULL
#' @export
designSafe1aHelper <- function(
    samplingResult, esMin, beta, ratio,
    testType=c("oneSample", "paired","twoSample")) {

  warning('The function designSafe1aHelper is deprecated;',
          'Please use designSavi1aHelper instead')

  designSavi1aHelper(
    samplingResult=samplingResult, esMin=esMin,
    beta=beta, ratio=ratio, testType=testType)
}

#' @describeIn designSavi2Helper Deprecated version of designSavi2Helper
#' @usage NULL
#' @export
designSafe2Helper <- function(
    samplingResult, esMin, nPlan, ratio,
    testType=c("oneSample", "paired","twoSample")) {

  warning('The function designSafe2Helper is deprecated;',
          'Please use designSavi2Helper instead')

  designSavi2Helper(
    samplingResult=samplingResult, esMin=esMin,
    nPlan=nPlan, ratio=ratio, testType=testType)
}

#' @describeIn constructSaviDesignObj Deprecated version of constructSaviDesignObj
#' @usage NULL
#' @export
constructSafeDesignObj <- function(testName) {

  warning('The function constructSafeDesignObj is deprecated;',
          'Please use constructSaviDesignObj instead')

  constructSaviDesignObj(testName=testName)
}

#' @describeIn constructSaviTestObj Deprecated version of constructSaviTestObj
#' @usage NULL
#' @export
constructSafeTestObj <- function(testName) {
  warning('The function constructSafeTestObj is deprecated;',
          'Please use constructSaviTestObj instead')

  constructSaviTestObj(testName=testName)
}

#' @describeIn print.saviTest Deprecated version of print.saviTest
#' @usage NULL
#' @export
print.safeTest <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  print.saviTest(x=x, digits=digits, prefix=prefix, ...)
}

#' @describeIn print.saviDesign Deprecated version of print.saviDesign
#' @usage NULL
#' @export
print.safeDesign <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  print.saviDesign(x=x, digits=digits, prefix=prefix, ...)
}


#' @describeIn plot.saviDesign Deprecated version of plot.saviDesign
#' @usage NULL
#' @export
plot.safeDesign <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                            xlim=NULL, ylim=NULL, maxNBins=35,
                            numSamplePaths=100, wantStepLines=FALSE,
                            wantQuantiles=NULL, border="#1F78B4E6",
                            breaks=NULL, lwd=2, pch=15,
                            histInnerColour="#A6CEE380", col="#DAA52066",
                            colQuant="#AA0000", cex=1.3, ...) {
  plot.saviDesign(x=x, main=main, xlab=xlab, ylab=ylab,
                  xlim=xlim, ylim=ylim, maxNBins=maxNBins,
                  numSamplePaths=numSamplePaths, wantStepLines=wantStepLines,
                  wantQuantiles=wantQuantiles, border=border,
                  breaks=breaks, lwd=lwd, pch=pch,
                  histInnerColour=histInnerColour, col=col,
                  colQuant=colQuant, cex=cex, ...)
}

#' @describeIn plot.saviTest Deprecated version of plot.saviTest
#' @usage NULL
#' @export
plot.safeTest <- function(x, main=NULL, xlab=NULL, ylab=NULL,
                          xlim=NULL, ylim=NULL, lwd=3, cex=1.3,
                          fillPlot=NULL, switchNFill=1e4,
                          logScale=NULL, switchNLog=30,
                          h0Colour="gold", lineColour="#1F78B4E6",
                          col="#A6CEE380", border="#1F78B4E6",
                          wantConfSeqPlot=FALSE, add=FALSE,
                          density=NULL, angle=45,
                          fillOddEven=FALSE, ...) {

  plot.saviTest(x=x, main=main, xlab=xlab, ylab=ylab,
                xlim=xlim, ylim=ylim, lwd=lwd, cex=cex,
                fillPlot=fillPlot, switchNFill=switchNFill,
                logScale=logScale, switchNLog=switchNLog,
                h0Colour=h0Colour, lineColour=lineColour,
                col=col, border=border,
                wantConfSeqPlot=wantConfSeqPlot, add=add,
                density=density, angle=angle,
                fillOddEven=fillOddEven, ...)
}

# T-TEST -------

#' @describeIn saviTTest Deprecated version of saviTTestStat
#' @usage NULL
#' @export
safeTTestStat <- function(
    t, n1, n2=NULL, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    ...) {

  warning('The function safeTTestStat is deprecated;',
          'Please use saviTTestStat instead')

  saviTTestStat(t=t, n1=n1, n2=n2, parameter=parameter,
                alternative=alternative, tDensity=tDensity,
                paired=paired, eType=eType, ...)
}

#' @describeIn saviTTest Deprecated version of saviTTestStatNEffNu
#' @usage NULL
#' @export
safeTTestStatNEffNu <- function(
    t, nEff, nu, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    ...) {

  warning('The function safeTTestStatNEffNu is deprecated;',
          'Please use saviTTestStatNEffNu instead')

  saviTTestStat(t=t, nEff=nEff, nu=nu, parameter=parameter,
                alternative=alternative, tDensity=tDensity,
                paired=paired, eType=eType, ...)
}

#' @describeIn saviTTest Deprecated version of saviTTestStatNEffNuMom
#' @usage NULL
#' @export
safeTTestStatNEffNuMom <- function(
    t, nEff, nu, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE,
    k=1,
    ...) {

  warning('The function safeTTestStatNEffNuMom is deprecated;',
          'Please use saviTTestStatNEffNuMom instead')

  saviTTestStatNEffNuMom(t=t, nEff=nEff, nu=nu, parameter=parameter,
                alternative=alternative, tDensity=tDensity,
                paired=paired, k=1, ...)
}

#' @describeIn saviTTest Deprecated version of saviTTestStatTDensity
#' @usage NULL
#' @export
safeTTestStatTDensity <- function(t, parameter, nu, nEff,
                                  alternative=c("twoSided", "less", "greater"),
                                  paired=FALSE, ...) {

  warning('The function safeTTestStatTDensity is deprecated;',
          'Please use saviTTestStatTDensity instead')

  saviTTestStatTDensity(t=t, nEff=nEff, nu=nu, parameter=parameter,
                        alternative=alternative, paired=paired, ...)
}

#' @describeIn saviTTest Deprecated version of saviTTest
#' @usage NULL
#' @export
safeTTest <- function(x, ...) {

  warning('The function safeTTest is deprecated;',
          'Please use saviTTest instead')

  saviTTest(x, ...)
}

#' @describeIn saviTTest Deprecated version of saviTTest.default
#' @usage NULL
#' @export
safeTTest.default <- function(
    x, y=NULL, designObj=NULL, paired=FALSE,
    varEqual=TRUE, ciValue=NULL,
    maxRoot=10, sequential=NULL, ...) {

  warning('The function safeTTest.default is deprecated;',
          'Please use saviTTest.default instead')

  saviTTest.default(x=x, y=y, designObj=designObj, paired=paired,
                    varEqual=varEqual, ciValue=ciValue,
                    maxRoot=maxRoot, sequential=sequential, ...)
}

#' @describeIn saviTTest Deprecated version of saviTTest.formula
#' @usage NULL
#' @export
safeTTest.formula <- function(
    formula, data, subset, na.action, ...) {

  warning('The function safeTTest.formula is deprecated;',
          'Please use saviTTest.formula instead')

  saviTTest.formula(formula=formula, data=data, subset=subset,
                    na.action=na.action, ...)
}


#' @describeIn saviTTest Deprecated version of savi.t.test
#' @usage NULL
#' @export
safe.t.test <- function(x, y=NULL, paired=FALSE, designObj=NULL, varEqual=TRUE,
                        ciValue=NULL, ...) {

  warning('The function safe.t.test is deprecated;',
          'Please use savi.t.test instead')

  savi.t.test(x=x, y=y, designObj=designObj, paired=paired,
              varEqual=varEqual, ciValue=ciValue, ...)
}

#' @describeIn designSaviT Deprecated version of designSaviT
#' @usage NULL
#' @export
designSafeT <- function(
    deltaMin=NULL, beta=NULL, nPlan=NULL,
    alpha=0.05, h0=0, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    wantSamplePaths=TRUE,
    lowEsTrue=0.01, highEsTrue=3,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  warning('The function designSafeT is deprecated;',
          'Please use designSaviT instead')

  designSaviT(
    deltaMin=deltaMin, beta=beta, nPlan=nPlan,
    alpha=alpha, h0=h0, alternative=alternative,
    testType=testType, ratio=ratio, parameter=parameter,
    eType=eType, wantSamplePaths=wantSamplePaths,
    lowEsTrue=lowEsTrue, highEsTrue=highEsTrue,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @describeIn designSaviT1aWantNPlan Deprecated version of designSaviT1aWantNPlan
#' @usage NULL
#' @export
designSafeT1aWantNPlan <- function(
    deltaMin, beta, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  warning('The function designSafeT1aWantNPlan is deprecated;',
          'Please use designSaviT1aWantNPlan instead')

  designSaviT1aWantNPlan(
    deltaMin=deltaMin, beta=beta,
    alpha=alpha, alternative=alternative,
    testType=testType, ratio=ratio, parameter=parameter,
    eType=eType, wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @describeIn designSaviT Deprecated version of designSaviT2WantBeta
#' @usage NULL
#' @export
designSafeT2WantBeta <- function(
    deltaMin, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  warning('The function designSafeT2WantBeta is deprecated;',
          'Please use designSaviT2WantBeta instead')

  designSaviT2WantBeta(
    deltaMin=deltaMin, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    testType=testType, ratio=ratio, parameter=parameter,
    eType=eType, wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @describeIn designSaviT Deprecated version of designSafeT3WantEsMin
#' @usage NULL
#' @export
designSafeT3WantEsMin <- function(
    beta, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  warning('The function designSafeT3WantEsMin is deprecated;',
          'Please use designSaviT3WantEsMin instead')

  designSaviT3WantEsMin(
    beta=beta, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    testType=testType, parameter=parameter,
    eType=eType,
    lowEsTrue=lowEsTrue, highEsTrue=highEsTrue,
    ...)
}

#' @describeIn designSaviT Deprecated version of designSaviT3bWantParameter
#' @usage NULL
#' @export
designSafeT3bWantParameter <- function(
    nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    ...) {

  warning('The function designSafeT3bWantParameter is deprecated;',
          'Please use designSaviT3bWantParameter instead')

  designSaviT3bWantParameter(
    nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    testType=testType, parameter=parameter,
    eType=eType, ...)
}

#' @describeIn designSaviT Deprecated version of computeNPlanBatchSaviT
#' @usage NULL
#' @export
computeNPlanBatchSafeT <- function(
    deltaTrue, alpha=0.05, beta=0.2,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    parameter=NULL, ratio=1) {

  warning('The function computeNPlanBatchSafeT is deprecated;',
          'Please use computeNPlanBatchSaviT instead')

  computeNPlanBatchSaviT(
    deltaTrue=deltaTrue, beta=beta,
    alpha=alpha, alternative=alternative,
    testType=testType, ratio=ratio, parameter=parameter,
    eType=eType)
}

#' @describeIn designSaviT Deprecated version of computeMinEsBatchSaviT
#' @usage NULL
#' @export
computeMinEsBatchSafeT <- function(
    nPlan, alpha=0.05, beta=0.2,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  warning('The function computeMinEsBatchSafeT is deprecated;',
          'Please use computeMinEsBatchSaviT instead')

  computeMinEsBatchSaviT(
    beta=beta, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    testType=testType, ratio=ratio, parameter=parameter,
    eType=eType,
    lowEsTrue=lowEsTrue, highEsTrue=highEsTrue, ...)
}

#' @describeIn sampleStoppingTimesSaviT Deprecated version of sampleStoppingTimesSaviT
#' @usage NULL
#' @export
sampleStoppingTimesSafeT <- function(
    deltaTrue, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, lowN=3L, nMax=1e8L,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    wantEValuesAtNMax=FALSE,
    wantSamplePaths=TRUE, wantSimData=FALSE,
    pb=TRUE, seed=NULL, nSim=1e3L, ...) {

  warning('The function sampleStoppingTimesSafeT is deprecated;',
          'Please use sampleStoppingTimesSaviT instead')

  sampleStoppingTimesSaviT(
    deltaTrue=deltaTrue, beta=beta,
    alpha=alpha, alternative=alternative,
    testType=testType,
    ratio=ratio, parameter=parameter, lowN=lowN, nMax=nMax,
    eType=eType,
    wantEValuesAtNMax=wantEValuesAtNMax,
    wantSamplePaths=wantSamplePaths, wantSimData=wantSimData,
    pb=pb, seed=seed, nSim=nSim, ...)
}


#' @describeIn computePowerSaviT Deprecated version of computeBetaSaviT
#' @usage NULL
#' @export
computeBetaSafeT <- function(
    deltaTrue, nPlan, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {


  warning('The function computeBetaSafeT is deprecated;',
          'Please use computePowerSaviT instead')

  computePowerSaviT(
    deltaTrue=deltaTrue, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    testType=testType,
    parameter=parameter,
    eType=eType,
    wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @describeIn computeNPlanSaviT Deprecated version of computeNPlanSaviT
#' @usage NULL
#' @export
computeNPlanSafeT <- function(
    deltaTrue, beta=0.2, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, nMax=1e8,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai", "bayarri"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  warning('The function computeNPlanSafeT is deprecated;',
          'Please use computeNPlanSaviT instead')

  computeNPlanSaviT(
    deltaTrue=deltaTrue, beta=beta,
    alpha=alpha, alternative=alternative,
    testType=testType,
    ratio=ratio, parameter=parameter, nMax=nMax,
    eType=eType,
    wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

# Z-TEST ------
#' @rdname saviZTestStat
#' @aliases saviZTestStat
#' @usage NULL
#' @export
safeZTestStat <- function(
    z, n1, n2=NULL, parameter,
    alternative=c("twoSided", "less", "greater"),
    paired=FALSE, sigma=1,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    ...) {

  warning('The function safeZTestStat is deprecated;',
          'Please use saviZTestStat instead')

  saviZTestStat(z=z, n1=n1, n2=n2, parameter=parameter,
                alternative=alternative,
                paired=paired, sigma=sigma,
                eType=eType, ...)
}

#' @describeIn saviZTest Deprecated version of saviZTest
#' @usage NULL
#' @export
safeZTest <- function(x, ...) {

  warning('The function safeZTest is deprecated;',
          'Please use saviZTest instead')

  UseMethod("saviZTest")
}

#' @describeIn saviZTest Deprecated version of saviZTest.default
#' @usage NULL
#' @export
safeZTest.default <- function(
    x, y=NULL, paired=FALSE, designObj=NULL,
    ciValue=NULL, maxRoot=10, sequential=NULL, ...) {

  warning('The function safeZTest.default is deprecated;',
          'Please use saviZTest.default instead')

  saviZTest.default(x=x, y=y, paired=paired, designObj=designObj,
                    ciValue=ciValue, maxRoot=maxRoot,
                    sequential=sequential, ...)

}

#' @rdname saviZTest
#' @aliases saviZTest
#' @usage NULL
#' @export
safeZTest.formula <- function(
    formula, data, subset, na.action, ...) {

  warning('The function safeZTest.formula is deprecated;',
          'Please use saviZTest.formula instead')

  saviZTest.formula(formula=formula,
                    data=data, subset=subset,
                    na.action=na.action, ...)
}

#' @describeIn saviZTest Deprecated version of saviZTest
#' @usage NULL
#' @export
safe.z.test <- function(x, y=NULL, paired=FALSE,
                        designObj=NULL, ...) {

  warning('The function safe.z.test is deprecated;',
          'Please use savi.z.test instead')

  savi.z.test(x=x, y=y, paired=paired,
              designObj=designObj, ...)
}

#' @describeIn designSaviZ Deprecated version of designSaviZ
#' @usage NULL
#' @export
designSafeZ <- function(
    meanDiffMin=NULL, beta=NULL, nPlan=NULL,
    alpha=0.05, h0=0, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    lowEsTrue=0.01, highEsTrue=3,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    minEffiTest=FALSE, growFutility=FALSE, ...) {

  warning('The function designSafeZ is deprecated;',
          'Please use designSaviZ instead')

  designSaviZ(
    meanDiffMin=meanDiffMin, beta=beta, nPlan=nPlan,
    alpha=alpha, h0=h0, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    ratio=ratio, parameter=parameter,
    eType=eType,
    wantSamplePaths=wantSamplePaths,
    lowEsTrue=lowEsTrue, highEsTrue=highEsTrue,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot,
    minEffiTest=futility, growFutility=growFutility, ...)
}

#' @describeIn designSaviZ1aWantNPlan Deprecated version of designSaviZ1aWantNPlan
#' @usage NULL
#' @export
designSafeZ1aWantNPlan <- function(
    meanDiffMin, beta, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    minEffiTest=FALSE, growFutility=FALSE, ...) {

  warning('The function designSafeZ1aWantNPlan is deprecated;',
          'Please use designSaviZ1aWantNPlan instead')

  designSaviZ1aWantNPlan(
    meanDiffMin=meanDiffMin, beta=beta,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    ratio=ratio, parameter=parameter,
    eType=eType,
    wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot,
    minEffiTest=minEffiTest, growFutility=growFutility, ...)
}

#' @describeIn designSaviZ2WantPower Deprecated version of designSaviZ2WantBeta
#' @usage NULL
#' @export
designSafeZ2WantBeta <- function(
    meanDiffMin, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  warning('The function designSafeZ2WantBeta is deprecated;',
          'Please use designSaviZ2WantBeta instead')

  designSaviZ2WantPower(
    meanDiffMin=meanDiffMin, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    ratio=ratio, parameter=parameter,
    eType=eType,
    wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @describeIn designSaviZ3WantEsMin Deprecated version of designSaviZ3WantEsMin
#' @usage NULL
#' @export
designSafeZ3WantEsMin <- function(
    beta, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  warning('The function designSafeZ3WantEsMin is deprecated;',
          'Please use designSaviZ3WantEsMin instead')

  designSaviZ3WantEsMin(
    beta=beta, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    parameter=parameter,
    eType=eType,
    lowEsTrue=lowEsTrue, highEsTrue=highEsTrue, ...)
}

#' @describeIn designSaviZ3bWantParameter Deprecated version of designSaviZ3bWantParameter
#' @usage NULL
#' @export
designSafeZ3bWantParameter <- function(
    nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"), ...) {

  warning('The function designSafeZ3bWantParameter is deprecated;',
          'Please use designSaviZ3bWantParameter instead')

  designSaviZ3bWantParameter(
    nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    parameter=parameter,
    eType=eType,
    ...)
}

#' @describeIn computeNPlanBatchSaviZ Deprecated version of computeNPlanBatchSaviZ
#' @usage NULL
#' @export
computeNPlanBatchSafeZ <- function(
    meanDiffTrue, alpha=0.05, beta=0.2, sigma=1, kappa=sigma,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    parameter=NULL,
    highN=1e6, ratio=1, ...) {

  warning('The function computeNPlanBatchSafeZ is deprecated;',
          'Please use computeNPlanBatchSaviZ instead')

  computeNPlanBatchSaviZ(
    meanDiffTrue=meanDiffTrue, beta=beta,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    ratio=ratio, parameter=parameter,
    eType=eType, highN=highN, ...)
}

#' @describeIn computeBetaBatchSaviZ Deprecated version of computeBetaBatchSaviZ
#' @usage NULL
#' @export
computeBetaBatchSafeZ <- function(
    meanDiffTrue, nPlan, alpha=0.05, sigma=1, kappa=sigma,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow")) {

  warning('The function computeBetaBatchSafeZ is deprecated;',
          'Please use computeBetaBatchSaviZ instead')

  computeBetaBatchSaviZ(
    meanDiffTrue=meanDiffTrue, nPlan=nPlan,
    alpha=alpha, sigma=sigma, kappa=kappa,
    alternative=alternative,
    testType=testType,
    parameter=parameter,
    eType=eType)
}

#' @describeIn computeMinEsBatchSaviZ Deprecated version of computeMinEsBatchSaviZ
#' @usage NULL
#' @export
computeMinEsBatchSafeZ <- function(
    nPlan, alpha=0.05, beta=0.2, sigma=1, kappa=sigma,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    lowEsTrue=0.01, highEsTrue=0.002, ...) {

  warning('The function computeMinEsBatchSafeZ is deprecated;',
          'Please use computeMinEsBatchSaviZ instead')

  computeMinEsBatchSaviZ(
    beta=beta, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa, testType=testType,
    parameter=parameter,
    eType=eType,
    lowEsTrue=lowEsTrue, highEsTrue=highEsTrue, ...)
}

#' @describeIn sampleStoppingTimesSaviZ Deprecated version of sampleStoppingTimesSaviZ
#' @usage NULL
#' @export
sampleStoppingTimesSafeZ <- function(
    meanDiffTrue, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    sigma=1, kappa=sigma,
    ratio=1, parameter=NULL, nMax=1e3L,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantEValuesAtNMax=FALSE,
    wantSamplePaths=TRUE, wantSimData=FALSE,
    pb=TRUE, seed=NULL, nSim=1e3L, minEffiTest=FALSE,
    growFutility=FALSE, beta=NULL, ...) {

  warning('The function sampleStoppingTimesSafeZ is deprecated;',
          'Please use sampleStoppingTimesSaviZ instead')

  sampleStoppingTimesSaviZ(
    meanDiffTrue=meanDiffTrue, beta=beta,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa,
    ratio=ratio, parameter=parameter,
    eType=eType,
    wantEValuesAtNMax=wantEValuesAtNMax,
    wantSamplePaths=wantSamplePaths,
    wantSimData=wantSimData,
    pb=pb, seed=seed, nSim=nSim, minEffiTest=minEffiTest,
    growFutility=growFutility, ...)
}

#' @describeIn computePowerSaviZ Deprecated version of computeBetaSaviZ
#' @usage NULL
#' @export
computeBetaSafeZ <- function(
    meanDiffTrue, nPlan, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    sigma=1, kappa=sigma,
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {


  warning('The function computeBetaSafeZ is deprecated;',
          'Please use computePowerSaviZ instead')

  computePowerSaviZ(
    meanDiffTrue=meanDiffTrue, nPlan=nPlan,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa,
    testType=testType, parameter=parameter,
    eType=eType,
    wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @describeIn computeNPlanSaviZ Deprecated version of computeNPlanSaviZ
#' @usage NULL
#' @export
computeNPlanSafeZ <- function(
    meanDiffTrue, beta=0.2, alpha=0.05,
    alternative=c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    sigma=1, kappa=sigma,
    ratio=1, parameter=NULL, nMax=1e8,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    minEffiTest=FALSE, growFutility=FALSE,
    ...) {

  warning('The function computeNPlanSafeZ is deprecated;',
          'Please use computeNPlanSaviZ instead')

  computeNPlanSaviZ(
    meanDiffTrue=meanDiffTrue, beta=beta,
    alpha=alpha, alternative=alternative,
    sigma=sigma, kappa=kappa,
    testType=testType, parameter=parameter,
    eType=eType, ratio=ratio,
    nMax=nMax,
    wantSamplePaths=wantSamplePaths,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot,
    minEffiTest=minEffiTest, growFutility=growFutility, ...)
}

# LOGRANK ------

#' @rdname saviLogrankTest
#' @aliases saviLogrankTest
#' @usage NULL
#' @export
safeLogrankTest <- function(formula, designObj=NULL, ciValue=NULL, data=NULL, survTime=NULL,
                            group=NULL, pilot=FALSE, exact=TRUE, computeZ=TRUE, ...) {

  warning('The function safeLogrankTest is deprecated;',
          'Please use saviLogrankTest instead')

  saviLogrankTest(formula=formula, designObj=designObj, ciValue=ciValue,
                  data=data, survTime=survTime, group=group,
                  pilot=pilot, exact=exact, computeZ=computeZ, ...)

}

#' @rdname saviLogrankTest
#' @aliases saviLogrankTest
#' @usage NULL
#' @export
safeLogrankTestStat <- function(z, nEvents, designObj, ciValue=NULL,
                                dataNull=1, sigma=1) {

  warning('The function safeLogrankTestStat is deprecated;',
          'Please use saviLogrankTestStat instead')

  saviLogrankTestStat(z=z, nEvents=nEvents, designObj=designObj,
                      ciValue=ciValue, dataNull=dataNull, sigma=sigma)
}

#' @rdname designSaviLogrank
#' @aliases designSaviLogrank
#' @usage NULL
#' @export
designSafeLogrank <- function(
    hrMin=NULL, beta=NULL, nEvents=NULL,
    alpha=0.05, h0=1, alternative=c("twoSided", "greater", "less"),
    m0=50000L, m1=50000L,
    testType=c("exactLogrank", "gaussianLogrank"),
    ratio=1, exact=TRUE, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    groupSizePerTimeFunction=returnOne,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {

  warning('The function designSafeLogrank is deprecated;',
          'Please use designSaviLogrank instead')

  designSaviLogrank(
    hrMin=hrMin, beta=beta, nEvents=nEvents,
    alpha=alpha, h0=h0, alternative=alternative,
    m0=m0, m1=m1,
    testType=testType,
    ratio=ratio, exact=exact, parameter=parameter,
    eType=eType, wantSamplePaths=wantSamplePaths,
    groupSizePerTimeFunction=groupSizePerTimeFunction,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}

#' @rdname designSaviLogrank2WantBeta
#' @aliases designSaviLogrank2WantBeta
#' @usage NULL
#' @export
designSafeLogrank2WantBeta <- function(
    hrMin, nEvents,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    m0=50000L, m1=50000L,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow"),
    wantSamplePaths=TRUE,
    groupSizePerTimeFunction=returnOne,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim, ...) {


  warning('The function designSafeLogrank2WantBeta is deprecated;',
          'Please use designSaviLogrank2WantBeta instead')

  designSaviLogrank2WantBeta(
    hrMin=hrMin, nEvents=nEvents,
    alpha=alpha, alternative=alternative,
    m0=m0, m1=m1,
    testType=testType,
    ratio=ratio, parameter=parameter,
    eType=eType, wantSamplePaths=wantSamplePaths,
    groupSizePerTimeFunction=groupSizePerTimeFunction,
    pb=pb, seed=seed, nSim=nSim, nBoot=nBoot, ...)
}


# 2x2 ------
#' @rdname designSaviTwoProportions
#' @aliases designSaviTwoProportions
#' @usage NULL
#' @export
designSafeTwoProportions <- function(na, nb,
                                     nBlocksPlan = NULL,
                                     beta = NULL,
                                     delta = NULL,
                                     alternativeRestriction = c("none", "difference", "logOddsRatio"),
                                     alpha = 0.05,
                                     pilot = "FALSE",
                                     hyperParameterValues = NULL,
                                     previousSafeTestResult = NULL,
                                     M = 1e3,
                                     simThetaAMin = NULL,
                                     simThetaAMax = NULL) {

  warning('The function designSafeTwoProportions is deprecated;',
          'Please use designSaviTwoProportions instead')

  designSaviTwoProportions(
    na=na, nb=nb, nBlocksPlan=nBlocksPlan,
    beta=beta, delta=delta, alternativeRestriction=alternativeRestriction,
    alpha=alpha, pilot=pilot, hyperParameterValues=hyperParameterValues,
    previousSaviTestResult=previousSafeTestResult, M=M,
    simThetaAMin=simThetaAMin, simThetaAMax=simThetaAMax)
}

#' @rdname saviTwoProportionsTest
#' @aliases saviTwoProportionsTest
#' @usage NULL
#' @export
safeTwoProportionsTest <- function(ya, yb, designObj = NULL, wantConfidenceSequence = FALSE, ciValue = NULL,
                                   confidenceBoundGridPrecision = 20, logOddsConfidenceSearchBounds = c(0.01, 5), pilot = FALSE) {

  warning('The function safeTwoProportionsTest is deprecated;',
          'Please use saviTwoProportionsTest instead')

  saviTwoProportionsTest(
    ya=ya, yb=yb, designObj=designObj,
    wantConfidenceSequence=wantConfidenceSequence,
    ciValue=ciValue, confidenceBoundGridPrecision=confidenceBoundGridPrecision,
    logOddsConfidenceSearchBounds=logOddsConfidenceSearchBounds,
    pilot=pilot)
}

#' @rdname saviTwoProportionsTest
#' @aliases savi.prop.test
#' @usage NULL
#' @export
safe.prop.test <- function(ya, yb, designObj = NULL, wantConfidenceSequence = FALSE, ciValue = NULL,
                           confidenceBoundGridPrecision = 20, logOddsConfidenceSearchBounds = c(0.01, 5), pilot = FALSE) {

  warning('The function safe.prop.test is deprecated;',
          'Please use saviTwoProportionsTest instead')

  saviTwoProportionsTest(
    ya=ya, yb=yb, designObj=designObj, wantConfidenceSequence=wantConfidenceSequence,
    ciValue=ciValue, confidenceBoundGridPrecision=confidenceBoundGridPrecision,
    logOddsConfidenceSearchBounds=logOddsConfidenceSearchBounds,
    pilot=pilot)
}

#' @rdname print.savi2x2Sim
#' @aliases print.savi2x2Sim
#' @usage NULL
#' @export
print.safe2x2Sim <- function(x, ...) {

  warning('The function print.safe2x2Sim is deprecated;',
          'Please use print.savi2x2Sim instead')

  print.savi2x2Sim(x, ...)
}


#' @rdname plot.savi2x2Sim
#' @aliases plot.savi2x2Sim
#' @usage NULL
#' @export
plot.safe2x2Sim <- function(x, ...) {

  warning('The function print.safe2x2Sim is deprecated;',
          'Please use print.savi2x2Sim instead')

  plot.savi2x2Sim(x, ...)
}
