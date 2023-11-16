bob <- designSafeZ(0.5, beta=0.2, eType="mom")
bobLarger <- designSafeZ(0.6, beta=0.2, eType="mom",
                         parameter=bob$parameter)
bobSmaller <- designSafeZ(0.1, beta=0.2, eType="mom",
                          parameter=bob$parameter)

bobSmallerFromTheStart <- designSafeZ(0.1, beta=0.2, eType="mom")

brie <- designSafeZ(0.5, beta=0.2, eType="eGauss")
brieLarger <- designSafeZ(0.6, beta=0.2, eType="eGauss",
                          parameter=brie$parameter)

brieSmaller <- designSafeZ(0.1, beta=0.2, eType="eGauss",
                           parameter=brie$parameter)
brieSmallerFromTheStart <- designSafeZ(0.1, beta=0.2, eType="eGauss")


bob$nPlan
bobLarger$nPlan

brie$nPlan
brieLarger$nPlan

bob3$nPlan
bob4$nPlan

bob2$nPlan
brie2$nPlan

brie3$nPlan
brie4$nPlan


alternative=alternative,
testType=testType, parameter=bob$parameter,

kaas <- sampleStoppingTimesSafeZ(meanDiffTrue, alpha=alpha,
                                 alternative=alternative,
                                 testType=testType, parameter=bob$parameter,
                                 nMax=bob$nPlan)


a1 <- safeZTestStat(1, n1=23, parameter=0.4, eType="imom",
                    alternative="less")$eValue
a2 <- safeZTestStat(1, n1=23, parameter=0.4, eType="imom",
                    alternative="greater")$eValue
b1 <- safeZTestStat(1, n1=23, parameter=0.4, eType="imom",
                    alternative="twoSided")$eValue
(a1+a2)
2*b1
