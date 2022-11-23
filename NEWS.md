# safestats 0.8.7
* Added:
  - safeDesign objects for z-, t-, and logrank tests now show nMean
* Bug fixes:
  - Removed browser statement in designSafeLogrank
  - Deprecated alternative="two.sided" in favour of alternative="twoSided"
  - checkAndReturnsNPlan now correctly show nPlan for testType="paired"

# safestats 0.8.6
* First release:
  - safe z-tests
	- safe t-tests
	- safe test for two proportions
  - safe logrank tests
