# safestats 0.8.9
* Added:
  - runningIntersection argument added to display the resulting running
  intersection confidence sequence if needed
  - added savi 2x2 test conditional
  - added savi one-way anova test functions (initial versions)
* Bug fixes:
  - Removed some html code from the help files 
  - Cleaned up vignette code
  
# safestats 0.8.8
* Added:
  - Generic plot functions for saviDesign and saviTest objects. 
  - Cleaned up print functions for saviDesigns.
  - Added relevance testing for z- and t-tests.
  - References in help files.
  - Possibility to design z and t-tests based on nPlan only.
* Bug fixes:
  - Deprecated "safeDesign"" and "safeTest" function identifiers and 
  changed to "saviDesign" and "saviTest" functions instead. 

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
