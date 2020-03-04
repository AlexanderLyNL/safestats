# Safestats: an R package for safe, anytime-valid inference.

```Safestats``` is for designing and applying **safe** hypothesis tests. It can be used for designing hypothesis tests in the prospective or randomised controlled trial (RCT) setting, where the tests can be used under optional stopping and experiments often can be stopped early. The 'pilot' functions in
the package also enable using the safe tests in the observational/ retrospective
setting. For examples and explanation about which test to choose for which study setup, our ```vignette``` can be used (see below). The current version includes safe t-tests and tests of two proportions. The initial paper on the theory of safe testing and a worked-out example for the t-test can be found [in this paper](https://arxiv.org/abs/1906.07801). More on the theory behind the development of the safe tests for proportions can be found [here](https://www.universiteitleiden.nl/binaries/content/assets/science/mi/scripties/statscience/2019-2020/thesis_rjturner_for_publication.pdf).
## Installation
The development version can be found on [GitHub](https://github.com/alexanderlynl/safestats), which can be installed with the `devtools` package from [CRAN](https://cran.r-project.org/package=devtools) by entering in `R`:

```{r devtools, eval=FALSE}
devtools::install_github("AlexanderLyNL/safestats", build_vignettes = TRUE)
```
Once a stable version is finished, you can install it from CRAN in R through:
```
install.packages("safestats")
```

## Usage
Load the safestats package:
```
library("safestats")
```
For an extensive overview of all functions and their usage, we recommend following the safestats vignette tutorial:
```
browseVignettes("safestats")
```
