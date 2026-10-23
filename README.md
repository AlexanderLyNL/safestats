# Safestats: an R package for safe, anytime-valid inference.

```Safestats``` is for designing and applying **safe** hypothesis tests. It can be used for designing hypothesis tests in the prospective or randomised controlled trial (RCT) setting, where the tests can be used under optional stopping allowing experiments to be stopped early without the loss of statistical validity. The 'pilot' functions in the package also enable using the safe tests in the observational/ retrospective setting. For examples and explanation about which test to choose for which study setup, our ```vignette``` can be used (see below). The current version includes safe t-tests and tests of two proportions. The initial paper on the theory of safe testing and a worked-out example for the t-test can be found [in this paper](https://arxiv.org/abs/1906.07801). More on the theory behind the development of the safe tests for proportions can be found [here](https://arxiv.org/abs/2106.02693).
## Installation
The main package can be found on CRAN. To install run the following in R:
```
install.packages("safestats")
```
The development version can be found on [GitHub]([https://github.com/alexanderlynl/safestats](https://github.com/AlexanderLyNL/safestats/tree/development)), which can be installed with the `pak` package from [CRAN](https://cran.r-project.org/package=pak) by entering in `R`:

```{r pak, eval=FALSE}
pak::pak("AlexanderLyNL/safestats@development")
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
