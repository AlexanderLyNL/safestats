
#' Print function for safe test designs for two proportions
#'
#' @param x Design for a safe test
#' @param ... Other arguments for the generic print function
#'
#' @return No return value
#' @export
#'
print.safe2x2_result <- function(x, ...){
  cat("\n")
  cat(paste("       ", "Safe Design for Test of Two Proportions", "\n"))
  cat("\n")
  if(x[["pilot"]] == FALSE){
    cat("Requires an experiment with sample sizes: ")
    cat("\n")
    cat(paste("    naPlan =", x[["na"]], "and nbPlan =", x[["nb"]]))
    cat("\n")
    cat("to find an effect size of at least: ")
    cat("\n")
    cat("    deltaMin =", round5(x[["deltaMin"]]))
    cat("\n")
    cat("with:")
    cat("\n")
    cat("    power = ", 1 - x[["beta"]], " (thus, beta = ", x[["beta"]], ")", sep="")
    cat("\n")
    cat("under the alternative:")
    cat("\n")
    cat("   ", getNameAlternative(x[["alternative"]], class(x)))
    cat("\n")
    cat("\n")
    cat("Based on the decision rule S > 1/alpha:")
    cat("\n")
    cat("    S > ", round5(1/x[["alpha"]]), sep="")
    cat("\n")
    cat("which occurs with chance less than:")
    cat("\n")
    cat("    alpha =", x[["alpha"]])
    cat("\n")
    cat("under iid Bernoulli distributed data with the same mean ( proportion) in group a and b.\n")
  } else {
    cat("The experiment is not planned.")
    cat("\n")
    cat("This design object only valid for experiments with:")
    cat("\n")
    cat("    na =", x[["na"]], "and nb =", x[["nb"]])
    cat("\n")
  }

}


#' Print function for objects of the \code{safe2x2_test} class.
#'
#' @param x Safe test result
#' @param ... Other arguments for the generic print function
#'
#' @return No return value
#' @export
#'
print.safe2x2_test <- function(x, ...){
  designObj <- x[["design"]]
  alternativeName <- getNameAlternative("alternative"=designObj[["alternative"]], "testType"=class(designObj))

  cat("\n")
  cat("    Safe test for 2x2 contingency tables \n")
  cat("\n")
  cat("data: \n")
  print(x[["data"]])

  if (designObj[["pilot"]]) {
    cat("The pilot test is based on an exploratory alpha =", designObj[["alpha"]])
    cat("\n")
    cat("and resulted in:  s-value =", round5(x[["s_value"]]))
    cat("\n")
    cat("Alternative hypothesis:")
  } else {
    cat("The test designed with alpha =", designObj[["alpha"]])
    cat("\n")
    cat("s-value =", round5(x[["s_value"]]), "> 1/alpha =", round5(1/designObj[["alpha"]]), ":",
        x[["s_value"]] > 1/designObj[["alpha"]])
    cat("\n")
    # Iets over n1Plan, n2Plan, etc
    cat("\n")
    cat(paste("Experiments required naPlan =", designObj[["na"]], "and nbPlan =",
              designObj[["nb"]], "samples,"))
    cat("\n")

    n1Diff <- designObj[["na"]] - x[["na"]]
    n2Diff <- designObj[["nb"]] - x[["nb"]]

    if (n1Diff > 0 || n2Diff > 0) {
      cat("    Note: ")
      if (n1Diff > 0) {
        cat("naPlan - na = ", n1Diff, ", ", sep="")
      }
      if (n2Diff > 0) {
        cat("nbPlan - nb =", n2Diff)
      }
      cat("\n")
    }

    cat("to guarantee a power of at least ", round5(1 - designObj[["beta"]]),
        " (beta =", round5(designObj[["beta"]]), "),", sep="")
    cat("\n")
    cat("under the alternative hypothesis: ")
  }
  cat("\n       ")
  cat(alternativeName)
  cat("\n")

  if (!designObj[["pilot"]]) {
    cat("and deltaMin =", designObj[["deltaMin"]])
  }
}
