
#' Print function for safe test designs for two proportions
#'
#' @param x Design for a safe test
#' @param ... Other arguments for the generic print function
#'
#' @return No return value
#' @export
#'
print.safe2x2_result <- function(x, ...){
  cat("Call: \n")
  cat(deparse(x[["call"]]))
  cat("\n"); cat("\n")
  cat("Sample sizes for this design:\na = ")
  cat(x[["na"]]); cat(" and b = ");
  cat(x[["nb"]]); cat("\n");
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
  cat("\n")
  cat("    Safe test for 2x2 contingency tables \n")
  cat("\n")
  cat("data: \n")
  print(x[["data"]])
  cat("\n")
  cat("s-value: "); cat(x[["s_value"]])
  cat("\n")
  cat("safe p-value: "); cat(min(1,x[["safe_p_value"]]))
}
