#' Perform a safe test for two proportions
#'
#' Perform a safe test for two proportions (a 2x2 contingency table test) with a
#' result object retrieved through one of design functions for two proportions
#' in this package, \code{\link{designPilotSafeTwoProportions}} or
#' \code{\link{designSafeTwoProportions}}.
#'
#' @param x either a two-dimensional contingency table in matrix form,
#' or a factor object indicating the groups (equivalent to \code{fisher.test}).
#' @param y a factor object indicating the observations; ignored if x is a matrix.
#' default \code{NULL}.
#' @param testDesign  a safe test design for two proportions retrieved through \code{\link{designPilotSafeTwoProportions}} or
#' \code{\link{designSafeTwoProportions}}.
#'
#' @return safe2x2 test result object with the data used for the test, the S-value and the corresponding
#' safe p-value.
#'
#' @export
#'
#' @examples
#' testDesign <- designPilotSafeTwoProportions(na = 10, nb = 10)
#' #test with two factor ojects
#' groups <- factor(c(rep("a", 10), rep("b", 10)))
#' observed_data <- factor(c(rep("a", 7), rep("b", 3),rep("a", 3), rep("b", 7)))
#' safeTwoProportionsTest(x = groups, y = observed_data, testDesign = testDesign)
#'
#' #test with a table
#' dataTable <- table(groups, observed_data)
#' safeTwoProportionsTest(x = dataTable, testDesign = testDesign)
safeTwoProportionsTest <- function(x, y = NULL, testDesign) {
  if (!(is.factor(x) & is.factor(y)) & !is.table(x)) {
    stop('provide a matrix or two factors representing the observed data,
         see ?S_test')
  }

  if (!all(c("w1", "point_h0", "H1set", "na", "nb") %in% names(testDesign))) {
    stop(
      "passed safe2x2 design object is incomplete:
         provide parameters, weights and group sizes"
    )
  }

  #retrieve designed model
  w1 <- testDesign[["w1"]]
  point_h0 <- testDesign[["point_h0"]]
  H1set <- testDesign[["H1set"]]
  na <- testDesign[["na"]]
  nb <- testDesign[["nb"]]

  #retrieve observed data from table
  if (is.factor(x)) {
    d.table <- table(x, y)
  } else {
    d.table <- x
  }
  na1 <- d.table[1, 2]
  nb1 <- d.table[2, 2]

  #check if group size observed in table matches group size S-test was designed for
  #maybe insert a warning here later
  if (all.equal(rowSums(d.table), c(na, nb), check.attributes = FALSE) != TRUE) {
    na <- rowSums(d.table)[1]
    nb <- rowSums(d.table)[2]
  }

  n1 <- na1 + nb1
  n <- na + nb

  # TODO(Rosanne): Dit heb ik eerder gezien, misschien hier een functie voor schrijven?

  s_value <- sum(exp(
    #Bayes marginal alternative
    na1 * log(H1set[, 1]) +
      (na - na1) * log(1 - H1set[, 1]) +
      nb1 * log(H1set[, 2]) +
      (nb - nb1) * log(1 - H1set[, 2]) +
      log(w1) -
      #Bayes marginal null
      (n1 * log(point_h0) + (n - n1) * log(1 - point_h0))
  ))
  safe_p_value <- 1 / s_value

  colnames(d.table) <- c(0, 1)

  testResult <- list('data' = d.table,
                     'na' = na,
                     'nb' = nb,
                     's_value' = s_value,
                     'design' = testDesign)
  class(testResult) <- "safe2x2_test"

  return(testResult)
}


#' Alias for \code{\link{safeTwoProportionsTest}}
#'
#' @inheritParams safeTwoProportionsTest
#' @return safe2x2 test result object with the data used for the test, the S-value and the corresponding
#' safe p-value.
#' @export
#'
#' @examples
#' testDesign <- designPilotSafeTwoProportions(na = 10, nb = 10)
#' groups <- factor(c(rep("a", 10), rep("b", 10)))
#' observed_data <- factor(c(rep("a", 7), rep("b", 3),rep("a", 3), rep("b", 7)))
#' safe.proportion.test(x = groups, y = observed_data, testDesign = testDesign)
safe.proportion.test <- function(x, y = NULL, testDesign) {
  return(safeTwoProportionsTest(x, y, testDesign))
}
