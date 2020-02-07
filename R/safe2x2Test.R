create_empty_safe_2x2_design <- function() {
  new_safe_2x2_design <- list(
    'delta.star' = NA,
    'na' = NA,
    'nb' = NA,
    'point_h0' = 0.5,
    'H1set' = NA,
    'w1' = c(0.5, 0.5),
    'call' = NA
  )

  class(new_safe_2x2_design) <- "safe2x2_result"
  return(new_safe_2x2_design)
}

add_attributes <- function(object, attributes_defined) {
  object[names(attributes_defined)] <- attributes_defined
  return(object)
}

create_safe_2x2_design <- function(attributes_defined) {
  stopifnot(is.list(attributes_defined))
  safe_2x2_design <- create_empty_safe_2x2_design()
  safe_2x2_design <- add_attributes(safe_2x2_design, attributes_defined)
  return(safe_2x2_design)
}

#finds the greatest common divisor of two integers a and b
#returns this as 'n.iter'
#together with a/n.iter and b/n.iter
# findGreatestCommonDivisor <- function(a, b) {
#   imax <- 100
#
#   for (i in 1:imax) {
#     ratio <- i * a / b
#     if (abs(round(ratio) - ratio) < 1e-8) {
#       a.iter <- ratio
#       b.iter <- i
#       n.iter <- a / ratio
#       break
#     }
#   }
#
#   if (!exists("a.iter")) {
#     #no greatest common divisor has been found
#     a.iter <- a
#     b.iter <- b
#     n.iter <- 1
#   }
#
#   return(
#     list(
#       'a.iter' = a.iter,
#       'b.iter' = b.iter,
#       'n.iter' = n.iter
#     )
#   )
# }

#helper function for finding adjustment phi to simple S
#for case of unequal group sizes
#calculates the expected capital growth of a simple S value adjusted with
#value 'transl'
GetExpectedCapitalGrowthSimpleSWithAdjustment <- function(transl, H1set.neutral, n.grid, na,
                                                          nb, binom.coef, pbar0) {
  H1set.transl <-
    H1set.neutral + rbind(rep(transl, 2), rep(-transl, 2))

  # TODO(Rosanne): Zou je dit kunnen controleren en een betere naam kunnen geven?
  helpFunc <- function(h) {
    result <- exp((n.grid[, 1]) * log(h[1]) + (na - n.grid[, 1]) * log(1 - h[1]) +
                    (n.grid[, 2]) * log(h[2]) + (nb - n.grid[, 2]) * log(1 - h[2])
    )
    return(result)
  }

  pbar1 <- c(0.5, 0.5) %*% t(apply(H1set.transl, 1, helpFunc))
  return(sum(binom.coef * pbar1 * log(pbar1 / pbar0)))
}

#helper function to find derivative of the inverse of S: T -> 1/alpha
calculate_derivative_wrt_delta_for_S_inverse <- function(d, a, n) {
  # TODO(Rosanne): Denk je dat we dit kunnen simplificeren?
  (1 / (1 / ((1 - d ^ 2) ^ (0.5 * n) * a) + sqrt(1 / ((
    1 - d ^ 2
  ) ^ n * a ^ 2) - 1)) * n * d * (a ^ (-1)) * (1 / ((1 - d ^ 2) ^ (0.5 * n + 1)) + (1 /
                                                                                      a) * ((1 / a ^ 2) * (1 / ((1 - d ^ 2) ^ n
                                                                                      )) - 1) ^ (-0.5) *
                                                 (1 - d ^ 2) ^ (-n - 1)) * log((1 + d) / (1 - d)) -
    2 / (1 - d ^ 2) * log(1 / ((1 - d ^ 2) ^ (0.5 * n) * a) + sqrt(1 /
                                                                     ((
                                                                       1 - d ^ 2
                                                                     ) ^ n * a ^ 2) - 1))) / (log((1 + d) / (1 - d))) ^ 2
}

#helper function to find a delta ump in the case of na == nb
#throws error if uniroot fails
#returns numeric lenght 1 if success
get_delta_ump_analytically <- function(n, alpha, lower = 0.01, upper = 0.9) {
  # TODO(Rosanne): Denk je dat je tryOrFailWithNA() hiervoor helpt?
  delta.ump <-
    tryCatch(
      stats::uniroot(
        f = calculate_derivative_wrt_delta_for_S_inverse,
        interval = c(lower, upper),
        a = alpha,
        n = n
      )$root,
      error = function(e) {
        NA
      }
    )
  #if ( is.na( delta.ump)) { stop("Choose different range for delta.stars")}
  return(delta.ump)
}

calculate_derivative_wrt_theta_a_for_exp_log_S <- function(theta_a, na, nb, delta) {
  na * log((theta_a * (1 - theta_a - (nb / (na + nb)) * delta)) /
             ((theta_a +  (nb / (na + nb)) * delta) * (1 - theta_a))) +
    nb * log(((theta_a + delta) * (1 - theta_a - (nb / (na + nb)) * delta)) /
               ((theta_a +  (nb / (na + nb)) * delta) * (1 - theta_a  - delta)))
}

get_GROW_theta_tilde_one_sided <- function(delta, na, nb) {
  #theta_a can take on values between 0 and 1-delta in the less case
  #and between delta and 1 in the greater case

  # TODO(Rosanne): Denk je dat je tryOrFailWithNA() hiervoor helpt?
  derivative_root <-
    tryCatch(
      stats::uniroot(
        f = calculate_derivative_wrt_theta_a_for_exp_log_S,
        interval = c(max(-delta, 0), min(1 - delta, 1)),
        na = na,
        nb = nb,
        delta = delta
      )$root,
      error = function(e) {
        return(NA)
      }
    )

  if (is.na(derivative_root))
    return(NA)

  theta_tilde <- derivative_root + (nb / (na + nb)) * delta
  return(theta_tilde)
}

create_point_h1_and_h0_one_sided <- function(delta, na, nb) {
  if (na == nb) {
    theta_tilde <- 0.5
    theta_a <- 0.5 - 0.5 * delta
    theta_b <- 0.5 + 0.5 * delta
  } else {
    theta_tilde <-
      suppressWarnings(get_GROW_theta_tilde_one_sided(delta = delta, na = na, nb = nb))

    if (is.na(theta_tilde))
      return(NA)

    theta_a <- theta_tilde - delta / (1 + na / nb)
    theta_b <- theta_a + delta
  }

  if (any(c(theta_a, theta_b) < 0))
    return(NA)

  result <- list('H1set' = matrix(c(theta_a, theta_b), nrow = 1, byrow = TRUE),
                 'point_h0' = theta_tilde)

  return(result)
}


create_h1_set_equal_group_sizes <- function(delta.ump) {
  result <- matrix(1/2*c(1 - delta.ump, 1 + delta.ump,
                         1 + delta.ump, 1- delta.ump),
                   byrow = TRUE, nrow = 2)
  return(result)
}

create_h1_set_unequal_group_sizes <- function(delta_ump, na, nb) {
  pbar0 <- 0.5 ^ (na + nb)
  n.grid <- expand.grid(0:na, 0:nb)

  # TODO(Rosanne): Misschien is het robuster als je lchoose gebruikt en de twee termen optelt.
  binom.coef <- apply(n.grid, 1, function(nc) {
    return(choose(na, nc[1]) * choose(nb, nc[2]))
  })

  H1set.neutral <- create_h1_set_equal_group_sizes(delta_ump)

  # TODO(Rosanne): Heeft optimize ook een error? Dat kunnen we dan ook rapporteren als het nodig is
  transl <- stats::optimize(
    GetExpectedCapitalGrowthSimpleSWithAdjustment,
    lower = -H1set.neutral[1, 1],
    upper = H1set.neutral[1, 1],
    H1set.neutral = H1set.neutral,
    n.grid = n.grid,
    na = na,
    nb = nb,
    binom.coef = binom.coef,
    pbar0 = pbar0
  )$minimum
  H1set <- H1set.neutral + rbind(rep(transl, 2), rep(-transl, 2))
  return(H1set)
}

create_data_generating_distributions <- function(deltaDesign, alternative, length.out) {
  data_generating_distributions <- rbind(
    cbind(seq(deltaDesign, 1, length.out = length.out),
          seq(0, 1 - deltaDesign, length.out = length.out)),
    cbind(seq(0, 1 - deltaDesign, length.out = length.out),
          seq(deltaDesign, 1, length.out = length.out))
  )

  if (alternative == "greater") {
    data_generating_distributions <-
      data_generating_distributions[data_generating_distributions[, 1] > data_generating_distributions[, 2], ]
  }

  if (alternative == "less") {
    data_generating_distributions <-
      data_generating_distributions[data_generating_distributions[, 1] < data_generating_distributions[, 2], ]
  }

  return(data_generating_distributions)
}

perform_MC_simulations_for_each_data_generating_distribution <- function(H1.deltamin, na, nb, seed.time, M,
                                                                         H1set, w1, point_h0, alpha) {
    # reserve memory space for the current delta star results
    percent.reject.S <- numeric(nrow(H1.deltamin))
    set.seed(seed.time)

    n <- na + nb

    #start sampling
    for (i in 1:nrow(H1.deltamin)) {
      mu.a <- H1.deltamin[i, 1]
      mu.b <- H1.deltamin[i, 2]
      reject.S <- numeric(M)

      for (m in 1:M) {
        na1 <- sum(stats::rbinom(na, 1, mu.a))
        nb1 <- sum(stats::rbinom(nb, 1, mu.b))
        n1 <- na1 + nb1

        #perform the S-test
        s_value <- sum(exp(
          #Bayes marginal alternative
          na1 * log(H1set[, 1]) + (na - na1) * log(1 - H1set[, 1]) +
            nb1 * log(H1set[, 2]) + (nb - nb1) * log(1 - H1set[, 2]) + log(w1) -
            #Bayes marginal null
            (n1 * log(point_h0) + (n - n1) * log(1 - point_h0))
        ))

        reject.S[m] <- s_value >= (1 / alpha)
      }
      percent.reject.S[i] <- mean(reject.S)
    }
    return(percent.reject.S)
}

simulate_maximin_delta_and_power <- function(delta.min, na, nb, delta.stars, seed.time, M, alpha, alternative) {
  H1.deltamin <- create_data_generating_distributions(deltaDesign = delta.min, alternative = alternative,
                                                      length.out = 5)

  #reserve memory space for the current delta.min results
  all.powers <- matrix(nrow = length(delta.stars), ncol = nrow(H1.deltamin))

  #try out different delta star candidates for this delta min
  for (j in 1:length(delta.stars)) {
    delta <- delta.stars[j]

    #construct thetas in de alternative of the simple S-value
    if (alternative == "two.sided") {
      H1set <- create_h1_set_unequal_group_sizes(delta_ump = delta, na =  na, nb =  nb)
      w1 <- c(0.5, 0.5)
      point_h0 <- 0.5
    } else {
      #if alternative is greater, feed negative delta to create point h1 and h0
      delta_design <- ifelse(alternative == "greater", -delta, delta)
      point_H1_and_H0 <- create_point_h1_and_h0_one_sided(delta = delta_design, na = na, nb = nb)
      # TODO(Rosanne): Misschien is het beter om in functies [["H1set"]] te gebruiken
      H1set <- point_H1_and_H0[["H1set"]]
      w1 <- 1
      point_h0 <- point_H1_and_H0[["point_h0"]]
    }

    if (any(is.na(c(H1set, point_h0)))) {
      next
    }

    percent.reject.S <-
      perform_MC_simulations_for_each_data_generating_distribution(H1.deltamin = H1.deltamin, na = na, nb = nb,
                                                                     seed.time = seed.time, M = M, H1set = H1set,
                                                                     w1 = w1, point_h0 = point_h0, alpha = alpha)

    all.powers[j,] <- percent.reject.S
  }

  #which delta star provided the best power in the worst case for this delta min?
  best.delta <- delta.stars[which.max(apply(all.powers, 1, min))]
  supinf.power <- max(apply(all.powers, 1, min))

  return(list('bestDelta' = best.delta, 'supinfPower' = supinf.power))
}

#-----------------------------------------------------------------------------
#' Design a the best S-value for a test of two proportions in terms of power for
#' fixed group sizes and a given significance level.
#'
#' \code{designPilotSafeTwoProportions} finds for given sample size and significance level the
#' 'most powerful' simple S-value formula: the S-value formula that provides the best
#' power in the worst case for all distributions in the alternative hypothesis.
#' In case of equal group sizes and two-sided testing, the best S-value can be found
#' analytically (working paper by Grunwald, Ly and Turner).
#' Otherwise, power simulations are used to find the S-value formula with the
#' most power in the worst case for a set of data generating distributions in the
#' alternative hypothesis.
#'
#' @param na Group size a
#' @param nb Group size b
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less"
#' @param alpha significance level that will be used for testing.
#' Default \code{0.05}.
#' @param M number of iterations used in the power simulations.
#' Default \code{1000}.
#' @param deltaMins vector of effect size values between 0 and 1 to find the S-value formula with the best
#' worst-case power for. For each value of delta.min, the power of several S-values defined by
#' the test deltas is simulated for several data generating distributions where the group mean
#' in group a differs \code{delta.min} with group mean b. It is advised to choose a wide range.
#' Default \code{seq(0.1, 0.9, 0.1)}.
#' @param lowDelta numeric that defines the smallest delta of our search space for the test-defining deltaS
#' @param highDelta numeric that defines the largest delta of our search space for the test-defining deltaS
#' @param tol a number that defines the stepsizes between the lowDelta and highDelta
#' @param numberForSeed Optionally provide your own seed for the simulations.
#' Default \code{NA}.
#'
#' @return Safe design object that can be used within \code{\link{safeTwoProportionsTest}}
#' @export
#'
#' @examples
#' safeDesign <- designPilotSafeTwoProportions(na = 10, nb = 10, alpha = 0.05)
#' \dontrun{
#' safeDesign <- designPilotSafeTwoProportions(na = 15, nb = 10, alpha = 0.05)
#' }
#'
designPilotSafeTwoProportions <- function(na, nb, alternative = c("two.sided", "greater", "less"),
                                          alpha = 0.05, M = 100, deltaMins = seq(0.1, 0.9, 0.1),
                                          lowDelta = 0.1, highDelta = 0.8, tol = 0.1,
                                          numberForSeed = NA) {
  alternative <- match.arg(alternative)
  safe_2x2_design <- create_safe_2x2_design(list(na = na, nb = nb))

  delta.stars <- seq(lowDelta, highDelta, by = tol)

  #check if delta ump can be found analytically
  if (na == nb & alternative == "two.sided") {
    delta.ump <- get_delta_ump_analytically(n = na + nb, alpha = alpha, lower = min(delta.stars),
                                            upper = max(delta.stars))

    safe_2x2_design <- add_attributes(safe_2x2_design, list('delta.star' = delta.ump,
                                                            'H1set' = create_h1_set_equal_group_sizes(delta.ump),
                                                            'call' = sys.call())
    )
    return(safe_2x2_design)
  }

  #--------------------------------------
  #na not nb case or one-sided: MC simulations

  #determine the way the seed is set for the MC simulations
  seed.time <-
    ifelse(is.na(numberForSeed), Sys.time(), numberForSeed)

  #reserve memory space for results
  best.delta <- supinf.power <- numeric(length(deltaMins))

  #start simulations for determining power
  for (d in 1:length(deltaMins)) {
    #define the grid of data-generating parameters to investigate power for
    delta.min <- deltaMins[d]
    simulation_result <- simulate_maximin_delta_and_power(delta.min = delta.min, na = na, nb = nb,
                                                          delta.stars = delta.stars, seed.time = seed.time,
                                                          M = M, alpha = alpha, alternative = alternative
    )
    best.delta[d] <- simulation_result[['bestDelta']]
    supinf.power[d] <- simulation_result[['supinfPower']]
  }

  # retrieve the mode in the table
  # as there is a 'most powerful region', due to simulation error,
  # different values might be found
  # we take the most prevalent value as a solution
  delta.ump <- as.numeric(names(sort(table(best.delta), decreasing = TRUE)[1]))

  if (alternative == "two.sided") {
    if (na == nb) {
      H1set <- create_h1_set_equal_group_sizes(delta.ump)
      w1 <- c(0.5, 0.5)
      point_h0 <- 0.5
    } else {
      H1set <- create_h1_set_unequal_group_sizes(delta_ump = delta.ump, na =  na, nb =  nb)
      w1 <- c(0.5, 0.5)
      point_h0 <- 0.5
    }
  } else {
    #if alternative is greater, feed negative delta to create point h1 and h0
    delta_design <- ifelse(alternative == "greater", -delta.ump, delta.ump)
    point_H1_and_H0 <- create_point_h1_and_h0_one_sided(delta = delta_design, na = na, nb = nb)
    H1set <- point_H1_and_H0[["H1set"]]
    w1 <- 1
    point_h0 <- point_H1_and_H0[["point_h0"]]
  }

  safe_2x2_design <- add_attributes(safe_2x2_design, list('delta.star' = delta.ump, 'H1set' = H1set,
                                                          'call' = sys.call(), 'w1' = w1,'point_h0' = point_h0))
  return(safe_2x2_design)
}


#-----------------------------------------------------------------------------

#' Design an S-value for a test of two proportions with a certain minimal power
#' for detecting at least a certain difference delta.min in proportionts between
#' group a and b.
#'
#' \code{designSafeTwoProportions} finds for a given delta.min indicating the minimal
#' mean difference one wants to detect between group a and b, and a given significance level alpha
#' and desired type II error guarantee beta the smallest sample size that is needed to detect
#' that difference. For this sample size, the S-value formula achieving the desired power is
#' found as well and can be used for performing the a safe test.
#'
#' As a side effect, prints each sample size being tried.
#'
#'
#' @param deltaMin numeric that defines the minimal absolute difference between the group means the
#' safe test should be designed for, i.e. the minimal relevant difference we want to be able to detect.
#' @param alpha numeric in (0, 1) that specifies the target type 1 error control --independent on n-- that the
#' designed test has to adhere to. Note that it also defines the rejection rule S10 > 1/alpha. Default
#' \code{0.05}.
#' @param beta numeric in (0, 1) that speficies the target type II error control necessary to calculate both "n"
#' and "deltaS". Note that 1-beta defines the power. Default \code{0.20}.
#' @param alternative a character string specifying the alternative hypothesis must be one of "two.sided" (default),
#' "greater" or "less"
#' @param M M Number of simulations used per sample size, per possible S-value, per
#' data-generating distribution. Default \code{100}.
#' @param lowN integer that defines the smallest n of our search space for n. Default 20.
#' @param highN integer that defines the largest n of our search space for n. This might be the largest n that we
#' are able to fund. Default 200.
#' @param sampleSizeRatio numeric representing nb/na. Default 1.
#' @param lowDelta numeric that defines the smallest delta of our search space for the test-defining deltaS
#' @param highDelta numeric that defines the largest delta of our search space for the test-defining deltaS
#' @param tol a number that defines the stepsizes between the lowDelta and highDelta
#' @param numberForSeed Optionally provide your own seed for the simulations.
#' Default \code{NA}.
#'
#' @return Result object that can be used within \code{\link{safeTwoProportionsTest}}
#' @export
#'
#' @examples
#' designSafeTwoProportions(deltaMin = 0.7)
#' designSafeTwoProportions(deltaMin = 0.7, sampleSizeRatio = 2, lowN = 30)
designSafeTwoProportions <- function(deltaMin, alpha = 0.05, beta = 0.20,
                                     alternative = c("two.sided", "greater", "less"), M = 100,
                                     lowN = 20, highN = 200, sampleSizeRatio = 1,
                                     lowDelta = 0.1, highDelta = 0.8, tol = 0.1, numberForSeed = NA) {

  alternative <- match.arg(alternative)

  #initialize starting parameters
  deltas <- seq(lowDelta, highDelta, by = tol)

  #calculate the number of complete groups a and b collected at low N
  #for example lowN = 20, sampleSizeRatio = 1, then we start with 10 groups
  #that have been collected
  nGroups <- lowN/(sampleSizeRatio+1)
  if (abs(round(nGroups) - nGroups) > 1e-8)
    stop("lowN should be multiple of a.iter + b.iter")

  #in the iteration we start directly with incrementing nGroups,
  #so now subtract 1
  nGroups <- nGroups - 1
  nb <- sampleSizeRatio*(nGroups)
  na <- nGroups
  n <- na + nb

  #for H1_delta.min, we want to test the power:
  H1DeltaMin <- create_data_generating_distributions(deltaDesign = deltaMin, alternative = alternative,
                                                      length.out = 5)

  notreached <- TRUE

  #save the start time to use as a seed inside the sampling process
  #so the same set of samples from H1 delta.circ will be presented
  #to each S-test constructed with delta.star

  seedTime <- ifelse(is.na(numberForSeed), Sys.time(), numberForSeed)

  cat("Trying n = ")

  while (notreached & (n < highN)) {
    #increment the group sizes
    nGroups <- nGroups + 1
    na <- nGroups
    nb <- sampleSizeRatio*nGroups
    if(nb %%1 != 0){
      #no whole group
      next
    }

    n <- na + nb
    cat(paste0(n, " "))

    #if na = nb and we are testing twosided we only need to test the delta ump,
    #since this is then known analytically

    if (na == nb & alternative == "two.sided") {
      delta.ump <- get_delta_ump_analytically(n = na + nb, alpha = alpha,
                                              lower = 0.1, upper = 0.9)

      if (!is.na(delta.ump)) {
        deltas <- delta.ump
      }
    }

    #start iterating over possible delta stars at this n
    for (j in 1:length(deltas)) {
      delta <- deltas[j]

      #construct parameters to use in the S-test
      if (alternative == "two.sided") {
        if (na == nb) {
          H1set <- create_h1_set_equal_group_sizes(delta)
          w1 <- c(0.5, 0.5)
          point_h0 <- 0.5
        } else {
          H1set <- create_h1_set_unequal_group_sizes(delta, na = na, nb = nb)
          w1 <- c(0.5, 0.5)
          point_h0 <- 0.5
        }
      } else {
        #if alternative is greater, feed negative delta to create point h1 and h0
        delta_design <- ifelse(alternative == "greater", -delta, delta)
        point_H1_and_H0 <- create_point_h1_and_h0_one_sided(delta = delta_design, na = na, nb = nb)

        if (any(is.na(point_H1_and_H0))) {
          next
        }

        H1set <- point_H1_and_H0$H1set
        w1 <- 1
        point_h0 <- point_H1_and_H0$point_h0
      }

      percent.reject.S <-
        perform_MC_simulations_for_each_data_generating_distribution(H1.deltamin = H1DeltaMin, na = na, nb = nb,
                                                                     seed.time = seedTime, M = M, H1set = H1set,
                                                                     w1 = w1, point_h0 = point_h0, alpha = alpha)

      #check if desired power has been reached at desired precision
      if (min(percent.reject.S) >= (1 - beta)) {
        notreached <- FALSE

        cat("\nFor all p1, power above desired level. Worst case: ")
        cat(min(percent.reject.S))
        cat(" with data generated from: ")
        cat(H1DeltaMin[which.min(percent.reject.S),]); cat(".\n")
        break
      }
    }
  }

  if (notreached)
    stop("Desired power for this alternative hypothesis not reached below maximal n value")

  safe_design <- create_safe_2x2_design(list('n.star' = n, 'delta.star' = delta, 'na' = na,
                                             'nb' = nb, 'w1' = w1, 'point_h0' = point_h0,
                                             'H1set' = H1set, 'call' = sys.call(), 'alpha' = alpha,
                                             'beta' = beta))

  return(safe_design)
}

#------------------------------------------------------------------------------
# frequentist experiment: when will we stop?

#' Simulate the spread of realized sample sizes during an optional stopping design setup.
#'
#' In simulations, data are collected in groups of 2 (1 from group a, and 1 from group b)
#' and s-values are evaluated after each group of 2 is collected. When an S-value
#' exceeds the critical value as defined by the alpha specified in the design,
#' that simulation is stopped, the null hypothesis is recorded to be rejected,
#' and the realized sample size is saved. One can then calculate the expected
#' sample size collected in the optional stopping setting, and the power of the S-value
#' in the optional stopping setting (see examples).
#'
#' Side effects: plots a histogram of realized sample sizes in the optional
#' stopping setting.
#'
#' @param safeDesign a safe design for a test of two proportions retrieved
#' through \code{\link{designSafeTwoProportions}}
#' @param M number of simulations to carry out
#' @param parametersDataGeneratingDistribution group means a and b in the data
#' generating distribution to simulate from.
#'
#' @return \code{n.final}, vector of realized sample sizes, and \code{rejected},
#' logical vector indicating whether the null hypothesis was rejected in that simulation.
#' @export
#'
#' @examples
#' safeDesignProportionTest <- designSafeTwoProportions(deltaMin = 0.4, lowN = 80)
#' #sample size planned, according to design
#' safeDesignProportionTest$n.star
#'
#' simulationResult <- simulateSpreadSampleSizeTwoProportions(safeDesignProportionTest,
#' M = 1000, parametersDataGeneratingDistribution = c(0.2, 0.6))
#' #what was the power?
#' mean(simulationResult$rejected)
#' #what sample size can be expected to be collected with optional stopping?
#' mean(simulationResult$actually_collected)
#'
simulateSpreadSampleSizeTwoProportions <- function(safeDesign, M, parametersDataGeneratingDistribution) {

  H1set <- safeDesign$H1set
  na.max <- safeDesign$na
  nb.max <- safeDesign$nb
  critical_value_s <- 1 / safeDesign$alpha

  if (na.max != nb.max)
    stop("optional stopping simulation only available for equal group sizes")

  if (nrow(H1set) != 2)
    stop("optional stopping simulation only available for two-sided test")

  #for H1_delta.min, we want to retrieve E(stopping time) IN THE WORST CASE
  H1.deltamin <- parametersDataGeneratingDistribution

  theta.a <- H1.deltamin[1]
  theta.b <- H1.deltamin[2]
  rejected <- logical(M)
  n.final <- numeric(M)

  for (i in 1:M) {
    #perform the "optional stopping" experiment: collect data sequentially
    a.sample <- stats::rbinom(n = na.max, size = 1, prob = theta.a)
    b.sample <- stats::rbinom(n = nb.max, size = 1, prob = theta.b)

    #retrieve na1 and nb1 at n=2, n=4, n=6, etc
    na <- 1:na.max
    nb <- 1:nb.max

    na1 <- sapply(1:na.max, function(n.cur) {sum(a.sample[1:n.cur])})

    nb1 <- sapply(1:nb.max, function(n.cur) {sum(b.sample[1:n.cur])})

    calculateLikelihoodUnderH <- function(h) {
      exp(na1 * log(h[1]) + (na - na1) * log(1 - h[1]) + nb1 * log(h[2]) +
            (nb - nb1) * log(1 - h[2]))
    }
    #retrieve pbar1 at each time point (is a vector now!)
    pbar1 <- rowSums(0.5 * apply(H1set, 1, calculateLikelihoodUnderH))

    #n is a vector now! And so is pbar0!
    n <- na + nb
    pbar0 <- 0.5 ^ n

    if (sum(pbar1 / pbar0 >= critical_value_s) > 0) {
      #at some time point, we rejected: record the first time point
      rejected[i] <- TRUE
      n.final[i] <- n[min(which(pbar1 / pbar0 >= critical_value_s))]
    } else {
      #we have not rejected and stopped collecting according to our stopping rule
      rejected[i] <- FALSE
      n.final[i] <- na.max + nb.max
    }
  }

  graphics::hist(n.final, breaks = na.max+nb.max, xlim = c(0, max(n.final)), xlab = "n collected",
                 main = bquote(~"Mean 1" == .(parametersDataGeneratingDistribution[1]) ~"," ~"Mean 2" == .(parametersDataGeneratingDistribution[2])  ~"," ~"alpha" == ~.(safeDesign$alpha)))
  return(list("actually_collected" = n.final, "rejected" = rejected))
}

#' Function that can be used to illustrate that H0 is falsely rejected too often when using
#' Fisher's exact test in the optional stopping setting.
#'
#' Simulates what would happen if we would use Fisher's exact test in the optional
#' stopping setting. First, determines the sample size we should use according to design with
#' simulations. Can be omitted if planned sample size is already known and passed through
#' Then, starts sampling from the passed distribution up to this sample size in each simulation, and stops
#' if a significant Fisher's exact test is recorded during sampling.
#'
#' @param deltaDesign minimal difference in mean proportions between group a and b
#' we want to be able to detect with a certain power.
#' @param alpha significance level used for testing.
#' @param power minimal power we want to achieve
#' @param M number of simulations to carry out.
#' @param parametersDataGeneratingDistribution Vector with actual values of the means in groups a and
#' b in the data generating distribution to simulate for. Values between 0 and 1.
#' @param nDesign if n has been designed/ determined previously, it can be passed
#' here and simulations to find the design sample size are omitted. This saves
#' time. Default \code{NA}.
#' @param highN maximal nDesign: if no nDesign is found to achieve the desired power
#' for the \code{deltaDesign} below \code{highN}, simulations are stopped.
#'
#' @return see \code{\link{simulateSpreadSampleSizeTwoProportions}}
#' @export
#'
#' @examples
#' \dontrun{
#' simulationResult <- simulateFisherSpreadSampleSizeOptionalStopping(deltaDesign = 0.5,
#' alpha = 0.05, power = 0.8, M = 100,
#' parametersDataGeneratingDistribution = c(0.4, 0.4))
#' #what was the power?
#' mean(simulationResult$rejected)
#' #what sample size can be expected to be collected with optional stopping?
#' mean(simulationResult$actually_collected)
#'
#' #nDesign already known
#' simulationResult <- simulateFisherSpreadSampleSizeOptionalStopping(deltaDesign = 0.5,
#' alpha = 0.05, power = 0.8, M = 1000, nDesign = 40,
#' parametersDataGeneratingDistribution = c(0.4, 0.4))
#' }
#'
simulateFisherSpreadSampleSizeOptionalStopping <- function(deltaDesign, alpha, power,
                                                           M, parametersDataGeneratingDistribution,
                                                           nDesign = NA, highN = 200) {

  if (is.na(nDesign)) {
    cat("nDesign not known, determining sample size through power simulations.\nTrying n: ")
    #determine design na and nb for fisher
    #for H1_delta.min, we want to test the power:
    H1.deltamin <- create_data_generating_distributions(deltaDesign, alternative = "two.sided", length.out = 5)

    #condition for looping
    # TODO(Rosanne): Hier heb TRUE van gemaakt. Klopt dat?

    not.satisfied <- TRUE
    n <- 20
    seed.time <- Sys.time()

    #case na=nb now
    while (not.satisfied) {
      n <- n + 2
      na <- nb <- n / 2

      cat(n); cat(" ")

      #reserve memory space for power estimate for each distribution in
      #the alternative of interest
      percent.reject.F <- numeric(nrow(H1.deltamin))

      #start MC simulation
      set.seed(seed.time)
      for (i in 1:nrow(H1.deltamin)) {
        #data are generated by a distribution in H1.deltamin
        mu.a <- H1.deltamin[i, 1]
        mu.b <- H1.deltamin[i, 2]
        reject.F <- numeric(M)

        for (m in 1:M) {
          na1 <- sum(stats::rbinom(na, 1, mu.a))
          nb1 <- sum(stats::rbinom(nb, 1, mu.b))

          # TODO(Rosanne): Dit heb ik uit elkaar getrokken. Klopt dit nog?
          somePValue <- stats::fisher.test(x = matrix(c(na - na1, na1, nb - nb1, nb1),
                                               byrow = TRUE, nrow = 2))$p.value

          reject.F[m] <- somePValue <= alpha
        }

        percent.reject.F[i] <- mean(reject.F)
      }

      if (min(percent.reject.F) >= power) {
        # TODO(Rosanne): Hier heb ik FALSE van gemaakt klopt dat?
        not.satisfied <- F
      }

      if (n > highN) {
        stop("highN too low or deltaDesign too low; desired power not achieved with these settings.")
      }

    }

    na.max <- nb.max <- na
    nDesign <- n
    cat("\nFound Fisher sample design size: ")
    cat(nDesign); cat(".\n")
  } else {
    na.max <- nb.max <- nDesign / 2
  }

  #for H1_delta.min, we want to retrieve E(stopping time) IN THE WORST CASE
  H1.deltamin <- parametersDataGeneratingDistribution

  theta.a <- H1.deltamin[1]
  theta.b <- H1.deltamin[2]

  rejected <- logical(M)
  n.final <- rep(nDesign, times = M)

  cat("Starting optional stopping simulations:\n")
  pbar <- utils::txtProgressBar(min = 1, max = M)


  for (i in 1:M) {
    utils::setTxtProgressBar(pbar, i)

    #perform the "optional stopping" experiment: collect data sequentially
    a.sample <- stats::rbinom(n = na.max, size = 1, prob = theta.a)
    b.sample <- stats::rbinom(n = nb.max, size = 1, prob = theta.b)

    #retrieve na1 and nb1 at n=2, n=4, n=6, etc
    na <- 1:na.max

    # TODO(Rosanne): Ik dacht er aan om hier een generieke functie te maken die een functie uit spuugt
    # Omdat je dit ook eerder al doet. Zie "function factory"

    na1 <- sapply(1:na.max, function(n.cur) {sum(a.sample[1:n.cur])})
    nb1 <- sapply(1:nb.max, function(n.cur) {sum(b.sample[1:n.cur])})

    for (j in 1:length(na)) {
      na.cur <- nb.cur <- na[j]
      na1.cur <- na1[j]
      nb1.cur <- nb1[j]

      somePValue <- stats::fisher.test(x = matrix(c(na.cur - na1.cur, na1.cur, nb.cur - nb1.cur, nb1.cur),
                                           byrow = TRUE,nrow = 2))$p.value

      if (somePValue <= alpha) {
        n.final[i] <- na.cur + nb.cur
        rejected[i] <- TRUE
        break
      }
    } # TODO(Rosanne): Wat denk je van de volgende comment: End for loop over the first sample
  } # TODO(Rosanne): Wat denk je van de volgende comment: End the number of iterations
  close(pbar)

  graphics::hist(n.final, breaks = nDesign, xlim = c(0, max(n.final)), xlab = "n collected",
                 main = bquote(~"Mean 1" == .(parametersDataGeneratingDistribution[1]) ~","
                               ~"Mean 2" == .(parametersDataGeneratingDistribution[2])  ~"," ~"alpha" == ~.(alpha)))
  return(list("actually_collected" = n.final, "rejected" = rejected))
}



#' Function that illustrates that experiments can, on average, be stopped sooner
#' than planned in the optional stopping setting with tests of two proportions.
#'
#' For each effect size (absolute difference between mean proportions in group a
#' and b), the 'design' sample size and S-value to reach at least a power are found
#' through simulations. Then, for each found combination of delta min, n design and the S-value,
#' the optional stopping setting is simulated \code{M} times and the average number
#' of samples collected during optional stopping is recorded.
#'
#' A plot is then generated, illustrating the difference between the designed sample
#' size and the sample size collected during optional stopping.
#'
#' Note that this function only demonstrates the two-sided setting, with equal
#' group sizes.
#'
#' @param alpha Significance level used for testing.
#' @param beta Maximally acceptable type-II error.
#' @param maxN numeric, the maximum number of samples one has budget for to collect data.
#' Default 100.
#' @param deltaMinVec Vector with effect sizes one would want to detect with an
#' experiment.
#' @param numberForSeed Optionally, a seed can be set to make results replicable.
#' @param M Number of simulations to be carried out per data generating distribution.
#' Default 200.
#' @param highN Maximal n to sample while determining the design n for each. Default 200.
#' Increasing this will increase runtime of function significantly.
#' \code{delta.min}.
#'
#' @return The data used to generate the plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plotSafeTwoProportionsSampleSizeProfile(alpha = 0.05, beta = 0.20, highN = 100)
#' }
plotSafeTwoProportionsSampleSizeProfile <- function(alpha, beta, maxN = 100, deltaMinVec = seq(0.9, 0.2, -0.1),
                                                    numberForSeed = NA, M = 200, highN = 200) {
  #determine seed
  seed.time <- ifelse(is.na(numberForSeed), Sys.time(), numberForSeed)

  #initialize parameters
  deltaMinVec <- sort(deltaMinVec, decreasing = TRUE)
  n.min.vec <- numeric(length(deltaMinVec))
  n <- 10
  power <- 1 - beta

  #progress bar
  cat("Simulation progress part 1, determining n design:\n")
  pbar <- utils::txtProgressBar(min = 1, max = length(deltaMinVec))

  #---------------------------------------------------
  #start filling the plot data for each delta min
  #to retrieve n 'standard/ frequentist':
  for (d in 1:length(deltaMinVec)) {
    #continue with the previous n: lowN monotonically decreasing in deltamin
    utils::setTxtProgressBar(pbar, d)

    delta.min <- deltaMinVec[d]
    #for H1_delta.min, we want to test the power:
    H1.deltamin <- rbind(
      cbind(seq(delta.min, 1, length.out = 5),
            seq(0, 1 - delta.min, length.out = 5)),
      cbind(seq(0, 1 - delta.min, length.out = 5),
            seq(delta.min, 1, length.out = 5))
    )

    #condition for looping
    not.satisfied <- TRUE

    #case na=nb now
    while (not.satisfied) {
      n <- n + 2
      na <- nb <- n / 2

      #design the simple S for this n
      pbar0 <- 0.5 ^ n
      delta <- get_delta_ump_analytically(n, alpha)

      if (is.na(delta))
        next

      H1set <- matrix(1/2*c(1- delta, 1 + delta,
                            1 + delta, 1 - delta),
                      byrow = TRUE, nrow = 2)

      #reserve memory space for power estimate for each distribution in
      #the alternative of interest
      percent.reject.S <- numeric(nrow(H1.deltamin))

      #start MC simulation
      set.seed(seed.time)

      for (i in 1:nrow(H1.deltamin)) {
        #data are generated by a distribution in H1.deltamin
        mu.a <- H1.deltamin[i, 1]
        mu.b <- H1.deltamin[i, 2]
        reject.S <- numeric(M)

        for (m in 1:M) {
          na1 <- sum(stats::rbinom(na, 1, mu.a))
          nb1 <- sum(stats::rbinom(nb, 1, mu.b))

          #perform the safe test designed with delta star on the sample

          helpFunc <- function(h) {
            result <- exp(na1 * log(h[1]) + (na - na1) * log(1 - h[1]) +
                            nb1 * log(h[2]) + (nb - nb1) * log(1 - h[2]))
            return(result)
          }

          pbar1 <- sum(0.5 * apply(H1set, 1, helpFunc))

          reject.S[m] <- (pbar1 / pbar0) >= (1 / alpha)
        } # End iterations

        percent.reject.S[i] <- mean(reject.S)
      }

      if (min(percent.reject.S) >= power)
        not.satisfied <- FALSE

      if (n > highN)
        break

    } # TODO(Rosanne): End while satisfied loop for fixed deltamin

    if (n > highN)
      break

    n.min.vec[d] <- n
  } # TODO(Rosanne): End looping over all deltaMin
  close(pbar)

  #---------------------------------------------------
  #to retrieve n optional stop
  n.s.opt <- numeric(length(deltaMinVec))
  cat("Simulation progress part 2, determining n optional stop:\n")
  pbar <- utils::txtProgressBar(min = 1, max = length(deltaMinVec))

  for (d in 1:length(deltaMinVec)) {
    utils::setTxtProgressBar(pbar, d)
    #maximal n: the 'design' n retrieved in the previous simulations
    n.min <- n.min.vec[d]
    na.max <- nb.max <- n.min / 2

    delta <- get_delta_ump_analytically(n.min, alpha)

    if (is.na(delta))
      next

    H1 <- matrix(1/2*c(1 - delta, 1 + delta,
                       1 + delta, 1 - delta),
                 byrow = TRUE, nrow = 2)

    delta.min <- deltaMinVec[d]

    # TODO(Rosanne): Dit zie ik nu een paar keer, maar ik begrijp het niet helemaal
    #
    #for H1_delta.min, we want to retrieve E(stopping time) IN THE WORST CASE
    H1.deltamin <- rbind(
      cbind(seq(delta.min, 1, length.out = 8),
            seq(0, 1 - delta.min, length.out = 8)),
      cbind(seq(0, 1 - delta.min, length.out = 8),
            seq(delta.min, 1, length.out = 8))
    )

    mean.s.os <- numeric(nrow(H1.deltamin))

    set.seed(seed.time)

    for (h in 1:nrow(H1.deltamin)) {
      theta.a <- H1.deltamin[h, 1]
      theta.b <- H1.deltamin[h, 2]

      rejected <- logical(M)
      n.final <- numeric(M)


      for (i in 1:M) {
        #perform the "optional stopping" experiment: collect data sequentially

        a.sample <- stats::rbinom(n = na.max, size = 1, prob = theta.a)
        b.sample <- stats::rbinom(n = nb.max, size = 1, prob = theta.b)

        #retrieve na1 and nb1 at n=2, n=4, n=6, etc
        na <- 1:na.max
        nb <- 1:nb.max

        # TODO(Rosanne): Function factory zoals boven
        na1 <- sapply(1:na.max, function(n.cur) {sum(a.sample[1:n.cur])})
        nb1 <- sapply(1:nb.max, function(n.cur) {sum(b.sample[1:n.cur])})

        # TODO(Rosanne): Is dit nog okay? Ken je een betere naam hiervoor?
        helpFunc <-  function(h) {
          exp(na1 * log(h[1]) + (na - na1) * log(1 - h[1]) + nb1 * log(h[2]) +
                (nb -nb1) * log(1 - h[2]))
        }

        #retrieve pbar1 at each time point (is a vector now!)
        pbar1 <- rowSums(0.5 * apply(H1, 1, helpFunc))

        #n is a vector now! And so is pbar0!
        n <- na + nb
        pbar0 <- 0.5 ^ n

        if (sum(pbar1 / pbar0 >= 20) > 0) {
          #at some time point, we rejected: record the first time point
          rejected[i] <- TRUE
          n.final[i] <- n[min(which(pbar1 / pbar0 >= 20))]
        } else {
          #we have not rejected and stopped collecting according to our stopping rule
          rejected[i] <- FALSE
          n.final[i] <- na.max + nb.max
        }
      }

      mean.s.os[h] <- mean(n.final)
    }
    #what was the worst case?
    n.s.opt[d] <- max(mean.s.os)
  }
  close(pbar)

  #--------------------------------
  #save the results and plot
  resultForPlot <- list()
  resultForPlot[["deltaMinVec"]] <- deltaMinVec[n.min.vec > 0]
  resultForPlot[["n.min.vec"]] <- n.min.vec[n.min.vec > 0]
  resultForPlot[["n.os.vec"]] <- n.s.opt[n.min.vec > 0]

  graphics::par(cex.main=1.5, mar=c(5, 6, 4, 7)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
                font.lab=2, cex.axis=1.3, bty="n", las=1)

  graphics::plot(resultForPlot[["deltaMinVec"]], resultForPlot[["n.min.vec"]],
                 ylim = c(0, highN), type="l", col="blue", lty=2, lwd=2,
                 ylab="n1", xlab=expression(delta["min"]),
                 main=bquote(~alpha == ~.(alpha) ~ "and" ~beta== ~.(beta)))

  graphics::lines(resultForPlot[["deltaMinVec"]], resultForPlot[["n.os.vec"]],  col="black", lwd=2, lty=1)

  graphics::abline(h=maxN, col="red", lty=2)

  legendName <- c("Average n", "Safe design", "max n")
  legendCol <- c("black", "blue", "red")
  legendLty <- c(1, 2, 2)

  graphics::legend("topright", legend = legendName, col = legendCol, lty=legendLty, bty="n")

  return(resultForPlot)
}
