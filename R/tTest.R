# Testing fnts -------

#' Computes E-Values Based on the T-Statistic
#'
#' A summary stats version of \code{\link{saviTTest}()} with
#' the data replaced by t, n1 and n2, and the design object by parameter
#'
#' @inheritParams designSaviT
#' @inheritParams saviTTest
#' @param t numeric that represents the observed t-statistic.
#' @param parameter numeric > 0, the savi test defining parameter,
#' see \code{\link{matchEParameterWith}} for details.
#' @param n1 integer that represents the size of the (first) sample.
#' Default n2=\code{NULL} implies a one-sample T-test.
#' @param n2 an optional integer that represents the size of the second sample,
#' which implies a two-sample t-statistic
#'
#' @return Returns a numeric that represent the e10, that is, the e-value in favour of the alternative over the null
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'   `r addCite(perez2024estatistics)`
#'   `r addCite(wang2025anytime)`
#'
#' @export
#'
#' @examples
#' saviTTestStat(t=1, n1=100, parameter=0.4)
#' saviTTestStat(t=3, n1=100, parameter=0.3)
saviTTestStat <- function(
    t, n1, n2=NULL, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE, nuMin=2,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    ...) {

  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  nEff <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1 else (1/n1+1/n2)^(-1)
  nu <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1-1 else n1+n2-2

  result <- suppressWarnings(
    saviTTestStatNEffNu("t"=t, "nEff"=nEff, "nu"=nu,
                        "parameter"=parameter, "alternative"=alternative,
                        "tDensity"=tDensity, "paired"=paired, "eType"=eType,
                        "nuMin"=nuMin, ...)
  )

  return(result)
}


#' SaviTTestStat based on the t-statistic, nEff and nu
#'
#' @rdname saviTTestStat
#' @inheritParams computeConfidenceIntervalT
#'
#' @param nEff numeric > 0, the effective sample size. For one sample tests,
#' this is just n.
#' @param nu numeric > 0, the degrees of freedom.
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'   `r addCite(perez2024estatistics)`
#'   `r addCite(wang2025anytime)`
#'
#' @export
saviTTestStatNEffNu <- function(
    t, nEff, nu, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE, nuMin=2,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    ...) {

  if (nu < 0 || nEff < 0)
    return(list(eValue=1))

  # stopifnot(nEff >= 0, nu >= 0)

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  # Note(Alexander): For "greater" and "less" this might not be ideal
  if (nu <= 1)
    return(list("eValue"=1))

  if (is.infinite(t) || is.na(t) && nu <= nuMin)
    return(list("eValue"=1))

  if (eType=="grow") {

    tempResult <- saviTTestStatNEffNuGrow(
      t=t, nEff=nEff, nu=nu, parameter=parameter,
      alternative=alternative, tDensity=tDensity,
      paired=paired, ...)

    return(tempResult)
  } else if (eType=="eCauchy") {
    kappaG <- parameter

    if (alternative=="twoSided") {
      twoSidedCauchyIntegrand <- function(g) {
        exp(-1/2*log(1+nEff*g)+(nu+1)/2*(log(1+t^2/nu)-log(1+t^2/(nu*(1+nEff*g))))
            - 2*log(g) + stats::dgamma(x=1/g, shape=1/2, rate=kappaG^2/2, log=TRUE))
      }

      tempResult <- stats::integrate(twoSidedCauchyIntegrand, 0, Inf)
    } else {
      oneSidedCauchyIntegrand <- function(delta) {
        2/kappaG*exp(
          stats::dt(t, df=nu, ncp=sqrt(nEff)*delta, log=TRUE)
          -stats::dt(t, df=nu, ncp=0, log=TRUE)
          +stats::dt(delta/kappaG, df=1, ncp=0, log=TRUE)
        )
      }

      if (alternative=="greater") {
        tempResult <- stats::integrate(oneSidedCauchyIntegrand, 0, Inf)
      } else if (alternative=="less") {
        tempResult <- stats::integrate(oneSidedCauchyIntegrand, -Inf, 0)
      }
    }

    result <- list("eValue"=tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
    return(result)
  } else if (eType=="eGauss") {
    g <- parameter

    if (alternative=="twoSided") {
      logResult <- -1/2*log(1+nEff*g)+((nu+1)/2)*(log((1+t^2/nu))-log(1+t^2/(nu*(1+nEff*g))))

      return(list("eValue"=exp(logResult)))
    } else {
      oneSidedGaussIntegrand <- function(delta) {
        2*exp(
          stats::dt(t, df=nu, ncp=sqrt(nEff)*delta, log=TRUE)
          -stats::dt(t, df=nu, ncp=0, log=TRUE)
          +stats::dnorm(delta, mean=0, sd=sqrt(g), log=TRUE)
        )
      }

      if (alternative=="greater") {
        tempResult <- stats::integrate(oneSidedGaussIntegrand, 0, Inf)
      } else if (alternative=="less") {
        tempResult <- stats::integrate(oneSidedGaussIntegrand, -Inf, 0)
      }
    }

    result <- list("eValue"=tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
    return(result)
  } else if (eType=="lai") {

    if (alternative!="twoSided")
      stop("Lai's e-variable currently only for twoSided tests")

    eValue <- (nu+1)/2*log(1+t^2/nu)+1/2*log(2*pi/nEff)
    return(list(eValue=eValue))
  } else if (eType=="mom") {
    tempResult <- saviTTestStatNEffNuMom(
      t=t, nEff=nEff, nu=nu, parameter=parameter,
      alternative=alternative, tDensity=tDensity,
      paired=paired)

    return(tempResult)
  }
}



#' SaviTTestStat based on the t-statistic, nEff and nu and the non-local moment prior
#'
#' @rdname saviTTestStat
#' @inheritParams saviTTestStatNEffNu
#'
#' @param k the moment used for the non-local moment prior. Default 1
#'
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'   `r addCite(perez2024estatistics)`
#'   `r addCite(wang2025anytime)`
#'
#' @export
#'
saviTTestStatNEffNuMom <- function(
    t, nEff, nu, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE,
    k=1,
    ...) {

  g <- parameter
  alternative <- match.arg(alternative)

  # Do t-density approximation numeric
  if (tDensity && is.finite(t)) {
    momIntegrand <- function(delta) {
      exp(
        stats::dt(t, df=nu, ncp=sqrt(nEff)*delta, log=TRUE)
        -stats::dt(t, df=nu, ncp=0, log=TRUE)
        +stats::dnorm(delta, mean=0, sd=sqrt(g), log=TRUE)
      )*delta^2/g
    }

    lowerBound <- switch(alternative,
                         "twoSided"=-Inf,
                         "greater"=0,
                         "less"=-Inf)

    upperBound <- switch(alternative,
                         "twoSided"=Inf,
                         "greater"=Inf,
                         "less"=0)

    tempResult <- stats::integrate(momIntegrand, lowerBound, upperBound)

    sidedConstant <- if (alternative=="twoSided") 1 else 2

    result <- list("eValue"=sidedConstant*tempResult[["value"]],
                   "eValueApproxError"=tempResult[["abs.error"]])
    return(result)
  }

  # For t infinite simplify some ratios
  if (is.infinite(t)) {
    logTerm <- nu/2*log(1+nEff*g)

    if (k==1) {
      hypTerm <- (1+nEff*g*(nu+1))/(1+nEff*g)
    } else {
      hypTerm <- suppressWarnings(
        try(
          hypergeo::f15.3.3(
            A=-k, B=-nu/2, C=1/2,
            z=nEff*g/(1+nEff*g))
        )
      )

      if (isTryError(hypTerm) || is.na(hypTerm))
        hypTerm <- Re(hypergeo::hypergeo(A=-k, B=-nu/2, C=1/2, z=nEff*g/(1+nEff*g)))
    }

    aPart <- exp(logTerm)*hypTerm

    if (alternative=="twoSided") {
      return(list(eValue=aPart))
    } else {
      logBPart <-
        lgamma(k+1)-lgamma(k+1/2)+lgamma(nu/2+1)-lgamma((nu+1)/2)+
        1/2*(log(nEff*g)-log(1+nEff*g))+log(2)

      hypTermB <- suppressWarnings(
        try(
          hypergeo::f15.3.3(
            A=1/2-k, B=1-(nu+1)/2, C=3/2,
            z=nEff*g/(1+nEff*g))
        )
      )

      if (isTryError(hypTermB) || is.na(hypTermB)) {
        hypTermB <- Re(hypergeo::hypergeo(A=1/2-k, B=1-(nu+1)/2, C=3/2,
                                          z=nEff*g/(1+nEff*g)))
      }

      signAlternative <- switch(alternative,
                                "greater"=1,
                                "less"=-1)

      eValueAlmost <- hypTerm+signAlternative*sign(t)*exp(logBPart)*hypTermB

      if (eValueAlmost < 0)
        stop("Overflow: eValue should be positive")

      if (eValueAlmost==0)
        eValueAlmost <- .Machine$double.xmin

      eValue <- exp(logTerm)*eValueAlmost

      if (alternative=="greater") {
        if (t < 0 && eValue > 1)
          warning("Overflow: Subtraction of two large numbers: eValue should be smaller than one")
      } else if (alternative=="less") {
        if (t > 0 && eValue > 1)
          warning("Overflow: Subtraction of two large numbers: eValue should be smaller than one")
      }

      return(list(eValue=eValue))
    }
  }

  # Normal run here
  #
  logTerm <- -(k+1/2)*log(1+nEff*g) +
    (k+(nu+1)/2)*(log((1+t^2/nu))-log(1+t^2/(nu*(1+nEff*g))))

  if (k==1) {
    hypTerm <- exp(log(1+(1+nEff*g*(nu+1))/(nu*(1+nEff*g))*t^2)-log(1+t^2/nu))
  } else {
    hypTerm <- suppressWarnings(
      try(
        hypergeo::f15.3.3(
          A=-k, B=-nu/2, C=1/2,
          z=nEff*g/(1+nEff*g)*t^2/(nu+t^2))
      )
    )

    if (isTryError(hypTerm) || is.na(hypTerm))
      hypTerm <- Re(hypergeo::hypergeo(A=-k, B=-nu/2, C=1/2, z=nEff*g/(1+nEff*g)*t^2/(nu+t^2)))
  }

  aPart <- exp(logTerm)*hypTerm

  if (alternative=="twoSided") {
    return(list(eValue=aPart))
  } else {
    logBPart <-
      lgamma(k+1)-lgamma(k+1/2)+lgamma(nu/2+1)-lgamma((nu+1)/2)+
      log(2)+1/2*(log(nEff*g)-log(1+nEff*g)-log(nu+t^2))

    hypTermB <- suppressWarnings(
      try(
        hypergeo::f15.3.3(
          A=1/2-k, B=1-(nu+1)/2, C=3/2,
          z=nEff*g/(1+nEff*g)*t^2/(nu+t^2))
      )
    )

    if (isTryError(hypTermB) || is.na(hypTermB)) {
      hypTermB <- Re(hypergeo::hypergeo(A=1/2-k, B=1-(nu+1)/2, C=3/2,
                                        z=nEff*g/(1+nEff*g)*t^2/(nu+t^2)))
    }

    signAlternative <- switch(alternative,
                              "greater"=1,
                              "less"=-1)

    eValueAlmost <- hypTerm+signAlternative*t*exp(logBPart)*hypTermB
    eValueAlmost <- round(eValueAlmost, 9)

    if (eValueAlmost < 0)
      stop("Overflow: eValue should be positive")

    if (eValueAlmost==0)
      eValueAlmost <- .Machine$double.xmin

    eValue <- exp(logTerm)*eValueAlmost

    if (alternative=="greater") {
      if (t < 0 && eValue > 1)
        warning("Overflow: Subtraction of two large numbers: eValue should be smaller than one")
    } else if (alternative=="less") {
      if (t > 0 && eValue > 1)
        warning("Overflow: Subtraction of two large numbers: eValue should be smaller than one")
    }

    return(list(eValue=eValue))
  }
}

#' SaviTTestStat based on the t-statistic, nEff and nu and the grow prior
#'
#' @rdname saviTTestStat
#' @inheritParams saviTTestStatNEffNu
#'
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'   `r addCite(perez2024estatistics)`
#'   `r addCite(wang2025anytime)`
#'
#' @export
#'
saviTTestStatNEffNuGrow <- function(
    t, nEff, nu, parameter,
    alternative=c("twoSided", "less", "greater"),
    tDensity=FALSE,
    paired=FALSE, ...) {

  # TODO(Alexander):
  #   One-sided not as stable as two-sided due to hypergeo::genhypergeo for the odd component
  #   1. Use Kummer's transform again (??)
  #   2. Switch to numerical integration. Boundary case
  #

  deltaS <- parameter

  alternative <- match.arg(alternative)

  if (isTRUE(tDensity)) {
    eValue <- saviTTestStatTDensity(
      "t"=t, "parameter"=parameter,
      "nu"=nu, "nEff"=nEff, "alternative"=alternative)
    return(list("eValue"=eValue))
  }

  if (is.finite(t)) {
    a <- t^2/(nu+t^2)
    expTerm <- exp((a-1)*nEff*deltaS^2/2)

    zArg <- (-1)*a*nEff*deltaS^2/2

    compMode <- 1

    aKummerFunction <- Re(hypergeo::genhypergeo(U=-nu/2, L=1/2, zArg))
    aKummerFunction2 <- aKummerFunction

    if (aKummerFunction > 0) {
      aPart <- expTerm * aKummerFunction
    } else if (aKummerFunction <= 0) {
      compMode <- 2

      aKummerFunction2 <- max(compute1F1AllVersions(U=-nu/2, L=1/2, z=zArg))
      aPart <- expTerm * aKummerFunction2
    }

    if (alternative=="twoSided") {
      eValue <- aPart
    } else {
      bKummerFunction <- exp(lgamma(nu/2+1)-lgamma((nu+1)/2))*sqrt(2*nEff)*deltaS*t/sqrt(t^2+nu) *
        Re(hypergeo::genhypergeo(U=(1-nu)/2, L=3/2, zArg))

      signAlt <- if (alternative=="greater") 1 else -1

      bPart <- expTerm*signAlt*bKummerFunction

      eValue <- aPart + bPart

      if (eValue <= 0) {
        signBPart <- sign(bPart)

        if (compMode==1) {
          compMode <- 2
          aKummerFunction2 <- max(compute1F1AllVersions(U=-nu/2, L=1/2, z=zArg))
        }

        bKummerFunction2 <- exp(lgamma(nu/2+1)-lgamma((nu+1)/2))*sqrt(2*nEff)*deltaS*t/sqrt(t^2+nu) *
          min(compute1F1AllVersions(U=(1-nu)/2, L=3/2, z=zArg))

        eValue <- expTerm*(aKummerFunction2 + signBPart*bKummerFunction2)
      }
    }

    if (is.null(eValue) || is.na(eValue) || eValue <= 0) {
      eValue <- saviTTestStatTDensity(
        "t"=t, "parameter"=parameter,
        "nu"=nu, "nEff"=nEff, "alternative"=alternative)
    }
  } else if (is.infinite(t)) {
    expTerm <- exp(-nEff*parameter^2/2)
    aHypTerm <- hypergeo::genhypergeo(U=(nu+1)/2, L=1/2, z=nEff*parameter^2/2)

    if (expTerm==0 && is.infinite(aHypTerm)) {
      x <- nEff*parameter^2/2

      expTerm <- exp(lgamma(1/2)-lgamma((nu+1)/2)+nu/2*log(x))
      seriesTerm <- 1 + nu^2/4*x^(-1)
      aPart <- expTerm*seriesTerm

      if (alternative=="twoSided") {
        eValue <- aPart
      } else {
        expTerm <- 2*exp(lgamma(3/2)-lgamma((nu+1)/2)+nu/2*log(x))
        seriesTerm <- 1 + nu*(nu-1)/4*x^(-1)

        bPart <- expTerm*seriesTerm

        if (alternative=="greater") {
          eValue <- aPart+sign(t)*bPart
        } else if (alternative=="less") {
          eValue <- aPart-sign(t)*bPart
        }
      }
    } else {
      if (alternative=="twoSided") {
        eValue <- expTerm*aHypTerm
      } else {
        bHypTerm <- sqrt(2*nEff)*exp(lgamma((nu+2)/2)-lgamma((nu+1)/2))*parameter*
          hypergeo::genhypergeo(U=nu/2+1, L=3/2, z=nEff*parameter^2/2)

        if (alternative=="greater") {
          eValue <- expTerm*(aHypTerm+sign(t)*bHypTerm)
        } else if (alternative=="less") {
          eValue <- expTerm*(aHypTerm-sign(t)*bHypTerm)
        }
      }
    }

    if (eValue < 0)
      eValue <- 1e-270
  }

  result <- list("eValue"=eValue)
  return(result)
}

#' saviTTestStat() based on t-densities
#'
#' This is \code{\link{saviTTestStat}()} based on t-densities instead of
#' hypergeometric functions.
#'
#' @inheritParams saviTTest
#' @inheritParams saviTTestStatNEffNu
#' @rdname saviTTestStat
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'   `r addCite(perez2024estatistics)`
#'   `r addCite(wang2025anytime)`
#'
#'
saviTTestStatTDensity <- function(t, parameter, nu, nEff,
                                  alternative=c("twoSided", "less", "greater"),
                                  paired=FALSE, ...) {
  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  deltaS <- parameter

  if (alternative=="twoSided") {
    logTerm1 <- stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)
    logTerm2 <- stats::dt(t, df=nu, ncp=-sqrt(nEff)*deltaS, log=TRUE)-stats::dt(t, df=nu, ncp=0, log=TRUE)

    term1 <- if (is.infinite(logTerm1)) 0 else exp(logTerm1)
    term2 <- if (is.infinite(logTerm2)) 0 else exp(logTerm2)

    result <- try((term1+term2)/2)
  } else {
    result <- try(exp(
      stats::dt(t, df=nu, ncp=sqrt(nEff)*deltaS, log=TRUE) -
        stats::dt(t, df=nu, ncp=0, log=TRUE)))
  }

  if (result < 0) {
    warning("Numerical overflow: E-value is essentially zero")
    return(2^(-25))
  }


  return(result)
}


#' Safe Anytime-Valid Student's T-Test.
#'
#' Savi one- and two-sample T-tests. Takes as input vector(s) of data and a designObj from
#' \code{\link{designSaviT}}. The function is modelled after \code{\link[stats]{t.test}()}.
#'
#' @aliases savi.t.test
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param designObj an object obtained from \code{\link{designSaviT}}.
#' @param paired a logical, if \code{TRUE} then pair the data.
#' @param ciValue numeric representing the confidence level.
#' Default ciValue=NULL yields ciValue = 1 - alpha
#' @param varEqual a logical variable indicating whether to treat
#' the two variances as being equal. Default varEqual=TRUE.
#' @param maxRoot Used to bound the candidate set of width of the
#' confidence interval/
#' @param formula a formula of the form lhs ~ rhs where lhs
#' is a numeric variable giving the data values and rhs
#' either 1 for a one-sample or paired test or a factor
#' with two levels giving the corresponding groups. If lhs
#' is of class "Pair" and rhs is 1, a paired test is done
#' @param data an optional matrix or data frame (or similar:
#' see \code{\link[stats]{model.frame}()}) containing the variables in
#' the formula. By default the variables are taken from
#' environment(formula).
#' @param subset an optional vector specifying a subset of
#' observations to be used.
#' @param na.action a function which indicates what should
#' happen when the data contain \code{NA}s. Defaults to
#' getOption("na.action")..
#' @param sequential a logical indicating whether a sequential
#' analysis should be performed.
#' @param tDensity Uses the the representation of the savi T-test as the likelihood ratio of t densities.
#' @param wantCi default TRUE
#' @param ... further arguments to be passed to or from methods.
#'
#' @return Returns an object of class 'saviTest'. An object of class 'saviTest'
#' is a list containing at least the following components:
#'
#' \describe{
#'   \item{statistic}{the value of the t-statistic.}
#'   \item{n}{The realised sample size(s).}
#'   \item{eValue}{the realised e-value from the savi test.}
#'   \item{confSeq}{A savi confidence interval for the mean (difference)}
#'   \item{estimate}{the estimated means or mean (difference) depending on
#'   whether it was a one-sample test or a two-sample test.}
#'   \item{stderr}{the standard error of the mean (difference), used as
#'   denominator in the t-statistic formula.}
#'   \item{dataName}{a character string giving the name(s) of the data.}
#'   \item{designObj}{an object of class "saviDesign" obtained from \code{\link{designSaviT}()}.}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'   `r addCite(perez2024estatistics)`
#'   `r addCite(wang2025anytime)`
#'
#' @export
#'
#' @examples
#' ## Examples with simulated data ----
#' set.seed(1)
#' x <- rnorm(30, mean=0)
#' y <- rnorm(40, mean=0)
#'
#' # Because no designObj is specified, a default
#' # designObj is used with deltaMin = 1/2,
#' # which can be thought of as a medium effect size
#' # according to Cohen.
#'
#' res <- saviTTest(x=x, y=y)
#'
#' # By default sequential=TRUE, because length(x) <= 200.
#' # This allows us to visualise the e-value as a function of
#' # the n1 and associated n2, where the ratio of sample sizes,
#' # ratio=n2/n1 is maintained. Here the e-value at n1=6 uses data
#' # x[1:6] and y[1:ceil(ratio*6)]
#' plot(res)
#'
#' # Plots the confidence sequence
#' plot(res, wantConfSeqPlot=TRUE)
#'
#' # See ?designSaviT for more info
#' # This designObj also allows for
#' # evidence quantification that the
#' # mean difference is minimal clinically
#' # relevant, here, larger than
#' # relevanceSize=meanDiffMin
#' designObj <- designSaviT(deltaMin=0.7, alpha=0.05,
#'                          alternative="twoSided",
#'                          testType="twoSample",
#'                          relevanceTest=TRUE)
#'
#' res <- saviTTest(x=x, y=y, designObj=designObj)
#'
#' plot(res)
#'
#' # Note that the e-value against relevance falls below alphaRelevance
#' # We can reject the hypothesis that the effect is relevantly large,
#' # larger than designObj$relevanceTestSim$parameter, after a sample size of
#' min(which(res$eRelevanceVec <= designObj$relevanceTestSim$alpha))
#'
#' set.seed(2)
#' x <- rnorm(30, mean=0.6)
#' y <- rnorm(40, mean=0)
#'
#' res <- saviTTest(x, y, designObj=designObj)
#'
#' plot(res)
#' # We could have stopped sampling after
#' min(which(res$eValueVec >= 1/designObj$alpha))
#' # the yellow curve crosses 1/alpha sooner,
#' # but we do **not** compare eRelevance >= 1/alpha, only eRelevance <= alphaRelevance.
#'
#' ## Classical example: Student's sleep data -----
#' plot(extra ~ group, data = sleep)
#'
#' designObj <- designSaviT(deltaMin=0.6, testType="twoSample")
#'
#' ## Traditional interface
#' with(sleep, saviTTest(extra[group == 1], extra[group == 2],
#'                       designObj=designObj))
#'
#' ## Formula interface
#' saviTTest(extra ~ group, data = sleep, designObj=designObj)
#'
#' ## Formula interface to one-sample test
#' designObj1 <- designSaviT(deltaMin=0.6,
#'                           testType="oneSample",
#'                           sigma=2)
#'
#' saviTTest(extra ~ 1, data = sleep, designObj=designObj1)
#'
#' ## Formula interface to paired test
#' ## The sleep data are actually paired, so could have been in wide format:
#' designObjPaired <- designSaviT(deltaMin=0.6,
#'                                testType="paired",
#'                                sigma=1.4)
#' sleep2 <- reshape(sleep, direction = "wide",
#'                   idvar = "ID", timevar = "group")
#' saviTTest(Pair(extra.1, extra.2) ~ 1, data = sleep2,
#'           designObj=designObjPaired)
saviTTest <- function(x, ...) {
  UseMethod("saviTTest")
}


#' @rdname saviTTest
#' @param nuMin numeric > 0, the minimum degrees of freedom under which the results
#' are trivial, thus, 1.
#' @export
saviTTest.default <- function(
    x, y=NULL, designObj=NULL, paired=FALSE,
    varEqual=TRUE, ciValue=NULL,
    maxRoot=10, sequential=NULL,
    tDensity=FALSE, nuMin=2, wantCi=TRUE, ...) {

  result <- constructSaviTestObj("T-Test")

  eRelevance <- NULL

  if (is.null(varEqual))
    varEqual <- designObj[["varEqual"]]

  # Vars for sequential analysis
  eValueVec <- NULL
  eRelevanceVec <- NULL
  confSeqMatrix <- NULL
  n1Vec <- NULL
  n2Vec <- NULL

  fpt <- NULL
  fptRelevance <- NULL

  ## Def: test type -------

  if (is.null(y)) {
    testType <- "oneSample"

    if (paired)
      stop("Data error: Paired analysis requested without specifying the second variable")

    dataName <- deparse1(substitute(x))
  } else {
    testType <- if (paired) "paired" else "twoSample"

    dataName <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  }

  ## Check: designObj ----
  if (is.null(designObj)) {
    designObj <- designSaviT(0.5, "eType"="mom",
                             "testType"=testType)
    designObj[["pilot"]] <- TRUE

    warningMessage <- paste("No designObj given. Default test computed based",
                            "on a non-local prior at +1/2 and -1/2.")
    warning(warningMessage)
  }

  if (designObj[["testName"]] != "T-Test")
    warning("The provided design is not constructed for the t-test,",
            "please use designSaviT() instead. The test results might be invalid.")

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  ## Check: Data -----
  #
  sumStats <- computeZTSumStats(
    "x"=x, "y"=y, "sequential"=sequential,
    "varEqual"=varEqual, "paired"=paired,
    "testType"=testType)

  # Dummies that get filled by list2env
  nu <- sdObs <- nEff <- meanObs <- nEffVec <- meanObsVec <-
    sdObsVec <- nuVec <- estimate <- n <- nEff <- meanObs <- n1 <-
    n2 <- nEffVec <- meanObsVec <- estimate <- n <- sdObsVec <- nuVec <- NULL

  list2env(sumStats, envir=environment())

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]

  if (is.null(ciValue))
    ciValue <- 1-alpha

  if (wantCi && ciValue < 0 || ciValue > 1)
    stop("Can't make a confidence sequence with ciValue < 0 or ciValue > 1, or alpha < 0 or alpha > 1")

  if (nu <= nuMin || is.na(sdObs)) {
    tStat <- 0
  } else {
    tStat <- try(sqrt(nEff)*(meanObs - h0)/sdObs)
  }

  # tStat <- if (nu <= nuMin) 0 else tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sdObs)
  #
  # if (meanObs==0 && sdObs==0) {
  #   tStat <- 0
  # } else {
  #   tStat <- if (nu <= nuMin) 0 else tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sdObs)
  # }

  if (is.na(tStat) && sdObs==0 && meanObs-h0==0)
    tStat <- 0

  if (is.na(tStat))
    stop("Data error: Could not compute the t-statistic")

  names(tStat) <- "t"

  ### Compute: eValue ----
  #
  testResult <- suppressWarnings(
    saviTTestStatNEffNu("t"=tStat, nEff=nEff, nu=nu,
                        "parameter"=designObj[["parameter"]],
                        "alternative"=alternative, "paired"=paired,
                        "tDensity"=tDensity,
                        "nuMin"=nuMin, "eType"=designObj[["eType"]])
  )

  if (designObj[["relevanceTest"]]) {
    relevanceRes <- suppressWarnings(
      saviRelevanceTStatNEffNu("t"=tStat, "nEff"=nEff, "nu"=nu,
                               "parameter"=designObj[["relevanceTestSim"]][["parameter"]],
                               "alternative"=designObj[["alternative"]], "paired"=paired,
                               "nuMin"=nuMin)
    )
    eRelevance <- unname(relevanceRes[["eValue"]])
  }

  ### Compute: confSeq ----
  #
  if (wantCi) {
    result[["confSeq"]] <- computeConfidenceIntervalT(
      "meanObs"=meanObs, "sdObs"=sdObs,
      "nEff"=nEff, "nu"=nu,
      "parameter"=designObj[["parameter"]],
      "eType"=designObj[["eType"]], "ciValue"=ciValue, "maxRoot"=maxRoot)
  }

  ## Compute: Sequential ----
  #
  # TODO(Alexander):
  #
  if (sequential) {

    tStatVec <- sqrt(nEffVec)*(meanObsVec-h0)/sdObsVec

    mIter <- length(n1Vec)

    eRelevanceVec <- eValueVec <- numeric(mIter)
    confSeqMatrix <- matrix(nrow=mIter, ncol=2)

    for (i in seq_along(n1Vec)) {
      res <- suppressWarnings(
        saviTTestStatNEffNu(
          "t"=tStatVec[i], "nEff"=nEffVec[i], "nu"=nuVec[i],
          "parameter"=designObj[["parameter"]], "alternative"=alternative,
          "paired"=paired, "tDensity"=tDensity, "nuMin"=nuMin,
          "eType"=designObj[["eType"]]
        )
      )

      eValueVec[i] <- unname(res[["eValue"]])

      if (designObj[["relevanceTest"]]) {
        relevanceRes <- suppressWarnings(
          saviRelevanceTStatNEffNu("t"=tStatVec[i],
                                   "nEff"=nEffVec[i], "nu"=nuVec[i],
                                   "parameter"=designObj[["relevanceTestSim"]][["parameter"]],
                                   "alternative"=designObj[["alternative"]], "paired"=paired,
                                   "nuMin"=nuMin)
        )
        eRelevanceVec[i] <- unname(relevanceRes[["eValue"]])
      }

      if (wantCi) {
        kaas <- computeConfidenceIntervalT("meanObs"=meanObsVec[i], "sdObs"=sdObsVec[i],
                                           "nEff"=nEffVec[i], "nu"=nuVec[i],
                                           "parameter"=designObj[["parameter"]],
                                           "eType"=designObj[["eType"]], "ciValue"=ciValue,
                                           "maxRoot"=maxRoot)

        confSeqMatrix[i, ] <- kaas
      }
    }

    tempConfSeq <- c(max(confSeqMatrix[, 1]), min(confSeqMatrix[, 2]))

    if (tempConfSeq[1] >= tempConfSeq[2]) {
      warning("Possible high degree of heterogeneity",
              "leading to an empty running intersection confidence sequence")
    } else if (tempConfSeq[1] < tempConfSeq[2]) {
      result[["confSeq"]] <- tempConfSeq
    }

    fpt <- suppressWarnings(
      min(which(eValueVec >= 1/designObj[["alpha"]]))
    )

    if (designObj[["relevanceTest"]])
      fptRelevance <- suppressWarnings(
        min(which(eRelevanceVec <= designObj[["relevanceTestSim"]][["alpha"]]))
      )
  }

  ### Fill: Result -----
  #
  result[["statistic"]] <- tStat
  result[["estimate"]] <- estimate
  result[["stderr"]] <- sdObs/sqrt(nEff)
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["testType"]] <- testType
  result[["n"]] <- n
  result[["ciValue"]] <- ciValue

  result[["eRelevance"]] <- eRelevance

  result[["eValueVec"]] <- eValueVec
  result[["eRelevanceVec"]] <- eRelevanceVec
  result[["confSeqMatrix"]] <- confSeqMatrix
  result[["n1Vec"]] <- n1Vec
  result[["n2Vec"]] <- n2Vec

  result[["eValue"]] <- testResult[["eValue"]]
  result[["eValueApproxError"]] <- testResult[["eValueApproxError"]]

  # Sumstats
  result[["meanObs"]] <- meanObs
  result[["sdObs"]] <- sdObs
  result[["meanObsVec"]] <- meanObsVec
  result[["sdObsVec"]] <- sdObsVec
  result[["nEffVec"]] <- nEffVec
  result[["nuVec"]] <- nuVec
  result[["fpt"]] <- fpt
  result[["fptRelevance"]] <- fptRelevance

  return(result)
}


#' @rdname saviTTest
#' @aliases saviTTest
#' @export
#'
saviTTest.formula <- function(
    formula, data, subset, na.action, ...) {

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  wantTwoSample <- TRUE

  if (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L)
    if (formula[[3L]] == 1L)
      wantTwoSample <- FALSE
  else
    stop("'formula' missing or incorrect")

  matchedCall <- match.call(expand.dots = FALSE)

  if (is.matrix(eval(matchedCall[["data"]], parent.frame())))
    matchedCall[["data"]] <- as.data.frame(data)

  # Note: Prepare calling stats::model.frame instead of saviTTest
  #
  matchedCall[[1L]] <- quote(stats::model.frame)
  matchedCall[["..."]] <- NULL

  # Call: stats::model.frame
  #
  modelFrame <- eval(matchedCall,
                     parent.frame())

  # Naming
  dataName <- paste(names(modelFrame),
                    collapse=" by ")

  names(modelFrame) <- NULL
  response <- attr(attr(modelFrame, "terms"),
                   "response")

  if (isTRUE(wantTwoSample)) {
    groupingFactor <- factor(modelFrame[[-response]])

    if (nlevels(groupingFactor) != 2L)
      stop("grouping factor must have exactly 2 levels")

    dataList <- split(modelFrame[[response]], groupingFactor)

    tResult <- saviTTest("x"=dataList[[1L]], "y"=dataList[[2L]], ...)

    if (length(tResult[["estimate"]]) == 2L) {
      names(tResult[["estimate"]]) <- paste("mean in group", levels(groupingFactor))
      names(tResult[["designObj"]][["h0"]]) <-
        paste("true difference in means between",
              paste("group", levels(groupingFactor), collapse = " and "))
    }
  } else {
    respVar <- modelFrame[[response]]

    if (inherits(respVar, "Pair")) {
      tResult <- saviTTest("x"=respVar[, 1L], "y"=respVar[, 2L],
                           paired=TRUE, ...)
      firstVar <- substring(dataName,
                            first=6,
                            last=regexpr(",", dataName)-1)
      secondVar <- substring(dataName,
                             first=regexpr(",", dataName)+2,
                             last=regexpr(")", dataName)-1)
      names(tResult[["estimate"]]) <-
        paste("mean difference between", firstVar, "and", secondVar)
      names(tResult[["designObj"]][["h0"]]) <-
        paste("true mean difference between",
              paste(c(firstVar, secondVar), collapse = " and "))
    } else {
      tResult <- saviTTest("x"=respVar, "y"=NULL, ...)
    }
  }

  tResult[["dataName"]] <- dataName
  return(tResult)
}


#' @rdname saviTTest
#' @aliases saviTTest
#' @export
savi.t.test <- function(x, y=NULL, paired=FALSE, designObj=NULL, varEqual=TRUE,
                        ciValue=NULL, ...) {
  result <- saviTTest("x"=x, "y"=y, "designObj"=designObj,
                      "paired"=paired, "varEqual"=varEqual,
                      ...)

  argumentNames <- getArgs()
  xLabel <- extractNameFromArgs(argumentNames, "x")

  if (is.null(y)) {
    dataName <- xLabel
  } else {
    yLabel <- extractNameFromArgs(argumentNames, "y")
    dataName <- paste(xLabel, "and", yLabel)
  }

  result[["dataName"]] <- dataName
  return(result)
}


#' Computes E-Relevance Values Based on the T-Statistic in favour of the alternative over practical equivalence
#'
#' Evidence for practical equivalence requires the e-value to be small, i.e.
#' smaller than alphaRelevance.
#' If the alternative holds true, i.e. deltaTrue >= relevanceSize, then
#' there is no more than alphaRelevance probability of ever seeing
#' eRelevance <= alphaRelevance.
#'
#' @rdname saviTTestStat
#' @inheritParams designSaviT
#'
#' @export
#'
#' @examples
#' # evidence for the alternative over minimal efficacy
#' saviRelevanceTStatNEffNu(t=3, nEff=100, nu=60, parameter=0.4)
#' # evidence for minimal efficacy over the alternative
#' saviRelevanceTStatNEffNu(t=0.35, nEff=100, nu=60, parameter=0.4)
saviRelevanceTStatNEffNu <- function(
    t, nEff, nu, parameter=NULL,
    alternative = c("twoSided", "less", "greater"), eType="grow",
    tDensity = FALSE, paired = FALSE, relevanceSize, nuMin=2, ...) {
  # Note overflow for t big
  #
  #     saviRelevanceTStatNEffNu(t=40.017, nEff=10000, nu=6000, parameter=0.4)
  #     saviRelevanceTStatNEffNu(t=40.018, nEff=10000, nu=6000, parameter=0.4)


  if (is.null(parameter) && is.null(relevanceSize))
    stop("No parameter, nor minimal clinically relevant effect size given")

  if (is.null(parameter) && !is.null(relevanceSize))
    parameter <- abs(relevanceSize)

  if (nu < nuMin)
    return(list("eValue"=1))

  alternative <- match.arg(alternative)

  if (alternative %in% c("twoSided", "greater"))
    sPlus0 <- suppressWarnings(
      saviTTestStatNEffNuGrow(
        "t"=t, "nEff"=nEff, "nu"=nu, "parameter"=parameter,
        "alternative"="greater", "paired"=paired)[["eValue"]]
    )

  if (alternative %in% c("twoSided", "less"))
    sMin0 <- suppressWarnings(
      saviTTestStatNEffNuGrow(
        "t"=t, "nEff"=nEff, "nu"=nu, "parameter"=parameter,
        "alternative"="less", "paired"=paired)[["eValue"]]
    )

  eValue <- switch(alternative,
                   "twoSided"=max(sPlus0, sMin0),
                   "greater"=sPlus0,
                   "less"=sMin0)

  if (eValue < 0) {
    warning("Overflow: e-value smaller than 0")
    eValue <- 1e-270
  }

  result <- list("eValue"=eValue)

  return(result)
}


#' Helper function: Computes the savi confidence sequence for the mean in a T-test
#'
#' @inheritParams saviTTestStatNEffNu
#' @inheritParams saviTTest
#' @inheritParams designSaviT
#'
#' @param meanObs numeric, the observed mean. For two sample tests this
#' is difference of the means.
#' @param sdObs numeric, the observed standard deviation. For a two-sample
#' test this is the root of the pooled variance.
#'
#' @return numeric vector that contains the upper and lower bound of the savi confidence sequence
#' @export
#'
#' @examples
#' computeConfidenceIntervalT(meanObs=0.3, sdObs=2, nEff=12, nu=11, parameter=0.4)
computeConfidenceIntervalT <- function(
    meanObs, sdObs, nEff, nu, parameter,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    alternative=c("twoSided", "greater", "less"),
    ciValue=0.95, maxRoot=11, nuMin=2) {

  eType <- match.arg(eType)
  alternative <- match.arg(alternative)

  trivialConfInt <- switch(alternative,
                           "twoSided"=c(-Inf, Inf),
                           "greater"=c(0, Inf),
                           "less"=c(-Inf, 0))

  g <- parameter

  if (nu <= 1) return(trivialConfInt)

  if (sdObs==0 && nu <= nuMin) return(trivialConfInt)

  alpha <- 1-ciValue

  numeratorW <- nu*(((1+nEff*g)/alpha^2)^(1/(nu+1))-1)
  denominatorW <- 1-((1+nEff*g)/alpha^2)^(1/(nu+1))/(1+nEff*g)

  W <- numeratorW/denominatorW

  if (eType=="eGauss" && alternative=="twoSided") {
    if (W < 0) return(trivialConfInt)

    width <- sdObs/sqrt(nEff)*sqrt(W)
  } else {

    if (alternative=="twoSided") {
      ciLogPenaltyFunc <- function(ciValue) 1/(1-ciValue)
    } else if (alternative %in% c("greater", "less")) {
      ciLogPenaltyFunc <- function(ciValue) 1/(2*(1-ciValue))
    }

    targetFunction <- function(t) {
      saviTTestStatNEffNu("t"=t, "nEff"=nEff, "nu"=nu, "parameter"=parameter,
                          "eType"=eType)$eValue-ciLogPenaltyFunc(ciValue)
    }

    if (W > 0) {
      lowerB <- max(sqrt(W)-5,0)
      upperB <- max(sqrt(W)+5, 5)
    } else {
      lowerB <- 0
      upperB <- maxRoot
      maxRoot <- 2*maxRoot
    }

    targetFunction <- Vectorize(targetFunction)

    tempResult <- suppressWarnings(
      tryCatch(stats::uniroot(targetFunction, c(lowerB, upperB)),
               error=identity)
    )

    iterationN <- 1

    while ( (inherits(tempResult, "simpleError") || is.null(tempResult)) && iterationN <= 15) {
      iterationN <- iterationN+1

      tempResult <- suppressWarnings(
        tryCatch(stats::uniroot(targetFunction, c(0, maxRoot)),
                 error=identity)
      )

      maxRoot <- maxRoot*2
    }

    if (inherits(tempResult, "simpleError")) {
      return(trivialConfInt)
      # stop("Can't compute the width of the interval")
    }

    width <- sdObs/sqrt(nEff)*tempResult$root
  }


  if (alternative=="twoSided") {
    lowerCS <- meanObs - width
    upperCS <- meanObs + width
  } else if (alternative=="greater") {
    lowerCS <- meanObs + width
    upperCS <- Inf
  } else if (alternative=="less") {
    lowerCS <- -Inf
    upperCS <- meanObs - width
  }

  return(unname(c(lowerCS, upperCS)))
}


# Design fnts -------

#' Design a Frequentist T-Test
#'
#' Computes the number of samples necessary to reach a tolerable type I and
#' desired power for the frequentist T-test.
#'
#' @inheritParams designSaviT
#'
#' @return Returns an object of class 'freqTDesign'. An object of class
#' 'freqTDesign' is a list containing at least the following components:
#' \describe{
#'   \item{nPlan}{the planned sample size(s).}
#'   \item{esMin}{the minimal clinically relevant standardised effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
#'   \item{power}{the desired power provided by the user.}
#'   \item{lowN}{the smallest n of the search space for n provided by the user.}
#'   \item{highN}{the largest n of the search space for n provided by the user.}
#'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
#'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
#' }
#' @export
#'
#' @examples
#' designFreqT(0.5)
designFreqT <- function(deltaMin, alpha=0.05, power=0.8,
                        alternative=c("twoSided", "greater", "less"),
                        h0=0, testType=c("oneSample", "paired", "twoSample"), ...) {
  stopifnot(alpha > 0, power > 0, power < 1, alpha < 1, power < 1)

  beta <- 1-power

  testType <- match.arg(testType)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)

  alternativeFreq <- switch(alternative,
                            "greater"="one.sided",
                            "less"="one.sided",
                            "twoSided"="two.sided")

  testTypeFreq <- switch(testType,
                         "twoSample"="two.sample",
                         "oneSample"="one.sample",
                         "paired"="paired")

  tempResult <- stats::power.t.test("delta"=deltaMin, "power"=1-beta, "type"=testTypeFreq,
                                    "alternative"=alternativeFreq)

  n1Plan <- ceil(tempResult[["n"]])
  n2Plan <- NULL

  if (testType!="oneSample") n2Plan <- n1Plan

  if (is.null(n2Plan)) {
    nPlan <- n1Plan
    names(nPlan) <- "n1Plan"
  } else {
    nPlan <- c(n1Plan, n2Plan)
    names(nPlan) <- c("n1Plan", "n2Plan")
  }

  result <- list(nPlan=nPlan, "esMin"=deltaMin, "alpha"=alpha, "power"=1-beta,
                 "testType"=testType, "alternative"=alternative, "ratio"=1, "h0"=h0)
  class(result) <- "freqTDesign"
  return(result)
}

#' Design a Safe Anytime-Valid Experiment to Test Means with a Z Test
#'
#' A designed experiment requires (1) a sample size nPlan to plan for, and
#' (2) a savi test defining parameter. The design involves alpha and the
#' three quantities: (1) nPlan, (2) power, and (3) a minimal
#' clinically relevant standarised mean difference deltaMin.
#' \describe{
#'   \item{Scenario 1.a}{Goal: "nPlan" and optimal E-variable. Given: deltaMin and power.}
#'   \item{Scenario 1.b}{Goal: an optimal E-variable. Given: deltaMin only.}
#'   \item{Scenario 2}{Goal: "power" and optimal E-variable. Given: deltaMin and nPlan.}
#'   \item{Scenario 3.a}{Goal: "deltaMin" and optimal E-variable. Given: power and nPlan.}
#'   \item{Scenario 3.b}{Goal: an optimal E-variable. Given: nPlan only.}
#'}
#'
#' Every scenario returns an E-variable adapted to the input. Scenario 1.a,
#' for instance, outputs the parameter of the provided eType (default mom)
#' savi test, see \code{\link{matchEParameterWith}} for details, and nPlan.
#' The nPlan is based on samples paths drawn under deltaTrue (if not specified,
#' then deltaTrue=deltaMin by default). The resulting nPlan corresponds to the
#' power (say 80%) quantile of the first-passage time distribution associated
#' with E crossing threshold 1/alpha.
#'
#' @param deltaMin numeric that defines the minimal relevant standardised mean difference,
#' the smallest population effect size that we would like to detect (with sufficient power).
#' @param power numeric in (0, 1) that specifies the desired power, that is, the targetted
#' chance to stop in favour of the alternative over the null hypothesis, when the alternative
#' holds true. Note that prior to version 0.8.8 power <- 1-beta. The "beta" argument does not
#' need to be specified anymore.
#' @param nPlan optional numeric vector of length at most 2, see scenario 2 and 3 above.
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error and the null
#' rejection rule e >= 1/alpha.
#' @param h0 numeric, representing the null value, default h0=0.
#' @param alternative a character string specifying the alternative hypothesis. Must be one
#' of "twoSided" (default), "greater" or "less".
#' @param sigma numeric > 0 representing the population standard deviation used for the test.
#' @param sigma2 numeric > 0 representing the population standard deviation used for the test,
#' for the second group in a two-sample t-test
#' @param deltaTrue numeric, data governing effect size used for simulations. Default
#' deltaTrue=deltaMin.
#' @param beta numerical in (0,1). Old parameter now replaced by the power parameter
#' @param testType either one of "oneSample", "paired", "twoSample".
#' @param ratio numeric > 0 representing the randomisation ratio of condition 2 over condition 1.
#' If testType is not equal to "twoSample", or if nPlan is of length(1) then ratio=1.
#' @param parameter numeric, an optional savi test defining parameter. Default set to \code{NULL}.
#' and adapts to meanDiffMin and eType, see \code{\link{matchEParameterWith}} for details.
#' @param eType character one of "mom", "grow", "eGauss", and "eCauchy". "mom" is default
#' and uses a non-local moment prior with bump(s) at meanDiffMin, "grow" uses point prior(s) at
#' meanDiffMin, "eGauss" a zero-centred normal prior, "eCauchy" a zero centred Cauchy prior.
#' @param wantSamplePaths logical, if \code{TRUE} then also outputs the sample paths.
#' @param wantSimData logical, if \code{TRUE} then also output the simulated data
#' @param lowEsTrue numeric, lower bound for the candidate set of the
#' targeted minimal clinically relevant effect size for scenario 3.a.
#' @param highEsTrue numeric, upper bound for the candidate set of the
#' targeted minimal clinically relevant effect size for scenario 3.a.
#' @param pb logical, if \code{TRUE}, then show progress bar.
#' @param seed integer, seed number.
#' @param nSim integer > 0, the number of simulations needed to compute power or the number of
#' samples paths for the savi t test under continuous monitoring.
#' @param nBoot integer > 0 representing the number of bootstrap samples to assess the accuracy
#' of the approximations of the power, the number of samples for the savi t test under continuous
#' monitoring,or for the computation of the logarithm of the implied target.
#' @param varEqual a logical variable indicating whether to treat
#' the two variances as being equal. Default varEqual=TRUE.
#' @param relevanceTest logical, if \code{TRUE} then impose rule to stop
#' for minimal efficiency if e <= alphaRelevance. Default \code{FALSE}.
#' @param relevanceTest logical, if \code{TRUE} then impose a rule to stop
#' for minimal efficiency if e <= alphaRelevance. Default \code{FALSE}.
#' @param relevanceSize numeric, the minimal clinical relevant mean
#' difference that we do not want to miss under the alternative.
#' Default relevanceSize=NULL implies relevanceSize=abs(meanDiffMin)
#' @param alphaRelevance numeric, the threshold for relevance test. Taken to be minimum of
#' alpha and 1-power.
#' @param betaDefault numeric, defaulting value for 1-power and alphaRelevance
#' @param highN integer, largest possible sampling horizon. This might be the
#' largest n that we are able to fund, which by default is set to 1e4L.
#' Typically, highN is not used, as the function
#' `computeNPlanBatchSaviT()` tries to find the sampling horizon.
#' If all fails, then use highN as the sampling horizon.
#' @param wantSampling logical, default TRUE so sampling paths are drawn.
#' For instance, if meanDiffMin and power, are given, then nPlan
#' (scenario 1a) is derived by sampling. Set this to FALSE, whenever we
#' want to run a minimal efficacy test without needing to know nPlan
#' @param nuMin numeric > 0, the minimum degrees of freedom under which the results
#' are trivial, thus, 1.
#' @param ... further arguments to be passed to or from methods, but mainly to perform do.calls.
#'
#' @return Returns an object of class 'saviDesign'. An object of class 'saviDesign' is a list containing at least the
#' following components:
#'
#' \describe{
# #'   \item{nPlan}{the planned sample size(s).}
#'   \item{parameter}{the savi test defining parameter, see \code{\link{matchEParameterWith}}.}
# #'   \item{esMin}{the minimal clinically relevant standardised effect size provided by the user.}
#'   \item{alpha}{the tolerable type I error provided by the user.}
# #'   \item{power}{the desired power provided by the user.}
# #'   \item{alternative}{any of "twoSided", "greater", "less" provided by the user.}
# #'   \item{testType}{any of "oneSample", "paired", "twoSample" provided by the user.}
# #'   \item{paired}{logical, \code{TRUE} if "paired", \code{FALSE} otherwise.}
# #'   \item{h0}{the specified hypothesised value of the mean or mean difference depending on
# #'   whether it was a one-sample or a two-sample test.}
# #'   \item{ratio}{default is 1. Different from 1, whenever testType equals "twoSample", then it defines
# #'   ratio between the planned randomisation of condition 2 over condition 1.}
#'   \item{pilot}{logical, specifying whether it's a pilot design, which occurs
#'   when saviTTest is called without a designObj.}
#'   \item{testName}{"T-Test".}
#'   \item{call}{the expression with which this function is called.}
#' }
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' # Scenario 1.b: Goal: an E-variable
#' designObj <- designSaviT(deltaMin=0.8)
#'
#' # Scenario 1.a: Goal: "nPlan" and optimal E-variable.
#' designObj <- designSaviT(deltaMin=0.8, power=0.6, alpha=0.2,
#'                          alternative="greater", nSim=100)
#'
#' plot(designObj)
#'
#' # Scenario 1a. with relevance testing, also stopping for practically null
#' designObj <- designSaviT(deltaMin=0.8, power=0.6, alpha=0.2,
#'                          alternative="greater", nSim=100,
#'                          relevanceTest=TRUE)
#'
#' plot(designObj)
#'
#' # Scenario 2: Goal: "power" and optimal E-variable
#' designObj <- designSaviT(deltaMin=0.8, nPlan=16, nSim=100)
#'
#' # Scenario 3.a: Goal: "meanDiffMin" and optimal E-variable
#' designObj <- designSaviT(power=0.7, nPlan=16)
#'
#' # Scenario 3.b: Goal: an optimal E-variable. Given: nPlan only.
#' designObj <- designSaviT(nPlan=16)
#'
designSaviT <- function(
    deltaMin=NULL, power=NULL, nPlan=NULL,
    alpha=0.05, h0=0, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, beta=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    deltaTrue=NULL, sigma=1, sigma2=sigma,
    lowEsTrue=0.01, highEsTrue=3, varEqual=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    relevanceTest=FALSE, relevanceSize=NULL,
    alphaRelevance=NULL, betaDefault=0.2,
    highN=1e4L, wantSampling=TRUE, nuMin=2, ...) {

  stopifnot(alpha > 0, alpha < 1)

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  result <- constructSaviDesignObj("T-Test")

  if (!is.null(parameter)) {
    if (eType=="grow") {
      parameter <- checkAndReturnEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="deltaS",
        "alternative"=alternative)
    } else if (eType %in% "eGauss") {
      parameter <- checkAndReturnEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="g",
        "alternative"=alternative)
    } else if (eType=="eCauchy") {
      parameter <- checkAndReturnEsMinParameterSide(
        "paramToCheck"=parameter, "esMinName"="kappaG",
        "alternative"=alternative)
    }
  }

  if (!is.null(deltaMin)) {
    deltaMin <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=deltaMin, "esMinName"="deltaMin",
      "alternative"=alternative)

    if (is.null(deltaTrue))
      deltaTrue <- deltaMin

    parameter <- matchEParameterWith(
      "parameter"=parameter, "analysisType"="t",
      "esMin"=deltaMin, "alternative"=alternative,
      "eType"=eType)

  }

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (is.null(beta) && !is.null(power))
    beta <- 1-power

  if (relevanceTest) {
    relevanceSize <- matchRelevanceParameterWith(
      "relevanceSize"=relevanceSize, "esMin"=deltaMin,
      "esTrue"=deltaTrue)

    if (is.null(relevanceSize))
      stop("Can't run a minimal efficacy analysis without relevanceSize or deltaMin")

    alphaRelevance <- matchAlphaRelevanceWith(
      "alphaRelevance"=alphaRelevance, "alpha"=alpha,
      "power"=power, "beta"=beta, "betaDefault"=betaDefault)
  }

  designScenario <- NULL

  tempResult <- list()

  if (!is.null(deltaMin) && !is.null(power) && is.null(nPlan) && wantSampling) {
    designScenario <- "1a"

    tempResult <- designSaviT1aWantNPlan(
      "deltaMin"=deltaMin, "power"=power, "deltaTrue"=deltaTrue,
      "alpha"=alpha, "alternative"=alternative,
      "ratio"=ratio, "parameter"=parameter, "testType"=testType,
      "eType"=eType, "wantSamplePaths"=wantSamplePaths,
      "wantSimData"=wantSimData, "pb"=pb, "seed"=seed, "nSim"=nSim,
      "nBoot"=nBoot, "relevanceTest"=relevanceTest,
      "sigma"=sigma, "sigma2"=sigma2,
      "relevanceSize"=relevanceSize, "alphaRelevance"=alphaRelevance,
      "nuMin"=nuMin, ...)

  } else if (!is.null(deltaMin) && !is.null(power) && is.null(nPlan) && isFALSE(wantSampling) ||
             !is.null(deltaMin) && is.null(power) && is.null(nPlan)) {
    designScenario <- "1b"

    tempResult <- list("parameter"=parameter,
                       "esMin"=deltaMin, "relevanceTest"=relevanceTest)

    if (relevanceTest) {
      relevanceParameter <- matchRelevanceParameterWith(
        "relevanceSize"=relevanceSize,
        "esMin"=deltaMin, "esTrue"=deltaTrue)

      relevanceTestSim <- list("parameter"=relevanceParameter,
                             "alpha"=alphaRelevance)

      tempResult[["relevanceTestSim"]] <- relevanceTestSim
    }

  } else if (!is.null(deltaMin) && is.null(power) && !is.null(nPlan)) {
    # scenario 2: given effect size and nPlan, calculate power and implied target
    designScenario <- "2"

    tempResult <- designSaviT2WantPower(
      "deltaTrue"=deltaTrue, "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative, "testType"=testType,
      "ratio"=ratio, "parameter"=parameter, "eType"=eType,
      "wantSamplePaths"=wantSamplePaths, "deltaMin"=deltaMin,
      "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot,
      #
      #
      #   TODO(Alexander): Check if necessary
      #
      #
      "sigma"=sigma, "sigma2"=sigma2,
      "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
      "alphaRelevance"=alphaRelevance, "nuMin"=nuMin, ...)
  } else if (is.null(deltaMin) && !is.null(power) && !is.null(nPlan)) {
    designScenario <- "3"

    tempResult <- designSaviT3WantEsMin(
      "power"=power, "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative, "testType"=testType,
      "parameter"=parameter, "eType"=eType,
      #
      #
      #   TODO(Alexander): Check if necessary
      #
      #
      "lowEsTrue"=lowEsTrue, "highEsTrue"=highEsTrue,
      ...)
  } else if (is.null(deltaMin) && is.null(power) && !is.null(nPlan)) {
    #scenario 3b: only nPlan known, find the parameter at which the confidence interval
    # is the most narrow at nPlan

    designScenario <- "3b"

    tempResult <- designSaviT3bWantParameter(
      "nPlan"=nPlan, "alpha"=alpha,
      "alternative"=alternative, "testType"=testType,
      #
      #
      #   TODO(Alexander): Check if necessary
      #
      #
      "parameter"=parameter, "eType"=eType, ...)
  }

  if (is.null(designScenario)) {
    stop("Can't design: Please provide this function with either: \n",
         "(1.a) non-null deltaMin, non-null power and NULL nPlan, or \n",
         "(1.b) non-null deltaMin, NULL power, and NULL nPlan, or \n",
         "(1.c) NULL deltaMin, NULL power, non-null nPlan, or \n",
         "(2) non-null deltaMin, NULL power and non-null nPlan, or \n",
         "(3) NULL deltaMin, non-null power, and non-null nPlan.")
  }

  # Fill and name ----
  result <- utils::modifyList(result, tempResult)

  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["testType"]] <- testType
  result[["ratio"]] <- ratio
  result[["eType"]] <- eType
  result[["designScenario"]] <- designScenario

  ## Name esMin ----
  esMin <- result[["esMin"]]

  if (is.na(esMin))
    esMin <- NULL

  if (!is.null(esMin))
    names(esMin) <- "standardised mean difference"

  result[["esMin"]] <- esMin

  if (!is.null(deltaTrue))
    result[["esTrue"]] <- deltaTrue

  ## Name nPlan ----
  nPlan <- result[["nPlan"]]

  if (!is.null(nPlan)) {
    if (designScenario %in% c(2:3, "3b")) {
      n2Plan <- nPlan[2]

      names(nPlan) <- if (is.na(n2Plan)) "n1Plan" else c("n1Plan", "n2Plan")
    }
  }

  result[["nPlan"]] <- nPlan

  ## Name parameter ----
  parameter <- result[["parameter"]]

  if (!is.null(parameter) && is.null(names(parameter))) {
    names(parameter) <- switch(eType,
                               "mom"="gMom",
                               "eGauss"="g",
                               "imom"="tau",
                               "eCauchy"="kappaG",
                               "grow"="deltaS")
  }

  result[["parameter"]] <- parameter

  ## Name h0 -----
  names(h0) <- "mu"
  result[["h0"]] <- h0
  result[["relevanceTest"]] <- relevanceTest

  result[["varEqual"]] <- varEqual

  result[["call"]] <- sys.call()

  result <- Filter(Negate(is.null), result)
  class(result) <- "saviDesign"

  return(result)
}


#' Helper function to designing a savi T-test (output nPlan)
#'
#' Finds the parameter and power when provided with only alpha, esMin, and nPlan
#'
#' @inheritParams designSaviT
#'
#' @return A list with the parameter and the targeted nPlan amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviT1aWantNPlan(deltaMin=0.9, power=0.7, nSim=10)
designSaviT1aWantNPlan <- function(
    deltaMin, power, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL, deltaTrue=NULL, beta=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    relevanceTest=FALSE, relevanceSize=NULL,
    sigma=1, sigma2=1,
    alphaRelevance=NULL, nuMin=2, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  if (is.null(deltaTrue)) deltaTrue <- deltaMin

  power <- matchPowerWith("power"=power, "beta"=beta)

  samplingResult <- computeNPlanSaviT(
    "deltaTrue"=deltaTrue, "power"=power, "alpha"=alpha,
    "alternative"=alternative, "ratio"=ratio,
    "parameter"=parameter, "testType"=testType, "eType"=eType,
    "wantSamplePaths"=wantSamplePaths,
    "deltaMin"=deltaMin, "wantSimData"=wantSimData,
    "pb"=pb, "seed"=seed, "nSim"=nSim, "nBoot"=nBoot,
    "sigma"=sigma, "sigma2"=sigma2,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "highN"=NULL, "alphaRelevance"=alphaRelevance,
    "nuMin"=nuMin, ...)


  result <- designSavi1aHelper("samplingResult"=samplingResult,
                               "esMin"=deltaMin, "power"=power,
                               "beta"=NULL, "ratio"=ratio, "testType"=testType)
  return(result)
}

#' Helper function to designing a savi T-test (output power)
#'
#' Finds the parameter and power when provided with only alpha, esMin, and nPlan
#'
#' @inheritParams designSaviT
#'
#' @return A list with the parameter and beta amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviT2WantPower(deltaTrue=0.9, nPlan=7, nSim=10)
designSaviT2WantPower <- function(
    deltaTrue, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    deltaMin=deltaTrue,
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    wantSamplePaths=TRUE, wantSimData=TRUE,
    relevanceTest=FALSE, relevanceSize=NULL,
    alphaRelevance=NULL,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    nuMin=2, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1

  nPlan <- checkAndReturnNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  samplingResult <- computePowerSaviT(
    "deltaTrue"=deltaTrue, "nPlan"=nPlan, "alpha"=alpha,
    "alternative"=alternative,
    "testType"=testType, "parameter"=parameter,
    "eType"=eType, "wantSamplePaths"=wantSamplePaths,
    "wantSimData"=wantSimData, "deltaMin"=deltaMin,
    "seed"=seed, "nSim"=nSim, "nBoot"=nBoot, "pb"=pb,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "alphaRelevance"=alphaRelevance, "nuMin"=nuMin, ...)

  result <- designSavi2Helper("samplingResult"=samplingResult,
                              "esMin"=deltaMin, "nPlan"=nPlan, "ratio"=ratio,
                              "testType"=testType)
  return(result)
}

#' Helper function to designing a Savi T-test (output esMin)
#'
#' Finds the parameter and esMin when provided with only alpha, power, and nPlan
#'
#' @inheritParams designSaviT
#'
#' @return A list with the parameter and the targeted esMin amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviT3WantEsMin(power=0.7, nPlan=10)
designSaviT3WantEsMin <- function(
    power, nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL, beta=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1
  nPlan <- checkAndReturnNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  result <- list("parameter"=NULL, "esMin"=NULL,
                 "nPlan"=nPlan, "beta"=beta, "ratio"=ratio,
                 "note"=NULL)

  deltaMin <- tryOrFailWithNA(
    computeMinEsBatchSaviT(
      "nPlan"=nPlan, "alpha"=alpha, "power"=power, "beta"=NULL,
      "alternative"=alternative, "testType"=testType,
      "parameter"=parameter, "eType"=eType,
      "lowEsTrue"=lowEsTrue, "highEsTrue"=highEsTrue)
  )

  if (is.null(parameter)) {
    parameter <- switch(eType,
                        "mom"=deltaMin^2/2,
                        "eGauss"=deltaMin^2,
                        "imom"=abs(deltaMin),
                        "eCauchy"=abs(deltaMin),
                        "grow"=deltaMin)
  }

  result[["parameter"]] <- parameter
  result[["esMin"]] <- deltaMin

  if (is.na(deltaMin))
    result[["note"]] <- "No deltaMin found for the provided beta and nPlan"
  else
    result[["note"]] <- "The reported deltaMin is based on the batch analysis."

  return(result)
}

#' Helper function to designing a savi T-test (output deltaMin based on the shortest interval at nPlan)
#'
#' Finds the parameter and deltaMin when provided with only alpha and nPlan
#'
#' @inheritParams designSaviT
#'
#' @return A list with the parameter and the parameter amongst other items
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' designSaviT3bWantParameter(nPlan=20)
designSaviT3bWantParameter <- function(
    nPlan,
    alpha=0.05, alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    ...) {
  # TODO(Alexander): Two-sample and imom don't play well


  defaultErrorText <-
    "Can't compute a design based on alpha and the provided planned sample size(s) alone."

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1
  nPlan <- checkAndReturnNPlan("nPlan"=nPlan, "ratio"=ratio, "testType"=testType)

  n1 <- nPlan[1]
  n2 <- nPlan[2]

  paired <- if (testType=="paired") TRUE else FALSE

  nEff <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1 else (1/n1+1/n2)^(-1)
  nu <- if (is.null(n2) || is.na(n2) || paired==TRUE) n1-1 else n1+n2-2

  if (nu <= 0)
    stop(defaultErrorText)

  result <- list("parameter"=NULL, "esMin"=NULL,
                 "nPlan"=nPlan, "ratio"=ratio,
                 "note"=NULL)

  minG <- (alpha^(-2/nu)-1)/nEff

  someTargetFunction <- function(g)
    tTestWidthDerivative(g, nEff=nEff, nu=nu, alpha=alpha)

  someTargetFunction <- Vectorize(someTargetFunction)

  tempResult <- stats::uniroot(
    someTargetFunction,
    c(minG, max(exp(-log(alpha))*minG, 1e6)),
    tol=min(.Machine$double.eps^0.25, 1/nEff))

  gCandidate <- tempResult[["root"]]
  deltaMinCandidate <- sqrt(gCandidate)

  if (eType=="eGauss") {
    parameter <- gCandidate
    deltaMin <- deltaMinCandidate
  } else {
    upperDelta <- if (eType %in% c("mom", "grow")) 2*deltaMinCandidate else max(2*deltaMinCandidate, 0.02)

    deltaDomain <- seq(deltaMinCandidate/4, upperDelta, length.out=1e3)
    ciWidths <- rep(Inf, length(deltaDomain))

    parameterDomain <- if (eType=="mom") deltaDomain^2/2 else deltaDomain

    for (i in seq_along(ciWidths)) {
      tempRes <- computeConfidenceIntervalT(meanObs=0, sdObs=1,
                                            nEff=nEff, nu=nu,
                                            parameter=parameterDomain[i],
                                            eType=eType, ciValue=1-alpha,
                                            alternative="twoSided")
      ciWidths[i] <- tempRes[2]
    }

    if (sum(is.infinite(ciWidths))==length(ciWidths))
      stop(defaultErrorText)

    minIndex <- which.min(ciWidths)

    deltaMin <- deltaDomain[minIndex]
    parameter <- parameterDomain[minIndex]

    if (minIndex!=1 && is.infinite(ciWidths[minIndex-1])) {
      result[["note"]] <- "Unstable design based on alpha and nPlan alone."
      warning("Unstable: The parameter corresponds to the smallest parameter value",
              "for which the ci width can be calculated. Another eType might yield more stable designs.")
    }
  }

  if (eType=="grow" && alternative=="less") {
    parameter <- -parameter
    deltaMin <- -deltaMin
  }

  result[["parameter"]] <- parameter
  result[["esMin"]] <- deltaMin

  return(result)
}

# Batch design fnts ------

#' Helper function: Computes the planned sample size for the savi T-test
#' based deltaMin, alpha and power
#'
#' @inheritParams designSaviT
#' @inheritParams sampleStoppingTimesSaviT
#'
#' @return a list which contains at least nPlan and the savi test defining parameter
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
computeNPlanBatchSaviT <- function(
    deltaTrue, alpha=0.05, power=0.8,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    sigma=1, sigma2=1,
    parameter=NULL, beta=NULL, ratio=1, deltaMin=NULL, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  result <- list(nPlan=NULL, "parameter"=parameter)

  n1Plan <- NULL
  n2Plan <- NULL

  n1OverNEffRatio <- if (testType=="twoSample") (1+ratio)/ratio else 1

  if (is.null(deltaMin))
    deltaMin <- deltaTrue

  deltaMin <- suppressWarnings(
    checkAndReturnEsMinParameterSide(
      "paramToCheck"=deltaMin, "alternative"=alternative,
      "esMinName"="deltaMin")
  )

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="z",
    "esMin"=deltaMin, "alternative"=alternative, "eType"=eType
  )

  deltaTrue <- abs(deltaTrue)

  # Sample size of greater sided Z-test as lower/upper bound for
  # the candidate set of nEff
  qB <- stats::qnorm(1-power)

  nTemp <- exp(-2*log(deltaTrue))*
    (2*qB^2 - 2*qB*sqrt(qB^2+2*log(1/alpha))
     +2*log(1/alpha))

  tempAlternative <- switch(alternative,
                            "twoSided"="twoSided",
                            "greater"="greater",
                            "less"="greater")

  if (testType=="twoSample") {
    n1Func <- function(nEff) (1+ratio)/ratio*nEff
    n2Func <- function(nEff) (1+ratio)*nEff
    nuFunc <- function(nEff) (1+ratio)^2/ratio*nEff-2
  } else if (testType %in% c("oneSample", "paired")) {
    n1Func <- function(nEff) nEff
    n2Func <- function(nEff) NULL
    nuFunc <- function(nEff) nEff-1
  }

  # TODO(Alexander):
  #
  #       THIS IS PROBABLY THE PLACE TO FIGURE OUT how Welch works
  #
  targetFunction <- function(nEff) {
    saviTTestStat(
      stats::qt("p"=1-power, "df"=nuFunc(nEff), "ncp"=sqrt(nEff)*deltaTrue),
      "n1"=n1Func(nEff), "n2"=n2Func(nEff), "parameter"=parameter, "alternative"=tempAlternative,
      "eType"=eType)$eValue-1/alpha
  }

  targetFunction <- Vectorize(targetFunction)

  tempResult <- suppressWarnings(
    tryCatch(stats::uniroot(targetFunction, interval=c(nTemp/2, 2*nTemp)),
             error=identity)
  )


  if (inherits(tempResult, "simpleError")) {
    tempResult <- suppressWarnings(
      tryCatch(stats::uniroot(targetFunction, interval=c(10*nTemp, 50*nTemp)),
               error=identity)
    )

    # if (eType=="bayarri") {
    #   tempResult <- suppressWarnings(
    #     tryCatch(stats::uniroot(targetFunction, interval=c(nTemp/2, 50*nTemp)),
    #              error=identity)
    #   )
    # }
  }

  if (inherits(tempResult, "simpleError"))
    stop("Can't compute the batched planned sample size")

  nEff <- tempResult[["root"]]

  # Two-sample/paired stuff
  if (testType == "twoSample") {
    n1Plan <- ceil(nEff * n1OverNEffRatio)
    n2Plan <- ceil(nEff * n1OverNEffRatio * ratio)
  } else {
    n1Plan <- ceil(nEff)
    n2Plan <- if (testType == "paired") n1Plan else NULL
  }

  if (is.null(n2Plan)) {
    result[["nPlan"]] <- n1Plan
    names(result[["nPlan"]]) <- "n1Plan"
  } else {
    result[["nPlan"]] <- c(n1Plan, n2Plan)
    names(result[["nPlan"]]) <- c("n1Plan", "n2Plan")
  }

  if (eType=="grow" && alternative=="less" && parameter > 0)
    parameter <- -parameter

  names(parameter) <- switch(eType,
                             "mom"="gMom",
                             "eGauss"="g",
                             "imom"="tau",
                             "eCauchy"="kappaG",
                             "grow"="deltaS",
                             "lai"="")

  result[["parameter"]] <- parameter

  return(result)
}

#' Computes the smallest detectable deltaMin with power probability, for the provided
#' sample size
#'
#' @inheritParams  designSaviT
#'
#' @return numeric > 0 that represents the minimal detectable effect size
#' @export
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @examples
#' computeMinEsBatchSaviT(27)
computeMinEsBatchSaviT <- function(
    nPlan, alpha=0.05, power=0.8,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL, beta=NULL,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    lowEsTrue=0.01, highEsTrue=3, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)
  testType <- match.arg(testType)

  nEff <- computeNEff("n"=nPlan, "testType" = testType)

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (eType=="mom") {
    paramFunc <- function(deltaTrue) deltaTrue^2/2
  } else if (eType=="eGauss") {
    paramFunc <- function(deltaTrue) deltaTrue^2
  } else if (eType=="imom") {
    paramFunc <- function(deltaTrue) abs(deltaTrue)
  } else if (eType=="eCauchy") {
    paramFunc <- function(deltaTrue) abs(deltaTrue)
  } else if (eType=="grow") {
    paramFunc <- function(deltaTrue) deltaTrue
  }

  ratio <- if (length(nPlan)==2) nPlan[2]/nPlan[1] else 1

  tempAlternative <- switch(alternative,
                            "twoSided"="twoSided",
                            "greater"="greater",
                            "less"="greater")

  if (testType %in% c("oneSample", "paired")) {
    n1 <- nPlan[1]
    n2 <- NULL
    nu <- n1-1
  } else if (testType=="twoSample") {
    n1 <- nPlan[1]
    n2 <- nPlan[2]
    nu <- n1+n2-2
  }

  targetFunction <- function(deltaTrue) {
    saviTTestStat(
      stats::qt("p"=1-power, "df"=nu, "ncp"=sqrt(nEff)*deltaTrue),
      "n1"=n1, "n2"=n2, "parameter"=paramFunc(deltaTrue),
      "alternative"=tempAlternative, "eType"=eType)$eValue-1/alpha
  }

  targetFunction <- Vectorize(targetFunction)

  if (eType=="grow")  {
    gaussResult <- computeMinEsBatchSaviT(
      "nPlan"=nPlan, "alpha"=alpha, "power"=power,
      "beta"=NULL, "alternative"=tempAlternative,
      testType=testType, eType="eGauss")
  }

  tempResult <- try(stats::uniroot(targetFunction, interval=c(lowEsTrue, highEsTrue)))

  result <- tempResult[["root"]]

  if (alternative=="less")
    result <- -result

  return(result)
}

# Sampling functions for design ----

#' Simulate stopping times for the savi T-test
#'
#' @inheritParams designSaviT
#' @inheritParams generateNormalData
#' @param nMax integer > 0, maximum sample size of the (first) sample in each sample path.
#' @param wantEValuesAtNMax logical. If \code{TRUE} then compute eValues at nMax. Default \code{FALSE}.
#' @param wantSamplePaths logical. If \code{TRUE} then output the (stopped) sample paths. Default \code{TRUE}.
#' @param wantSimData logical. If \code{TRUE}, then output the simulated data.
#' @param lowN integer, smallest sample size (of the first group).
#'
#' @return a list with stoppingTimes and breakVector. Entries of breakVector are 0, 1. A 1 represents stopping
#' due to exceeding nMax, and 0 due to 1/alpha threshold crossing, which implies that in corresponding stopping
#' time is Inf.
#'
#' @export
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @examples
#' sampleStoppingTimesSaviT(0.7, nSim=10, nMax=20)
sampleStoppingTimesSaviT <- function(
    deltaTrue, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    ratio=1, deltaMin=NULL, parameter=NULL,
    lowN=3L, nMax=1e8L,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    wantEValuesAtNMax=FALSE, nuMin=2, power=NULL,
    wantSamplePaths=TRUE, wantSimData=TRUE,
    sigma=1, sigma2=1,
    pb=TRUE, seed=NULL, nSim=1e3L, relevanceTest=FALSE,
    relevanceSize=NULL, beta=NULL, alphaRelevance=NULL, ...) {

  stopifnot(alpha > 0, alpha <= 1,
            is.finite(nMax),
            is.finite(deltaTrue))

  power <- matchPowerWith("power"=power, "beta"=beta)

  if (relevanceTest) {
    relevanceSize <- matchRelevanceParameterWith(
      "relevanceSize"=relevanceSize, "esMin"=deltaMin,
      "esTrue"=deltaTrue
    )

    alphaRelevance <- matchAlphaRelevanceWith(
      "alphaRelevance"=alphaRelevance, "alpha"=alpha, "power"=power, "beta"=1-power
    )

    stopifnot(alphaRelevance > 0, alphaRelevance < 1)
  }

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  result <- constructSampleStoppingTimesList(
    "nSim"=nSim, "nMax"=nMax,
    "wantEValuesAtNMax"=wantEValuesAtNMax,
    "wantSamplePaths"=wantSamplePaths)

  deltaTrue <- abs(deltaTrue)

  deltaMin <- if (is.null(deltaMin)) abs(deltaTrue) else abs(deltaMin)

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="t",
    "esMin"=deltaMin, "alternative"=alternative,
    "eType"=eType
  )

  relevanceTestSim <- NULL

  if (relevanceTest) {
    relevanceParameter <- matchRelevanceParameterWith(
      "relevanceSize"=relevanceSize,
      "esMin"=deltaMin, "esTrue"=deltaTrue)

    relevanceTestSim <-
      list("eValuesStopped"=Matrix::sparseVector(x=0, i=1, length=nSim),
           "samplePaths"=result[["samplePaths"]],
           "stoppingTimes"=Matrix::sparseVector(x=0, i=1, length=nSim),
           "parameter"=relevanceParameter, "alpha"=alphaRelevance)
  }

  if (testType=="twoSample" && length(nMax)==1) {
    nMax <- c(nMax, ceil(ratio*nMax))
    n1Max <- nMax[1]
    n2Max <- nMax[2]
    ratio <- nMax[2]/nMax[1]
  } else if (testType %in% c("paired", "oneSample")){
    n1Max <- nMax[1]
    n2Max <- NULL
    nMax <- n1Max
    ratio <- 1
  }

  if (pb)
    pbSavi <- utils::txtProgressBar(style=3, title="Savi test threshold crossing")

  tempN <- defineTTestN("lowN"=1, "highN"=nMax[1], "ratio"=ratio, "testType"=testType)

  n1Vector <- tempN[["n1"]]
  n2Vector <- tempN[["n2"]]
  nEffVector <- tempN[["nEff"]]
  nuVector <- tempN[["nu"]]

  # TODO(Alexander): Here check WElch
  #
  #
  simData <- generateNormalData("nPlan"=nMax, "nSim"=nSim,
                                "deltaTrue"=deltaTrue,
                                "sigma"=sigma, "sigma2"=sigma,
                                "paired"=FALSE, "seed"=seed)

  for (sim in seq_along(result[["stoppingTimes"]])) {
    if (testType %in% c("oneSample", "paired")) {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/n1Vector*cumsum(x1)
      x1SquareVector <- cumsum(x1^2)
      sX1Vector <- sqrt(1/(n1Vector-1)*(x1SquareVector - n1Vector*x1BarVector^2))

      badIndeces <- which(n1Vector-1 <= 0)
      sX1Vector[badIndeces] <- 1

      tValues <- sqrt(nEffVector)*x1BarVector/sX1Vector
    } else {
      x1 <- simData[["dataGroup1"]][sim, ]
      x1BarVector <- 1/(n1Vector)*cumsum(x1)
      x1BarVector <- x1BarVector[n1Vector]
      x1SquareVector <- cumsum(x1^2)[n1Vector]

      x2 <- simData[["dataGroup2"]][sim, ]
      x2CumSum <- cumsum(x2)[n2Vector]
      x2BarVector <- 1/(n2Vector)*x2CumSum
      x2SquareVector <- cumsum(x2^2)[n2Vector]

      sPVector <- sqrt(1/(n1Vector+n2Vector-2)*
                         (x1SquareVector-n1Vector*x1BarVector^2 + x2SquareVector - n2Vector*x2BarVector^2))

      badIndeces <- which(n1Vector+n2Vector-2 <= 0)
      sPVector[badIndeces] <- 1

      tValues <- sqrt(nEffVector)*(x1BarVector-x2BarVector)/sPVector
    }

    if (wantEValuesAtNMax) {
      tempResult <- saviTTestStat("t"=tValues[length(tValues)],
                                  "parameter"=parameter,
                                  "n1"=nMax[1], n2=nMax[2],
                                  "alternative"=alternative, "eType"=eType,
                                  "nuMin"=nuMin)
      result[["eValuesAtNMax"]][sim] <- tempResult[["eValue"]]
    }

    for (j in seq_along(n1Vector)) {
      tempResult <- suppressWarnings(
        saviTTestStat("t"=tValues[j], "parameter"=parameter,
                      "n1"=n1Vector[j], "n2"=n2Vector[j],
                      "alternative"=alternative,
                      "eType"=eType, "nuMin"=nuMin)
      )

      evidenceNow <- tempResult[["eValue"]]

      if (wantSamplePaths)
        result[["samplePaths"]][sim, j] <- evidenceNow

      if (evidenceNow >= 1/alpha) {
        result[["stoppingTimes"]][sim] <- n1Vector[j]
        result[["eValuesStopped"]][sim] <- evidenceNow

        if (wantSamplePaths) {
          result[["samplePaths"]][sim, j:nMax[1]] <- evidenceNow
        }
        break()
      }

      if (relevanceTest) {
        relevanceRes <- saviRelevanceTStatNEffNu(
          "t"=tValues[j], "nEff"=nEffVector[j], "nu"=nuVector[j],
          "parameter"=relevanceParameter,
          "alternative"=alternative,
          "tDensity"=FALSE,
          "paired"=ifelse(testType=="paired", TRUE, FALSE),
          "nuMin"=nuMin)

        relevanceEValue <- relevanceRes[["eValue"]]

        if (wantSamplePaths)
          relevanceTestSim[["samplePaths"]][sim, j] <- relevanceEValue

        if (relevanceEValue < alphaRelevance) {
          result[["breakVector"]][sim] <- -1
          result[["stoppingTimes"]][sim] <- nMax

          relevanceTestSim[["stoppingTimes"]][sim] <- n1Vector[j]
          relevanceTestSim[["eValuesStopped"]][sim] <- relevanceEValue

          break()
        }
      }

      # Note(Alexander): If passed maximum nPlan[1] stop.
      #   For power calculations if beyond nPlan[1], then set to Inf, doesn't matter for the quantile
      #
      if (n1Vector[j] >= nMax[1]) {
        result[["stoppingTimes"]][sim] <- n1Vector[j]
        result[["breakVector"]][sim] <- 1
        result[["eValuesStopped"]][sim] <- evidenceNow
        break()
      }
    }

    if (pb)
      utils::setTxtProgressBar(pbSavi, "value"=sim/nSim, "title"="Trials")
  }

  if (pb)
    close(pbSavi)


  result[["parameter"]] <- parameter
  result[["n1Vector"]] <- n1Vector
  result[["ratio"]] <- ratio

  if (isTRUE(wantSimData))
    result[["simData"]] <- simData


  if (wantSamplePaths) {
    samplePaths <- result[["samplePaths"]]
    samplePaths[is.na(samplePaths)] <- 0
    result[["samplePaths"]] <- Matrix::Matrix(samplePaths, sparse=TRUE)

    if (relevanceTest) {
      relevanceSamplePaths <- relevanceTestSim[["samplePaths"]]
      relevanceSamplePaths[is.na(relevanceSamplePaths)] <- 0
      relevanceTestSim[["samplePaths"]] <- Matrix::Matrix(relevanceSamplePaths, sparse=TRUE)
    }
  }

  result[["relevanceTestSim"]] <- relevanceTestSim

  return(result)
}


#' Helper function: Computes the power of the saviTTest based on deltaMin and nPlan
#'
#' @inheritParams designSaviT
#' @inheritParams saviTTestStat
#' @inheritParams sampleStoppingTimesSaviT
#'
#' @return a list which contains at least power and an adapted bootObject of class
#' \code{\link[boot]{boot}()}.
#' @export
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @examples
#' computePowerSaviT(deltaTrue=0.7, 27, nSim=10)
computePowerSaviT <- function(
    deltaTrue, nPlan, alpha=0.05,
    alternative=c("twoSided", "greater", "less"),
    testType=c("oneSample", "paired", "twoSample"),
    parameter=NULL, deltaMin=deltaTrue,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    wantSamplePaths=TRUE, nuMin=2, wantSimData=TRUE,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    relevanceTest=FALSE, relevanceSize=NULL,
    alphaRelevance=NULL, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  ratio <- if (length(nPlan) == 2) nPlan[2]/nPlan[1] else 1

  if (testType=="twoSample" && length(nPlan)==1) {
    nPlan <- c(nPlan, nPlan)
    warning('testType=="twoSample" specified, but nPlan[2] not provided. nPlan[2] is set to ratio = ', ratio,
            'times nPlan[1] = ', nPlan[2])
  }

  deltaTrue <- checkAndReturnEsMinParameterSide(
    "paramToCheck"=deltaTrue, "alternative"=alternative,
    "esMinName"="deltaTrue")

  parameter <- matchEParameterWith(
    "parameter"=parameter, "analysisType"="t",
    "esMin"=deltaMin,
    "alternative"=alternative, "eType"=eType
  )

  samplingResult <- sampleStoppingTimesSaviT(
    "deltaTrue"=deltaTrue, "alpha"=alpha, "power"=NULL,
    "alternative" = alternative, "testType"=testType,
    "ratio"=ratio, "parameter"=parameter, "nMax"=nPlan,
    "eType"=eType, "nuMin"=nuMin, "wantSimData"=wantSimData,
    "wantEValuesAtNMax"=TRUE, "wantSamplePaths"=wantSamplePaths,
    "pb"=pb, "seed"=seed, "nSim"=nSim, "beta"=NULL,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "alphaRelevance"=alphaRelevance, ...)

  result <- computePowerBootstrapper(
    "samplingResult"=samplingResult, "parameter"=parameter,
    "nPlan"=nPlan, "nBoot"=nBoot)

  return(result)
}


#' Helper function: Computes nPlan based on deltaMin and power
#'
#'
#' @inheritParams designSaviT
#' @inheritParams sampleStoppingTimesSaviT
#'
#' @return a list which contains at least nPlan and an adapted bootObject of class  \code{\link[boot]{boot}()}.
#'
#' @references
#'   `r addCite(grunwald2024safe)`
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' computeNPlanSaviT(0.7, 0.2, nSim=10)
computeNPlanSaviT <- function(
    deltaTrue, power=0.8, alpha=0.05,
    alternative = c("twoSided", "less", "greater"),
    testType=c("oneSample", "paired", "twoSample"),
    deltaMin=NULL, beta=NULL,
    ratio=1, parameter=NULL, nMax=1e8,
    eType=c("mom", "eGauss", "imom", "eCauchy", "grow", "lai"),
    wantSamplePaths=TRUE, wantSimData=TRUE, nuMin=2,
    pb=TRUE, seed=NULL, nSim=1e3L, nBoot=nSim,
    sigma=1, sigma2=1,
    relevanceTest=FALSE, relevanceSize=NULL,
    alphaRelevance=NULL, ...) {

  # TODO(Alexander): Remove in v0.9.0
  #
  if (length(alternative)==1 && alternative=="two.sided") {
    warning('The option alternative="two.sided" is deprecated;',
            'Please use alternative="twoSided" instead')
    alternative <- "twoSided"
  }

  alternative <- match.arg(alternative)
  testType <- match.arg(testType)
  eType <- match.arg(eType)

  if (is.null(deltaMin)) {
    deltaTrue <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=deltaTrue, "alternative"=alternative,
      "esMinName"="deltaTrue")

    deltaMin <- deltaTrue
  } else {
    deltaMin <- checkAndReturnEsMinParameterSide(
      "paramToCheck"=deltaMin, "alternative"=alternative,
      "esMinName"="deltaMin")
  }

  power <- matchPowerWith("power"=power, "beta"=beta)

  tempObj <- computeNPlanBatchSaviT(
    "deltaTrue"=deltaMin, "alpha"=alpha, "power"=power,
    "alternative"=alternative, "testType"=testType,
    "parameter"=parameter, "ratio"=ratio, "eType"=eType,
    "deltaMin"=deltaMin, #"highN"=highN,
    "sigma"=sigma, "sigma2"=sigma2,
    "wantSimData"=wantSimData, ...)

  nPlanBatch <- tempObj[["nPlan"]]

  if (is.null(parameter))
    parameter <- tempObj[["parameter"]]

  samplingResult <- sampleStoppingTimesSaviT(
    "deltaTrue"=deltaTrue, "alpha"=alpha, "power"=power,
    "alternative"=alternative, "testType"=testType,
    "ratio"=ratio, "parameter"=parameter, "nMax"=nPlanBatch,
    "eType"=eType, "deltaMin"=deltaMin, "nuMin"=nuMin,
    "wantSamplePaths"=wantSamplePaths, "wantSimData"=wantSimData,
    "pb"=pb, "seed"=seed, "nSim"=nSim, "beta"=NULL,
    "sigma"=sigma, "sigma2"=sigma2,
    "relevanceTest"=relevanceTest, "relevanceSize"=relevanceSize,
    "alphaRelevance"=alphaRelevance, ...)

  result <- computeNPlanBootstrapper(
    "samplingResult"=samplingResult, "parameter"=parameter,
    "beta"=NULL, "power"=power, "nPlanBatch"=nPlanBatch, "nBoot"=nBoot
  )
  return(result)
}


# Helper fnts ------

#' Computes a Sequence of (Effective) Sample Sizes of Z- and T-Tests
#'
#' Helper function that outputs a sequence of sample sizes, effective sample sizes,
#' and the degrees of freedom depending on the type of T-test. Also used for Z-tests.
#'
#' @inheritParams designSaviT
#'
#' @param lowN integer that defines the smallest n of our search space for n.
#' @param highN integer largest sample size of the (first) sample. Default set to 100.
#'
#' @return Returns the sample sizes and degrees of freedom.
defineTTestN <- function(lowN=3, highN=100, ratio=1,
                         testType=c("oneSample", "paired", "twoSample")) {
  testType <- match.arg(testType)

  if (testType %in% c("twoSample")) {
    n1 <- lowN:highN
    n2 <- ceil(ratio*n1)
    nEff <- (1/n1+1/n2)^(-1)
    nu <- n1+n2-2
  } else if (testType %in% c("oneSample", "paired")) {
    n1 <- lowN:highN
    n2 <- NULL
    nEff <- n1
    nu <- nEff-1
  }
  result <- list("n1"=n1, "n2"=n2, "nEff"=nEff, "nu"=nu)
  return(result)
}

#' Helper function to compute the relevant summary statistics in z and t-tests
#'
#' @inheritParams saviTTest
#' @inheritParams designSaviT
#'
#' @returns a list of summary statistics
#' @export
#'
#' @examples
#' n1 <- 23
#' n2 <- 39
#'
#' x <- rnorm(n1, sd=7)
#' y <- rnorm(n2, sd=3)
#'
#' computeZTSumStats(x=x, y=y)
#'
computeZTSumStats <- function(
    x, y=NULL, sequential=NULL,
    paired=FALSE, varEqual=TRUE,
    testType=NULL) {

  # check data
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (paired) {
    xGoodIndeces <- yGoodIndeces  <-
      stats::complete.cases(x, y)
    x <- x[xGoodIndeces]
    y <- y[yGoodIndeces]
  }

  n1 <- length(x)
  n2 <- length(y)

  if (is.null(testType)) {
    if (is.null(y)) {
      testType <- "oneSample"
    } else if (paired) {
      testType <- "paired"
    } else {
      testType <- "twoSample"
    }
  }

  if (is.null(sequential))
    sequential <- if (n1 <= 200) TRUE else FALSE

  if (sequential) {
    tempN <- defineTTestN("lowN"=1, "highN"=n1, "ratio"=n2/n1, "testType"=testType)
  } else {
    tempN <- list("meanObsVec"=NULL, "sdObsVec"=NULL, "nEffVec"=NULL,
                  "n1Vec"=NULL, "n2Vec"=NULL, "nuVec"=NULL)
  }

  meanObsVec <- NULL
  sdObsVec <- NULL

  n1Vec <- tempN[["n1"]]
  n2Vec <- tempN[["n2"]]
  nuVec <- tempN[["nu"]]
  nEffVec <- tempN[["nEff"]]

  meanObsVec <- NULL
  sdObsVec <- NULL

  if (testType=="oneSample") {
    n <- nEff <- n1 <- length(x)
    n2 <- NULL
    nu <- n-1

    meanObs <- estimate <- mean(x)
    sdObs <- stats::sd(x)

    names(estimate) <- "mean of x"
    names(n) <- "n1"

    if (sequential) {
      meanObsVec <- 1/nEffVec*cumsum(x)
      sdObsVec <- sqrt(1/nuVec*(cumsum(x^2)-nEffVec*meanObsVec^2))
    }
  } else if (testType=="paired") {
    if (n1 != n2)
      stop("Data error: Error in complete.cases(x, y): Paired analysis requested, ",
           "but the two samples are not of the same size.")

    nEff <- n1
    nu <- n1-1

    meanObs <- estimate <- mean(x-y)
    sdObs <- stats::sd(x-y)
    names(estimate) <- "mean of the differences"

    if (sequential) {
      meanObsVec <- 1/nEffVec*cumsum(x-y)
      sdObsVec <- sqrt(1/nuVec*(cumsum((x-y)^2)-nEffVec*meanObsVec^2))
    }
  } else if (testType=="twoSample") {
    nEff <- (1/n1+1/n2)^(-1)

    varX <- stats::var(x)
    varY <- stats::var(y)

    if (is.na(varX))
      varX <- 0

    if (is.na(varY))
      varY <- 0

    if (varEqual) {
      nu <- n1+n2-2
    } else {
      nu <- (varX/n1+varY/n2)^2/
        ((varX/n1)^2/(n1-1)+(varY/n2)^2/(n2-1))

      if (is.na(nu))
        nu <- 0

    }

    if (varEqual) {
      sPooledSquared <- ((n1-1)*varX+(n2-1)*varY)/nu
      sdObs <- sqrt(sPooledSquared)
    } else {
      sdObs <- sqrt((1/n1*varX+1/n2*varY)*nEff)
    }

    estimate <- c(mean(x), mean(y))
    names(estimate) <- c("mean of x", "mean of y")
    meanObs <- estimate[1]-estimate[2]

    if (sequential) {
      xMeanObsRaw <- 1/(1:n1)*cumsum(x)
      yMeanObsRaw <- 1/(1:n2)*cumsum(y)

      xSumsOfSquaresRaw <- (cumsum(x^2)-(1:n1)*xMeanObsRaw^2)
      ySumsOfSquaresRaw <- (cumsum(y^2)-(1:n2)*yMeanObsRaw^2)

      if (n2/n1==1) {
        xMeanObsVec <- xMeanObsRaw
        yMeanObsVec <- yMeanObsRaw
        xSumsOfSquaresVec <- xSumsOfSquaresRaw
        ySumsOfSquaresVec <- ySumsOfSquaresRaw
      } else {
        vecLength <- length(n1Vec)

        xMeanObsVec <- yMeanObsVec <-
          xSumsOfSquaresVec <- ySumsOfSquaresVec <- numeric(vecLength)

        for (j in 1:vecLength) {
          nowN1 <- n1Vec[j]
          nowN2 <- n2Vec[j]

          xMeanObsVec[j] <- xMeanObsRaw[nowN1]
          yMeanObsVec[j] <- yMeanObsRaw[nowN2]
          xSumsOfSquaresVec[j] <- xSumsOfSquaresRaw[nowN1]
          ySumsOfSquaresVec[j] <- ySumsOfSquaresRaw[nowN2]
        }
      }

      meanObsVec <- xMeanObsVec-yMeanObsVec

      if (varEqual==TRUE) {
        sPooledSquaredVec <- (xSumsOfSquaresVec+ySumsOfSquaresVec)/nuVec
        sdObsVec <- sqrt(sPooledSquaredVec)
      } else {
        varXVec <- xSumsOfSquaresVec/(n1Vec-1)
        varYVec <- ySumsOfSquaresVec/(n2Vec-1)

        logNumeratorVec <- 2*log(varXVec/n1Vec+varYVec/n2Vec)

        logDenominatorVec <- log(
          exp(2*log(varXVec)-2*log(n1Vec)-log(n1Vec-1)) +
            exp(2*log(varYVec)-2*log(n2Vec)-log(n2Vec-1))
        )

        nuVec <- exp(logNumeratorVec-logDenominatorVec)

        # Alexander: Remove bad data
        # Zero sums-of-squares means Inf t, but nu = 0
        badIndeces <- which(is.na(nuVec))

        nuVec[badIndeces] <- 0
        varXVec[badIndeces] <- 0
        varYVec[badIndeces] <- 0

        sSquaredVec <- 1/n1Vec*varXVec + 1/n2Vec*varYVec
        sdObsVec <- sqrt(sSquaredVec*nEff)
      }
    }
  }

  n <- if (testType=="oneSample") n1 else c(n1, n2)
  names(n) <- if (testType=="oneSample") "n1" else c("n1", "n2")

  res <- list("n"=n, "nEff"=nEff, "n1"=n1, "n2"=n2, "nEff"=nEff,
              "nu"=nu, "meanObs"=meanObs,"sdObs"=sdObs,
              "meanObsVec"=meanObsVec, "sdObsVec"=sdObsVec,
              "estimate"=estimate, "nEffVec"=nEffVec, "n1Vec"=n1Vec,
              "n2Vec"=n2Vec, "nuVec"=nuVec, "sequential"=sequential)

  return(res)
}

#' Computes 6 equivalent forms of the confluent hypergeometric function 1f1
#'
#' Repeated calls of \code{\link[hypergeo]{genhypergeo}},
#' which provides further details.
#'
#' @param U Upper arguments respectively (real or complex)
#' @param L Lower arguments respectively (real or complex)
#' @param z Primary complex argument
#' @param tol tolerance with default zero meaning to iterate
#' until additional terms to not change the partial sum
#' @param maxiter Maximum number of iterations to perform
#'
#' @returns a vector
#' @export
#'
#' @examples
#' compute1F1AllVersions(U=-359, L=1/2, z=-0.1891234)
compute1F1AllVersions <- function(U, L, z, tol=0, maxiter=2000) {
  a1 <- Re(hypergeo::genhypergeo_contfrac_single(
    "U"=U, "L"=L, "z"=z, "tol"=tol, "maxiter"=maxiter))
  a2 <- Re(hypergeo::genhypergeo_series(
    "U"=U, "L"=L, "z"=z, "tol"=tol, "maxiter"=maxiter))
  a3 <- Re(hypergeo::genhypergeo_shanks(
    "U"=U, "L"=L, "z"=z, "maxiter"=maxiter))

  a4 <- Re(hypergeo::genhypergeo_contfrac_single(
    "U"=L-U, "L"=L, "z"=-z, "tol"=tol, "maxiter"=maxiter))*exp(z)
  a5 <- Re(hypergeo::genhypergeo_series(
    "U"=L-U, "L"=L, "z"=-z, "tol"=tol, "maxiter"=maxiter))*exp(z)
  a6 <- Re(hypergeo::genhypergeo_shanks(
    "U"=L-U, "L"=L, "z"=-z, "maxiter"=maxiter))*exp(z)

  # res <- list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5, a6=a6)
  res <- c(a1, a2, a3, a4, a5, a6)
  return(res)
}

# Data generating fnt ------

#' Generates Normally Distributed Data Depending on the Design
#'
#' The designs supported are "oneSample", "paired", "twoSample".
#'
#' @inheritParams designSaviT
#' @inheritParams saviTTest
#'
#' @param meanDiffTrue numeric representing the true mean for
#' simulations with a Z-test. Default \code{NULL}
#' @param muGlobal numeric, population grand mean
#' @param sigma numeric > 0, population standard deviation
#' @param meanDiffTrue numeric, data governing parameter value
#' @param deltaTrue numeric, the value of the true standardised effect size (test-relevant parameter).
#' This argument is used by `designSaviT()` with `deltaTrue <- deltaMin`
#'
#' @return Returns a list of two data matrices contains at least the following components:
#'
#' \describe{
#'   \item{dataGroup1}{a matrix of data dimension nSim by \code{nPlan[1]}.}
#'   \item{dataGroup2}{a matrix of data dimension nSim by \code{nPlan[2]}.}
#' }
#' @export
#'
#' @examples
#' generateNormalData(20, 15, deltaTrue=0.3)
generateNormalData <- function(nPlan, nSim=1000L,
                               deltaTrue=NULL, muGlobal=0, sigma=1,
                               sigma2=1, paired=FALSE,
                               seed=NULL, meanDiffTrue=NULL) {
  stopifnot(all(nPlan > 0))

  if ((is.null(deltaTrue) && is.null(meanDiffTrue)) || !is.null(deltaTrue) && !is.null(meanDiffTrue))
    stop("Please provide either deltaTrue (T-test), or meanDiffTrue (Z-test).")

  result <- list("dataGroup1"=NULL, "dataGroup2"=NULL)
  set.seed(seed)

  # TODO(Alexander): vector("mode"="list", length=length(nPlan))

  n1Plan <- nPlan[1]

  # TODO(Alexander): Figure out here
  #
  if (is.null(meanDiffTrue))
    meanDiffTrue <- deltaTrue*sigma

  if (length(nPlan)==1) {
    dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=meanDiffTrue, "sd"=sigma)
    dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
    dataGroup2 <- NULL
  } else {
    n2Plan <- nPlan[2]

    if (paired) {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muGlobal + meanDiffTrue/sqrt(2), "sd"=sigma)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nSim, "mean"=muGlobal - meanDiffTrue/sqrt(2), "sd"=sigma)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nSim)
    } else {
      dataGroup1 <- stats::rnorm("n"=n1Plan*nSim, "mean"=muGlobal + meanDiffTrue/2, "sd"=sigma)
      dataGroup1 <- matrix(dataGroup1, "ncol"=n1Plan, "nrow"=nSim)
      dataGroup2 <- stats::rnorm("n"=n2Plan*nSim, "mean"=muGlobal - meanDiffTrue/2, "sd"=sigma)
      dataGroup2 <- matrix(dataGroup2, "ncol"=n2Plan, "nrow"=nSim)
    }
  }

  return(list("dataGroup1"=dataGroup1, "dataGroup2"=dataGroup2))
}

# Workshop functions ---------

#' A "subjective" Bayes factor for the two-sample T-test
#'
#' Based on conjugate priors with a total of 8 hyperparameters.
#'
#' @param x1 numeric, sample mean of group 1
#' @param sdObs1 numeric, the observed standard deviation of the first group.
#' @param n1 integer sample size of group 1
#' @param x2 numeric, sample mean of group 2
#' @param sdObs2 numeric, the observed standard deviation of the second group.
#' @param n2 integer sample size of group 2
#' @param a1 numeric, prior mean of the population mean mu1 of group 1
#' @param g1 numeric > 0, conditional prior variance of the population mean
#' mu1 of group 1 is given by \code{g1*sigma^2}
#' @param a2 numeric, prior mean of the population mean mu2 of group 2
#' @param g2 numeric > 0, conditional prior variance of the population mean
#' mu2 of group 2 is given by \code{g2*sigma^2}
#' @param a0 numeric, prior mean of the overall population mean mu0 of both groups
#' @param g0 numeric > 0, conditional prior variance of the population mean
#' mu0 of both groups is given by \code{g1*sigma^2}
#' @param aGamma numeric > 0, shape parameter of the prior on the
#' standard deviation sigma
#' @param bGamma numeric > 0, rate parameter of the prior on the
#' standard deviation sigma
#' @param log logical, default FALSE, if TRUE then return logarithm of the subjective Bayes factor outcome
#'
#' @return numeric > 0 representing the subjective Bayes factor outcome in favour of the alternative over the null
#'
#' @references
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' conjugateBfTStat(5.2, 2, 3, 3.4, 2, 12)
#'
conjugateBfTStat <- function(
    x1, sdObs1, n1, x2, sdObs2, n2,
    # a1=3.98, g1=0.05, a2=4.02, g2=0.02,
    # a0=4, g0=2,
    # a1=3.98, g1=0.5, a2=4.02, g2=0.1,
    # a0=4, g0=2,
    a1=3.98, g1=0.03,
    a2=4.02, g2=0.05,
    a0=4, g0=2,
    # a1=3.7, g1=0.1, a2=4.3, g2=0.3,
    # a0=4, g0=1e3,
    aGamma=2, bGamma=1/2, log=FALSE) {

  nu1 <- n1-1
  nu2 <- n2-1

  if (n1 <= 1) {
    sdObs1 <- 0
    nu1 <- 0
  }

  if (n2 <= 1) {
    sdObs2 <- 0
    nu2 <- 0
  }

  ssTerm <- nu1*sdObs1^2+nu2*sdObs2^2+2*bGamma
  nP <- n1+n2

  logBf10 <- 1/2*(log(1+g0*nP)-log(1+g1*n1)-log(1+g2*n2))+
    (nP+2*aGamma)/2*
    (log(n1*n2/nP*(x1-x2)^2+nP/(1+g0*nP)*(n1/nP*x1+n2/nP*x2-a0)^2+ssTerm) -
       log(n1/(1+g1*n1)*(x1-a1)^2+n2/(1+g2*n2)*(x2-a2)^2+ssTerm))

  if (isTRUE(log))
    return(logBf10)
  else
    return(exp(logBf10))
}

#' A "subjective" Bayes factor for the two-sample T-test
#'
#' Based on conjugate priors with a total of 8 hyperparameters.
#'
#' @param x1 numeric, sample mean of group 1
#' @param sdObs1 numeric, the observed standard deviation of the first group.
#' @param n1 integer sample size of group 1
#' @param x2 numeric, sample mean of group 2
#' @param sdObs2 numeric, the observed standard deviation of the second group.
#' @param n2 integer sample size of group 2
#' @param a1 numeric, prior mean of the population mean mu1 of group 1
#' @param g1 numeric > 0, conditional prior variance of the population mean
#' mu1 of group 1 is given by \code{g1*sigma^2}
#' @param a2 numeric, prior mean of the population mean mu2 of group 2
#' @param g2 numeric > 0, conditional prior variance of the population mean
#' mu2 of group 2 is given by \code{g2*sigma^2}
#' @param a0 numeric, prior mean of the overall population mean mu0 of both groups
#' @param g0 numeric > 0, conditional prior variance of the population mean
#' mu0 of both groups is given by \code{g1*sigma^2}
#' @param aGamma numeric > 0, shape parameter of the prior on the
#' standard deviation sigma
#' @param bGamma numeric > 0, rate parameter of the prior on the
#' standard deviation sigma
#' @param log logical, default FALSE, if TRUE then return logarithm of the subjective Bayes factor outcome
#'
#' @return numeric > 0 representing the subjective Bayes factor outcome in favour of the alternative over the null
#'
#' @references
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' conjugateBfTStat(5.2, 2, 3, 3.4, 2, 12)
#'
conjugateBfTStatOld <- function(
    x1, sdObs1, n1, x2, sdObs2, n2,
    # a1=3.98, g1=0.05, a2=4.02, g2=0.02,
    # a0=4, g0=2,
    a1=3.98, g1=0.5, a2=4.02, g2=0.1,
    a0=4, g0=2,
    # a1=3.7, g1=0.1, a2=4.3, g2=0.3,
    # a0=4, g0=1e3,
    aGamma=2, bGamma=1/2, log=FALSE) {

  nPlus <- n1+n2
  q1 <- n1/nPlus
  q2 <- n2/nPlus
  # nuCombined <- nCombined-1
  xGlobal <- q1*x1+q2*x2
  nu1 <- n1-1
  nu2 <- n2-1

  if (n1 <= 1 && n2 <= 1) {
    sdObs1 <- 0
    sdObs2 <- 0
    nu1 <- 0
    nu2 <- 0
  }

  ssTerm <- nu1*sdObs1^2+nu2*sdObs2^2+2*bGamma

  logBf10 <- 1/2*(log(1+nPlus*g0)-log(1+n1*g1)-log(1+n2*g2)) +
    (nPlus/2+bGamma)*(
      log(ssTerm+nPlus/(1+nPlus*g0)*(q1*x1+q2*x2-a0)^2 +nPlus*q1*q2*(x1-x2)) -
        log(ssTerm+n1/(1+n1*g1)*(x1-a1)^2+n1/(1+n2*g2)*(x2-a2)^2)
    )

  if (isTRUE(log))
    return(logBf10)
  else
    return(exp(logBf10))
}

#   {
#
#   nCombined <- n1+n2
#   nuCombined <- nCombined-1
#   xCombined <- (n1*x1+n2*x2)/nCombined
#   nu1 <- n1-1
#   nu2 <- n2-1
#
#   if (n1 <= 1 && n2 <= 1) {
#     sdObsCombined <- stats::sd(c(x1, x2))
#     sdObs1 <- 0
#     sdObs2 <- 0
#     nu1 <- 0
#     nu2 <- 0
#   } else {
#     nuCombined <- nCombined - 1
#     sdObsCombined <- sqrt(
#       ((n1-1)*sdObs1^2+n1*x1^2+(n2-1)*sdObs2^2+n2*x2^2-nCombined*xCombined^2)/nuCombined
#     )
#   }
#
#   logBf10 <- 1/2*log((1+nCombined*g0)/((1+n1*g1)*(1+n2*g2))) +
#     (nCombined/2+aGamma)*log(
#       (2*bGamma + nuCombined*sdObsCombined^2+nCombined/(1+nCombined*g0)*(xCombined-a0)^2)/
#         (2*bGamma + nu1*sdObs1^2+n1/(1+n1*g1)*(x1-a1)^2 + nu2*sdObs2^2+n2/(1+n2*g2)*(x2-a2)^2)
#     )
#
#   if (isTRUE(log))
#     return(logBf10)
#   else
#     return(exp(logBf10))
#
# }

#' Computes the credible interval of a two-sample t-test based on conjugate priors
#'
#' @inheritParams conjugateBfTStat
#' @inheritParams computeConfidenceIntervalT
#'
#' @return a vector of length two representing the credible interval
#'
#' @references
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' computeConjugateCredibleIntervalTwoSampleT(1, 1, 3, 1, 1, 3)
computeConjugateCredibleIntervalTwoSampleT <- function(
    x1, sdObs1, n1, x2, sdObs2, n2,
    # a1=3.98, g1=0.05, a2=4.02, g2=0.02,
    # a1=3.98, g1=0.5, a2=4.02, g2=0.1,
    a1=3.98, g1=0.03,
    a2=4.02, g2=0.05,
    aGamma=2, bGamma=1/2, ciValue=0.95) {

  # posterior mean conditional on sigma
  u <- (n1*g1*x1+a1)/(1+n1*g1)-(n2*g2*x2+a2)/(1+n2*g2)
  nP <- n1+n2

  w <- sqrt(
    (g1+g2+g1*g2*nP)/((1+g1*n1)*(1+g2*n2))*
    (2*bGamma+(n1-1)*sdObs1^2+(n2-1)*sdObs2^2)/(nP+2*aGamma))

  rightQuantile <- abs(
    stats::qt((1-ciValue)/2, df=nP+2*aGamma, ncp=0, lower.tail=TRUE))

  lowerCS <- u-w*rightQuantile
  upperCS <- u+w*rightQuantile

  return(unname(c(lowerCS, upperCS)))
}


#' Internal function to solve the smallest width of an eGauss t-test
#'
#' @inheritParams designSaviT
#' @inheritParams saviTTestStat
#' @param g prior variance of the eGauss t-test
#'
#' @return a number that should be zero when g is optimal
#'
#' @references
#'   `r addCite(ly2024safe)`
#'
#' @export
#'
#' @examples
#' tTestWidthDerivative(1, 1, 1)
tTestWidthDerivative <- function(g, nEff, nu, alpha=0.05) {
  nEff*nu*(alpha^(2/(nu+1))*(1+g*nEff)^(-1/(nu+1))*(1+g*nEff+nu)-1-nu)
}
