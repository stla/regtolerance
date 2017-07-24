kLM <- function(p, alpha, ell, d){
  e <- (1+d^2)^2/d^4
  f <- d^4/(1+d^2) # faux dans K&M
  delta <- d^2*((3*d^2+sqrt(9*d^4+6*d^2+3))/(2*d^2+1))
  sqrt(e*f/(1+delta)*rcpp_qchisq(p, 1, delta)*qf(1-alpha, e, ell)) # manque sqrt dans K&M
}

solveK <- function(p, alpha, ell, d, ...){
  k_low <- sqrt(ell*qchisq(p,1)/qchisq(alpha, ell))
  k_upp <- kLM(p, alpha, ell, d) + 1
  uniroot(function(k) integral(ell, p, k, d) - (1-alpha), 
          lower=k_low, upper=k_upp, ...)
}

#' Tolerance factor for linear regression
#'
#' @param X model matrix
#' @param xnew new values of the predictors
#' @param p proportion to be covered by the tolerance interval
#' @param alpha significance level
#' @param side either \code{"one-sided"} or \code{"two-sided"}
#' @param ... arguments passed to \code{\link[base]{uniroot}} if \code{side="two-sided"}, 
#' otherwise ignored
#'
#' @return A number, the one-sided tolerance factor, if \code{side="one-sided"}, 
#' otherwise the output of \code{\link[base]{uniroot}}; the component \code{root} is 
#' the two-sided tolerance factor
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ disp+drat, data = mtcars)
#' tolfactor(model.matrix(fit), xnew=c(1, 160, 4), side="two-sided") 
#' data(KM31)
#' fit <- lm(y ~ x1+x2, data = KM31)
#' tolfactor(model.matrix(fit), xnew=c(1, 88, 9), side="one-sided")
#' tolfactor(model.matrix(fit), xnew=c(1, 88, 9), side="two-sided")$root
tolfactor <- function(X, xnew, p=0.9, alpha=0.05, side="one-sided", ...){
  side <- match.arg(side, choices = c("one-sided", "two-sided"))
  H <- chol2inv(chol(t(X) %*% X))
  d <- sqrt(c(t(xnew) %*% H %*% matrix(xnew)))
  if(side == "one-sided"){
    return(d * rcpp_qt(1-alpha, nrow(X)-ncol(X), qnorm(p)/d))
  }
  solveK(p=p, alpha=alpha, ell=nrow(X)-ncol(X), d=d, ...)
}

#' @title Tolerance interval for linear regression
#' @description Tolerance interval for a Gaussian linear regression model.
#' @param fit a \code{\link[stats]{lm}} object
#' @param xnew new values of the predictors for which to get the tolerance interval
#' @param p proportion to be covered by the tolerance interval
#' @param alpha significance level
#' @param side either \code{"left-sided"}, \code{"right-sided"} or \code{"two-sided"}
#'
#' @return The tolerance interval, a vector of length two.
#' @export
#' @importFrom stats sigma model.matrix
#'
#' @examples
#' data(KM31)
#' fit <- lm(y ~ x1+x2, data = KM31)
#' tolinterval(fit, xnew=c(1, 88, 9), side="two-sided")
tolinterval <- function(fit, xnew, p=0.9, alpha=0.05, side="left-sided"){
  side <- match.arg(side, choices = c("left-sided", "right-sided", "two-sided"))
  k <- tolfactor(model.matrix(fit), xnew=xnew, p=p, alpha=alpha, side=side)
  estimates <- fit[["coefficients"]]
  yhat <- c(t(xnew) %*% matrix(estimates))
  if(side == "two-sided"){
    k <- k$root
    return(yhat + c(-1,1)*k*sigma(fit))
  }
  if(side == "left-sided"){
    return(c(-Inf, yhat + k*sigma(fit)))
  }
  c(yhat - k*sigma(fit), Inf)
}