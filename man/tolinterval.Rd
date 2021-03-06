% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tolerance.R
\name{tolinterval}
\alias{tolinterval}
\title{Tolerance interval for linear regression}
\usage{
tolinterval(fit, xnew, p = 0.9, alpha = 0.05, side = "left-sided")
}
\arguments{
\item{fit}{a \code{\link[stats]{lm}} object}

\item{xnew}{new values of the predictors for which to get the tolerance interval}

\item{p}{proportion to be covered by the tolerance interval}

\item{alpha}{significance level}

\item{side}{either \code{"left-sided"}, \code{"right-sided"} or \code{"two-sided"}}
}
\value{
The tolerance interval, a vector of length two.
}
\description{
Tolerance interval for a Gaussian linear regression model.
}
\examples{
data(KM31)
fit <- lm(y ~ x1+x2, data = KM31)
tolinterval(fit, xnew=c(1, 88, 9), side="two-sided")
}
