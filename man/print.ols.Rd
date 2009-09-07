\name{print.ols}
\alias{print.ols}
\title{Print ols}
\description{
formatted printing of an object of class \code{ols} using methods taken from
\code{print.lm} and \code{summary.lm}. Prints R-squared, adjusted R-squared,
parameter estimates,
standard errors, and t-statistics (Z statistics if penalized estimation was used).
For penalized
estimation, prints the maximum penalized likelihood estimate of the residual
standard deviation (\code{Sigma}) instead of the usual root mean squared error.
}
\usage{
\method{print}{ols}(x, digits=4, long=FALSE, \dots)
}
\arguments{
\item{x}{fit object}
\item{digits}{number of significant digits to print}
\item{long}{set to \code{TRUE} to print the correlation matrix of
  parameter estimates}
\item{\dots}{other parameters to pass to \code{print} or \code{format}}
}
\seealso{
\code{\link{ols}}, \code{\link{print.lm}}, \code{\link{summary.lm}}
}
\keyword{print}
