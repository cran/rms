\name{print.cph}
\alias{print.cph}
\title{Print cph Results}

\description{
Formatted printing of an object of class \code{cph}. Prints strata
frequencies, parameter estimates, standard errors, z-statistics, numbers
of missing values, VIFs.
}
\usage{
\method{print}{cph}(x, long=FALSE, digits=3, conf.int=FALSE, table=TRUE, \dots)
}
\arguments{
\item{x}{
fit object
}
\item{long}{
set to \code{TRUE} to print the centering constant
}
\item{digits}{
number of significant digits to print
}
\item{conf.int}{
set to e.g. .95 to print 0.95 confidence intervals on simple hazard ratios
}
\item{table}{
set tp \code{FALSE} to suppress event frequency statistics
}
\item{...}{
arguments passed to \code{print.cphfit} : \code{coef}, \code{scale}
}}
\seealso{
\code{\link[survival]{coxph}}
}
\keyword{print}

