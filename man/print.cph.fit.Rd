\name{print.cph.fit}
\alias{print.cph.fit}
\title{Print cph.fit}

\description{
Formatted printing of an object of class \code{cph.fit} created
by \code{cph.fit} (which is usually called by \code{cph}). Most of the
logic for \code{print.cph.fit} came from Therneau's \code{print.coxreg}.
}
\usage{
\method{print}{cph.fit}(x, table=TRUE, coef=TRUE, conf.int=FALSE, scale=1, digits=NULL,\dots)
}
\arguments{
\item{x}{
object created by \code{cph.fit}
}
\item{table}{
print table of frequencies of events, by strata (if any)
}
\item{coef}{
print coefficient estimates, standard errors, and z-statistics
}
\item{conf.int}{
set to e.g. .95 to print 0.95 confidence intervals. Default is \code{FALSE}
to suppress confidence intervals
}
\item{scale}{
constant by which to multiply coefficients and standard errors if
printing confidence intervals.
}
\item{digits}{
number of significant digits to print
}
\item{\dots}{ignored}
}
\seealso{
\code{\link[survival]{coxph}}
}
\keyword{print}
