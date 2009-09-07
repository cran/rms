\name{Rq}
\Rdversion{1.1}
\alias{Rq}
\alias{RqFit}
\alias{print.Rq}
\alias{latex.Rq}
\title{rms Package Interface to quantreg Package}
\description{
The \code{Rq} function is the \code{rms} front-end to the
\code{quantreg} package's \code{rq} function.  \code{print} and
\code{latex} methods are also provided, and a fitting function
\code{RqFit} is defined for use in bootstrapping, etc.  Its result is a
function definition.
}
\usage{
Rq(formula, tau = 0.5, data, subset, weights, na.action=na.delete,
   method = "br", model = FALSE, contrasts = NULL,
   se = "nid", hs = TRUE, x = FALSE, y = FALSE, ...)

\method{print}{Rq}(x, digits=4, \dots)

\method{latex}{Rq}(object,
           file = paste(first.word(deparse(substitute(object))),
             ".tex", sep = ""), append=FALSE,
           which, varnames, columns=65, inline=FALSE, caption=NULL,
           ...)

RqFit(fit, wallow=TRUE, passdots=FALSE)
}
\arguments{
  \item{formula}{model formula}
  \item{tau}{
	the single quantile to estimate.  Unlike \code{rq} you cannot estimate
	more than one quantile at one model fitting.
  }
  \item{data}{}
  \item{subset}{}
  \item{weights}{}
  \item{na.action}{}
  \item{method}{}
  \item{model}{}
  \item{contrasts}{}
  \item{se}{}
  \item{hs}{see \code{\link[quantreg]{rq}}}
  \item{x}{set to \code{TRUE} to store the design matrix with the fit.
	For \code{print} is an \code{Rq} object.}
  \item{y}{set to \code{TRUE} to store the response vector with the fit}
  \item{\dots}{
	other arguments passed to one of the \code{rq} fitting routines.
	For \code{latex.Rq} these are optional arguments passed to
	\code{latexrms}.  Ignored for \code{print.Rq}.
  }
  \item{digits}{
	number of significant digits used in formatting results in
	\code{print.Rq}.
  }
  \item{object}{an object created by \code{Rq}}
  \item{file}{}
  \item{append}{}
  \item{which}{}
  \item{varnames}{}
  \item{columns}{}
  \item{inline}{}
  \item{caption}{see \code{\link{latexrms}}}
  \item{fit}{an object created by \code{Rq}}
  \item{wallow}{
	set to \code{TRUE} if \code{weights} are allowed in the
	current context.
  }
  \item{passdots}{
	set to \code{TRUE} if \dots may be passed to the fitter}
}
\value{
  \code{Rq} returns a list of class \code{"rms", "lassorq"} or \code{"scadrq",
	"Rq"}, and \code{"rq"}.  \code{RqFit} returns a function
  definition.  \code{latex.Rq} returns an object of class \code{"latex"}.
}
\author{
  Frank Harrell
}
\note{
  The author and developer of methodology in the \code{quantreg} package
  is Roger Koenker.
}
\seealso{
  \code{\link[quantreg]{rq}}
}
\examples{
\dontrun{
set.seed(1)
n     <- 30
y     <- rnorm(n)
x1    <- rnorm(n)
w     <- runif(n)
x1[2] <- NA
dd <- datadist(x1); options(datadist='dd')
f <- Rq(y ~ pol(x1,2))
anova(f)
plot(Predict(f, x1=.))
}
}
\keyword{models}
\keyword{nonparametric}
