% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Glm.r
\name{print.Glm}
\alias{print.Glm}
\title{print.glm}
\usage{
\method{print}{Glm}(x, digits = 4, coefs = TRUE, title = "General Linear Model", ...)
}
\arguments{
\item{x}{`Glm` object}

\item{digits}{number of significant digits to print}

\item{coefs}{specify `coefs=FALSE` to suppress printing the table of
model coefficients, standard errors, etc.  Specify `coefs=n` to print
only the first `n` regression coefficients in the model.}

\item{title}{a character string title to be passed to `prModFit`}

\item{...}{ignored}
}
\description{
Print a `Glm` Object
}
\details{
Prints a `Glm` object, optionally in LaTeX or html
}
\author{
Frank Harrell
}
