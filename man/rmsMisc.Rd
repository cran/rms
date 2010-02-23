\name{rmsMisc}
\alias{rmsMisc}
\alias{calibrate.rms}
\alias{vcov.cph}
\alias{vcov.Glm}
\alias{vcov.Gls}
\alias{vcov.lrm}
\alias{vcov.ols}
\alias{vcov.psm}
\alias{oos.loglik}
\alias{oos.loglik.ols}
\alias{oos.loglik.lrm}
\alias{oos.loglik.cph}
\alias{oos.loglik.psm}
\alias{oos.loglik.Glm}
\alias{num.intercepts}
\alias{Getlim}
\alias{Getlimi}
\alias{related.predictors}
\alias{interactions.containing}
\alias{param.order}
\alias{Penalty.matrix}
\alias{Penalty.setup}
\alias{lrtest}
\alias{univarLR}
\alias{Newlabels}
\alias{Newlevels}
\alias{Newlabels.rms}
\alias{Newlevels.rms}
\alias{rmsArgs}
\alias{rmsFit}
\alias{print.rms}
\alias{residuals.rms}
\alias{print.lrtest}
\alias{survest.rms}
\title{Miscellaneous Design Attributes and Utility Functions}
\description{
These functions are used internally to \code{anova.rms},
\code{fastbw}, etc., to retrieve various attributes of a design.  These
functions allow some fitting functions not in the \code{rms} series
(e.g,, \code{lm}, \code{glm}) to be used with \code{rms.Design},
\code{fastbw}, and similar functions. 


For \code{vcov}, there are these functions: \code{vcov.default},
\code{vcov.lm}, \code{vcov.glm}.  The \code{oos.loglik} function for
each type of model implemented computes the -2 log likelihood for
out-of-sample data (i.e., data not necessarily used to fit the model)
evaluated at the parameter estimates from a model fit.  Vectors for the
model's linear predictors and response variable must be given.
\code{oos.loglik} is used primarily by \code{bootcov}.


The \code{Getlim} function retrieves distribution summaries
from the fit or from a \code{datadist} object.  It handles getting summaries
from both sources to fill in characteristics for variables that were not
defined during the model fit.  \code{Getlimi} returns the summary
for an individual model variable.  


The \code{related.predictors} function
returns a list containing variable numbers that are directly or
indirectly related to each predictor.  The \code{interactions.containing}
function returns indexes of interaction effects containing a given
predictor.  The \code{param.order} function returns a vector of logical
indicators for whether parameters are associated with certain types of
effects (nonlinear, interaction, nonlinear interaction).


The \code{Penalty.matrix} function builds a default penalty matrix for
non-intercept term(s) for use in penalized maximum likelihood
estimation.  The \code{Penalty.setup} function takes a constant or list
describing penalty factors for each type of term in the model and
generates the proper vector of penalty multipliers for the current model.


The \code{lrtest} function does likelihood ratio tests for
two nested models, from fits that have \code{stats} components with
\code{"Model L.R."} values.  For models such as \code{psm, survreg, ols, lm} which have
scale parameters, it is assumed that scale parameter for the smaller model
is fixed at the estimate from the larger model (see the example).

\code{univarLR} takes a multivariable model fit object from
\code{rms} and re-fits a sequence of models containing one predictor
at a time.  It prints a table of likelihood ratio \eqn{chi^2} statistics
from these fits.

The \code{Newlabels} function is used to override the variable labels in a
fit object.  Likewise, \code{Newlevels} can be used to create a new fit object
with levels of categorical predictors changed.  These two functions are
 especially useful when constructing nomograms.


\code{rmsFit} is used to convert a fit from non-rms functions (e.g.,
\code{glm}) that were invoked with rms in effect to rms functions so
that \code{anova.rms} will be called by \code{anova()}, etc.  So that the
original fit's \code{residuals} and \code{print} methods, if they exist, will be
called, there are functions \code{print.rms} and \code{residuals.rms} to
dispatch them.

\code{rmsArgs} handles \dots arguments to functions such as
\code{Predict}, \code{summary.rms}, \code{nomogram} so that variables to
vary may be specified without values (after an equals sign).
}
\usage{
\method{vcov}{cph}(object, regcoef.only=TRUE, \dots)
\method{vcov}{Glm}(object, \dots)
\method{vcov}{Gls}(object, \dots)
\method{vcov}{lrm}(object, regcoef.only=TRUE, \dots)
\method{vcov}{ols}(object, regcoef.only=TRUE, \dots)
\method{vcov}{psm}(object, regcoef.only=TRUE, \dots)

oos.loglik(fit, \dots)

\method{oos.loglik}{ols}(fit, lp, y, \dots)
\method{oos.loglik}{lrm}(fit, lp, y, \dots)
\method{oos.loglik}{cph}(fit, lp, y, \dots)
\method{oos.loglik}{psm}(fit, lp, y, \dots)
\method{oos.loglik}{Glm}(fit, lp, y, \dots)

num.intercepts(fit)

Getlim(at, allow.null=FALSE, need.all=TRUE)
Getlimi(name, Limval, need.all=TRUE)

related.predictors(at, type=c("all","direct"))
interactions.containing(at, pred)
param.order(at, term.order)

Penalty.matrix(at, X)
Penalty.setup(at, penalty)

lrtest(fit1, fit2)
\method{print}{lrtest}(x, \dots)

univarLR(fit)

Newlabels(fit, \dots)
Newlevels(fit, \dots)
\method{Newlabels}{rms}(fit, labels, \dots)
\method{Newlevels}{rms}(fit, levels, \dots)

rmsFit(fit)  # fit from glm, lm, etc.,then use anova etc. on result
}
\arguments{
\item{fit}{result of a fitting function}
\item{object}{result of a fitting function}
\item{regcoef.only}{For fits such as parametric survival models
  which have a final row and column of the covariance matrix for a
  non-regression parameter such as a log(scale) parameter, setting
  \code{regcoef.only=TRUE} causes only the first 
  \code{p} rows and columns of the covariance matrix to be returned,
  where \code{p} is the length of \code{object$coef}.
}
\item{at}{
\code{Design} element of a fit
}
\item{pred}{
index of a predictor variable (main effect)
}
\item{fit1}{
}
\item{fit2}{
fit objects from \code{lrm,ols,psm,cph} etc.  It doesn't matter which
fit object is the sub-model.
}
\item{lp}{
linear predictor vector for \code{oos.loglik}.  For proportional odds
ordinal logistic models, this should have used the first intercept
only.  If \code{lp} and \code{y} are omitted, the -2 log likelihood for the
original fit are returned.
}
\item{y}{
values of a new vector of responses passed to \code{oos.loglik}.
}
\item{name}{
the name of a variable in the model
}
\item{Limval}{
an object returned by \code{Getlim}
}
\item{allow.null}{
prevents \code{Getlim} from issuing an error message if no limits are found
in the fit or in the object pointed to by \code{options(datadist=)}
}
\item{need.all}{
set to \code{FALSE} to prevent \code{Getlim} or \code{Getlimi} from issuing an error message
if data for a variable are not found
}
\item{type}{
set to \code{"direct"} to return lists of indexes of directly related
factors only (those in interactions with the predictor)
}
\item{term.order}{
1 for all parameters, 2 for all parameters associated with either nonlinear
or interaction effects, 3 for nonlinear effects (main or interaction),
4 for interaction effects, 5 for nonlinear interaction effects.
}
\item{X}{
a design matrix, not including columns for intercepts
}
\item{penalty}{
a vector or list specifying penalty multipliers for types of model terms
}
\item{x}{a result of \code{lrtest}}
\item{labels}{
a character vector specifying new labels for variables in a fit.
To give new labels for all variables, you can specify \code{labels} of the
form \code{labels=c("Age in Years","Cholesterol")}, where the list of new labels is
assumed to be the length of all main effect-type variables in the fit and
in their original order in the model formula.  You may specify a named
vector to give new labels in random order or for a subset of the 
variables, e.g., \code{labels=c(age="Age in Years",chol="Cholesterol")}.
}
\item{levels}{
a list of named vectors specifying new level labels for categorical
predictors.  This will override \code{parms} as well as \code{datadist} information
(if available) that were stored with the fit.  
}
\item{\dots}{other arguments}
}
\value{
\code{vcov} returns a variance-covariance matrix, and \code{num.intercepts}
returns an integer with the number of intercepts in the model.
\code{oos.loglik} returns a scalar -2 log likelihood value.
\code{Getlim} returns a list with components \code{limits} and \code{values}, either
stored in \code{fit} or retrieved from the object created by \code{datadist} and
pointed to in \code{options(datadist=)}.
\code{related.predictors} returns a list of vectors, and \code{interactions.containing}
returns a vector.  \code{param.order} returns a logical vector corresponding
to non-strata terms in the model.
\code{Penalty.matrix} returns a symmetric matrix with dimension equal to the
number of slopes in the model.  For all but categorical predictor main
effect elements, the matrix is diagonal with values equal to the variances
of the columns of \code{X}.  For segments corresponding to \code{c-1} dummy variables
for \code{c}-category predictors,  puts a \code{c-1} x \code{c-1} sub-matrix in
\code{Penalty.matrix} that is constructed so that a quadratic form with 
\code{Penalty.matrix} in the middle computes the sum of squared differences
in parameter values about the mean, including a portion for the reference
cell in which the parameter is by definition zero.
\code{Newlabels} returns a new fit object with the labels adjusted.
\code{rmsFit} returns the original object but with \code{oldClass} of
\code{"rms"} and with a new attribute \code{"fitFunction"} containing the
original vector of classes.
}
\seealso{
\code{\link{rms}}, \code{\link{fastbw}}, \code{\link{anova.rms}}, \code{\link{summary.lm}}, \code{\link{summary.glm}}, \code{\link{datadist}}, \code{\link{vif}}, \code{\link{bootcov}}
}
\examples{
\dontrun{
f <- psm(S ~ x1 + x2 + sex + race, dist='gau')
g <- psm(S ~ x1 + sex + race, dist='gau', 
         fixed=list(scale=exp(f$parms)))
lrtest(f, g)


g <- Newlabels(f, c(x2='Label for x2'))
g <- Newlevels(g, list(sex=c('Male','Female'),race=c('B','W')))
nomogram(g)
}
}
\keyword{models}
\keyword{methods}






