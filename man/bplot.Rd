\name{bplot}
\alias{bplot}
\alias{perimeter}
\alias{lines.perimeter}
\alias{iLegend}
\title{
3-D Plots Showing Effects of Two Continuous Predictors in a Regression
Model Fit}
\description{
Uses base graphics and the output from \code{Predict} to plot image,
contour, or perspective plots showing the simultaneous effects of two
continuous predictor variables.  The \eqn{x}-axis is constructed from
the first variable listed in the call to \code{Predict} and the
\eqn{y}-axis variable comes from the second.

The \code{perimeter} function is used to generate the boundary of data
to plot when a 3-d plot is made.  It finds the area where there are
sufficient data to generate believable interaction fits.

\code{iLegend} is function for adding legends to an existing graph
produced with \code{bplot} using \code{method="image"}.
}
\usage{
bplot(x, xlab, ylab, zlab,
      adj.subtitle=TRUE, cex.adj, perim,
      method=c("image", "persp", "contour"),
      zlim=range(yhat, na.rm=TRUE), nlevels=10, \dots)

perimeter(x, y, xinc=diff(range(x))/10, n=10, lowess.=TRUE)

\method{lines}{perimeter}(x, \dots)

iLegend(object, x, y, size=c(1,1), horizontal=TRUE,
       nint=50, fun.=NULL, at=NULL, zlab, zlim, par.=NULL, \dots)
}
\arguments{
  \item{x}{
	for \code{bplot}, an object created by \code{Predict} for which
    exactly two numeric predictors varied.
	For \code{iLegend}, is either a vector of 1 or 2
	\eqn{x}-coordinates or a list with elements \code{x} and \code{y} each
	with 1 or 2 coordinates.  A single coordinate represents the upper
    left corner of the legend box.  \code{x} and \code{y} are optional.  If
	omitted, \code{locator} is used to position legends with the mouse.
	For \code{perim} is
	the first variable of a pair of predictors forming a 3-d plot.
	For \code{lines.perimeter}, \code{x} is the result of \code{perimeter}.
  }
\item{xlab}{
  Label for \eqn{x}-axis. Default is given by \code{Predict}.
}
\item{ylab}{
  Label for \eqn{y}-axis
}
\item{zlab}{
  \eqn{z}-axis label for perspective plots or for \code{iLegend}.
       Default comes from \code{Predict}.  \code{zlab} will often be
       specified if \code{fun} was specified to \code{Predict} or
       \code{iLegend}.
  }
\item{adj.subtitle}{
Set to \code{FALSE} to suppress subtitling the graph with the list of
       settings of non-graphed adjustment values. Default is \code{TRUE}
       if there are non-plotted adjustment variables.
}
\item{cex.adj}{
\code{cex} parameter for size of adjustment settings in subtitles.  Default is
0.75 times \code{par("cex")}.
}
\item{perim}{
names a matrix created by \code{perimeter} when used for 3-d plots of
two continuous predictors.  When the combination of variables is outside
the range in \code{perim}, that section of the plot is suppressed.  If
\code{perim} 
is omitted, 3-d plotting will use the marginal distributions of the
two predictors to determine the plotting region, when the grid is
not specified explicitly in \code{variables}.  When instead a series of
curves is being plotted, \code{perim} specifies a function having two
arguments.  The first is the vector of values of the first variable that
is about to be plotted on the \eqn{x}-axis.  The second argument is the single
value of the variable representing different curves, for the current
curve being plotted.  The function's returned value must be a logical
vector whose length is the same as that of the first argument, with
values \code{TRUE} if the corresponding point should be plotted for the
current curve, \code{FALSE} otherwise.  See one of the latter examples.
}
\item{method}{
Use \code{method="persp"} for perspective plots (\code{persp()}), 
  \code{method="contour"} to use \code{contour()}, or \code{method="image"}
  to use \code{image()} (the default).
}
\item{zlim}{
If 'type="persp"' controls the range for plotting in the
\eqn{z}-axis. Computed by default.  Also used by \code{iLegend}.
}
\item{nlevels}{passed to \code{contour}}
\item{par.}{a list of arguments to pass to \code{par}, such as
       \code{mgp}.  The user can use the argument name \code{par},
       omitting the period.}
\item{\dots}{
  optional arguments to pass to \code{persp}, \code{image}, or
  \code{contour}.   Specify a vector to assign different
  colors to different curves or to components of 3-d plots,
  e.g. \code{col=rainbow(30)} when \code{method='image'}. 
}
\item{y}{
  second variable of the pair for \code{perim}, or \eqn{y}-coordinates for
  \code{iLegend}.  If omitted, \code{x} is assumed to be a list with both
  \code{x} and \code{y} components.
}
\item{xinc}{
  increment in \code{x} over which to examine the density of \code{y} in
       \code{perimeter} 
}
\item{n}{
  within intervals of \code{x} for \code{perimeter}, takes the
  informative range of \code{y} to be the \eqn{n}th smallest to the
  \eqn{n}th largest values of \code{y}.  If there aren't 
  at least 2\eqn{n} \code{y} values in the \code{x} interval, no
  \code{y} ranges are used for that interval.
}
\item{lowess.}{
set to \code{FALSE} to not have \code{lowess} smooth the data perimeters
}
\item{object}{object created by \code{Predict}}
\item{size}{size of legend in inches; omitted with 2 points are given in
       \code{x} and \code{y}}
\item{horizontal}{set to \code{FALSE} to make the labeled legend axis vertical}
\item{nint}{number of intervals of predicted values to use in legend}
\item{fun.}{
  For \code{iLegend}, \code{fun} is
  a function for transforming tick mark labels for color or gray scale
  legends for \code{method="image"}.   For example, if \code{bplot} is used
  to make an image plot of log odds ratios, specifying \code{fun=plogis} will
  cause the color legend to be labeled with probability values rather
  than log odds.  The user will usually omit the . from the argument
       name \code{fun.}.
}
\item{at}{
If \code{fun} is specified to \code{iLegend}, \code{at} may be given.
\code{at} is a vector 
of values at which to evaluate \code{fun} for drawing tick marks in the
color legend.  For example, if you want to show the median survival time
for a log-normal survival model whereas the linear predictor (log median)
was used in constructing the image plot, and if you want to place tick
marks at nice median values, specify \code{fun=exp, at=log(c(1,10,100,1000))}.
}
}
\value{
\code{perimeter} returns a matrix of class \code{perimeter}.  This
       outline can be conveniently plotted by \code{lines.perimeter}.
\code{iLegend} invisibly returns the position of the legend.
}
\details{
\code{perimeter} is a kind of generalization of \code{datadist} for 2
continuous variables.  First, the \code{n} smallest and largest \code{x}
values are determined.  These form the lowest and highest possible
\code{x}s to display.  Then \code{x} is grouped into intervals bounded
by these two numbers, with the interval widths defined by \code{xinc}.
Within each interval, \code{y} is sorted and the \eqn{n}th smallest and
largest \code{y} are taken as the interval containing sufficient data
density to plot interaction surfaces.  The interval is ignored when
there are insufficient \code{y} values.  When the data are being
readied for \code{persp}, \code{bplot} uses the \code{approx} function to do
linear interpolation of the \code{y}-boundaries as a function of the
\code{x} values actually used in forming the grid (the values of the
first variable specified to \code{Predict}).  To make the perimeter smooth,
specify \code{lowess.=TRUE} to \code{perimeter}.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{datadist}}, \code{\link{Predict}},
\code{\link{rms}}, \code{\link{rmsMisc}}, \code{\link{image}},
       \code{\link{contour}}, \code{\link{persp}},
       \code{\link[Hmisc]{subplot}} 
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'   # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'

# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
p <- Predict(fit, age=., cholesterol=., np=100)
bplot(p)                 # image plot for age, cholesterol with color
                         # coming from yhat; use default ranges for
                         # both continuous predictors
bplot(p, method='persp', theta=30, ticktype='detailed')
                         # 3-d perspective plot
bplot(p, method='contour')       # contour plot
boundaries <- perimeter(age, cholesterol, lowess=TRUE)
plot(age, cholesterol)   # show bivariate data density
lines(boundaries)        # and perimeter that will be used for 3-D plot
p <- Predict(fit, age=., cholesterol=.)
par(mgp=c(1.7, .35, 0))
bplot(p, perim=boundaries)
                         # draws image() plot
                         # don't show estimates where data are sparse
                         # doesn't make sense here since vars don't interact
iLegend(p, x=c(30,40), y=c(230, 245))   # original logit scale
iLegend(p, x=c(65,75), y=c(230, 245), fun=plogis,
        at=qlogis(c(.1,.25,.5,.75,.9)),
        zlab='Probability')             # probability scale
options(datadist=NULL)
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
