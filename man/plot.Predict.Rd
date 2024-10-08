\name{plot.Predict}
\alias{plot.Predict}
\alias{pantext}
\title{Plot Effects of Variables Estimated by a Regression Model Fit}

\description{
  Uses \code{lattice} graphics to plot the effect of one or two predictors
  on the linear predictor or X beta scale, or on some transformation of
  that scale.  The first argument specifies the result of the
  \code{Predict} function.  The predictor is always plotted in its
  original coding.  \code{plot.Predict} uses the
  \code{xYplot} function unless \code{formula} is omitted and the x-axis
  variable is a factor, in which case it reverses the x- and y-axes and
  uses the \code{Dotplot} function.

  If \code{data} is given, a rug plot is drawn showing
  the location/density of data values for the \eqn{x}-axis variable.  If
  there is a \code{groups} (superposition) variable that generated separate
  curves, the data density specific to each class of points is shown.
  This assumes that the second variable was a factor variable.  The rug plots
  are drawn by \code{scat1d}.  When the same predictor is used on all
  \eqn{x}-axes, and multiple panels are drawn, you can use
  \code{subdata} to specify an expression to subset according to other
  criteria in addition.

  To plot effects instead of estimates (e.g., treatment differences as a
  function of interacting factors) see \code{contrast.rms} and
  \code{summary.rms}.

  \code{pantext} creates a \code{lattice} panel function for including
  text such as that produced by \code{print.anova.rms} inside a panel or
  in a base graphic.
}
\usage{
\method{plot}{Predict}(x, formula, groups=NULL,
     cond=NULL, varypred=FALSE, subset,
     xlim, ylim, xlab, ylab, 
     data=NULL, subdata, anova=NULL, pval=FALSE, cex.anova=.85,
     col.fill=gray(seq(.825, .55, length=5)),
     adj.subtitle, cex.adj, cex.axis, perim=NULL, digits=4, nlevels=3,
     nlines=FALSE, addpanel, scat1d.opts=list(frac=0.025, lwd=0.3),
     type=NULL, yscale=NULL, scaletrans=function(z) z, ...)

pantext(object, x, y, cex=.5, adj=0, fontfamily="Courier", lattice=TRUE)

}
\arguments{
  \item{x}{a data frame created by \code{Predict}, or for \code{pantext}
  the x-coordinate for text}
\item{formula}{
  the right hand side of a \code{lattice} formula reference variables in
  data frame \code{x}.  You may not specify \code{formula} if you varied
  multiple predictors separately when calling \code{Predict}.
  Otherwise, when \code{formula} is not given, \code{plot.Predict}
  constructs one from information in \code{x}.
}
\item{groups}{an optional name of one of the variables in \code{x} that
  is to be used as a grouping (superpositioning) variable.  Note that
  \code{groups} does not contain the groups data as is customary in
  \code{lattice}; it is only a single character string specifying the
  name of the grouping variable.}
\item{cond}{when plotting effects of different predictors, \code{cond}
  is a character string that specifies a single variable name in
  \code{x} that can be used to form panels.  Only applies if using
  \code{rbind} to combine several \code{Predict} results.}
\item{varypred}{set to \code{TRUE} if \code{x} is the result of
  passing multiple \code{Predict} results, that represent different
  predictors, to \code{rbind.Predict}.  This will cause the \code{.set.}
  variable created by \code{rbind} to be copied to the
  \code{.predictor.} variable.}
\item{subset}{a subsetting expression for restricting the rows of
  \code{x} that are used in plotting.  For example, predictions may have
  been requested for males and females but one wants to plot only females.}
\item{xlim}{
This parameter is seldom used, as limits are usually controlled with
\code{Predict}.  One reason to use \code{xlim} is to plot a
\code{factor} variable on the x-axis that was created with the \code{cut2} function
with the \code{levels.mean} option, with \code{val.lev=TRUE} specified to \code{plot.Predict}. 
In this case you may want the axis to
have the range of the original variable values given to \code{cut2} rather
than the range of the means within quantile groups.
}
\item{ylim}{
Range for plotting on response variable axis. Computed by default.
}
\item{xlab}{
Label for \code{x}-axis. Default is one given to \code{asis, rcs}, etc.,
which may have been the \code{"label"} attribute of the variable.
}
\item{ylab}{
Label for \code{y}-axis.  If \code{fun} is not given,
default is \code{"log Odds"} for
\code{lrm}, \code{"log Relative Hazard"} for \code{cph}, name of the response
variable for \code{ols}, \code{TRUE} or \code{log(TRUE)} for \code{psm}, or \code{"X * Beta"} otherwise.
}
\item{data}{a data frame containing the original raw data on which the
  regression model were based, or at least containing the \eqn{x}-axis
  and grouping variable.  If \code{data} is present and contains the
  needed variables, the original data are added to the graph in the form
  of a rug plot using \code{scat1d}.
}
\item{subdata}{if \code{data} is specified, an expression to be
  evaluated in the \code{data} environment that evaluates to a logical
  vector specifying which observations in \code{data} to keep.  This
  will be intersected with the criterion for the \code{groups}
  variable.  Example: if conditioning on two paneling variables using
  \code{|a*b} you can specify
  \code{subdata=b==levels(b)[which.packet()[2]]}, where the \code{2}
  comes from the fact that \code{b} was listed second after the
  vertical bar (this assumes \code{b} is a \code{factor} in
  \code{data}.  Another example:
  \code{subdata=sex==c('male','female')[current.row()]}.}
\item{anova}{an object returned by \code{\link{anova.rms}}.  If
	\code{anova} is specified, the overall test of association for
	predictor plotted is added as text to each panel, located at the spot
	at which the panel is most empty unless there is significant empty
	space at the top or bottom of the panel; these areas are given preference.}
\item{pval}{specify \code{pval=TRUE} for \code{anova} to include not
	only the test statistic but also the P-value}
\item{cex.anova}{character size for the test statistic printed on the panel}
\item{col.fill}{
  a vector of colors used to fill confidence bands for successive
  superposed groups.  Default is inceasingly dark gray scale.
  }
\item{adj.subtitle}{
Set to \code{FALSE} to suppress subtitling the graph with the list of
settings of non-graphed adjustment values.
}
\item{cex.adj}{
\code{cex} parameter for size of adjustment settings in subtitles.  Default is
0.75 times \code{par("cex")}.
}
\item{cex.axis}{
  \code{cex} parameter for x-axis tick labels
  }
\item{perim}{
\code{perim} specifies a function having two
arguments.  The first is the vector of values of the first variable that
is about to be plotted on the x-axis.  The second argument is the single
value of the variable representing different curves, for the current
curve being plotted.  The function's returned value must be a logical
vector whose length is the same as that of the first argument, with
values \code{TRUE} if the corresponding point should be plotted for the
current curve, \code{FALSE} otherwise.  See one of the latter examples.
If a predictor is not specified to \code{plot}, \code{NULL} is passed as
the second argument to \code{perim}, although it makes little sense to
use \code{perim} when the same \code{perim} is used for multiple predictors.
}
\item{digits}{
Controls how numeric variables used for panel labels are formatted. The
default is 4 significant digits.
}
\item{nlevels}{
  when \code{groups} and \code{formula} are not specified, if any panel
  variable has \code{nlevels} or fewer values, that variable is
  converted to a \code{groups} (superpositioning) variable.  Set
  \code{nlevels=0} to prevent this behavior.  For other situations, a
  numeric x-axis variable with \code{nlevels} or fewer unique values
  will cause a dot plot to be drawn instead of an x-y plot.
}
\item{nlines}{If \code{formula} is given, you can set \code{nlines} to
  \code{TRUE} to convert the x-axis variable to a factor and then to an
  integer.  Points are plotted at integer values on the x-axis but
  labeled with category levels.  Points are connected by lines.}
\item{addpanel}{an additional panel function to call along with panel
  functions used for \code{xYplot} and \code{Dotplot} displays}
\item{scat1d.opts}{a list containing named elements that specifies
  parameters to \code{\link[Hmisc]{scat1d}} when \code{data} is given.  The
  \code{col} parameter is usually derived from other plotting
  information and not specified by the user.}
\item{type}{a value (\code{"l","p","b"}) to override default choices
  related to showing or connecting points.  Especially  useful for
  discrete x coordinate variables.}
\item{yscale}{a \code{lattice} scale \code{list} for the \code{y}-axis
	to be added to what is automatically generated for the \code{x}-axis.
	Example:
	\code{yscale=list(at=c(.005,.01,.05),labels=format(c(.005,.01,.05)))}.
  See \link[lattice]{xyplot}}
\item{scaletrans}{a function that operates on the \code{scale} object
	created by \code{plot.Predict} to produce a modified \code{scale}
	object that is passed to the lattice graphics function.  This is
	useful for adding other \code{scales} options or for changing the
	\code{x}-axis limits for one predictor.}
\item{\dots}{
  extra arguments to pass to \code{xYplot} or \code{Dotplot}.  Some
  useful ones are \code{label.curves} and \code{abline}.
  Set \code{label.curves} to \code{FALSE} to suppress labeling of
  separate curves. Default is \code{TRUE}, which
  causes \code{labcurve} to be invoked to place labels at positions where the
  curves are most separated, labeling each curve with the full curve label.
  Set \code{label.curves} to a \code{list} to specify options to
  \code{labcurve}, e.g., \code{label.curves=} \code{list(method="arrow",
	cex=.8)}. 
  These option names may be abbreviated in the usual way arguments
  are abbreviated.  Use for example \code{label.curves=list(keys=letters[1:5])}
  to draw single lower case letters on 5 curves where they are most
  separated, and automatically position a legend
  in the most empty part of the plot.  The \code{col}, \code{lty}, and
  \code{lwd} parameters are passed automatically to \code{labcurve}
  although they may be overridden here.
  It is also useful to use \dots to pass \code{lattice} graphics parameters, e.g.
  \code{par.settings=list(axis.text=list(cex=1.2), par.ylab.text=list(col='blue',cex=.9),par.xlab.text=list(cex=1))}.
}
\item{object}{an object having a \code{print} method}
\item{y}{y-coordinate for placing text in a \code{lattice} panel
or on a base graphics plot}
\item{cex}{character expansion size for \code{pantext}}
\item{adj}{text justification.  Default is left justified.}
\item{fontfamily}{
font family for \code{pantext}.  Default is \code{"Courier"} which
will line up columns of a table.
}
\item{lattice}{set to \code{FALSE} to use \code{text} instead of
  \code{ltext} in the function generated by \code{pantext}, to use base
  graphics}
}
\value{
  a \code{lattice} object ready to \code{print} for rendering.
}
\details{
When a \code{groups} (superpositioning) variable was used, you can issue
the command \code{Key(\dots)} after printing the result of
\code{plot.Predict}, to draw a key for the groups.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\references{
  Fox J, Hong J (2009): Effect displays in R for multinomial and
  proportional-odds logit models: Extensions to the effects package.  J
  Stat Software 32 No. 1.
}
\note{If plotting the effects of all predictors you can reorder the
  panels using for example \code{p <- Predict(fit); p$.predictor. <-
	factor(p$.predictor., v)} where \code{v} is a vector of predictor
  names specified in the desired order.
  }
\seealso{
  \code{\link{Predict}}, \code{\link{ggplot.Predict}},
	\code{link{plotp.Predict}}, \code{\link{rbind.Predict}},
  \code{\link{datadist}}, \code{\link{predictrms}}, \code{\link{anova.rms}},
  \code{\link{contrast.rms}}, \code{\link{summary.rms}},
  \code{\link{rms}}, \code{\link{rmsMisc}}, 
  \code{\link[Hmisc]{labcurve}}, \code{\link[Hmisc]{scat1d}},
  \code{\link[Hmisc]{xYplot}}, \code{\link[Hmisc]{Overview}}
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
an <- anova(fit)
# Plot effects of all 4 predictors with test statistics from anova, and P
plot(Predict(fit), anova=an, pval=TRUE)
plot(Predict(fit), data=llist(blood.pressure,age))
                         # rug plot for two of the predictors

p <- Predict(fit, name=c('age','cholesterol'))   # Make 2 plots
plot(p)

p <- Predict(fit, age=seq(20,80,length=100), sex, conf.int=FALSE)
                         # Plot relationship between age and log
                         # odds, separate curve for each sex,
plot(p, subset=sex=='female' | age > 30)
# No confidence interval, suppress estimates for males <= 30

p <- Predict(fit, age, sex)
plot(p, label.curves=FALSE, data=llist(age,sex))
                         # use label.curves=list(keys=c('a','b'))'
                         # to use 1-letter abbreviations
                         # data= allows rug plots (1-dimensional scatterplots)
                         # on each sex's curve, with sex-
                         # specific density of age
                         # If data were in data frame could have used that
p <- Predict(fit, age=seq(20,80,length=100), sex='male', fun=plogis)
                         # works if datadist not used
plot(p, ylab=expression(hat(P)))
                         # plot predicted probability in place of log odds

per <- function(x, y) x >= 30
plot(p, perim=per)       # suppress output for age < 30 but leave scale alone

# Take charge of the plot setup by specifying a lattice formula
p <- Predict(fit, age, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215), sex)
plot(p, ~ age | blood.pressure*cholesterol, subset=sex=='male')
# plot(p, ~ age | cholesterol*blood.pressure, subset=sex=='female')
# plot(p, ~ blood.pressure|cholesterol*round(age,-1), subset=sex=='male')
plot(p)

# Plot the age effect as an odds ratio
# comparing the age shown on the x-axis to age=30 years

ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
fit <- update(fit)   # make new reference value take effect
p <- Predict(fit, age, ref.zero=TRUE, fun=exp)
plot(p, ylab='Age=x:Age=30 Odds Ratio',
     abline=list(list(h=1, lty=2, col=2), list(v=30, lty=2, col=2)))

# Compute predictions for three predictors, with superpositioning or
# conditioning on sex, combined into one graph

p1 <- Predict(fit, age, sex)
p2 <- Predict(fit, cholesterol, sex)
p3 <- Predict(fit, blood.pressure, sex)
p <- rbind(age=p1, cholesterol=p2, blood.pressure=p3)
plot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE)
plot(p, cond='sex', varypred=TRUE, adj.subtitle=FALSE)

\dontrun{
# For males at the median blood pressure and cholesterol, plot 3 types
# of confidence intervals for the probability on one plot, for varying age
ages <- seq(20, 80, length=100)
p1 <- Predict(fit, age=ages, sex='male', fun=plogis)  # standard pointwise
p2 <- Predict(fit, age=ages, sex='male', fun=plogis,
              conf.type='simultaneous')               # simultaneous
p3 <- Predict(fit, age=c(60,65,70), sex='male', fun=plogis,
              conf.type='simultaneous')               # simultaneous 3 pts
# The previous only adjusts for a multiplicity of 3 points instead of 100
f <- update(fit, x=TRUE, y=TRUE)
g <- bootcov(f, B=500, coef.reps=TRUE)
p4 <- Predict(g, age=ages, sex='male', fun=plogis)    # bootstrap percentile
p <- rbind(Pointwise=p1, 'Simultaneous 100 ages'=p2,
           'Simultaneous     3 ages'=p3, 'Bootstrap nonparametric'=p4)
xYplot(Cbind(yhat, lower, upper) ~ age, groups=.set.,
       data=p, type='l', method='bands', label.curve=list(keys='lines'))
}

# Plots for a parametric survival model
require(survival)
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
              rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
label(t) <- 'Follow-up Time'
e <- ifelse(t<=cens,1,0)
t <- pmin(t, cens)
units(t) <- "Year"
ddist <- datadist(age, sex)
Srv <- Surv(t,e)


# Fit log-normal survival model and plot median survival time vs. age
f <- psm(Srv ~ rcs(age), dist='lognormal')
med <- Quantile(f)       # Creates function to compute quantiles
                         # (median by default)
p <- Predict(f, age, fun=function(x) med(lp=x))
plot(p, ylab="Median Survival Time")
# Note: confidence intervals from this method are approximate since
# they don't take into account estimation of scale parameter


# Fit an ols model to log(y) and plot the relationship between x1
# and the predicted mean(y) on the original scale without assuming
# normality of residuals; use the smearing estimator
# See help file for rbind.Predict for a method of showing two
# types of confidence intervals simultaneously.
set.seed(1)
x1 <- runif(300)
x2 <- runif(300)
ddist <- datadist(x1,x2)
y  <- exp(x1+x2-1+rnorm(300))
f <- ols(log(y) ~ pol(x1,2)+x2)
r <- resid(f)
smean <- function(yhat)smearingEst(yhat, exp, res, statistic='mean')
formals(smean) <- list(yhat=numeric(0), res=r[!is.na(r)])
#smean$res <- r[!is.na(r)]   # define default res argument to function
plot(Predict(f, x1, fun=smean), ylab='Predicted Mean on y-scale')

# Make an 'interaction plot', forcing the x-axis variable to be
# plotted at integer values but labeled with category levels
n <- 100
set.seed(1)
gender <- c(rep('male', n), rep('female',n))
m <- sample(c('a','b'), 2*n, TRUE)
d <-  datadist(gender, m); options(datadist='d')
anxiety <- runif(2*n) + .2*(gender=='female') + .4*(gender=='female' & m=='b')
tapply(anxiety, llist(gender,m), mean)
f <- ols(anxiety ~ gender*m)
p <- Predict(f, gender, m)
plot(p)     # horizontal dot chart; usually preferred for categorical predictors
Key(.5, .5)
plot(p, ~gender, groups='m', nlines=TRUE)
plot(p, ~m, groups='gender', nlines=TRUE)
plot(p, ~gender|m, nlines=TRUE)

options(datadist=NULL)

\dontrun{
# Example in which separate curves are shown for 4 income values
# For each curve the estimated percentage of voters voting for
# the democratic party is plotted against the percent of voters
# who graduated from college.  Data are county-level percents.

incomes <- seq(22900, 32800, length=4)  
# equally spaced to outer quintiles
p <- Predict(f, college, income=incomes, conf.int=FALSE)
plot(p, xlim=c(0,35), ylim=c(30,55))

# Erase end portions of each curve where there are fewer than 10 counties having
# percent of college graduates to the left of the x-coordinate being plotted,
# for the subset of counties having median family income with 1650
# of the target income for the curve

show.pts <- function(college.pts, income.pt) {
  s <- abs(income - income.pt) < 1650  #assumes income known to top frame
  x <- college[s]
  x <- sort(x[!is.na(x)])
  n <- length(x)
  low <- x[10]; high <- x[n-9]
  college.pts >= low & college.pts <= high
}

plot(p, xlim=c(0,35), ylim=c(30,55), perim=show.pts)

# Rename variables for better plotting of a long list of predictors
f <- ...
p <- Predict(f)
re <- c(trt='treatment', diabet='diabetes', sbp='systolic blood pressure')

for(n in names(re)) {
  names(p)[names(p)==n] <- re[n]
  p$.predictor.[p$.predictor.==n] <- re[n]
  }
plot(p)
}
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
