survfit.formula <- function(formula, data, ...)
{
  sff <- survival:::survfit.formula
  ## Can't use as.name('survival:::survfit.formula')
  mc <- match.call()
  mc[[1L]] <- as.name('sff')
  f <- eval(mc, list(sff=survival:::survfit.formula), parent.frame())
  
  f$maxtime <- max(f$time)

  g <- if(missing(data)) model.frame(formula)
  else model.frame(formula, data=data)
  Y <- model.extract(g, 'response')
  f$units <- valueUnit(Y)
  f$time.label <- attr(Y, "time.label")
  f$call <- match.call()
  f
}
