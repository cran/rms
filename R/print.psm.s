print.psm <- function(x, correlation = FALSE, ...)
{
  if (length(z <- x$fail) && z)
    {
      warning(" psm failed.", x$fail, "   No summary provided\n")
      return(invisible(x))
    }

  dist <- x$dist
  name <- survreg.distributions[[dist]]$name
  cat("Parametric Survival Model:",name,"Distribution\n\n")
  dput(x$call)
  cat("\n")
  if(length(x$nmiss))
    {   #backward compatibility
      cat("Frequencies of Missing Values Due to Each Variable\n")
      print(x$nmiss)
      cat("\n")
    }
  else if(length(z <- x$na.action)) naprint(z)
  stats <- x$stats
  stats[3] <- round(stats[3],2)
  stats[5] <- round(stats[5],4)
  stats[6] <- round(stats[6],2)
  print(stats)
  cat("\n")

  summary.survreg <- survival:::summary.survreg
  if(!x$fail) x$fail <- NULL    # summary.survreg uses NULL for OK
  s <- summary.survreg(x, correlation=correlation)
  print.summary.survreg2(s, correlation=correlation)
  return(invisible())

    
  wt <- x$weights
  fparms <- x$fixed
  coef <- c(x$coef, x$parms[!fparms])
  resid <- x$residuals
  dresid <- x$dresiduals
  n <- length(resid)
  p <- x$rank
  if(!length(p)) p <- sum(!is.na(coef))
  if(!p)
    {
      warning("This model has zero rank --- no summary is provided")
      return(x)
    }
  nsingular <- length(coef) - p
  rdf <- x$df.resid
  if(!length(rdf))
    rdf <- n - p
  R <- x$R   #check for rank deficiencies
  if(p < max(dim(R)))
    R <- R[1:p,     #coded by pivoting
           1:p]
  if(length(wt))
    {
      wt <- wt^0.5
      resid <- resid * wt
      excl <- wt == 0
      if(any(excl))
        {
          warning(paste(sum(excl), 
                        "rows with zero weights not counted"))
          resid <- resid[!excl]
          if(!length(x$df.residual))
            rdf <- rdf - sum(excl)
        }
    }
  famname <- x$family["name"]
  if(!length(famname)) famname <- "Gaussian"
  scale <- x$fparms
  nas <- is.na(coef)
  cnames <- names(coef[!nas])
  coef <- matrix(rep(coef[!nas], 4), ncol = 4)
  dimnames(coef) <- list(cnames, c("Value", "Std. Error", "z value", "p"))
  stds <- sqrt(diag(x$var[!nas,!nas,drop=FALSE]))
  coef[, 2] <- stds
  coef[, 3] <- coef[, 1]/stds
  coef[, 4] <- 2*pnorm(-abs(coef[,3]))
  if(correlation)
    {
      if(sum(nas)==1) ss <- 1/stds else ss <- diag(1/stds)
      correl <- ss %*% x$var[!nas, !nas, drop=FALSE] %*% ss
      dimnames(correl) <- list(cnames, cnames)
    }
  else
    correl <- NULL
  ocall <- x$call
  if(length(form <- x$formula))
    {
      if(!length(ocall$formula))
        ocall <- match.call(get("survreg"), ocall)
      ocall$formula <- form
    }
  dig <- .Options$digits
  survival:::print.summary.survreg(
                        list(call = ocall, terms = x$terms, coefficients = coef,
                             df = c(p, rdf), deviance.resid = dresid,
                             var=x$var, correlation = correl, deviance = deviance(x),
                             null.deviance = x$null.deviance, loglik=x$loglik,
                             iter = x$iter,
                             nas = nas))
  options(digits=dig)   #recovers from bug in print.summary.survreg
  invisible()
}

## Mod of print.summary.survreg from survival5 - suppresses printing a
## few things, added correlation arg

print.summary.survreg2 <-
  function (x, digits = max(options()$digits - 4, 3),
            correlation=FALSE, ...) 
  {
    correl <- x$correl
    n <- x$n
    if (is.null(digits)) 
      digits <- options()$digits
    print(x$table, digits = digits)
    if (nrow(x$var) == length(x$coefficients)) 
      cat("\nScale fixed at", format(x$scale, digits = digits), 
          "\n")
    else
      if (length(x$scale) == 1) 
        cat("\nScale=", format(x$scale, digits = digits), "\n")
      else
        {
          cat("\nScale:\n")
          print(x$scale, digits = digits, ...)
        }

    if (correlation && length(correl))
      {
        p <- dim(correl)[2]
        if (p > 1)
          {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits = digits))
            correl[!ll] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
          }
      }
    cat("\n")
    invisible(NULL)
  }
