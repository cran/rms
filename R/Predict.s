Predict <-
  function(x, ...,
           fun, type=c("predictions","model.frame","x"),
           np=200,
           conf.int=.95, conf.type=c('mean','individual','simultaneous'),
           adj.zero=FALSE, ref.zero=FALSE,
           non.slopes, time=NULL, loglog=FALSE, digits=4, name, factors=NULL)
{

  fit       <- x
  type      <- match.arg(type)
  conf.type <- match.arg(conf.type)

  oldopt <- options(digits=digits)
  on.exit(options(oldopt))

  dotlist <- if(length(factors)) factors else rmsArgs(substitute(list(...)))
  fname   <- if(missing(name)) '' else name
  at      <- fit$Design
  assume  <- at$assume.code
  name    <- at$name	##interactions are placed at end by design

  if(any(name == 'time'))
    {
      dotlist$time <- time
      time <- NULL
    }

  if(length(fname)>1 || (length(dotlist)==0 && fname==''))
    {
      m <- match.call(expand=FALSE)
      m[[1]] <- as.name('Predict')
      nams <- if(length(fname)>1) fname else name[assume!=9]
      res <- vector('list', length(nams))
      names(res) <- nams

      i <- 0
      info <- NULL # handles case where nams is empty, when no predictors
      for(nam in nams)
        {
          i <- i + 1
          m$name <- nam
          lv     <- eval(m)
          j <- attr(lv, 'info')
          if(i==1) info <- j
          else
            {
              info$varying <- c(info$varying, j$varying)
              info$adjust  <- c(info$adjust, j$adjust)
            }
          attr(lv, 'info') <- NULL
          lv$.predictor. <- nam
          res[[nam]] <- lv
        }
      lv <- do.call('rbind.data.frame', res)
      class(lv) <- c('Predict', 'data.frame')
      attr(lv, 'info') <- info
      return(lv)
    }

  f <- sum(assume!=9)	##limit list to main effects factors
  parms  <- at$parms
  label  <- at$label
  values <- at$values
  yunits <- fit$units
  units  <- at$units
  scale  <- fit$scale.pred
  if(!length(scale)) scale <- "X * Beta"

  if(missing(fun))
    {
      ylab <- scale[1]
      if(length(time))
        ylab <- ylabPlotmath <-
          if(loglog) paste("log[-log S(",format(time),")]",sep="")
          else paste(format(time),yunits,"Survival Probability")
      else if(scale[1] == 'X * Beta')
        ylabPlotmath <- expression(X*hat(beta))
      else ylabPlotmath <- ylab
    }
  else ylab <- ylabPlotmath <- ''

  if(ref.zero & length(time))
    stop("ref.zero=TRUE does not make sense with time given")

  if(fname=='') factors <- dotlist ## name not an argument
  else 
    {
      factors <- list()
      for(g in fname) factors[[g]] <- NA
    }
  nf <- length(factors)
  fnam <- names(factors)
  
  if(nf < 1) stop("must specify predictors to vary")

  which <- charmatch(fnam, name, 0)
  if(any(which==0))
    stop(paste("predictors(s) not in model:",
               paste(names(factors)[which==0],collapse=" ")))

  if(any(assume[which]==9))
    stop("cannot plot interaction terms")

  lim <- Getlim(at, allow.null=TRUE, need.all=FALSE)

  fnam   <- names(factors)
  nf     <- length(factors)
  xadjdf <- lim$limits[2,,drop=FALSE]
  xadj   <- unclass(xadjdf)
  varying <- NULL

  if(nf==0) return(as.data.frame(xadj))
  if(nf < f)  ## number of arguments < total number of predictors
    {
      ## Some non-varying predictors
      settings <- xadj
      if(adj.zero) for(x in names(settings))
        {
          i <- match(x, name)
          settings[[x]] <- if(assume[i] %in% c(5,8)) parms[[i]][1]
          else if(length(V <- lim$values[[x]]) & is.character(V)) V[1]
          else 0
        }
      for(n in fnam) settings[[n]] <- factors[[n]]
    }
  else settings <- factors

  for(i in 1:nf)
    {
      n <- fnam[i]
      v <- settings[[n]]
      lv <- length(v)
      if(lv == 0) stop('a predictor has zero length')
      if(lv == 1 && is.na(v))
        settings[[n]] <- value.chk(at, which(name==n), NA, np, lim)
      if(length(settings[[n]]) > 1) varying <- c(varying, n)
    }
  if(prod(sapply(settings,length)) > 1e5)
    stop('it is not sensible to predict more than 100,000 combinations')
  settings <- expand.grid(settings)
  adjust <- NULL
  for(n in name[assume != 9 & name %nin% fnam])
    adjust <- paste(adjust, n, "=", 
                    if(is.factor(xadj[[n]])) as.character(xadj[[n]])
                    else format(xadj[[n]])," ",sep="")

  j <- assume != 9
  label <- label[j]
  units <- units[j]
  assume <- assume[j]
  names(label) <- names(units) <- names(assume) <- name[j]
  at <- list(label=label, units=units, assume.code=assume)
  
  info <- list(varying=varying, adjust=adjust, Design=at,
               ylabPlotmath=ylabPlotmath, ylab=ylab, yunits=yunits,
               ref.zero=ref.zero, adj.zero=adj.zero, time=time,
               conf.int=conf.int)
  if(type=='model.frame')
    {
      attr(settings, 'info') <- info
      return(settings)
    }
  ##Number of non-slopes
  nrp <- num.intercepts(fit)
  if(missing(non.slopes))
    {
      non.slopes <- rep(0, nrp)
      if(!adj.zero) non.slopes[1] <- 1
    }
  if(nrp>0 && length(non.slopes) != nrp)
    stop("wrong # values in non.slopes")

  beta <- fit$coefficients
  bootdone <- length(boot.Coef <- fit$boot.Coef)
  if(bootdone && (conf.type == 'individual'))
    stop('conf.type="individual" not compatible with bootcov with coef.reps=TRUE')
  
  if(!length(time))
    {
      xx <- predictrms(fit, settings, non.slopes=non.slopes,
                       conf.int=conf.int, conf.type=conf.type,
                       ref.zero=ref.zero)
      if(length(attr(xx,"strata")) && any(is.na(attr(xx,"strata"))))
        warning("Computed stratum NA.  Requested stratum may not\nexist or reference values may be illegal strata combination\n")

      if(length(xx)==0)
        stop("model has no covariables and survival not plotted")
      xb <- if(is.list(xx)) xx$linear.predictors else xx
      if(bootdone)
        {
          X <- predictrms(fit, settings, non.slopes=non.slopes,
                          ref.zero=ref.zero, type='x')
          pred <- X %*% t(boot.Coef)
          lim  <- apply(pred, 1, quantile,
                        probs=c((1-conf.int)/2,
                              1-(1-conf.int)/2), na.rm=TRUE)
          xx$lower <- lim[1,]
          xx$upper <- lim[2,]
        }
    }
  else   ## time specified
    {
      if(bootdone) stop('time may not be specified if bootcov was used with coef.reps=TRUE')
      xx <- survest(fit, settings, times=time, loglog=loglog, 
                    conf.int=conf.int)
      xb <- as.vector(xx$surv)
    }
  if(conf.int > 0)
    {
      lower <- as.vector(xx$lower)
      upper <- as.vector(xx$upper)
    }
      
  if(!missing(fun))
    {
      xb <- fun(xb)
      if(conf.int > 0)
        {
          lower <- fun(lower)
          upper <- fun(upper)
        }	
    }

  
  settings$yhat   <- xb
  if(conf.int > 0)
    {
      settings$lower  <- lower
      settings$upper  <- upper
    }
  class(settings) <- c('Predict', 'data.frame')
  attr(settings, 'info') <- info
  settings
}

print.Predict <- function(x, ...)
{
  print.data.frame(x)
  info <- attr(x, 'info')
  cat('\nResponse variable (y):', info$ylab,'\n')
  if(length(info$adjust)==1)
    cat('\nAdjust to:',info$adjust,'\n')
  ci <- info$conf.int
  if(ci > 0)
    cat('\nLimits are', ci, 'confidence limits\n')
  invisible()
}

perimeter <- function(x, y, xinc=diff(range(x))/10, n=10,
                      lowess.=TRUE)
{

  s <- !is.na(x+y)
  x <- x[s]
  y <- y[s]
  m <- length(x)
  if(m<n)
    stop("number of non-NA x must be >= n")

  i <- order(x)
  x <- x[i]
  y <- y[i]
  s <- n:(m-n+1)
  x <- x[s]
  y <- y[s]

  x <- round(x/xinc)*xinc

  g <- function(y, n)
    {
      y <- sort(y)
      m <- length(y)
      if(n>(m-n+1)) c(NA,NA)
      else c(y[n], y[m-n+1])
    }

  r <- unlist(tapply(y, x, g, n=n))
  i <- seq(1, length(r), by=2)
  rlo <- r[i]
  rhi <- r[-i]
  s <- !is.na(rlo+rhi)
  if(!any(s))
    stop("no intervals had sufficient y observations")

  x <- sort(unique(x))[s]
  rlo <- rlo[s]
  rhi <- rhi[s]
  if(lowess.)
    {
      rlo <- lowess(x, rlo)$y
      rhi <- lowess(x, rhi)$y
    }

  structure(cbind(x, rlo, rhi),
            dimnames=list(NULL, c("x","ymin","ymax")),
            class='perimeter')
}

rbind.Predict <- function(..., rename)
  {
    d <- list(...)
    ns <- length(d)
    if(ns==1) return(d[[1]])
    
    info <- attr(d[[1]], 'info')

    if(!missing(rename))
      {
        trans <- function(input, rename)
          {
            k <- input %in% names(rename)
            if(any(k)) input[k] <- rename[input[k]]
            input
          }
        info$varying <- trans(info$varying, rename)
        names(info$Design$label) <- trans(names(info$Design$label), rename)
        names(info$Design$units) <- trans(names(info$Design$units), rename)
        names(info$Design$assume.code) <-
          trans(names(info$Design$assume.code), rename)
      }
    
    info$Design$label <- c(info$Design$label, .set.='Set')
    info$Design$units <- c(info$Design$units, .set.='')
    info$varying <- c(info$varying, '.set.')

    sets <- names(d)
    if(!length(sets)) sets <- paste('Set', 1:ns)
    obs.each.set <- sapply(d, function(x) length(x[[1]]))
    .set. <- rep(sets, obs.each.set)
    .set. <- factor(.set., levels=unique(.set.))

    if(!missing(rename)) for(i in 1:ns)
      names(d[[i]]) <- trans(names(d[[i]]), rename)

    result <- do.call('rbind.data.frame', d)
    result$.set. <- .set.
    attr(result, 'info') <- info
    class(result) <- c('Predict', 'data.frame')
    result
  }
