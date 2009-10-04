cph <- function(formula=formula(data),
                data=parent.frame(),
                weights,
                subset,
                na.action=na.delete, 
                method=c("efron","breslow","exact",
                  "model.frame", "model.matrix"),
                singular.ok=FALSE,
                robust=FALSE,
                model=FALSE,
                x=FALSE,
                y=FALSE,
                se.fit=FALSE,
                eps=1e-4,
                init,
                iter.max=10,
                tol=1e-9,
                surv=FALSE,
                time.inc,
                type,
                vartype,
                ...)
{
  require(survival)
  method <- match.arg(method)
  call <- match.call()
  m <- match.call(expand=FALSE)
  mc <- match(c("formula", "data", "subset", "weights", "na.action"), 
              names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action

  m$drop.unused.levels <- TRUE
  
  m[[1]] <- as.name("model.frame")

  if (!inherits(formula,"formula"))
    {
      ## I allow a formula with no right hand side
      ## The dummy function stops an annoying warning message "Looking for
      ##  'formula' of mode function, ignored one of mode ..."
      if (inherits(formula,"Surv")) {
        xx <- function(x) formula(x)
        
        formula <- xx(paste(deparse(substitute(formula)), 1, sep="~"))
      }
      else stop("Invalid formula")
    }

  m$formula <- formula

  nstrata <- 0
  Strata <- NULL

  if(!missing(data) ||
     (length(z <- attr(terms(formula, allowDotAsName=TRUE),"term.labels"))>0 &&
                        any(z!=".")))
    { #X's present
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul)
        {
          on.exit(options(drop.unused.levels=dul))
          options(drop.unused.levels=FALSE)
        }

      X <- Design(eval.parent(m))
      atrx <- attributes(X)
      atr  <- atrx$Design
      nact <- atrx$na.action
      if(method=="model.frame") return(X)

      Terms <- if(missing(data))
        terms(formula, specials=c("strat","cluster"))
      else
        terms(formula, specials=c("strat","cluster"), data=data)

      asm   <- atr$assume.code
      name  <- atr$name

      cluster <- attr(Terms, "specials")$cluster
      stra    <- attr(Terms, "specials")$strat

      if(length(cluster))
        {
          if(missing(robust)) robust <- TRUE
          Terms <- Terms[-(cluster - 1)]
          cluster <- attr(X, 'cluster')
          attr(X, 'cluster') <- NULL
        }

      Terms.ns <- Terms
      if(length(stra))
        {
          temp <- untangle.specials(Terms.ns, "strat", 1)
          Terms.ns <- Terms.ns[-temp$terms]	#uses [.terms function
          ##  Set all factors=2
          ## (-> interaction effect not appearing in main effect
          ##  that was deleted strata effect)
          
          Strata <- list()
          
          for(i in (1:length(asm))[asm==8])
            {
              nstrata <- nstrata+1
              xi <- X[[i+1]]
              levels(xi) <- paste(name[i],"=",levels(xi),sep="")
              Strata[[nstrata]] <- xi
            }

          names(Strata) <- paste("S",1:nstrata,sep="")
          Strata <- interaction(as.data.frame(Strata),drop=TRUE)
        }

      offs <- offset <- attr(Terms, "offset")
      
      xpres <- length(asm) && any(asm!=8)
      Y <- model.extract(X, 'response')
      if(!inherits(Y,"Surv"))
        stop("response variable should be a Surv object")
    
      weights <- model.extract(X, 'weights')
      tt <- length(offset)
      offset <- if(tt == 0) rep(0, nrow(Y))
      else if(tt == 1) X[[offset]]
      else
        {
          ff <- X[[offset[1]]]
          for(i in 2:tt)   # for case with multiple offset terms
            ff <- ff + X[[offset[i]]]
          ff
        }
    
      if(model) m <- X
      
      ##No mf if only strata factors
      if(!xpres)
        {
          X <- matrix(nrow=0,ncol=0)
          assign <- NULL
        }
      else
        {
          X <- model.matrix(Terms.ns, X)[,-1,drop=FALSE]
          assign <- attr(X, "assign")
          assign[[1]] <- NULL  # remove intercept position, renumber
        }
      
      nullmod <- FALSE
    }
  else
    {	# model with no right-hand side
      X <- NULL
      Terms <- terms(formula)
      yy <- attr(terms(formula),"variables")[1]
      Y <- eval(yy, data)
      if(!inherits(Y,"Surv"))
        stop("response variable should be a Surv object")
    
      Y <- Y[!is.na(Y)]
      assign <- NULL
      xpres <- FALSE
      nullmod <- TRUE
      nact <- NULL
    }

  ny <- ncol(Y)
  time.units <- attr(Y, "units")
  maxtime <- max(Y[,ny-1])

  rnam <- dimnames(Y)[[1]]
  if(xpres) dimnames(X) <- list(rnam, atr$colnames)

  if(method=="model.matrix") return(X)

  if(!length(time.units)) time.units <- "Day"
  
  if(missing(time.inc))
    {
      time.inc <- switch(time.units,
                         Day=30,
                         Month=1,
                         Year=1,
                         maxtime/10)
      
      if(time.inc >= maxtime | maxtime/time.inc>25)
        time.inc <- max(pretty(c(0,maxtime)))/10
    }

  if(nullmod) f <- NULL
  else
    {
      ytype <- attr(Y, "type")
      if( method=="breslow" || method =="efron")
        {
          if (ytype== 'right')
            fitter <- survival:::coxph.fit
          else if (ytype=='counting')
            fitter <- survival:::agreg.fit
          else
            stop(paste("Cox model doesn't support \"", ytype,
                       "\" survival data", sep=''))
        }
      else if (method=='exact')
        fitter <- survival:::agexact.fit
      else
        stop(paste ("Unknown method", method))

      if (missing(init)) init <- NULL

      f <- fitter(X, Y, strata=Strata, offset=offset,
                  weights=weights, init=init,
                  method=method, rownames=rnam,
                  control=survival:::coxph.control(eps=eps, toler.chol=tol,
                    toler.inf=1, iter.max=iter.max))
    }
  if (is.character(f))
    {
      cat("Failure in cph:\n",f,"\n")
      return(structure(list(fail=TRUE), class="cph"))
    }
  else
    {
      if(length(f$coefficients) && any(is.na(f$coefficients))) {
        vars <- names(f$coefficients)[is.na(f$coefficients)]
        msg <- paste("X matrix deemed to be singular; variable",
                     paste(vars, collapse=" "))
        if(singular.ok) warning(msg)
        else
          {
            cat(msg,"\n")
            return(structure(list(fail=TRUE),class="cph"))
          }
      }
    }
  f$terms <- Terms
  
  if(robust)
    {
      f$naive.var <- f$var
      ## Terry gets a little tricky here, calling resid before adding
      ## na.action method to avoid re-inserting NAs.  Also makes sure
      ## X and Y are there
      if(!length(cluster)) cluster <- FALSE
    
      fit2 <- c(f, list(x=X, y=Y, method=method))
      if(length(stra)) fit2$strata <- Strata
    
      temp <- survival:::residuals.coxph(fit2, type='dfbeta', collapse=cluster)
      f$var <- t(temp) %*% temp
    }
  
  if(length(weights) && any(weights!=1)) f$weights <- weights

  nvar <- length(f$coefficients)

  temp <- factor(Y[,ny], levels=0:1, labels=c("No Event","Event"))
  n.table <- {
    if(!length(Strata)) table(temp,dnn='Status')
    else table(Strata, temp, dnn=c('Stratum','Status'))
  }
  f$n <- n.table
  nnn <- nrow(Y)
  nevent <- sum(Y[,ny])
  if(xpres)
    {
      logtest <- -2 * (f$loglik[1] - f$loglik[2])
      R2.max <- 1 - exp(2*f$loglik[1]/nnn)
      R2 <- (1 - exp(-logtest/nnn))/R2.max
      P <- 1-pchisq(logtest,nvar)
      stats <- c(nnn, nevent, logtest, nvar, P, f$score, 
                 1-pchisq(f$score,nvar), R2)
      names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P", 
                        "Score", "Score P","R2")
    }
  else
    {
      stats <- c(nnn, nevent)
      names(stats) <- c("Obs","Events")
    }

  f$method <- NULL
  if(xpres)
    dimnames(f$var) <- list(atr$colnames, atr$colnames)

  f <- c(f, list(call=call, Design=atr,
                 assign=DesignAssign(atr, 0, atrx$terms),
                 na.action=nact,
                 fail = FALSE, non.slopes = 0, stats = stats, method=method,
                 maxtime = maxtime, time.inc = time.inc,
                 units = time.units, fitFunction=c('cph','coxph')))

  if(xpres)
    {
      f$center <- sum(f$means*f$coefficients)
      f$scale.pred <- c("log Relative Hazard","Hazard Ratio")
      attr(f$linear.predictors,"strata") <- Strata
      names(f$linear.predictors) <- rnam
      if(se.fit)
        {
          XX <- X - rep(f$means, rep.int(nnn, nvar))   # see scale() function
          ##  XX <- sweep(X, 2, f$means)	# center   (slower)
          se.fit <- drop(((XX %*% f$var) * XX) %*% rep(1,ncol(XX)))^.5
          names(se.fit) <- rnam
          f$se.fit <- se.fit  	
        }
    }
  if(model) f$model <- m
  
  if(nstrata > 0)
    {
      attr(X, "strata") <- attr(Y, "strata") <- Strata
      f$strata <- levels(Strata)
    }
  
  if(x) f$x <- X
  if(y) f$y <- Y
  
  if(is.character(surv) || surv)
    {
      Strata <- if(!length(Strata)) rep(1, nnn) else unclass(Strata)
      
      nstr <- max(Strata, na.rm=TRUE)
      srv  <- NULL
      tim  <- NULL
      s.e. <- NULL
      timepts <- seq(0, maxtime, by=time.inc)
      s.sum <- array(double(1),
                     c(length(timepts),nstr,3),
                     list(format(timepts),paste("Stratum",1:nstr),
                          c("Survival","n.risk","std.err")))
      
      attr(Terms,'specials')$strata <- attr(Terms,'specials')$strat
      g <- list(formula=list(n=sum(f$n), coefficients=f$coefficients,
                  linear.predictors=f$linear.predictors,
                  method=f$method, means=f$means, var=f$var,
                  x=X, y=Y, strata=Strata, terms=Terms))
      class(g$formula) <- if(xpres) 'coxph' else 'coxph.null'

      if(!missing(type))      g$type      <- type
      if(!missing(vartype))   g$vartype   <- vartype

      if(length(Strata))
        {
          n.all <- table(Strata) # patch bug in survfit.coxph.null
          storeTemp(n.all)       # survfit.coxph.null can't find here
                                 # (why? namespace?)
        }
      g <- do.call('survfit', g)

      stemp <- if(nstr==1) rep(1, length(g$time)) else rep(1:nstr, g$strata)

      i <- 0
      for(k in 1:nstr)
        {
          j    <- stemp==k
          i    <- i+1
          yy   <- Y[Strata==i, ny-1]
          maxt <- max(yy)
          ##n.risk from surv.fit does not have usual meaning if not Kaplan-Meier
      
          tt <- c(0,  g$time[j])
          su <- c(1,  g$surv[j])
          se <- c(NA, g$std.err[j])
          
          if(maxt > tt[length(tt)])
            {
              tt <- c(tt, maxt)
              su <- c(su, su[length(su)])
              se <- c(se, NA)
            }

          kk <- 0
          for(tp in timepts)
            {
              kk <- kk + 1
              t.choice <- max((1:length(tt))[tt <= tp+1e-6])
              if(tp > max(tt)+1e-6 & su[length(su)] > 0)
                {
                  Su <- NA
                  Se <- NA
                }
              else
                {
                  Su <- su[t.choice]
                  Se <- se[t.choice]
                }
              
              n.risk <- sum(yy >= tp)
              s.sum[kk, i, 1:3] <- c(Su, n.risk, Se)
            }
          
          if(!is.character(surv))
            {
              if(nstr==1)
                {
                  tim  <- tt
                  srv  <- su
                  s.e. <- se
                }
              else
                {
                  tim  <- c(tim,  list(tt))
                  srv  <- c(srv,  list(su))
                  s.e. <- c(s.e., list(se))
                }
            }
        }

      if(is.character(surv)) f$surv.summary <- s.sum
      else
        {
          attr(srv, "type") <- if(missing(type)) method
          else type

          if(nstr>1)
            {
              names(srv) <- names(tim) <- names(s.e.) <- f$strata
            }

          f <- c(f, list(time=tim, surv=srv,
                         std.err=s.e., surv.summary=s.sum))		
        }
    }
  
  class(f) <- c("cph", "rms", "coxph")
  
  f
}

coxphFit <- function(..., method, strata=NULL, rownames=NULL, offset=NULL,
                     init=NULL, toler.chol=1e-9, eps=.0001, iter.max=10,
                     type)
{
  if( method == "breslow" || method == "efron")
    {
      fitter <- if (type == 'right')
        getFromNamespace('coxph.fit', 'survival')
      else getFromNamespace('agreg.fit', 'survival')
    }
  else if (method == 'exact')
    fitter <- getFromNamespace('agexact.fit', 'survival')
  else stop("Unkown method ", method)
  
  if(!existsFunction('coxph.control'))
    coxph.control <- getFromNamespace('coxph.control', 'survival')

  res <- fitter(..., strata=strata, rownames=rownames,
                offset=offset, init=init, method=method,
                control=coxph.control(toler.chol=toler.chol, toler.inf=1,
                  eps=eps, iter.max=iter.max))
  
  if(is.character(res)) return(list(fail=TRUE))
  
  if(iter.max > 1 && res$iter >= iter.max) return(list(fail=TRUE))

  res$fail <- FALSE
  res
}

Survival.cph <- function(object, ...)
{
  if(!length(object$time) || !length(object$surv))
    stop("did not specify surv=T with cph")
  f <- function(times, lp=0, stratum=1, type=c("step","polygon"),
                time, surv)
    {
      type <- match.arg(type)
      if(length(stratum)>1) stop("does not handle vector stratum")
      if(length(times)==0)
        {
          if(length(lp)>1) stop("lp must be of length 1 if times=NULL")
          return(surv[[stratum]]^exp(lp))
        }
      s <- matrix(NA, nrow=length(lp), ncol=length(times),
                  dimnames=list(names(lp), format(times)))
      if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
      if(type=="polygon")
        {
          if(length(lp)>1 && length(times)>1)
            stop('may not have length(lp)>1 & length(times>1) when type="polygon"')
          su <- approx(time, surv, times, ties=mean)$y
          return(su ^ exp(lp))
        }
      for(i in 1:length(times))
        {
          tm <- max((1:length(time))[time <= times[i]+1e-6])
          su <- surv[tm]
          if(times[i] > max(time)+1e-6) su <- NA
          s[,i] <- su^exp(lp)
        }
      drop(s)
    }
  formals(f) <- list(times=NULL, lp=0, stratum=1,
                     type=c("step","polygon"),
                     time=object$time, surv=object$surv)
  f
}

Quantile.cph <- function(object, ...)
{
  if(!length(object$time) || !length(object$surv))
    stop("did not specify surv=T with cph")
  f <- function(q=.5, lp=0, stratum=1, type=c("step","polygon"), time, surv)
    {
      type <- match.arg(type)
      if(length(stratum)>1) stop("does not handle vector stratum")
      if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
      Q <- matrix(NA, nrow=length(lp), ncol=length(q),
                  dimnames=list(names(lp), format(q)))
      for(j in 1:length(lp))
        {
          s <- surv^exp(lp[j])
          if(type=="polygon") Q[j,] <- approx(s, time, q, ties=mean)$y
          else for(i in 1:length(q))
            if(any(s <= q[i])) Q[j,i] <- min(time[s<=q[i]])  #is NA if none
        }
      drop(Q)
    }
  formals(f) <- list(q=.5, lp=0, stratum=1,
                     type=c('step','polygon'),
                     time=object$time, surv=object$surv)
  f
}


Mean.cph <- function(object, method=c("exact","approximate"),
                     type=c("step","polygon"), n=75, tmax, ...)
{
  method <- match.arg(method)
  type   <- match.arg(type)

  if(!length(object$time) || !length(object$surv))
    stop("did not specify surv=T with cph")
  
  if(method=="exact")
    {
      f <- function(lp=0, stratum=1, type=c("step","polygon"),
                    tmax=NULL, time, surv)
        {
          type <- match.arg(type)
          if(length(stratum)>1) stop("does not handle vector stratum")
          if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
          Q <- lp
          if(!length(tmax))
            {
              if(min(surv)>1e-3)
                warning(paste("Computing mean when survival curve only defined down to",
                              format(min(surv)),"\n Mean is only a lower limit"))
              k <- rep(TRUE,length(time))
            }
          else
            {
              if(tmax>max(time)) stop(paste("tmax=",format(tmax),
                                            "> max follow-up time=",
                                            format(max(time))))
              k <- (1:length(time))[time<=tmax]
            }
          for(j in 1:length(lp))
            {
              s <- surv^exp(lp[j])
              Q[j] <- if(type=="step") sum(c(diff(time[k]),0) * s[k]) else 
              trap.rule(time[k], s[k])
            }
          Q
        }
      formals(f) <- alist(lp=0, stratum=1,
                          type=if(!missing(type))type else c("step","polygon"),
                          tmax=tmax,
                          time=object$time, surv=object$surv)
    }
  else
    {
      lp     <- object$linear.predictors
      lp.seq <- if(length(lp)) lp.seq <- seq(min(lp), max(lp), length=n) else 0
  
      time   <- object$time
      surv   <- object$surv
      nstrat <- if(is.list(time)) length(time) else 1
      areas  <- list()

      for(is in 1:nstrat)
        {
          tim <- if(nstrat==1) time else time[[is]]
          srv <- if(nstrat==1) surv else surv[[is]]
          if(!length(tmax))
            {
              if(min(srv)>1e-3)
                warning(paste("Computing mean when survival curve only defined down to",
                              format(min(srv)),
                              "\n Mean is only a lower limit"))
              k <- rep(TRUE,length(tim))
            }
          else
            {
              if(tmax>max(tim)) stop(paste("tmax=",format(tmax),
                                           "> max follow-up time=",
                                           format(max(tim))))
              k <- (1:length(tim))[tim<=tmax]
            }
          ymean <- lp.seq
          for(j in 1:length(lp.seq))
            {
              s <- srv^exp(lp.seq[j])
              ymean[j] <- if(type=="step") sum(c(diff(tim[k]),0) * s[k]) else 
              trap.rule(tim[k], s[k])
            }
          areas[[is]] <- ymean
        }
      if(nstrat>1) names(areas) <- names(time)

      f <- function(lp=0, stratum=1, lp.seq, areas)
        {

          if(length(stratum)>1) stop("does not handle vector stratum")
          area <- areas[[stratum]]
          if(length(lp.seq)==1 && all(lp==lp.seq))
            ymean <- rep(area,length(lp)) else
          ymean <- approx(lp.seq, area, xout=lp, ties=mean)$y
          if(any(is.na(ymean)))
            warning("means requested for linear predictor values outside range of linear\npredictor values in original fit")
          names(ymean) <- names(lp)
          ymean
        }
      formals(f) <- list(lp=0, stratum=1, lp.seq=lp.seq, areas=areas)
    }
  eval(f)
}

predict.cph <- function(object, newdata=NULL,
                        type=c("lp", "x", "data.frame", "terms", "adjto",
                          "adjto.data.frame", "model.frame"),
                        se.fit=FALSE, conf.int=FALSE,
                        conf.type=c('mean','individual'),
                        incl.non.slopes=NULL, non.slopes=NULL, kint=1,
                        na.action=na.keep, expand.na=TRUE,
                        center.terms=TRUE, ...)
  predictrms(object, newdata, type, se.fit, conf.int, conf.type,
             incl.non.slopes, non.slopes, kint,
             na.action, expand.na, center.terms, ...)

#This is Terry Therneau's old print.coxreg with conf.int default to F
#Add Nagelkerke R2 9Jun92
#Remove printing hazard ratios 17Jun92

print.cph.fit <- 
  function(x, table = TRUE, coef = TRUE, conf.int = FALSE, scale = 1,
           digits = NULL, ...)
{
  if(table && !is.null(x$n) && is.matrix(x$n))
    print(x$n)
  if(is.null(digits))
    digits <- 3
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  beta <- x$coef
  se <- sqrt(diag(x$var))
  if(is.null(beta) | is.null(se))
    stop("Input is not valid")
  if(coef)
    {
      tmp <- cbind(beta, se, beta/se, 1 - pchisq((beta/
                                                  se)^2, 1))
      dimnames(tmp) <- list(names(beta), c("coef", 
                                           "se(coef)", "z", "p"))
      cat("\n")
      prmatrix(tmp)
	}
  if(conf.int)
    {
      z <- qnorm((1 + conf.int)/2, 0, 1)
      beta <- beta * scale
      se <- se * scale
      tmp <- cbind(exp(beta), exp( - beta), exp(beta - z * se),
                   exp(beta + z * se))
      dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
			paste("upper .", round(100 * conf.int, 2), sep = "")))
      cat("\n")
      prmatrix(tmp)
	}
  invisible(x)
}

print.cph <- function(x, long=FALSE, digits=3, conf.int=FALSE,
                      table=TRUE,  ...)
{ 

  cat("\n")
  if(x$fail) stop("Model Did Not Converge")


  cat("Cox Proportional Hazards Model\n\n")
  dput(x$call)
  cat("\n")
  if(length(z <- x$na.action)) naprint(z)
  if(length(x$coef))
    {
      stats <- x$stats
      stats[3] <- round(stats[3],2)
      stats[5] <- round(stats[5],4)
      stats[6] <- round(stats[6],2)
      stats[7] <- round(stats[7],4)
      stats[8] <- round(stats[8],3)
      print(format.sep(stats), quote=FALSE)
      cat("\n")
      print.cph.fit(x, digits=digits, conf.int=conf.int, table=table, ...)
      if(long)cat("Centering constant:",format(x$center),"\n")
    }
  else if(table) print(x$n)
  invisible()
}
