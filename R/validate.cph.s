validate.cph <- function(fit, method="boot",
                         B=40, bw=FALSE, rule="aic", type="residual",
                         sls=.05, aics=0, pr=FALSE,
                         dxy=FALSE, u, tol=1e-9, ...)
{
  atr <- fit$Design

  need.surv <- dxy & any(atr$assume.code==8)
  
  if(need.surv & missing(u))
    stop("Presence of strata -> survival estimates needed for dxy; u omitted")
  
  modtype <- fit$method
  
  discrim <- function(x, y, strata, fit, iter, evalfit=FALSE, dxy=FALSE,
                      need.surv=FALSE, u, modtype, pr=FALSE, ...)
    {
      n <- nrow(y)
      if(!length(x) || length(unique(x))==1)
        {
          Dxy <- 0
          slope <- 1
          D <- 0
          U <- 0
          R2 <- 0
        }
      else
        {
          x <- as.matrix(x)
          dimnames(x) <- list(as.character(1:nrow(x)),as.character(1:ncol(x)))
          if(evalfit)
            {	#Fit was for training sample
              lr  <- -2 * (fit$loglik[1] - fit$loglik[2])
              ll0 <- -2 * fit$loglik[1]
              slope <- 1
              D <- (lr - 1)/ll0
              U <- -2/ll0
              R2.max <- 1 - exp(-ll0/n)
              R2 <- (1 - exp(-lr/n))/R2.max
              g  <- GiniMd(x)
            }
          else
            {
              type <- attr(y, "type")
              storage.mode(x) <- "double"
              f <- coxphFit(x=x, y=y, strata=strata, iter.max=10, eps=.0001,
                            method=modtype, type=type)
              if(f$fail)
                stop('fit failure in discrim,coxphFit')

              ##x is x*beta from training sample
              lr <- -2 * (f$loglik[1]-f$loglik[2])
              ll0 <- -2 * f$loglik[1]
              slope <- f$coef[1]
              D <- (lr - 1)/ll0
              R2.max <- 1 - exp(-ll0/n)
              R2 <- (1 - exp(-lr/n))/R2.max
              f.frozen <- coxphFit(x=x, y=y, strata=strata,
                                   iter.max=0, eps=.0001,
                                   method=modtype, init=1, type=type)
              if(f.frozen$fail) stop('fit failure in discrim for f.frozen')
              U <- -2 * (f.frozen$loglik[2] - f$loglik[2]) / ll0
              g <- GiniMd(slope*x)
            }
        }
      
      Q <- D - U
      z   <- c(R2,  slope,    D,  U,   Q,   g)
      nam <- c("R2","Slope", "D", "U", "Q", "g")
      if(dxy)
        {
          if(need.surv)
            {
              attr(x, "strata") <- strata
              x <- survest(fit, linear.predictors=x, times=u,
                           conf.int=FALSE)$surv
            }
          Dxy <- rcorr.cens(x, y)["Dxy"]
          z <- c(Dxy, z)
          nam <- c("Dxy", nam)
        }
      names(z) <- nam
      z
    }
  
  cox.fit <- function(x, y, strata, u, need.surv=FALSE, modtype, tol=1e-9,
                      ...)
    {
      if(!length(x))
        return(list(fail=FALSE,coefficients=numeric(0)))
      
      if(!need.surv)
        u <- 0
      
      ##	coxph(x,y,e,pr=F,surv=need.surv)
      if(!need.surv)
        {
          type <- attr(y, 'type')
          storage.mode(x) <- "double"
          x <- as.matrix(x)
          dimnames(x) <- list(as.character(1:nrow(x)),as.character(1:ncol(x)))
          
          f <- coxphFit(x=x, y=y, strata=strata, iter.max=10, eps=.0001,
                        method=modtype, toler.chol=tol, type=type)
          
          if(f$fail) return(f)
          
          if(any(is.na(f$coef)))
            {
              cat('Singularity in coxph.fit. Coefficients:\n'); print(f$coef)
              return(list(fail=TRUE))
            }
          return(f)
        }

      x <- x      #Don't want lazy evaluation of complex expression
      f <- if(length(strata))
        cph(y ~ x + strat(strata), surv=TRUE, method=modtype)
      else
        cph(y ~ x, surv=TRUE, method=modtype)
      f$non.slopes <- f$assume.code <- f$assign <- f$name <- f$assume <- NULL
      ##Don't fool fastbw called from predab.resample
      f
    }
  
  predab.resample(fit, method=method, fit=cox.fit, measure=discrim,
                  pr=pr, B=B, bw=bw, rule=rule, type=type, sls=sls,
                  aics=aics, dxy=dxy, u=u, need.surv=need.surv,
                  modtype=modtype,tol=tol, ...)
}
