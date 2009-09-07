calibrate.psm <- function(fit,method="boot",u,m=150,cuts,B=40,
		bw=FALSE,rule="aic",
		type="residual",sls=.05,aics=0,
		pr=FALSE,what="observed-predicted",tol=1e-12, maxiter=15, 
	    rel.tolerance=1e-5, ...)
{
  call <- match.call()
  if(!length(fit$y)) stop("fit did not store y")
  oldopt <- options(digits=3)
  on.exit(options(oldopt))
  unit <- fit$units
  if(unit=="") unit <- "Day"
  ny <- dim(fit$y)
  nevents <- sum(fit$y[,ny[2]])

  ##Note: fit$y already has been transformed by the link function by psm

  if(missing(cuts))
    {
      g <- max(1,floor(ny[1]/m))
      survival <- survest.psm(fit,times=u,conf.int=FALSE)$surv
      cuts <- quantile(c(0,1,survival), seq(0,1,length=g+1),na.rm=TRUE)
    }
  
  dist <- fit$dist
  inverse <- survreg.distributions[[dist]]$itrans
  if(!length(inverse)) inverse <- function(x) x
  parms <- fit$parms
  
  distance <- function(x,y,fit,iter,u,fit.orig,what="observed",inverse,
                       orig.cuts, ...)
    {
      ##Assumes y is matrix with 1st col=time, 2nd=event indicator
      if(sum(y[,2])<5) return(NA)
      oldClass(fit) <- 'psm'   # for survest.psm which uses Survival.psm
      fit$dist <- fit.orig$dist
      psurv <- survest.psm(fit, linear.predictors=x,
                           times=u, conf.int=FALSE)$surv
      ##Assumes x really= x * beta
      pred.obs <- 
        groupkm(psurv,Surv(inverse(y[,1]),y[,2]),u=u,cuts=orig.cuts)
      if(what=="observed") dist <- pred.obs[,"KM"]	else
      dist <- pred.obs[,"KM"] - pred.obs[,"x"]
      if(iter==0) storeTemp(pred.obs)
      dist
    }

  b <- min(10,B)
  overall.reps <- max(1,round(B/b))
  ## Bug in S prevents>10 loops in predab.resample
  if(pr)
    cat("\nAveraging ", overall.reps," repetitions of B=",b,"\n\n")
  rel <- 0
  opt <- 0
  nrel <- 0
  B <- 0
  
  for(i in 1:overall.reps)
    {
      reliability <- predab.resample(fit, method=method,
                                     fit=survreg.fit2,measure=distance,
                                     pr=pr, B=b, bw=bw, rule=rule, type=type,  
                                     u=u, m=m, what=what, 
                                     dist=dist, inverse=inverse, parms=parms,
                                     fixed=fixed, family=family,
                                     sls=sls, aics=aics, strata=FALSE,
                                     tol=tol, orig.cuts=cuts, maxiter=maxiter,
                                     rel.tolerance=rel.tolerance, ...)
      n <- reliability[,"n"]
      rel <- rel + n * reliability[,"index.corrected"]
      opt <- opt + n * reliability[,"optimism"]
      nrel <- nrel + n
      B <- B + max(n)	
      if(pr) print(reliability)
    }

  mean.corrected <- rel/nrel
  mean.opt <- opt/nrel
  rel <- cbind(mean.optimism=mean.opt,mean.corrected=mean.corrected,n=nrel)
  if(pr)
    {
      cat("\nMean over ",overall.reps," overall replications\n\n")
      print(rel)
    }
  
  pred <- pred.obs[,"x"]
  KM <- pred.obs[,"KM"]
  se <- pred.obs[,"std.err"]
  obs.corrected <- KM - mean.opt
  
  structure(cbind(reliability[,c("index.orig","training","test"),drop=FALSE],
                  rel,mean.predicted=pred,KM=KM,
                  KM.corrected=obs.corrected,std.err=se), 
            class="calibrate", u=u, units=unit, n=ny[1], d=nevents, 
            p=length(fit$coef)-1, m=m, B=B, what=what, call=call)
}
