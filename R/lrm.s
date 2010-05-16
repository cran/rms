lrm <- function(formula,data,subset,na.action=na.delete,
				method="lrm.fit",model=FALSE, x=FALSE, y=FALSE, 
				linear.predictors=TRUE, se.fit=FALSE, 
				penalty=0, penalty.matrix, tol=1e-7, strata.penalty=0,
                var.penalty=c('simple','sandwich'),
                weights, normwt=FALSE, ...)
{
  call <- match.call()
  var.penalty <- match.arg(var.penalty)
  m <- match.call(expand=FALSE)
  mc <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action
  if(.R.) m$drop.unused.levels <- TRUE
  
  m[[1]] <- as.name("model.frame")
  nact <- NULL
  if(missing(data)) data <- NULL

  tform <- terms(formula, specials='strat', data=data)
  offs <- attr(tform, "offset")
  nstrata <- 1
  if(!missing(data) || (
						length(atl <- attr(tform,"term.labels")) && 
						any(atl!=".")))	{ ##X's present

    dul <- .Options$drop.unused.levels
    if(!length(dul) || dul)
      {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }

    X <- Design(eval.parent(m))
    atrx <- attributes(X)
    nact <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms <- atrx$terms
    attr(Terms, "formula") <- formula
    atr <- atrx$Design

    Y <- model.extract(X, 'response')
    weights <- wt <- model.extract(X, 'weights')
    if(length(weights))
      warning('currently weights are ignored in model validation and bootstrapping lrm fits')
    ofpres <- length(offs) > 0
    if(ofpres) offs <- X[[offs]]  else offs <- 0
    if(model) m <- X
    stra <- attr(tform,'specials')$strat
    Strata <- NULL
    Terms.ns <- Terms
    if(length(stra))
      {
        temp <- untangle.specials(Terms.ns, 'strat', 1)
        Terms.ns <- Terms.ns[-temp$terms]
        attr(Terms,   "factors") <- pmin(attr(Terms,"factors"),1)
        attr(Terms.ns,"factors") <- pmin(attr(Terms.ns,"factors"),1)
        Strata <- X[[stra]]
        nstrata <- length(levels(Strata))
      }
    X <- model.matrix(Terms.ns, X)
    X <- X[,-1,drop=FALSE]
    dimnames(X)[[2]] <- atr$colnames
    xpres <- length(X)>0

    p <- length(atr$colnames)
    n <- length(Y)

    penpres <- !(missing(penalty) && missing(penalty.matrix))
    if(penpres && missing(var.penalty))
      warning('default for var.penalty has changed to "simple"')

    if(!penpres) penalty.matrix <- matrix(0,ncol=p,nrow=p) else { 
      if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr, X) else
      if(nrow(penalty.matrix)!=p || ncol(penalty.matrix)!=p) stop(
             paste("penalty.matrix does not have",p,"rows and columns"))
      psetup <- Penalty.setup(atr, penalty)
      penalty <- psetup$penalty
      multiplier <- psetup$multiplier
      if(length(multiplier)==1)
        penalty.matrix <- multiplier*penalty.matrix
      else
        {
          a <- diag(sqrt(multiplier))
          penalty.matrix <- a %*% penalty.matrix %*% a
        }
    }
  }
  else
    {
      X <- eval.parent(m)
      ofpres <- FALSE
      if(length(offs))
        {
          ofpres <- TRUE
          offs <- X[[offs]]
        }
      Y <- model.extract(X, 'response')
      Y <- Y[!is.na(Y)]
      Terms <- X <- NULL
      xpres <- FALSE
      penpres <- FALSE
      penalty.matrix <- NULL
    }  ##Model: y~. without data= -> no predictors
  
  if(method=="model.matrix") return(X)

  if(nstrata > 1)
    {
      f <- if(ofpres)
        lrm.fit.strat(X,Y,Strata,offset=offs,
                      penalty.matrix=penalty.matrix,
                      strata.penalty=strata.penalty,
                      tol=tol,
                      weights=weights,normwt=normwt, ...)
      else
        lrm.fit.strat(X,Y,Strata,penalty.matrix=penalty.matrix,
                      strata.penalty=strata.penalty,tol=tol,
                      weights=weights, normwt=normwt, ...)
    }
  else
    {
      if(existsFunction(method))
        {
          fitter <- getFunction(method)
          if(ofpres) f <- fitter(X,Y,offset=offs,
                                 penalty.matrix=penalty.matrix,tol=tol,
                                 weights=weights, normwt=normwt, ...)
          else f <- fitter(X,Y,
                           penalty.matrix=penalty.matrix,tol=tol,
                           weights=weights, normwt=normwt, ...)
        }
	  else stop(paste("unimplemented method:", method))
    }
  
  if(f$fail) stop("Unable to fit model using ", dQuote(method))
  
  f$call <- NULL
  if(model) f$model <- m
  if(x)
    {
      f$x <- X
      f$strata <- Strata
	}
  if(y) f$y <- Y
  nrp <- f$non.slopes
  if(penpres)
    {
      f$penalty <- penalty
      if(nstrata==1) {
        ## Get improved covariance matrix
        v <- f$var
        if(var.penalty=='sandwich') f$var.from.info.matrix <- v
        f.nopenalty <- if(ofpres) fitter(X,Y,offset=offs,initial=f$coef, maxit=1, tol=tol) else
        fitter(X,Y,initial=f$coef, maxit=1, tol=tol)
        ##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)
        info.matrix.unpenalized <- f.nopenalty$info.matrix
        dag <- diag(info.matrix.unpenalized %*% v)
        f$effective.df.diagonal <- dag
        f$var <- if(var.penalty=='simple')v else
        v %*% info.matrix.unpenalized %*% v
        df <- sum(dag[-(1:nrp)])
        lr <- f.nopenalty$stats["Model L.R."]
        pval <- 1-pchisq(lr, df)
        f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)  
      }
    }
  ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()
  
  if(xpres)
    {
      if(linear.predictors) names(f$linear.predictors) <- names(Y)
      else
        f$linear.predictors <- NULL

      if(se.fit)
        {
          if(nstrata > 1) stop('se.fit=T not implemented for strat')
          xint <- matrix(0, nrow=length(Y), ncol=f$non.slopes)
          xint[,1] <- 1
          X <- cbind(xint, X)
          se <- drop((((X %*% f$var) * X) %*% rep(1, ncol(X)))^.5)
          names(se) <- names(Y)
          f$se.fit <- se
        }
    }
  f <- c(f, list(call=call, Design=if(xpres)atr,
                 scale.pred=c("log odds","Odds Ratio"),
				 terms=Terms, assign=ass, na.action=nact, fail=FALSE,
                 nstrata=nstrata, fitFunction=c('lrm','glm')))
  
  oldClass(f) <- c("lrm","rms","glm")
  f
}

print.lrm <- function(x, digits=4, strata.coefs=FALSE, ...)
{
  sg <- function(x,d)
    {
      oldopt <- options(digits=d)
      on.exit(options(oldopt))
      format(x)
    }
  rn <- function(x,d) format(round(as.single(x),d))

  cat("\n")
  if(x$fail)
    {
      cat("Model Did Not Converge\n")
      return()
    }

  cat("Logistic Regression Model\n\n")
  dput(x$call)
  cat("\n\nFrequencies of Responses\n")
  print(x$freq)
  if(length(x$sumwty))
    {
      cat('\n\nSum of Weights by Response Category\n')
      print(x$sumwty)
    }
  cat("\n")
  if(!is.null(x$nmiss))
    {  #for backward compatibility
      cat("Frequencies of Missing Values Due to Each Variable\n")
      print(x$nmiss)
      cat("\n")
    }
  else if(length(x$na.action)) naprint(x$na.action)

  ns <- x$non.slopes
  nstrata <- x$nstrata
  if(!length(nstrata)) nstrata <- 1

  pm <- x$penalty.matrix
  if(length(pm))
    {
      psc <- if(length(pm)==1) sqrt(pm)
      else
        sqrt(diag(pm))
      penalty.scale <- c(rep(0,ns),psc)
      cof <- matrix(x$coef[-(1:ns)], ncol=1)
      cat("Penalty factors:\n\n")
      print(as.data.frame(x$penalty, row.names=''))
      cat("\nFinal penalty on -2 log L:",
          rn(t(cof) %*% pm %*% cof,2),"\n\n")
    }

  ## ?ok to have uncommented next 3 lines?
  est.exp <- 1:ns
  if(length(x$est)) est.exp <-
    c(est.exp, ns + x$est[x$est+ns <= length(x$coefficients)])
  vv <- diag(x$var)
  cof <- x$coef
  if(strata.coefs)
    {
      cof <- c(cof, x$strata.coef)
      vv  <- c(vv,  x$Varcov(x,which='strata.var.diag'))
      if(length(pm)) penalty.scale <- c(penalty.scale,rep(NA,x$nstrat-1))
    }
  score.there <- nstrata==1 && (length(x$est) < length(x$coef)-ns)
  stats <- x$stats
  stats[2] <- signif(stats[2],1)
  stats[3] <- round(stats[3],2)
  stats[4] <- round(stats[4],2)
  stats[5] <- round(stats[5],4)
  stats[6] <- round(stats[6],3)
  stats[7] <- round(stats[7],3)
  if(nstrata==1)
    {
      for(j in 8:13) stats[j] <- round(stats[j], 3)
      if(length(stats) > 13)
        {
          stats[14] <- round(stats[14], 3)
          if(length(x$weights)) stats[15] <- round(stats[15], 3)
        }
    }
  else
    stats <- c(stats, Strata=x$nstrat)

  nst <- length(stats)
  cstats <- character(nst)
  names(cstats) <- names(stats)
  for(i in 1:nst) cstats[i] <- format(stats[i])
  print(cstats, quote=FALSE)

  cat("\n")

  z <- cof/sqrt(vv)
  stats <- cbind(sg(cof,digits), sg(sqrt(vv),digits), 
                 rn(cof/sqrt(vv),2))
  stats <- cbind(stats, rn(1-pchisq(z^2,1),4))
  dimnames(stats) <- list(names(cof),
                          c("Coef","S.E.","Wald Z","P"))
  if(length(pm))
    stats <- cbind(stats, "Penalty Scale"=sg(penalty.scale,digits))
  print(stats, quote=FALSE)
  cat("\n")


  if(score.there)
    {
      q <- (1:length(cof))[-est.exp]
      if(length(q)==1) vv <- x$var[q,q] else vv <- diag(x$var[q,q])
      z <- x$u[q]/sqrt(vv)
      stats <- cbind(rn(z,2), rn(1-pchisq(z^2,1),4))
      dimnames(stats) <- list(names(cof[q]),c("Score Z","P"))
      print(stats,quote=FALSE)
      cat("\n")
    }
  invisible()

}


