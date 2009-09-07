Gls <-
  function (model, data = sys.frame(sys.parent()), correlation = NULL, 
    weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail, 
    control = list(), verbose = FALSE, B=0, dupCluster=FALSE,
            pr=FALSE, opmeth=c('optimize','optim'))
            
{
    require(nlme)
    glsEstimate <- nlme:::glsEstimate
    
    Call <- match.call()
    opmeth <- match.arg(opmeth)
    controlvals <- glsControl()
    if (!missing(control))
      {
        if (!is.null(control$nlmStepMax) && control$nlmStepMax < 0) {
          warning("Negative control$nlmStepMax - using default value")
          control$nlmStepMax <- NULL
        }
        controlvals[names(control)] <- control
      }
    if (!inherits(model, "formula") || length(model) != 3)
        stop("\nModel must be a formula of the form \"resp ~ pred\"")

    method <- match.arg(method)
    REML <- method == "REML"
    if (length(correlation))
      groups <- getGroupsFormula(correlation)
    else groups <- NULL
    glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
    mfArgs <- list(formula = asOneFormula(formula(glsSt), model, 
                     groups), data = data, na.action = na.action)
    if (!missing(subset)) 
      mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]

    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call("model.frame", mfArgs)
    rn <- origOrder <- row.names(dataMod)  ## rn FEH 6apr03
    if (length(groups))
      {
        groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), 
                               sep = "|")))
        grps <- getGroups(dataMod, groups,
                          level = length(getGroupsFormula(groups, 
                            asList = TRUE)))
        ord <- order(grps)
        grps <- grps[ord]
        dataMod <- dataMod[ord, , drop = FALSE]
        rn <- rn[ord]
        revOrder <- match(origOrder, rn)
      }
    else grps <- NULL
    
    X <- model.frame(model, dataMod)
    dul <- .Options$drop.unused.levels
    if(!length(dul) || dul)
      {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    X <- Design(X)
    atrx <- attributes(X)
    desatr <- atrx$Design
    mt <- atrx$terms
    attr(X,'Design') <- NULL
    
    contr <- lapply(X, function(el) if (inherits(el, "factor")) 
                    contrasts(el))
    contr <- contr[!unlist(lapply(contr, is.null))]
    X <- model.matrix(model, X)
    dimnames(X)[[2]] <- cn <- c('Intercept',desatr$colnames)
    y <- eval(model[[2]], dataMod)
    N <- nrow(X)
    p <- ncol(X)
    parAssign <- attr(X, "assign")
    fTerms <- terms(as.formula(model))
    namTerms <- attr(fTerms, "term.labels")
    if (attr(fTerms, "intercept") > 0)
      namTerms <- c("Intercept", namTerms)
    namTerms <- factor(parAssign, labels = namTerms)
    parAssign <- split(order(parAssign), namTerms)
    ## Start FEH 4apr03
    if(B > 0)
      {
        bootcoef <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL,cn))
        bootcorr <- numeric(B)
        Nboot    <- integer(B)
        if(length(grps))
          {
            obsno    <- split(1:N,grps)
            levg     <- levels(grps)
            ng       <- length(levg)
            if(!length(levg)) stop('program logic error')
          }
        else
          {
            obsno <- 1:N
            levg  <- NULL
            ng    <- N
          }
      }
    for(j in 0:B)
      {
        if(j == 0) s <- 1:N else
        {
          if(ng == N) s <- sample(1:N, N, replace=TRUE)
          else
            {
              grps.sampled <- sample(levg, ng, replace=TRUE)
              s <- unlist(obsno[grps.sampled])
              dataMods <- dataMod[s,]
              if(!dupCluster)
                {
                  grp.freqs <- table(grps)
                  newgrps <- factor(rep(paste('C',1:ng,sep=''),
                                        table(grps)[grps.sampled]))
                  dataMods$id <- newgrps
                }
            }
          Nboot[j] <- Nb <- length(s)
          if(pr) cat(j,'\r')
        }
        attr(glsSt, "conLin") <-
          if(j==0)
            list(Xy = array(c(X, y), c(N, p + 1),
                   list(rn, c(cn, deparse(model[[2]])))), 
                 dims = list(N = N, p = p, REML = as.integer(REML)),
                 logLik = 0)
          else
            list(Xy = array(c(X[s,,drop=FALSE], y[s]), c(Nb, p + 1),
                   list(rn[s], c(cn, deparse(model[[2]])))), 
                 dims = list(N = Nb, p = p, REML = as.integer(REML)),
                 logLik = 0)
        ## FEH colnames(X) -> cn, ncol(X) -> p, j>0 case above
        glsEstControl <- controlvals[c("singular.ok", "qrTol")]
        glsSt <- Initialize(glsSt, if(j==0) dataMod else dataMods,
                            glsEstControl)
        parMap <- attr(glsSt, "pmap")
        numIter <- numIter0 <- 0
        repeat
          {
            oldPars <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
            if (length(coef(glsSt)))
              {
                co <- c(coef(glsSt))  ## FEH
                if(opmeth == 'optimize' && co > 1) opmeth <- 'optim'
                best <- if(opmeth=='optimize')
                  optimize(function(z)
                           -logLik(glsSt,z), lower=-12, upper=12)$minimum else
                optim(fn = function(glsPars)
                      -logLik(glsSt, glsPars),
                      par = co,  ## FEH
                      method = "BFGS", 
                      control = list(trace = controlvals$msVerbose, 
                        reltol = if (numIter == 0)
                        controlvals$msTol else 100 * 
                        .Machine$double.eps,
                        maxit = controlvals$msMaxIter))$par
                coef(glsSt) <- best   ## FEH
              }
            attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
            if (!needUpdate(glsSt)) 
              break
            numIter <- numIter + 1
            glsSt <- update(glsSt, if(j==0) dataMod else dataMods)  ## FEH
            aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
            conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
            aConv <- c(beta = max(conv[1:p]))
            conv <- conv[-(1:p)]
            for (i in names(glsSt))
              {
                if (any(parMap[, i]))
                  {
                    aConv <- c(aConv, max(conv[parMap[, i]]))
                    names(aConv)[length(aConv)] <- i
                  }
              }
            if (verbose)
              {
                cat("\nIteration:", numIter)
                ## cat("\nObjective:", format(aNlm$value), "\n")
                ## ERROR: aNlm doesn't exist.  Need to fix.
                print(glsSt)
                cat("\nConvergence:\n")
                print(aConv)
              }
            if (max(aConv) <= controlvals$tolerance) break
            if (numIter > controlvals$maxIter)
              stop("Maximum number of iterations reached without convergence.")
      }
      if(j > 0)
        {
          bootcoef[j,] <- attr(glsSt, "glsFit")[["beta"]]
          bootcorr[j]  <- coef(glsSt$corStruct, unconstrained=FALSE)
        }
        if(j==0) glsSt0 <- glsSt  ## FEH 4apr03
      }  ## end bootstrap reps
    if(pr && B > 0) cat('\n')
    glsSt <- glsSt0   ## FEH
        
    glsFit <- attr(glsSt, "glsFit")
    namBeta <- names(glsFit$beta)
    attr(parAssign, "varBetaFact") <- varBeta <- glsFit$sigma * 
      glsFit$varBeta * sqrt((N - REML * p)/(N - p))
    varBeta <- crossprod(varBeta)
    dimnames(varBeta) <- list(namBeta, namBeta)
    Fitted <- fitted(glsSt)
    if (length(grps))
      {
        grps <- grps[revOrder]
        Fitted <- Fitted[revOrder]
        Resid <- y[revOrder] - Fitted
        attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt)[revOrder])
      }
    else
      {
        Resid <- y - Fitted
        attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt))
      }
    if (controlvals$apVar && FALSE) ## FEH 3apr03
      apVar <-
        nlme:::glsApVar(glsSt, glsFit$sigma,
                        .relStep = controlvals[[".relStep"]],  
                        minAbsPar = controlvals[["minAbsParApVar"]],
                        natural = controlvals[["natural"]])
    else
      apVar <- "Approximate variance-covariance matrix not available"
    dims <- attr(glsSt, "conLin")[["dims"]]
    dims[["p"]] <- p
    attr(glsSt, "conLin") <- NULL
    attr(glsSt, "glsFit") <- NULL
    estOut <- list(modelStruct = glsSt, dims = dims, contrasts = contr, 
                   coefficients = glsFit[["beta"]], varBeta = varBeta,
                   sigma = glsFit$sigma, 
                   apVar = apVar, logLik = glsFit$logLik,
                   numIter = if(needUpdate(glsSt)) numIter else numIter0,  
                   groups = grps, call = Call, method = method,
                   fitted = Fitted,  
                   residuals = Resid, parAssign = parAssign,
                   Design=desatr,assign=DesignAssign(desatr,1,mt),
                   formula=model, opmeth=opmeth,
                   B=B, boot.Coef=if(B > 0) bootcoef,
                   boot.Corr=if(B > 0) bootcorr,
                   Nboot=if(B > 0) Nboot,
                   var=if(B > 0) var(bootcoef))
    
    ## Last 2 lines FEH 29mar03
    if (inherits(data, "groupedData"))
      {
        attr(estOut, "units") <- attr(data, "units")
        attr(estOut, "labels") <- attr(data, "labels")
      }
    attr(estOut, "namBetaFull") <- colnames(X)
    estOut$fitFunction <- 'Gls'
    class(estOut) <- c('Gls','rms','gls')
    estOut
  }


print.Gls <- function(x, digits=4, ...)
{
  ## Following taken from print.gls with changes marked FEH

  summary.gls <- nlme:::summary.gls

  dd <- x$dims
  mCall <- x$call
  if (inherits(x, "gnls"))
    cat("Generalized nonlinear least squares fit\n")
  else
    {
      cat("Generalized least squares fit by ")
      cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    }
  cat("  Model:", deparse(mCall$model), "\n")
  cat("  Data:", deparse(mCall$data), "\n")
  if (length(mCall$subset))
    {
      cat("  Subset:", deparse(asOneSidedFormula(mCall$subset)[[2]]), 
          "\n")
    }
  if (inherits(x, "gnls"))
    cat("  Log-likelihood: ", format(x$logLik), "\n", sep = "")
  else
    {
      cat("  Log-", ifelse(x$method == "REML", "restricted-", 
                           ""),
          "likelihood: ", format(x$logLik), "\n", sep = "")
  }
  cat('\n')
  if(any(names(x)=='var') && length(x$var))
    {
      cat('Using bootstrap variance estimates\n\n')
      se <- sqrt(diag(x$var))
      beta <- coef(x)
      zTable <- cbind(Coef=format(beta,digits=digits),
                      'S.E'=format(se, digits=digits),
                      Z   =format(beta/se, digits=2),
                      'Pr(>|Z|)'=format.pval(2*pnorm(-abs(beta/se)),digits=4))
      print(zTable, quote=FALSE)
    }
  else
    print(summary.gls(x)$tTable)
  
  cat("\n")
  if (length(x$modelStruct) > 0)
    print(summary(x$modelStruct))

  cat("Degrees of freedom:", dd[["N"]], "total;", dd[["N"]] - 
      dd[["p"]], "residual\n")
  cat("Residual standard error:", format(x$sigma), "\n")
  
  cat('Clusters:',length(unique(x$groups)),'\n')
  if(x$B > 0)
    {
      cat('Bootstrap repetitions:',x$B,'\n')
      tn <- table(x$Nboot)
      if(length(tn) > 1)
        {
          cat('Table of Sample Sizes used in Bootstraps\n')
          print(tn)
        }
      else
        cat('Bootstraps were all balanced with respect to clusters\n')
      dr <- diag(x$varBeta)/diag(x$var)
      cat('Ratio of Original Variances to Bootstrap Variances\n')
      print(round(dr,2))
      cat('Bootstrap Nonparametric 0.95 Confidence Limits for Correlation Parameter\n')
      r <- round(quantile(x$boot.Corr, c(.025,.975)),3)
      names(r) <- c('Lower','Upper')
      print(r)
    }
  invisible()
}

vcov.Gls <- function(object, ...)
  if(any(names(object)=='var') && length(object$var))
  object$var else object$varBeta

predict.Gls <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictrms(object, newdata, type, se.fit, conf.int, conf.type,
             incl.non.slopes, non.slopes, kint,
             na.action, expand.na, center.terms, ...)


latex.Gls <- function(...) latexrms(...)
