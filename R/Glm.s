Glm <- 
  function(formula, family = gaussian, data = list(), weights = NULL,
           subset = NULL, na.action = na.delete, start = NULL, offset = NULL,
           control = glm.control(...), model = TRUE, method = "glm.fit",
           x = FALSE, y = TRUE, contrasts = NULL, ...)
{
  call <- match.call()
  if (is.character(family)) family <- get(family)
  if (is.function(family))  family <- family()
  if (!length(family$family))
    {
      print(family)
      stop("`family' not recognized")
    }
  mt <- terms(formula, data = data)
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
  mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$... <- NULL
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  mf[[1]] <- as.name("model.frame")

  dul <- .Options$drop.unused.levels
  if(!length(dul) || dul)
    {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }

  mf <- Design(eval(mf, parent.frame()))
  desatr <- attr(mf,'Design')
  attr(mf,'Design') <- NULL
  nact <- attr(mf,'na.action')
    
  switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1,
         stop(paste("invalid `method':", method)))
  xvars <- as.character(attr(mt, "variables"))[-1]
  if ((yvar <- attr(mt, "response")) > 0)
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0)
    {
      xlev <- lapply(mf[xvars], levels)
      xlev[!sapply(xlev, is.null)]
    }
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
  Y <- model.response(mf, "numeric")
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if (length(weights) && any(weights < 0))
    stop("Negative wts not allowed")
  if (length(offset) && length(offset) != NROW(Y))
    stop(paste("Number of offsets is", length(offset), ", should equal",
               NROW(Y), "(number of observations)"))
  fit <- (if (is.empty.model(mt)) glm.fit.null
  else glm.fit)(x = X, y = Y, weights = weights, start = start,
                offset = offset, family = family, control = control,
                intercept = attr(mt, "intercept") > 0)
  if (length(offset) && attr(mt, "intercept") > 0)
    {
      fit$null.deviance <- if (is.empty.model(mt))
        fit$deviance
      else glm.fit(x = X[, "(Intercept)", drop = FALSE], y = Y,
                   weights = weights, start = start, offset = offset,
                   family = family, control = control,
                   intercept = TRUE)$deviance
    }
  if (model) fit$model <- mf
  if (x)  fit$x <- X
  if (!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
                     data = data, offset = offset, control = control,
                     method = method,
                     contrasts = attr(X, "contrasts"), xlevels = xlev,
                     Design=desatr, na.action=nact, fitFunction='Glm',
                     assign=DesignAssign(desatr,1,mt),
                     g=GiniMd(fit$linear.predictors)))
  class(fit) <- c('Glm', 'rms',
                  if (is.empty.model(mt)) "glm.null", "glm",
                  "lm")
  fit
}

print.Glm <- function(x, digits=4, ...)
{
  cat('General Linear Model\n\n')
  dput(x$call); cat('\n\n')
  if(length(z <- x$na.action)) naprint(z)
  cof <- coef(x)
  lr <- x$null.deviance - x$deviance
  names(cof) <- ifelse(names(cof)=='(Intercept)','Intercept',names(cof))
  dof <- x$rank - (names(cof)[1]=='Intercept')
  pval <- 1 - pchisq(lr, dof)
  print(c('Model L.R.'=format(lr,digits=2), 'd.f.'=format(dof),
          'P'=format(pval,digits=4), 'g'=round(x$g,3)), quote=FALSE)
  cat('\n')
  se <- sqrt(diag(vcov(x)))
  z <- cof/se
  p <- 1 - pchisq(z^2, 1)
  w <- cbind(format(cof, digits=digits),
             format(se,  digits=digits),
             format(z,   digits=2),
             format(p,   digits=4))
  dimnames(w) <- list(names(cof), c('Coef','S.E.','Wald Z','P'))
  print(w, quote=FALSE)
  invisible()
}

summary.Glm <- function(...) summary.rms(...)

vcov.Glm <- function(object, ...) stats:::vcov.glm(object, ...)

## Varcov.Glm <- function(object, ...) vcov.
#  Varcov.glm(object, regcoef.only, ...)
# Varcov.glm <- function(object, ...)
#{
#  if(length(object$var))
#    return(object$var)  ## for Glm
#  
#  s <- summary.glm(object)
#  s$cov.unscaled * s$dispersion
#}

predict.Glm <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms", "adjto",
             "adjto.data.frame", "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...)
  {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               incl.non.slopes, non.slopes, kint,
               na.action, expand.na, center.terms, ...)
  }

latex.Glm <- function(...) latexrms(...)
