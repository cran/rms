Rq <- function (formula, tau = 0.5, data, subset, weights, na.action=na.delete, 
                method = "br", model = FALSE, contrasts = NULL,
                se='nid', hs=TRUE, x=FALSE, y=FALSE, ...) 
{
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  mf[[1]] <- as.name("model.frame")
  mf <- Design(eval.parent(mf))

  at <- attributes(mf)
  desatr <- at$Design
  attr(mf,'Design') <- NULL
  if (method == "model.frame") return(mf)
  mt <- at$terms
  weights <- model.weights(mf)
  Y <- model.response(mf)
  X <- model.matrix(mt, mf, contrasts)
  if (length(desatr$colnames)) colnames(X) <- c("Intercept", desatr$colnames)

  eps <- .Machine$double.eps^(2/3)
  Rho <- function(u, tau) u * (tau - (u < 0))
  if (length(tau) > 1)
    stop('does not allow more than one quantile to be estimated simultaneously')
  require(quantreg)
  fit <- if (length(weights)) 
    rq.wfit(X, Y, tau = tau, weights, method, ...)
  else rq.fit(X, Y, tau = tau, method, ...)
  rownames(fit$residuals) <- rownames(dimnames(X)[[1]])
  rho <- sum(Rho(fit$residuals, tau))
  
  stats <- c(n=length(fit$residuals),
             p=length(fit$coefficients),
             g=GiniMd(fit$fitted.values))
  
  fit <- c(fit,
           list(
                na.action = at$na.action,
                formula   = formula,
                terms     = mt,
                xlevels   = .getXlevels(mt, mf),
                call      = call,
                tau       = tau,
                method    = method,
                weights   = weights,
                residuals = drop(fit$residuals),
                rho       = rho,
                fitted.values = drop(fit$fitted.values),
                model     = mf,
                Design    = desatr,
                assign    = DesignAssign(desatr, 1, mt),
                fitFunction=c("Rq", "rq"),
                stats     = stats))
  attr(fit, "na.message") <- attr(m, "na.message")
  
  s <- summary.rq(fit, cov=TRUE, se=se, hs=hs)
  k <- s$coefficients
  nam <- names(fit$coefficients)
  rownames(k) <- nam
  fit$summary <- k
  cov <- s$cov
  dimnames(cov) <- list(nam, nam)
  fit$var <- cov
  
  ## Remove the following since summary.rq has done its job
  if(!model) fit$model <- NULL
  if(!x) fit$x <- NULL
  if(!y) fit$y <- NULL
  class(fit) <- c('rms',
                  if (method == "lasso") "lassorq"
                  else if (method == "scad") "scadrq",
                  "Rq", "rq")
  fit
}

## Thanks to Duncan Murdoch for the formals alist substitute technique
RqFit <- function(fit, wallow=TRUE, passdots=FALSE)
  {
    w <- fit$weights
    if(length(w))
      {
        if(!wallow) stop('weights not implemented')
        g <- if(passdots) function(x, y, weights, tau, method, ...)
          rq.wfit(x, y, tau = tau, weights=weights, method=method, ...)
        else function(x, y, weights, tau, method, ...)
          rq.wfit(x, y, tau = tau, weights=weights, method=method)
        formals(g) <- eval(substitute(
                       alist(x=,y=, weights=,tau=deftau,method=defmethod,...=),
                       list(deftau=fit$tau, defmethod=fit$method)))
      }
    else
      {
        g <- if(passdots) function(x, y, tau, method, ...)
          rq.fit(x, y, tau = tau, method=method, ...)
        else
          function(x, y, tau, method, ...)
            rq.fit(x, y, tau = tau, method=method)
        formals(g) <-
          eval(substitute(alist(x=,y=, tau=deftau, method=defmethod,...=),
                          list(deftau=fit$tau, defmethod=fit$method)))
      }
    g
  }

print.Rq <- function(x, digits=4, ...)
  {
    cat('Quantile Regression\t\ttau:', format(x$tau), '\n\n')
    dput(x$call); cat('\n')
    if(length(z <- x$na.action)) naprint(z)
    s <- x$stats
    n <- s['n']; p <- s['p']; g <- s['g']
    cat('n=', n, '  p=', p, '  residual d.f.=', n-p,
        '  g=', round(g,3), '\n\n')
    print(x$summary)
    if (length(attr(x, "na.message"))) cat(attr(x, "na.message"), "\n")
  }

latex.Rq <-
  function(object,
           file = paste(first.word(deparse(substitute(object))),
             ".tex", sep = ""), append=FALSE,
           which, varnames, columns=65, inline=FALSE, caption=NULL,
           ...)
  {
    f   <- object
    tau <- f$tau
    at  <- f$Design
    
    w <- if (length(caption)) 
        paste("\\begin{center} \\bf", caption, "\\end{center}")
    if (missing(which) & !inline)
      {
        Y <- paste("{\\rm ", as.character(formula(f))[2], 
                   "}", sep = "")
        w <- c(w, paste("\\[", Y, "_{", tau, "} = X\\beta, {\\rm \\ \\ where} \\\\ \\]", 
                        sep = ""))
      }
    if(missing(which)) which <- 1:length(at$name)
    if(missing(varnames)) varnames <- at$name
    
    cat(w, file = file, sep = if (length(w)) "\n"  else "", append = append)
    latexrms(f, file=file, append=TRUE, which=which, inline=inline,
             varnames=varnames, columns=columns, caption, ...)
  }
