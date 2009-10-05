bplot <-
  function(x, xlab, ylab, zlab,
           adj.subtitle=TRUE, cex.adj, 
           perim,  method=c("image","persp","contour"),
           zlim=range(yhat, na.rm=TRUE), nlevels=10, ...)
{
  method  <- match.arg(method)
  fit     <- x
  info    <- attr(x, 'info')
  varying <- info$varying
  if(length(varying) != 2) stop('two variables should be varying')
  nx      <- varying[1]
  ny      <- varying[2]
  
  yhat    <- x$yhat
  y       <- x[[ny]]
  x       <- x[[nx]]
  xu      <- sort(unique(x))
  yu      <- sort(unique(y))

  at      <- info$Design
  label   <- at$label
  units   <- at$units
  npersp  <- method != 'persp'

  if(missing(xlab))
    xlab  <- labelPlotmath(label[nx], units[nx], plotmath=npersp)
  if(missing(ylab))
    ylab  <- labelPlotmath(label[ny], units[ny], plotmath=npersp)
  if(missing(zlab))
    zlab  <- if(npersp) info$ylabPlotmath else info$ylab
  
  adjust  <- info$adjust
  if(!adj.subtitle) adjust <- NULL
  
  cex <- par('cex')
  if(missing(cex.adj)) cex.adj <- .75*cex

  if(!missing(perim))
    {
      Ylo <- approx(perim[,1], perim[,2], x, ties=mean)$y
      Yhi <- approx(perim[,1], perim[,3], x, ties=mean)$y
      Ylo[is.na(Ylo)] <-  1e30
      Yhi[is.na(Yhi)] <- -1e30
      yhat[y < Ylo] <- NA
      yhat[y > Yhi] <- NA
    }
      
  zmat <- matrix(pmin(zlim[2], pmax(zlim[1], yhat)),
                 nrow=length(xu),
                 ncol=length(yu), byrow=TRUE)

  switch(method,
         contour = contour(xu, yu, zmat,
           xlab=xlab, ylab=ylab, nlevels=nlevels, ...),
         persp   = persp(xu, yu, zmat, zlim=zlim, xlab=xlab, ylab=ylab,
           zlab=zlab, box=TRUE, ...),
         image   = image(xu, yu, zmat, xlab=xlab, ylab=ylab, ...))

  if(length(adjust)) title(sub=paste('Adjusted to:', adjust),
                           cex.sub=cex.adj, adj=0)
}

iLegend <- function(object, x, y, size = c(1, 1), 
                    horizontal = TRUE, nint = 50, fun.=NULL, at=NULL, 
                    zlab, zlim, par.=NULL, ...)
{
  ## Note: fun. is used instead of fun because subplot has arg fun
  if(missing(x)) 
    if(missing(size))
      {
        cat("Using function \"locator(2)\" to place opposite corners of legend\n")
        x   <- locator(2)
        x$x <- sort(x$x)
        x$y <- sort(x$y)
      }
    else
      {
        cat("Using function \"locator(1)\" to place upper left corner of legend\n")
        x <- locator(1)
      }

  if(missing(zlab)) zlab <- attr(object, 'info')$ylabPlotmath
    
  if(!missing(y)) x <- list(x=x,y=y)
  z <- object$yhat

  if(missing(zlim)) zlim <- range(z, na.rm=TRUE)
  
  irgz  <- seq(zlim[1], zlim[2], length = nint)
  lirgz <- length(irgz)

  dotlist <- list(...)
  
  if(horizontal)
    {
      f <- function()
        {
          if(length(par.))
            {
              opar <- do.call('par', par.)
              on.exit(par(opar))
            }
          ##axis() does not respect mgp

          image(x=irgz, y=1:lirgz, z=matrix(irgz, lirgz, lirgz),
                xlab=zlab, ylab='',
                yaxt="n", xaxt=if(!length(fun.))"s" else "n", ...)
        
          if(length(fun.))
            mgp.axis(1,
                     if(!length(at)) pretty(irgz)
                     else at,
                     labels=format(fun.(if(!length(at)) pretty(irgz)
                     else at)))
          
          title(xlab=zlab)
        }
      subplot(x = x$x, y = x$y, size = size, fun = f, hadj=0, vadj=1)
    }
  else
    {
      f <- function()
        {
          if(length(par.))
            {
              opar <- do.call('par', par.)
              on.exit(par(opar))
            }
          image(x = 1:lirgz, y = irgz,
                z = matrix(irgz, lirgz, lirgz,byrow=TRUE),
                xlab='', ylab=zlab,
                xaxt = "n", yaxt=if(!length(fun.))"s" else "n", ...)
          
          if(length(fun.))
            mgp.axis(2, if(!length(at)) pretty(irgz) else at,
                     labels=format(fun.(if(!length(at)) pretty(irgz) else at)))
          
          title(ylab=zlab)
        }
  subplot(x = x$x, y = x$y, size = size, fun = f, hadj=0, vadj=1)
}
  
invisible(x)
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
      if(n > (m - n + 1)) c(NA, NA)
      else c(y[n], y[m-n+1])
    }

  r <- unlist(tapply(y, x, g, n=n))
  i <- seq(1, length(r), by=2)
  rlo <- r[i]
  rhi <- r[-i]
  s <- !is.na(rlo + rhi)
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
            dimnames=list(NULL,
              c("x","ymin","ymax")), class='perimeter')
}

lines.perimeter <- function(x, ...)
{
  lines(x[,'x'], x[,'ymin'],...)
  lines(x[,'x'], x[,'ymax'],...)
  invisible()
}

