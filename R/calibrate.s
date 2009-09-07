calibrate <- function(fit, ...) UseMethod("calibrate")

print.calibrate <- function(x, ...)
{
  at <- attributes(x)
  dput(at$call); cat("\n")
  print.default(x)
  invisible()
}

plot.calibrate <- function(x, xlab, ylab, subtitles=TRUE,
                           conf.int=TRUE, cex.subtitles=.75, ...)
{
  at    <- attributes(x)
  u     <- at$u
  units <- at$units
  pred  <- x[,"mean.predicted"]
  KM    <- x[,"KM"]
  obs.corrected <- x[,"KM.corrected"]
  se <- x[,"std.err"]
  if(missing(xlab)) xlab <- paste("Predicted ",format(u),units,"Survival")
  if(missing(ylab)) ylab <- paste("Fraction Surviving ",format(u)," ",units,
                                  "s",sep="")

  ##Remember that groupkm stored the se of the log(-log) survival

  if(conf.int)
    errbar(pred, KM,
           ifelse(KM==0 | KM==1, NA, exp(-exp(logb(-logb(KM))-1.96*se))),
           ifelse(KM==0 | KM==1, NA, exp(-exp(logb(-logb(KM))+1.96*se))),
           xlab=xlab, ylab=ylab, type="b", axes=TRUE, ...)
  else plot(pred, KM, xlab=xlab, ylab=ylab, type="b", ...)

  if(subtitles)
    {
      title(sub=paste("n=",at$n," d=",at$d," p=",at$p,
              ", ",at$m," subjects per group",sep=""),adj=0,
            cex.sub=cex.subtitles)
      title(sub=paste("X - resampling optimism added, B=",at$B,
              "\nBased on ",at$what,sep=""),
            adj=1,cex.sub=cex.subtitles)
    }
  
  abline(0,1,lty=2)	#ideal
  points(pred, obs.corrected, pch=4)
  
  invisible()
}
