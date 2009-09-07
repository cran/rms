gendata <- function(fit, ...) UseMethod("gendata")

gendata.default <- function(fit, ...) gendata.rms(obj, ...)

gendata.rms <- function(fit, nobs, viewvals=FALSE,
	editor=.Options$editor, ..., factors)
{
  at <- fit$Design

  nam <- at$name[at$assume!="interaction"]

  if(!missing(nobs) && !is.logical(nobs))
    {
      df <- predictrms(fit, type="adjto.data.frame")
      df[1:nobs,] <- df
      cat("Edit the list of variables you would like to vary.\nVariables not listed will be set to reference values.\n")
      if(editor=="xedit")
        cat("To delete an individual variable, type Cntl-k\nTo delete blocks of variables, highlight the block by holding down the left\nmouse button, then type Cntl-w.\n")
      nam.sub <- edit(nam, editor=editor)
      if(!all(nam.sub %in% nam)) stop("misspelled a variable name")
      df.sub <- as.data.frame(df[,nam.sub])  #df[,] was returning list (?)
      cat("Edit the predictor settings to use.\n")
      if(viewvals && 
         length(val <- Getlim(at, allow.null=TRUE,
                              need.all=FALSE)$values[nam.sub]))
        {
          cat("A window is being opened to list the valid values of discrete variables.\n")
          sink(tf <- tempfile())
          print.datadist(list(values=val))
          sink()
          file.show(tf)
        }
##      if(existsFunction('Edit.data'))
##        {
##          stop('use of S-PLUS 4.x GUI not yet implemented for gendata')
##          assign('.df.sub.', df.sub, where=1)
##          Edit.data(.df.sub., '.df.sub.')
##          df.sub <- get('.df.sub.', where=1)
##          remove('.df.sub.', where=1)
##        }
##      else
##      if(existsFunction('data.ed')) df.sub <- data.ed(df.sub)
##      else
      if(existsFunction('data.entry')) df.sub <- data.entry(df.sub)
      df[nam.sub] <- df.sub
      return(structure(df, names.subset=nam.sub))
    }

  factors <- if(missing(factors)) list(...) else factors
  fnam <- names(factors)
  nf <- length(factors)
  if(nf==0) return(predictrms(fit, type="adjto.data.frame"))
  which <- charmatch(fnam, nam, 0)
  if(any(which==0)) stop(paste("factor(s) not in design:",
           paste(names(factors)[which==0],collapse=" ")))
  settings <- if(nf<length(nam)) predictrms(fit, type="adjto.data.frame")
  else
    list()
  settings <- oldUnclass(settings)
  if(nf>0) for(i in 1:nf) settings[[fnam[i]]] <- factors[[i]]
  if(nf==0) return(as.data.frame(settings))
  expand.grid(settings)
}
