latex.lrm <-
  function(object, title, 
           file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='', ...)
{
  f <- object
  
  if(missing(which) & !inline)
    {
      Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")
      lev <- names(f$freq)
      nrp <- f$non.slopes
      
      w <- '\\['
      
      j <- if(lev[2]=="TRUE") "" else paste("=",lev[2],sep="")
      if(nrp==1) w <- paste(w,"{\\rm Prob}\\{",Y, j,
           "\\} = \\frac{1}{1+\\exp(-X\\beta)}", sep="")

      else
        w <- paste(w,"{\\rm Prob}\\{", Y, 
                   "\\geq j\\} = \\frac{1}{1+\\exp(-\\alpha_{j}-X\\beta)}",
                   sep="")

      w <- paste(w, ", {\\rm \\ \\ where} \\\\ \\]", sep="")

      if(length(caption)) w <- c(paste('\\begin{center} \\bf',caption,
                                       '\\end{center}'), w)
      
      if(nrp>1)
        {
          w <- c(w,"\\begin{eqnarray*}")
          cof <- format(f$coef[1:nrp])
          for(i in 1:nrp)
            w <- c(w, paste("\\hat{\\alpha}_{\\rm ",
                            lev[i+1],"} &=&",cof[i],"\\\\",sep=""))
          w <- c(w,"\\end{eqnarray*}",sep="")
        }
							}
  else w <- NULL
  if(missing(which) | missing(varnames)) at <- f$Design

  if(missing(which)) which <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  cat(w, file=file, append=append, sep=if(length(w))"\n" else "")
  latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
           columns=columns, 
           before=before, after=after, prefix="X\\hat{\\beta}",
           inline=inline, pretrans=pretrans, digits=digits, size=size)
}


latex.orm <-
  function(object, title, 
           file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           intercepts=nrp < 10, ...)
{
  f <- object
  
  if(missing(which) & !inline)
    {
      Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")
      lev <- names(f$freq)
      nrp <- f$non.slopes

      z <- '\\alpha_{y} + X\\beta'
      zm <- '- \\alpha_{y} - X\\beta'
      dist <-
        switch(f$family,
               logistic = paste('\\frac{1}{1+\\exp(', zm, ')}', sep=''),
               probit   = paste('\\Phi(', z, ')', sep=''),
               cauchit  = paste('\\frac{1}{\\pi}\\tan^{-1}(', z,
                 ') + \\frac{1}{2}', sep=''),
               loglog   = paste('\\exp(-\\exp(', zm, '))', sep=''),
               cloglog  = paste('1 - \\exp(-\\exp(', z, ')', sep=''))
                     
      w <- '\\['
      
      w <- paste(w, "{\\rm Prob}\\{", Y, 
                   "\\geq y | X\\} = ", dist, sep='')

      w <- paste(w, ", {\\rm \\ \\ where} \\\\ \\]", sep="")

      if(length(caption)) w <- c(paste('\\begin{center} \\bf',caption,
                                       '\\end{center}'), w)
      
      if(intercepts) {
        nl <- as.numeric(lev)
        if(!any(is.na(nl))) lev <- format(nl, digits=digits)
          w <- c(w,"\\begin{eqnarray*}")
          cof <- format(f$coef[1:nrp], digits=digits)
          for(i in 1:nrp)
            w <- c(w, paste("\\hat{\\alpha}_{\\rm ",
                            lev[i+1],"} &=&",cof[i],"\\\\",sep=""))
          w <- c(w, "\\end{eqnarray*}", sep="")
      }
    }
  else w <- NULL
  if(missing(which) | missing(varnames)) at <- f$Design

  if(missing(which)) which <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  cat(w, file=file, append=append, sep=if(length(w))"\n" else "")
  latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
           columns=columns, 
           before=before, after=after, prefix="X\\hat{\\beta}",
           inline=inline, pretrans=pretrans, digits=digits, size=size)
}
