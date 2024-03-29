latex.pphsm <-
  function(object, title,
           file='',
           append=FALSE, which=NULL, varnames, 
           columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           ...)
{
  md <- prType() %in% c('html', 'md', 'markdown')
  
  whichThere <- length(which)
  w <- if(length(caption)) {
         if(md) paste('<div align=center><strong>', caption,
                      '</strong></div>', sep='')
         else
           paste('\\begin{center} \\bf',caption,'\\end{center}')
       }

  sc <- object$scale
  at <- object$Design

  if(!whichThere & !inline)
    {
      dist <- paste("\\exp(-t^{", format(1 / sc, digits=digits),
                    "} \\exp(X\\hat{\\beta}))")
      w <- c(w,paste("$$\\Pr(T\\geq t) = ", dist,
                     "~\\mathrm{where}~~$$",sep=""))
    }				
  if(!whichThere) which <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code != 9]
  if(file != '') cat(w, file=file, sep=if(length(w))"\n" else "",
                     append=append)

  ltx <- latexrms(object, file=file, append=TRUE, which=which,
                  varnames=varnames, 
                  columns=columns, 
                  before=before, after=after,
                  prefix=if(!whichThere)"X\\hat{\\beta}" else NULL, 
                  inline=inline,pretrans=pretrans, digits=digits,
                  size=size)
  if(inline) return(ltx)
  z <- c(w, ltx)
  if(file == '' && prType() != 'plain') return(rendHTML(z, html=FALSE))
  cat(z, file=file, append=append, sep='\n')
  invisible()
}


