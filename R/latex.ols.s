latex.ols <-
  function(object, title,
           file='',
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           ...)
{
  f <- object

  md <- prType() %in% c('html', 'md', 'markdown')
  
  w <- if(length(caption)) {
         if(md) paste('<div align=center><strong>', caption,
                      '</strong></div>', sep='')
         else
           paste('\\begin{center} \\bf',
                 caption,'\\end{center}')
         }
  
  if(missing(which) & !inline)
    {
      Y <- paste("\\mathrm{",as.character(attr(f$terms,"formula")[2]), "}",
                 sep="")
      
      w <- c(w, paste("$$\\mathrm{E}(", Y,
                      ") = X\\beta,~~\\mathrm{where}$$", sep=""))
    }
  at <- f$Design
  
  if(missing(which)) which <- 1:length(at$name)
  
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  if(file != '') cat(w, file=file, sep=if(length(w)) "\n" else "",
                     append=append)
  ltx <- latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
                  columns=columns, 
                  before=before, after=after, prefix="X\\hat{\\beta}",
                  inline=inline, 
                  pretrans=pretrans, digits=digits, size=size)
  if(inline) return(ltx)
  z <- c(w, ltx)
  if(file == '' && prType() != 'plain') return(rendHTML(z, html=FALSE))
  cat(z, file=file, append=append, sep='\n')
  invisible()
}
