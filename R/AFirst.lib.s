.noGenenerics <- TRUE  # faster loading as new methods not used

. <- NA   ## used by Predict etc.

.First.lib <- function(lib, pkg, ...)
{
  library.dynam("rms", pkg, lib)
  require(Hmisc) || stop('Hmisc package not available')
  invisible()
}
