##  spatstat.core/R/First.R

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.core"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatCoreVersion", vs)
  packageStartupMessage(paste("spatstat.core", vs))
  return(invisible(NULL))
}

  
