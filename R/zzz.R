########################################################
## S+WMTSA library attachment/detachment functions
##
##  .First.lib
##  .Last.lib
##
########################################################

if (!is.R()) {
	.First.lib <- function(libname, pkgname)
	{
    ifultools.pos <- attached.where("ifultools", nomatch=0)
    if (!ifultools.pos) {
      cat("Loading required library: ifultools.\n")
      library("ifultools", first=TRUE)
      assign("ifultools.previously.wmtsa.loaded", TRUE, frame=0)
    }
    invisible(NULL)
	}
	
  .Last.lib <- function(library , section , .data , where)
  {
    # ifultools library
    ifultools.id <- "ifultools.previously.wmtsa.loaded"
    if (exists(ifultools.id, frame=0)){
      was.loaded <- get(ifultools.id, frame=0)
      if (was.loaded && attached.where("ifultools", nomatch=0)){
         cat("Detaching library: ifultools.\n")
         detach("ifultools")
      }
      remove(ifultools.id, frame=0)
    }
    invisible(NULL)
  }
}
