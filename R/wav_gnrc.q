################################################
## S+WMTSA generic functionality
##
##  Functions:
##
##    crystal.names
##    eda.plot
##    reconstruct
##    stack.plot
##
################################################

###
# crystal.names
##

"crystal.names" <- function(x, ...)
  UseMethod("crystal.names")

###
# eda.plot
##

"eda.plot" <- function(x, ...)
  UseMethod("eda.plot")

###
# reconstruct
###

"reconstruct" <- function(x,...)
  UseMethod("reconstruct")

###
# stack.plot
##

"stack.plot" <- function(x, ...)
  UseMethod("stack.plot")