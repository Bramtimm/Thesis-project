############## Modeling with LNR ###############################################

# Author: Bram Timmers
# use: adding a lognormal model to depmixS4
# Dependencies: dependencies.R <- need gamlss and depmixS4 at minimum
################################################################################

source('dependencies.R')

# define a response class which only contains the standard slots, no additional slots
setClass("lnorm", contains="response") 

# note that we name our distribution now as lnorm.

# define a generic for the method defining the response class
setGeneric("lnorm", function(y, pstart=NULL, fixed=NULL,...) standardGeneric("lnorm"))

# define the method that creates the response class
setMethod("lnorm",
          signature(y="ANY"),
          function(y,pstart=NULL,fixed=NULL,...){
            y <- matrix(y, length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 2
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
            if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
              parameters$mu <- pstart[1]
              parameters$sigma <- pstart[2]
            }
            mod <- new("lnorm", parameters = parameters, fixed = fixed, x=x,y=y,npar=npar)
            mod
          }
)