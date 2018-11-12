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
            npar <- 3
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
            if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar,"\n","The third parameter here relates to nu, set it to zero for a regular lognormal distribution")
              parameters$mu <- pstart[1]
              parameters$sigma <- log(pstart[2])
              parameters$nu <- pstart[3]
            }
            mod <- new("lnorm", parameters = parameters, fixed = fixed, x=x,y=y,npar=npar)
            mod
          }
)

setMethod("show","lnorm",
          function(object) {
            cat("Model of type lnorm/lognormal (see ?gamlss for details) \n")
            cat("Parameters: \n")
            cat("mu: ", object@parameters$mu, "\n")
            cat("sigma: ", object@parameters$sigma, "\n")
            cat("nu: ", object@parameters$nu, "\n")
          }
)

setMethod("dens","lnorm",
          function(object,log=FALSE) {
            dLNO(object@y, mu = predict(object), sigma=exp(object@parameters$sigma), nu=object@parameters$nu, log=log)
          }
)

setMethod("getpars","response",
          function(object,which="pars",...) {
            switch(which,
                   "pars"= {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   })
            return(pars)
          }
)

setMethod("setpars","lnorm",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$mu <- values[1]
                     object@parameters$sigma <- values[2]
                     object@parameters$nu <- values[3]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

setMethod("fit","lnorm",
          function(object,w) {
            if(missing(w)) w <- NULL
            y <- object@y
            fit <- gamlss(y~1,weights=w,family=LNO(),
                          control=gamlss.control(n.cyc=100,trace=FALSE),
                          mu.start=object@parameters$mu,
                          sigma.start=exp(object@parameters$sigma),
                          nu.start=object@parameters$nu)
            pars <- c(fit$mu.coefficients,fit$sigma.coefficients,fit$nu.coefficients)
            object <- setpars(object,pars)
            object
          }
)

setMethod("predict","exgaus", 
          function(object) {
            ret <- object@parameters$mu
            return(ret)
          }
)