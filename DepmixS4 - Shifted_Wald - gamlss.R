setClass("SWald", contains="response") 

setGeneric("SWald", function(y, pstart=NULL, fixed=NULL,...) standardGeneric("SWald"))

setMethod("SWald",
          signature(y="ANY"),
          function(y,pstart=NULL,fixed=NULL,...){
            y <- matrix(y, length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 3
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar) #"\n","The third parameter here relates to nu, set it to zero for a regular lognormal distribution")
              parameters$mu <- log(pstart[1])
              parameters$sigma <- log(pstart[2])
              parameters$nu <- log(pstart[3])
            }
            mod <- new("SWald", parameters=parameters, fixed=fixed, x=x,y=y,npar=npar)
            mod
          }
)

setMethod("show","SWald",
          function(object) {
            cat("Model of type Shifted Wald (see ?statmod::dinvgauss for details) \n")
            cat("Parameters: \n")
            cat("mu: ", object@parameters$mu, "\n")
            cat("sigma: ", object@parameters$sigma, "\n")
            cat("nu: ", object@parameters$nu, "\n")
          }
)

setMethod("dens","SWald",
          function(object,log=FALSE) {
            dSwald(object@y,mu=predict(object),sigma=exp(object@parameters$sigma),nu=exp(object@parameters$nu),log=log)  }
)

setMethod("getpars","response",
          function(object,which="pars",...) {
            switch(which,
                   "pars" = {
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

setMethod("setpars","SWald",
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

setMethod("fit","SWald",
          function(object,w) {
            if(missing(w)) w <- NULL
            y <- object@y
            #y <- as.vector(y)
            #fit <- optim(par=c(exp(object@parameters$alpha),exp(object@parameters$gamma),exp(object@parameters$theta)),fn=obj.seq,w=w, y=y)
            #pars <- c(fit$par[1],fit$par[2],fit$par[3])#,fit$nu.coefficients)
            fit <- gamlss(y~1,weights=w,family=SWald(),
                                 control=gamlss.control(c.crit=1e-5,n.cyc=100,trace=FALSE),
                                 mu.start=exp(object@parameters$mu),
                                 sigma.start=exp(object@parameters$sigma),
                                 nu.start=exp(object@parameters$nu))
            pars <- c(fit$mu.coefficients,fit$sigma.coefficients,fit$nu.coefficients)
            object <- setpars(object,pars)
            object
          }
)

setMethod("predict","SWald", 
          function(object) {
            ret <- object@parameters$mu
            return(ret)
          }
)
