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
              parameters$alpha <- pstart[1]
              parameters$gamma <- log(pstart[2])
              parameters$theta <- log(pstart[3])
            }
            mod <- new("SWald", parameters=parameters, fixed=fixed, x=x,y=y,npar=npar)
            mod
          }
)

setMethod("show","SWald",
          function(object) {
            cat("Model of type Shifted Wald (see ?statmod::dinvgauss for details) \n")
            cat("Parameters: \n")
            cat("alpha: ", object@parameters$alpha, "\n")
            cat("gamma: ", object@parameters$gamma, "\n")
            cat("theta: ", object@parameters$theta, "\n")
          }
)

setMethod("dens","SWald",
          function(object,log=FALSE) {
            dSwald(object@y,mu=predict(object),sigma=exp(object@parameters$gamma),nu=exp(object@parameters$theta),log=log)
          }
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
                     object@parameters$alpha <- values[1]
                     object@parameters$gamma <- values[2]
                     object@parameters$theta <- values[3]
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
            fit <- optim(par=c(object@parameters$alpha,exp(object@parameters$gamma),exp(object@parameters$theta)),fn=obj.seq,w=w, y=y)
            
            ######### ps optim
            #lb <- rep(0,3)
            #ub <- rep(20,3)
            #fit <- psoptim(par=c(exp(object@parameters$alpha),exp(object@parameters$gamma),exp(object@parameters$theta)),obj.seq,lower=lb,upper=ub,y=y,w=w,control=list(trace=T,abstol=1e-8,hybrid="off",hybrid.control=list(pgtol=1e-8)))
            pars <- c(fit$par[1],fit$par[2],fit$par[3])#,fit$nu.coefficients)
            object <- setpars(object,pars)
            object
          }
)

setMethod("predict","SWald", 
          function(object) {
            ret <- object@parameters$alpha
            return(ret)
          }
)
