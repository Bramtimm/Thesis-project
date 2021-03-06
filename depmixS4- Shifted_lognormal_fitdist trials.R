dSLOGNO <- function(x, meanlog=meanlog, sdlog = sdlog, shift=shift, log = FALSE) 
{
  if (any(sdlog <= 0)) 
   stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(x < 0)){ 
    stop(paste("x must be >=0", "\n", ""))}
  fy <- ifelse(shift>=x,0,dshifted_lnorm(x = x, meanlog = meanlog, sdlog = sdlog, shift=shift , log = log))
  fy
}

pSLOGNO <- function(q, meanlog = meanlog, sdlog = sdlog, shift=shift, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sdlog <= 0)) 
    stop(paste("sdlog must be greater than 0 ", "\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be >=0", "\n", ""))
  if(any(shift<0)){
    stop(paste("Gamma must be smaller than x and higher than zero"))
  }
  cdf <- pshifted_lnorm(q,meanlog = meanlog, sdlog = sdlog, shift=shift,lower.tail=lower.tail,
                        log.p = log.p)
  cdf
}

qSLOGNO <- function(p,meanlog = meanlog, sdlog = sdlog, shift=shift, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sdlog <= 0)) 
    stop(paste("sdlog must be greater than 0 ", "\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  q <- qshifted_lnorm(p = p, meanlog = meanlog, sdlog = sdlog, shift=shift, lower.tail = lower.tail, 
                      log.p = log.p)
  q
}

rSLOGNO <- function(n,meanlog = meanlog, sdlog = sdlog, shift=shift) 
{
  if (any(sdlog <= 0)) 
    stop(paste("sdlog must be greater than 0 ", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  r <- rshifted_lnorm(n = n, meanlog = meanlog, sdlog = sdlog, shift = shift)
}

setClass("s_lnorm_fitdist", contains="response") 

setGeneric("s_lnorm_fitdist", function(y, pstart=NULL, fixed=NULL,...) standardGeneric("s_lnorm_fitdist"))

setMethod("s_lnorm_fitdist",
          signature(y="ANY"),
          function(y,pstart=NULL,fixed=NULL,...){
            y <- matrix(y, length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 3
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar) #"\n","The third parameter here relates to nu, set it to zero for a regular lognormal distribution")
              parameters$meanlog <- pstart[1]
              parameters$sdlog <- log(pstart[2])
              parameters$shift <- log(pstart[3])
            }
            mod <- new("s_lnorm_fitdist", parameters=parameters, fixed=fixed, x=x,y=y,npar=npar)
            mod
          }
)

setMethod("show","s_lnorm_fitdist",
          function(object) {
            cat("Model of type shifted lnorm/lognormal \n")
            cat("Parameters: \n")
            cat("meanlog: ", object@parameters$meanlog, "\n")
            cat("sdlog: ", object@parameters$sdlog, "\n")
            cat("shift: ", object@parameters$shift, "\n")
          }
)

setMethod("dens","s_lnorm_fitdist",
          function(object,log=FALSE) {
            dshifted_lnorm(object@y, meanlog = predict(object), sdlog=exp(object@parameters$sdlog),shift=exp(object@parameters$shift), log=log)
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

setMethod("setpars","s_lnorm_fitdist",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$meanlog <- values[1]
                     object@parameters$sdlog <- values[2]
                     object@parameters$shift <- values[3]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

setMethod("fit","s_lnorm_fitdist",
          function(object,w) {
            if(missing(w)) w <- NULL
            y <- object@y
            
            ###################### with fitdistr (does often not work)
            #fit <- fitdistr(y, dens=dSLOGNO,lower=c(-Inf,1e-10,0),upper=c(Inf,Inf,max(y)),
             #               method="L-BFGS-B",control=list(maxit=500,pgtol=1e-20,trace=TRUE),
              #                                 start=list(meanlog=object@parameters$meanlog,
               #                                           sdlog=object@parameters$sdlog,
                #                                          shift=object@parameters$shift)) 
            #pars <- c(fit$estimate[1],fit$estimate[2],fit$estimate[3])
            
            #######################with psoptim
            lb <- c(-100,0.01,0.01)
            ub <- c(100,100,100)
            fit <- psoptim(c(object@parameters$meanlog,exp(object@parameters$sdlog),exp(object@parameters$shift)),obj,lower=lb,upper=ub,y=y,w=w,control=list(trace=F,abstol=1e-8,hybrid="improved",hybrid.control=list(pgtol=1e-8)))
            pars <- c(fit$par[1],fit$par[2],fit$par[3])
            #######################with DEoptim
            #lb <- c(-100,0,-1)
            #ub <- c(100,100,max(y))
            #fit <- DEoptim(obj, lower=lb, upper=ub, control=DEoptim.control(
            #       itermax=500, trace=F, NP=1000), rt=y,w=w)
            #pars <- fit$optim$bestmem
            #pars <- unname(pars)
            object <- setpars(object,pars)
            object
          }
)

setMethod("predict","s_lnorm_fitdist", 
          function(object) {
            ret <- object@parameters$meanlog
            return(ret)
          }
)
