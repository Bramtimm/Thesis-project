# GAMLSS object for WALD

dSwald<-function(x, mu=mu, sigma=sigma,nu=nu, log=FALSE){
  t=x
  x=t
  fy <-seqmodels::dinvgauss(t=t,kappa=mu,xi=sigma,tau=nu,sigma=1,ln=log)
  if(is.vector(x)==TRUE){ as.vector(fy)
    }else if(is.matrix(x)==TRUE){
      fy <- matrix(fy,nrow=nrow(x),ncol=ncol(x))
    }else if(is.array(x)==TRUE){
      fy <- array(fy,dim=dim(x))
    }
  fy
  }


pSwald <- function(q, mu=mu, sigma=sigma,nu=nu, lower.tail = TRUE, log.p = FALSE){
  t=q
  q=t
  cdf <- seqmodels::pinvgauss(t=t, kappa=mu, xi=sigma, tau = nu, sigma = as.numeric(c(1)),
                              ln = log.p, lower_tail = lower.tail)
  if(is.vector(q)==TRUE){ as.vector(cdf)
  }else if(is.matrix(q)==TRUE){
    cdf <- matrix(cdf,nrow=nrow(q),ncol=ncol(q))
    }else if(is.array(x)==TRUE){
      cdf <- array(fy,dim=dim(q))
    }
    cdf
  }

qSwald <- function(p, mu=mu,sigma=sigma,nu=nu,lower.tail=TRUE,log.p=FALSE){
  q <- seqmodels::qinvgauss(p=p, kappa=mu, xi=sigma, tau = nu, sigma = as.numeric(c(1)),
                         bounds = 3, em_stop = 20, err = 1e-08)
  if(is.vector(p)==TRUE){ as.vector(q)
  }else if(is.matrix(p)==TRUE){
    q <- matrix(q,nrow=nrow(p),ncol=ncol(p))
  }else if(is.array(p)==TRUE){
    q <- array(fy,dim=dim(p))
  }
  q
}

rSwald <- function(n, mu=mu, sigma=sigma, nu=nu){
  r<-seqmodels::rinvgauss(n, kappa=mu, xi=sigma, tau = nu, sigma = as.numeric(c(1)))
  if(is.vector(n)==TRUE){ as.vector(r)
  }else if(is.matrix(n)==TRUE){
    r <- matrix(r,nrow=nrow(n),ncol=ncol(n))
  }else if(is.array(n)==TRUE){
    r <- array(fy,dim=dim(n))
  }
  r
}
  
  
Swald <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link", "Shifted Wald", substitute(mu.link), 
                      c("1/mu^2", "inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Shifted Wald", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "Shifted Wald", substitute(nu.link),
                      c("inverse", "log", "identity", "own"))

structure(list(family = c("Swald", "Shifted wald"),
               parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
               nopar = 3,
               type = "Continuous",
               mu.link = as.character(substitute(mu.link)),
               sigma.link = as.character(substitute(sigma.link)),
               nu.link = as.character(substitute(nu.link)),
               mu.linkfun = mstats$linkfun, 
               sigma.linkfun = dstats$linkfun,
               nu.linkfun = vstats$linkfun,
               mu.linkinv = mstats$linkinv, 
               sigma.linkinv = dstats$linkinv,
               nu.linkinv = vstats$linkinv,
               mu.dr = mstats$mu.eta, 
               sigma.dr = dstats$mu.eta,
               nu.dr = vstats$mu.eta,
               #step1 done
               dldm = function(y,mu,sigma,nu){
                 dldm <- sigma+1/mu-sigma*(1/(y-nu))
                 dldm
               },
               #step2 done
               d2ldm2 = function(y,mu,nu){
                 d2ldm2 <- -1/mu^2-1/(y-nu)
                 d2ldm2
               },
               #step3 done
               dldd = function(y,mu,sigma,nu){
                 dldd <- mu-sigma*(y-nu)
                 dldd
               },
               #step4 done
               d2ldd2 = function(y,nu){
                 d2ldd2 <- nu-y
                 d2ldd2
               },
               #step5 done
               dldv = function(y,mu,sigma,nu){
                 dldv <- (3/2)*(1/(y-nu))-(mu^2/2)*(1/(y-nu)^2)+(sigma^2/2)
                 d1dv
               },
               #step6 done
               d2ldv2 = function(y,mu,nu){
                 d2ldv2 <- ((3/2)-mu^2*(1/(y-nu)))*(1/(y-nu)^2)
                 d2ldv2
               },
               #step7
               d2ldmdd = function(y){rep(1,length(y))},
               d2ldmdv = function(y,mu,nu){-mu/(y-nu)^2},
               d2ldddv = function(y,sigma){rep(sigma,length(y))}, 
               G.dev.incr = function(y, mu, sigma, nu,...) -2 * dSwald(x=y, mu = mu, sigma = sigma, nu=nu,log = TRUE),
               rqres = expression(rqres(pfun = "pSwald",type = "Continuous", y = y, mu = mu,
                                        sigma = sigma, nu = nu)),
               mu.initial = expression(mu <- (y + mean(y))/2),
               sigma.initial = expression(sigma <- sd(y)/(mean(y))^1.5),
               nu.initial = expression({nu <- rep(min(y),length(y))}),
               mu.valid = function(mu) all(mu > 0),
               sigma.valid = function(sigma) all(sigma > 0),
               nu.valid = function(nu) all(nu > 0),
               y.valid = function(y) all(y > 0),
               mean = function(mu,sigma,nu) ((mu/sigma)+nu),
               variance = function(mu,sigma)((mu/sigma^3)),
               shift = function(nu) nu),
          class = c("gamlss.family", "family"))
}  
