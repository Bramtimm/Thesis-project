dSLOGNO <- function(x, mu = mu, sigma = sigma, nu=nu, log = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(x < 0)){ 
    stop(paste("x must be >=0", "\n", ""))}
  if(any(nu > x & nu<0)){
    stop(paste("nu must be smaller than x and higher than zero"))
  }
  fy <- dshifted_lnorm(x = x, meanlog = mu, sdlog = sigma, shift=nu , log = log)
  fy
}

pSLOGNO <- function(q, mu = mu, sigma = sigma,nu =nu, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be >=0", "\n", ""))
  if(any(nu<0)){
    stop(paste("nu must be smaller than x and higher than zero"))
  }
  cdf <- pshifted_lnorm(q, meanlog = mu, sdlog = sigma, shift = nu,lower.tail=lower.tail,
                        log.p = log.p)
  cdf
}

qSLOGNO <- function(p, mu = 0, sigma = 1,nu=0, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  q <- qshifted_lnorm(p = p, meanlog = mu, sdlog = sigma, shift = nu, lower.tail = lower.tail, 
              log.p = log.p)
  q
}

rSLOGNO <- function(n, mu = 0, sigma = 1,nu=0) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  r <- rshifted_lnorm(n = n, meanlog = mu, sdlog = sigma, shift = nu)
}


SLOGNO <- function(mu.link="identity",sigma.link="log", nu.link="identity")
{
  mstats <- checklink("mu.link", "Shifted Log Normal", substitute(mu.link), 
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Shifted Log Normal", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "Shifted Log Normal", substitute(nu.link),
                      c("inverse", "log", "identity", "own"))
  
structure(list(family = c("SLOGNO", "Shifted Log Normal"),
               parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), nopar = 3,
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
               dldm = function(y,mu,sigma,nu){
                 dldm <- (log(y-nu) - mu)*(1/(sigma^2))
                 dldm
               },
               d2ldm2 = function(sigma){
                 d2ldm2 <- -(1/sigma^2)
                 d2ldm2
               },
               dldd = function(y,mu,sigma,nu){
                 dldd <- (1/(sigma^3)) * ((log(y-nu) - mu)^2 - sigma^2)
                 dldd
               },
               d2ldd2 = function(y,mu,sigma,nu){
                 d2ldd2 <- (sigma^2-3*((log(y-nu)-mu)^2))*(1/sigma^4)
                 d2ldd2
               },
               dldv = function(y,mu,sigma,nu){
                 dldv <- (-mu+sigma^2+log(y-nu))*(1/(sigma^2*(y-nu))) 
                 dldv
               },
               d2ldv2 = function(y,mu,sigma,nu){
                 d2ldv2 <- (mu-sigma^2-log(y-nu)-1)*(1/(sigma^2*(y-nu)^2))
                 d2ldv2
               },
               d2ldmdd = function(y,mu,sigma,nu) (2*(mu-log(y-nu)))*(1/(sigma^3)),
               d2ldmdv = function(y,sigma,nu){-(1/(sigma^2*(y-nu)))},
               d2ldddv = function(y,mu,sigma,nu){(2*log(y-nu)-mu)*1/(sigma^3*(nu-y))}, 
               G.dev.incr = function(y, mu, sigma, nu,...) -2 * dSLOGNO(y, mu = mu, sigma = sigma, nu=nu,log = TRUE),
               rqres = expression(rqres(pfun = "pSLOGNO",type = "Continuous", y = y, mu = mu,
                                        sigma = sigma, nu = nu)),
               mu.initial = expression({mu <- (log(y) + mean(log(y)))/2}),
               sigma.initial = expression({sigma <- rep(sd(log(y)), length(y))
               }),
               nu.initial = expression({nu <- rep(0,length(y))}),
               mu.valid = function(mu) TRUE,
               sigma.valid = function(sigma) all(sigma > 0),
               nu.valid = function(nu,y) all(nu>=0 & nu<y),
               y.valid = function(y) all(y > 0),
               mean = function(mu,sigma,nu) exp((mu + sigma^2)/2)+nu,
               variance = function(mu,sigma) (exp(2 * mu + sigma^2) * (exp(sigma^2)-1)),
               shift = function(nu) nu),
               class = c("gamlss.family", "family"))
               }  
               
