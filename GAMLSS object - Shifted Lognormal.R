dSLOGNO <- function(x, mu = 0, sigma = 1, nu=0.1, log = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(x < 0)){ 
    stop(paste("x must be >=0", "\n", ""))}
  fy <- dshifted_lnorm(x = x, meanlog = mu, sdlog = sigma, shift=nu , log = log)
  fy
}

pSLOGNO <- function(q, mu = 0, sigma = 1,nu =0, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be >=0", "\n", ""))
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
  mstats <- checklink("mu.link",
                      "Shifted Log Normal", substitute(mu.link), 
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Shifted Log Normal", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "Shifted Log Normal", substitute(nu.link),
                      c("inverse", "log", "identity", "own"))
  
structure(list(family = c("SLOGNO", "Shifted Log Normal"),
               parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), nopar = 3,
               type = "continuous",
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
               d1dm = function(y,mu,sigma,nu){
                 d1dm <- (log(y-nu) - mu)/sigma^2
                 d1dm
               },
               d2ldm2 = function(sigma) -1/sigma^2,
               dldd = function(y,mu,sigma,nu){
                 dldd <- (1/(sigma^3)) * ((log(y-nu) - mu)^2 - sigma^2)
                 dldd
               },
               d2ldd2 = function(sigma) -2/sigma^2,
               d1dv = function(y,mu,sigma,nu){
                 d1dv <- (1/(y-nu)) + (1/(sigma^2)) * (1/(y-nu)) * (log(y-nu)-mu) 
                 d1dv
               },
               d2ldv2 = function(mu,y,nu,sigma){
                 d21dv2 <- (-mu+log(y-nu)+1)*(1/(sigma^2*(y-nu)^2))+(1/((y-nu)^2))},
               d2ldmdd = function(y) rep(0,length(y)),
               d2ldmdv = function(y,mu,sigma,nu){-1/(sigma^2*(y-nu))},
               d2ldddv = function(y,mu,sigma,nu){(1/(sigma^3*(nu-y)))*(2*(log(nu-y)-mu))}, 
               G.dev.incr = function(y, mu, sigma, nu,...) -2 * dSLOGNO(y, mu = mu, sigma = sigma, nu=nu,log = TRUE),
               rqres = expression(rqres(pfun = "pSLOGNO",type = "Continuous", y = y, mu = mu,
                                        sigma = sigma, nu = nu)),
               mu.initial = expression({mu <- (log(y) + mean(log(y)))/2}),
               sigma.initial = expression({sigma <- rep(sd(log(y)), length(y))
               }),
               nu.initial = expression({nu <- rep(0,length(y))}),
               mu.valid = function(mu) all(mu > 0),
               sigma.valid = function(sigma) all(sigma > 0),
               nu.valid = function(nu,y,mu) all(nu >= 0 & nu < y & nu < mu),
               y.valid = function(y) all(y > 0),
               mean = function(mu,sigma,nu) exp((mu + sigma^2)/2)+nu,
               variance = function(mu,sigma) exp(2 * mu + sigma^1) * (exp(sigma^2 + 2))),
               shift = function(nu) nu,
               class = c("gamlss.family", "family"))
               }  
               
