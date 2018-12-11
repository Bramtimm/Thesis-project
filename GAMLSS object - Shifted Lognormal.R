dSLOGNO <- function(x, mu = 0, sigma = 1, Ter=0.1, log = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(x < 0)){ 
    stop(paste("x must be >=0", "\n", ""))}
  fy <- dshifted_lnorm(x = x, meanlog = mu, sdlog = sigma, shift=Ter , log = log)
  fy
}

pSLOGNO <- function(q, mu = 0, sigma = 1,Ter=0, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be >=0", "\n", ""))
  cdf <- pshifted_lnorm(q, meanlog = mu, sdlog = sigma, shift = Ter,lower.tail=lower.tail,
                        log.p = log.p)
  cdf
}

qSLOGNO <- function (p, mu = 0, sigma = 1,Ter=0, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  q <- qshifted_lnorm(p = p, meanlog = mu, sdlog = sigma, shift = Ter, lower.tail = lower.tail, 
              log.p = log.p)
  q
}

rSLOGNO <- function (n, mu = 0, sigma = 1,Ter=0) 
{
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  r <- rshifted_lnorm(n = n, meanlog = mu, sdlog = sigma, shift = Ter)
}


SLOGNO <- function(mu.link="identity",sigma.link="log", Ter.link="identity")
{
  mstats <- checklink("mu.link",
                      "Shifted Log Normal", substitute(mu.link), 
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Shifted Log Normal", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  tstats <- checklink("Ter.link", "Shifted Log Normal", substitute(Ter.link),
                      c("inverse", "log", "identity", "own"))
  
structure(list(family = c("SLOGNO", "Shifted Log Normal"),
               parameters = list(mu=TRUE, sigma=TRUE, Ter=TRUE), nopar = 3,
               type = "continuous",
               mu.link = as.character(substitute(mu.link)),
               sigma.link = as.character(substitute(sigma.link)),
               Ter.link = as.character(substitute(Ter.link)),
               mu.linkfun = mstats$linkfun,
               sigma.linkfun = dstats$linkfun,
               Ter.linkfun = tstats$linkfun,
               mu.linkinv = mstats$linkinv,
               sigma.linkinv = dstats$linkinv,
               Ter.linkinv = tstats$linkinv,
               mu.dr = mstats$mu.eta,
               sigma.dr = dstats$mu.eta,
               Ter.dr = tstats$mu.eta,
               d1dm = function(y, mu, sigma,Ter){
                 d1dm <- (log(y-Ter) - mu)/sigma^2
                 d1dm
               },
               d2ldm2 = function(sigma) -1/sigma^2,
               dldd = function(y, mu, sigma,Ter) {
                 dldd <- (1/(sigma^3)) * ((log(y-Ter) - mu)^2 - sigma^2)
                 dldd
               },
               d2ldd2 = function(sigma) -2/sigma^2,
               d1dt = function(y,mu,sigma,Ter){
                 d1dt <- (1/(y-Ter)) + (1/(sigma^2)) * (1/(y-Ter)) * (log(y-Ter)-mu) 
                 d1dt
               },
               d2ldt2 = function(Ter)(-mu+log(y-Ter)+1)*(1/(sigma^2*(y-Ter)^2))+(1/((y-Ter)^2)),
               d2ldmdd = function(y) rep(0,length(y)),
               d2ldmdt = function(y,mu,sigma,Ter){-1/(sigma^2*(y-Ter))},
               d2ldddt = function(y,mu,sigma){(1/(sigma^3*(Ter-y)))*(2*(log(Ter-y)-mu))}, 
               G.dev.incr = function(y, mu, sigma, Ter,...) -2 * dSLOGNO(x = y, mu = mu, sigma = sigma, Ter=Ter,log = TRUE),
               rqres = expression(rqres(pfun = "pSLOGNO",type = "Continuous", y = y, mu = mu,
                                        sigma = sigma, Ter = Ter)),
               mu.initial = expression({mu <- (log(y) + mean(log(y)))/2}),
               sigma.initial = expression({sigma <- rep(sd(log(y)), length(y))
               }),
               Ter.initial = expression({Ter <- rep(0,length(y))}),
               mu.valid = function(mu) all(mu > 0),
               sigma.valid = function(sigma) all(sigma > 0),
               Ter.valid = function(Ter) all(Ter >= 0),
               y.valid = function(y) all(y > 0),
               mean = function(mu,sigma,Ter) exp((mu + sigma^2)/2)+Ter,
               variance = function(mu,sigma) exp(2 * mu + sigma^1) * (exp(sigma^2 + 2))),
               shift = function(Ter) Ter,
               class = c("gamlss.family", "family"))
               }  
               
