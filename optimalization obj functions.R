obj.own <- function(par,y){
  par <- par
  y <- as.vector(y)
  if(par[3]<0){stop('threshold lower than y')}
  res <- par[1]/sqrt(2*pi*(y-par[3])^3)*exp(-(par[1]-par[2]*(y-par[3]))^2/(2*(y-par[3])))
  if(any(is.na(res))) res <- rep(0,length(y))
  -sum(log(res))
}

obj.seq <- function(par,y,w){
  par<- par
  y <- y
  w <- w
  res <- dSwald(x=y,mu=par[1],sigma=par[2],nu=par[3],log=FALSE)
  res.test <- w*res
  res <- res.test
  LL.1<- -sum(log(res))
  LL.1
}

obj <- function(par,y){
  par <- par
  y <- as.vector(y)
  # shifted wald (inverse Gauss) Likelihood function using supdist
  res <- statmod::dinvgauss(y-par[3], mean=par[1],shape=par[2])
  if(any(is.na(res))) res <- rep(10,length(y))
  -sum(log(res))
}
1