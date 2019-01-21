obj <- function(par,y,w){
  par <- par
  y <- y
  w <- w
  # shifted wald (inverse Gauss) Likelihood function using supdist
  res <- w*dSLOGNO(y, meanlog=par[1], sdlog = par[2], shift=par[3], log = FALSE)
  if(any(is.na(res))) res <- rep(0,length(y))
  -sum(log(res))
}