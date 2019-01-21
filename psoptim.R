obj <- function(pars, rt) {
  
  meanlog <- pars[1] 
  sdlog <- pars[2]
  shift <- pars[3]

  densities <- tryCatch(
    dSLOGNO(rt, meanlog=meanlog, 
            sdlog=sdlog, shift=shift))
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

obj(c(-1,3,0),RT[1:10])

dSLOGNO(RT,-1,3,0)==0
debug(obj)
